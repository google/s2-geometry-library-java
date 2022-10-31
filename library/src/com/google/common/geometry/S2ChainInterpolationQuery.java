/*
 * Copyright 2021 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.google.common.geometry;

import static com.google.common.base.Preconditions.checkArgument;

import com.google.common.collect.Iterables;
import java.util.ArrayList;

/**
 * S2ChainInterpolationQuery is a helper class for querying points along an S2Shape's edges (chains
 * of vertices) by spherical angular distances.  The distance, measured in radians, is computed
 * accumulatively along the edges contained in the shape, using the order in which the edges are
 * stored by the S2Shape object.
 *
 * <p>If a particular edge chain is specified at the query initialization, then the distances are
 * along that single chain, which allows per-chain operations.  If no chain is specified, then the
 * interpolated points as a function of distance will have discontinuities at chain boundaries.
 * This makes it easier to implement some algorithms, such as random sampling along the total length
 * of a multiline.
 *
 * <p>Once the query object is initialized, the complexity of each subsequent query is O(log(m)),
 * where m is the number of edges.  The complexity of the constructor and the memory footprint of
 * the query object are both O(m).
 */
public class S2ChainInterpolationQuery {
  private final S2Shape shape;
  private final ArrayList<S1Angle> cumulativeValues;
  private final int firstEdgeId;
  private final int lastEdgeId;

  private S2Point resultPoint;
  private int resultEdgeId;
  private S1Angle resultDistance;

  /**
   * Constructs an S2ChainInterpolationQuery for all the edge chains of the given shape.
   */
  public S2ChainInterpolationQuery(S2Shape shape) {
    this(shape, -1);
  }

  /**
   * Constructs an S2ChainInterpolationQuery for the given shape.
   *
   * <p>If a non-negative chainId is supplied by the caller, then only the edges belonging to the
   * chain with that chainId are used. Otherwise the edges from all chains are used in the order
   * they are stored in the shape object.
   */
  public S2ChainInterpolationQuery(S2Shape shape, int chainId) {
    checkArgument(chainId < shape.numChains());

    this.shape = shape;
    if (chainId >= 0) {
      // If a valid chain id was provided, then the range of edge ids is defined by the start and
      // the length of the chain.
      firstEdgeId = shape.getChainStart(chainId);
      lastEdgeId = firstEdgeId + shape.getChainLength(chainId) - 1;
    } else {
      // If no valid chain id was provided then we use the whole range of the shape's edge ids.
      firstEdgeId = 0;
      lastEdgeId = shape.numEdges() - 1;
    }

    cumulativeValues = new ArrayList<>();
    S1Angle cumulativeAngle = S1Angle.ZERO;
    for (int i = firstEdgeId; i <= lastEdgeId; ++i) {
      cumulativeValues.add(cumulativeAngle);
      S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
      shape.getEdge(i, edge);
      cumulativeAngle = cumulativeAngle.add(new S1Angle(edge.a, edge.b));
    }
    if (!cumulativeValues.isEmpty()) {
      cumulativeValues.add(cumulativeAngle);
    }
  }


  /**
   * Gets the maximum absolute accumulated length, which corresponds to the end vertex of the last
   * selected edge.  Returns zero for shapes containing no edges.
   */
  public S1Angle getLength() {
    // The total length equals the cumulative value at the end of the last edge, iff there is at
    // least one edge in the shape.
    return cumulativeValues.isEmpty() ? S1Angle.ZERO : Iterables.getLast(cumulativeValues);
  }

  /**
   * Computes the S2Point located at the given distance along the edges from the first vertex of the
   * first edge.  Also computes the edge id and the actual normalized distance corresponding to the
   * resulting point.  The values computed during the last findPoint() call are accessible via the
   * respective result*() methods, and are valid iff the last findPoint() call returned true.
   *
   * <p>This method returns true iff the query has been initialized with at least one edge.
   *
   * <p>If the input distance exceeds the total length, then the resulting point is the end vertex
   * of the last edge, and the resulting distance is set to the total length. If the input distance
   * is negative, then the resulting point is the first vertex of the first edge, and the resulting
   * distance is set to 0.
   *
   * <p>If there are one or more degenerate (zero-length) edges corresponding to the given distance,
   * then the resulting point is located on the first of these edges.
   */
  public boolean findPoint(S1Angle distance) {
    // Return false for uninitialized queries or for shapes containing no edges.
    if (cumulativeValues.isEmpty()) {
      return false;
    }

    // Binary search in the list of cumulative values, which by construction are sorted in ascending
    // order, to obtain the lowest cumulative value that is not smaller than the target distance.
    // The input distance must then be on the edge ending at the corresponding vertex.  If every
    // cumulative value is smaller, lowerBound returns cumulativeValues.size().
    final int lowerBound = S2ShapeUtil.lowerBound(0, cumulativeValues.size(),
        i -> cumulativeValues.get(i).compareTo(distance) < 0);

    S2Shape.MutableEdge shapeEdge = new S2Shape.MutableEdge();
    if (lowerBound == 0) {
      // Corner case: the first vertex of the shape at distance = 0.
      shape.getEdge(firstEdgeId, shapeEdge);
      resultPoint = shapeEdge.getStart();
      resultEdgeId = firstEdgeId;
      resultDistance = cumulativeValues.get(0);
    } else if (lowerBound == cumulativeValues.size()) {
      // Corner case: the input distance is greater than the total length, hence we snap the result
      // to the last vertex of the shape at distance = total length.
      shape.getEdge(lastEdgeId, shapeEdge);
      resultPoint = shapeEdge.getEnd();
      resultEdgeId = lastEdgeId;
      resultDistance = Iterables.getLast(cumulativeValues);
    } else {
      // Obtain the edge index and compute the interpolated result from the edge vertices.
      resultEdgeId = lowerBound + firstEdgeId - 1;
      shape.getEdge(resultEdgeId, shapeEdge);
      resultDistance = distance;
      // Interpolate along shapeEdge by the distance beyond the cumulative distance at the starting
      // point of shapeEdge.
      resultPoint =
          S2EdgeUtil.interpolateAtDistance(
              distance.sub(cumulativeValues.get(lowerBound - 1)),
              shapeEdge.getStart(),
              shapeEdge.getEnd());
    }

    return true;
  }

  /**
   * Similar to the above function.  Takes a normalized fraction of the distance as input, with
   * fraction = 0 corresponding to the beginning of the shape or chain and fraction = 1 to the end.
   * Forwards the call to findPoint(S1Angle distance).  A small precision loss may occur due to
   * converting the fraction to a distance by multiplying it by the total length.
   */
  public boolean findPointAtFraction(double fraction) {
    return findPoint(getLength().mul(fraction));
  }

  /**
   * Returns the point which is the result of the last query.  The point is valid iff the last
   * findPoint*() call returned true.
   */
  public S2Point resultPoint() {
    return resultPoint;
  }

  /**
   * Returns the index of the edge on which the point from the last query is located.  The edge id
   * is valid iff the last findPoint*() call returned true.
   */
  public int resultEdgeId() {
    return resultEdgeId;
  }

  /**
   * Returns the actual distance of the resulting point from the last query.  It may differ from the
   * input distance if the latter is negative or exceeds the total length of the shape's chain(s),
   * in which case the resulting distance is snapped to 0 or total length, respectively.  The
   * distance value is valid iff the last findPoint*() call returned true.
   */
  public S1Angle resultDistance() {
    return resultDistance;
  }

}
