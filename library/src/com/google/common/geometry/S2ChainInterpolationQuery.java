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

import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * S2ChainInterpolationQuery is a helper class for querying points along an S2Shape's edges (chains
 * of vertices) by spherical angular distances. The distance, measured in radians, is computed by
 * accumulating the lengths of the edges of the shape, in the order the edges are stored by the
 * S2Shape object.
 *
 * <p>If a particular edge chain is specified at the query initialization, then the distances are
 * along that single chain, which allows per-chain operations. If no chain is specified, then the
 * interpolated points as a function of distance will have discontinuities at chain boundaries. This
 * makes it easier to implement some algorithms, such as random sampling along the total length of a
 * multiline.
 *
 * <p>Once the query object is initialized, the complexity of each subsequent query is O(log(m)),
 * where m is the number of edges. The complexity of the constructor and the memory footprint of the
 * query object are both O(m).
 */
public class S2ChainInterpolationQuery {
  /** The shape being interpolated. */
  private final S2Shape shape;

  /**
   * Cumulative value 'n' is the total length of the edges from the beginning of the first edge to
   * the end of edge 'n', in radians.
   */
  private final ArrayList<S1Angle> cumulativeValues;

  /** The first edge id of the chain being interpolated. */
  private final int firstEdgeId;

  /** The last edge id of the chain being interpolated. */
  private final int lastEdgeId;

  // resultPoint, resultEdgeId, and resultDistance are set by findPoint() or findPointAtFraction().
  private S2Point resultPoint;
  private int resultEdgeId;
  private S1Angle resultDistance;

  /** Constructs an S2ChainInterpolationQuery for all the edge chains of the given shape. */
  public S2ChainInterpolationQuery(S2Shape shape) {
    this(shape, -1);
  }

  /**
   * Constructs an S2ChainInterpolationQuery for a specific chain of the given shape (or all edges).
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
    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
    for (int i = firstEdgeId; i <= lastEdgeId; ++i) {
      cumulativeValues.add(cumulativeAngle);
      shape.getEdge(i, edge);
      cumulativeAngle = cumulativeAngle.add(new S1Angle(edge.a, edge.b));
    }
    if (!cumulativeValues.isEmpty()) {
      cumulativeValues.add(cumulativeAngle);
    }
  }

  /**
   * Gets the maximum absolute accumulated length, from the beginning of the first selected edge to
   * the end of the last selected edge. Returns zero for shapes containing no edges, or zero-length
   * chains.
   */
  public S1Angle getLength() {
    // The total length equals the cumulative value at the end of the last edge, iff there is at
    // least one edge in the shape.
    return cumulativeValues.isEmpty() ? S1Angle.ZERO : Iterables.getLast(cumulativeValues);
  }

  /**
   * Returns the cumulative length along the edges being interpolated, up to the end of the given
   * edge id. Returns S1Angle.INFINITY if the given edge id does not lie within the set of edges
   * being interpolated. Returns S1Angle.ZERO if the S2ChainInterpolationQuery is empty.
   */
  public S1Angle getLengthAtEdgeEnd(int edgeId) {
    if (cumulativeValues.isEmpty()) {
      return S1Angle.ZERO;
    }

    if (edgeId < firstEdgeId || edgeId > lastEdgeId) {
      return S1Angle.INFINITY;
    }

    return cumulativeValues.get(edgeId - firstEdgeId + 1);
  }

  /** Returns a slice of the chain from 'a' to 'b', in reverse order if {@code b > a}. */
  public List<S2Point> slice(double a, double b) {
    List<S2Point> result = new ArrayList<>();
    addSlice(result, a, b);
    return result;
  }

  /**
   * Appends the chain sliced from 'a' to 'b' to 'result', in reverse order if {@code b > a}. If the
   * query is empty, 'result' is unchanged and 'a' and 'b' are unused. Otherwise, 'a' and 'b' must
   * be between 0 and 1, inclusive.
   */
  public void addSlice(List<S2Point> result, double a, double b) {
    if (cumulativeValues.isEmpty()) {
      return;
    }

    int start = result.size();
    boolean reverse = a > b;
    if (reverse) {
      double t = a;
      a = b;
      b = t;
    }
    Preconditions.checkArgument(findPointAtFraction(a), "Invalid value of A %s", a);
    int startEdge = resultEdgeId();
    S2Point last = resultPoint();
    result.add(last);
    Preconditions.checkArgument(findPointAtFraction(b), "Invalid value of B %s", b);
    int endEdge = resultEdgeId();
    MutableEdge edge = new MutableEdge();
    for (int id = startEdge; id < endEdge; id++) {
      shape.getEdge(id, edge);
      if (!last.equalsPoint(edge.b)) {
        last = edge.b;
        result.add(edge.b);
      }
    }
    result.add(resultPoint());
    if (reverse) {
      Collections.reverse(result.subList(start, result.size()));
    }
  }

  /**
   * Computes the S2Point located at the given distance along the edges from the first vertex of the
   * first edge. Also computes the edge id and the actual normalized distance corresponding to the
   * resulting point. The values computed during the last findPoint() or findPointAtFraction() call
   * are accessible via the respective result*() methods, and are valid iff the last such call
   * returned true.
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
    // The input distance must then be on the edge ending at the corresponding vertex. If every
    // cumulative value is smaller, lowerBound returns cumulativeValues.size().
    final int lowerBound =
        S2ShapeUtil.lowerBound(
            0, cumulativeValues.size(), i -> cumulativeValues.get(i).compareTo(distance) < 0);

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
          S2EdgeUtil.getPointOnLine(
              shapeEdge.getStart(),
              shapeEdge.getEnd(),
              distance.sub(cumulativeValues.get(lowerBound - 1)));
    }

    return true;
  }

  /**
   * Computes the S2Point located at the given normalized fraction along the edges from the first
   * vertex of the first edge. A fraction = 0 corresponds to the beginning of the shape or chain,
   * and fraction = 1 to the end. Also computes the edge id and the actual normalized distance
   * corresponding to the resulting point. The values computed during the last findPoint() or
   * findPointAtFraction() call are accessible via the respective result*() methods, and are valid
   * iff the last such call returned true.
   *
   * <p>This method returns true iff the query has been initialized with at least one edge.
   *
   * <p>Forwards the call to findPoint(S1Angle distance). A small precision loss may occur due to
   * converting the fraction to a distance by multiplying it by the total length.
   */
  public boolean findPointAtFraction(double fraction) {
    return findPoint(getLength().mul(fraction));
  }

  /**
   * Returns the point which is the result of the last query. The point is valid iff the last
   * findPoint*() call returned true.
   */
  public S2Point resultPoint() {
    return resultPoint;
  }

  /**
   * Returns the index of the edge on which the point from the last query is located. The edge id is
   * valid iff the last findPoint*() call returned true.
   */
  public int resultEdgeId() {
    return resultEdgeId;
  }

  /**
   * Returns the actual distance of the resulting point from the last query. It may differ from the
   * input distance if the latter is negative or exceeds the total length of the shape's chain(s),
   * in which case the resulting distance is snapped to 0 or total length, respectively. The
   * distance value is valid iff the last findPoint*() call returned true.
   */
  public S1Angle resultDistance() {
    return resultDistance;
  }
}
