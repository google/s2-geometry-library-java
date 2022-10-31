/*
 * Copyright 2022 Google Inc.
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

import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.List;
import java.util.Optional;

/**
 * S2HausdorffDistanceQuery is a helper class for computing discrete Hausdorff distances between two
 * geometries. This can be useful for e.g. computing how far apart a highway and a semi-parallel
 * frontage road ever get.
 *
 * <p>Both geometries are provided as S2ShapeIndex objects. A S2ShapeIndex is a collection of
 * points, polylines, and/or polygons. See S2ShapeIndex.java for details.
 *
 * <p>Discrete directed Hausdorff distance from target geometry A to source geometry B is defined as
 * the maximum, over all vertices of A, of the closest edge distance from the vertex to geometry B.
 * It is called _discrete_ because the maximum is computed over all _vertices_ of A, rather than
 * over other positions such as midpoints of the edges. The current implementation computes the
 * discrete Hausdorff distance instead of exact Hausdorff distance because the latter incurs
 * significantly larger computational expenses, while the former is suitable for most practical use
 * cases.
 *
 * <p>The undirected Hausdorff distance (usually referred to more simply as just Hausdorff distance)
 * between geometries A and B is the maximum of the directed Hausdorff distances from A to B, and
 * from B to A.
 *
 * <p>The difference between directed and undirected Hausdorff distances can be illustrated by the
 * following example. Let's say we have two polygonal geometries, one representing the continental
 * contiguous United States, and the other, Catalina island off the coast of California. The
 * directed Hausdorff distance from the island to the continental US is going to be about 60 km -
 * this is how far the furthest point on the island gets from the continent. At the same time, the
 * directed Hausdorff distance from the continental US to Catalina island is several thousand
 * kilometers - this is how far from the island one can get at the North-East corner of the US. The
 * undirected Hausdorff distance between these two entities is the maximum of the two directed
 * distances, that is a few thousand kilometers.
 *
 * <p>For example, given two geometries A and B in the form of S2 shape indexes, the following code
 * finds the directed Hausdorff distance from A to B and the [undirected] Hausdorff distance between
 * A and B:
 *
 * <pre>{@code
 * boolean computeHausdorffDistances(
 *     S2ShapeIndex A,
 *     S2ShapeIndex B,
 *     S1ChordAngle directedDistance,
 *     S1ChordAngle& undirectedDistance) {
 *   S2HausdorffDistanceQuery query = S2HausdorffDistanceQuery.builder().build();
 *
 *   Optional<DirectedResult> directedResult = query.getDirectedHausdorffDistance(A, B);
 *   if (!directedResult.isPresent()) {
 *     return false;
 *   }
 *   directedDistance = directedResult.get().distance();
 *
 *   Optional<Result> undirectedResult = query.getHausdorffDistance(A, B);
 *   undirectedDistance = undirectedResult.get().distance();
 *
 *   return true;
 * }
 * }</pre>
 *
 * <p>For the definition of Hausdorff distance and other details see
 * https://en.wikipedia.org/wiki/HausdorffDistance.
 */
@CheckReturnValue
public final class S2HausdorffDistanceQuery {
  // Options are constructed from a Builder.
  private final Options options;

  // maxDistance is a temporary value, repeatedly updated as a query runs with the largest target-
  // to-source distance found so far.
  private final DistanceCollector<S1ChordAngle> maxDistance = S1ChordAngle.maxCollector();

  // sourcePoint is a temporary value. When maxDistance is updated, sourcePoint is also updated to
  // the corresponding point on the source.
  private S2Point sourcePoint = null;

  // targetPoint is a temporary value. When maxDistance is updated, targetPoint is also updated to
  // the corresponding point on the target.
  private S2Point targetPoint = null;

  /** Options control the set of edges returned by S2HausdorffDistanceQuery. */
  public static class Options {
    // Are polygon interiors included when considering distance?
    protected final boolean includeInteriors;

    /** Internal constructor from a Builder. */
    Options(Builder b) {
      includeInteriors = b.includeInteriors();
    }

    /** Returns a new Builder with values copied from these Options. */
    public Builder toBuilder() {
      Builder b = new Builder();
      b.setIncludeInteriors(includeInteriors());
      return b;
    }

    /** True if polygon interiors in the queried index are considered when computing distance. */
    public boolean includeInteriors() {
      return includeInteriors;
    }
  }

  /** The Builder for S2HausdorffDistanceQuery and its Options. */
  public static class Builder {
    private boolean includeInteriors;

    /** Constructs a new Builder with default values. */
    public Builder() {
      includeInteriors = true;
    }

    /**
     * The includeInteriors flag (true by default) indicates that the distance should be computed
     * not just to polygon boundaries of the source index, but to polygon interiors as well. For
     * example, if target shape A is fully contained inside the source shape B, and
     * includeInteriors is set to true, then the directed Hausdorff distance from A to B is zero.
     */
    @CanIgnoreReturnValue
    public Builder setIncludeInteriors(boolean includeInteriors) {
      this.includeInteriors = includeInteriors;
      return this;
    }

    /** Returns the current value of the includeInteriors option. */
    public boolean includeInteriors() {
      return includeInteriors;
    }

    /**
     * Returns a new S2HausdorffDistanceQuery with options set from current values of this
     * Builder.
     */
    public S2HausdorffDistanceQuery build() {
      return new S2HausdorffDistanceQuery(new Options(this));
    }
  }

  /** DirectedResult stores the results of directed Hausdorff distance queries. */
  public static class DirectedResult {
    private final S1ChordAngle distance;
    private final S2Point targetPoint;

    /* Constructs a new DirectedResult with the given values. */
    public DirectedResult(S1ChordAngle distance, S2Point targetPoint) {
      this.distance = distance;
      this.targetPoint = targetPoint;
    }

    /** Returns the resulting directed Hausdorff distance value. */
    public S1ChordAngle distance() {
      return distance;
    }

    /**
     * Returns the point on the target index on which the directed Hausdorff distance is achieved.
     */
    public S2Point targetPoint() {
      return targetPoint;
    }
  }

  /**
   * Result stores the output of an undirected Hausdorff distance query. It consists of two directed
   * query results, forward and reverse.
   */
  public static class Result {
    private final DirectedResult targetToSource;
    private final DirectedResult sourceToTarget;

    /**
     * Constructs a new Result with the directed result of the forward query (target index to source
     * index) and that of the reverse (source index to target index) query.
     */
    public Result(DirectedResult targetToSource, DirectedResult sourceToTarget) {
      this.targetToSource = targetToSource;
      this.sourceToTarget = sourceToTarget;
    }

    /**
     * Returns the actual (undirected) Hausdorff distance, which is the maximum of the two directed
     * query results.
     */
    public S1ChordAngle distance() {
      return S1ChordAngle.max(targetToSource.distance(), sourceToTarget.distance());
    }

    /** Returns the result for the target-to-source directed Hausdorff distance call. */
    public DirectedResult targetToSource() {
      return targetToSource;
    }

    /** Returns the result for the source-to-target directed Hausdorff distance call. */
    public DirectedResult sourceToTarget() {
      return sourceToTarget;
    }
  }

  /**
   * Internal S2HausdorffDistanceQuery constructor from provided Options. Clients should use
   * {@code S2HausdorffDistanceQuery query = S2HausdorffDistanceQuery.builder().build();}
   */
  S2HausdorffDistanceQuery(Options options) {
    this.options = options;
  }

  /** Returns a new Builder for S2HausdorffDistanceQuery with default options. */
  public static Builder builder() {
    return new Builder();
  }

  /** Returns the immutable current Options on this S2HausdorffDistanceQuery. */
  public Options options() {
    return options;
  }

  /**
   * Compute directed Hausdorff distance from the target index to the source index. Returns
   * Optional.empty() iff at least one of the shape indexes is empty.
   *
   * <p>Note that directed Hausdorff distance from geometry A (as target) to geometry B (as source)
   * is not (in general case) equal to that from B (as target) to A (as source).
   */
  public Optional<DirectedResult> getDirectedResult(S2ShapeIndex target, S2ShapeIndex source) {
    S2ClosestEdgeQuery.Query closestEdgeQuery = S2ClosestEdgeQuery.builder()
        .setMaxResults(1)
        .setIncludeInteriors(options.includeInteriors())
        .build(source);

    // maxDistance is repeatedly updated with the largest target-to-source distance found so far.
    maxDistance.reset();
    // When maxDistance is updated, sourcePoint is updated to the corresponding point on the source.
    sourcePoint = null;
    // When maxDistance is updated, targetPoint is updated to the corresponding point on the target.
    targetPoint = null;

    // This approximation of Haussdorff distance is based on computing closest point distances from
    // the _vertices_ of the target index to the _edges_ of the source index. Hence we iterate over
    // all shapes in the target index, then over all chains in those shapes, then over all edges in
    // those chains, and then over the edges' vertices.
    for (S2Shape shape : target.getShapes()) {
      for (List<S2Point> chain : shape.chains()) {
        for (S2Point vertex : chain) {
          updateMaxDistance(vertex, closestEdgeQuery);
        }
      }
    }

    if (maxDistance.distance().lessThan(S1ChordAngle.ZERO)) {
      // maxDistance was never updated.
      return Optional.empty();
    } else {
      return Optional.of(new DirectedResult(maxDistance.distance(), targetPoint));
    }
  }

  /**
   * Same as {@link #getDirectedResult(S2ShapeIndex, S2ShapeIndex)}, except returns the actual
   * distance, or S1ChordAngle.INFINITY iff at least one of the shape indexes is empty.
   */
  public S1ChordAngle getDirectedDistance(S2ShapeIndex target, S2ShapeIndex source) {
    Optional<DirectedResult> directedResult = getDirectedResult(target, source);
    return directedResult.isPresent() ? directedResult.get().distance() : S1ChordAngle.INFINITY;
  }

  /**
   * Compute the undirected Hausdorff distance between the target index and the source index.
   * Returns Optional.empty() iff at least one of the shape indexes is empty.
   *
   * <p>Note that the result of this query is symmetrical with respect to target vs. source, i.e. if
   * target and source indices are swapped, the resulting Hausdorff distance remains unchanged.
   */
  public Optional<Result> getResult(S2ShapeIndex target, S2ShapeIndex source) {
    Optional<DirectedResult> targetToSource = getDirectedResult(target, source);
    if (targetToSource.isPresent()) {
      return Optional.of(new Result(targetToSource.get(), getDirectedResult(source, target).get()));
    } else {
      return Optional.empty();
    }
  }

  /**
   * Same as {@link #getResult(S2ShapeIndex, S2ShapeIndex)}, but returns the maximum of forward and
   * reverse distances, or S1ChordAngle.Infinity() iff at least one of the shape indexes is empty.
   */
  public S1ChordAngle getDistance(S2ShapeIndex target, S2ShapeIndex source) {
    Optional<Result> result = getResult(target, source);
    return result.isPresent() ? result.get().distance() : S1ChordAngle.INFINITY;
  }

  /**
   * This internally used function computes the closest edge distance from the given target "point"
   * to the source index using the given S2ClosestEdgeQuery. If necessary, updates the maxDistance,
   * the targetPoint, and the sourcePoint.
   */
  private void updateMaxDistance(S2Point point, S2ClosestEdgeQuery.Query closestEdgeQuery) {
    // In case we already have a valid result, it can be used as the lower bound estimate for the
    // final Hausdorff distance. Therefore, if the distance between the current target point and the
    // last source point does not exceed this lower bound, we can safely skip this target point, not
    // updating the maximum distance.
    if (maxDistance.distance().greaterOrEquals(S1ChordAngle.ZERO)
        && S2Predicates.compareDistance(
                point, sourcePoint, maxDistance.distance().getLength2()) <= 0) {
      return;
    }

    // Otherwise, find the single closest edge and the corresponding closest point in the source
    // geometry to the target point.
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.PointTarget<>(point);
    closestEdgeQuery.findClosestEdges(target, (distance, shape, edgeId) -> {
      // If this closest distance is greater than the previous maximum, update it and the
      // corresponding source and target points.
      if (maxDistance.update(distance)) {
        targetPoint = point;

        // If the edgeId is -1, the target point is contained by the source shape. Otherwise,
        // find the closest point on the closest source (shape, edge) to the target point.
        if (edgeId == -1) {
          sourcePoint = targetPoint;
        } else {
          S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
          shape.getEdge(edgeId, edge);
          sourcePoint = S2EdgeUtil.getClosestPoint(targetPoint, edge.getStart(), edge.getEnd());
        }
      }
      return true;
    });
  }
}
