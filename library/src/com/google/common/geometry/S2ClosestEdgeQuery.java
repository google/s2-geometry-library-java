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

import static com.google.common.base.Preconditions.checkArgument;
import static java.util.Comparator.naturalOrder;

import com.google.common.annotations.VisibleForTesting;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

/**
 * {@link S2ClosestEdgeQuery} is a helper class for searching within an {@link S2ShapeIndex} to find
 * the closest edge(s) to a given target point, edge, s2 cell, or {@link S2ShapeIndex} geometry
 * collection. Closest edges are those with minimum distance from any point on that edge to any
 * point on the target geometry.
 *
 * <p>For example, given a set of polylines, the following code efficiently finds up to the closest
 * 100 edges of those polylines to a query target point, and at most 5km from it. Then, the single
 * polyline edge closest to any edge in a second target ShapeIndex is found.
 *
 * <p>This example code demonstrates the use of the provided convenience methods for working with
 * S2ClosestEdgeQuery using an S1Distance abstract type that is specifically S1ChordAngle. This is
 * recommended for most clients.
 *
 * <pre>{@code

class ClosestEdgeDemo {
  // Recommended: S2Earth.toAngle(new Length.Kilometers(5));
  static final S1Angle MAX_DISTANCE = S1Angle.fromEarthDistance(5000);

  public static void example(
      List<S2Polyline> polylines, S2ShapeIndex targetIndex, S2Point targetPoint) {
    // Add the given polylines to a shape index.
    S2ShapeIndex polylinesIndex = new S2ShapeIndex();
    for (S2Polyline polyline : polylines) {
      polylinesIndex.add(S2LaxPolylineShape.create(polyline));
    }

    // Create a Builder, modifying the default values with a maximum distance and number of
    // results.
    S2ClosestEdgeQuery.Builder builder =
        S2ClosestEdgeQuery.builder().setInclusiveMaxDistance(MAX_DISTANCE).setMaxResults(100);

    // Construct the query.
    S2ClosestEdgeQuery.Query query = builder.build(polylinesIndex);

    // Create a target that is a single point.
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> pointTarget =
        new S2ClosestEdgeQuery.PointTarget<>(targetPoint);

    // Find and visit up to 100 polyline edges in polylinesIndex that are the closest to the
    // target point, and at most 5 km from it.
    query.findClosestEdges(
        pointTarget,
        (S1ChordAngle distance, S2Shape shape, int edgeId) -> {
          // Do something with each of these closest edge results.
          return true;
        });

    // Create a target that is another shapeIndex.
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> shapeIndexTarget =
        S2ClosestEdgeQuery.createShapeIndexTarget(targetIndex);
    shapeIndexTarget.setIncludeInteriors(true);

    // Find the single closest polyline edge to any edge or polygon interior in the
    // shapeIndexTarget, if any matches the query options, i.e. at most 5000 meters.
    Optional<S2BestEdgesQueryBase.Result<S1ChordAngle>> result =
        query.findClosestEdge(shapeIndexTarget);

    // If the result.isPresent(), it contains distance(), shape(), and edgeId().
    if (result.isPresent()) {
      // Get the actual closest polyline edge endpoints from the Result shape and edge id.
      S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
      S2ClosestEdgeQuery.getEdge(result.get(), edge);
      // Use the edge endpoints.
    }
  }
}

 * }</pre>
 * Query options are immutable and are set on the Builder, which is then used to construct queries
 * as shown in the example above. To find the closest edges to a query edge rather than a point,
 * use:
 *
 * <pre>{@code
 * S2ClosestEdgeQuery.EdgeTarget<S1ChordAngle> target = new S2ClosestEdgeQuery.EdgeTarget<>(v0, v1);
 * query.findClosestEdges(target);
 * }</pre>
 *
 * <p>Similarly you can find the closest edges to an S2Cell by using a
 * {@link S2ClosestEdgeQuery.CellTarget}, and you can find the closest edges to an arbitrary
 * collection of points, polylines, and polygons by using an
 * {@link S2ClosestEdgeQuery.ShapeIndexTarget}.
 *
 * <p>There are two overloads of the findClosestEdges() method:
 * <ol>
 *   <li> Users may simply get a returned list of Result objects. However, for cases where many
 *        calls to findClosestEdges will be made for a single query, the alternative can be much
 *        more efficient.
 *   <li> Users may instead provide a ResultVisitor that accepts (distance, shape, edge id) tuples.
 *        This avoids boxing values, but it does fully compute the set of results before any are
 *        visited.
 * </ol>
 *
 * The {@link Result} object contains the following accessors:
 *
 * <ul>
 *   <li>{@link Result#distance()} is the distance to the edge.
 *   <li>{@link Result#shape()} is the S2Shape containing the edge.
 *   <li>{@link Result#edgeId()} identifies the edge with the given shape.
 *   <li>{@link Result#isInterior()} indicates that the result is an interior point.
 * </ul>
 *
 * <p>The following convenience methods may also be useful:
 *
 * <ul>
 *   <li>{@link S2ClosestEdgeQuery#getEdge(Result, MutableEdge)} fills the MutableEdge with the
 *       endpoints of the Result edge.
 *   <li>{@link S2ClosestEdgeQuery#project(S2Point, Result)} computes the closest point on the
 *       Result edge to the given S2Point.
 * </ul>
 *
 * <p>You can find either the k closest edges, or all edges within a given radius, or both (i.e.,
 * the k closest edges up to a given maximum radius). E.g. to find all the edges within 5
 * kilometers, call
 *
 * <pre>{@code
 * options.setMaxDistance(S2Earth.toAngle(new Length.Kilometers(5)));
 * }</pre>
 *
 * <p>By default, *all* edges are returned, so you should always specify either maxResults() or
 * maxDistance() or both. There is also a findClosestEdge() convenience method that returns only the
 * single closest edge.
 *
 * <p>Note that by default, distances are measured to the boundary and interior of polygons. For
 * example, if a point is inside a polygon then its distance is zero. To change this behavior,
 * setIncludeInteriors(false) when building the query options.
 *
 * <p>If you only need to test whether the distance is below a given threshold (e.g., 10 km), you
 * can use the isDistanceLess() method. This is much faster than actually calculating the distance
 * with findClosestEdge(), since the implementation can stop as soon as it can prove that the
 * minimum distance is either above or below the threshold.
 *
 * <p>The implementation is designed to be fast for both simple and complex geometric objects.
 */
@CheckReturnValue
public abstract class S2ClosestEdgeQuery<D extends S1Distance<D>> extends S2BestEdgesQueryBase<D> {

  /** Target is the base interface of all the S2ClosestEdgeQuery-specific targets. */
  public interface Target<D extends S1Distance<D>> extends S2BestDistanceTarget<D> {}

  /** A target for finding the closest edges to a point. */
  public static class PointTarget<D extends S1Distance<D>>
      extends S2MinDistanceTargets.PointTarget<D> implements Target<D> {
    public PointTarget(S2Point p) {
      super(p);
    }

    /** See {@link S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }
  }

  /** A target for finding the closest edges to an edge. */
  public static class EdgeTarget<D extends S1Distance<D>> extends S2MinDistanceTargets.EdgeTarget<D>
      implements Target<D> {
    public EdgeTarget(S2Point a, S2Point b) {
      super(a, b);
    }

    /** See {@link S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }
  }

  /** A target for finding the closest edges to an S2Cell. */
  public static class CellTarget<D extends S1Distance<D>> extends S2MinDistanceTargets.CellTarget<D>
      implements Target<D> {
    public CellTarget(S2Cell c) {
      super(c);
    }

    /** See {@link S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 16;
    }
  }

  /** A target for finding the closest edges to an S2ShapeIndex. */
  public static class ShapeIndexTarget<D extends S1Distance<D>>
      extends S2MinDistanceTargets.ShapeIndexTarget<D>
      implements Target<D> {
    public ShapeIndexTarget(S2ShapeIndex index, S2BestEdgesQueryBase.Builder<D> queryBuilder) {
      super(index, queryBuilder);
    }

    /** See {@link S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }
  }

  /** The constructor is internal, as clients should use the Builder to construct queries. */
  S2ClosestEdgeQuery(Options<D> options) {
    super(options);
  }

  /**
   * Returns the closest edges to the given target that satisfy the current options. This method may
   * be called multiple times, with the same or different targets, and the same or different
   * options. However, if many calls will be made, the two alternative findClosestEdges() methods
   * below are more efficient.
   *
   * <p>Note that if {@code options().includeInteriors()} is true, the result vector may include
   * some entries with edgeId == -1. This indicates that the target is contained by or intersects
   * the indexed polygon with the given shape. Such results may be identifed by calling
   * {@code result.isInterior()}.
   */
  public List<Result<D>> findClosestEdges(S2BestDistanceTarget<D> target) {
    return findBestEdges(target);
  }

  /**
   * This version calls the provided ResultVisitor with the closest (shape, edge) pairs and their
   * distances to the target that satisfy the current options. This is the most efficient approach.
   *
   * <p>If no shape edges satisfy the query options, the visitor will not be called.
   */
  public void findClosestEdges(S2BestDistanceTarget<D> target, ResultVisitor<D> visitor) {
    findBestEdges(target, visitor);
  }

  /**
   * Returns the closest edge to the target, if one satisfies the search options. Otherwise, if none
   * satisfies the options, the returned Optional will not be present.
   *
   * <p>Note that if options.includeInteriors() is true, result.isInterior() should be called to
   * check whether the result represents an interior point (in which case edgeId() == -1).
   */
  public Optional<Result<D>> findClosestEdge(S2BestDistanceTarget<D> target) {
    return findBestEdge(target);
  }

  /**
   * Returns the minimum distance to the target. If the index or target is empty, returns a distance
   * greater than any valid distance.
   *
   * <p>Use isDistanceLess() if you only want to compare the distance against a threshold value,
   * since it is often much faster.
   */
  public D getDistance(S2BestDistanceTarget<D> target) {
    Optional<Result<D>> result = findBestEdge(target);
    return result.isPresent() ? result.get().distance() : beyondWorstDistance();
  }

  /**
   * Returns true if the distance to "target" is less than "limit".
   *
   * <p>This method is usually much faster than getDistance(), since it is less work to determine
   * whether the minimum distance is above or below a threshold than it is to calculate the actual
   * minimum distance.
   */
  public boolean isDistanceLess(S2BestDistanceTarget<D> target, D limit) {
    maxResults = 1;
    maxError = worstDistance();
    distanceLimit = limit;

    // Determine if the closest edge is within the limit.
    findBestEdgesInternal(target);
    boolean result = !resultQueue.isEmpty();
    resultQueue.clear();
    return result;
  }

  /** Gets the edge endpoints for the given result, which must not be interior. */
  public static void getEdge(Result<?> result, S2Shape.MutableEdge edge) {
    result.shape().getEdge(result.edgeId(), edge);
  }

  /** Returns the point on the given Result edge that is closest to the given targetPoint. */
  public S2Point project(S2Point targetPoint, Result<D> result) {
    if (result.edgeId() < 0) {
      return targetPoint;
    }
    S2Shape.MutableEdge resultEdge = new S2Shape.MutableEdge();
    getEdge(result, resultEdge);
    return S2EdgeUtil.getClosestPoint(targetPoint, resultEdge.getStart(), resultEdge.getEnd());
  }

  /**
   * Visits shapes in the index with interiors containing a point of a connected component of the
   * given target. Returns true if all such shapes were visited, or false if the {@link
   * Options#maxResults()} limit was reached. Note that the visited shapes may either intersect or
   * completely contain a connected component of the target.
   *
   * <p>This is a low-level method, visible for testing. Clients should use {@link
   * #findClosestEdges} with a maxDistance of S1ChordAngle.ZERO and check {@link
   * Result#isInterior()}.
   */
  @VisibleForTesting
  @CanIgnoreReturnValue
  public boolean visitContainingShapes(
      Target<D> target, S2ContainsPointQuery.ShapeVisitor visitor) {
    return visitBestDistanceContainingShapes(target, visitor);
  }

  @Override
  @CanIgnoreReturnValue
  protected boolean visitBestDistanceContainingShapes(
      S2BestDistanceTarget<D> target, S2ContainsPointQuery.ShapeVisitor visitor) {
    S2ContainsPointQuery containsPointQuery = new S2ContainsPointQuery(index);
    return target.visitConnectedComponentPoints(
        targetPoint -> containsPointQuery.visitContainingShapes(targetPoint, visitor));
  }

  /**
   * Subclasses {@code S2BestEdgesQueryBase.Builder<S1ChordAngle>} for finding closest edges using
   * S1ChordAngle as the distance type.
   *
   * <p>Provides the additional convenience methods setMaxDistance(S1Angle),
   * setInclusiveMaxDistance(S1ChordAngle), and setConservativeMaxDistance(S1ChordAngle).
   */
  public static class Builder extends S2BestEdgesQueryBase.Builder<S1ChordAngle> {
    /**
     * Constructs a new Builder with inclusive default values: unlimited results, no maximum
     * distance, maxError of zero, includeInteriors true, brute force false.
     */
    public Builder() {
      super(S1ChordAngle.INFINITY, S1ChordAngle.ZERO);
    }

    /** Constructs a new Builder with values copied from the given Options. */
    public Builder(S2BestEdgesQueryBase.Options<S1ChordAngle> options) {
      super(options);
    }

    /**
     * Builds a new Query, which is a {@code S2ClosestEdgeQuery<S1ChordAngle>} using the current
     * Options of this Builder.
     */
    @Override
    public Query build() {
      return new Query(new Options<>(this));
    }

    /**
     * Builds a new Query, which is a {@code S2ClosestEdgeQuery<S1ChordAngle>} for the provided
     * S2ShapeIndex using the current Options of this Builder.
     */
    @Override
    public Query build(S2ShapeIndex index) {
      return new Query(new Options<>(this), index);
    }

    /** Specifies the maximum number of results to be returned. */
    @CanIgnoreReturnValue
    @Override
    public Builder setMaxResults(int maxResults) {
      this.maxResults = maxResults;
      return this;
    }

    /**
     * Specifies that only edges whose distance to the target is less than maxDistance should be
     * returned.
     *
     * <p>Note that edges whose distance is exactly equal to "maxDistance" are not returned. In
     * most cases this doesn't matter, since distances are not computed exactly in the first
     * place, but if such edges are needed then you can use {@link
     * Builder#setInclusiveMaxDistance(S1ChordAngle)}.
     */
    @CanIgnoreReturnValue
    public Builder setMaxDistance(S1ChordAngle maxDistance) {
      this.distanceLimit = maxDistance;
      return this;
    }

    /**
     * Like {@link Builder#setMaxDistance(S1ChordAngle)}, but maxDistance is provided as an
     * S1Angle.
     */
    @CanIgnoreReturnValue
    public Builder setMaxDistance(S1Angle maxDistance) {
      this.distanceLimit = S1ChordAngle.fromS1Angle(maxDistance);
      return this;
    }

    /**
     * Specifies that only edges whose distance to the target is less than or equal to maxDistance
     * should be returned.
     */
    @CanIgnoreReturnValue
    public Builder setInclusiveMaxDistance(S1ChordAngle maxDistance) {
      this.distanceLimit = maxDistance.successor();
      return this;
    }

    /** Like {@link setInclusiveMaxDistance(S1ChordAngle)} but takes an S1Angle for convenience. */
    @CanIgnoreReturnValue
    public Builder setInclusiveMaxDistance(S1Angle maxDistance) {
      setInclusiveMaxDistance(S1ChordAngle.fromS1Angle(maxDistance));
      return this;
    }

    /**
     * Like setInclusiveMaxDistance(), except that "maxDistance" is also increased by the maximum
     * error in distance calculations. This ensures that all edges whose true distance is less than
     * or equal to "maxDistance" will be returned (along with some edges whose true distance is
     * slightly greater).
     *
     * <p>Algorithms that need to do exact distance comparisons can use this option to find a set of
     * candidate edges that can then be filtered further (e.g., using {@code
     * S2Predicates.compareDistance}).
     */
    @CanIgnoreReturnValue
    public Builder setConservativeMaxDistance(S1ChordAngle maxDistance) {
      this.distanceLimit =
          maxDistance.plusError(S2EdgeUtil.getMinDistanceMaxError(maxDistance)).successor();
      return this;
    }

    /**
     * Like {@link setConservativeMaxDistance(S1ChordAngle)} but takes an S1Angle for convenience.
     */
    @CanIgnoreReturnValue
    public Builder setConservativeMaxDistance(S1Angle maxDistance) {
      setConservativeMaxDistance(S1ChordAngle.fromS1Angle(maxDistance));
      return this;
    }

    /** All results will have distance less than maxDistance() from the target. */
    public S1ChordAngle maxDistance() {
      return distanceLimit;
    }

    /**
     * A non-zero maxError specifies that edges with distance up to maxError further than the true
     * closest edges may be substituted in the result set, as long as such edges satisfy all the
     * the remaining search criteria. This option only has an effect if {@link maxResults()} is
     * also specified; otherwise all edges that satisfy the maxDistance() will be returned.
     *
     * <p>Note that this does not affect how the distance between edges is computed; it simply
     * gives the algorithm permission to stop the search early, as soon as the best possible
     * improvement drops below {@link maxError()}.
     *
     * <p>This can be used to implement distance predicates efficiently. For example, to determine
     * whether the minimum distance is less than 'threshold', set {@code maxResults == 1} and {@code
     * maxDistance() == maxError() == threshold}. This causes the algorithm to terminate as soon as
     * it finds any edge whose distance is less than 'threshold', rather than continuing to search
     * for an edge that is even closer.
     */
    @CanIgnoreReturnValue
    @Override
    public Builder setMaxError(S1ChordAngle maxError) {
      this.maxError = maxError;
      return this;
    }

    /** Like {@link setMaxError(S1ChordAngle)}, but maxError is provided as an S1Angle. */
    @CanIgnoreReturnValue
    public Builder setMaxError(S1Angle maxError) {
      this.maxError = S1ChordAngle.fromS1Angle(maxError);
      return this;
    }

    /**
     * True if polygon interiors in the queried index are considered when computing distance.
     *
     * <p>When true, polygons that contain some geometry have a distance of zero to it. For targets
     * consisting of multiple connected components, this occurs if any component is contained.
     * This is indicated in the results by returning a (shape, edgeId) pair with {@code edgeId ==
     * -1}, i.e. this value denotes the polygons's interior.
     *
     * <p>Note that this does not control if the interiors of _target_ polygons should be
     * included; that is a separate option on targets.
     *
     * <p>Note that for efficiency, any polygon that intersects the target may or may not have an
     * (edgeId == -1) result. Such results are optional because in that case the distance from the
     * polygon is already zero.
     */
    @CanIgnoreReturnValue
    @Override
    public Builder setIncludeInteriors(boolean includeInteriors) {
      this.includeInteriors = includeInteriors;
      return this;
    }

    /**
     * If true, distances should be computed by examining every edge rather than using the
     * S2ShapeIndex. This is useful for testing, benchmarking, and debugging.
     */
    @CanIgnoreReturnValue
    @Override
    public Builder setUseBruteForce(boolean useBruteForce) {
      this.useBruteForce = useBruteForce;
      return this;
    }
  }

  /**
   * Subclasses {@code S2ClosestEdgeQuery<S1ChordAngle>} for finding closest edges using
   * S1ChordAngle as the distance type.
   *
   * <p>Provides the additional convenience methods isDistanceLessOrEqual() and
   * isConservativeDistanceLessOrEqual().
   */
  public static class Query extends S2ClosestEdgeQuery<S1ChordAngle> {
    /** Constructor for internal use. Clients should use the Builder to create queries. */
    Query(Options<S1ChordAngle> options) {
      super(options);
    }

    /** Constructor for internal use. Clients should use the Builder to create queries. */
    Query(Options<S1ChordAngle> options, S2ShapeIndex index) {
      super(options);
      init(index);
    }

    /** Convenience method to get a Builder from this Query's current options. */
    public Builder toBuilder() {
      return new Builder(options());
    }

    @Override
    protected DistanceCollector<S1ChordAngle> newDistanceCollector() {
      return S1ChordAngle.minCollector();
    }

    @Override
    protected boolean atBestLimit(DistanceCollector<S1ChordAngle> distanceCollector) {
      return distanceCollector.distance().getLength2() <= 0;
    }

    @Override
    protected Comparator<S1ChordAngle> distanceComparator() {
      return naturalOrder();
    }

    @Override
    protected S1ChordAngle zeroDistance() {
      return S1ChordAngle.ZERO;
    }

    @Override
    protected S1ChordAngle bestDistance() {
      return S1ChordAngle.ZERO;
    }

    @Override
    protected S1ChordAngle worstDistance() {
      return S1ChordAngle.STRAIGHT;
    }

    @Override
    protected S1ChordAngle beyondWorstDistance() {
      return S1ChordAngle.INFINITY;
    }

    /** Adjust the given 'value' towards zero by the maximum error allowed by options. */
    @Override
    protected S1ChordAngle errorBoundedDistance(S1ChordAngle value) {
      return S1ChordAngle.sub(value, options().maxError());
    }

    /**
     * For closest edges, the search cap radius is the sum of the target cap radius bounding the
     * area with distance zero from the target, the maximum distance from the target, and the
     * potential error in computing the max distance from an angle.
     */
    @Override
    protected S1ChordAngle searchCapRadius(S1ChordAngle targetCapRadius, S1ChordAngle maxDistance) {
      checkArgument(!targetCapRadius.isNegative());
      checkArgument(!targetCapRadius.isInfinity());
      checkArgument(!maxDistance.isNegative());
      checkArgument(!maxDistance.isInfinity());
      return S1ChordAngle.add(
          targetCapRadius, maxDistance.plusError(maxDistance.getS1AngleConstructorMaxError()));
    }

    /**
     * Like {@link #isDistanceLess(S2BestDistanceTarget, S1Distance)}, but also returns true if
     * the distance to "target" is exactly equal to "limit".
     */
    public boolean isDistanceLessOrEqual(Target<S1ChordAngle> target, S1ChordAngle limit) {
      // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
      // not the same-named Options fields.
      distanceLimit = limit.successor();
      maxResults = 1;
      maxError = worstDistance();

      findBestEdgesInternal(target);
      boolean result = !resultQueue.isEmpty();
      resultQueue.clear();

      return result;
    }

    /**
     * Like {@link isDistanceLessOrEqual(Target, S1ChordAngle)}, except that "limit" is increased
     * by the maximum error in the distance calculation. This ensures that this function returns
     * true whenever the true, exact distance is less than or equal to "limit".
     *
     * <p>For example, suppose that we want to test whether two geometries might intersect each
     * other after they are snapped together using S2Builder (using the IdentitySnapFunction with a
     * given "snapRadius"). Since S2Builder uses exact distance predicates, we need to measure the
     * distance between the two geometries conservatively. If the distance is definitely greater
     * than "snapRadius", then the geometries are guaranteed to not intersect after snapping.
     */
    public boolean isConservativeDistanceLessOrEqual(
        Target<S1ChordAngle> target, S1ChordAngle limit) {
      // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
      // not the same-named Options fields.
      distanceLimit = limit.plusError(S2EdgeUtil.getMinDistanceMaxError(limit)).successor();
      maxResults = 1;
      maxError = worstDistance();

      findBestEdgesInternal(target);
      boolean result = !resultQueue.isEmpty();
      resultQueue.clear();

      return result;
    }
  }

  /**
   * Convenience method to create a new Builder, which extends {@link S2ClosestEdgeQuery.Builder}
   * using S1ChordAngle as the distance type.
   *
   * <p>The initial options are the most permissive default values. There is no limit on the maximum
   * distance or number of results returned. The maxError option is zero. The interiors of polygons
   * in the index are included, i.e. a target contained by an index polygon is at distance zero.
   */
  public static Builder builder() {
    return new Builder();
  }

  /**
   * Convenience method to create a ShapeIndexTarget using S1ChordAngle as the distance type.
   *
   * <p>This is used for finding closest edges from a query ShapeIndex to a second ShapeIndex
   * wrapped by the ShapeIndexTarget.
   */
  public static ShapeIndexTarget<S1ChordAngle> createShapeIndexTarget(S2ShapeIndex index) {
    return new ShapeIndexTarget<>(index, new Builder());
  }
}
