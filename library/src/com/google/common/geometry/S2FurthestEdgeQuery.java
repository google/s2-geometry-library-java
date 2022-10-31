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

import com.google.common.annotations.VisibleForTesting;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

/**
 * {@link S2FurthestEdgeQuery} is a helper class for searching within an {@link S2ShapeIndex} to
 * find the furthest edge(s) to a given target point, edge, S2Cell, or {@link ShapeIndex geometry
 * collection}. Furthest edges are those with maximum distance from any point on that edge to any
 * point on the target geometry.
 *
 * <p>For example, given a set of polylines, the following code efficiently finds the furthest 100
 * polyline edges from a query target point that are at least 10 kilometers from it, and then finds
 * the single furthest polyline edge from a second target ShapeIndex.
 *
 * <p>This example code demonstrates the use of the provided convenience methods for working with
 * S2FurthestEdgeQuery using an S1Distance abstract type that is specifically S1ChordAngle. This is
 * recommended for most clients.
 *
 * <pre>{@code

class FurthestEdgeDemo {
  // Recommended: S2Earth.toAngle(new Length.Kilometers(10))
  static final S1Angle MIN_DISTANCE = S1Angle.fromEarthDistance(10000);

  public static void example(
      List<S2Polyline> polylines, S2ShapeIndex targetIndex, S2Point targetPoint) {
    // Add the given polylines to a shape index.
    S2ShapeIndex polylinesIndex = new S2ShapeIndex();
    for (S2Polyline polyline : polylines) {
      polylinesIndex.add(S2LaxPolylineShape.create(polyline));
    }

    // Create a Builder, modifying the default values with a minimum distance and maximum number
    // of results.
    S2FurthestEdgeQuery.Builder builder =
        S2FurthestEdgeQuery.builder().setInclusiveMinDistance(MIN_DISTANCE).setMaxResults(100);

    // Construct the query. Note that "Query" extends S2FurthestEdgeQuery<S1ChordAngle>.
    S2FurthestEdgeQuery.Query query = builder.build(polylinesIndex);

    // Create a target that is a single point.
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> pointTarget =
        new S2FurthestEdgeQuery.PointTarget<>(targetPoint);

    // Find and visit up to 100 polyline edges in polylinesIndex that are the farthest from the
    // target point, and at least 10 km from it.
    query.findFurthestEdges(
        pointTarget,
        (S1ChordAngle distance, S2Shape shape, int edgeId) -> {
          // Do something with each of these furthest edge results.
          return true;
        });

    // Create a target that is another shapeIndex.
    S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> shapeIndexTarget =
        S2FurthestEdgeQuery.createShapeIndexTarget(targetIndex);
    shapeIndexTarget.setIncludeInteriors(true);

    // Find the single farthest polyline edge to any edge or polygon interior in the
    // ShapeIndexTarget, if any matches the query options, i.e. at least 10000 meters.
    Optional<S2BestEdgesQueryBase.Result<S1ChordAngle>> result =
        query.findFurthestEdge(shapeIndexTarget);

    // If the result.isPresent() it contains distance(), shape(), and edgeId().
    if (result.isPresent()) {
      // Get the actual furthest polyline edge endpoints from the Result shape and edge id.
      S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
      S2ClosestEdgeQuery.getEdge(result.get(), edge);
      // Do something with the edge endpoints.
    }
  }
}

 * }</pre>
 *
 * Query options are immutable and are set on the Builder, which is then used to construct queries
 * as shown in the example above. To find the furthest edges from a query edge rather than a point,
 * use:
 *
 * <pre>{@code
 * S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target
 *     = new S2FurthestEdgeQuery.EdgeTarget<>(v0, v1);
 * query.findFurthestEdges(target);
 * }</pre>
 *
 * <p>Similarly you can find the furthest edges from an S2Cell by using a
 * {@link S2FurthestEdgeQuery.CellTarget}, and you can find the furthest edges from an arbitrary
 * collection of points, polylines, and polygons by using an
 * {@link S2FurthestEdgeQuery.ShapeIndexTarget}.
 *
 * <p>There are two overloads of the findFurthestEdges() method:
 * <ol>
 *   <li> Users may simply get a returned list of Result objects. However, for cases where many
 *        calls to findFurthestEdges will be made for a single query, the alternative can be much
 *        more efficient.
 *   <li> Users may instead provide a ResultVisitor that accepts (distance, shape, edge id) tuples
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
 * <p>The following convenience method may also be useful:
 * {@link S2FurthestEdgeQuery#getEdge(Result, MutableEdge)} fills the MutableEdge with the endpoints
 * of the Result edge.
 *
 * <p>You can find either the k furthest edges, or all edges no closer than a given radius, or both
 * (i.e., the k furthest edges no closer than a given minimum radius). E.g. to find all the edges
 * further than 5 kilometers, call:
 *
 * <pre>{@code
 * options.setMinDistance(S2Earth.toAngle(new Length.Kilometers(5)));
 * }</pre
 *
 * <p>
 *
 * By default, *all* edges are returned, so you should always specify either maxResults() or
 * minDistance() or both. Setting min distance may not be very restrictive, so strongly consider
 * using maxResults().There is also a findFurthestEdge() convenience method that returns only the
 * single furthest edge.
 *
 * <p>Note that by default, distances are measured to the boundary and interior of polygons. For
 * example, if a point is inside a polygon then its distance is zero. To change this behavior,
 * setIncludeInteriors(false) when building the query options.
 *
 * <p>If you only need to test whether the distance is above a given threshold (e.g., 10 km), you
 * can use the isDistanceGreater() method. This is much faster than actually calculating the
 * distance with findFurthestEdge(), since the implementation can stop as soon as it can prove that
 * the maximum distance is either above or below the threshold.
 *
 * <p>The implementation is designed to be fast for both simple and complex geometric objects.
 */
@CheckReturnValue
public abstract class S2FurthestEdgeQuery<D extends S1Distance<D>> extends S2BestEdgesQueryBase<D> {

  /** Target is the base interface of all the S2FurthestEdgeQuery-specific targets. */
  public interface Target<D extends S1Distance<D>> extends S2BestDistanceTarget<D> {}

  /** A target for finding the furthest edges from a point. */
  public static class PointTarget<D extends S1Distance<D>>
      extends S2MaxDistanceTargets.PointTarget<D> implements Target<D> {
    public PointTarget(S2Point p) {
      super(p);
    }

    /** See {@link S2FurthestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 64;
    }
  }

  /** A target for finding the furthest edges from an edge. */
  public static class EdgeTarget<D extends S1Distance<D>> extends S2MaxDistanceTargets.EdgeTarget<D>
      implements Target<D> {
    public EdgeTarget(S2Point a, S2Point b) {
      super(a, b);
    }

    /** See {@link S2FurthestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }
  }

  /** A target for finding the furthest edges from an S2Cell. */
  public static class CellTarget<D extends S1Distance<D>> extends S2MaxDistanceTargets.CellTarget<D>
      implements Target<D> {
    public CellTarget(S2Cell c) {
      super(c);
    }

    /** See {@link S2FurthestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 96;
    }
  }

  /** A target for finding the furthest edges from an S2Cell. */
  public static class ShapeIndexTarget<D extends S1Distance<D>>
      extends S2MaxDistanceTargets.ShapeIndexTarget<D>
      implements Target<D> {
    public ShapeIndexTarget(S2ShapeIndex index, S2BestEdgesQueryBase.Builder<D> queryBuilder) {
      super(index, queryBuilder);
    }

    /** See {@link S2FurthestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 32;
    }
  }

  /** The constructor is internal, as clients should use the Builder to construct queries. */
  S2FurthestEdgeQuery(Options<D> options) {
    super(options);
  }

  /**
   * Returns the furthest edges to the given target that satisfy the current options. This method
   * may be called multiple times, with the same or different targets, and the same or different
   * options. However, if many calls will be made, the two alternative findFurthestEdges() methods
   * below are more efficient.
   *
   * <p>Note that if options().includeInteriors() is true, the result list may include some entries
   * with edgeId == -1. This indicates that the furthest distance is attained at a point in the
   * interior of the indexed polygon with the given shape. Such results may be identifed by calling
   * result.isInterior().
   */
  public List<Result<D>> findFurthestEdges(S2BestDistanceTarget<D> target) {
    return findBestEdges(target);
  }

  /**
   * This version calls the provided ResultVisitor with the furthest (shape, edge) pairs and their
   * distances to the target that satisfy the current options. This is the most efficient approach.
   *
   * <p>If no shape edges satisfy the query options, the visitor will not be called.
   */
  public void findFurthestEdges(S2BestDistanceTarget<D> target, ResultVisitor<D> visitor) {
    findBestEdges(target, visitor);
  }

  /**
   * Returns the furthest edge from the target. If no edge satisfies the search criteria, then the
   * Optional will not be present.
   *
   * <p>Note that if options.includeInteriors() is true, result.isInterior() should be called to
   * check whether the result represents an interior point (in which case edgeId() == -1).
   */
  public Optional<Result<D>> findFurthestEdge(S2BestDistanceTarget<D> target) {
    return findBestEdge(target);
  }

  /**
   * Returns the maximum distance to the target. If the index or target is empty, returns a distance
   * less than any valid distance.
   *
   * <p>Use isDistanceGreater() if you only want to compare the distance against a threshold value,
   * since it is often much faster.
   */
  public D getDistance(S2BestDistanceTarget<D> target) {
    Optional<Result<D>> result = findBestEdge(target);
    return result.isPresent() ? result.get().distance() : beyondWorstDistance();
  }

  /**
   * Returns true if the distance to "target" is greater than "limit".
   *
   * <p>This method is usually much faster than getDistance(), since it is less work to determine
   * whether the maximum distance is above or below a threshold than it is to calculate the actual
   * maximum distance.
   */
  public boolean isDistanceGreater(S2BestDistanceTarget<D> target, D limit) {
    maxResults = 1;
    maxError = worstDistance();
    distanceLimit = limit;

    // Determine if the furthest edge is beyond the limit.
    findBestEdgesInternal(target);
    boolean result = !resultQueue.isEmpty();
    resultQueue.clear();
    return result;
  }

  /** Gets the edge endpoints for the given result, which must not be interior. */
  public static void getEdge(Result<?> result, S2Shape.MutableEdge edge) {
    result.shape().getEdge(result.edgeId(), edge);
  }

  /**
   * Visits shapes in the index with interiors antipodal to a point of a connected component of the
   * given target. Returns true if all such shapes were visited, or false if the {@link
   * Options#maxShapes()} limit was reached. Note that the visited shapes may either intersect or
   * completely contain a connected component of the target.
   *
   * <p>This is a low-level method, visible for testing. Clients should use {@link
   * findFurthestEdges} with a minDistance of S1ChordAngle.STRAIGHT and then check {@link
   * Result#isInterior()}.
   */
  @VisibleForTesting
  @CanIgnoreReturnValue
  public boolean visitAntipodalShapes(Target<D> target, S2ContainsPointQuery.ShapeVisitor visitor) {
    return visitBestDistanceContainingShapes(target, visitor);
  }

  @Override
  @CanIgnoreReturnValue
  protected boolean visitBestDistanceContainingShapes(
      S2BestDistanceTarget<D> target, S2ContainsPointQuery.ShapeVisitor visitor) {
    S2ContainsPointQuery containsPointQuery = new S2ContainsPointQuery(index);
    return target.visitConnectedComponentPoints(
        targetPoint -> containsPointQuery.visitContainingShapes(targetPoint.neg(), visitor));
  }

  /**
   * Subclasses {@code S2BestEdgesQueryBase.Builder<S1ChordAngle>} for finding furthest edges using
   * S1ChordAngle as the distance type.
   *
   * <p>Provides the additional convenience methods setMinDistance(S1Angle),
   * setInclusiveMinDistance(S1ChordAngle), and setConservativeMinDistance(S1ChordAngle).
   */
  public static class Builder extends S2BestEdgesQueryBase.Builder<S1ChordAngle> {
    /**
     * Constructs a new Builder with inclusive default values: unlimited results, no minimum
     * distance, maxError of zero, includeInteriors true, brute force false.
     */
    public Builder() {
      super(S1ChordAngle.ZERO, S1ChordAngle.ZERO);
    }

    /** Constructs a new Builder with values copied from the given Options. */
    public Builder(Options<S1ChordAngle> options) {
      super(options);
    }

    /**
     * Builds a new Query, which is a {@code S2FurthestEdgeQuery<S1ChordAngle>} using the current
     * Options of this Builder.
     */
    @Override
    public Query build() {
      return new Query(new Options<>(this));
    }

    /**
     * Builds a new Query, which is a {@code S2FurthestEdgeQuery<S1ChordAngle>} for the provided
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
     * Specifies that only edges whose distance to the target is greater than minDistance should
     * be returned.
     *
     * <p>Note that edges whose distance is exactly equal to "minDistance" are not returned. In
     * most cases this doesn't matter, since distances are not computed exactly in the first
     * place, but if such edges are needed then you can use {@link
     * Builder#setInclusiveMinDistance(S1ChordAngle)}.
     */
    @CanIgnoreReturnValue
    public Builder setMinDistance(S1ChordAngle minDistance) {
      this.distanceLimit = minDistance;
      return this;
    }

    /**
     * Like {@link Options.Builder#setMinDistance(S1ChordAngle)}, but minDistance is provided as an
     * S1Angle.
     */
    @CanIgnoreReturnValue
    public Builder setMinDistance(S1Angle minDistance) {
      this.distanceLimit = S1ChordAngle.fromS1Angle(minDistance);
      return this;
    }

    /**
     * Specifies that only edges whose distance to the target is greater than or equal to
     * minDistance should be returned.
     */
    @CanIgnoreReturnValue
    public Builder setInclusiveMinDistance(S1ChordAngle minDistance) {
      this.distanceLimit = minDistance.predecessor();
      return this;
    }

    /** Like {@link setInclusiveMinDistance(S1ChordAngle)} but takes an S1Angle for convenience. */
    @CanIgnoreReturnValue
    public Builder setInclusiveMinDistance(S1Angle minDistance) {
      setInclusiveMinDistance(S1ChordAngle.fromS1Angle(minDistance));
      return this;
    }

    /**
     * Like setInclusiveMinDistance(), except that "minDistance" is also decreased by the maximum
     * error in the distance calculation. This ensures that all edges whose true distance is more
     * than or equal to "minDistance" will be returned (along with some edges whose true distance is
     * slightly less).
     *
     * <p>Algorithms that need to do exact distance comparisons can use this option to find a set of
     * candidate edges that can then be filtered further (e.g., using {@code
     * S2Predicates.compareDistance}).
     */
    @CanIgnoreReturnValue
    public Builder setConservativeMinDistance(S1ChordAngle minDistance) {
      this.distanceLimit =
          minDistance.plusError(-S2EdgeUtil.getMinDistanceMaxError(minDistance)).predecessor();
      return this;
    }

    /**
     * Like {@link setConservativeMinDistance(S1ChordAngle)} but takes an S1Angle for convenience.
     */
    @CanIgnoreReturnValue
    public Builder setConservativeMinDistance(S1Angle minDistance) {
      setConservativeMinDistance(S1ChordAngle.fromS1Angle(minDistance));
      return this;
    }

    /** All results will have distance more than minDistance() from the target. */
    public S1ChordAngle minDistance() {
      return distanceLimit;
    }

    /**
     * A non-zero maxError specifies that edges with distance up to maxError closer than the true
     * furthest edges may be substituted in the result set, as long as such edges satisfy all the
     * the remaining search criteria. This option only has an effect if {@link maxResults()} is
     * also specified; otherwise all edges that satisfy the minDistance() will be returned.
     *
     * <p>Note that this does not affect how the distance between edges is computed; it simply
     * gives the algorithm permission to stop the search early, as soon as the best possible
     * improvement drops below {@link maxError()}.
     *
     * <p>This can be used to implement distance predicates efficiently. For example, to determine
     * whether the maximum distance is more than 'threshold', set {@code maxResults == 1} and {@code
     * minDistance() == maxError() == threshold}. This causes the algorithm to terminate as soon as
     * it finds any edge whose distance is greater than 'threshold', rather than continuing to
     * search for an edge that is even further.
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
     * <p>When true, polygons that contain some geometry have a distance of zero to it, and polygons
     * with interiors that are antipodal to some geometry have the maximum possible distance to it.
     * For targets consisting of multiple connected components, this occurs if any component is
     * contained. This is indicated in the results by returning a (shape, edgeId) pair with
     * {@code edgeId == -1}, i.e. this value denotes the polygons's interior.
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
   * Subclasses {@code S2FurthestEdgeQuery<S1ChordAngle>} for finding furthest edges using
   * S1ChordAngle as the distance type.
   *
   * <p>Provides the additional convenience methods isDistanceGreaterOrEqual() and
   * isConservativeDistanceGreaterOrEqual().
   */
  public static class Query extends S2FurthestEdgeQuery<S1ChordAngle> {
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
      return S1ChordAngle.maxCollector();
    }

    @Override
    protected boolean atBestLimit(DistanceCollector<S1ChordAngle> distanceCollector) {
      return distanceCollector.distance().getLength2() >= S1ChordAngle.MAX_LENGTH2;
    }

    @Override
    protected Comparator<S1ChordAngle> distanceComparator() {
      return Comparator.reverseOrder();
    }

    @Override
    protected S1ChordAngle zeroDistance() {
      return S1ChordAngle.ZERO;
    }

    @Override
    protected S1ChordAngle bestDistance() {
      return S1ChordAngle.STRAIGHT;
    }

    @Override
    protected S1ChordAngle worstDistance() {
      return S1ChordAngle.ZERO;
    }

    @Override
    protected S1ChordAngle beyondWorstDistance() {
      return S1ChordAngle.NEGATIVE;
    }

    /** Adjust the given 'value' towards the maximum by the maximum error allowed by options. */
    @Override
    protected S1ChordAngle errorBoundedDistance(S1ChordAngle value) {
      return S1ChordAngle.add(value, options().maxError());
    }

    /**
     * For furthest edges, the search cap radius is the sum of radius of the cap bounding points
     * that are antipodal (at maximum distance) relative to the target, and 180 degrees minus the
     * minimum distance. As the minimum distance from the target increases, the search cap gets
     * smaller, until at the limit, when the minimum distance from the target is 180 degrees, only
     * the antipodal cap must be searched. At the other limit, with no minimum distance, the whole
     * sphere must be searched for results, so the search cap radius is STRAIGHT.
     */
    @Override
    protected S1ChordAngle searchCapRadius(
        S1ChordAngle antipodalCapRadius, S1ChordAngle minDistance) {
      checkArgument(!antipodalCapRadius.isNegative());
      checkArgument(!antipodalCapRadius.isInfinity());
      checkArgument(!minDistance.isNegative());
      checkArgument(!minDistance.isInfinity());
      return S1ChordAngle.add(
          antipodalCapRadius, S1ChordAngle.sub(S1ChordAngle.STRAIGHT, minDistance));
    }

    /**
     * Like {@link S2FurthestEdgeQuery#isDistanceGreater(Target, D)}, but also returns true if the
     * distance to "target" is exactly equal to "limit".
     */
    public boolean isDistanceGreaterOrEqual(Target<S1ChordAngle> target, S1ChordAngle limit) {
      // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
      // not the same-named Options fields.
      distanceLimit = limit.predecessor();
      maxResults = 1;
      maxError = worstDistance();

      findBestEdgesInternal(target);
      boolean result = !resultQueue.isEmpty();
      resultQueue.clear();

      return result;
    }

    /**
     * Like {@link isDistanceGreaterOrEqual(Target, S1ChordAngle)}, except that "limit" is
     * decreased by the maximum error in the distance calculation. This ensures that this function
     * returns true whenever the true, exact distance is greater than or equal to "limit".
     */
    public boolean isConservativeDistanceGreaterOrEqual(
        Target<S1ChordAngle> target, S1ChordAngle limit) {
      // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
      // not the same-named Options fields.
      distanceLimit = limit.plusError(-S2EdgeUtil.getMinDistanceMaxError(limit)).predecessor();
      maxResults = 1;
      maxError = worstDistance();

      findBestEdgesInternal(target);
      boolean result = !resultQueue.isEmpty();
      resultQueue.clear();

      return result;
    }
  }

  /**
   * Convenience method to create a new Builder, which extends {@link S2FurthestEdgeQuery.Builder}
   * using S1ChordAngle as the distance type.
   *
   * <p>The initial options are the most permissive default values. There is no limit on the minimum
   * distance or number of results returned. The maxError option is zero. The interiors of polygons
   * in the index are included, i.e. a target antipodal to the interior of an index polygon is at
   * maximum distance.
   */
  public static Builder builder() {
    return new Builder();
  }

  /**
   * Convenience method to create a ShapeIndexTarget using S1ChordAngle as the distance type.
   *
   * <p>This is used for finding farthest edges from a query ShapeIndex to a second ShapeIndex
   * wrapped by the ShapeIndexTarget.
   */
  public static ShapeIndexTarget<S1ChordAngle> createShapeIndexTarget(S2ShapeIndex index) {
    return new ShapeIndexTarget<>(index, new Builder());
  }
}
