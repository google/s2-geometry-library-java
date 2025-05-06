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
import org.jspecify.annotations.Nullable;

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
        (S1ChordAngle distance, int shapeId, int edgeId) -> {
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

    // If the result.isPresent(), it contains distance(), shapeId(), and edgeId().
    if (result.isPresent() && !result.get().isInterior()) {
      // Get the actual closest polyline edge endpoints from the Result shape and edge id.
      S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
      query.getEdge(result.get(), edge);
      // Use the edge endpoints.
    }
  }
}

 * }</pre>
 * Query options are immutable and are set on the Builder, which is then used to construct queries
 * as shown in the example above. To find the closest edges to a query edge rather than a point,
 * use:
 *
 * {@snippet :
 * S2ClosestEdgeQuery.EdgeTarget<S1ChordAngle> target = new S2ClosestEdgeQuery.EdgeTarget<>(v0, v1);
 * query.findClosestEdges(target);
 * }
 *
 * <p>Similarly you can find the closest edges to an S2Cell by using a {@link
 * S2ClosestEdgeQuery.CellTarget}, and you can find the closest edges to an arbitrary collection of
 * points, polylines, and polygons by using an {@link S2ClosestEdgeQuery.ShapeIndexTarget}.
 *
 * <p>There are two overloads of the findClosestEdges() method:
 *
 * <ol>
 *   <li>Users may simply get a returned list of Result objects. However, for cases where many calls
 *       to findClosestEdges will be made for a single query, the alternative can be much more
 *       efficient.
 *   <li>Users may instead provide a ResultVisitor that accepts (distance, shape, edge id) tuples.
 *       This avoids boxing values, but it does fully compute the set of results before any are
 *       visited.
 * </ol>
 *
 * The {@link Result} object contains the following accessors:
 *
 * <ul>
 *   <li>{@link Result#distance()} is the distance to the edge.
 *   <li>{@link Result#shapeId()} is the index of the S2Shape containing the edge in the index.
 *   <li>{@link Result#edgeId()} identifies the edge with the given shape.
 *   <li>{@link Result#isInterior()} indicates that the result is an interior point.
 * </ul>
 *
 * <p>{@link S2ClosestEdgeQuery#project(S2Point, Result)} computes the closest point on the Result
 * edge to the given S2Point.
 *
 * <p>You can find either the k closest edges, or all edges within a given radius, or both (i.e.,
 * the k closest edges up to a given maximum radius). E.g. to find all the edges within 5
 * kilometers, call
 *
 * {@snippet :
 * options.setMaxDistance(S2Earth.toAngle(new Length.Kilometers(5)));
 * }
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
  public static class PointTarget<D2 extends S1Distance<D2>>
      extends S2BestEdgesQueryBase.PointTarget<D2> implements Target<D2> {
    public PointTarget(S2Point p) {
      super(p);
    }

    /** See {@code S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }

    @Override
    public S2Cap getCapBound() {
      return S2Cap.fromAxisChord(point, S1ChordAngle.ZERO);
    }
  }

  /** A target for finding the closest edges to an edge. */
  public static class EdgeTarget<D2 extends S1Distance<D2>>
      extends S2BestEdgesQueryBase.EdgeTarget<D2> implements Target<D2> {
    public EdgeTarget(S2Point a, S2Point b) {
      super(a, b);
    }

    /** See {@code S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }

    @Override
    public S2Cap getCapBound() {
      double r2 = getHalfEdgeLength2();
      return S2Cap.fromAxisChord(a.add(b).normalize(), S1ChordAngle.fromLength2(r2));
    }
  }

  /** A target for finding the closest edges to an S2Cell. */
  public static class CellTarget<D2 extends S1Distance<D2>>
      extends S2BestEdgesQueryBase.CellTarget<D2> implements Target<D2> {
    public static final int MAX_BRUTE_FORCE_INDEX_SIZE = 16;

    public CellTarget(S2Cell c) {
      super(c);
    }

    /** See {@code S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return MAX_BRUTE_FORCE_INDEX_SIZE;
    }

    @Override
    public S2Cap getCapBound() {
      return cell.getCapBound();
    }
  }

  /**
   * A target for finding the closest edges to an S2ShapeIndex. Shape interiors are included by
   * default.
   */
  public static class ShapeIndexTarget<D extends S1Distance<D>>
      extends S2BestEdgesQueryBase.ShapeIndexTarget<D> implements Target<D> {

    /**
     * Clients using S1ChordAngle as their S1Distance type may find it convenient to use
     * {@link S2ClosestEdgeQuery#createShapeIndexTarget(S2ShapeIndex)}.
     *
     * <p>Otherwise, constructing a ShapeIndexTarget for a specific S1Distance type requires
     * providing a S2BestEdgesQueryBase.Builder for closest edges and the templated S1Distance
     * type.
     */
    public ShapeIndexTarget(S2ShapeIndex index, S2BestEdgesQueryBase.Builder<D> queryBuilder) {
      super(index, queryBuilder);
    }

    /** See {@code S2ClosestEdgeQueryBenchmark.DetermineGoodBruteForceIndexSizes}. */
    @Override
    public int maxBruteForceIndexSize() {
      return 40;
    }

    @Override
    @CanIgnoreReturnValue
    public boolean updateBestDistance(S2Point p, DistanceCollector<D> collector) {
      return updateBestDistance(new PointTarget<D>(p), collector);
    }

    @Override
    @CanIgnoreReturnValue
    public boolean updateBestDistance(S2Point v0, S2Point v1, DistanceCollector<D> collector) {
      return updateBestDistance(new EdgeTarget<D>(v0, v1), collector);
    }

    @Override
    @CanIgnoreReturnValue
    public boolean updateBestDistance(S2Cell cell, DistanceCollector<D> collector) {
      return updateBestDistance(new CellTarget<D>(cell), collector);
    }

    @Override
    public S2Cap getCapBound() {
      // TODO(torrey): This is only called once, from S2BestEdgesQueryBase.initQueue. It might be
      // better to rework that algorithm to use a covering of the target directly, rather than
      // getting a cap, enlarging it, and covering that, particularly here where the cap is over
      // a covering.

      // If the index hasn't been built and the number of edges is sufficiently small, we avoid
      // building the index just to compute the cap bound here. If the index has more than
      // CellTarget.MAX_BRUTE_FORCE_INDEX_SIZE edges, it will be built later anyway.
      int maxSize = CellTarget.MAX_BRUTE_FORCE_INDEX_SIZE;
      if (!index.isFresh() && S2ShapeUtil.countEdgesUpTo(index, maxSize) < maxSize) {
        S2Cap.Builder builder = new S2Cap.Builder();
        S2Shape.MutableEdge e = new S2Shape.MutableEdge();
        for (S2Shape shape : index.getShapes()) {
          for (int i = 0; i < shape.numEdges(); i++) {
            shape.getEdge(i, e);
            builder.add(e);
          }
        }
        return builder.build();
      }

      return new S2ShapeIndexRegion(index).getCapBound();
    }
  }

  /** The constructor is internal, as clients should use the Builder to construct queries. */
  S2ClosestEdgeQuery(Options<D> options) {
    super(options);
  }

  /**
   * Returns the closest (shape, edge, distance) Results to the given target that satisfy the
   * current options. The returned list is sorted by increasing edge distance from the target. This
   * method may be called multiple times, with the same or different targets and options. Note that
   * the findClosestEdges() methods below that take a visitor are more efficient, particularly if
   * many calls are made.
   *
   * <p>See {@link #findClosestEdges(S2BestDistanceTarget, ShapeFilter, ResultVisitor)} for more
   * details.
   */
  public List<Result<D>> findClosestEdges(S2BestDistanceTarget<D> target) {
    return findBestEdges(target);
  }

  /**
   * Returns the closest (shape, edge, distance) Results to the given target that satisfy the
   * current options, and with shapes that pass the given ShapeFilter. The returned list is sorted
   * by increasing distance from the target. This method may be called multiple times, with the same
   * or different targets, shapeFilter, and options. Note that the findClosestEdges() methods below
   * that take a visitor are more efficient, particularly if many calls are made.
   *
   * <p>See {@link #findClosestEdges(S2BestDistanceTarget, ShapeFilter, ResultVisitor)} for more
   * details.
   */
  public List<Result<D>> findClosestEdges(
      S2BestDistanceTarget<D> target, @Nullable ShapeFilter shapeFilter) {
    return findBestEdges(target, shapeFilter);
  }

  /**
   * Visits the closest (shape, edge, distance) Results to the given target that satisfy the current
   * options. This method may be called multiple times, with the same or different targets
   * and options. Results are visited in order of increasing distance from the target.
   *
   * <p>See {@link #findClosestEdges(S2BestDistanceTarget, ShapeFilter, ResultVisitor)} for more
   * details.
   */
  public void findClosestEdges(S2BestDistanceTarget<D> target, ResultVisitor<D> visitor) {
    findBestEdges(target, visitor);
  }

  /**
   * Visits the closest (shape, edge, distance) Results to the given target that satisfy the current
   * options, and with shapes that pass the given ShapeFilter. This method may be called multiple
   * times, with the same or different targets, options, and ShapeFilters. Results are visited in
   * order of increasing distance from the target.
   *
   * <p>Note that if {@code options().includeInteriors()} is true, the results may include some
   * entries with edgeId == -1. This indicates that the target is contained by or intersects the
   * indexed polygon with the returned shapeId. Such results may also be identified by calling
   * {@link Result#isInterior()}.
   *
   * <p>The ShapeFilter may return different values for the same shapeId as it is repeatedly called.
   * For instance, the filter could accept a given shapeId the first time it is seen, and false
   * afterwards. If you only need to find the closest shapes to a target, this technique can speed
   * up the query significantly by returning fewer edges, especially for large indexes.
   *
   * <p>However, the ShapeFilter isn't called for every edge: a single passed filter test may
   * produce Results for many or even all the edges of that shape. Also, the Results produced would
   * not necessarily be the closest possible Results, as edges are not discovered in distance order
   * and the ShapeFilter is not called for edges in distance order. Finally, note that filtering
   * shapes and visiting results are not interleaved in the current implementation. First, edges and
   * shapes are discovered and filtered, and Results collected. Then the remaining Results are
   * visited in distance order.
   */
  public void findClosestEdges(
      S2BestDistanceTarget<D> target, @Nullable ShapeFilter shapeFilter, ResultVisitor<D> visitor) {
    findBestEdges(target, shapeFilter, visitor);
  }

  /**
   * Returns the closest (shape, edge, distance) Result to the target, if one satisfies the current
   * options. Otherwise, the returned Optional will not be present.
   *
   * <p>See {@link #findClosestEdges(S2BestDistanceTarget, ShapeFilter, ResultVisitor)} for more
   * details.
   */
  public Optional<Result<D>> findClosestEdge(S2BestDistanceTarget<D> target) {
    return findBestEdge(target);
  }

  /**
   * Returns the closest (shape, edge, distance) Result to the target, if one satisfies the current
   * options with a shape id that passes the given ShapeFilter. Otherwise, the returned Optional
   * will not be present.
   *
   * <p>See {@link #findClosestEdges(S2BestDistanceTarget, ShapeFilter, ResultVisitor)} for more
   * details.
   */
  public Optional<Result<D>> findClosestEdge(
      S2BestDistanceTarget<D> target, @Nullable ShapeFilter shapeFilter) {
    return findBestEdge(target, shapeFilter);
  }

  /**
   * Returns the minimum distance to the target. If the index or target is empty, or no edges
   * satisfy the current Options, returns a distance greater than any valid distance.
   *
   * <p>Use isDistanceLess() if you only want to compare the distance against a threshold value,
   * since it is often much faster.
   */
  public D getDistance(S2BestDistanceTarget<D> target) {
    return getDistance(target, null);
  }

  /**
   * Returns the minimum distance to the target from edges in the index satisfying the current
   * Options, and from shapes accepted by the given ShapeFilter. If the index or target is empty, or
   * no edges satisfy the current Options or ShapeFilter, returns a distance greater than any valid
   * distance.
   *
   * <p>Use isDistanceLess() if you only want to compare the distance against a threshold value,
   * since it is often much faster.
   */
  public D getDistance(S2BestDistanceTarget<D> target, @Nullable ShapeFilter shapeFilter) {
    Optional<Result<D>> result = findBestEdge(target, shapeFilter);
    return result.isPresent() ? result.get().distance() : beyondWorstDistance();
  }

  /**
   * Returns true if the distance to "target" from any edge in the index is less than "limit".
   *
   * <p>This method is usually much faster than getDistance(), since it is less work to determine
   * whether the minimum distance is above or below a threshold than it is to calculate the actual
   * minimum distance.
   */
  public boolean isDistanceLess(S2BestDistanceTarget<D> target, D limit) {
    return isDistanceLess(target, limit, null);
  }

  /**
   * Returns true if the distance to "target" from shapes accepted by the given ShapeFilter is less
   * than "limit".
   *
   * <p>This method is usually much faster than getDistance(), since it is less work to determine
   * whether the minimum distance is above or below a threshold than it is to calculate the actual
   * minimum distance.
   */
  public boolean isDistanceLess(S2BestDistanceTarget<D> target, D limit, ShapeFilter shapeFilter) {
    maxResults = 1;
    maxError = worstDistance();
    distanceLimit = limit;
    this.shapeFilter = shapeFilter;

    // Determine if any edge is found within the limit.
    findBestEdgesInternal(target);
    this.shapeFilter = null;
    boolean result = bestResult != null;
    bestResult = null;
    return result;
  }

  /** Returns the point on the given Result edge that is closest to the given targetPoint. */
  public S2Point project(S2Point targetPoint, Result<D> result) {
    if (result.edgeId() < 0) {
      return targetPoint;
    }
    S2Shape.MutableEdge resultEdge = new S2Shape.MutableEdge();
    index.getShapes().get(result.shapeId()).getEdge(result.edgeId(), resultEdge);
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

    /** Like {@link #setInclusiveMaxDistance(S1ChordAngle)} but takes an S1Angle for convenience. */
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
     * Like {@link #setConservativeMaxDistance(S1ChordAngle)} but takes an S1Angle for convenience.
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
     * the remaining search criteria. This option only has an effect if {@link #maxResults()} is
     * also specified; otherwise all edges that satisfy the maxDistance() will be returned.
     *
     * <p>Note that this does not affect how the distance between edges is computed; it simply
     * gives the algorithm permission to stop the search early, as soon as the best possible
     * improvement drops below {@link #maxError()}.
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

    /** Like {@link #setMaxError(S1ChordAngle)}, but maxError is provided as an S1Angle. */
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
    public S2ClosestEdgeQuery.Builder toBuilder() {
      return new S2ClosestEdgeQuery.Builder(options());
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
      // TODO(torrey): When Android supports Comparator.naturalOrder(), just return that.
      return (S1ChordAngle a, S1ChordAngle b) -> {
        return a.compareTo(b);
      };
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

    /**
     * Adjust the given 'value' towards zero by the current maximum error.
     *
     * <p>For {@link #isDistanceLess(S2BestDistanceTarget, S1Distance)}, maxError is STRAIGHT, so
     * the errorBoundedDistance will be zero for any given 'value'.
    */
    @Override
    protected S1ChordAngle errorBoundedDistance(S1ChordAngle value) {
      return S1ChordAngle.sub(value, maxError);
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
     return isDistanceLessOrEqual(target, limit, null);
    }

    /**
     * Like {@link #isDistanceLess(S2BestDistanceTarget, S1Distance, ShapeFilter)}, but also returns
     * true if the distance to "target" is exactly equal to "limit".
     */
    public boolean isDistanceLessOrEqual(
        Target<S1ChordAngle> target, S1ChordAngle limit, @Nullable ShapeFilter shapeFilter) {
      // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
      // not the same-named Options fields.
      distanceLimit = limit.successor();
      maxResults = 1;
      maxError = worstDistance();
      this.shapeFilter = shapeFilter;

      findBestEdgesInternal(target);
      this.shapeFilter = null;
      boolean result = bestResult != null;
      bestResult = null;

      return result;
    }

    /**
     * Like {@link #isDistanceLessOrEqual(Target, S1ChordAngle)}, except that "limit" is increased
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
      return isConservativeDistanceLessOrEqual(target, limit, null);
    }

    /**
     * Like {@link #isDistanceLessOrEqual(Target, S1ChordAngle, ShapeFilter)}, except that "limit"
     * is increased by the maximum error in the distance calculation. This ensures that this
     * function returns true whenever the true, exact distance is less than or equal to "limit".
     *
     * <p>For example, suppose that we want to test whether two geometries might intersect each
     * other after they are snapped together using S2Builder (using the IdentitySnapFunction with a
     * given "snapRadius"). Since S2Builder uses exact distance predicates, we need to measure the
     * distance between the two geometries conservatively. If the distance is definitely greater
     * than "snapRadius", then the geometries are guaranteed to not intersect after snapping.
     */
    public boolean isConservativeDistanceLessOrEqual(
        Target<S1ChordAngle> target, S1ChordAngle limit, @Nullable ShapeFilter shapeFilter) {
      // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
      // not the same-named Options fields.
      distanceLimit = limit.plusError(S2EdgeUtil.getMinDistanceMaxError(limit)).successor();
      maxResults = 1;
      maxError = worstDistance();
      this.shapeFilter = shapeFilter;

      findBestEdgesInternal(target);
      this.shapeFilter = null;
      boolean result = bestResult != null;
      bestResult = null;

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
   * Convenience method to create a ShapeIndexTarget using S1ChordAngle as the distance type. Shape
   * interiors are included, by default.
   *
   * <p>This is used for finding closest edges from a query ShapeIndex to a second ShapeIndex
   * wrapped by the ShapeIndexTarget.
   */
  public static ShapeIndexTarget<S1ChordAngle> createShapeIndexTarget(S2ShapeIndex index) {
    return new ShapeIndexTarget<>(index, new Builder());
  }
}
