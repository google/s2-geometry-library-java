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

import static java.lang.Math.sqrt;

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeUtil.PointVisitor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Optional;
import java.util.PriorityQueue;
import java.util.logging.Logger;
import javax.annotation.Nullable;

/**
 * S2BestEdgesQueryBase is an abstract class for finding the best edge(s) between two geometries. It
 * is not intended to be used directly, but rather to serve as the core implementation for various
 * specialized classes with more convenient APIs, such as S2ClosestEdgeQuery, which define what
 * "Best" means, and what kind of S1Distance is used for measuring distance.
 *
 * <p>It is flexible enough to be adapted to compute maximum distances and even potentially
 * Hausdorff distances. By using the appropriate options, implementations can answer questions such
 * as:
 *
 * <ul>
 *   <li>Find the minimum or maximum distance between two geometries A and B.
 *   <li>Find all edges of geometry A that are within a given distance 'd' of geometry B.
 *   <li>Find all edges of geometry A that are more than a given distance 'd' from geometry B.
 *   <li>Find the 'k' edges of geometry A that are closest to, or farthest from a given point P.
 * </ul>
 *
 * <p>You can also specify whether polygons should include their interiors (i.e., if a target point
 * is contained by an indexed polygon, should the distance be zero or should it be measured to the
 * polygon boundary?)
 *
 * <p>The input geometries may consist of any number of points, polylines, and polygons
 * (collectively referred to as "shapes"). Shapes do not need to be disjoint; they may overlap or
 * intersect arbitrarily. The implementation is designed to be fast for both simple and complex
 * geometries.
 *
 * <p>The S1Distance template argument defines the representation of distances. Usually it is
 * S1ChordAngle, but another implementation of the S1Distance interface may be used.
 *
 * <p>There are two versions of the findBestEdges() method. Users can simply get a List of Result
 * objects, but for cases where many calls to findBestEdges will be made for a single query, the
 * alternative can be much more efficient, as Results will be reused internally. Instead of
 * obtaining Result objects, clients provide a ResultVisitor that accepts (distance, shape, edge id)
 * tuples.
 */
@CheckReturnValue
public abstract class S2BestEdgesQueryBase<D extends S1Distance<D>> {
  private static final Logger log = Platform.getLoggerForClass(S2BestEdgesQueryBase.class);

  /**
   * When considering enqueing an index cell, if it has a sufficiently small number of edges, then
   * it is faster to check them directly rather than computing the minimum distance to the S2Cell,
   * enqueing it, and possibly later dequeing it.
   */
  // TODO(user): Adjust this using benchmarks.
  private static final int MIN_EDGES_TO_ENQUEUE = 10;

  /** Compares results by distance only, using the abstract distanceComparator(). */
  private final Comparator<Result<D>> resultDistanceComparator =
      Comparator.comparing((Result<D> result) -> result.distance, distanceComparator());

  /** Compares results by shape identity, breaking ties by edge id. */
  private final Comparator<Result<D>> resultShapeEdgeComparator =
      Comparator.comparingInt((Result<D> result) -> System.identityHashCode(result.shape))
          .thenComparingInt((Result<D> result) -> result.edgeId);

  /**
   * Compares Results by distance, then shape, then edge id. This breaks ties for distance in an
   * unpredictable but stable (per JVM instance) way.
   */
  private final Comparator<Result<D>> resultDistanceStableComparator =
      resultDistanceComparator.thenComparing(resultShapeEdgeComparator);

  // A queue of unprocessed S2CellIds, ordered by best distance to the target.
  private final PriorityQueue<QueueEntry<D>> queue =
      new PriorityQueue<>(
          Comparator.comparing((QueueEntry<D> arg) -> arg.distance, distanceComparator()));

  /**
   * Initial set of S2 cells to be searched for results, obtained by intersecting a covering of the
   * cap covering the area of acceptable distance to the target with the underlying S2ShapeIndex
   * cells.
   */
  private final ArrayList<S2CellId> initialCells = new ArrayList<>();

  /** The comparator of distances of a concrete type, providing the specific "best" ordering. */
  private final Comparator<D> distanceComparator = distanceComparator();

  /** The query options. */
  protected final Options<D> options;

  /** The geometry from which Results are found. */
  protected S2ShapeIndex index;

  /** The geometry we are measuring distance to. */
  protected S2BestDistanceTarget<D> target;

  /** Polygon shapes with interiors at the best possible distance to the target. */
  private final IdentityHashMap<S2Shape, Boolean> shapeInteriors = new IdentityHashMap<>();

  /**
   * This is a temporary value, containing the maximum number of results that should be returned for
   * the current query in progress. It is set to 1 for findBestEdge, otherwise is set from {@link
   * Options#maxResults()}.
   */
  protected int maxResults = Integer.MAX_VALUE;

  /**
   * This is a temporary value, containing the maximum error to be used for the current query in
   * progress. It is set to worstDistance() for isDistanceLess(), otherwise is set from {@link
   * Options#maxError()}.
   */
  protected D maxError;


  /**
   * This is a temporary value, computed from the target and current options. It is true if priority
   * queue cell distances must be adjusted by maxError in order to ensure that such distances are
   * measured conservatively. This is the case only if the target takes advantage of maxError in
   * order to return faster results, and maxError is between bestDistance() and the current
   * distanceLimit.
   */
  private boolean useConservativeCellDistance;

  /**
   * For the optimized algorithm we precompute the top-level S2CellIds that will be added to the
   * priority queue. There can be at most 6 of these cells. Essentially this is just a covering of
   * the indexed edges, except that we also store references to the corresponding Cell to reduce the
   * number of index seeks required. The covering needs to be stored in a List so that we can use
   * S2CellUnion.getIntersection().
   */
  private final ArrayList<S2CellId> indexCovering = new ArrayList<>(6);

  /**
   * These Cells correspond to the S2CellIds of indexCovering just above. If not null, these are
   * Cells in the underlying index that correspond to the top-level cells added to the priority
   * queue.
   */
  private final ArrayList<Cell> indexCells = new ArrayList<>(6);

  /**
   * The decision about whether to use the brute force algorithm is based on counting the total
   * number of edges in the index. However, if the index contains a large number of shapes, this in
   * itself might take too long. So we only count indexNumEdges up to maxBruteForceIndexSize() + 1
   * for the current target type, stored as indexNumEdgesLimit.
   */
  private int indexNumEdges;

  private int indexNumEdgesLimit;

  /**
   * The distance beyond which we can safely ignore further candidate edges: only edges with
   * distance better than this need to be considered. Candidates that are exactly at the limit are
   * ignored; this is more efficient for updateBestDistance() and should not affect clients since
   * distance measurements have a small amount of error anyway.
   *
   * <p>Initially this is the same as the minimum or maximum distance specified by the user, but it
   * can also be updated by the algorithm. See {@code maybeAddResult}.
   */
  protected D distanceLimit;

  /**
   * The accumulating result set is stored in a queue ordered by resultDistanceComparator reversed,
   * which places the worst distance result at the head of the queue.
   */
  protected final PriorityQueue<Result<D>> resultQueue =
      new PriorityQueue<>(resultDistanceComparator.reversed());

  /**
   * Some (shape, edge) pairs may be discovered as potential results multiple times, with the same
   * or different distances. As each potential Result is found, we check if it has already been
   * added to the resultQueue before computing the distance to it and adding it to the queue.
   */
  private final HashSet<ShapeEdgeId> testedEdges = new HashSet<>();

  /** A pool of Result objects that can be reused. */
  private final List<Result<D>> resultPool = new ArrayList<>();

  // S2BestEdgesQueryBase deals with "best" edges, e.g. edges with minimum or maximum distances to a
  // target. The following abstract methods define the meaning of "best" in their implementation:

  /** Returns a new DistanceCollector for best distances. */
  protected abstract DistanceCollector<D> newDistanceCollector();

  /**
   * Returns true if the given {@link DistanceCollector#distance()} has reached the best
   * possible distance, and therefore no further calls to update the DistanceCollector will change
   * the value.
   */
  protected abstract boolean atBestLimit(DistanceCollector<D> distanceCollector);

  /** Returns the distance zero. */
  protected abstract D zeroDistance();

  /** Returns the best valid distance a result can have. */
  protected abstract D bestDistance();

  /** Returns the worst valid distance. */
  protected abstract D worstDistance();

  /** Returns an invalid distance, worse than any valid result. */
  protected abstract D beyondWorstDistance();

  /**
   * Given a "bestCapRadius" of an S2Cap bounding the area of best possible distance to the target,
   * returns an enlarged (if necessary) cap radius such that all valid results with distance better
   * than the given "distanceLimit" can be found in the resulting search cap.
   */
  protected abstract S1ChordAngle searchCapRadius(S1ChordAngle bestCapRadius, D distanceLimit);

  /** Comparator for the concrete distance type and desired ordering. */
  protected abstract Comparator<D> distanceComparator();

  /**
   * Returns the "distance to beat" for new results by adjusting the given 'value' towards
   * bestDistance() by {@link maxError}.
   */
  protected abstract D errorBoundedDistance(D value);

  /**
   * Visits index shapes that are at the best distance to some point on a connected component of the
   * target. If the visitor returns false, no more shapes will be visited and this method returns
   * false. Otherwise all the best distance shapes will be visited for all connected components of
   * the target and this method returns true.
   */
  @CanIgnoreReturnValue
  protected abstract boolean visitBestDistanceContainingShapes(
      S2BestDistanceTarget<D> target, S2ContainsPointQuery.ShapeVisitor visitor);

  /**
   * A Result represents an edge of an S2Shape in the underlying index and its distance to a query
   * target. Note the following special case:
   *
   * <p>An {@code edgeId() < 0} represents the interior of a shape. Such results may be returned
   * when {@link Options#includeInteriors()} is true. Such results can be identified using the
   * {@link isInterior()} method.
   */
  public static class Result<D extends S1Distance<D>> {
    private D distance;
    private S2Shape shape;
    private int edgeId;

    /**
     * Constructs a Result for the given arguments.
     *
     * @param distance The distance from the indexed shape to the target.
     * @param shape The S2Shape in the underlying S2ShapeIndex of the query.
     * @param edgeId The edge of the indexed shape, or -1 if this Result represents the shape
     *     interior.
     */
    public Result(D distance, S2Shape shape, int edgeId) {
      this.distance = distance;
      this.shape = shape;
      this.edgeId = edgeId;
    }

    /** Updates this Result to have the given values. */
    public void set(D distance, S2Shape shape, int edgeId) {
      this.distance = distance;
      this.shape = shape;
      this.edgeId = edgeId;
    }

    /**
     * Returns true if this Result represents the interior of a shape. (Such results may be returned
     * when options.includeInteriors() is true.)
     */
    public boolean isInterior() {
      return shape != null && edgeId < 0;
    }

    /** The distance to the target from this shape edge. */
    public D distance() {
      return distance;
    }

    /** Returns the indexed shape. */
    public S2Shape shape() {
      return shape;
    }

    /** Identifies an edge within the shape. */
    public int edgeId() {
      return edgeId;
    }

    /**
     * Returns true if this Result is equal to the other, where equality means the same S2Shape (by
     * identity), the same edge id, and an exactly equal S1Distance of the same type.
     */
    public boolean equalsResult(Result<D> other) {
      return (this.distance.equals(other.distance)
          && this.shape == other.shape
          && this.edgeId == other.edgeId);
    }
  }

  /**
   * A ResultVisitor implements an accept() method with parameters of distance, index shape, and
   * shape edge id.
   */
  public interface ResultVisitor<D> {
    /**
     * Each call to accept() indicates that the given 'edge' of the given 'shape' in the query index
     * is at the given 'distance' to the query target.
     */
    boolean accept(D distance, S2Shape shape, int edgeId);
  }

  /**
   * The base class Builder. Subclasses provide additional constructors and setters.
   */
  public abstract static class Builder<D extends S1Distance<D>> {
    protected int maxResults;
    protected D distanceLimit;
    protected D maxError;
    protected boolean includeInteriors;
    protected boolean useBruteForce;

    /**
     * Constructs a Builder with default options. The default distance limit and value of zero must
     * be supplied by the subclass.
     */
    protected Builder(D defaultDistanceLimit, D zero) {
      distanceLimit = defaultDistanceLimit;
      maxError = zero;
      maxResults = Integer.MAX_VALUE;
      includeInteriors = true;
      useBruteForce = false;
    }

    /** Constructs a Builder with option values copied from the given Options. */
    public Builder(Options<D> options) {
      distanceLimit = options.distanceLimit;
      maxError = options.maxError;
      maxResults = options.maxResults;
      includeInteriors = options.includeInteriors;
      useBruteForce = options.useBruteForce;
    }

    /** Copy constructor. */
    public Builder(Builder<D> other) {
      distanceLimit = other.distanceLimit;
      maxError = other.maxError;
      maxResults = other.maxResults;
      includeInteriors = other.includeInteriors;
      useBruteForce = other.useBruteForce;
    }

    /** Builds a new query using the current option values. */
    public abstract S2BestEdgesQueryBase<D> build();

    /** Create new query on the given index, using the current option values. */
    public S2BestEdgesQueryBase<D> build(S2ShapeIndex index) {
      S2BestEdgesQueryBase<D> query = build();
      query.init(index);
      return query;
    }

    /** Specifies the maximum number of results to be returned. */
    @CanIgnoreReturnValue
    public Builder<D> setMaxResults(int maxResults) {
      this.maxResults = maxResults;
      return this;
    }

    /** Returns the current value of the maxResults option in this Builder. */
    public int maxResults() {
      return maxResults;
    }

    /** All returned results must have a distance better than the distanceLimit. */
    @CanIgnoreReturnValue
    public Builder<D> setDistanceLimit(D distanceLimit) {
      this.distanceLimit = distanceLimit;
      return this;
    }

    /** Returns the current value of the distanceLimit option in this Builder. */
    public D distanceLimit() {
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
     * whether the minimum distance is less than minDist, set {@code maxResults == 1} and {@code
     * maxDistance() == maxError() == minDist}. This causes the algorithm to terminate as soon as
     * it finds any edge whose distance is less than minDist, rather than continuing to search for
     * an edge that is even closer.
     */
    @CanIgnoreReturnValue
    public Builder<D> setMaxError(D maxError) {
      this.maxError = maxError;
      return this;
    }

    /** Returns the current value of the maxError option in this Builder. */
    public D maxError() {
      return maxError;
    }

    /**
     * True if polygon interiors in the queried index are considered when computing distance.
     *
     * <p>When true, polygons that contain the target have a distance of zero. (For targets
     * consisting of multiple connected components, the distance is zero if any component is
     * contained.) This is indicated in the results by returning a (shape, edgeId) pair with
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
    public Builder<D> setIncludeInteriors(boolean includeInteriors) {
      this.includeInteriors = includeInteriors;
      return this;
    }

    /** Returns the current value of the includeInteriors option in this Builder. */
    public boolean includeInteriors() {
      return includeInteriors;
    }

    /**
     * If true, distances should be computed by examining every edge rather than using the
     * S2ShapeIndex. This is useful for testing, benchmarking, and debugging.
     */
    @CanIgnoreReturnValue
    public Builder<D> setUseBruteForce(boolean useBruteForce) {
      this.useBruteForce = useBruteForce;
      return this;
    }

    /** Returns the current value of the useBruteForce option in this Builder. */
    public boolean useBruteForce() {
      return useBruteForce;
    }
  }

  /**
   * Options control the set of edges returned by S2BestEdgesQueryBase. You will usually want to
   * provide Options that set a {@link distanceLimit()}, or limit {@link maxResults()} or both,
   * otherwise all edges will be returned.
   *
   * <p>Options are intended to be immutable, and created from a Builder.
   */
  public static final class Options<D extends S1Distance<D>> {
    final int maxResults;
    final D distanceLimit;
    final D maxError;
    final boolean includeInteriors;
    final boolean useBruteForce;

    /**
     * Options may only be constructed from a Builder. Clients should use {@link Builder#build()}.
     */
    Options(Builder<D> b) {
      maxResults = b.maxResults();
      distanceLimit = b.distanceLimit();
      maxError = b.maxError();
      includeInteriors = b.includeInteriors();
      useBruteForce = b.useBruteForce();
    }

    /**
     * maxResults specifies the maximum number of results to be returned by a query.
     *
     * <p>maxResults must be greater than or equal to 1, and the default value set by
     * implementations should normally be Integer.MAX_VALUE.
     */
    public int maxResults() {
      return maxResults;
    }

    /**
     * All returned results must have a distance better than the distanceLimit.
     *
     * <p>The default value set by implementations should normally be the most permissive possible
     * distance.
     *
     * <p>Note that edges whose distance is exactly equal to "distanceLimit" are not returned. In
     * most cases this doesn't matter, since distances are not computed exactly in the first place,
     * but if such edges are needed then you can retrieve them by specifying "distanceLimit" as the
     * next largest or smallest representable S1Distance.
     */
    public D distanceLimit() {
      return distanceLimit;
    }

    /**
     * A non-zero maxError specifies that edges with distance up to maxError worse than the true
     * best-distance edges may be substituted in the result set, as long as such edges satisfy all
     * the remaining search criteria. This option only has an effect if {@link maxResults()} is also
     * specified; otherwise all edges that satisfy the distanceLimit() will be returned.
     *
     * <p>Note that this does not affect how the distance between edges is computed; it simply gives
     * the algorithm permission to stop the search early, as soon as the best possible improvement
     * drops below {@link maxError()}.
     *
     * <p>This can be used to implement distance predicates efficiently. For example, to determine
     * whether the minimum distance is less than minDist, set {@code maxResults == 1} and {@code
     * maxDistance() == maxError() == minDist}. This causes the algorithm to terminate as soon as it
     * finds any edge whose distance is less than minDist, rather than continuing to search for an
     * edge that is even closer.
     *
     * <p>The default set by implementations should normally be zero of type D.
     */
    public D maxError() {
      return maxError;
    }

    /**
     * True if polygon interiors in the queried index are considered when computing distance.
     *
     * <p>When true, polygons that contain the target have a distance of zero. (For targets
     * consisting of multiple connected components, the distance is zero if any component is
     * contained.) This is indicated in the results by returning a (shape, edgeId) pair with {@code
     * edgeId == -1}, i.e. this value denotes the polygons's interior.
     *
     * <p>Note that this does not control if the interiors of _target_ polygons should be included;
     * that is a separate option on targets.
     *
     * <p>Note that for efficiency, any polygon that intersects the target may or may not have an
     * (edgeId == -1) result. Such results are optional because in that case the distance from the
     * polygon is already zero.
     *
     * <p>The default set by implementations should normally be true.
     */
    public boolean includeInteriors() {
      return includeInteriors;
    }

    /**
     * If true, distances should be computed by examining every edge rather than using the
     * S2ShapeIndex. This is useful for testing, benchmarking, and debugging.
     *
     * <p>The default set by implementations should normally be false.
     */
    public boolean useBruteForce() {
      return useBruteForce;
    }
  }

  /** A base class for targets that are points. */
  protected abstract static class PointTarget<D extends S1Distance<D>>
      implements S2BestDistanceTarget<D> {
    protected final S2Point point;

    /**
     * Constructs a new PointTarget.
     *
     * @param point the {@link S2Point} defining this target.
     */
    public PointTarget(S2Point point) {
      this.point = point;
    }

    /**
     * If the distance to this target point from the given point {@code p} is better than the
     * current value in {@code collector}, update {@code collector} and return true. Otherwise
     * return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point p, DistanceCollector<D> collector) {
      return collector.update(point, p);
    }

    /**
     * If the distance to this target point from the given edge {@code (v0, v1)} is better than the
     * current value in {@code collector}, update {@code collector} and return true. Otherwise
     * return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point v0, S2Point v1, DistanceCollector<D> collector) {
      return collector.update(point, v0, v1);
    }

    /**
     * If the distance to this target point from the given {@link S2Cell} {@code cell} is better
     * than the current value in {@code collector}, update {@code collector} and return true.
     * Otherwise return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Cell cell, DistanceCollector<D> collector) {
      return collector.update(point, cell);
    }

    @CanIgnoreReturnValue
    @Override
    public boolean visitConnectedComponentPoints(PointVisitor visitor) {
      return visitor.apply(point);
    }
  }

  /** A base class for targets that are edges. */
  protected abstract static class EdgeTarget<D extends S1Distance<D>>
      implements S2BestDistanceTarget<D> {
    protected final S2Point a;
    protected final S2Point b;

    /**
     * Constructs a new EdgeTarget.
     *
     * @param a One endpoint of the edge defining this target.
     * @param b The other endpoint of the edge defining this target.
     */
    public EdgeTarget(S2Point a, S2Point b) {
      this.a = a.normalize();
      this.b = b.normalize();
    }

    /**
     * If the distance to this target edge from the given point {@code p} is better than the current
     * value in {@code collector}, update {@code collector} and return true. Otherwise return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point p, DistanceCollector<D> collector) {
      return collector.update(p, a, b);
    }

    /**
     * If the distance to this target edge from the given edge {@code (v0, v1)} is better than the
     * current value in {@code collector}, update {@code collector} and return true. Otherwise
     * return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point v0, S2Point v1, DistanceCollector<D> collector) {
      return collector.update(v0, v1, a, b);
    }

    /**
     * If the distance to this target edge from the given {@link S2Cell} {@code cell} is better than
     * the current value in {@code collector}, update {@code collector} and return true. Otherwise
     * return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Cell cell, DistanceCollector<D> collector) {
      return collector.update(a, b, cell);
    }

    /**
     * Computes the square of half the edge length in an efficient and numerically stable way. For
     * constructing bounding caps of the edge.
     */
    protected double getHalfEdgeLength2() {
      double d2 = new S1ChordAngle(a, b).getLength2();
      return (0.5 * d2) / (1 + sqrt(1 - 0.25 * d2));
    }

    @CanIgnoreReturnValue
    @Override
    public boolean visitConnectedComponentPoints(PointVisitor visitor) {
      // We visit the center of the edge in order to ensure that edge targets AB and BA yield
      // identical results (which is not guaranteed by the API but users might expect).
      S2Point edgeCenter = a.add(b).normalize();
      return visitor.apply(edgeCenter);
    }
  }

  /** A base class for targets that are S2 cells. */
  protected abstract static class CellTarget<D extends S1Distance<D>>
      implements S2BestDistanceTarget<D> {
    protected final S2Cell cell;

    /**
     * Constructs a new S2DistanceCellTarget.
     *
     * @param cell The {@link S2CellId} defining this target.
     */
    public CellTarget(S2Cell cell) {
      this.cell = cell;
    }

    /**
     * If the distance to this target cell from the given point {@code p} is better than the current
     * value in {@code collector}, update {@code collector} and return true. Otherwise return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point p, DistanceCollector<D> collector) {
      return collector.update(p, cell);
    }

    /**
     * If the distance to this target cell from the given edge {@code (v0, v1)} is better than the
     * current value in {@code collector}, update {@code collector} and return true. Otherwise
     * return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point v0, S2Point v1, DistanceCollector<D> collector) {
      return collector.update(v0, v1, cell);
    }

    /**
     * If the distance to this target cell from the given {@link S2Cell} {@code c} is better than
     * the current value in {@code collector}, update {@code collector} and return true. Otherwise
     * return false.
     */
    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Cell c, DistanceCollector<D> collector) {
      return collector.update(cell, c);
    }

    @CanIgnoreReturnValue
    @Override
    public boolean visitConnectedComponentPoints(PointVisitor visitor) {
      return visitor.apply(cell.getCenter());
    }
  }

  /**
   * The algorithm maintains a priority queue of unprocessed S2CellIds, sorted in order with best
   * cells first.
   */
  protected static class QueueEntry<D extends S1Distance<D>> {
    /**
     * A bound on the distance from the cell "id" to the target. This is the key of the priority
     * queue.
     */
    final D distance;

    /** The cell being queued. */
    final S2CellId id;

    /**
     * If "id" belongs to the index, this field stores the corresponding S2ShapeIndex.Cell.
     * Otherwise "id" is a proper ancestor of one or more S2ShapeIndexCells and this field is null.
     * The purpose of this field is to avoid an extra seek() when the queue entry is processed.
     */
    final Cell indexCell;

    QueueEntry(D distance, S2CellId id, Cell indexCell) {
      Preconditions.checkNotNull(distance);
      this.distance = distance;
      this.id = id;
      this.indexCell = indexCell;
    }
  }

  /**
   * ShapeEdgeId is a unique identifier for an edge within an {@link S2ShapeIndex}, consisting of a
   * shape and edgeId pair. These are stored in the testedEdges HashSet.
   *
   * TODO(torrey): Could this be optimized by encoding the shape and edge id into a Long? This may
   * require ShapeIndex to provide shape ids.
   */
  private static class ShapeEdgeId {
    private final S2Shape shape;
    private final int edgeId;

    /** Constructs a new ShapeEdgeId for the given shape and edge id. */
    public ShapeEdgeId(S2Shape shape, int edgeId) {
      this.shape = shape;
      this.edgeId = edgeId;
    }

    @Override
    public int hashCode() {
      return System.identityHashCode(shape) + (17 * edgeId);
    }

    @Override
    public boolean equals(Object obj) {
      if (obj == this) {
        return true;
      }
      if (!(obj instanceof ShapeEdgeId)) {
        return false;
      }
      return (shape.equals(((ShapeEdgeId) obj).shape) && edgeId == ((ShapeEdgeId) obj).edgeId);
    }
  }

  /** The constructor is internal, to be used by the Builder. Requires init() to be called. */
  S2BestEdgesQueryBase(Options<D> options) {
    this.options = options;
  }

  /** Initializes the query. Note that reInit() must be called if {@code index} is modified. */
  public void init(S2ShapeIndex index) {
    this.index = index;
    reInit();
  }

  /** Reinitializes the query. This method must be called whenever the {@code index} is modified. */
  public void reInit() {
    indexNumEdges = 0;
    indexNumEdgesLimit = 0;
    indexCells.clear();
  }

  /** Returns the underlying S2ShapeIndex. */
  public S2ShapeIndex index() {
    return index;
  }

  /** Get a reference to the current options. Options are immutable. */
  public Options<D> options() {
    return options;
  }

  /** True if the given distance has reached or passed the current distance limit. */
  protected boolean distanceAtLimit(D distance) {
    return distanceComparator.compare(distance, distanceLimit) >= 0;
  }

  /**
   * Convenience method that returns the one best edge, if one is found that satisfies the current
   * search options. Otherwise, if none satisfies the options, the {@code Optional<Result>} will be
   * empty.
   */
  public Optional<Result<D>> findBestEdge(S2BestDistanceTarget<D> target) {
    // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
    // not the same-named Options fields.
    distanceLimit = options.distanceLimit();
    maxResults = 1;
    maxError = options.maxError;

    // Find best edge and put it in resultQueue.
    findBestEdgesInternal(target);

    Optional<Result<D>> result =
        resultQueue.isEmpty() ? Optional.empty() : Optional.of(resultQueue.poll());
    resultQueue.clear();
    return result;
  }

  /**
   * Returns true if the distance to "target" is better than the current value in the given
   * DistanceCollector, and updates the DistanceCollector with that smaller distance. Otherwise
   * returns false.
   */
  public boolean updateBestDistance(
      S2BestDistanceTarget<D> target, DistanceCollector<D> collector) {
    // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
    // not the same-named Options fields.
    maxResults = 1;
    maxError = options.maxError;
    distanceLimit = collector.distance();

    findBestEdgesInternal(target);
    if (resultQueue.isEmpty()) {
      resultQueue.clear();
      return false;
    }

    Result<D> result = resultQueue.peek();
    resultPool.add(result);
    resultQueue.clear();
    return collector.update(result.distance());
  }

  /**
   * Returns the best-distance shape edges from the query index to the given target that satisfy the
   * current options. This method may be called multiple times. If no shape edges satisfy the
   * current options, an empty list is returned.
   *
   * <p>Note that if {@link Options#includeInteriors()} is true, the results may include some
   * entries with {@code edgeId == -1}. This indicates that some point on a connected component of
   * the target is at the best possible distance from the interior of the {@link Result#shape()}.
   *
   * <p>If the target has s2 cells, or polygons with
   * {@link S2BestDistanceTarget#includeInteriors()}, then the results will include entries for any
   * indexed shape edges that have the best possible distance target polygon or s2 cell interiors.
   */
  public List<Result<D>> findBestEdges(S2BestDistanceTarget<D> target) {
    // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
    // not the same-named Options fields.
    distanceLimit = options.distanceLimit();
    maxResults = options.maxResults();
    maxError = options.maxError;

    // Find best edges and put them in resultQueue, reusing Results from the pool.
    findBestEdgesInternal(target);

    // Copy the resultQueue into a list and sort it to return.
    List<Result<D>> results = new ArrayList<>(resultQueue);
    Collections.sort(results, resultDistanceStableComparator);
    resultQueue.clear();
    return results;
  }

  /**
   * This version of findBestEdges calls the provided ResultVisitor with the best (shape, edge)
   * pairs and their distances to the target.
   */
  public void findBestEdges(S2BestDistanceTarget<D> target, ResultVisitor<D> visitor) {
    // Note that from here on down, the distanceLimit, maxResults, and maxError fields are used,
    // not the same-named Options fields.
    distanceLimit = options.distanceLimit();
    maxResults = options.maxResults();
    maxError = options.maxError;

    // Find best edges and put them in resultQueue, reusing Results from the pool.
    findBestEdgesInternal(target);

    // Copy the resultQueue into a list and sort it before visiting.
    List<Result<D>> results = new ArrayList<>(resultQueue);
    Collections.sort(results, resultDistanceStableComparator);
    resultQueue.clear();

    for (Result<D> result : results) {
      if (!visitor.accept(result.distance, result.shape, result.edgeId)) {
        break;
      }
    }
    // Return Results to the pool for reuse.
    resultPool.addAll(results);
  }

  /**
   * Finds the best distance shape edges from the query index to the given target that satisfy the
   * given values of the temporary maxResults, maxError, and distanceLimit fields. Results matching
   * those (and options.includeInteriors and options.useBruteForce) are added to the resultQueue.
   */
  protected void findBestEdgesInternal(S2BestDistanceTarget<D> target) {
    // assert resultQueue.isEmpty();
    // assert target.maxBruteForceIndexSize() >= 0;

    this.target = target;
    testedEdges.clear();

    // If the distance limit is the best possible distance, no results can be added.
    if (distanceLimit.equals(bestDistance())) {
      log.fine("No possible results, cannot be better than the best possible distance.");
      return;
    }

    if (maxResults == Integer.MAX_VALUE && distanceLimit.equals(beyondWorstDistance())) {
      log.fine("Returning all edges! maxResults and distanceLimit set no limits.");
    }

    // If polygon interiors in the query index are included, then we collect indexed polygons with
    // interiors at the best possible distance to a point on a connected component of the target.
    if (options.includeInteriors()) {
      shapeInteriors.clear();
      visitBestDistanceContainingShapes(
          target,
          (shape) -> {
            if (!shapeInteriors.containsKey(shape)) {
              shapeInteriors.put(shape, true);
              addResult(bestDistance(), shape, -1);

              // If the distance limit is the best possible, no more results can be added.
              if (distanceLimit.equals(bestDistance())) {
                return false;
              }
            }
            return true;
          });
    }

    // If maxError > 0 and the target takes advantage of this, then we may need to adjust the
    // distance estimates to the priority queue cells to ensure that they are always a bound on the
    // true distance in the "best" direction. For example, suppose we are finding closest edges
    // with distanceLimit == 70, maxError == 30, and we compute the distance to the target from
    // some cell C0 as d(C0) == 80. Because the target takes advantage of maxError, the true
    // distance could be as low as 50. In order not to miss edges contained by such cells, we need
    // to adjust the distance estimates by maxError. This behavior is controlled by the
    // useConservativeCellDistance value.
    //
    // However there is one important case where this adjustment is not necessary. (The following
    // explanation assumes a search for closest edges, but the situation is analogous for furthest
    // edges). This case is when distanceLimit < maxError. This is because maxError only affects the
    // algorithm once at least maxResults edges have been found that satisfy the given distance
    // limit. At that point, maxError is subtracted from distanceLimit in order to ensure that any
    // further matches are closer by at least that amount. But when distanceLimit < maxError, this
    // reduces the distance limit to 0, i.e. all remaining candidate cells and edges can safely be
    // discarded. (Note that this is how isDistanceLess() and friends are implemented.)
    //
    // Note that S2BestDistanceTarget.setMaxError() returns true if the implementation takes
    // advantage of maxError, and false otherwise.
    boolean targetUsesMaxError = !maxError.isZero() && target.setMaxError(maxError);

    useConservativeCellDistance =
        targetUsesMaxError
            && (distanceLimit.equals(beyondWorstDistance())
                || maxError.lessThan(distanceLimit));

    // Use the brute force algorithm if the index is small enough. To avoid spending too much time
    // counting edges and points when there are many shapes, we stop counting once there are too
    // many. We may need to recount if we later see a target with a larger brute force threshold.
    // Note that points are represented and counted as degenerate edges.
    int minOptimizedEdges = target.maxBruteForceIndexSize() + 1;
    if (minOptimizedEdges > indexNumEdgesLimit && indexNumEdges >= indexNumEdgesLimit) {
      indexNumEdges = S2ShapeUtil.countEdgesUpTo(index, minOptimizedEdges);
      indexNumEdgesLimit = minOptimizedEdges;
    }

    if (options.useBruteForce() || indexNumEdges < minOptimizedEdges) {
      findBestEdgesBruteForce();
    } else {
      findBestEdgesOptimized();
    }
  }

  private void findBestEdgesBruteForce() {
    for (S2Shape shape : index.getShapes()) {
      if (shape == null) {
        continue;
      }
      int numEdges = shape.numEdges();
      for (int e = 0; e < numEdges; ++e) {
        maybeAddResult(shape, e);
      }
    }
  }

  private void findBestEdgesOptimized() {
    S2Iterator<Cell> iter = initQueue();

    // Repeatedly find the best-distance S2Cell to "target" and either split it into its four
    // children, or process all of its edges.
    while (!queue.isEmpty()) {
      QueueEntry<D> entry = queue.poll();

      // If the best (nearest / furthest) cell to the target is past the distance limit, we're done.
      if (distanceAtLimit(entry.distance)) {
        queue.clear(); // Clear any remaining entries.
        break;
      }

      // If this is already known to be an index cell, just process it.
      if (entry.indexCell != null) {
        processEdges(entry.indexCell);
        continue;
      }

      // Otherwise, split the cell into its four children. Before adding a child back to the queue,
      // we first check whether it is empty. We do this in two seek operations rather than four by
      // seeking to the key between children 0 and 1, and to the key between children 2 and 3.
      S2CellId id = entry.id;
      iter.seek(id.child(1).rangeMin());
      if (!iter.done() && iter.id().lessOrEquals(id.child(1).rangeMax())) {
        processOrEnqueueChild(id.child(1), iter);
      }
      if (iter.prev() && iter.id().greaterOrEquals(id.rangeMin())) {
        processOrEnqueueChild(id.child(0), iter);
      }
      iter.seek(id.child(3).rangeMin());
      if (!iter.done() && iter.id().lessOrEquals(id.child(3).rangeMax())) {
        processOrEnqueueChild(id.child(3), iter);
      }
      if (iter.prev() && iter.id().greaterOrEquals(id.child(2).rangeMin())) {
        processOrEnqueueChild(id.child(2), iter);
      }
    }
  }

  /** Enqueue the given cell id. Requires that iter is positioned at a cell contained by "id". */
  private void processOrEnqueueChild(S2CellId id, S2Iterator<Cell> iter) {
    // assert id.contains(iter.id());
    if (iter.id().equals(id)) {
      processOrEnqueue(id, iter.entry());
    } else {
      processOrEnqueue(id, null);
    }
  }

  private S2Iterator<Cell> initQueue() {
    // assert queue.isEmpty();
    S2Iterator<Cell> iter = index.iterator();

    // Optimization: if the user is searching for just one best edge, and the center of the target's
    // bounding cap happens to intersect an index cell, then we try to limit the search region to a
    // small disc by first processing the edges in that cell. This sets distanceLimit based on the
    // best edge in that cell, which we can then use to limit the search area. This means that the
    // cell containing "target" will be processed twice, but in general this is still faster.
    //
    // TODO(user): Even if the cap center is not contained, we could still process one or
    // both of the adjacent index cells in S2CellId order, provided that those cells are better than
    // distanceLimit.

    // Obtain a bound on points with the best possible distance to the target.
    S2Cap targetCap = target.getCapBound();
    if (targetCap.isEmpty()) {
      // If the target has no points, no results are possible.
      return iter;
    }
    if (maxResults == 1 && iter.locate(targetCap.axis())) {
      processEdges(iter.entry());
      // Skip the rest of the algorithm if we found an intersecting edge and only need one result.
      if (distanceLimit.equals(bestDistance())) {
        return iter;
      }
    }

    if (indexCovering.isEmpty()) {
      initCovering();
    }

    // If the distance limit allows for results at any distance, use the precomputed covering of
    // the whole index.
    if (distanceLimit == beyondWorstDistance()) {
      for (int i = 0; i < indexCovering.size(); ++i) {
        processOrEnqueue(indexCovering.get(i), indexCells.get(i));
      }
    } else {
      // Otherwise, we compute a covering of the search disc and intersect it with the precomputed
      // index covering.
      ArrayList<S2CellId> searchCapCovering = new ArrayList<>();
      initialCells.clear();
      S1ChordAngle radius = searchCapRadius(targetCap.radius(), distanceLimit);
      S2Cap searchCap = S2Cap.fromAxisChord(targetCap.axis(), radius);
      S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(4).build();
      coverer.getFastCovering(searchCap, searchCapCovering);
      S2CellUnion.getIntersection(indexCovering, searchCapCovering, initialCells);

      // Now we need to clean up the initial cells to ensure that they all contain at least one cell
      // of the S2ShapeIndex. (Some may not intersect the index at all, while others may be
      // descendants of an index cell.)
      for (int i = 0, j = 0; i < initialCells.size(); ) {
        S2CellId initialCellId = initialCells.get(i);
        // Find the top-level cell that contains this initial cell.
        while (indexCovering.get(j).rangeMax().lessThan(initialCellId)) {
          ++j;
        }
        S2CellId coveringCellId = indexCovering.get(j);

        if (initialCellId.equals(coveringCellId)) {
          // This initial cell is one of the top-level cells. Use the precomputed S2ShapeIndex.Cell
          // to avoid an index seek.
          processOrEnqueue(coveringCellId, indexCells.get(j));
          ++i;
          ++j;
        } else {
          // This initial cell is a proper descendant of a top-level cell. Check how it is related
          // to the cells of the S2ShapeIndex.
          S2ShapeIndex.CellRelation r = iter.locate(initialCellId);
          if (r == S2ShapeIndex.CellRelation.INDEXED) {
            // This cell is a descendant of an index cell. Enqueue it and skip any other initial
            // cells that are also descendants of this cell.
            processOrEnqueue(iter.id(), iter.entry());
            S2CellId lastId = iter.id().rangeMax();
            while (++i < initialCells.size() && initialCells.get(i).lessOrEquals(lastId)) {
              continue;
            }
          } else {
            // Enqueue the cell only if it contains at least one index cell.
            if (r == S2ShapeIndex.CellRelation.SUBDIVIDED) {
              processOrEnqueue(initialCellId, null);
            }
            ++i;
          }
        }
      }
    }
    return iter;
  }

  /**
   * Find the range of S2Cells spanned by the index and choose a level such that the entire index
   * can be covered with just a few cells. These are the "top-level" cells. There are two cases:
   *
   * <p>If the index spans more than one face, then there is one top-level cell per spanned face,
   * just big enough to cover the index cells on that face.
   *
   * <p>If the index spans only one face, then we find the smallest cell "C" that covers the index
   * cells on that face (just like the case above). Then for each of the 4 children of "C", if the
   * child contains any index cells then we create a top-level cell that is big enough to just fit
   * those index cells (i.e., shrinking the child as much as possible to fit its contents). This
   * essentially replicates what would happen if we started with "C" as the top-level cell, since
   * "C" would immediately be split, except that we take the time to prune the children further
   * since this will save work on every subsequent query.
   */
  private void initCovering() {
    indexCovering.clear();
    indexCells.clear();

    // TODO(user): Have a single implementation of the algorithm for covering a shape index.

    // TODO(user): Use a single iterator below and save position information using a
    // Pair<S2CellId, S2ShapeIndexCell> type?

    S2Iterator<Cell> next = index.iterator();
    if (next.done()) {
      // The index is empty.
      return;
    }

    // TODO(torrey): index should have an iteratorLast() or equivalent, this is clumsy.
    S2Iterator<Cell> last = index.iterator();
    last.finish();
    last.prev();

    if (!next.id().equals(last.id())) {
      // The index has at least two cells. Choose a level such that the entire index can be spanned
      // with at most 6 cells (if the index spans multiple faces) or 4 cells (if the index spans a
      // single face).
      int level = next.id().getCommonAncestorLevel(last.id()) + 1;

      // Visit each potential top-level cell except the last (handled below).
      S2CellId lastId = last.id().parent(level);
      for (S2CellId id = next.id().parent(level); !id.equals(lastId); id = id.next()) {
        // Skip any top-level cells that don't contain any index cells.
        if (id.rangeMax().lessThan(next.id())) {
          continue;
        }
        // Find the range of index cells contained by this top-level cell and then shrink the cell
        // if necessary so that it just covers them.
        S2Iterator<Cell> cellFirst = S2Iterator.copy(next);
        next.seek(id.rangeMax().next());
        S2Iterator<Cell> cellLast = S2Iterator.copy(next);
        cellLast.prev();
        addInitialRange(cellFirst, cellLast);
      }
    }
    addInitialRange(next, last);
  }

  /**
   * Add an entry to indexCovering and indexCells that covers the given inclusive range of cells.
   * Requires that "first" and "last" have a common ancestor.
   */
  private void addInitialRange(S2Iterator<Cell> first, S2Iterator<Cell> last) {
    if (first.id().equals(last.id())) {
      // The range consists of a single index cell.
      indexCovering.add(first.id());
      indexCells.add(first.entry());
    } else {
      // Add the lowest common ancestor of the given range.
      int level = first.id().getCommonAncestorLevel(last.id());
      // assert level >= 0; // Ensure that they do in fact have a common ancestor.
      indexCovering.add(first.id().parent(level));
      indexCells.add(null);
    }
  }

  private void maybeAddResult(S2Shape shape, int edgeId) {
    if (!testedEdges.add(new ShapeEdgeId(shape, edgeId))) {
      return;
    }

    // TODO(user): If we knew when the shape edges are actually points, we could call
    // updateBestDistance(point, distance) instead, which might be meaningfully more efficient.
    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
    shape.getEdge(edgeId, edge);

    // Determine if the distance to the edge is better than the limit beyond which edges may be
    // ignored. If so, add (or visit) a new result with that edge and distance.
    DistanceCollector<D> collector = newDistanceCollector();
    collector.set(distanceLimit);
    if (target.updateBestDistance(edge.a, edge.b, collector)) {
      addResult(collector.distance(), shape, edgeId);
    }
  }

  // Add the new result to the queue, and update the distance limit based on the worst-distance
  // result remaining in the queue.
  private void addResult(D distance, S2Shape shape, int edgeId) {
    Result<D> result;
    if (resultPool.isEmpty()) {
      result = new Result<>(distance, shape, edgeId);
    } else {
      // Remove and reuse a result object from the pool.
      result = resultPool.remove(resultPool.size() - 1);
      result.set(distance, shape, edgeId);
    }
    resultQueue.add(result);

    int size = resultQueue.size();
    if (size >= maxResults) {
      // If we have more results than needed, remove the worst-distance one from the queue and
      // return it to the pool.
      if (size > maxResults) {
        resultPool.add(resultQueue.poll());
      }
      // Once we've reached the required number of results, keep updating the distance beyond
      // which candidate results may be ignored.
      distanceLimit = errorBoundedDistance(resultQueue.peek().distance());
    }
  }

  // Return the number of edges in the given index cell.
  // TODO(torrey): This should be a method on S2ShapeIndex.Cell. Same for C++.
  private static int countEdges(Cell cell) {
    int count = 0;
    for (int s = 0; s < cell.numShapes(); ++s) {
      count += cell.clipped(s).numEdges();
    }
    return count;
  }

  // Process all the edges of the given index cell.
  private void processEdges(Cell shapeIndexCell) {
    for (int s = 0; s < shapeIndexCell.numShapes(); ++s) {
      S2ShapeIndex.S2ClippedShape clipped = shapeIndexCell.clipped(s);
      S2Shape shape = clipped.shape();
      for (int j = 0; j < clipped.numEdges(); ++j) {
        maybeAddResult(shape, clipped.edge(j));
      }
    }
  }

  /**
   * Add the given cell id to the queue. "indexCell" is the corresponding S2ShapeIndexCell, or null
   * if "id" is not an index cell.
   *
   * <p>This version is called directly only by initQueue().
   */
  private void processOrEnqueue(S2CellId id, @Nullable Cell indexCell) {
    if (indexCell != null) {
      // If this index cell has only a few edges, then it is faster to check them directly rather
      // than computing the minimum distance to the S2Cell and inserting it into the queue.
      int numEdges = countEdges(indexCell);
      if (numEdges == 0) {
        return;
      }
      if (numEdges < MIN_EDGES_TO_ENQUEUE) {
        processEdges(indexCell);
        return;
      }
    }

    // Otherwise, determine if the distance to the cell is better than the limit at which we can
    // ignore candidate edges. If so, the cell may have viable candidates, so add the cell to the
    // priority queue.
    S2Cell cell = new S2Cell(id);
    DistanceCollector<D> collector = newDistanceCollector();
    collector.set(distanceLimit);
    if (!target.updateBestDistance(cell, collector)) {
      return;
    }
    // Adjust the distance used to enqueue the cell by so it is a lower or upper bound on the
    // true distance to the cell, for minimum or maximum distance queries respectively.
    D distance =
        useConservativeCellDistance
            ? errorBoundedDistance(collector.distance())
            : collector.distance();
    queue.add(new QueueEntry<D>(distance, id, indexCell));
  }
}
