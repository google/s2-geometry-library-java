/*
 * Copyright 2015 Google Inc.
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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2PointIndex.Entry;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

/**
 * Given a set of points stored in an S2PointIndex, S2ClosestPointQuery provides methods that find
 * the closest point(s) to a given query point.
 *
 * <p>Example usage:
 *
 * <pre>{@code
 * void test(List<S2Point> points, List<S2Point> targets) {
 *   // The template argument allows auxiliary data to be attached to each point (in this case, the
 *   // array index).
 *   S2PointIndex<Integer> index = new S2PointIndex<>();
 *   for (int i = 0; i < points.size(); i++) {
 *     index.add(points.get(i), i);
 *   }
 *   S2ClosestPointQuery<Integer> query = new S2ClosestPointQuery<>(index);
 *   query.setMaxPoints(15);
 *   for (S2Point target : targets) {
 *     for (Result<Integer> result : query.findClosestPoints(target)) {
 *       // result.entry().point() is one of the found closest points.
 *       // result.entry().data() is the auxiliary data (the "points" array index).
 *       // result.distance() is the distance to the target point.
 *       doSomething(target, result.entry().point(), result.entry().data(), result.distance());
 *     }
 *   }
 * }
 * }</pre>
 *
 * <p>You can find either the k closest points, or all points within a given radius, or both (i.e.,
 * the k closest points up to a given maximum radius). E.g. to find all the points within 5
 * kilometers, call {@code query.setMaxDistance(S2Earth.toAngle(new Kilometers(5)));}.
 *
 * <p>You can also restrict the results to an arbitrary S2Region via {@link #setRegion(S2Region)}.
 *
 * <p>The implementation is designed to be very fast for both small and large point sets.
 *
 * <p>This class is not thread-safe, since setters and find* methods may mutate local state. This
 * class can however share an existing index with other queries, and is cheap to build so each
 * thread may simply create its own.
 */
public final class S2ClosestPointQuery<T> {
  // TODO(user): retune the constants.

  /** The maximum number of points to process by brute force. */
  private static final int MAX_BRUTE_FORCE_POINTS = 150;

  /** The maximum number of points to process without subdividing further. */
  private static final int MAX_LEAF_POINTS = 12;

  /** The index being queried. */
  private final S2PointIndex<T> index;

  /** The max number of closest points to find. */
  private int maxPoints;

  /** The max distance to search for points. */
  private S1Angle maxDistance;

  /** The region to restrict closest point search to. */
  private S2Region region;

  /** Whether to use brute force, which is cheaper when the index has few edges. */
  private boolean useBruteForce;

  /** A small (<6) cell covering of the indexed points. */
  private final List<S2CellId> indexCovering = Lists.newArrayList();

  /** Unprocessed cells for the current query being processed. */
  private final PriorityQueue<QueueEntry> queue = new PriorityQueue<>();

  /** The iterator for the last-known state of the index. New instance built by {@link #reset()}. */
  private S2Iterator<Entry<T>> iter;

  /** The covering of {@link #indexCovering}. Type is ArrayList due to {@link S2RegionCoverer}. */
  private final ArrayList<S2CellId> regionCovering = Lists.newArrayList();

  /** The covering of {@link #maxDistance}. Type is ArrayList due to {@link S2RegionCoverer}. */
  private final ArrayList<S2CellId> maxDistanceCovering = Lists.newArrayList();

  /** The intersection between the index and {@link #regionCovering}. */
  private final List<S2CellId> intersectionWithRegion = Lists.newArrayList();

  /** The intersection between the index and {@link #maxDistance}. */
  private final List<S2CellId> intersectionWithMaxDistance = Lists.newArrayList();

  /** Temporary storage for index entries that are of interest during query processing. */
  @SuppressWarnings({"rawtypes", "unchecked"})
  private final Entry<T>[] tmpPoints = new Entry[MAX_LEAF_POINTS];

  /** Temporary queue of results sorted in descending order. */
  private final PriorityQueue<Result<T>> results = new PriorityQueue<>();

  /**
   * Temporary distance to continue searching during a query, generally the distance of the furthest
   * point in the results found so far. Beyond this distance, we can safely ignore further candidate
   * points. Candidates that are exactly at the limit are ignored; this makes things easier in the
   * case of S2ClosestEdgeQuery and should not affect clients since distance measurements have a
   * small amount of error anyway.
   *
   * <p>Initially this is the same as the maximum distance specified by the user, but it can also be
   * updated by the algorithm (see maybeAddResult).
   */
  private S1ChordAngle maxDistanceLimit;

  /**
   * Construct a new query for the given index. Must call reset() before using the query, if the
   * index has been modified since the query was constructed.
   */
  public S2ClosestPointQuery(S2PointIndex<T> index) {
    this.index = index;
    maxPoints = Integer.MAX_VALUE;
    maxDistance = S1Angle.INFINITY;
    region = null;
    reset();
  }

  /** Resets the query state. This method must be called after modifying the underlying index. */
  public void reset() {
    iter = index.iterator();
    useBruteForce(index.numPoints() <= MAX_BRUTE_FORCE_POINTS);
  }

  /** Returns the underlying S2PointIndex. */
  public S2PointIndex<T> index() {
    return index;
  }

  /** Returns the max number of closest points to find. */
  public int getMaxPoints() {
    return maxPoints;
  }

  /** Sets a new max number of closest points to find. */
  public void setMaxPoints(int maxPoints) {
    Preconditions.checkArgument(maxPoints >= 1, "Must be at least 1.");
    this.maxPoints = maxPoints;
  }

  /** Returns the max distance between returned points and the given target. Default is +inf. */
  public S1Angle getMaxDistance() {
    return maxDistance;
  }

  /** Sets a new max distance to search for points. */
  public void setMaxDistance(S1Angle maxDistance) {
    this.maxDistance = maxDistance;
  }

  /** Returns the region in which point searches will be done. */
  public S2Region getRegion() {
    return region;
  }

  /*
   * Sets the region in which point searches will be done, or clears the region if {@code region} is
   * null.
   *
   * <p>Note that if you want to set the region to a disc around the target point, it is faster to
   * use setMaxDistance() instead. You can also call both methods, e.g. if you want to limit the
   * maximum distance to the target and also require that points lie within a given rectangle.
   */
  public void setRegion(S2Region region) {
    this.region = region;
  }

  /**
   * Sets whether distances are computed using "brute force" (i.e., by examining every point) rather
   * than using the S2PointIndex.
   *
   * <p>This is intended only for testing, benchmarking, and debugging.
   *
   * <p>Do not call before reset().
   */
  @VisibleForTesting
  public void useBruteForce(boolean useBruteForce) {
    this.useBruteForce = useBruteForce;
    if (!useBruteForce) {
      initIndexCovering();
    }
  }

  /**
   * Creates an empty list if 'list' is null, otherwise uses 'list'. Polls all results out of
   * {@link #results} and appends them to 'list' in reverse order. Returns the resulting list.
   */
  @CanIgnoreReturnValue
  private List<Result<T>> toList(List<Result<T>> list) {
    int size = results.size();
    int index = size;
    if (list == null) {
      list = Lists.newArrayListWithCapacity(size);
    } else {
      index += list.size();
    }
    // Allocate 'size' elements at the end of the list, and fill the items in reverse.
    list.addAll(Collections.<Result<T>>nCopies(size, null));
    while (size-- > 0) {
      list.set(--index, results.poll());
    }
    return list;
  }

  /**
   * Returns the closest points to {@code target} that satisfy the {@link #getMaxDistance()}, {@link
   * #getMaxPoints()}, and {@link #getRegion()} criteria, ordered by increasing distance. If there
   * are no criteria set, then all points are returned.
   *
   * <p>This class, including this method, is not thread-safe.
   */
  public List<Result<T>> findClosestPoints(S2Point target) {
    findClosestPointsToTarget(new PointTarget(target));
    return toList(null);
  }

  /**
   * As {@link #findClosestPoints(S2Point)}, but sorts the results and adds them at the end of the
   * given list.
   *
   * <p>This class, including this method, is not thread-safe.
   */
  public void findClosestPoints(List<Result<T>> results, S2Point target) {
    findClosestPointsToTarget(new PointTarget(target));
    toList(results);
  }

  /**
   * Convenience method that returns the closest point to the given target point, or null if no
   * points satisfy the {@link #getMaxDistance()} and {@link #getRegion()} criteria.
   *
   * <p>This class, including this method, is not thread-safe.
   */
  public Result<T> findClosestPoint(S2Point target) {
    setMaxPoints(1);
    return Iterables.getOnlyElement(findClosestPoints(target), null);
  }

  /**
   * Returns the closest points to the given edge AB. Otherwise similar to {@link
   * #findClosestPoints(S2Point)}.
   *
   * <p>This class, including this method, is not thread-safe.
   */
  public List<Result<T>> findClosestPointsToEdge(S2Point a, S2Point b) {
    findClosestPointsToTarget(new EdgeTarget(a, b));
    return toList(null);
  }

  /**
   * As {@link #findClosestPointsToEdge(S2Point, S2Point)}, but adds results to the given list.
   *
   * <p>This class, including this method, is not thread-safe.
   */
  public void findClosestPointsToEdge(List<Result<T>> results, S2Point a, S2Point b) {
    findClosestPointsToTarget(new EdgeTarget(a, b));
    toList(results);
  }

  /** A kind of query target. */
  private interface Target {
    /** Returns the approximate center of the target. */
    S2Point center();
    /** Returns the distance between this target and the given cell. */
    S1ChordAngle getDistance(S2Cell cell);
    /** Returns the radian radius of an angular cap that encloses this target. */
    double radius();
    /** Returns the smaller of {@code distance} and a new distance from target to {@code point}. */
    S1ChordAngle getMinDistance(S2Point point, S1ChordAngle distance);
  }

  /** A point query, used to find the closest points to a query point. */
  private static class PointTarget implements Target {
    private final S2Point point;

    public PointTarget(S2Point point) {
      this.point = point;
    }

    @Override
    public S2Point center() {
      return point;
    }

    @Override
    public double radius() {
      return 0;
    }

    @Override
    public S1ChordAngle getMinDistance(S2Point x, S1ChordAngle minDist) {
      S1ChordAngle angle = new S1ChordAngle(x, point);
      // See comment regarding ">=" in the findClosestPoints() main loop.
      return angle.compareTo(minDist) > 0 ? minDist : angle;
    }

    @Override
    public S1ChordAngle getDistance(S2Cell cell) {
      return cell.getDistance(point);
    }
  }

  /** An edge query, used to find the closest points to a query edge. */
  private static class EdgeTarget implements Target {
    private final S2Point a;
    private final S2Point b;

    public EdgeTarget(S2Point a, S2Point b) {
      this.a = a;
      this.b = b;
    }

    @Override
    public S2Point center() {
      return a.add(b).normalize();
    }

    @Override
    public double radius() {
      return 0.5 * a.angle(b);
    }

    @Override
    public S1ChordAngle getMinDistance(S2Point x, S1ChordAngle minDist) {
      return S2EdgeUtil.updateMinDistance(x, a, b, minDist);
    }

    @Override
    public S1ChordAngle getDistance(S2Cell cell) {
      return cell.getDistanceToEdge(a, b);
    }
  }

  /**
   * Computes the "index covering", which is a small number of S2CellIds that cover the indexed
   * points.
   */
  private void initIndexCovering() {
    //  There are two cases:
    // - If the index spans more than one face, then there is one covering cell per spanned face,
    //   just big enough to cover the index cells on that face.
    // - If the index spans only one face, then we find the smallest cell "C" that covers the index
    //   cells on that face (just like the case above). Then for each of the 4 children of "C", if
    //   the child contains any index cells then we create a covering cell that is big enough to
    //   just fit those index cells (i.e., shrinking the child as much as possible to fit its
    //   contents). This essentially replicates what would happen if we started with "C" as the
    //   covering cell, since "C" would immediately be split, except that we take the time to prune
    //   the children further since this will save work on every subsequent query.
    indexCovering.clear();
    iter.restart();
    if (iter.done()) {
      // Empty index.
      return;
    }

    S2Iterator<Entry<T>> nextIt = iter.copy();
    S2CellId indexNext = nextIt.id();
    S2Iterator<Entry<T>> lastIt = iter.copy();
    lastIt.finish();
    lastIt.prev();
    S2CellId indexLast = lastIt.id();
    if (!nextIt.equalIterators(lastIt)) {
      // The index has at least two cells. Choose a level such that the entire index can be spanned
      // with at most 6 cells (if the index spans multiple faces) or 4 cells (if the index spans a
      // single face).
      int level = indexNext.getCommonAncestorLevel(indexLast) + 1;

      // Visit each potential covering cell except the last (handled below).
      S2CellId coverLast = indexLast.parent(level);
      for (S2CellId cover = indexNext.parent(level);
          !cover.equals(coverLast) && !nextIt.done();
          cover = cover.next()) {
        // Skip any covering cells that don't contain any index cells.
        S2CellId coverMax = cover.rangeMax();
        if (nextIt.compareTo(coverMax) <= 0) {
          // Find the range of index cells contained by this covering cell and then shrink the cell
          // if necessary so that it just covers them.
          S2CellId prevId = indexNext;
          nextIt.seek(coverMax.next());
          indexNext = nextIt.id();
          S2Iterator<Entry<T>> cellLast = nextIt.copy();
          cellLast.prev();
          coverRange(prevId, cellLast.id());
        }
      }
    }
    coverRange(indexNext, indexLast);
  }

  /**
   * Adds a cell to indexCovering that covers the given inclusive range. This is done with {@link
   * S2CellId#getCommonAncestorLevel(S2CellId)}, which requires the cells have a common ancestor.
   */
  private void coverRange(S2CellId firstId, S2CellId lastId) {
    int level = firstId.getCommonAncestorLevel(lastId);
    indexCovering.add(firstId.parent(level));
  }

  private void findClosestPointsToTarget(Target target) {
    maxDistanceLimit = S1ChordAngle.fromS1Angle(maxDistance);
    if (useBruteForce) {
      findClosestPointsBruteForce(target);
    } else {
      findClosestPointsOptimized(target);
    }
  }

  private void findClosestPointsBruteForce(Target target) {
    for (iter.restart(); !iter.done(); iter.next()) {
      maybeAddResult(iter.entry(), target);
    }
  }

  private void findClosestPointsOptimized(Target target) {
    initQueue(target);
    while (!queue.isEmpty()) {
      QueueEntry entry = queue.poll();
      if (entry.distance().compareTo(maxDistanceLimit) >= 0) {
        queue.clear();
        break;
      }
      S2CellId child = entry.id.childBegin();
      // We already know that it has too many points, so process its children. Each child may either
      // be processed directly or enqueued again. The loop is optimized so that we don't seek
      // unnecessarily.
      boolean seek = true;
      for (int i = 0; i < 4; i++, child = child.next()) {
        seek = addCell(child, iter, seek, target);
      }
    }
  }

  private void maybeAddResult(Entry<T> entry, Target target) {
    S1ChordAngle distance = target.getMinDistance(entry.point(), maxDistanceLimit);
    if (distance == maxDistanceLimit) {
      // The previous 'max' reference is returned in this case, so only check object identity.
      return;
    }
    if (region != null && !region.contains(entry.point())) {
      return;
    }

    // Add this point to results.
    if (results.size() >= maxPoints) {
      // Replace the furthest result point.
      results.poll();
    }
    results.add(new Result<>(distance, entry));
    if (results.size() >= maxPoints) {
      maxDistanceLimit = results.peek().distance();
    }
  }

  private void initQueue(Target target) {
    // assert queue.isEmpty();

    // Optimization: rather than starting with the entire index, see if we can limit the search
    // region to a small disc. Then we can find a covering for that disc and intersect it with the
    // covering for the index. This can save a lot of work when the search region is small.

    if (maxPoints == 1) {
      // If the user is searching for just the closest point, we can compute an upper bound on
      // search radius by seeking to the target point in the index and looking at the adjacent index
      // points (in S2CellId order). The minimum distance to either of these points is an upper
      // bound on the search radius.
      //
      // TODO(user): The same strategy would also work for small values of maxPoints() > 1,
      // e.g. maxPoints() == 20, except that we would need to examine more neighbors (at least 20,
      // and preferably 20 in each direction). It's not clear whether this is a common case,
      // though, and also this would require extending maybeAddResult() so that it can remove
      // duplicate entries. (The points added here may be re-added by addCell(), but this is okay
      // when maxPoints() == 1.)
      iter.seek(S2CellId.fromPoint(target.center()));
      if (!iter.done()) {
        maybeAddResult(iter.entry(), target);
      }
      if (!iter.atBegin()) {
        iter.prev();
        maybeAddResult(iter.entry(), target);
      }
    }

    // We start with a covering of the set of indexed points, then intersect it with the given
    // region (if any) and maximum search radius disc (if any).
    List<S2CellId> initialCells = indexCovering;
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(4).build();
    if (region != null) {
      coverer.getCovering(region, regionCovering);
      S2CellUnion.getIntersection(indexCovering, regionCovering, intersectionWithRegion);
      initialCells = intersectionWithRegion;
    }
    if (!maxDistanceLimit.isInfinity()) {
      S2Cap searchCap =
          S2Cap.fromAxisAngle(
              target.center(),
              S1Angle.radians(target.radius() + maxDistanceLimit.toAngle().radians()));
      coverer.getFastCovering(searchCap, maxDistanceCovering);
      S2CellUnion.getIntersection(initialCells, maxDistanceCovering, intersectionWithMaxDistance);
      initialCells = intersectionWithMaxDistance;
    }
    iter.restart();
    for (int i = 0; i < initialCells.size() && !iter.done(); i++) {
      S2CellId id = initialCells.get(i);
      boolean seek = iter.compareTo(id.rangeMin()) <= 0;
      addCell(id, iter, seek, target);
    }
  }

  /**
   * Processes the cell at {@code id}, adding the contents of the cell immediately, or if there are
   * too many points, adding it to the queue to be subdivided. If {@code seek} is false, then
   * {@code iter} must already be positioned at the first indexed point within this cell.
   *
   * @return true if the cell was added to the queue, and false if it was processed immediately (in
   *     which case {@code iter} is left positioned at the next cell in S2CellId order.
   */
  @CanIgnoreReturnValue
  private boolean addCell(S2CellId id, S2Iterator<Entry<T>> iter, boolean seek, Target target) {
    if (seek) {
      iter.seek(id.rangeMin());
    }
    if (id.isLeaf()) {
      // Leaf cells can't be subdivided.
      for (; !iter.done() && iter.compareTo(id) == 0; iter.next()) {
        maybeAddResult(iter.entry(), target);
      }
      // No need to seek to next child.
      return false;
    }
    S2CellId last = id.rangeMax();
    int numPoints = 0;
    for (; !iter.done() && iter.compareTo(last) <= 0; iter.next()) {
      if (numPoints == MAX_LEAF_POINTS) {
        // This child cell has too many points, so enqueue it.
        S2Cell cell = new S2Cell(id);
        S1ChordAngle distance = target.getDistance(cell);
        if (distance.compareTo(maxDistanceLimit) < 0) {
          // We delay checking "region" as long as possible because it may be relatively expensive.
          if (region == null || region.mayIntersect(cell)) {
            queue.add(new QueueEntry(distance, id));
          }
        }
        // Seek to next child.
        return true;
      }
      tmpPoints[numPoints++] = iter.entry();
    }
    // There were few enough points that we might as well process them now.
    for (int i = 0; i < numPoints; i++) {
      maybeAddResult(tmpPoints[i], target);
    }
    // No need to seek to next child.
    return false;
  }

  /** A type that is comparable on distance only. */
  public abstract static class ChordComparable implements Comparable<ChordComparable> {
    protected final S1ChordAngle distance;

    ChordComparable(S1ChordAngle distance) {
      this.distance = distance;
    }

    public final S1ChordAngle distance() {
      return distance;
    }
  }

  /**
   * A query result paired with the distance to the query target. Natural order is by distance in
   * descending order.
   */
  public static class Result<T> extends ChordComparable {
    private final Entry<T> pointData;

    private Result(S1ChordAngle distance, Entry<T> pointData) {
      super(distance);
      this.pointData = pointData;
    }

    public Entry<T> entry() {
      return pointData;
    }

    @Override
    public int hashCode() {
      return pointData.hashCode();
    }

    @Override
    public boolean equals(Object o) {
      if (o instanceof Result) {
        // Don't need to test distance, it's derived from point data.
        Result<?> other = (Result<?>) o;
        return pointData.equals(other.pointData);
      } else {
        return false;
      }
    }

    @Override
    public String toString() {
      return distance().toAngle().degrees() + ": " + pointData;
    }

    @Override
    public int compareTo(ChordComparable other) {
      // The algorithm works by replacing the result whose distance is largest when a better
      // candidate is found, so we keep the entries sorted such that the largest distance is at the
      // top of the heap.
      return other.distance.compareTo(distance);
    }
  }

  /**
   * A queued cell waiting to be processed by the current query, ordered by distance to any point in
   * the cell in ascending order.
   */
  private static class QueueEntry extends ChordComparable {
    private final S2CellId id;

    QueueEntry(S1ChordAngle distance, S2CellId id) {
      super(distance);
      this.id = id;
    }

    @Override
    public int hashCode() {
      return id.hashCode() * 31 + distance().hashCode();
    }

    @Override
    public boolean equals(Object o) {
      if (o instanceof QueueEntry) {
        QueueEntry q = (QueueEntry) o;
        return id.equals(q.id) && distance.equals(q.distance);
      } else {
        return false;
      }
    }

    @Override
    public int compareTo(ChordComparable other) {
      // The algorithm works by replacing the result whose distance is largest when a better
      // candidate is found, so we keep the entries sorted such that the largest distance is at the
      // top of the heap.
      return distance.compareTo(other.distance);
    }
  }
}
