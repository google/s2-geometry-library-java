// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package com.google.common.geometry;

import static com.google.common.geometry.S2EdgeUtil.getUpdateMinDistanceMaxError;
import static java.util.Collections.reverseOrder;

import com.google.common.geometry.S2CellIndex.ContentsIterator;
import com.google.common.geometry.S2CellIndex.NonEmptyRangeIterator;
import com.google.common.geometry.S2CellIndex.RangeIterator;
import com.google.common.geometry.S2ShapeUtil.PointVisitor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import org.jspecify.annotations.Nullable;

/**
 * S2ClosestCellQuery is a helper class for finding the closest cell(s) to a given point, edge,
 * S2Cell, S2CellUnion, or geometry collection. A typical use case would be to add a collection of
 * S2Cell coverings to an S2CellIndex (representing a collection of original geometry), and then use
 * S2ClosestCellQuery to find all coverings that are within a given distance of some target geometry
 * (which could be represented exactly, or could also be a covering). The distance to the original
 * geometry corresponding to each covering could then be measured more precisely if desired.
 *
 * <p>For example, here is how to find all cells that are closer than "distanceLimit" to a given
 * target point:
 * <pre>
 *   S2ClosestCellQuery query = new S2ClosestCellQuery(cellIndex);
 *   query.setMaxDistance(distanceLimit);
 *   S2ClosestCellQuery.PointTarget target = new S2ClosestCellQuery.PointTarget(targetPoint);
 *   for (S2ClosestCellQuery.Result result : query.findClosestCells(target)) {
 *     // result.distance() is the distance to the target.
 *     // result.cellId() is the indexed S2CellId.
 *     // result.label() is the integer label associated with the S2CellId.
 *     doSomething(targetPoint, result);
 *   }
 * </pre>
 *
 * <p>You can find either the k closest cells, or all cells within a given radius, or both (i.e.,
 * the k closest cells up to a given maximum radius). By default *all* cells are returned, so you
 * should always specify either maxResults() or maxDistance() or both. You can also restrict the
 * results to cells that intersect a given S2Region; for example:
 * <pre>
 *   S2LatLngRect rect = ...;
 *   query.setRegion(rect);
 * </pre>
 *
 * <p>There is a findClosestCell() convenience method that returns the closest cell. However, if you
 * only need to test whether the distance is above or below a given threshold (e.g., 10 km), it is
 * typically much faster to use the isDistanceLess() method instead. Unlike findClosestCell(), this
 * method stops as soon as it can prove that the minimum distance is either above or below the
 * threshold. Example usage:
 * <pre>
 *   if (query.isDistanceLess(target, limitDistance)) ...
 * </pre>
 *
 * <p>To find the closest cells to a query edge rather than a point, use:
 * <pre>
 *   S2ClosestCellQuery.EdgeTarget target = new S2ClosestCellQuery.EdgeTarget(v0, v1);
 *   query.findClosestCells(target);
 * </pre>
 *
 * <p>Similarly you can find the closest cells to an S2Cell using an S2ClosestCellQuery.CellTarget,
 * you can find the closest cells to an S2CellUnion using an S2ClosestCellQuery.CellUnionTarget, and
 * you can find the closest cells to an arbitrary collection of points, polylines, and polygons by
 * using an S2ClosestCellQuery.ShapeIndexTarget.
 *
 * <p>The implementation is designed to be fast for both simple and complex geometric objects.
 */
@CheckReturnValue
public class S2ClosestCellQuery {

  /**
   * The minimum number of ranges that a cell must contain to enqueue it rather than processing its
   * contents immediately.
   */
  private static final int MIN_RANGES_TO_ENQUEUE = 6;

  private final S2CellIndex index;
  private int maxResults = Integer.MAX_VALUE;
  private S1ChordAngle maxDistance = S1ChordAngle.INFINITY;
  private S1ChordAngle maxError = S1ChordAngle.ZERO;
  private @Nullable S2Region region = null;
  private boolean useBruteForce = false;

  /**
   * The distance beyond which we can safely ignore further candidate cells. (Candidates that are
   * exactly at the limit are ignored; this is more efficient for updateBestDistance() and should
   * not affect clients since distance measurements have a small amount of error anyway.)
   *
   * <p>Initially this is the same as the maximum distance specified by the user, but it can also be
   * updated by the algorithm (see maybeAddResult).
   */
  private S1ChordAngle distanceLimit;

  /**
   * True if maxError() must be subtracted from priority queue cell distances in order to ensure
   * that such distances are measured conservatively. This is true only if the target takes
   * advantage of maxError() in order to return faster results, and
   * {@code 0 < maxError() < distanceLimit}.
   */
  private boolean useConservativeCellDistance;

  private Target target;

  /**
   * For the optimized algorithm we precompute the top-level S2CellIds that will be added to the
   * priority queue. There can be at most 6 of these cells. Essentially this is just a covering of
   * the indexed cells.
   */
  private final List<S2CellId> indexCovering = new ArrayList<>();

  /** Temporary, defined here to avoid multiple allocations / initializations. */
  private final List<S2CellId> maxDistanceCovering = new ArrayList<>();

  /** Temporary, defined here to avoid multiple allocations / initializations. */
  private final List<S2CellId> intersectionWithMaxDistance = new ArrayList<>();

  /**
   * The current result set is stored in one of three ways:
   *
   * <ul>
   * <li>If maxResults == 1, the best result is kept in resultSingleton.
   * <li>If maxResults == Integer.MAX_VALUE, results are appended to resultVector and sorted/uniqued
   *     at the end.
   * <li>Otherwise results are kept in a PriorityQueue so that we can progressively reduce the
   *     distance limit once maxResults results have been found.
   * </ul>
   */
  private Result resultSingleton = new Result(S1ChordAngle.INFINITY, S2CellId.none(), -1);
  private final List<Result> resultVector = new ArrayList<>();
  private final PriorityQueue<Result> resultSet = new PriorityQueue<>(8, reverseOrder());

  /**
   * When the results are stored in a PriorityQueue, duplicates can usually be removed simply by
   * inserting candidate cells into the current result set. However this is not true if
   * {@code maxError() > 0} and the Target subtype takes advantage of this by returning suboptimal
   * distances. This is because when updateBestDistance() is called with different "minDist"
   * parameters, the implementation may return a different distance for the same cell. Since the
   * result queue is keyed by distance this can create duplicate results.
   *
   * <p>The flag below is true when duplicates must be avoided explicitly. This is achieved by
   * maintaining a separate set keyed by (cellId, label) only, and checking whether each edge is in
   * that set before computing the distance to it.
   */
  private boolean avoidDuplicates;
  private final HashSet<LabelledCell> testedCellsSet = new HashSet<>();

  /**
   * The algorithm maintains a priority queue of unprocessed S2CellIds, sorted in increasing order
   * of distance from the target.
   */
  private final PriorityQueue<QueueEntry> queue = new PriorityQueue<>();

  /**
   * Used to iterate over the contents of an S2CellIndex range. It is defined here to take advantage
   * of the fact that when multiple ranges are visited in increasing order, duplicates can
   * automatically be eliminated.
   */
  private final ContentsIterator contentsIt;

  /**
   * Initializes the query.
   *
   * <p>REQUIRES: "index" must persist for the lifetime of this object. REQUIRES: reset() must be
   * called if "index" is modified.
   */
  public S2ClosestCellQuery(S2CellIndex index) {
    this.index = index;
    this.contentsIt = index.contents();
    reset();
  }

  /**
   * Reinitializes the query. This method must be called if the underlying S2CellIndex is modified.
   */
  public void reset() {
    indexCovering.clear();
  }

  /** Returns the underlying S2CellIndex. */
  public S2CellIndex index() {
    return index;
  }

  /** Specifies that at most "maxResults" cells should be returned. */
  public void setMaxResults(int maxResults) {
    this.maxResults = maxResults;
  }

  /** Returns the maximum number of results. */
  public int getMaxResults() {
    return maxResults;
  }

  /**
   * Specifies that only cells whose distance to the target is less than "maxDistance" should be
   * returned.
   *
   * <p>Note that cells whose distance is exactly equal to "maxDistance" are not returned. Normally
   * this doesn't matter, because distances are not computed exactly in the first place, but if such
   * cells are needed then see setInclusiveMaxDistance() below.
   */
  public void setMaxDistance(S1ChordAngle maxDistance) {
    this.maxDistance = maxDistance;
  }

  /**
   * Versions of setMaxDistance that take an S1Angle argument. (Note that these functions require a
   * conversion, and that the S1ChordAngle versions are preferred.)
   */
  public void setMaxDistance(S1Angle maxDistance) {
    this.maxDistance = S1ChordAngle.fromS1Angle(maxDistance);
  }

  /**
   * Like setMaxDistance(), except that cells whose distance is exactly equal to "maxDistance" are
   * also returned. Equivalent to calling setMaxDistance(maxDistance.successor()).
   */
  public void setInclusiveMaxDistance(S1ChordAngle maxDistance) {
    setMaxDistance(maxDistance.successor());
  }

  /** Version of setInclusiveMaxDistance that takes an S1Angle argument. */
  public void setInclusiveMaxDistance(S1Angle maxDistance) {
    setInclusiveMaxDistance(S1ChordAngle.fromS1Angle(maxDistance));
  }

  /**
   * Like setInclusiveMaxDistance(), except that "maxDistance" is also increased by the maximum
   * error in the distance calculation. This ensures that all cells whose true distance is less than
   * or equal to "maxDistance" will be returned (along with some cells whose true distance is
   * slightly greater).
   *
   * <p>Algorithms that need to do exact distance comparisons can use this option to find a set of
   * candidate cells that can then be filtered further.
   */
  public void setConservativeMaxDistance(S1ChordAngle maxDistance) {
    setMaxDistance(maxDistance.plusError(getUpdateMinDistanceMaxError(maxDistance)).successor());
  }

  /** Version of setConservativeMaxDistance that takes an S1Angle argument. */
  public void setConservativeMaxDistance(S1Angle maxDistance) {
    setConservativeMaxDistance(S1ChordAngle.fromS1Angle(maxDistance));
  }

  /**
   * Specifies that cells up to maxError() further away than the true closest cells may be
   * substituted in the result set, as long as such cells satisfy all the remaining search criteria
   * (such as maxDistance). This option only has an effect if maxResults() is also specified;
   * otherwise all cells closer than maxDistance() will always be returned.
   *
   * <p>Note that this does not affect how the distance between cells is computed; it simply gives
   * the algorithm permission to stop the search early as soon as the best possible improvement
   * drops below maxError().
   *
   * <p>This can be used to implement distance predicates efficiently. For example, to determine
   * whether the minimum distance is less than D, the isDistanceLess() method sets maxResults() == 1
   * and maxDistance() == maxError() == D. This causes the algorithm to terminate as soon as it
   * finds any cell whose distance is less than D, rather than continuing to search for a cell that
   * is even closer.
   */
  public void setMaxError(S1ChordAngle maxError) {
    this.maxError = maxError;
  }

  /** Version of setMaxError that takes an S1Angle argument. */
  public void setMaxError(S1Angle maxError) {
    this.maxError = S1ChordAngle.fromS1Angle(maxError);
  }

  /** Returns the maximum error in the distance calculation. */
  public S1ChordAngle getMaxError() {
    return maxError;
  }

  /**
   * Specifies that cells must intersect the given S2Region. "region" is owned by the caller and
   * must persist during the lifetime of this object. The value may be changed between calls to
   * findClosestCells(), or reset by calling setRegion(null).
   *
   * <p>Note that if you want to set the region to a disc around a target point, it is faster to use
   * a PointTarget with setMaxDistance() instead. You can also call both methods, e.g. to set a
   * maximum distance and also require that cells lie within a given rectangle.
   */
  public void setRegion(@Nullable S2Region region) {
    this.region = region;
  }

  /** Returns the region in which to restrict the results. */
  public @Nullable S2Region getRegion() {
    return region;
  }

  /**
   * Specifies that distances should be computed by examining every cell rather than using the
   * S2ShapeIndex. This is useful for testing, benchmarking, and debugging.
   */
  public void useBruteForce(boolean useBruteForce) {
    this.useBruteForce = useBruteForce;
  }

  /** Returns whether the brute force algorithm is used. */
  public boolean isUseBruteForce() {
    return useBruteForce;
  }

  /**
   * Returns the closest cells to the given target that satisfy the current options. This method may
   * be called multiple times.
   */
  public List<Result> findClosestCells(Target target) {
    List<Result> results = new ArrayList<>();
    findClosestCells(target, results);
    return results;
  }

  /**
   * This version can be more efficient when this method is called many times, since it does not
   * require allocating a new vector on each call.
   */
  public void findClosestCells(Target target, List<Result> results) {
    findClosestCellsInternal(target);
    results.clear();
    if (maxResults == 1) {
      if (!resultSingleton.isEmpty()) {
        results.add(resultSingleton);
      }
    } else if (maxResults == Integer.MAX_VALUE) {
      Collections.sort(resultVector);
      // Remove duplicates
      Result last = null;
      for (Result r : resultVector) {
        if (!r.equals(last)) {
          results.add(r);
          last = r;
        }
      }
      resultVector.clear();
    } else {
      results.addAll(resultSet);
      Collections.sort(results); // Sort from closest to furthest
      resultSet.clear();
    }
  }

  /**
   * Returns the closest cell to the target. If no cell satisfies the search criteria, then the
   * Result object will have distance == Infinity() and isEmpty() == true.
   */
  public Result findClosestCell(Target target) {
    int oldMaxResults = maxResults;
    maxResults = 1;
    findClosestCellsInternal(target);
    maxResults = oldMaxResults;
    return resultSingleton;
  }

  /**
   * Returns the minimum distance to the target. If the index or target is empty, returns
   * S1ChordAngle.INFINITY.
   *
   * <p>Use isDistanceLess() if you only want to compare the distance against a threshold value,
   * since it is often much faster.
   */
  public S1ChordAngle getDistance(Target target) {
    return findClosestCell(target).distance();
  }

  /**
   * Returns true if the distance to "target" is less than "limit".
   *
   * <p>This method is usually much faster than getDistance(), since it is much less work to
   * determine whether the minimum distance is above or below a threshold than it is to calculate
   * the actual minimum distance.
   */
  public boolean isDistanceLess(Target target, S1ChordAngle limit) {
    int oldMaxResults = maxResults;
    S1ChordAngle oldMaxDistance = maxDistance;
    S1ChordAngle oldMaxError = maxError;
    maxResults = 1;
    maxDistance = limit;
    maxError = S1ChordAngle.STRAIGHT;
    findClosestCellsInternal(target);
    maxResults = oldMaxResults;
    maxDistance = oldMaxDistance;
    maxError = oldMaxError;
    return !resultSingleton.isEmpty();
  }

  /**
   * Like isDistanceLess(), but also returns true if the distance to "target" is exactly equal to
   * "limit".
   */
  public boolean isDistanceLessOrEqual(Target target, S1ChordAngle limit) {
    return isDistanceLess(target, limit.successor());
  }

  /**
   * Like isDistanceLessOrEqual(), except that "limit" is increased by the maximum error in the
   * distance calculation. This ensures that this function returns true whenever the true, exact
   * distance is less than or equal to "limit".
   */
  public boolean isConservativeDistanceLessOrEqual(Target target, S1ChordAngle limit) {
    return isDistanceLess(target, limit.plusError(getUpdateMinDistanceMaxError(limit)).successor());
  }

  private void findClosestCellsInternal(Target target) {
    this.target = target;
    testedCellsSet.clear();
    distanceLimit = maxDistance;
    resultSingleton = new Result(S1ChordAngle.INFINITY, S2CellId.none(), -1);
    resultVector.clear();
    resultSet.clear();

    if (distanceLimit.equals(S1ChordAngle.ZERO)) {
      return;
    }

    boolean targetUsesMaxError =
        !maxError.equals(S1ChordAngle.ZERO) && target.setMaxError(maxError);
    useConservativeCellDistance =
        targetUsesMaxError
            && (distanceLimit.equals(S1ChordAngle.INFINITY)
                || distanceLimit.compareTo(maxError) > 0);

    if (useBruteForce || index.numCells() <= target.maxBruteForceIndexSize()) {
      avoidDuplicates = false;
      findClosestCellsBruteForce();
    } else {
      avoidDuplicates = targetUsesMaxError && maxResults > 1;
      findClosestCellsOptimized();
    }
  }

  private void findClosestCellsBruteForce() {
    contentsIt.clear();
    RangeIterator range = index.ranges();
    for (range.begin(); !range.done(); range.next()) {
      addRange(range);
    }
  }

  private void findClosestCellsOptimized() {
    initQueue();
    while (!queue.isEmpty()) {
      QueueEntry entry = queue.poll();
      if (entry.distance.compareTo(distanceLimit) >= 0) {
        queue.clear();
        break;
      }
      S2CellId child = entry.id.childBegin();
      boolean seek = true;
      NonEmptyRangeIterator range = index.nonEmptyRanges();
      for (int i = 0; i < 4; ++i, child = child.next()) {
        seek = processOrEnqueue(child, range, seek);
      }
    }
  }

  private void initQueue() {
    queue.clear();
    S2Cap cap = target.getCapBound();
    if (cap.isEmpty()) {
      return;
    }

    if (maxResults == 1) {
      NonEmptyRangeIterator range = index.nonEmptyRanges();
      S2CellId targetId = S2CellId.fromPoint(cap.axis());
      range.seek(targetId);
      addRange(range);
      if (distanceLimit.equals(S1ChordAngle.ZERO)) {
        return;
      }

      if (range.startId().greaterThan(targetId) && range.prev()) {
        addRange(range);
        if (distanceLimit.equals(S1ChordAngle.ZERO)) {
          return;
        }
      }
    }

    if (indexCovering.isEmpty()) {
      initCovering();
    }
    List<S2CellId> initialCells = indexCovering;
    if (distanceLimit.compareTo(S1ChordAngle.INFINITY) < 0) {
      S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(4).build();
      S1ChordAngle radius = S1ChordAngle.add(cap.radius(), distanceLimit);
      S2Cap searchCap = S2Cap.fromAxisChord(cap.axis(), radius);
      maxDistanceCovering.clear();
      coverer.getFastCovering(searchCap, (ArrayList<S2CellId>) maxDistanceCovering);
      intersectionWithMaxDistance.clear();
      S2CellUnion.getIntersection(initialCells, maxDistanceCovering, intersectionWithMaxDistance);
      initialCells = intersectionWithMaxDistance;
    }

    NonEmptyRangeIterator range = index.nonEmptyRanges();
    for (int i = 0; i < initialCells.size(); ++i) {
      S2CellId id = initialCells.get(i);
      boolean seek = (i == 0) || id.rangeMin().greaterOrEquals(range.limitId());
      var unused = processOrEnqueue(id, range, seek);
      if (range.done()) {
        break;
      }
    }
  }

  private void initCovering() {
    // Compute the "index covering", which is a small number of S2CellIds that cover the indexed
    // cells. There are two cases:
    //
    //  - If the index spans more than one face, then there is one covering cell per spanned face,
    //    just big enough to cover the indexed cells on that face.
    //
    //  - If the index spans only one face, then we find the smallest cell "C" that covers the
    //    indexed cells on that face (just like the case above). Then for each of the 4 children of
    //    "C", if the child contains any index cells then we create a covering cell that is big
    //    enough to just fit those indexed cells (i.e., shrinking the child as much as possible to
    //    fit its contents). This essentially replicates what would happen if we started with "C" as
    //    the covering cell, since "C" would immediately be split, except that we take the time to
    //    prune the children further since this will save work on every subsequent query.
    indexCovering.clear();
    NonEmptyRangeIterator it = index.nonEmptyRanges();
    NonEmptyRangeIterator last = index.nonEmptyRanges();
    it.begin();
    last.finish();
    if (!last.prev()) {
      return;
    }

    S2CellId indexLastId = last.limitId().prev();
    if (!it.startId().equals(last.startId())) {
      // The index contains at least two distinct S2CellIds (because otherwise there would only be
      // one non-empty range). Choose a level such that the entire index can be spanned with at most
      // 6 cells (if the index spans multiple faces) or 4 cells (it the index spans a single face).
      int level = it.startId().getCommonAncestorLevel(indexLastId) + 1;
      S2CellId startId = it.startId().parent(level);
      S2CellId lastId = indexLastId.parent(level);
      for (S2CellId id = startId; !id.equals(lastId); id = id.next()) {
        // Skip any covering cells that don't contain an indexed range.
        if (id.rangeMax().lessThan(it.startId())) {
          continue;
        }

        // Find the indexed range contained by this covering cell and then shrink the cell if
        // necessary so that it just covers this range.
        S2CellId cellFirstId = it.startId();
        it.seek(id.rangeMax().next());
        // Find the last leaf cell covered by the previous non-empty range.
        last.seek(it.startId());
        var unused = last.prev();
        addInitialRange(cellFirstId, last.limitId().prev());
      }
    }
    addInitialRange(it.startId(), indexLastId);
  }

  /** Adds a cell to indexCovering that covers the given inclusive range. */
  private void addInitialRange(S2CellId firstId, S2CellId lastId) {
    // Add the lowest common ancestor of the given range.
    int level = firstId.getCommonAncestorLevel(lastId);
    indexCovering.add(firstId.parent(level));
  }

  private void maybeAddResult(S2CellId cellId, int label) {
    if (avoidDuplicates && !testedCellsSet.add(new LabelledCell(cellId, label))) {
      return;
    }

    S2Cell cell = new S2Cell(cellId);
    DistanceCollector<S1ChordAngle> collector = S1ChordAngle.minCollector();
    collector.set(distanceLimit);

    if (!target.updateBestDistance(cell, collector)) {
      return;
    }

    if (region != null && !region.mayIntersect(cell)) {
      return;
    }

    Result result = new Result(collector.distance(), cellId, label);
    if (maxResults == 1) {
      resultSingleton = result;
      distanceLimit = S1ChordAngle.sub(result.distance(), maxError);
    } else if (maxResults == Integer.MAX_VALUE) {
      resultVector.add(result);
    } else {
      resultSet.add(result);
      if (resultSet.size() >= maxResults) {
        if (resultSet.size() > maxResults) {
          resultSet.poll();
        }
        distanceLimit = S1ChordAngle.sub(resultSet.peek().distance(), maxError);
      }
    }
  }

  /**
   * Either process the contents of the given cell immediately, or add it to the queue to be
   * subdivided. If "seek" is false, then "iter" must be positioned at the first non-empty range
   * (if any) with startId >= id.rangeMin().
   *
   * @return true if the cell was added to the queue, and false if it was processed immediately,
   *     in which case "iter" is positioned at the first non-empty range (if any) with startId >
   *     id.rangeMax().
   */
  @CanIgnoreReturnValue
  private boolean processOrEnqueue(S2CellId id, NonEmptyRangeIterator iter, boolean seek) {
    if (seek) {
      iter.seek(id.rangeMin());
    }
    S2CellId last = id.rangeMax();
    if (iter.startId().greaterThan(last)) {
      return false;
    }

    RangeIterator maxIt = index.ranges();
    maxIt.seek(iter.startId());
    if (maxIt.advance(MIN_RANGES_TO_ENQUEUE - 1) && maxIt.startId().lessOrEquals(last)) {
      S2Cell cell = new S2Cell(id);
      DistanceCollector<S1ChordAngle> collector = S1ChordAngle.minCollector();
      collector.set(distanceLimit);

      if (target.updateBestDistance(cell, collector)
          && (region == null || region.mayIntersect(cell))) {
        S1ChordAngle dist = collector.distance();
        if (useConservativeCellDistance) {
          dist = S1ChordAngle.sub(dist, maxError);
        }
        queue.add(new QueueEntry(dist, id));
      }
      return true;
    }
    for (; iter.startId().lessOrEquals(last); iter.next()) {
      addRange(iter);
      if (iter.done()) {
        break;
      }
    }
    return false;
  }

  private void addRange(RangeIterator range) {
    for (contentsIt.startUnion(range); !contentsIt.done(); contentsIt.next()) {
      maybeAddResult(contentsIt.cellId(), contentsIt.label());
    }
  }

  /** Each "Result" object represents a closest (s2cellId, label) pair. */
  public static class Result implements Comparable<Result> {
    private final S1ChordAngle distance;
    private final S2CellId cellId;
    private final int label;

    /** Constructs a Result object for the given (cellId, label) pair. */
    public Result(S1ChordAngle distance, S2CellId cellId, int label) {
      this.distance = distance;
      this.cellId = cellId;
      this.label = label;
    }

    /** The distance from the target to this cell. */
    public S1ChordAngle distance() {
      return distance;
    }

    /** The cell itself. */
    public S2CellId cellId() {
      return cellId;
    }

    /** The label associated with this S2CellId. */
    public int label() {
      return label;
    }

    /**
     * Returns true if this Result object does not refer to any cell. (The only case where an empty
     * Result is returned is when the findClosestCell() method does not find any cells that meet the
     * specified criteria.)
     */
    public boolean isEmpty() {
      return cellId.equals(S2CellId.none());
    }

    @Override public int compareTo(Result other) {
      int c = distance.compareTo(other.distance);
      if (c != 0) {
        return c;
      }
      c = cellId.compareTo(other.cellId);
      if (c != 0) {
        return c;
      }
      return Integer.compare(label, other.label);
    }

    @Override public boolean equals(Object o) {
      return o instanceof Result r
          && distance.equals(r.distance)
          && cellId.equals(r.cellId)
          && label == r.label;
    }

    @Override public int hashCode() {
      return cellId.hashCode() * 31 + label;
    }
  }

  /**
   * "Target" represents the geometry to which the distance is measured. There are subtypes for
   * measuring the distance to a point, an edge, an S2Cell, or an S2ShapeIndex (an arbitrary
   * collection of geometry).
   */
  public interface Target extends S2BestDistanceTarget<S1ChordAngle> {
    @Override int maxBruteForceIndexSize();
  }

  /** Target subtype that computes the closest distance to a point. */
  public static class PointTarget extends S2BestEdgesQueryBase.PointTarget<S1ChordAngle>
      implements Target {
    /** Constructs a PointTarget for the given point. */
    public PointTarget(S2Point p) {
      super(p);
    }
    @Override public int maxBruteForceIndexSize() {
      return 9;
    }
    @Override public S2Cap getCapBound() {
      return S2Cap.fromAxisChord(point, S1ChordAngle.ZERO);
    }
  }

  /** Target subtype that computes the closest distance to an edge. */
  public static class EdgeTarget extends S2BestEdgesQueryBase.EdgeTarget<S1ChordAngle>
      implements Target {
    /** Constructs an EdgeTarget for the given edge. */
    public EdgeTarget(S2Point a, S2Point b) {
      super(a, b);
    }
    @Override public int maxBruteForceIndexSize() {
      return 5;
    }
    @Override public S2Cap getCapBound() {
      double r2 = getHalfEdgeLength2();
      return S2Cap.fromAxisChord(a.add(b).normalize(), S1ChordAngle.fromLength2(r2));
    }
  }

  /**
   * Target subtype that computes the closest distance to an S2Cell (including the interior of the
   * cell).
   */
  public static class CellTarget extends S2BestEdgesQueryBase.CellTarget<S1ChordAngle>
      implements Target {
    /** Constructs a CellTarget for the given S2Cell. */
    public CellTarget(S2Cell cell) {
      super(cell);
    }
    @Override public int maxBruteForceIndexSize() {
      return 6;
    }
    @Override public S2Cap getCapBound() {
      return cell.getCapBound();
    }
  }

  /** Target subtype that computes the closest distance to an S2CellUnion. */
  public static class CellUnionTarget implements Target {
    private final S2CellUnion cellUnion;

    /** Constructs a CellUnionTarget for the given S2CellUnion. */
    public CellUnionTarget(S2CellUnion cellUnion) {
      this.cellUnion = cellUnion;
    }

    @Override public boolean updateBestDistance(
        S2Point p, DistanceCollector<S1ChordAngle> collector) {
      boolean updated = false;
      for (int i = 0; i < cellUnion.size(); i++) {
        if (collector.update(p, new S2Cell(cellUnion.cellId(i)))) {
          updated = true;
        }
      }
      return updated;
    }

    @Override public boolean updateBestDistance(
        S2Point v0, S2Point v1, DistanceCollector<S1ChordAngle> collector) {
      boolean updated = false;
      for (int i = 0; i < cellUnion.size(); i++) {
        if (collector.update(v0, v1, new S2Cell(cellUnion.cellId(i)))) {
          updated = true;
        }
      }
      return updated;
    }

    @Override public boolean updateBestDistance(
        S2Cell cell, DistanceCollector<S1ChordAngle> collector) {
      boolean updated = false;
      for (int i = 0; i < cellUnion.size(); i++) {
        if (collector.update(new S2Cell(cellUnion.cellId(i)), cell)) {
          updated = true;
        }
      }
      return updated;
    }

    @Override public boolean visitConnectedComponentPoints(PointVisitor visitor) {
      for (int i = 0; i < cellUnion.size(); i++) {
        if (!visitor.apply(cellUnion.cellId(i).toPoint())) {
          return false;
        }
      }
      return true;
    }

    @Override public S2Cap getCapBound() {
      return cellUnion.getCapBound();
    }
    @Override public int maxBruteForceIndexSize() {
      return 8;
    }
  }

  /**
   * Target subtype that computes the closest distance to an S2ShapeIndex (an arbitrary collection
   * of points, polylines, and/or polygons).
   *
   * <p>By default, distances are measured to the boundary and interior of polygons in the
   * S2ShapeIndex rather than to polygon boundaries only. If you wish to change this behavior, you
   * may call
   *
   * <pre>
   *   target.setIncludeInteriors(false);
   * </pre>
   */
  public static class ShapeIndexTarget extends S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle>
      implements Target {
    /** Constructs a ShapeIndexTarget for the given S2ShapeIndex. */
    public ShapeIndexTarget(S2ShapeIndex index) {
      super(index, S2ClosestEdgeQuery.builder());
    }

    @Override public int maxBruteForceIndexSize() {
      return 7;
    }
  }

  private static class LabelledCell {
    final S2CellId cellId;
    final int label;

    LabelledCell(S2CellId cellId, int label) {
      this.cellId = cellId;
      this.label = label;
    }

    @Override public int hashCode() {
      return cellId.hashCode() * 31 + label;
    }

    @Override public boolean equals(Object o) {
      return o instanceof LabelledCell c && label == c.label && cellId.equals(c.cellId);
    }
  }

  private static class QueueEntry implements Comparable<QueueEntry> {
    /** A lower bound on the distance from the target to "id". This is the key of the queue. */
    final S1ChordAngle distance;
    /** The cell being queued. */
    final S2CellId id;

    QueueEntry(S1ChordAngle distance, S2CellId id) {
      this.distance = distance;
      this.id = id;
    }

    @Override public int compareTo(QueueEntry other) {
      return distance.compareTo(other.distance);
    }
  }
}
