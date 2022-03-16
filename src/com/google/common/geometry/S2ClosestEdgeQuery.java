package com.google.common.geometry;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Objects;
import com.google.common.base.Preconditions;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import javax.annotation.concurrent.NotThreadSafe;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

/**
 * S2ClosestEdgeQuery is a helper class for searching within an S2ShapeIndex to find the closest edge(s)
 * to a given point, edge, S2Cell, or geometry collection.  For example, given a set of polylines,
 * the following code efficiently finds the closest 5 edges to a query point:
 *
 * <pre>
 * public void test(List<S2Polyline> polylines, S2Point point) {
 *     S2ShapeIndex index = new S2ShapeIndex();
 *     for (S2Polyline polyline : polylines) {
 *         index.add(polyline);
 *     }
 *     S2ClosestEdgeQuery.Options options = new S2ClosestEdgeQuery.Options().setMaxResults(5);
 *     S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index, options);
 *     S2ClosestEdgeQuery.PointTarget target = new S2ClosestEdgeQuery.PointTarget(point);
 *     S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
 *     for (Result result : query.findClosestEdges(target)) {
 *         // The Result object contains the following accessors:
 *         //   distance is the distance to the edge.
 *         //   shape the S2Shape containing the edge.
 *         //   edgeId identifies the edge with the given shape.
 *         //   isInterior indicates that the result is an interior point.
 *         //
 *         // The following convenience methods may also be useful:
 *         //   query.getEdge(result, edge) returns the endpoints of the edge.
 *         //   query.project(point, result) computes the closest point on the
 *         //       result edge to the given target point.
 *         S2Shape polyline = result.shape;
 *         int edgeIndex = result.edgeId;
 *         S1ChordAngle distance = result.distance;  // Can convert to S1Angle.
 *         query.GetEdge(result, edge);
 *         S2Point closestPoint = query.project(point, result);
 *         ...
 *     }
 * }
 * </pre>
 * <p>
 * You can find either the k closest edges, or all edges within a given radius, or both (i.e., the k
 * closest edges up to a given maximum radius).
 * <p>
 * By default *all* edges are returned, so you should always specify either setMaxResults() or setMaxDistance()
 * or both.  There is also a FindClosestEdge() convenience method that returns only the closest edge.
 * <p>
 * Note that by default, distances are measured to the boundary and interior of polygons. For example, if a point
 * is inside a polygon then its distance is zero. To change this behavior, call setIncludeInteriors(false).
 * <p>
 * If you only need to test whether the distance is above or below a given threshold (e.g., 10 km), you can use
 * the isDistanceLess() method.  This is much faster than actually calculating the distance with findClosestEdge(),
 * since the implementation can stop as soon as it can prove that the minimum distance is either above or below
 * the threshold.
 * <p>
 * To find the closest edges to a query edge rather than a point, use:
 *
 * <pre>
 *     S2ClosestEdgeQuery.EdgeTarget target = new S2ClosestEdgeQuery.EdgeTarget(v0, v1);
 *     query.findClosestEdges(target);
 * </pre>
 * <p>
 * The implementation is designed to be fast for both simple and complex geometric objects.
 */
@NotThreadSafe
public class S2ClosestEdgeQuery {

    /** Logger. */
    private static final Logger log = Platform.getLoggerForClass(S2ClosestEdgeQuery.class);

    /** The index being queried. */
    private final S2ShapeIndex index;

    /////////////// Options /////////////////

    /** The query options. */
    private Options options;

    /** Whether to use brute force, which is cheaper when the index has few edges. */
    private boolean useBruteForce;

    /** The query target. */
    private S2MinDistanceTarget target;

    /////////////// Internal state /////////////////

    // Temporaries, defined here to avoid multiple allocations / initializations.

    private S2Iterator<S2ShapeIndex.Cell> iter;
    private final ArrayList<S2CellId> maxDistanceCovering = Lists.newArrayList();
    private final ArrayList<S2CellId> initialCells = Lists.newArrayList();

    // For the optimized algorihm we precompute the top-level S2CellIds that will be added to the priority queue.
    // There can be at most 6 of these cells. Essentially this is just a covering of the indexed edges, except
    // that we also store pointers to the corresponding S2ShapeIndexCells to reduce the number of index seeks required.
    private final ArrayList<S2CellId> indexCovering = Lists.newArrayListWithCapacity(6);
    private final ArrayList<S2ShapeIndex.Cell> indexCells = Lists.newArrayListWithCapacity(6);

    // The decision about whether to use the brute force algorithm is based on counting the total number of edges in
    // the index.  However if the index contains a large number of shapes, this in itself might take too long.
    // So instead we only count edges up to (maxBruteForceIndexSize() + 1) for the current target type
    // (stored as indexNumEdgesLimit).
    private int indexNumEdges = 0;
    private int indexNumEdgesLimit = 0;

    /**
     * The distance beyond which we can safely ignore further candidate edges. (Candidates that are exactly at the limit
     * are ignored; this is more efficient for updateMinDistance() and should not affect clients since distance
     * measurements have a small amount of errornyway.)
     * <p>
     * Initially this is the same as the maximum distance specified by the user, but it can also be updated by the
     * algorithm (see maybeAddResult).
     */
    private S1ChordAngle distanceLimit;

    /**
     * The algorithm maintains a priority queue of unprocessed S2CellIds, sorted in increasing order of distance from
     * the target.
     */
    private final Queue<QueueEntry> queue = new PriorityQueue<>(16);

    /////////////// Results /////////////////

    // The current result set is stored in one of three ways:
    //
    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.

    /** If maxResults == 1, keep the best result */
    private Result resultSingleton;

    /** If maxResults() == "infinity", results are appended to this list and sorted/uniqued at the end. */
    private final ArrayList<Result> resultVector = Lists.newArrayList();

    /**
     * Otherwise results are kept in a TreeSet so that we can progressively reduce the distance limit once maxResults
     * results have been found. (A priority queue is not sufficient because we need to be able to check whether a
     * candidate edge is already in the result set.)
     */
    private final TreeSet<Result> resultSet = Sets.newTreeSet();

    /////////////////////////////////////////

    /**
     * Construct a new query for the given index. Must call reset() before using the query, if the
     * index has been modified since the query was constructed.
     */
    public S2ClosestEdgeQuery(S2ShapeIndex index, Options options) {
        this.index = index;
        this.options = options;
        reset();
    }

    public S2ClosestEdgeQuery(S2ShapeIndex index) {
        this(index, new Options());
    }

    /** Reinitializes the query. This method must be called whenever the underlying index is modified. */
    public void reset() {
        indexNumEdges = 0;
        indexNumEdgesLimit = 0;
        indexCovering.clear();
        indexCells.clear();
        // We don't initialize iter here to make queries on small indexes a bit faster (i.e., where brute force is used).
    }

    public S2ShapeIndex index() {
        return index;
    }

    public Options getOptions() {
        return options;
    }

    /**
     * Sets whether distances are computed using "brute force" (i.e., by examining every point) rather
     * than using the S2PointIndex.
     *
     * <p>This is package private, as it is intended only for testing, benchmarking, and debugging.
     *
     * <p>Do not call before init().
     */
    @VisibleForTesting
    void useBruteForce(boolean useBruteForce) {
        this.useBruteForce = useBruteForce;
    }

    /**
     * Returns the closest edges to the given target that satisfy the current options.
     * This method may be called multiple times.
     * <p>
     * Note that if includeInteriors() is true, the result vector may include some entries with edgeId == -1.
     * This indicates that the target intersects the indexed polygon.
     */
    public List<Result> findClosestEdges(S2Point target) {
        return findClosestEdges(new PointTarget(target));
    }

    /**
     * This version can be more efficient when this method is called many times, since it does not require
     * allocating a new list on each call.
     */
    public void findClosestEdges(List<Result> results, S2Point target) {
        findClosestEdges(results, new PointTarget(target));
    }

    /**
     * Returns the closest edges to the given target that satisfy the current options.
     * This method may be called multiple times.
     * <p>
     * Note that if includeInteriors() is true, the result vector may include some entries with edgeId == -1.
     * This indicates that the target intersects the indexed polygon with the given shape_id.
     */
    public List<Result> findClosestEdgesToEdge(S2Point a, S2Point b) {
        return findClosestEdges(new EdgeTarget(a, b));
    }

    /**
     * This version can be more efficient when this method is called many times,
     * since it does not require allocating a new vector on each call.
     */
    public void findClosestEdgesToEdge(List<Result> results, S2Point a, S2Point b) {
        findClosestEdges(results, new EdgeTarget(a, b));
    }

    /**
     * Returns the closest edge to the target point. If no edge satisfies the search criteria, then the Result
     * object will have distance == S1ChordAngle.INFINITY, isEmpty() == true, shape == null and edgeId == -1.
     *
     * <p>
     * Note that if includeInteriors is true, edgeId == -1 is also used to indicate that the target intersects
     * an indexed polygon (but in that case distance == S1ChordAngle.ZERO and shape != null).
     */
    public Result findClosestEdge(S2Point target) {
        return findClosestEdge(new PointTarget(target));
    }



    /**
     * Returns the closest edge to the target edge. If no edge satisfies the search criteria, then the Result
     * object will have distance == S1ChordAngle.INFINITY, isEmpty() == true, shape == null and edgeId == -1.
     *
     * <p>
     * Note that if includeInteriors is true, edgeId == -1 is also used to indicate that the target intersects
     * an indexed polygon (but in that case distance == S1ChordAngle.ZERO and shape != null).
     */
    public Result findClosestEdgeToEdge(S2Point a, S2Point b) {
        return findClosestEdge(new EdgeTarget(a, b));
    }

    /**
     * Finds the closest edges to the given target that satisfy the given options. This method may be called multiple
     * times.
     * <p>
     * Note that if options.includeInteriors is true, the result vector may include some entries with edge_id == -1.
     * This indicates that the target intersects the indexed polygon with the given shapeId.
     *
     * @param target The query target.
     * @return The query result ordered by increase distance.
     */
    private List<Result> findClosestEdges(S2MinDistanceTarget target) {
        ArrayList<Result> results = new ArrayList<>();
        findClosestEdges(results, target);
        return results;
    }

    /**
     * Finds the closest edges to the given target that satisfy the given options. This version can be more efficient
     * when this method is called many times, since it does not require allocating a new list on each call.
     *
     * @param target  The query target.
     * @param results The list to populate with the query result ordered by increase distance.
     */
    void findClosestEdges(List<Result> results, S2MinDistanceTarget target) {
        findClosestEdgesInternal(target, options);
        results.clear();
        if (options.maxResults == 1) {
            if (resultSingleton.shape != null) {
                results.add(resultSingleton);
            }
        } else if (options.maxResults == Integer.MAX_VALUE) {
            results.addAll(resultVector);
            sortAndRemoveDuplicates(results);
            resultVector.clear();
        } else {
            results.addAll(resultSet);
            resultSet.clear();
        }
    }

    /**
     * Convenience method that returns exactly one edge. If no edges satisfy the given search criteria, then a
     * Result with distance == infinity() and shapeId == edgeId == -1 is returned.
     * <p>
     * Note that if options.includeInteriors is true, edge_id == -1 is also used to indicate that the target intersects
     * an indexed polygon (but in that case distance == Zero() and shapeId >= 0).
     * <p>
     * REQUIRES: options.maxResults == 1
     *
     * @param target The query target.
     * @return The closest edge result.
     */
    Result findClosestEdge(S2MinDistanceTarget target, Options options) {
        Options currentOptions = this.options;
        Options singleResultOptions = options.copy().setMaxResults(1);
        findClosestEdgesInternal(target, singleResultOptions);
        this.options = currentOptions;
        return resultSingleton;
    }

    Result findClosestEdge(S2MinDistanceTarget target) {
        return findClosestEdge(target, options);
    }

    S1ChordAngle getDistance(S2MinDistanceTarget target) {
        return findClosestEdge(target).distance;
    }

    /**
     * Returns the minimum distance to the target point.  If the index or target is empty, returns
     * S1ChordAngle.INFINITY.
     *
     * Use IsDistanceLess() if you only want to compare the distance against a threshold value,
     * since it is often much faster.
     */
    public S1ChordAngle getDistance(S2Point target) {
        return getDistance(new PointTarget(target));
    }

    /**
     * Returns the minimum distance to the target edge.  If the index or target is empty, returns
     * S1ChordAngle.INFINITY.
     *
     * Use IsDistanceLess() if you only want to compare the distance against a threshold value,
     * since it is often much faster.
     */
    public S1ChordAngle getDistance(S2Point a, S2Point b) {
        return findClosestEdgeToEdge(a, b).distance;
    }

    /**
     * Returns true if the distance to "target" is less than "limit".
     *
     * This method is usually much faster than GetDistance(), since it is much less work to determine whether
     * the minimum distance is above or below a threshold than it is to calculate the actual minimum distance.
     */
    public boolean isDistanceLess(S2Point target, S1ChordAngle limit) {
        return isDistanceLess(new PointTarget(target), limit);
    }

    boolean isDistanceLess(S2MinDistanceTarget target, S1ChordAngle limit) {
        Options currentOptions = options;
        Options distanceLessOptions = options.copy()
                .setMaxResults(1)
                .setMaxDistance(limit)
                .setMaxError(S1ChordAngle.STRAIGHT);
        Result result = findClosestEdge(target, distanceLessOptions);
        options = currentOptions;
        return result.isNotEmpty();
    }

    /**
     * Like IsDistanceLess(), but also returns true if the distance to "target"is exactly equal to "limit".
     */
    public boolean isDistanceLessOrEqual(S2Point target, S1ChordAngle limit) {
        Options currentOptions = options;
        Options distanceLessOrEqualOptions = options.copy()
                .setMaxResults(1)
                .setInclusiveMaxDistance(limit)
                .setMaxError(S1ChordAngle.STRAIGHT);
        Result result = findClosestEdge(new PointTarget(target), distanceLessOrEqualOptions);
       options = currentOptions;
        return result.isNotEmpty();
    }

    /**
     * Like IsDistanceLessOrEqual(), except that "limit" is increased by the maximum error in the
     * distance calculation.  This ensures that this function returns true whenever the true, exact
     * distance is less than or equal to "limit".
     */
    public boolean isConservativeDistanceLessOrEqual(S2Point target, S1ChordAngle limit) {
        return isConservativeDistanceLessOrEqual(new PointTarget(target), limit);
    }

    boolean isConservativeDistanceLessOrEqual(S2MinDistanceTarget target, S1ChordAngle limit) {
        Options currentOptions = options;
        Options conservativeDistanceLessOrEqualOptions = options.copy()
                .setMaxResults(1)
                .setConservativeMaxDistance(limit)
                .setMaxError(S1ChordAngle.STRAIGHT);
        Result result = findClosestEdge(target, conservativeDistanceLessOrEqualOptions);
        options = currentOptions;
        return result.isNotEmpty();
    }

    /**
     * Returns the endpoints of the given result edge.
     *
     * CAVEAT: If includeInteriors is true, then clients must not pass this method any Result objects
     * that correspond to shape interiors, i.e. those where result.edgeId < 0.
     *
     * REQUIRES: result.edgeId >= 0
     */
    public void getEdge(Result result, S2Shape.MutableEdge edge) {
        Preconditions.checkArgument(result.edgeId >= 0);
        result.shape.getEdge(result.edgeId, edge);
    }

    /** Returns the point on given result edge that is closest to "point". */
    public S2Point project(S2Point point, Result result) {
        if (result.edgeId < 0) return point;
        S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
        getEdge(result, edge);
        return S2EdgeUtil.getClosestPoint(point, edge.a, edge.b);
    }


    ///////////////////////////////// Internal ////////////////////////////

    private void findClosestEdgesInternal(S2MinDistanceTarget target, Options options) {
        this.target = target;
        this.options = options;

        distanceLimit = options.maxDistance;
        resultSingleton = new Result(S1ChordAngle.INFINITY);
        Preconditions.checkState(resultVector.isEmpty());
        Preconditions.checkState(resultSet.isEmpty());
        Preconditions.checkState(target.maxBruteForceIndexSize() >= 0);
        if (distanceLimit == S1ChordAngle.ZERO) return;

        if (options.maxResults == Integer.MAX_VALUE && options. maxDistance == S1ChordAngle.INFINITY) {
            log.warning("Returning all edges (maxResults/maxDistance not set)");
        }

        // If includeInteriors == true. Process all the shapes contained by the target.
        if (options.includeInteriors) {
            Set<S2Shape> shapes = new HashSet<>();
            target.visitContainingShapes(index, input -> {
                shapes.add(input);
                return shapes.size() < options.maxResults;
            });
            for (S2Shape shape : shapes) {
                addResult(new Result(S1ChordAngle.ZERO, shape, -1));
            }
            if (distanceLimit == S1ChordAngle.ZERO) return;
        }

        // Use the brute force algorithm if the index is small enough.  To avoid
        // spending too much time counting edges when there are many shapes, we stop
        // counting once there are too many edges.  We may need to recount the edges
        // if we later see a target with a larger brute force edge threshold.
        int minOptimizedEdges = target.maxBruteForceIndexSize() + 1;
        if (minOptimizedEdges > indexNumEdgesLimit && indexNumEdges >= indexNumEdgesLimit) {
            indexNumEdges = countEdgesUpTo(index, minOptimizedEdges);
            indexNumEdgesLimit = minOptimizedEdges;
        }

        if (useBruteForce || indexNumEdges < minOptimizedEdges) {
            findClosestEdgesBruteForce();
        } else {
            findClosestEdgesOptimized();
        }
    }

    private void findClosestEdgesBruteForce() {
        for (S2Shape shape : index.shapes) {
            if (shape == null) continue;
            int numEdges = shape.numEdges();
            for (int e = 0; e < numEdges; ++e) {
                maybeAddResult(shape, e);
            }
        }
    }

    private void findClosestEdgesOptimized() {
        initQueue();
        // Repeatedly find the closest S2Cell to "target" and either split it into
        // its four children or process all of its edges.
        while (!queue.isEmpty()) {
            // We need to copy the top entry before removing it, and we need to
            // remove it before adding any new entries to the queue.
            QueueEntry entry = queue.poll();

            // Queue is sorted by distance. So if the current entry distance is greater to distanceLimit
            // we have collect all the matching edges.
            if (entry.distance.compareTo(distanceLimit) >= 0) {
                queue.clear();  // Clear any remaining entries.
                break;
            }

            // If this is already known to be an index cell, just process it.
            if (entry.indexCell != null) {
                processEdges(entry);
                continue;
            }

            // Otherwise split the cell into its four children.  Before adding a
            // child back to the queue, we first check whether it is empty.  We do
            // this in two seek operations rather than four by seeking to the key
            // between children 0 and 1 and to the key between children 2 and 3.
            S2CellId id = entry.id;
            iter.seek(id.child(1).rangeMin());
            if (!iter.done() && iter.id().lessOrEquals(id.child(1).rangeMax())) {
                processOrEnqueue(id.child(1));
            }
            if (!iter.atBegin()) {
                iter.prev();
                if (iter.id().greaterOrEquals(id.rangeMin())) {
                    processOrEnqueue(id.child(0));
                }
            }
            iter.seek(id.child(3).rangeMin());
            if (!iter.done() && iter.id().lessOrEquals(id.rangeMax())) {
                processOrEnqueue(id.child(3));
            }
            if (!iter.atBegin()) {
                iter.prev();
                if (iter.id().greaterOrEquals(id.child(2).rangeMin())) {
                    processOrEnqueue(id.child(2));
                }
            }
        }
    }

    private void initQueue() {
        Preconditions.checkState(queue.isEmpty());
        if (indexCovering.isEmpty()) {
            // We delay iterator initialization until now to make queries on very
            // small indexes a bit faster (i.e., where brute force is used).
            iter = index.iterator();
        }

        // Optimization: if the user is searching for just the closest edge, and the
        // center of the target's bounding cap happens to intersect an index cell,
        // then we try to limit the search region to a small disc by first
        // processing the edges in that cell.  This sets distance_limit_ based on
        // the closest edge in that cell, which we can then use to limit the search
        // area.  This means that the cell containing "target" will be processed
        // twice, but in general this is still faster.
        //
        // TODO(ericv): Even if the cap center is not contained, we could still
        // process one or both of the adjacent index cells in S2CellId order,
        // provided that those cells are closer than distance_limit_.
        S2Cap cap = target.getCapBound();
        if (cap.isEmpty()) return;  // Empty target.
        if (options.maxResults == 1 && iter.locate(cap.axis())) {
            processEdges(new QueueEntry(S1ChordAngle.ZERO, iter.id(), iter.entry()));
            // Skip the rest of the algorithm if we found an intersecting edge.
            if (distanceLimit == S1ChordAngle.ZERO) return;
        }
        if (indexCovering.isEmpty()) initCovering();
        if (distanceLimit == S1ChordAngle.INFINITY) {
            // Start with the precomputed index covering.
            for (int i = 0; i < indexCovering.size(); ++i) {
                processOrEnqueue(indexCovering.get(i), indexCells.get(i));
            }
        } else {
            // Compute a covering of the search disc and intersect it with the
            // precomputed index covering.
            S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(4).build();
            S1ChordAngle radius = S1ChordAngle.fromLength2(cap.radius().getLength2() + distanceLimit.plusError(distanceLimit.getS1AngleConstructorMaxError()).getLength2());
            S2Cap searchCap = S2Cap.fromAxisChord(cap.axis(), radius);
            coverer.getFastCovering(searchCap, maxDistanceCovering);
            S2CellUnion.getIntersection(indexCovering, maxDistanceCovering, initialCells);

            // Now we need to clean up the initial cells to ensure that they all
            // contain at least one cell of the S2ShapeIndex.  (Some may not intersect
            // the index at all, while other may be descendants of an index cell.)
            int i = 0;
            int j = 0;
            while (i < initialCells.size()) {
                S2CellId idi = initialCells.get(i);
                // Find the top-level cell that contains this initial cell.
                while (indexCovering.get(j).rangeMax().lessThan(idi)) {
                    ++j;
                }
                S2CellId idj = indexCovering.get(j);
                if (idi.equals(idj)) {
                    // This initial cell is one of the top-level cells.  Use the
                    // precomputed S2ShapeIndexCell pointer to avoid an index seek.
                    processOrEnqueue(idj, indexCells.get(j));
                    ++i;
                    ++j;
                } else {
                    // This initial cell is a proper descendant of a top-level cell.
                    // Check how it is related to the cells of the S2ShapeIndex.
                    S2ShapeIndex.CellRelation r = iter.locate(idi);
                    if (r == S2ShapeIndex.CellRelation.INDEXED) {
                        // This cell is a descendant of an index cell.  Enqueue it and skip
                        // any other initial cells that are also descendants of this cell.
                        processOrEnqueue(iter.id(), iter.entry());
                        S2CellId lastid = iter.id().rangeMax();
                        do {
                            ++i;
                        } while (i < initialCells.size() && initialCells.get(i).lessOrEquals(lastid));
                    } else {
                        // Enqueue the cell only if it contains at least one index cell.
                        if (r == S2ShapeIndex.CellRelation.SUBDIVIDED) processOrEnqueue(idi, null);
                        ++i;
                    }
                }
            }
        }
    }

    private void initCovering() {
        // Find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  These are
        // the "top-level" cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one top-level cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a top-level cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the top-level cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.

        // Don't need to reserve index_cells_ since it is an InlinedVector.
        indexCovering.ensureCapacity(6);

        // TODO(ericv): Use a single iterator (iter_) below and save position
        // information using pair<S2CellId, const S2ShapeIndexCell*> type.
        S2Iterator<S2ShapeIndex.Cell> next = index.iterator();
        S2Iterator<S2ShapeIndex.Cell> last = index.iterator();
        last.finish(); last.prev();

        if (!next.id().equals(last.id())) {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            int level = next.id().getCommonAncestorLevel(last.id()) + 1;

            // Visit each potential top-level cell except the last (handled below).
            S2CellId lastId = last.id().parent(level);
            for(S2CellId id = next.id().parent(level); !id.equals(lastId); id = id.next()) {
                // Skip any top-level cells that don't contain any index cells.
                if (id.rangeMax().lessThan(next.id())) {
                    continue;
                }

                // Find the range of index cells contained by this top-level cell and
                // then shrink the cell if necessary so that it just covers them.
                S2CellId fid = next.id();
                S2ShapeIndex.Cell fCell = next.entry();
                next.seek(id.rangeMax().next());
                next.prev();
                S2CellId lid = next.id();
                S2ShapeIndex.Cell lCell = next.entry();
                addInitialRange(fid, fCell, lid, lCell);
                next.next();
            }
        }
        addInitialRange(next.id(), next.entry(), last.id(), last.entry());
    }

    // Add an entry to index_covering_ and index_cells_ that covers the given
    // inclusive range of cells.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private void addInitialRange(S2CellId firstId, S2ShapeIndex.Cell firstCell, S2CellId lastId, S2ShapeIndex.Cell lastCell) {
        if (firstId.equals(lastId)) {
            // The range consists of a single index cell.
            indexCovering.add(firstId);
            indexCells.add(firstCell);
        } else {
            // Add the lowest common ancestor of the given range.
            int level = firstId.getCommonAncestorLevel(lastId);
            Preconditions.checkState(level >= 0);
            indexCovering.add(firstId.parent(level));
            indexCells.add(null);
        }
    }

    private void maybeAddResult(S2Shape shape, int edge_id) {
        S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
        shape.getEdge(edge_id, edge);

        S1ChordAngle distance = target.getMinDistance(edge.a, edge.b, distanceLimit);
        if (distance.compareTo(distanceLimit) < 0) {
            addResult(new Result(distance, shape, edge_id));
        }
    }

    private void addResult(Result result) {
        if (options.maxResults == 1) {
            // Optimization for the common case where only the closest edge is wanted.
            resultSingleton = result;
            distanceLimit = S1ChordAngle.sub(result.distance, options.maxError);
        } else if (options.maxResults == Integer.MAX_VALUE) {
            resultVector.add(result);  // Sort/unique at end.
        } else {
            // Add this edge to result_set_.  Note that even if we already have enough
            // edges, we can't erase an element before insertion because the "new"
            // edge might in fact be a duplicate.
            resultSet.add(result);
            int size = resultSet.size();
            if (size >= options.maxResults) {
                if (size > options.maxResults) {
                    resultSet.pollLast();
                }
                distanceLimit = S1ChordAngle.sub(resultSet.last().distance, options.maxError);
            }
        }
    }

    private void processEdges(QueueEntry entry) {
        S2ShapeIndex.Cell indexCell = entry.indexCell;
        for (int s = 0; s < indexCell.numShapes(); ++s) {
            S2ShapeIndex.S2ClippedShape clipped = indexCell.clipped(s);
            S2Shape shape = clipped.shape();
            if (shape != null) {
                for (int j = 0; j < clipped.numEdges(); ++j) {
                    maybeAddResult(shape, clipped.edge(j));
                }
            }
        }
    }

    // Enqueue the given cell id.
    // REQUIRES: iter_ is positioned at a cell contained by "id".
    private void processOrEnqueue(S2CellId id) {
        Preconditions.checkArgument(id.contains(iter.id()));
        if (iter.id().equals(id)) {
            processOrEnqueue(id, iter.entry());
        } else {
            processOrEnqueue(id, null);
        }
    }

    // Add the given cell id to the queue.  "index_cell" is the corresponding
    // S2ShapeIndexCell, or nullptr if "id" is not an index cell.
    //
    // This version is called directly only by InitQueue().
    private void processOrEnqueue(S2CellId id, S2ShapeIndex.Cell index_cell) {
        if (index_cell != null) {
            // If this index cell has only a few edges, then it is faster to check
            // them directly rather than computing the minimum distance to the S2Cell
            // and inserting it into the queue.
            int kMinEdgesToEnqueue = 10;
            int num_edges = countEdges(index_cell);
            if (num_edges == 0) return;
            if (num_edges < kMinEdgesToEnqueue) {
                // Set "distance" to zero to avoid the expense of computing it.
                processEdges(new QueueEntry(S1ChordAngle.ZERO, id, index_cell));
                return;
            }
        }
        // Otherwise compute the minimum distance to any point in the cell and add
        // it to the priority queue.
        S2Cell cell = new S2Cell(id);
        S1ChordAngle distance = target.getMinDistance(cell, distanceLimit);
        if (distance.compareTo(distanceLimit) >= 0) {
            return;
        }
        queue.offer(new QueueEntry(distance, id, index_cell));
    }

    /** Sort in place the given list in natural order and remove duplicates values. */
    private static <T extends Comparable<T>> void sortAndRemoveDuplicates(List<T> list) {
        list.sort(Comparator.naturalOrder());
        int idx = list.size() - 1;
        while (idx >= 1) {
            if (list.get(idx).equals(list.get(idx - 1))) {
                list.remove(idx);
            }
            --idx;
        }
    }

    /** Return the number of edges in the given index cell. */
    private static int countEdges(S2ShapeIndex.Cell cell) {
        int count = 0;
        for (int s = 0; s < cell.numShapes(); ++s) {
            count += cell.clipped(s).numEdges();
        }
        return count;
    }

    @VisibleForTesting
    static int countEdges(S2ShapeIndex index) {
        int numEdges = 0;
        for (S2Shape shape : index.getShapes()) {
            numEdges += shape.numEdges();
        }
        return numEdges;
    }

    /**
     * Return the number of edges in the given index, but stops once "maxEdges" edges have been found (in
     * which case the current running total is returned).
     */
    private static int countEdgesUpTo(S2ShapeIndex index, int maxEdges) {
        int numEdges = 0;
        for (S2Shape shape : index.getShapes()) {
            numEdges += shape.numEdges();
            if (numEdges >= maxEdges) return numEdges;
        }
        return numEdges;
    }

    // Target subtype that computes the closest distance to a point.
    static class PointTarget extends S2MinDistancePointTarget {

        public PointTarget(S2Point point) {
            super(point);
        }

        @Override
        public int maxBruteForceIndexSize() {
            // Using BM_FindClosest (which finds the single closest edge), the
            // break-even points are approximately 80, 100, and 250 edges for point
            // cloud, fractal, and regular loop geometry respectively.
            return 120;
        }
    }

    // Target subtype that computes the closest distance to an edge.
    static class EdgeTarget extends S2MinDistanceEdgeTarget {

        public EdgeTarget(S2Point a, S2Point b) {
            super(a, b);
        }

        @Override
        public int maxBruteForceIndexSize() {
            // Using BM_FindClosestToEdge (which finds the single closest edge), the
            // break-even points are approximately 40, 50, and 100 edges for point
            // cloud, fractal, and regular loop geometry respectively.
            return 60;
        }
    }

    public static class Options {

        /** The max distance to search for points. */
        private S1ChordAngle maxDistance;

        /** Max distance error. */
        private S1ChordAngle maxError;

        /** Indicates if polygon interiors should be included when measuring distances. */
        private boolean includeInteriors;

        /** The max number of returned results. */
        private int maxResults = Integer.MAX_VALUE;

        public Options() {
            maxDistance = S1ChordAngle.INFINITY;
            maxError = S1ChordAngle.ZERO;
            includeInteriors = true;
        }

        /**
         * Returns the max distance between returned points and the given target. Default is +inf.
         */
        public S1ChordAngle getMaxDistance() {
            return maxDistance;
        }

        /**
         * Sets a new max distance to search for shapes.
         */
        public Options setMaxDistance(S1ChordAngle maxDistance) {
            this.maxDistance = maxDistance;
            return this;
        }

        /**
         * Sets a new max distance to search for shapes.
         */
        public Options setMaxDistance(S1Angle maxDistance) {
            this.maxDistance = S1ChordAngle.fromS1Angle(maxDistance);
            return this;
        }

        /**
         * Like setMaxDistance(), except that edges whose distance is exactly equal to "maxDistance" are also returned.
         * Equivalent to calling setMaxDistance(maxDistance.successor()).
         */
        public Options setInclusiveMaxDistance(S1ChordAngle maxDistance) {
            setMaxDistance(maxDistance.successor());
            return this;
        }

        /**
         * Like setInclusiveMaxDistance(), except that "maxDistance" is also increased by the maximum error in
         * the distance calculation.  This ensures that all edges whose true distance is less than or equal to
         * "maxFistance" will be returned (along with some edges whose true distance is slightly greater).
         * <p>
         * Algorithms that need to do exact distance comparisons can use this option to find a set of candidate
         * edges that can then be filtered further (e.g., using s2pred::CompareDistance).
         */
        public Options setConservativeMaxDistance(S1ChordAngle maxDistance) {
            this.maxDistance = maxDistance.plusError(S2EdgeUtil.getMinDistanceMaxError(maxDistance)).successor();
            return this;
        }

        /**
         * Set inclusive max distance given as a {@link S1Angle}.
         */
        public Options setInclusiveMaxDistance(S1Angle maxDistance) {
            setInclusiveMaxDistance(S1ChordAngle.fromS1Angle(maxDistance));
            return this;
        }

        /**
         * Set conservative max distance given as a {@link S1Angle}.
         */
        public Options setConservativeMaxDistance(S1Angle maxDistance) {
            setConservativeMaxDistance(S1ChordAngle.fromS1Angle(maxDistance));
            return this;
        }

        /**
         * Gets the max error distance.
         */
        public S1ChordAngle getMaxError() {
            return maxError;
        }

        /**
         * Specifies that edges up to maxError further away than the true closest edges may be substituted in the
         * result set, as long as such edges satisfy all the remaining search criteria (such as max_distance).
         * This option only has an effect if max_results() is also specified; otherwise all edges closer than
         * maxDistance will always be returned.
         * <p>
         * Note that this does not affect how the distance between edges is computed; it simply gives the algorithm
         * permission to stop the search early as soon as the best possible improvement drops below maxError.
         * <p>
         * This can be used to implement distance predicates efficiently.  For example, to determine whether the minimum
         * distance is less than D, set maxResults == 1 and maxDistance == maxError == D.  This causes the algorithm to
         * terminate as soon as it finds any edge whose distance is less than D, rather than continuing to search for
         * an edge that is even closer.
         */
        public Options setMaxError(S1ChordAngle maxError) {
            this.maxError = maxError;
            return this;
        }

        /** Sets the maxError given as {@link S1Angle}. */
        public Options setMaxError(S1Angle maxError) {
            this.maxError = S1ChordAngle.fromS1Angle(maxError);
            return this;
        }

        /**
         * Indicates if polygon interiors must be included.
         */
        public boolean isIncludeInteriors() {
            return includeInteriors;
        }

        /**
         * Sets the include interiors option. Specifies that polygon interiors should be included when measuring
         * distances. In other words, polygons that contain the target should have a distance of zero.
         * (For targets consisting of multiple connected components, the distance is zero if any component is
         * contained.)
         * This is indicated in the results by returning a (shapeId, edgeId) pair with edgeId == -1, i.e.
         * this value denotes the polygons's interior.
         * <p>
         * Note that for efficiency, any polygon that intersects the target may or may not have an (edgeId == -1)
         * result.  Such results are optional because in that case the distance to the polygon is already zero.
         */
        public Options setIncludeInteriors(boolean includeInteriors) {
            this.includeInteriors = includeInteriors;
            return this;
        }

        /**
         * Specifies that at most "maxResults" edges should be returned.
         */
        public Options setMaxResults(int maxResults) {
            this.maxResults = maxResults;
            return this;
        }

        /** Get the current max results. */
        public int getMaxResults() {
            return maxResults;
        }

        public Options copy() {
            return new Options()
                    .setMaxDistance(maxDistance)
                    .setIncludeInteriors(includeInteriors)
                    .setMaxError(maxError)
                    .setMaxResults(maxResults);
        }
    }

    /**
     * Each "Result" object represents a closest edge. Note the following special cases:
     * <p>
     * - (shapeId >= 0) && (edgeId < 0) represents the interior of a shape.
     * Such results may be returned when options.includeInteriors is true.
     * Such results can be identified using the isInterior() method.
     * <p>
     * - (shapeId < 0) && (edgeId < 0) is returned by `findClosestEdge`  to indicate that no edge satisfies the given query
     * options.  Such results can be identified using isEmpty() method.
     */
    public static class Result implements Comparable<Result> {
        public final S1ChordAngle distance;
        public final S2Shape shape;
        public final int edgeId;

        public Result(S1ChordAngle distance, S2Shape shape, int edgeId) {
            this.distance = distance;
            this.shape = shape;
            this.edgeId = edgeId;
        }

        public Result(S1ChordAngle distance) {
            this(distance, null, -1);
        }

        /**
         * Indicates if the result represents the interior of a shape.
         * (Such results may be returned when options.includeInteriors is true.)
         *
         * @return true if this Result object represents the interior of a shape.
         */
        public boolean isInterior() {
            return shape != null && edgeId < 0;
        }

        /**
         * Indicates if the result is empty.
         * (This result is only returned in one special case, namely when findClosestEdge() does not find any suitable edges.
         * It is never returned by methods that return a list of results.)
         *
         * @return true if this Result object indicates that no edge satisfies the given query options.
         */
        public boolean isEmpty() {
            return shape == null;
        }

        public boolean isNotEmpty() {
            return shape != null;
        }

        @Override
        public String toString() {
            return "Result{" +
                    "distance=" + distance +
                    ", edgeId=" + edgeId +
                    ", shape=" + shape +
                    '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Result)) return false;
            Result result = (Result) o;
            return edgeId == result.edgeId && Objects.equal(distance, result.distance) && Objects.equal(shape, result.shape);
        }

        @Override
        public int hashCode() {
            return Objects.hashCode(distance, shape, edgeId);
        }

        // Compares edges first by distance, then by (shape_id, edge_id).
        @Override
        public int compareTo(Result other) {
            return ComparisonChain.start()
                    .compare(distance, other.distance)
                    .compare(System.identityHashCode(shape), System.identityHashCode(other.shape))
                    .compare(edgeId, other.edgeId)
                    .result();
        }

    }

    /**
     * ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
     * consisting of a (shapeId, edgeId) pair.
     */
    static class ShapeEdgeId implements Comparable<ShapeEdgeId> {

        private final int shapeId;
        private final int edgeId;

        ShapeEdgeId(int shapeId, int edgeId) {
            this.shapeId = shapeId;
            this.edgeId = edgeId;
        }

        @Override
        public int compareTo(ShapeEdgeId other) {
            return ComparisonChain.start()
                    .compare(shapeId, other.shapeId)
                    .compare(edgeId, other.edgeId)
                    .result();
        }
    }

    private static class QueueEntry implements Comparable<QueueEntry> {

        // A lower bound on the distance from the target to "id".  This is the key
        // of the priority queue.
        private final S1ChordAngle distance;

        // The cell being queued.
        private final S2CellId id;

        // If "id" belongs to the index, this field stores the corresponding
        // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
        // S2ShapeIndexCells and this field stores nullptr.  The purpose of this
        // field is to avoid an extra Seek() when the queue entry is processed.
        private final S2ShapeIndex.Cell indexCell;

        private QueueEntry(S1ChordAngle distance, S2CellId id, S2ShapeIndex.Cell indexCell) {
            this.distance = distance;
            this.id = id;
            this.indexCell = indexCell;
        }

        // The priority queue returns the smallest elements first
        @Override
        public int compareTo(QueueEntry other) {
            return distance.compareTo(other.distance);
        }

    }

}
