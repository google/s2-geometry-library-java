/*
 * Copyright 2014 Google Inc.
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

import static java.lang.Math.abs;
import static java.lang.Math.max;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.RandomAccess;
import javax.annotation.Nullable;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsEnum;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * S2ShapeIndex indexes a set of "shapes", where a shape is a collection of edges that optionally
 * defines an interior. A shape can be as simple as a single edge, or as complex as a collection of
 * loops. For shapes that have interiors, the index makes it very fast to determine the shape(s)
 * that contain a given point or region.
 *
 * <p>It is designed to handle up to millions of shapes and edges. All data structures are designed
 * to be small, so the index is compact; generally it is a 20-30% of the size of the underlying data
 * being indexed. The index is also very fast to construct. The index size and construction time are
 * guaranteed to be linear in the number of input edges.
 *
 * <p>The index is not yet dynamic, but it is planned to allow shapes to be inserted or removed at
 * any time (although each individual shape will have to be immutable). Presently, shapes may be
 * inserted until the first call to {@link #iterator()}, after which no more changes may be made.
 *
 * <p>There are a number of built-in classes that work with S2ShapeIndex objects. Generally these
 * classes accept any collection of geometry that can be represented by an S2ShapeIndex, i.e. any
 * combination of points, polylines, and polygons. Such classes include:
 * <ul>
 * <li>S2ContainsPointQuery: returns the shape(s) that contain a given point.
 * <li>S2EdgeQuery: returns the edges or shapes that are crossed by a given edge.
 * </ul>
 *
 * <p>The shape index is essentially a map from S2CellId to a set of clipped shapes that intersect
 * that cell id. It is adaptively refined to ensure that no cell contains more than {@link
 * Options#getMaxEdgesPerCell()} edges.
 *
 * <p>Insertions and removals are queued and processed on the first subsequent query. There are
 * several reasons to do this:
 *
 * <ul>
 *   <li>It is significantly more efficient to process updates in batches.
 *   <li>Often the index will never be queried, in which case we can save both the time and memory
 *       required to build it.
 * </ul>
 *
 * Example cases where lazy construction works especially well:
 *
 * <ul>
 *   <li>S2Loops that are created simply to pass to an S2Polygon. (We don't need the S2Loop index,
 *       because S2Polygon builds its own index.)
 *   <li>Applications that load a database of geometry and then query only a small fraction of it.
 *   <li>Applications that only read and write geometry.
 * </ul>
 *
 * <p>The main drawback to lazily building the index is that we need to go to some extra work to
 * ensure that the index behaves as though it is an immutable object from the perspective of an
 * external observer. Index construction is done lazily, but guarded by a lock on the S2ShapeIndex
 * instance.
 *
 * <p>Shapes added to the index are expected to be deeply immutable, or else not mutated. The state
 * of the index is undefined if a shape is mutated after passing it to {@link #add}.
 */
@JsType
public strictfp class S2ShapeIndex implements Serializable {
  private static final long serialVersionUID = 1L;

  /**
   * The amount in UV coordinates by which cells are "padded" to compensate for numerical errors
   * when clipping line segments to cell boundaries. The total error when clipping an edge comes
   * from two sources:
   *
   * <ol>
   *   <li>Clipping the original spherical edge to a cube face (the "face edge"). The maximum error
   *       in this step is {@link S2EdgeUtil#FACE_CLIP_ERROR_UV_COORD}.
   *   <li>Clipping the face edge to the u- or v-coordinate of a cell boundary. The maximum error in
   *       this step is {@link S2EdgeUtil#EDGE_CLIP_ERROR_UV_COORD}.
   * </ol>
   *
   * <p>Finally, since we encounter the same errors when clipping query edges, we double the total
   * error so that we only need to pad edges during indexing and not at query time.
   */
  public static final double CELL_PADDING =
      2 * (S2EdgeUtil.FACE_CLIP_ERROR_UV_COORD + S2EdgeUtil.EDGE_CLIP_ERROR_UV_COORD);

  /**
   * Default maximum number of edges per cell (not counting 'long' edges). Reasonable values range
   * from 10 to 50. Small values makes queries faster, while large values make construction faster
   * and use less memory.
   */
  public static final int DEFAULT_MAX_EDGES_PER_CELL = 10;

  /**
   * Default maximum cell size, relative to an edge's length, for which that edge is considered
   * 'long'. Long edges are not counted towards {@link Options#maxEdgesPerCell}. The size and speed
   * of the index are typically not very sensitive to this parameter. Reasonable values range from
   * 0.1 to 10, with smaller values causing more aggressive subdivision of long edges grouped
   * closely together.
   */
  public static final double DEFAULT_CELL_SIZE_TO_LONG_EDGE_RATIO = 1.0;

  /**
   * Default minimum fraction of 'short' edges that must be present in a cell in order for it to be
   * subdivided. The default value of 0.2 was chosen to make index construction as fast as possible
   * while still protecting against possible quadratic space usage.
   */
  static final double DEFAULT_MIN_SHORT_EDGE_FRACTION = 0.2;

  /**
   * The current encoding version. When adding a new encoding, be aware that old binaries will not
   * be able to decode it.
   */
  public static final int CURRENT_ENCODING_VERSION = 0;

  /** The options supplied for this index. */
  protected final Options options;

  /** Shapes currently in the index. */
  protected List<S2Shape> shapes;

  /**
   * Essentially a map from each non-overlapping cell id to the shapes that intersect that cell,
   * clipped to include only the edges that intersect.
   *
   * <p>Note that this field is updated lazily. Accessing the iterator is the most common way to
   * construct the index.
   */
  private List<Cell> cells = ImmutableList.of();

  /** The index of the first shape that has been queued for insertion but not processed yet. */
  private int pendingInsertionsBegin = 0;

  /** The shapes that have been queued for removal but not processed yet (not yet used.) */
  private final List<S2Shape> pendingRemovals = Lists.newArrayList();

  /**
   * If true, the index is up to date. If false the index is updating or stale and requires an
   * update. This is marked volatile to avoid memory consistency errors with this field, which
   * allows us to avoid taking a lock when no update is required.
   */
  private volatile boolean isIndexFresh = true;

  /** Creates an S2ShapeIndex that uses the default options, {@link Options}. */
  @JsIgnore // Only one constructor is allowed for J2CL.
  public S2ShapeIndex() {
    this(new Options());
  }

  /** Creates an S2ShapeIndex with the given options. */
  @JsConstructor // JSInterop requires all other constructors to delegate to this one.
  public S2ShapeIndex(Options options) {
    this.options = options;
    shapes = new ArrayList<>();
  }

  /**
   * Creates an S2ShapeIndex that uses the default options, {@link Options}, and contains the
   * given shapes.
   */
  public static S2ShapeIndex fromShapes(S2Shape... shapes) {
    S2ShapeIndex index = new S2ShapeIndex();
    for (S2Shape shape : shapes) {
      index.add(shape);
    }
    return index;
  }

  /** Returns the options used for this index. */
  public Options options() {
    return options;
  }

  /**
   * Returns an immutable list view of shapes in the index. When shapes are added or removed, the
   * returned view is updated as well.
   */
  public List<S2Shape> getShapes() {
    return Collections.unmodifiableList(shapes);
  }

  /** Adds the given shape to this index. Invalidates all iterators and their associated data. */
  public void add(S2Shape shape) {
    // Insertions are processed lazily by applyUpdates(). All we do now is assign a unique id to the
    // shape, sequentially starting from 0 in the order shapes are inserted.
    shapes.add(shape);
    isIndexFresh = false;
  }

  // TODO(user): Implement remove() and incremental calls to add().

  /**
   * Currently not implemented. Will eventually remove the given shape from the index, and
   * invalidate all iterators and their associated data.
   *
   * @param shape the shape to remove
   */
  public void remove(S2Shape shape) {
    throw new UnsupportedOperationException("Not implemented yet");
  }

  /** Clears the contents of the index and resets it to its original state. */
  public void reset() {
    cells = ImmutableList.of();
    pendingRemovals.clear();
    shapes.clear();
    isIndexFresh = false;
    pendingInsertionsBegin = 0;
  }

  /**
   * Returns a new iterator over the cells of this index, positioned at the first cell in the index,
   * after initializing any pending updates.
   */
  @CanIgnoreReturnValue
  public S2Iterator<Cell> iterator() {
    // TODO(user): Make applyUpdates public, migrate clients who are currently using
    // iterator() for the side effect, and annotate iterator() with @CheckReturnValue.
    applyUpdates();
    return S2Iterator.create(cells);
  }

  /**
   * Returns true if there are no pending updates that need to be applied. This can be useful to
   * avoid building the index unnecessarily, or for choosing between two different algorithms
   * depending on whether the index is available.
   */
  public boolean isFresh() {
    return isIndexFresh;
  }

  /**
   * Ensures pending updates have been applied, returning immediately if the index is fresh as
   * reported by {@link #isFresh()}, and otherwise blocking while the index is built.
   *
   * <p>This operation is thread safe, guarded by 'this'.
   */
  void applyUpdates() {
    // The most common case is that the index is already fresh. So just return in that case without
    // taking a lock.
    if (isIndexFresh) {
      return;
    }

    // Otherwise an update is needed, so lock on 'this'. One of the contending threads will "win"
    // and update the index, and the others will just immediately exit the lock.
    synchronized (this) {
      if (!isIndexFresh) {
        // This thread won the race and must do the update.

        Preconditions.checkState(cells.isEmpty(), "Incremental updates not supported yet");

        int numEdges = 0;
        for (int i = pendingInsertionsBegin; i < shapes.size(); i++) {
          numEdges += shapes.get(i).numEdges();
        }

        // As a a pessimistic overestimate, assume we will have about 50% of the edges cross into
        // other cells, and cells end up about 50% full.
        cells = createList(3 * numEdges / options.maxEdgesPerCell / 4);

        // Create a list to hold edges that intersect each face, assuming the worse case scenario of
        // every edge intersecting every face. By far the most common case is that all edges
        // intersect one face, and the other lists are unused. The createList() method thus attempts
        // to avoid overhead for lists until they are actively used; i.e. numEdges edges are not
        // allocated right here, this is just advice on what approach to use.
        List<List<FaceEdge>> allEdges = createList(6);
        for (int face = 0; face < 6; face++) {
          List<FaceEdge> edges = createList(numEdges);
          allEdges.add(edges);
        }

        InteriorTracker tracker = new InteriorTracker(shapes.size() - pendingInsertionsBegin);
        for (int i = pendingInsertionsBegin; i < shapes.size(); i++) {
          addShapeEdges(i, allEdges, tracker);
        }
        for (int face = 0; face < 6; face++) {
          updateFaceEdges(face, allEdges.get(face), tracker);
          // Save memory by clearing each set of face edges after we are done with them.
          allEdges.set(face, null);
        }
        pendingInsertionsBegin = shapes.size();
        isIndexFresh = true;
      }
    }
  }

  /**
   * Reserves an appropriate amount of space for the top-level face edges. These lists are
   * responsible for most of the temporary memory usage during index construction. Furthermore, if
   * the arrays are grown via add() then a large fraction of the total run time consists of copying
   * data as these arrays grow, or handling GC events since reallocating a large array generates
   * objects that are difficult for the GC to reap cleanly without a stop-the-world collection.
   *
   * <p>However, creating the lists of edges for each face with this method is about 1% slower than
   * just using a List with smooth capacity increases. See {@link ShardedList} for details on the
   * implementation, but note also that allocating single large lists is hard on the Java garbage
   * collection environment, so not only is this method slower in the absence of GC measurements, it
   * is significantly faster when you also consider the impact of very large array allocations on
   * the garbage collector.
   */
  @SuppressWarnings("unused")
  private void reserveSpace(int numEdges, List<List<FaceEdge>> allEdges) {
    // TODO(user): reserveSpace is currently unused.

    // If the number of edges is relatively small, then the fastest approach is to simply reserve
    // space on every face for the maximum possible number of edges.
    final int maxCheapMemoryBytes = 10 << 20; // 10MB
    final int faceEdgeSize = 140; // bytes
    final int maxCheapNumEdges = maxCheapMemoryBytes / (6 * faceEdgeSize);
    if (numEdges <= maxCheapNumEdges) {
      for (int face = 0; face < 6; face++) {
        List<FaceEdge> edges = new SimpleList<>(numEdges);
        allEdges.add(edges);
      }
      return;
    }

    // Otherwise we estimate the number of edges on each face by taking a random sample. The goal is
    // to come up with an estimate that is fast and accurate for non-pathological geometry. If our
    // estimates happen to be wrong, the list will still grow automatically - the main side effects
    // are that memory usage will be larger (by up to a factor of 3), and constructing the index
    // will be about 10% slower.
    //
    // Given a desired sample size, we choose equally spaced edges from throughout the entire data
    // set. We use a Bresenham-type algorithm to choose the samples.
    final int desiredSampleSize = 10000;
    final int sampleInterval = max(1, numEdges / desiredSampleSize);

    // Initialize "edgeId" to be midway through the first sample interval. Because samples are
    // equally spaced the actual sample size may differ slightly from the desired sample size.
    int edgeId = sampleInterval / 2;
    MutableEdge edge = new MutableEdge();
    final int actualSampleSize = (numEdges + edgeId) / sampleInterval;
    int[] faceCount = {0, 0, 0, 0, 0, 0};
    for (int i = pendingInsertionsBegin; i < shapes.size(); i++) {
      S2Shape shape = shapes.get(i);
      edgeId += shape.numEdges();
      while (edgeId >= sampleInterval) {
        edgeId -= sampleInterval;
        shape.getEdge(edgeId, edge);
        // For speed, we only count the face containing one endpoint of the edge. In general the
        // edge could span all 6 faces (with padding), but it's not worth the expense to compute
        // this more accurately.
        faceCount[S2Projections.xyzToFace(edge.a)]++;
      }
    }

    // Now given the raw face counts, compute a confidence interval such that we will be unlikely to
    // allocate too little space. Computing accurate binomial confidence intervals is expensive and
    // not really necessary. Instead we use a simple approximation:
    //
    //  - For any face with at least 1 sample, we use at least a 4-sigma confidence interval. (The
    //    chosen width is adequate for the worst case accuracy, which occurs when the face contains
    //    approximately 50% of the edges.)  Assuming that our sample is representative, the
    //    probability of reserving too little space is approximately 1 in 30,000.
    //  - For faces with no samples at all, we don't bother reserving space. It is quite likely that
    //    such faces are truly empty, so we save time and memory this way. If the face does contain
    //    some edges, there will only be a few so it is fine to let the list grow automatically.
    //
    // On average, we reserve 2% extra space for each face that has geometry.

    // maxSemiWidth is the maximum semi-width over all probabilities p of a 4-sigma binomial
    // confidence interval with a sample size of 10,000.
    final double maxSemiWidth = 0.02;
    final double sampleRatio = 1.0 / actualSampleSize;
    for (int face = 0; face < 6; face++) {
      List<FaceEdge> edges;
      if (faceCount[face] > 0) {
        int fraction = (int) (sampleRatio * faceCount[face] + maxSemiWidth);
        edges = new SimpleList<>(1 + fraction * numEdges);
      } else {
        edges = new SimpleList<>(1);
      }
      allEdges.add(edges);
    }
  }

  /**
   * Clips all edges of the given shape to the six cube faces, and adds the clipped edges to {@code
   * allEdges}.
   */
  private void addShapeEdges(int shapeId, List<List<FaceEdge>> allEdges, InteriorTracker tracker) {
    S2Shape shape = shapes.get(shapeId);
    boolean hasInterior = shape.hasInterior();
    if (hasInterior) {
      tracker.addShape(shapeId, shape);
    }
    int numEdges = shape.numEdges();
    MutableEdge edge = new MutableEdge();
    R2Vector a = new R2Vector();
    R2Vector b = new R2Vector();
    double ratio = options.getCellSizeToLongEdgeRatio();

    for (int e = 0; e < numEdges; e++) {
      shape.getEdge(e, edge);
      if (hasInterior) {
        tracker.testEdge(shapeId, edge.a, edge.b);
      }

      // Fast path: both endpoints are on the same face, and are far enough from the edge of the
      // face that don't intersect any (padded) adjacent face.
      int aFace = S2Projections.xyzToFace(edge.a);
      if (aFace == S2Projections.xyzToFace(edge.b)) {
        S2Projections.validFaceXyzToUv(aFace, edge.a, a);
        S2Projections.validFaceXyzToUv(aFace, edge.b, b);
        final double kMaxUV = 1 - CELL_PADDING;
        if (abs(a.x) <= kMaxUV
            && abs(a.y) <= kMaxUV
            && abs(b.x) <= kMaxUV
            && abs(b.y) <= kMaxUV) {
          allEdges.get(aFace).add(new FaceEdge(shapeId, e, edge.a, edge.b, a, b, ratio));
          continue;
        }
      }

      // Otherwise we simply clip the edge to all six faces.
      for (int face = 0; face < 6; face++) {
        if (S2EdgeUtil.clipToPaddedFace(edge.a, edge.b, face, CELL_PADDING, a, b)) {
          allEdges.get(face).add(new FaceEdge(shapeId, e, edge.a, edge.b, a, b, ratio));
        }
      }
    }
  }

  /**
   * Given a face and a list of edges that intersect that face, insert or remove all the edges from
   * the index. (An edge is inserted if shape(id) is not null, and removed otherwise.)
   */
  private void updateFaceEdges(int face, List<FaceEdge> faceEdges, InteriorTracker tracker) {
    int numEdges = faceEdges.size();
    if (numEdges == 0 && tracker.focusCount == 0) {
      return;
    }

    // Create the initial ClippedEdge for each FaceEdge. Additional clipped edges are created when
    // edges are split between child cells.
    List<ClippedEdge> clippedEdges = createList(numEdges);
    R2Rect bound = R2Rect.empty();
    for (int i = 0; i < numEdges; i++) {
      FaceEdge edge = faceEdges.get(i);
      ClippedEdge clipped = new ClippedEdge();
      clipped.orig = edge;
      clipped.bound.x().initFromPointPair(edge.ax, edge.bx);
      clipped.bound.y().initFromPointPair(edge.ay, edge.by);
      clippedEdges.add(clipped);
      bound.addRect(clipped.bound);
    }

    // Construct the initial face cell containing all the edges, and then update all the edges in
    // the index recursively.
    S2CellId faceId = S2CellId.fromFace(face);
    S2PaddedCell pcell = new S2PaddedCell(faceId, CELL_PADDING);
    EdgeAllocator alloc = new EdgeAllocator(clippedEdges.size());
    if (numEdges > 0) {
      S2CellId shrunkId = pcell.shrinkToFit(bound);
      if (shrunkId.id() != pcell.id().id()) {
        // All the edges are contained by some descendant of the face cell. We can save a lot of
        // work by starting directly with that cell, but if we are in the interior of at least one
        // shape then we need to create index entries for the cells we are skipping over.
        skipCellRange(faceId.rangeMin(), shrunkId.rangeMin(), tracker, alloc);
        pcell = new S2PaddedCell(shrunkId, CELL_PADDING);
        updateEdges(pcell, clippedEdges, tracker, alloc);
        skipCellRange(shrunkId.rangeMax().next(), faceId.rangeMax().next(), tracker, alloc);
        return;
      }
    }
    // Otherwise (no edges, or no shrinking is possible), subdivide normally.
    updateEdges(pcell, clippedEdges, tracker, alloc);
  }

  /**
   * Skips over the cells in the given range, creating index cells if we are currently in the
   * interior of at least one shape.
   */
  private void skipCellRange(
      S2CellId begin, S2CellId end, InteriorTracker tracker, EdgeAllocator alloc) {
    if (tracker.focusCount > 0) {
      // If we are in the interior of at least one shape, then generate the list of cell ids that we
      // need to visit, and create an index entry with no edges, for each one.
      S2CellUnion skipped = new S2CellUnion();
      skipped.initFromBeginEnd(begin, end);
      List<ClippedEdge> clippedEdges = ImmutableList.of();
      for (int i = 0; i < skipped.size(); i++) {
        S2PaddedCell pcell = new S2PaddedCell(skipped.cellId(i), CELL_PADDING);
        updateEdges(pcell, clippedEdges, tracker, alloc);
      }
    }
  }

  /**
   * Given a cell and a set of ClippedEdges whose bounding boxes intersect that cell, insert or
   * remove all the edges from the index. Temporary space for edges that need to be subdivided is
   * allocated from the given EdgeAllocator.
   */
  boolean makeIndexCell(
      S2PaddedCell pcell, List<ClippedEdge> edges, InteriorTracker tracker) {
    if (edges.isEmpty() && tracker.focusCount == 0) {
      // No index cell is needed. In most cases this situation is detected before we get to this
      // point, but this can happen when all shapes in a cell are removed.
      return true;
    }

    // We can show using amortized analysis that the total index size is
    //
    //     O(c1 * n + c2 * (1 - f) / f * n)
    //
    // where n is the number of input edges (and where we also count an "edge" for each shape with
    // an interior but no edges), f is the value of Options.minShortEdgeFraction, and c1 and c2 are
    // constants where c2 is about 20 times larger than c1.
    //
    // First observe that the space used by an S2ShapeIndex is proportional to the space used by all
    // of its index cells, and the space used by an S2ShapeIndex.Cell is proportional to the number
    // of edges that intersect that cell plus the number of shapes that contain the entire cell
    // ("containing shapes"). Define an "index entry" as an intersecting edge or containing shape
    // stored by an index cell. Our goal is then to bound the number of index entries.
    //
    // We divide the index entries into two groups. An index entry is "short" if it represents an
    // edge that was considered short in that index cell's parent, and "long" otherwise. (Note that
    // the long index entries also include the containing shapes mentioned above.) We then bound
    // the maximum number of both types of index entries by associating them with edges that were
    // considered short in those index cells' parents.
    //
    // First consider the short index entries for a given edge E. Let S be the set of index cells
    // that intersect E and where E was considered short in those index cells' parents. Since E was
    // short in each parent cell, the width of those parent cells is at least some fraction "g" of
    // E's length (as controlled by FLAGS_s2shape_index_cell_size_to_long_edge_ratio). Therefore
    // the minimum width of each cell in S is also at least some fraction of E's length (i.e.,
    // g / 2). This implies that there are at most a constant number c1 of such cells, since they
    // all intersect E and do not overlap, which means that there are at most (c1 * n) short entries
    // in total.
    //
    // With cellSizeToLongEdgeRatio = 1.0 (the default value), it can be shown that c1 = 10. In
    // other words, it is not possible for a given edge to intersect more than 10 index cells where
    // it was considered short in those cells' parents. The value of c1 can be reduced as low as 4
    // by increasing cellSizeToLongEdgeRatio to about 3.1. (The reason the minimum value is 3.1
    // rather than 2.0 is that this ratio is defined in terms of the average edge length of cells at
    // a given level, rather than their minimum width, and 2 * (S2Projections.PROJ.avgEdge /
    // S2Projections.PROJ.minWidth) ~= 3.1.)
    //
    // Next we consider the long index entries. Let c2 be the maximum number of index cells where a
    // given edge E was considered short in those cells' parents. (Unlike the case above, we do not
    // require that these cells intersect E.)  Because the minimum width of each parent cell is at
    // least some fraction of E's length and the parent cells at a given level do not overlap, there
    // can be at most a small constant number of index cells at each level where E is considered
    // short in those cells' parents. For example, consider a very short edge E that intersects the
    // midpoint of a cell edge at level 0. There are 16 cells at level 30 where E was considered
    // short in the parent cell, 12 cells at each of levels 29..2, and 4 cells at levels 1 and 0
    // (pretending that all 6 face cells share a common "parent"). This yields a total of c2 = 360
    // index cells. This is actually the worst case for cellSizeToLongEdgeRatio >= 3.1; with the
    // default value of 1.0 it is possible to have a few more index cells at levels 29 and 30, for a
    // maximum of c2 = 366 index cells.
    //
    // The code below subdivides a given cell only if
    //
    //     s > f * (s + l)
    //
    // where "f" is the minShortEdgeFraction option, "s" is the number of short edges that intersect
    // the cell, and "l" is the number of long edges that intersect the cell plus an upper bound on
    // the number of shapes that contain the entire cell. (It is an upper bound rather than an
    // exact count because we use the number of shapes that contain an arbitrary vertex of the
    // cell.) Note that the number of long index entries in each child of this cell is at most "l"
    // because no child intersects more edges than its parent or is entirely contained by more
    // shapes than its parent.
    //
    // The inequality above can be rearranged to give
    //
    //    l < s * (1 - f) / f
    //
    // This says that each long index entry in a child cell can be associated with at most (1 - f) /
    // f edges that were considered short when the parent cell was subdivided. Furthermore we know
    // that there are at most c2 index cells where a given edge was considered short in the parent
    // cell. Since there are only n edges in total, this means that the maximum number of long
    // index entries is at most
    //
    //    c2 * (1 - f) / f * n
    //
    // and putting this together with the result for short index entries gives the desired bound.
    //
    // There are a variety of ways to make this bound tighter, e.g. when "n" is relatively small.
    // For example when the indexed geometry satisfies the requirements of S2BooleanOperation (i.e.,
    // shape interiors are disjoint) and the minShortEdgeFraction option is not too large, then the
    // constant c2 above is only about half as big (i.e., c2 ~= 180). This is because the worst
    // case under these circumstances requires having many shapes whose interiors overlap.

    // Continue subdividing if the proposed index cell would contain too many edges that are "short"
    // relative to its size (as controlled by the cellSizeToLongEdgeRatio option). Usually "too
    // many" means more than Options.maxEdgesPerCell(), but this value might be increased if the
    // cell has a lot of long edges and/or containing shapes. This strategy ensures that the total
    // index size is linear (see above).
    if (edges.size() > options.maxEdgesPerCell) {
      int maxShortEdges = max(
          options.maxEdgesPerCell,
          (int) (options.minShortEdgeFraction * (edges.size() + tracker.focusCount)));
      int count = 0;
      for (ClippedEdge edge : edges) {
        count += pcell.level() < edge.orig.maxLevel ? 1 : 0;
        if (count > maxShortEdges) {
          return false;
        }
      }
    }

    // There aren't too many edges, so build a leaf cell from these edges.

    // We update the InteriorTracker as follows. For every S2Cell in the index we construct two
    // edges: one edge from entry vertex of the cell to its center, and one from the cell center to
    // its exit vertex. Here "entry" and "exit" refer to the S2CellId ordering, i.e. the order in
    // which points are encountered along the S2 space-filling curve. The exit vertex then becomes
    // the entry vertex for the next cell in the index, unless there are one or more empty
    // intervening cells, in which case the InteriorTracker state is unchanged because the
    // intervening cells have no edges.

    // Shift the InteriorTracker focus point to the center of the current cell.
    int numEdges = edges.size();
    if (tracker.isActive() && numEdges > 0) {
      if (!tracker.atCellId(pcell.id())) {
        tracker.moveTo(pcell.getEntryVertex());
      }
      tracker.drawTo(pcell.getCenter());
      for (int i = 0; i < numEdges; i++) {
        ClippedEdge edge = edges.get(i);
        FaceEdge orig = edge.orig;
        if (shapes.get(orig.shapeId).hasInterior()) {
          tracker.testEdge(orig.shapeId, orig.va, orig.vb);
        }
      }
    }

    // Allocate and fill a new index cell. To get the total number of shapes we need to merge the
    // shapes associated with the intersecting edges together with the shapes that happen to contain
    // the cell center.

    // The first clipped shape will contain this cell. After its used this value will be null, so
    // we don't waste space storing the cell ID repeatedly.
    S2CellId cellId = pcell.id();

    // To fill the index cell we merge the two sources of shapes: "edge shapes" (those that have at
    // least one edge that intersects this cell), and "containing shapes" (those that contain the
    // cell center.) We keep track of the index of the next intersecting edge and the next
    // containing shape as we go along.
    int numShapes = 0;
    int edgesIndex = 0;
    int trackerIndex = 0;
    int nextShapeId = shapes.size();
    while (edgesIndex < numEdges || trackerIndex < tracker.focusCount) {
      int edgeId;
      if (edgesIndex < numEdges) {
        edgeId = edges.get(edgesIndex).orig.shapeId;
      } else {
        edgeId = nextShapeId;
      }
      int trackerId;
      if (trackerIndex < tracker.focusCount) {
        trackerId = tracker.focusedShapes[trackerIndex];
      } else {
        trackerId = nextShapeId;
      }
      S2ClippedShape clipped;
      if (trackerId < edgeId) {
        // The entire cell is in the shape interior.
        clipped = S2ClippedShape.Contained.create(cellId, shapes.get(trackerId));
        cellId = null;
        trackerIndex++;
      } else {
        // Count the number of edges for this shape and allocate space for them.
        int firstEdge = edgesIndex;
        while (edgesIndex < numEdges && edges.get(edgesIndex).orig.shapeId == edgeId) {
          edgesIndex++;
        }
        boolean containsCenter = trackerId == edgeId;
        clipped =
            S2ClippedShape.create(
                cellId, shapes.get(edgeId), containsCenter, edges, firstEdge, edgesIndex);
        cellId = null;
        if (containsCenter) {
          trackerIndex++;
        }
      }
      tracker.tempClippedShapes[numShapes++] = clipped;
    }

    // updateEdges() visits cells in increasing order of S2CellId, so during initial construction of
    // the index we can just append new cells at the end of the list. This is much faster than
    // sorting the cells afterward.
    cells.add(Cell.create(numShapes, tracker.tempClippedShapes));

    // Shift the InteriorTracker focus point to the exit vertex of this cell.
    if (tracker.isActive() && !edges.isEmpty()) {
      tracker.drawTo(pcell.getExitVertex());
      for (int i = 0; i < numEdges; i++) {
        ClippedEdge edge = edges.get(i);
        FaceEdge orig = edge.orig;
        if (shapes.get(orig.shapeId).hasInterior()) {
          tracker.testEdge(orig.shapeId, orig.va, orig.vb);
        }
      }
      tracker.doneCellId(pcell.id());
    }

    return true;
  }

  /**
   * Given a cell and a set of ClippedEdges whose bounding boxes intersect that cell, add all the
   * edges from the index. Temporary space for edges that need to be subdivided is allocated from
   * the given EdgeAllocator.
   */
  private void updateEdges(
      S2PaddedCell pcell, List<ClippedEdge> edges, InteriorTracker tracker, EdgeAllocator alloc) {
    // TODO(user): Implement the "disjointFromIndex" optimization hint indicating that
    // cells does not contain any entries that overlap the given cell.

    // Cases where an index cell is not needed should be detected before this.
    assert !edges.isEmpty() || tracker.focusCount > 0;

    // This function is recursive with a maximum recursion depth of 30 (S2CellId.MAX_LEVEL).

    if (makeIndexCell(pcell, edges, tracker)) {
      // Skip splitting if we made a cell for 'edges'.
      return;
    }

    // Reserve space for the edges that will be passed to each child. We select the kind of list to
    // use based on how large the child edges lists could possibly be.
    int numEdges = edges.size();
    List<ClippedEdge> edges00 = createList(numEdges);
    List<ClippedEdge> edges01 = createList(numEdges);
    List<ClippedEdge> edges10 = createList(numEdges);
    List<ClippedEdge> edges11 = createList(numEdges);
    List<List<ClippedEdge>> ijEdges = ImmutableList.of(edges00, edges01, edges10, edges11);

    // Remember the current size of the EdgeAllocator so that we can free any edges that are
    // allocated during edge splitting.
    int allocSize = alloc.size();

    // Compute the middle of the padded cell, defined as the rectangle in (u,v)-space that belongs
    // to all four (padded) children. By comparing against the four boundaries of "middle" we can
    // determine which children each edge needs to be propagated to.
    R2Rect middle = pcell.middle();

    // Build up a list edges to be passed to each child cell. The (i,j) directions are left (i=0),
    // right (i=1), lower (j=0), and upper (j=1). Note that the vast majority of edges are
    // propagated to a single child.
    for (int i = 0; i < numEdges; i++) {
      ClippedEdge edge = edges.get(i);
      if (edge.bound.x().hi() <= middle.x().lo()) {
        // Edge is entirely contained in the two left children.
        clipVAxis(edge, middle.y(), edges00, edges01, alloc);
      } else if (edge.bound.x().lo() >= middle.x().hi()) {
        // Edge is entirely contained in the two right children.
        clipVAxis(edge, middle.y(), edges10, edges11, alloc);
      } else if (edge.bound.y().hi() <= middle.y().lo()) {
        // Edge is entirely contained in the two lower children.
        edges00.add(clipUBound(edge, true, middle.x().hi(), alloc));
        edges10.add(clipUBound(edge, false, middle.x().lo(), alloc));
      } else if (edge.bound.y().lo() >= middle.y().hi()) {
        // Edge is entirely contained in the two upper children.
        edges01.add(clipUBound(edge, true, middle.x().hi(), alloc));
        edges11.add(clipUBound(edge, false, middle.x().lo(), alloc));
      } else {
        // The edge bound spans all four children. The edge itself intersects either three or four
        // (padded) children.
        ClippedEdge left = clipUBound(edge, true, middle.x().hi(), alloc);
        clipVAxis(left, middle.y(), edges00, edges01, alloc);
        ClippedEdge right = clipUBound(edge, false, middle.x().lo(), alloc);
        clipVAxis(right, middle.y(), edges10, edges11, alloc);
      }
    }

    // Now recursively update the edges in each child. We call the children in increasing order of
    // S2CellId so that when the index is first constructed, all insertions into 'cells' are at the
    // end (which is much faster than sorting by cell id afterward.)
    for (int pos = 0; pos < 4; pos++) {
      List<ClippedEdge> childEdges = ijEdges.get(S2.posToIJ(pcell.orientation(), pos));
      if (!childEdges.isEmpty() || tracker.focusCount > 0) {
        S2PaddedCell childCell = pcell.childAtPos(pos);
        updateEdges(childCell, childEdges, tracker, alloc);
      }
    }
    alloc.reset(allocSize);
  }

  /**
   * Given an edge and two bound endpoints that need to be updated, allocates and returns a new edge
   * with the updated bound.
   */
  private static ClippedEdge updateBound(
      ClippedEdge edge, boolean uEnd, double u, boolean vEnd, double v, EdgeAllocator alloc) {
    ClippedEdge clipped = alloc.create();
    clipped.orig = edge.orig;
    if (uEnd) {
      clipped.bound.x().set(edge.bound.x().lo(), u);
    } else {
      clipped.bound.x().set(u, edge.bound.x().hi());
    }
    if (vEnd) {
      clipped.bound.y().set(edge.bound.y().lo(), v);
    } else {
      clipped.bound.y().set(v, edge.bound.y().hi());
    }
    // assert !clipped.bound.isEmpty();
    // assert edge.bound.contains(clipped.bound);
    return clipped;
  }

  private static ClippedEdge clipUBound(
      ClippedEdge edge, boolean uEnd, double u, EdgeAllocator alloc) {
    // First check whether the edge actually requires any clipping. (Sometimes this method is
    // called when clipping is not necessary, e.g. when one edge endpoint is in the overlap area
    // between two padded child cells.)
    if (!uEnd) {
      if (edge.bound.x().lo() >= u) {
        return edge;
      }
    } else {
      if (edge.bound.x().hi() <= u) {
        return edge;
      }
    }

    // We interpolate the new v-value from the endpoints of the original edge. This has two
    // advantages:
    //
    // (1) We don't need to store the clipped endpoints at all, just their bounding box; and
    // (2) It avoids the accumulation of roundoff errors due to repeated interpolations.
    //
    // The result needs to be clamped to ensure that it is in the appropriate range.
    FaceEdge e = edge.orig;
    double v = edge.bound.y().clampPoint(S2EdgeUtil.interpolateDouble(u, e.ax, e.bx, e.ay, e.by));

    // Determine which endpoint of the v-axis bound to update. If the edge slope is positive we
    // update the same endpoint, otherwise we update the opposite endpoint.
    boolean vEnd = ((e.ax > e.bx) != (e.ay > e.by)) ^ uEnd;
    return updateBound(edge, uEnd, u, vEnd, v, alloc);
  }

  private static ClippedEdge clipVBound(
      ClippedEdge edge, boolean vEnd, double v, EdgeAllocator alloc) {
    // See comments in clipUBound.
    if (!vEnd) {
      if (edge.bound.y().lo() >= v) {
        return edge;
      }
    } else {
      if (edge.bound.y().hi() <= v) {
        return edge;
      }
    }
    FaceEdge e = edge.orig;
    double u = edge.bound.x().clampPoint(S2EdgeUtil.interpolateDouble(v, e.ay, e.by, e.ax, e.bx));
    boolean uEnd = ((e.ax > e.bx) != (e.ay > e.by)) ^ vEnd;
    return updateBound(edge, uEnd, u, vEnd, v, alloc);
  }

  private static void clipVAxis(
      ClippedEdge edge,
      R1Interval middle,
      List<ClippedEdge> edges0,
      List<ClippedEdge> edges1,
      EdgeAllocator alloc) {
    if (edge.bound.y().hi() <= middle.lo()) {
      // Edge is entirely contained in the lower child.
      edges0.add(edge);
    } else if (edge.bound.y().lo() >= middle.hi()) {
      // Edge is entirely contained in the upper child.
      edges1.add(edge);
    } else {
      // The edge bound spans both children.
      edges0.add(clipVBound(edge, true, middle.hi(), alloc));
      edges1.add(clipVBound(edge, false, middle.lo(), alloc));
    }
  }

  /** Options that affect construction of the S2ShapeIndex. */
  @JsType
  public static class Options implements Serializable {
    private static final long serialVersionUID = 1L;

    private int maxEdgesPerCell = DEFAULT_MAX_EDGES_PER_CELL;
    private double cellSizeToLongEdgeRatio = DEFAULT_CELL_SIZE_TO_LONG_EDGE_RATIO;
    private double minShortEdgeFraction = DEFAULT_MIN_SHORT_EDGE_FRACTION;

    /**
     * Returns the maximum number of edges per cell (default 10.) If a cell has more than this many
     * edges that are "long" relative to the cell size, and it is not a leaf cell, then it is
     * subdivided. Whether an edge is considered "long" is controlled by the value returned by
     * {@link #getCellSizeToLongEdgeRatio()}.
     */
    public int getMaxEdgesPerCell() {
      return maxEdgesPerCell;
    }

    /**
     * Sets the new number of max edges per cell. Only has an effect during index construction,
     * usually triggered by calling {@link S2ShapeIndex#iterator()}.
     */
    public void setMaxEdgesPerCell(int maxEdgesPerCell) {
      this.maxEdgesPerCell = maxEdgesPerCell;
    }

    /**
     * Returns the cell size relative to the length of an edge at which it is first considered to be
     * "long" (default is 1.0). Long edges do not contribute toward the decision to subdivide a cell
     * further. The size and speed of the index are typically not very sensitive to this parameter.
     * Reasonable values range from 0.1 to 10, with smaller values causing more aggressive
     * subdivision of long edges grouped closely together. For example, a value of 2.0 means that
     * the cell must be at least twice the size of the edge in order for that edge to be counted.
     * There are two reasons for not counting long edges:
     *
     * <ol>
     *   <li>Such edges typically need to be propagated to several children, which increases time
     *       and memory costs without much benefit, and
     *   <li>In pathological cases, many long edges close together could force subdivision to
     *       continue all the way to the leaf cell level.
     * </ol>
     */
    public double getCellSizeToLongEdgeRatio() {
      return cellSizeToLongEdgeRatio;
    }

    /**
     * Sets the new ratio of cell size to long edges. Only has an effect during index construction,
     * usually triggered by calling {@link S2ShapeIndex#iterator()}.
     */
    public void setCellSizeToLongEdgeRatio(double cellSizeToLongEdgeRatio) {
      this.cellSizeToLongEdgeRatio = cellSizeToLongEdgeRatio;
    }

    /**
     * Returns the minimum fraction of 'short' edges that must be present in a cell in order for it
     * to be subdivided.
     */
    public double getMinShortEdgeFraction() {
      return minShortEdgeFraction;
    }

    /**
     * Sets the new minimum fraction of 'short' edges that must be present in a cell in order for it
     * to be subdivided. Only has an effect during index construction. If this option is non-zero
     * then the total index size and construction time are guaranteed to be linear in the number of
     * input edges; this prevents the worst-case quadratic space and time usage that can otherwise
     * occur with certain input configurations. Specifically, the maximum index size is:
     *
     * <pre>O((c1 + c2 * (1 - f) / f) * n)</pre>
     *
     * <p>where n is the number of input edges, f is this parameter value, and constant c2 is
     * roughly 20 times larger than constant c1. (The exact values of c1 and c2 depend on the
     * cellSizeToLongEdgeRatio and maxEdgesPerCell options and certain properties of the input
     * geometry such as whether it consists of O(1) shapes, whether it includes polygons, and
     * whether the polygon interiors are disjoint.)
     *
     * <p>Reasonable values range from 0.1 to perhaps 0.95. The main factors to consider when
     * choosing this parameter are:
     *
     * <ul>
     *   <li>For pathological geometry, larger values result in indexes that are smaller and faster
     *       to construct but have worse query performance, due to having more edges per cell.
     *       However, note that even a setting of 0.1 reduces the worst case by 100x compared with a
     *       setting of 0.001.
     *   <li>For normal geometry, values up to about 0.8 result in indexes that are virtually
     *       unchanged except for a slight increase in index construction time, proportional to the
     *       parameter value f, for very large inputs. With millions of edges, indexing time
     *       increases by about (15% * f), e.g. a parameter value of 0.5 slows down indexing for
     *       very large inputs by about 7.5%. Indexing time for small inputs is unaffected.
     *   <li>Values larger than about 0.8 start to affect index construction even for normal
     *       geometry, resulting in smaller indexes and faster construction times but gradually
     *       worse query performance.
     * </ul>
     *
     * Essentially this parameter provides control over a space-time trade-off that largely affects
     * only pathological geometry.
     */
    public void setMinShortEdgeFraction(double minShortEdgeFraction) {
      this.minShortEdgeFraction = minShortEdgeFraction;
    }
  }

  /**
   * This class contains the set of clipped shapes within a particular index cell, sorted in
   * increasing order of shape id.
   *
   * <p>To be as memory efficient as possible, we specialize two very common cases.
   *
   * <ul>
   *   <li>The Cell class is extended by S2ClippedShape, and in the *very* common case of a cell
   *       containing just one clipped shape, we return the shape directly without wrapping it (this
   *       requires that the clipped shapes contain the cell IDs, rather than the Cell; more about
   *       that on {@link S2ClippedShape}.)
   *   <li>In the fairly common case of a cell intersecting two shapes, we have a BinaryCell
   *       implementation that is half the size of the general purpose MultiCell in that case.
   * </ul>
   */
  @JsType
  public abstract static class Cell implements S2Iterator.Entry, Serializable {
    private static final long serialVersionUID = 1L;

    /** Returns a Cell with a copy of the given shapes, specialized for the number of elements. */
    static Cell create(int size, S2ClippedShape[] tempClippedShapes) {
      switch (size) {
        case 1:
          return tempClippedShapes[0];
        case 2:
          return new BinaryCell(tempClippedShapes[0], tempClippedShapes[1]);
        default:
          return new MultiCell(Arrays.copyOf(tempClippedShapes, size));
      }
    }

    @Override
    public long id() {
      // We store the cell ID on the first clipped shape.
      return clipped(0).id();
    }

    /** Returns the number of clipped shapes in this cell. */
    public abstract int numShapes();

    /** Returns the number of clipped edges in this cell (shape edges that intersect the cell). */
    public abstract int numEdges();

    /**
     * Returns the clipped shape at the given index.
     *
     * @param i must be at least 0 and less than {@link #numShapes()}
     */
    public abstract S2ClippedShape clipped(int i);

    /**
     * Returns the clipped shape corresponding to the given shape, or null if the shape does not
     * intersect this cell.
     */
    @Nullable
    S2ClippedShape findClipped(S2Shape shape) {
      // Linear search is fine because the number of shapes per cell is typically very small (most
      // often 1), and is large only for pathological inputs (e.g. very deeply nested loops).
      for (int i = 0; i < numShapes(); i++) {
        S2ClippedShape clipped = clipped(i);
        if (clipped.shape == shape) {
          return clipped;
        }
      }
      return null;
    }

    /** A specialization of Cell for the case of two clipped shapes. Also very common. */
    private static final class BinaryCell extends Cell {
      private static final long serialVersionUID = 1L;
      private final S2ClippedShape shape1;
      private final S2ClippedShape shape2;

      @JsConstructor
      BinaryCell(S2ClippedShape shape1, S2ClippedShape shape2) {
        this.shape1 = shape1;
        this.shape2 = shape2;
      }

      @Override
      public int numShapes() {
        return 2;
      }

      @Override
      public int numEdges() {
        return shape1.numEdges() + shape2.numEdges();
      }

      @Override
      public S2ClippedShape clipped(int i) {
        switch (i) {
          case 0:
            return shape1;
          case 1:
            return shape2;
          default:
            throw new ArrayIndexOutOfBoundsException();
        }
      }
    }

    /**
     * A specialization of Cell for multiple shapes per cell. Last resort, largest-memory
     * implementation that is not used very often.
     */
    private static final class MultiCell extends Cell {
      private final S2ClippedShape[] clippedShapes;

      @JsConstructor
      MultiCell(S2ClippedShape[] shapes) {
        this.clippedShapes = shapes;
      }

      @Override
      public int numShapes() {
        return clippedShapes.length;
      }

      @Override
      public int numEdges() {
        int sum = 0;
        for (S2ClippedShape shape : clippedShapes) {
          sum += shape.numEdges();
        }
        return sum;
      }

      @Override
      public S2ClippedShape clipped(int i) {
        return clippedShapes[i];
      }
    }
  }

  /**
   * The possible relationships between a "target" cell and the cells of the S2ShapeIndex. If the
   * target is an index cell or is contained by an index cell, it is "INDEXED". If the target is
   * subdivided into one or more index cells, it is "SUBDIVIDED". Otherwise it is "DISJOINT".
   */
  @JsEnum
  public enum CellRelation {
    /** Target is contained by an index cell. */
    INDEXED,
    /** Target is subdivided into one or more index cells. */
    SUBDIVIDED,
    /** Target does not intersect any index cells. */
    DISJOINT
  }

  /**
   * S2ClippedShape represents the part of a shape that intersects an S2Cell. It consists of the set
   * of edge ids that intersect that cell, and a boolean indicating whether the center of the cell
   * is inside the shape (for shapes that have an interior). Edges themselves are not clipped; we
   * always use the original edges for intersection tests so that the results will be the same as
   * the original shape.
   *
   * <p>Most of the index memory is consumed here, so we have several strategies to be as efficient
   * as possible:
   *
   * <ul>
   *   <li>Since only the first clipped shape in a cell needs to store the cell ID, we have
   *       subclasses for both cases based on whether a cell ID was provided.
   *   <li>Since no edges, one edge, and one dense range of edges are very common cases, we have
   *       more efficient subclasses for those situations, in addition to the general purpose
   *       implementation.
   * </ul>
   */
  public abstract static class S2ClippedShape extends Cell {
    static S2ClippedShape create(
        S2CellId cellId,
        S2Shape shape,
        boolean containsCenter,
        List<ClippedEdge> edges,
        int start,
        int end) {
      int numEdges = end - start;
      if (numEdges == 1) {
        return OneEdge.create(cellId, shape, containsCenter, edges.get(start));
      }
      int edge = edges.get(start).orig.edgeId;
      for (int i = 1; i < numEdges; i++) {
        if (edge + i != edges.get(start + i).orig.edgeId) {
          return ManyEdges.create(cellId, shape, containsCenter, edges, start, end);
        }
      }
      return EdgeRange.create(cellId, shape, containsCenter, edge, numEdges);
    }

    static S2ClippedShape create(
        @Nullable S2CellId cellId, S2Shape shape, boolean containsCenter, int[] edges) {
      return ManyEdges.create(cellId, shape, containsCenter, edges);
    }

    static S2ClippedShape create(
        @Nullable S2CellId cellId, S2Shape shape, boolean containsCenter, int offset, int count) {
      return EdgeRange.create(cellId, shape, containsCenter, offset, count);
    }

    /** The original shape that this S2ClippedShape represents. */
    private final S2Shape shape;

    @JsConstructor
    private S2ClippedShape(S2Shape shape) {
      this.shape = shape;
    }

    /** Returns the original shape this clipped shape was clipped from. */
    public final S2Shape shape() {
      return shape;
    }

    /**
     * Returns whether the center of the S2CellId is inside the shape, and always returns false for
     * shapes that do not have an interior according to {@link S2Shape#hasInterior()}.
     */
    public abstract boolean containsCenter();

    /** Returns the number of edges that intersect the S2CellId. */
    @Override
    public abstract int numEdges();

    /**
     * Returns the {@code i}<sup>th</sup> edge ID of this clipped shape. Edges are sorted in
     * increasing order of edge ID. The edge IDs may be passed to the corresponding shape's {@link
     * S2Shape#getEdge(int, com.google.common.geometry.S2Shape.MutableEdge)} method.
     *
     * @param i must be at least 0 and less than {@link #numEdges}
     */
    public abstract int edge(int i);

    /** Returns whether the clipped shape contains the given edge id. */
    public final boolean containsEdge(int edgeId) {
      // Linear search is fast because the number of edges per shape is typically very small (less
      // than 10).
      for (int e = 0; e < numEdges(); e++) {
        if (edge(e) == edgeId) {
          return true;
        }
      }
      return false;
    }

    /** For implementing the Cell interface, this class contains just 1 shape (itself.) */
    @Override
    public final int numShapes() {
      return 1;
    }

    /** For implementing the Cell interface, this class contains just 1 shape (itself.) */
    @Override
    public final S2ClippedShape clipped(int i) {
      // assert i == 0;
      return this;
    }

    /**
     * An S2ClippedShape for a shape that completely contains the cell (no edge intersections and
     * containsCenter is true.)
     */
    private abstract static class Contained extends S2ClippedShape {
      static Contained create(@Nullable S2CellId cellId, S2Shape shape) {
        if (cellId != null) {
          final long id = cellId.id();
          return new Contained(shape) {
            @Override
            public long id() {
              return id;
            }
          };
        } else {
          return new Contained(shape) {
            @Override
            public long id() {
              throw new UnsupportedOperationException();
            }
          };
        }
      }

      @JsConstructor
      private Contained(S2Shape shape) {
        super(shape);
      }

      @Override
      public final boolean containsCenter() {
        return true;
      }

      @Override
      public final int numEdges() {
        return 0;
      }

      @Override
      public final int edge(int i) {
        throw new ArrayIndexOutOfBoundsException();
      }
    }

    /** An S2ClippedShape that contains a single edge from a given shape. Very common. */
    private abstract static class OneEdge extends S2ClippedShape {
      static final OneEdge create(
          @Nullable S2CellId cellId,
          S2Shape shape,
          boolean containsCenter,
          ClippedEdge clippedEdge) {
        if (cellId != null) {
          final long id = cellId.id();
          if (containsCenter) {
            return new OneEdge(shape, clippedEdge) {
              @Override
              public long id() {
                return id;
              }

              @Override
              public boolean containsCenter() {
                return true;
              }
            };
          } else {
            return new OneEdge(shape, clippedEdge) {
              @Override
              public long id() {
                return id;
              }

              @Override
              public boolean containsCenter() {
                return false;
              }
            };
          }
        } else {
          if (containsCenter) {
            return new OneEdge(shape, clippedEdge) {
              @Override
              public long id() {
                throw new UnsupportedOperationException();
              }

              @Override
              public boolean containsCenter() {
                return true;
              }
            };
          } else {
            return new OneEdge(shape, clippedEdge) {
              @Override
              public long id() {
                throw new UnsupportedOperationException();
              }

              @Override
              public boolean containsCenter() {
                return false;
              }
            };
          }
        }
      }

      private final int edge;

      @JsConstructor
      private OneEdge(S2Shape shape, ClippedEdge clippedEdge) {
        super(shape);
        this.edge = clippedEdge.orig.edgeId;
      }

      @Override
      public final int numEdges() {
        return 1;
      }

      @Override
      public final int edge(int i) {
        return edge;
      }
    }

    /**
     * An S2ClippedShape that contains the non-contiguous edges from {@code start} to {@code end} in
     * {@code edges}. A much larger object than the other subclasses of S2ClippedShape, but also the
     * rarest.
     */
    private abstract static class ManyEdges extends S2ClippedShape {
      static ManyEdges create(
          @Nullable S2CellId cellId,
          S2Shape shape,
          boolean containsCenter,
          List<ClippedEdge> edges,
          int start,
          int end) {
        int[] edgeArray = new int[end - start];
        // Arrays.setAll(edgeArray, i -> edges.get(i + start).orig.edgeId); is Android-incompatible.
        for (int i = 0; i < edgeArray.length; i++) {
          edgeArray[i] = edges.get(i + start).orig.edgeId;
        }
        return create(cellId, shape, containsCenter, edgeArray);
      }

      static ManyEdges create(
          @Nullable S2CellId cellId, S2Shape shape, boolean containsCenter, int[] edges) {
        if (cellId != null) {
          final long id = cellId.id();
          if (containsCenter) {
            return new ManyEdges(shape, edges) {
              @Override
              public long id() {
                return id;
              }

              @Override
              public boolean containsCenter() {
                return true;
              }
            };
          } else {
            return new ManyEdges(shape, edges) {
              @Override
              public long id() {
                return id;
              }

              @Override
              public boolean containsCenter() {
                return false;
              }
            };
          }
        } else {
          if (containsCenter) {
            return new ManyEdges(shape, edges) {
              @Override
              public long id() {
                throw new UnsupportedOperationException();
              }

              @Override
              public boolean containsCenter() {
                return true;
              }
            };
          } else {
            return new ManyEdges(shape, edges) {
              @Override
              public long id() {
                throw new UnsupportedOperationException();
              }

              @Override
              public boolean containsCenter() {
                return false;
              }
            };
          }
        }
      }

      private final int[] edges;

      @JsConstructor
      private ManyEdges(S2Shape shape, int[] edges) {
        super(shape);
        this.edges = edges;
      }

      @Override
      public final int numEdges() {
        return edges.length;
      }

      @Override
      public final int edge(int i) {
        return edges[i];
      }
    }

    /** An S2ClippedShape containing a single range of contiguous edge IDs. Very common. */
    private abstract static class EdgeRange extends S2ClippedShape {
      static EdgeRange create(
          @Nullable S2CellId cellId, S2Shape shape, boolean containsCenter, int offset, int count) {
        if (cellId != null) {
          final long id = cellId.id();
          if (containsCenter) {
            return new EdgeRange(shape, offset, count) {
              @Override
              public long id() {
                return id;
              }

              @Override
              public boolean containsCenter() {
                return true;
              }
            };
          } else {
            return new EdgeRange(shape, offset, count) {
              @Override
              public long id() {
                return id;
              }

              @Override
              public boolean containsCenter() {
                return false;
              }
            };
          }
        } else {
          if (containsCenter) {
            return new EdgeRange(shape, offset, count) {
              @Override
              public long id() {
                throw new UnsupportedOperationException();
              }

              @Override
              public boolean containsCenter() {
                return true;
              }
            };
          } else {
            return new EdgeRange(shape, offset, count) {
              @Override
              public long id() {
                throw new UnsupportedOperationException();
              }

              @Override
              public boolean containsCenter() {
                return false;
              }
            };
          }
        }
      }

      private final int offset;
      private final int count;

      @JsConstructor
      private EdgeRange(S2Shape shape, int offset, int count) {
        super(shape);
        this.offset = offset;
        this.count = count;
      }

      @Override
      public final int numEdges() {
        return count;
      }

      @Override
      public final int edge(int i) {
        return offset + i;
      }
    }
  }

  /**
   * FaceEdge stores temporary edge data while the index is being updated. FaceEdge represents an
   * edge in the UV coordinates of a specific face, without any clipping. ClippedEdge, by
   * comparison, contains the UV bound of the portion of the edge that intersects each cell.
   *
   * <p>While it would be possible to combine all the edge information into ClippedEdge, there will
   * be many clipped edges for each original face edge, and only the UV bound is different. Keeping
   * the shared fields on this separate class provides two advantages:
   *
   * <ul>
   *   <li>Memory usage. Separating the two classes means that we only need to store one copy of the
   *       per-face data no matter how many times an edge is subdivided, and it also lets us delay
   *       computing bounding boxes until they are needed for processing each face (when the dataset
   *       spans multiple faces).
   *   <li>Performance. UpdateEdges is significantly faster on large polygons when the data is
   *       separated, because it often only needs to access the data in ClippedEdge and this data is
   *       cached more successfully.
   * </ul>
   */
  private static final class FaceEdge {
    /** The shape that this edge belongs to. */
    private final int shapeId;

    /** Edge id within that shape. */
    private final int edgeId;

    /** Not desirable to subdivide this edge beyond this level. */
    private final int maxLevel;

    /** The edge endpoints, clipped to a given face. */
    private final double ax;

    private final double ay;
    private final double bx;
    private final double by;

    /**
     * The corresponding S2Points, cached here to avoid repeated calls to {@link
     * S2Shape#getEdge(int, com.google.common.geometry.S2Shape.MutableEdge)}.
     */
    private final S2Point va;

    private final S2Point vb;

    private FaceEdge(
        int shapeId,
        int edgeId,
        S2Point va,
        S2Point vb,
        R2Vector a,
        R2Vector b,
        double cellSizeToLongEdgeRatio) {
      this.shapeId = shapeId;
      this.edgeId = edgeId;
      this.ax = a.x;
      this.ay = a.y;
      this.bx = b.x;
      this.by = b.y;
      this.va = va;
      this.vb = vb;
      this.maxLevel = getEdgeMaxLevel(va, vb, cellSizeToLongEdgeRatio);
    }

    @Override
    public String toString() {
      return "shape " + shapeId + " edge " + edgeId;
    }
  }

  /**
   * Returns the first level for which the given edge will be considered "long", i.e. it will not
   * count towards the {@link Options#maxEdgesPerCell} limit.
   */
  @VisibleForTesting
  static final int getEdgeMaxLevel(S2Point va, S2Point vb, double cellSizeToLongEdgeRatio) {
    // Compute the maximum cell edge length for which this edge is considered "long". The
    // calculation does not need to be perfectly accurate, so avoid angle() since it's slower.
    double maxCellEdge = va.getDistance(vb) * cellSizeToLongEdgeRatio;
    // Now set the first level encountered during subdivision where the average cell edge length at
    // that level is at most "cellSize".
    return S2Projections.PROJ.avgEdge.getMinLevel(maxCellEdge);
  }

  /** ClippedEdge represents the portion of a FaceEdge that has been clipped to an S2Cell. */
  @JsType
  private static class ClippedEdge {
    /**
     * The original unclipped edge. This field is not final so we can reuse ClippedEdge instances in
     * an object pool, {@link EdgeAllocator}.
     */
    private FaceEdge orig;

    /** Bounding box for the clipped portion. */
    private final R2Rect bound = new R2Rect();
  }

  /**
   * This class provides temporary storage for new ClippedEdges that are created during indexing. It
   * is essentially a stack-based object pool, where edges are allocated as the recursion goes down
   * the first time, put back in the pool as recursion come back up, and reused when recursion goes
   * back down another branch of the S2Cell tree.
   */
  private static final class EdgeAllocator {
    private int size;
    private final List<ClippedEdge> edges;

    public EdgeAllocator(int maxEdges) {
      edges = createList(maxEdges);
    }

    /** Returns an edge. */
    public ClippedEdge create() {
      if (size == edges.size()) {
        edges.add(new ClippedEdge());
      }
      return edges.get(size++);
    }

    /**
     * Returns the number of allocated edges. Before a thread calls {@link #create()}, this method
     * should be called to assess the size of the stack, and after all created edges are no longer
     * needed, call {@link #reset(int)} with the previous size.
     */
    public int size() {
      return size;
    }

    /** Returns all edges after 'size' to the object pool to be reused by another thread. */
    public void reset(int size) {
      this.size = size;
    }
  }

  /**
   * Given a set of shapes, InteriorTracker keeps track of which shapes contain a particular point
   * (the "focus".) It provides an efficient way to move the focus from one point to another and
   * incrementally update the set of shapes which contain it. We use this to compute which shapes
   * contain the center of every S2CellId in the index, by advancing the focus from one cell center
   * to the next.
   *
   * <p>Initially the focus is S2.origin(), and therefore we can initialize the state of every shape
   * to its containsOrigin() value. Next we advance the focus to the start of the S2CellId space-
   * filling curve, by drawing a line segment between this point and S2.origin() and testing whether
   * every edge of every shape intersects it. Then we visit all the cells that are being added to
   * the S2ShapeIndex in increasing order of S2CellId.
   *
   * <p>For each cell, we draw two edges: one from the entry vertex to the center, and another from
   * the center to the exit vertex (where "entry" and "exit" refer to the points where the space-
   * filling curve enters and exits the cell). By counting edge crossings we can incrementally
   * compute which shapes contain the cell center.
   *
   * <p>Note that the same set of shapes will always contain the exit point of one cell and the
   * entry point of the next cell in the index, because either (a) these two points are actually the
   * same, or (b) the intervening cells in S2CellId order are all empty, and therefore there are no
   * edge crossings if we follow this path from one cell to the other.
   */
  static final class InteriorTracker {
    /** Whether any shapes have an interior. */
    private boolean isActive;

    /** The prior focus point. */
    private S2Point oldFocus;

    /** The new focus point. */
    private S2Point newFocus;

    /** The last exit vertex. */
    private S2Point lastEnd;

    /**
     * The ideal next cell ID such that the entry vertex of that cell would match the exit vertex of
     * this cell.
     */
    private S2CellId nextCellId;

    /** A temporary crosser from the old focus to the new focus. */
    private S2EdgeUtil.EdgeCrosser crosser;

    /**
     * The set of shape ids (the indices of each shape in the S2ShapeIndex.shapes field) that
     * contain the current focus. Sorted by ID in ascending order.
     */
    private final int[] focusedShapes;

    /** The number of elements in 'focusedShapes' that are valid and in use. */
    private int focusCount;

    /** A temporary array in which to accumulate the clipped shapes for each cell. */
    private final S2ClippedShape[] tempClippedShapes;

    /**
     * Initializes the InteriorTracker. You must call {@link #addShape(int, S2Shape)} for each shape
     * that will be tracked before calling {@link #moveTo(S2Point)} or {@link #drawTo(S2Point)}.
     */
    public InteriorTracker(int numShapes) {
      this.tempClippedShapes = new S2ClippedShape[numShapes];
      this.focusedShapes = new int[numShapes];
      this.isActive = false;
      this.nextCellId = S2CellId.begin(S2CellId.MAX_LEVEL);
      // Draw from S2.origin() to the entry vertex of the very first cell. When shapes are added,
      // they will be initialized based on whether they contain the origin, and the call to testEdge
      // will then move the focus to the entry vertex of the first cell.
      this.newFocus = S2.origin();
      drawTo(S2Point.normalize(S2Projections.faceUvToXyz(0, -1, -1))); // S2CellId curve start
    }

    /** Returns true if {@link #addShape(int, S2Shape)} has been called at least once. */
    public boolean isActive() {
      return isActive;
    }

    /**
     * Adds a shape whose interior should be tracked. This should be followed by calling {@link
     * #testEdge(int, S2Point, S2Point)} with every edge of the given shape.
     */
    public void addShape(int shapeId, S2Shape shape) {
      // assert shapeId == shapes.indexOf(shape);
      isActive = true;
      if (shape.containsOrigin()) {
        toggleShape(shapeId);
      }
    }

    /**
     * Moves the focus to the given point. This method should only be used when it is known that
     * there are no edge crossings between the old and new focus locations; otherwise use {@link
     * #drawTo(S2Point)}.
     */
    public void moveTo(S2Point b) {
      newFocus = b;
    }

    /**
     * Moves the focus to the given point. After this method is called, {@link #testEdge(int,
     * S2Point, S2Point)} should be called with all edges that may cross the line segment between
     * the old and new focus locations.
     */
    public void drawTo(S2Point focus) {
      oldFocus = newFocus;
      newFocus = focus;
      crosser = new S2EdgeUtil.EdgeCrosser(oldFocus, newFocus);
      lastEnd = null;
    }

    /**
     * Tests whether the given edge of the given shape may cross the line segment between the old
     * and new focus locations (see {@link #drawTo(S2Point)}), and if there is a crossing the
     * shape's containment of the focus is toggled.
     */
    public void testEdge(int shapeId, S2Point start, S2Point end) {
      // Just check that 'lastEnd' set by the last call is the new 'start', so comparison by
      // reference is okay.
      if (start != lastEnd) {
        crosser.restartAt(start);
      }
      if (crosser.edgeOrVertexCrossing(end)) {
        toggleShape(shapeId);
      }
      lastEnd = end;
    }

    /**
     * Indicates that the caller has finished processing the given S2CellId. By using this method
     * together with {@link #atCellId(S2CellId)}, the caller can avoid calling {@link
     * #moveTo(S2Point)} in cases where the exit vertex of the previous cell is the same as the
     * entry vertex of the current cell.
     */
    public void doneCellId(S2CellId cellid) {
      nextCellId = cellid.rangeMax().next();
    }

    /**
     * Returns true if the focus is already at the entry vertex of the given S2CellId (provided that
     * the caller calls {@link #doneCellId(S2CellId)} as each cell is processed).
     */
    public boolean atCellId(S2CellId cellid) {
      return cellid.rangeMin().id() == nextCellId.id();
    }

    /**
     * Toggles the given shape ID from the list of shapes that contain the current focus; if the
     * shape was not in the set, add it; if it was in the set, remove it.
     */
    private void toggleShape(int shapeId) {
      // Since focusCount is typically *very* small (0, 1, or 2), it turns out to be significantly
      // faster to maintain a sorted array rather than using a Set<Integer>.
      if (focusCount == 0) {
        focusedShapes[0] = shapeId;
        focusCount++;
      } else if (focusedShapes[0] == shapeId) {
        if (focusCount-- > 1) {
          System.arraycopy(focusedShapes, 1, focusedShapes, 0, focusCount);
        }
      } else {
        int pos = 0;
        while (focusedShapes[pos] < shapeId) {
          if (++pos == focusCount) {
            focusedShapes[focusCount++] = shapeId;
            return;
          }
        }
        if (focusedShapes[pos] == shapeId) {
          focusCount--;
          System.arraycopy(focusedShapes, pos + 1, focusedShapes, pos, focusCount - pos);
        } else {
          System.arraycopy(focusedShapes, pos, focusedShapes, pos + 1, focusCount - pos);
          focusedShapes[pos] = shapeId;
          focusCount++;
        }
      }
    }
  }

  /**
   * Creates a new list, using a SimpleList when the predicted maximum size is small, and a sharded
   * list when the predicted size is large enough to be worth it.
   */
  static final <T> List<T> createList(int maxSize) {
    if (maxSize < 256) {
      return new SimpleList<T>(maxSize);
    } else {
      return new ShardedList<T>(maxSize);
    }
  }

  /** A simple append-only RandomAccess List similar to (but about 10% faster than) ArrayList. */
  private static final class SimpleList<T> extends AbstractList<T>
      implements RandomAccess, Serializable {
    private static final long serialVersionUID = 1L;

    private Object[] elements;
    private int size;

    public SimpleList(int maxSize) {
      elements = new Object[max(1, maxSize)];
    }

    @Override
    public T get(int index) {
      @SuppressWarnings("unchecked")
      T result = (T) elements[index];
      return result;
    }

    @Override
    public T set(int index, T value) {
      @SuppressWarnings("unchecked")
      T old = (T) elements[index];
      elements[index] = value;
      return old;
    }

    @Override
    public int size() {
      return size;
    }

    @Override
    public boolean add(T item) {
      if (size == elements.length) {
        elements = Arrays.copyOf(elements, size * 2);
      }
      elements[size++] = item;
      return true;
    }
  }

  /**
   * A more complex append-only RandomAccess List that allocates space in shards of 256 elements
   * each, avoiding reallocation as the list grows, and avoiding single allocations larger than 2KB.
   * This is about 1% faster than SimpleList with reserveEdges in use, but when the heap is anywhere
   * near full, this approach is dramatically faster than any technique requiring a large contiguous
   * allocation, since it will probably require a GC to supply one.
   */
  private static final class ShardedList<T> extends AbstractList<T>
      implements RandomAccess, Serializable {
    private static final long serialVersionUID = 1L;

    private Object[][] elements;
    private int size;

    public ShardedList(int maxItems) {
      // Don't allocate the inner list until items are actually added.
      elements = new Object[1 + (maxItems >> 8)][];
    }

    @Override
    public int size() {
      return size;
    }

    @Override
    public boolean add(T item) {
      int shard = size >> 8;
      if (shard == elements.length) {
        elements = Arrays.copyOf(elements, shard * 2);
        elements[shard] = new Object[256];
      } else if (elements[shard] == null) {
        elements[shard] = new Object[256];
      }
      elements[shard][size & 0xFF] = item;
      size++;
      return true;
    }

    @Override
    public T get(int index) {
      @SuppressWarnings("unchecked")
      T result = (T) elements[index >> 8][index & 0xFF];
      return result;
    }

    @Override
    public T set(int index, T value) {
      @SuppressWarnings("unchecked")
      T result = (T) elements[index];
      elements[index >> 8][index & 0xFF] = value;
      return result;
    }
  }
}
