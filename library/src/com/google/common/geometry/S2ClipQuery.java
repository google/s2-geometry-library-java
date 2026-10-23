/*
 * Copyright 2026 Google LLC.
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
import static com.google.common.geometry.S2ShapeUtil.lowerBound;

import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Shape.ChainPosition;
import com.google.common.geometry.S2Shape.ReferencePoint;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.PrimitiveIterator;
import java.util.PriorityQueue;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.IntPredicate;

/**
 * S2ClipQuery efficiently selects geometry from an S2ShapeIndex that covers a given set of regions,
 * preserving the exact edges that intersect the regions. This stability is essential when clipping
 * abutting regions such as tiles from the same index. It is also considerably more efficient than
 * exact boolean operations.
 *
 * <p>Example of clipping an index by an {@link S2LatLngRect}:
 *
 * {@snippet :
 * S2ClipQuery query = new S2ClipQuery(index);
 * query.clip(
 *     S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(-10, -10), S2LatLng.fromDegrees(10, 10)));
 * for (int shapeId : query.shapes()) {
 *   for (S2ClipQuery.Chain chain : query.chains(shapeId)) {
 *     for (S2Point vertex : chain.vertices()) {
 *       // Process vertex...
 *     }
 *   }
 * }
 * }
 *
 * <p>It is possible to filter the shapes that are returned using {@link #setShapeFilter}. This can
 * be used to avoid processing shapes that are not relevant to the query:
 *
 * {@snippet :
 * S2ClipQuery query = new S2ClipQuery(index);
 * query.setShapeFilter(shapeId -> isShapeRelevant(shapeId));
 * query.clip(region);
 * }
 *
 * <p>It is also possible to accumulate clipped shapes from multiple sources using {@link #addClip}.
 * Where the various clip methods replace existing results, addClip methods add to existing results.
 * For example, to find shapes that intersect either a {@link S2CellUnion} or a complex polygon
 * indexed by an {@link S2ShapeIndex}:
 *
 * {@snippet :
 * S2ClipQuery query = new S2ClipQuery(index);
 * query.addClip(cellUnion);
 * query.addClip(polygonIndex);
 * for (int shapeId : query.shapes()) {
 *   // Process distinct shapes / edges found in either region...
 * }
 * }
 *
 * <p>For point shapes, it is more efficient to access edges directly:
 *
 * {@snippet :
 * for (int shapeId : query.shapes()) {
 *   S2Shape shape = index.getShapes().get(shapeId);
 *   if (shape.dimension() == 0) {
 *     PrimitiveIterator.OfInt it = query.edgeIds(shapeId);
 *     while (it.hasNext()) {
 *       int edgeId = it.nextInt();
 *       S2Point point = shape.getChainVertex(edgeId, 0);
 *       // Process point...
 *     }
 *   }
 * }
 * }
 *
 * <p>For polygon shapes, there may still be intersection even if the polygon has no edges or chains
 * when the entire clip area is contained by the polygon. In this case, the shape is still included
 * in the {@link #shapes()} iterator, but {@link #edgeIds(int)} will be empty. The {@link
 * #referencePoint(int)} method will return a contained reference point near the clip area:
 *
 * {@snippet :
 * List<S2Shape> containing = new ArrayList<>();
 * for (int shapeId : query.shapes()) {
 *   S2Shape shape = index.getShapes().get(shapeId);
 *   if (shape.dimension() == 2 && !query.edgeIds(shapeId).hasNext()) {
 *     // The polygon contains the clip region.
 *     containing.add(shape);
 *   }
 * }
 * }
 *
 * <p>This class is not thread-safe.
 */
public class S2ClipQuery {
  /** The S2ShapeIndex to query. */
  private S2ShapeIndex index;

  /** An iterator over the S2ShapeIndex. */
  private S2CellRangeIterator<Cell> iter;

  /** The list of clipped shapes found so far. */
  private final List<S2ClippedShape> clippedShapes = new ArrayList<>();

  /** A map from shape ID to reference point for polygon shapes. */
  private final Int2ObjectOpenHashMap<ReferencePoint> refPoints = new Int2ObjectOpenHashMap<>();

  /** A priority queue of EdgeIterators, used to merge edges from multiple S2ClippedShapes. */
  private PriorityQueue<EdgeIterator> pq;

  /** The index cells used by each {@code addClip} call. */
  private final SortedSet<S2CellId> cells = new TreeSet<>();

  private IntPredicate shapeFilter = id -> true;

  /** Creates a new S2ClipQuery for the given S2ShapeIndex. */
  public S2ClipQuery(S2ShapeIndex index) {
    init(index);
  }

  /** Creates a new S2ClipQuery with no index. Must call {@link #init} before clip methods. */
  public S2ClipQuery() {}

  /** Initializes this query with the given S2ShapeIndex. */
  public void init(S2ShapeIndex index) {
    clear();
    this.index = index;
    this.iter = S2CellRangeIterator.makeS2CellRangeIterator(index);
  }

  /** Returns the S2ShapeIndex currently set in this query or null if not yet initialized. */
  public S2ShapeIndex index() {
    return index;
  }

  /** Sets a filter to avoid processing shapes. */
  public void setShapeFilter(IntPredicate filter) {
    this.shapeFilter = filter;
  }

  //
  // These methods set the clip to the given argument. Each will replaces any previous results.
  //

  /** Sets this query to the portion of 'index' that intersects 'cellUnion'. */
  public void clip(S2CellUnion cellUnion) {
    clear();
    addClip(cellUnion);
  }

  /** Sets this query to the portion of 'index' that intersects 'cellId'. */
  public void clip(S2CellId cellId) {
    clear();
    addClip(cellId);
  }

  /** Sets this query to the portion of 'index' that intersects 'index'. */
  public void clip(S2ShapeIndex index) {
    clear();
    addClip(index);
  }

  /**
   * Sets this query to the portion of 'index' that intersects 'it'.
   *
   * <p>Note: if the cells of the iterator are a superset of the actual clip region, then the clip
   * operation may select far more data than desired. {@link #clip(S2Region)} may be used to refine
   * the selection as far as the cells of the index we are clipping from.
   *
   * @param it the iterator to clip shapes to
   */
  public void clip(S2Iterator<?> it) {
    clear();
    addClip(it);
  }

  /** Sets this query to the portion of 'index' that intersects 'region'. */
  public void clip(S2Region region) {
    clear();
    addClip(region);
  }

  //
  // These methods add the given argument to the clip result.
  //

  /** Clears the results of previous clip / addClip calls. */
  public void clear() {
    clippedShapes.clear();
    refPoints.clear();
    cells.clear();
  }

  /** Adds to this query the portion of 'index' that intersects 'cellUnion'. */
  public void addClip(S2CellUnion cellUnion) {
    addClip(S2Iterator.fromList(cellUnion.cellIds()));
  }

  /** Adds to this query the portion of 'index' that intersects 'cellId'. */
  public void addClip(S2CellId cellId) {
    addClip(S2Iterator.fromList(ImmutableList.of(cellId)));
  }

  /** Adds to this query the portion of 'this.index' that intersects 'index'. */
  public void addClip(S2ShapeIndex index) {
    addClip(new S2ShapeIndexRegion(index));
  }

  /**
   * Adds to this query the portion of 'index' that intersects 'it'.
   *
   * <p>Note: if the cells of the iterator are a superset of the actual clip region, then the clip
   * operation may select far more data than desired. {@link #addClip(S2Region)} may be used to
   * refine the selection as far as the cells of the index we are clipping from.
   */
  public void addClip(S2Iterator<?> it) {
    new S2CellIteratorJoin<>(iter, S2CellRangeIterator.makeS2CellRangeIterator(it))
        .simpleJoin(
            idxIter -> {
              if (cells.add(idxIter.id())) {
                addClippedShapes(idxIter.entry());
              }
              return true;
            });
    clippedShapes.sort(Comparator.comparingInt(S2ClippedShape::shapeId));
  }

  /** Adds to this query the portion of 'index' that intersects 'region'. */
  public void addClip(S2Region region) {
    S2CellUnion covering = new S2CellUnion();
    region.getCellUnionBound(covering.cellIds());
    covering.normalize();
    iter.restart();
    for (S2CellId target : covering.cellIds()) {
      addShapesInRegion(new S2Cell(target), region);
    }
    clippedShapes.sort(Comparator.comparingInt(S2ClippedShape::shapeId));
  }

  /** Adds shapes in the target cell to the query if it intersects the region. */
  private void addShapesInRegion(S2Cell target, S2Region region) {
    S2ShapeIndex.CellRelation relation = iter.locate(target.id());
    if (relation == S2ShapeIndex.CellRelation.DISJOINT) {
      return;
    }

    if (relation == S2ShapeIndex.CellRelation.INDEXED) {
      S2CellId id = iter.id();
      if (cells.add(id)) {
        addClippedShapes(iter.entry());
      }
      return;
    }

    // The target is subdivided into one or more index cells. Recurse on children that intersect
    // the region.
    S2Cell[] children = new S2Cell[4];
    for (int i = 0; i < 4; ++i) {
      children[i] = new S2Cell();
    }
    target.subdivide(children);
    for (int i = 0; i < 4; i++) {
      if (region.mayIntersect(children[i])) {
        addShapesInRegion(children[i], region);
      }
    }
  }

  /** Adds clipped shapes in the target cell to the query. */
  private void addClippedShapes(Cell cell) {
    for (int i = 0; i < cell.numShapes(); i++) {
      S2ClippedShape clipped = cell.clipped(i);
      int shapeId = clipped.shapeId();
      if (!shapeFilter.test(shapeId)) {
        continue;
      }
      if (clipped.numEdges() > 0) {
        clippedShapes.add(clipped);
      }
      if (index.getShapes().get(shapeId).dimension() == 2) {
        refPoints.computeIfAbsent(
            shapeId,
            unused -> {
              if (clipped.numEdges() == 0) {
                clippedShapes.add(clipped);
              }
              return ReferencePoint.create(iter.center(), clipped.containsCenter());
            });
      }
    }
  }

  /** Returns an iterator of the clipped shape IDs. Allowed unboxed iteration. */
  public PrimitiveIterator.OfInt shapesIterator() {
    return new ShapesIterator();
  }

  /** Returns an iterable of the clipped shape IDs. */
  public Iterable<Integer> shapes() {
    return this::shapesIterator;
  }

  /** Returns the reference point of the original shape, or null if the shape isn't a polygon. */
  public ReferencePoint referencePoint(int shapeId) {
    return refPoints.get(shapeId);
  }

  /** Returns an iterator of unique edge IDs clipped from the shape, from {@link #shapes}. */
  public PrimitiveIterator.OfInt edgeIds(int shapeId) {
    int low = lowerBound(0, clippedShapes.size(), i -> clippedShapes.get(i).shapeId() < shapeId);
    checkArgument(low < clippedShapes.size() && clippedShapes.get(low).shapeId() == shapeId);
    PriorityQueue<EdgeIterator> moved;
    if (pq == null) {
      // First time, or a concurrent iteration is not yet exhausted.
      moved = new PriorityQueue<>();
    } else {
      // Subsequent iteration. Reuse the pq and return it after exhaustion to avoid reallocating it.
      moved = pq;
      this.pq = null;
    }
    return new EdgeMergeIterator(moved, low);
  }

  /**
   * Returns an iterable of {@link Chain}s for the given {@code shapeId}, one for each maximal
   * contiguous run of clipped edges belonging to the same shape chain.
   *
   * <p>If a shape chain is clipped into multiple disjoint pieces (or if a closed polygon chain's
   * selected edges wrap around from the last edge {@code getChainLength(chainId) - 1} to edge
   * {@code 0}), each contiguous edge-offset range is returned as a separate {@link Chain} with the
   * same {@link Chain#chainId()}.
   */
  public Iterable<Chain> chains(int shapeId) {
    return () -> new ChainIterator(index.getShapes().get(shapeId), edgeIds(shapeId));
  }

  /**
   * Represents a single nonempty contiguous run of clipped edges within one chain of an {@link
   * S2Shape}.
   *
   * <p>Within chain {@link #chainId()}, the clipped edges form the half-open chain-offset interval
   * {@code [runStart, runLimit)} where {@code 0 <= runStart < runLimit <=
   * shape.getChainLength(chainId)}:
   *
   * <ul>
   *   <li>For point shapes ({@code dimension() == 0}), each point is a single-edge chain with
   *       {@code runStart == 0} and {@code runLimit == 1}. {@link #offsets()} yields {@code 0} and
   *       {@link #vertices()} yields the point.
   *   <li>For polyline and polygon shapes ({@code dimension() > 0}), the {@code m = runLimit -
   *       runStart} contiguous edges at chain offsets {@code [runStart, runLimit)} have {@code m +
   *       1} endpoints at chain vertex offsets {@code [runStart, runLimit]} (inclusive). Offset
   *       {@code k} in {@code [runStart, runLimit)} is the start vertex ({@code v0}) of edge {@code
   *       k}, and offset {@code runLimit} is the end vertex ({@code v1}) of the last edge {@code
   *       runLimit - 1} (which wraps to vertex {@code 0} for a closed polygon loop when {@code
   *       runLimit == shape.getChainLength(chainId)}).
   * </ul>
   */
  public static class Chain {
    private final S2Shape shape;

    /** Index of the parent chain within {@link #shape}, in {@code [0, shape.numChains())}. */
    private final int chainId;

    /**
     * Inclusive starting edge offset of this run within {@link #chainId}, in {@code [0, runLimit)}.
     */
    private final int runStart;

    /**
     * Exclusive ending edge offset of this run within {@link #chainId}, in {@code (runStart,
     * shape.getChainLength(chainId)]}. For non-point shapes ({@code dimension() > 0}), {@code
     * runLimit} is also the chain vertex offset of the final endpoint ({@code v1} of edge {@code
     * runLimit - 1}) passed to {@link S2Shape#getChainVertex(int, int)}.
     */
    private final int runLimit;

    Chain(S2Shape shape, int chainId, int runStart, int runLimit) {
      this.shape = shape;
      this.chainId = chainId;
      this.runStart = runStart;
      this.runLimit = runLimit;
    }

    /** Returns the chain ID within the shape, in {@code [0, shape.numChains())}. */
    public int chainId() {
      return chainId;
    }

    /**
     * Returns the inclusive starting edge offset of this run within {@link #chainId}, in {@code [0,
     * limitOffset())}.
     */
    public int startOffset() {
      return runStart;
    }

    /**
     * Returns the exclusive ending edge offset of this run within {@link #chainId}, in {@code
     * (startOffset(), shape.getChainLength(chainId)]}.
     */
    public int limitOffset() {
      return runLimit;
    }

    /**
     * Returns the number of contiguous edges in this clipped run ({@code limitOffset() -
     * startOffset()}).
     */
    public int numEdges() {
      return runLimit - runStart;
    }

    /**
     * Returns an iterator of chain vertex offsets suitable for {@link S2Shape#getChainVertex(int,
     * int)}.
     *
     * <p>For point shapes ({@code dimension() == 0}), this yields {@code [runStart, runLimit)}
     * (i.e. {@code 0}). For polyline and polygon shapes ({@code dimension() > 0}), this yields
     * {@code [runStart, runLimit]} (inclusive), where {@code runLimit} is the trailing vertex
     * offset for the end vertex ({@code v1}) of edge {@code runLimit - 1}. Note that when {@code
     * runLimit == shape.getChainLength(chainId)}, the trailing offset {@code runLimit} is a valid
     * vertex offset for {@link S2Shape#getChainVertex(int, int)}, not an edge offset for {@link
     * S2Shape#getChainEdge}.
     */
    public PrimitiveIterator.OfInt offsets() {
      return new PrimitiveIterator.OfInt() {
        private final int limit = shape.dimension() > 0 ? runLimit + 1 : runLimit;
        private int nextOffset = runStart;

        @Override
        public boolean hasNext() {
          return nextOffset < limit;
        }

        @Override
        public int nextInt() {
          if (nextOffset < limit) {
            return nextOffset++;
          }
          throw new NoSuchElementException();
        }
      };
    }

    /**
     * Returns an iterable of the vertices of the clipped edge run in chain order: 1 vertex for a
     * point chain ({@code dimension() == 0}), or {@code runLimit - runStart + 1} vertices for a
     * polyline or polygon chain ({@code dimension() > 0}) of {@code runLimit - runStart} contiguous
     * edges.
     */
    public Iterable<S2Point> vertices() {
      return () ->
          new Iterator<S2Point>() {
            private final int limit = shape.dimension() > 0 ? runLimit + 1 : runLimit;
            private int nextOffset = runStart;

            @Override
            public boolean hasNext() {
              return nextOffset < limit;
            }

            @Override
            public S2Point next() {
              if (nextOffset < limit) {
                return shape.getChainVertex(chainId, nextOffset++);
              }
              throw new NoSuchElementException();
            }
          };
    }
  }

  /** An iterator over the clipped chains of a single shape. */
  private static class ChainIterator extends AbstractIterator<Chain> {
    private final S2Shape shape;
    private final PrimitiveIterator.OfInt edgeIter;
    private final ChainPosition chainPos = new ChainPosition();
    private int nextEdge;
    private boolean hasNextEdge;

    ChainIterator(S2Shape shape, PrimitiveIterator.OfInt edgeIter) {
      this.shape = shape;
      this.edgeIter = edgeIter;
      if (edgeIter.hasNext()) {
        nextEdge = edgeIter.nextInt();
        hasNextEdge = true;
      }
    }

    @Override
    protected Chain computeNext() {
      if (!hasNextEdge) {
        return endOfData();
      }

      int edgeStart = nextEdge;
      shape.getChainPosition(edgeStart, chainPos);
      int chainId = chainPos.getChainId();
      int chainStart = shape.getChainStart(chainId);
      int chainLimit = chainStart + shape.getChainLength(chainId);

      // Extend the half-open shape edge-ID range [edgeStart, edgeLimit) as long as edges are
      // contiguous and belong to the same shape chain.
      int edgeLimit = edgeStart;
      do {
        edgeLimit++;
        if (edgeIter.hasNext()) {
          nextEdge = edgeIter.nextInt();
        } else {
          hasNextEdge = false;
          break;
        }
      } while (nextEdge == edgeLimit && nextEdge < chainLimit);

      return new Chain(shape, chainId, edgeStart - chainStart, edgeLimit - chainStart);
    }
  }

  /** An iterator that merges the already-sorted edge IDs from multiple S2ClippedShapes. */
  private class EdgeMergeIterator implements PrimitiveIterator.OfInt {
    private final PriorityQueue<EdgeIterator> pq;
    private EdgeIterator active;

    EdgeMergeIterator(PriorityQueue<EdgeIterator> pq, int startIndex) {
      this.pq = pq;
      if (startIndex < clippedShapes.size()) {
        int shapeId = clippedShapes.get(startIndex).shapeId();
        for (int i = startIndex; i < clippedShapes.size(); i++) {
          S2ClippedShape s = clippedShapes.get(i);
          if (s.shapeId() != shapeId) {
            break;
          }
          if (s.numEdges() > 0) {
            pq.add(new EdgeIterator(s));
          }
        }
      }
      active = pq.poll();
    }

    @Override
    public boolean hasNext() {
      return active != null;
    }

    @Override
    public int nextInt() {
      // Fast path: if the next edge in the active iterator is contiguous with the last edge, we can
      // yield it without consulting the queue.
      int result = active.next();
      if (active.hasNext() && active.currentEdgeId == result + 1) {
        return result;
      }

      // Slow path: put the active iterator back, then poll the best available.
      if (active.hasNext()) {
        pq.add(active);
      }
      while ((active = pq.poll()) != null) {
        if (active.currentEdgeId > result) {
          // The smallest element in the queue is strictly greater than `result`.
          // This is our next valid edge.
          break;
        }
        // The iterator points to a duplicate of `result` (or less). Advance it.
        if (active.advanceTo(result + 1)) {
          // Re-enqueue the advanced iterator so the queue can correctly sort its new edge.
          pq.add(active);
        }
      }

      // Move queue back to S2ClipQuery to reuse it after we're done.
      if (active == null) {
        S2ClipQuery.this.pq = pq;
      }
      return result;
    }
  }

  /** An iterator over the edges of a single S2ClippedShape. */
  private static class EdgeIterator implements Comparable<EdgeIterator> {
    private final S2ClippedShape shape;
    private int index = 0;
    private int currentEdgeId;

    EdgeIterator(S2ClippedShape shape) {
      this.shape = shape;
      this.currentEdgeId = shape.edge(0);
    }

    /** Returns true if there are more edges to iterate. */
    boolean hasNext() {
      return index < shape.numEdges();
    }

    /** Returns the next edge ID. */
    int next() {
      int result = currentEdgeId;
      if (++index < shape.numEdges()) {
        currentEdgeId = shape.edge(index);
      }
      return result;
    }

    /** Advances past the given edge ID and returns {@link #hasNext} as a result. */
    boolean advanceTo(int id) {
      while (++index < shape.numEdges()) {
        currentEdgeId = shape.edge(index);
        if (currentEdgeId >= id) {
          return true;
        }
      }
      return false;
    }

    @Override
    public int compareTo(EdgeIterator other) {
      return Integer.compare(this.currentEdgeId, other.currentEdgeId);
    }
  }

  /** An iterator over the shape IDs found in clip calls. */
  private class ShapesIterator implements PrimitiveIterator.OfInt {
    private int offset = 0;

    @Override
    public boolean hasNext() {
      return offset < clippedShapes.size();
    }

    @Override
    public int nextInt() {
      if (!hasNext()) {
        throw new NoSuchElementException();
      }
      int shapeId = clippedShapes.get(offset).shapeId();

      // Find the end of this shape's run.
      // Use linear scan for small runs, then switch to binary search.
      int limit = offset + 1;
      int size = clippedShapes.size();
      int step = 1;
      while (limit < size && step <= 4) {
        if (clippedShapes.get(limit).shapeId() != shapeId) {
          break;
        }
        limit++;
        step++;
      }

      if (limit < size && clippedShapes.get(limit).shapeId() == shapeId) {
        // Switch to binary search for the upper bound.
        limit = lowerBound(limit + 1, size, i -> clippedShapes.get(i).shapeId() <= shapeId);
      }

      offset = limit;
      return shapeId;
    }
  }
}
