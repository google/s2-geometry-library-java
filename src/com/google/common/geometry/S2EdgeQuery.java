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

import com.google.common.annotations.GwtCompatible;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;

/**
 * S2EdgeQuery is used to find edges or shapes that are crossed by an edge. If you need to query
 * many edges, it is more efficient to declare a single S2EdgeQuery object and reuse it so that
 * temporary storage does not need to be reallocated each time.
 * 
 * <p>This class is not thread-safe.
 */
@GwtCompatible
public class S2EdgeQuery {
  // TODO(eengle): Make this class public once getCandidates has been optimized to avoid boxing edge
  // IDs.
  private final S2ShapeIndex index;
  /** Temporary list of cells that intersect the query edge AB. Used while processing a query. */
  private final List<S2ShapeIndex.Cell> cells;
  /** The following vectors are temporary storage used while processing a query. */
  private final S2ShapeIndex.CellIterator iter;
  /** An {@code Edges} implementation that contains no edges. */
  private static final Edges EMPTY_EDGE_LIST = new Edges() {
    @Override
    public int nextEdge() {
      return -1;
    }
    @Override
    public boolean isEmpty() {
      return true;
    }
  };
  
  /** Constructor from an {@link S2ShapeIndex}. */
  public S2EdgeQuery(S2ShapeIndex index) {
    this.index = index;
    this.iter = index.iterator();
    cells = Lists.newArrayList();
  }
  
  /**
   * Given a query edge AB and a shape {@code shape}, returns a superset of the edges of {@code
   * shape} that intersect AB. Consider using {@link ShapeEdges} instead, if the shape has few
   * enough edges.
   */
  public Edges getCandidates(S2Point a, S2Point b, S2Shape shape) {
    // For small loops it is faster to use brute force. The threshold below was determined using the
    // benchmarks in the C++ S2Loop unit test.
    // TODO(eengle) Update this value based on benchmarking in Java.
    int maxBruteForceEdges = 27;
    int maxEdges = shape.numEdges();
    if (maxEdges <= maxBruteForceEdges) {
      return new ShapeEdges(shape.numEdges());
    }

    getCells(a, b);

    // Compute and return the 'Edges', using different 'Edges' implementations based on how many
    // cells the query edge covers.
    if (cells.isEmpty()) {
      return EMPTY_EDGE_LIST;
    } else if (cells.size() == 1) {
      S2ClippedShape clippedShape = cells.get(0).findClipped(shape);
      if (clippedShape == null || clippedShape.numEdges() == 0) {
        return EMPTY_EDGE_LIST;
      } else {
        return new SimpleEdges(clippedShape);
      }
    } else {
      MergedEdges edges = new MergedEdges();
      for (int i = 0; i < cells.size(); ++i) {
        S2ClippedShape clippedShape = cells.get(i).findClipped(shape);
        if (clippedShape != null && clippedShape.numEdges() != 0) {
          edges.add(clippedShape);
        }
      }
      return edges;
    }
  }

  /**
   * Given a query edge AB, returns a map from the indexed shapes to a superset of the edges for
   * each shape that intersect AB. Consider using {@link ShapeEdges} instead, if there is just one
   * indexed shape with few enough edges.
   */
  public Map<S2Shape, Edges> getCandidates(S2Point a, S2Point b) {
    // If there are only a few edges then it's faster to use brute force. We only bother with this
    // optimization when there is a single shape, since then we can also use some tricks to avoid
    // reallocating the edge map.
    if (index.shapes.size() == 1) {
      S2Shape shape = index.shapes.get(0);
      Edges edges = getCandidates(a, b, shape);
      if (edges.isEmpty()) {
        return Collections.<S2Shape, Edges>emptyMap();
      } else {
        return Collections.<S2Shape, Edges>singletonMap(shape, edges);
      }
    }

    getCells(a, b);

    // Compute and return the map.  If the map is nonempty, use different 'Edges' implementations
    // based on how many cells the query edge covers.
    if (cells.isEmpty()) {
      return Collections.emptyMap();
    } else if (cells.size() == 1) {
      S2ShapeIndex.Cell cell = cells.get(0);
      if (cell.numShapes() == 1) {
        S2ClippedShape clippedShape = cell.clipped(0);
        if (clippedShape.numEdges() == 0) {
          return Collections.<S2Shape, Edges>emptyMap();
        } else {
          S2Shape shape = cell.clipped(0).shape();
          return Collections.<S2Shape, Edges>singletonMap(shape, new SimpleEdges(clippedShape));
        }
      }
      Map<S2Shape, Edges> edgeMap = Maps.newIdentityHashMap();
      for (int j = 0; j < cell.numShapes(); ++j) {
        S2ClippedShape clippedShape = cell.clipped(j);
        if (clippedShape.numEdges() > 0) {
          S2Shape shape = clippedShape.shape();
          edgeMap.put(shape, new SimpleEdges(clippedShape));
        }
      }
      return edgeMap;
    } else {
      Map<S2Shape, Edges> edgeMap = Maps.newIdentityHashMap();
      for (int i = 0; i < cells.size(); ++i) {
        S2ShapeIndex.Cell cell = cells.get(i);
        for (int j = 0; j < cell.numShapes(); ++j) {
          S2ClippedShape clippedShape = cell.clipped(j);
          if (clippedShape.numEdges() == 0) {
            continue;
          }
          S2Shape shape = clippedShape.shape();
          MergedEdges edges = (MergedEdges) edgeMap.get(shape);
          if (edges == null) {
            edges = new MergedEdges();
            edgeMap.put(shape, edges);
          }
          edges.add(clippedShape);
        }
      }
      return edgeMap;
    }
  }

  /** Sets cells to the set of index cells intersected by an edge AB. */
  private void getCells(S2Point a, S2Point b) {
    cells.clear();
    S2EdgeUtil.FaceSegment[] segments = new S2EdgeUtil.FaceSegment[6];
    int numSegments = S2EdgeUtil.getFaceSegments(a, b, segments);
    for (int i = 0; i < numSegments; ++i) {
      // Optimization: rather than always starting the recursive subdivision at the top level face
      // cell, we start at the smallest S2CellId that contains the edge (the "edge root cell"). This
      // typically lets us skip quite a few levels of recursion, since most edges are short.
      R2Rect edgeBound = R2Rect.fromPointPair(segments[i].a, segments[i].b);
      S2PaddedCell pCell = new S2PaddedCell(S2CellId.fromFace(segments[i].face), 0);
      S2CellId edgeRoot = pCell.shrinkToFit(edgeBound);

      // Now we need to determine how the edge root cell is related to the cells in the spatial
      // index ('cells' in S2ShapeIndex.java). There are three cases:
      //
      // 1. edgeRoot is an index cell or is contained within an index cell. In this case, we only
      //    need to look at the contents of that cell.
      // 2. edgeRoot is subdivided into one or more index cells. In this case we recursively
      //    subdivide to find the cells intersected by AB.
      // 3. edgeRoot does not intersect any index cells. In this case there is nothing to do.
      S2ShapeIndex.CellRelation relation = iter.locate(edgeRoot);
      if (relation == S2ShapeIndex.CellRelation.INDEXED) {
        // edgeRoot is an index cell or is contained by an index cell (case 1).
        // assert (iter.id().contains(edgeRoot));
        cells.add(iter.cell());
      } else if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) {
        // edgeRoot is subdivided into one or more index cells (case 2). We find the cells
        // intersected by AB using recursive subdivision.
        if (!edgeRoot.isFace()) {
          pCell = new S2PaddedCell(edgeRoot, 0);
        }
        getCells(pCell, edgeBound, segments[i].a, segments[i].b);
      }
    }
  }
  
  /**
   * Convenience method for calling {@link #getCells(S2Point, R2Vector, S2Point, R2Vector,
   * S2PaddedCell, List)}.
   */
  public boolean getCells(S2Point a, S2Point b, S2PaddedCell root, List<S2ShapeIndex.Cell> cells) {
    R2Vector aVector = new R2Vector();
    R2Vector bVector = new R2Vector();
    return getCells(a, aVector, b, bVector, root, cells);
  }

  /**
   * Adds all cells to {@code cells} that might intersect the query edge from {@code a} to {@code b}
   * and the cell {@code root}. The {@code aVector} and {@code bVector} parameters are cached R2
   * versions of the [A, B] edge projected onto the same cube face as {@code root}.
   */
  @VisibleForTesting
  boolean getCells(S2Point a, R2Vector aVector, S2Point b, R2Vector bVector, S2PaddedCell root, 
      List<S2ShapeIndex.Cell> cells) {
    this.cells.clear();
    if (S2EdgeUtil.clipToFace(a, b, root.id().face(), aVector, bVector)) {
      R2Rect edgeBound = R2Rect.fromPointPair(aVector, bVector);
      if (root.bound().intersects(edgeBound)) {
        getCells(root, edgeBound, aVector, bVector);
      }
    }
    if (this.cells.isEmpty()) {
      return false;
    }
    cells.addAll(this.cells);
    return true;
  }
  
  /**
   * Computes the index cells intersected by the current edge that are descendants of {@code pCell},
   * and adds them to {@code cells}.
   * 
   * <p>WARNING: This function is recursive with a maximum depth of 30.
   */
  private void getCells(S2PaddedCell pCell, R2Rect edgeBound, R2Vector aVector, R2Vector bVector) {
    iter.seek(pCell.id().rangeMin());
    if (iter.done() || iter.id().greaterThan(pCell.id().rangeMax())) {
      // The index does not contain 'pCell' or any of its descendants.
      return;
    }
    if (iter.id().equals(pCell.id())) {
      // The index contains this cell exactly.
      cells.add(iter.cell());
      return;
    }
    
    // Otherwise, split the edge among the four children of 'pCell'.
    R2Vector center = pCell.middle().lo();
    if (edgeBound.x().hi() < center.x()) {
      // Edge is entirely contained in the two left children.
      clipVAxis(edgeBound, center.y(), 0, pCell, aVector, bVector);
    } else if (edgeBound.x().lo() >= center.x()) {
      // Edge is entirely contained in the two right children.
      clipVAxis(edgeBound, center.y(), 1, pCell, aVector, bVector);
    } else {
      R2Rect[] childBounds = new R2Rect[2];
      splitUBound(edgeBound, center.x(), childBounds, aVector, bVector);
      if (edgeBound.y().hi() < center.y()) {
        // Edge is entirely contained in the two lower children.
        getCells(pCell.childAtIJ(0, 0), childBounds[0], aVector, bVector);
        getCells(pCell.childAtIJ(1, 0), childBounds[1], aVector, bVector);
      } else if (edgeBound.y().lo() >= center.y()) {
        // Edge is entirely contained in the two upper children.
        getCells(pCell.childAtIJ(0, 1), childBounds[0], aVector, bVector);
        getCells(pCell.childAtIJ(1, 1), childBounds[1], aVector, bVector);
      } else {
        // The edge bound spans all four children. The edge itself intersects at most three children
        // (since no padding is being used).
        clipVAxis(childBounds[0], center.y(), 0, pCell, aVector, bVector);
        clipVAxis(childBounds[1], center.y(), 1, pCell, aVector, bVector);
      }
    }
  }
  
  /**
   * Given either the left (i = 0) or right (i = 1) side of a padded cell {@code pCell}, determines
   * whether the current edge intersects the lower child, upper child, or both children, and calls
   * getCells() recursively on those children. {@code center} is the v-coordinate at the center of
   * {@code pCell}.
   */
  private void clipVAxis(R2Rect edgeBound, double center, int i, S2PaddedCell pCell,
      R2Vector aVector, R2Vector bVector) {
    if (edgeBound.y().hi() < center) {
      // Edge is entirely contained in the lower child.
      getCells(pCell.childAtIJ(i, 0), edgeBound, aVector, bVector);
    } else if (edgeBound.y().lo() >= center) {
      // Edge is entirely contained in the upper child.
      getCells(pCell.childAtIJ(i, 1), edgeBound, aVector, bVector);
    } else {
      // The edge intersects both children.
      R2Rect[] childBounds = new R2Rect[2];
      splitVBound(edgeBound, center, childBounds, aVector, bVector);
      getCells(pCell.childAtIJ(i, 0), childBounds[0], aVector, bVector);
      getCells(pCell.childAtIJ(i, 1), childBounds[1], aVector, bVector);
    }
  }

  /**
   * Splits the current edge into two child edges at {@code u} and returns the bound for each child.
   */  
  private void splitUBound(R2Rect edgeBound, double u, R2Rect[] childBounds, R2Vector aVector,
      R2Vector bVector) {
    // See comments in S2ShapeIndex.clipUBound.
    double v = edgeBound.y().
        clampPoint(S2EdgeUtil.interpolateDouble(u, aVector.x, bVector.x, aVector.y, bVector.y));
    
    // 'diag' indicates which diagonal of the bounding box is spanned by AB: It is 0 if AB has
    // positive slope, and 1 if AB has negative slope.
    int diag = ((aVector.x > bVector.x) != (aVector.y > bVector.y)) ? 1 : 0;
    splitBound(edgeBound, 0, u, diag, v, childBounds);
  }
  
  /**
   * Splits the current edge into two child edges at {@code v} and returns the bound for each child.
   */
  private void splitVBound(R2Rect edgeBound, double v, R2Rect[] childBounds, R2Vector aVector,
      R2Vector bVector) {
    double u = edgeBound.x().
        clampPoint(S2EdgeUtil.interpolateDouble(v, aVector.y, bVector.y, aVector.x, bVector.x));
    int diag = ((aVector.x > bVector.x) != (aVector.y > bVector.y)) ? 1 : 0;
    splitBound(edgeBound, diag, u, 0, v, childBounds);
  }
  
  /**
   * Splits the current edge into two child edges at the given point (u, v) and returns the bound
   * for each child. {@code uEnd} and {@code vEnd} indicate which bound endpoints of child 1 will
   * be updated.
   */
  private void splitBound(R2Rect edgeBound, int uEnd, double u, int vEnd, double v,
      R2Rect[] childBounds) {
    childBounds[0] = new R2Rect(edgeBound);
    childBounds[1] = new R2Rect(edgeBound);
    if (uEnd == 0) {
      childBounds[0].x().setHi(u);
      childBounds[1].x().setLo(u);
    } else {
      childBounds[0].x().setLo(u);
      childBounds[1].x().setHi(u);
    }
    if (vEnd == 0) {
      childBounds[0].y().setHi(v);
      childBounds[1].y().setLo(v);
    } else {
      childBounds[0].y().setLo(v);
      childBounds[1].y().setHi(v);
    }
    // assert (!childBounds[0].isEmpty());
    // assert (edgeBound.contains(childBounds[0]));
    // assert (!childBounds[1].isEmpty());
    // assert (edgeBound.contains(childBounds[1]));
  }
  
  /** An iterator over the sorted unique edge IDs of a shape that may intersect some query edge. */
  public interface Edges {
    /** Returns the next edge ID, or throws an exception if empty. */
    int nextEdge();
    
    /** Returns true if there are no more edges. */
    boolean isEmpty();
  }
  
  /** An {@code Edges} implementation that includes all the edges of a clipped shape. */
  private static final class SimpleEdges implements Edges {
    int index;
    final S2ClippedShape shape;
    
    SimpleEdges(S2ClippedShape shape) {
      index = 0;
      this.shape = shape;
    }
    
    @Override
    public int nextEdge() {
      Preconditions.checkState(!isEmpty(), "Cannot call nextEdge() on empty Edges.");
      if (index == shape.numEdges()) {
        return -1;
      } else {
        return shape.edge(index++);
      }
    }

    @Override
    public boolean isEmpty() {
      return index == shape.numEdges();
    }
  }
  
  /**
   * An {@code Edges} implementation optimized for merging edges from multiple S2ClippedShapes
   * already in sorted order.
   */
  private static final class MergedEdges implements Edges {
    final PriorityQueue<Stepper> steppers = new PriorityQueue<Stepper>();
    /**
     * The top of the priority queue (the stepper which currently has the least value for {@code
     * currentEdge}). It is stored separately as an optimization, to avoid repeatedly adding and
     * polling it from the top of the queue.
     */
    Stepper top;
    
    /** Note: {@code shape} should have at least one edge. */
    public void add(S2ClippedShape shape) {
      Stepper stepper = new Stepper(shape);

      if (top == null) {
        top = stepper;
      } else if (top.currentEdge() <= stepper.currentEdge()) {
        steppers.add(stepper);
      } else {
        steppers.add(top);
        top = stepper;
      }
    }
    
    @Override
    public int nextEdge() {
      Preconditions.checkState(!isEmpty(), "Cannot call nextEdge() on empty Edges.");

      // Store the value to be returned.
      int nextEdge = top.currentEdge();
      removeFromPriorityQueue(nextEdge);
      
      top.index++;
      if (top.index == top.clipped.numEdges()) {
        // Exhausted current stepper. Get the next one.
        top = steppers.isEmpty() ? null : steppers.poll();
      } else if (!steppers.isEmpty() && top.currentEdge() > steppers.peek().currentEdge()) {
        // New top.currentEdge() is no longer the next sorted edge, so swap top with the head of the
        // queue.
        steppers.add(top);
        top = steppers.poll();
      }
      
      return nextEdge;
    }
    
    /**
     * Updates the priority queue {@code steppers} so that no stepper in the queue will return
     * {@code n} if {@code currentEdge()} is called on it.
     */
    private void removeFromPriorityQueue(int n) {
      while (!steppers.isEmpty() && steppers.peek().currentEdge() == n) {
        Stepper stepper = steppers.poll();
        stepper.index++;
        if (stepper.index != stepper.clipped.numEdges()) {
          steppers.add(stepper);
        }
      }
    }

    @Override
    public boolean isEmpty() {
      return top == null;
    }
  }
  
  /** An {@code Edges} that contains all the edges of a shape with the given number of edges. */
  public static final class ShapeEdges implements Edges {
    private int edgeIndex = 0;
    private final int numEdges;
    
    ShapeEdges(int numEdges) {
      this.numEdges = numEdges;
    }
    
    @Override
    public int nextEdge() {
      Preconditions.checkState(!isEmpty(), "Cannot call nextEdge() on empty Edges.");
      return edgeIndex < numEdges ? edgeIndex++ : -1;
    }
    
    @Override
    public boolean isEmpty() {
      return edgeIndex == numEdges;
    }
  }
  
  /** Tracks the current edge index within a clipped shape. */
  private static final class Stepper implements Comparable<Stepper> {
    int index;
    final S2ClippedShape clipped;
   
    Stepper(S2ClippedShape shape) {
      this.index = 0;
      this.clipped = shape;
    }
    
    int currentEdge() {
      return clipped.edge(index);
    }
    
    @Override
    public int compareTo(Stepper that) {
      return Integer.compare(this.currentEdge(), that.currentEdge());
    }
  }
}