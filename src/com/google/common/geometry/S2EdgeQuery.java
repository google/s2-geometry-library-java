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
import com.google.common.collect.Lists;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * S2EdgeQuery is used to find edges or shapes that are crossed by an edge. If you need to query
 * many edges, it is more efficient to declare a single S2EdgeQuery object and reuse it so that
 * temporary storage does not need to be reallocated each time.
 * 
 * <p>This class is not thread-safe.
 */
@GwtCompatible
class S2EdgeQuery {
  // TODO(eengle): Make this class public once getCandidates has been optimized to avoid boxing edge
  // IDs.
  private final S2ShapeIndex index;
  /** Temporary list of cells that intersect the query edge AB. Used while processing a query. */
  private final List<S2ShapeIndex.Cell> cells;
  /** The following vectors are temporary storage used while processing a query. */
  private final S2ShapeIndex.CellIterator iter;
  /** Constructor from an {@link S2ShapeIndex}. */
  public S2EdgeQuery(S2ShapeIndex index) {
    this.index = index;
    this.iter = index.iterator();
    cells = Lists.newArrayList();
  }
  
  /**
   * Given a query edge AB and a shape {@code shape}, sets {@code edges} equal to a superset of the
   * edges of {@code shape} that intersect AB. Returns false if {@code edges} is empty.
   */
  public boolean getCandidates(S2Point a, S2Point b, S2ShapeIndex.IdShape shape,
      ArrayList<Integer> edges) {
    // For small loops it is faster to use brute force. The threshold below was determined using the
    // benchmarks in the C++ S2Loop unit test.
    // TODO(eengle) Update this value based on benchmarking in Java.
    int maxBruteForceEdges = 27;
    edges.clear();
    int maxEdges = shape.shape.numEdges();
    if (maxEdges <= maxBruteForceEdges) {
      edges.ensureCapacity(maxEdges);
      for (int i = 0; i < maxEdges; ++i) {
        edges.add(i);
      }
      return true;
    }
    
    // Compute the set of index cells intersected by the query edge.
    getCells(a, b);
    if (cells.isEmpty()) {
      return false;
    }

    // Gather all the edges that intersect these cells and sort them.
    int shapeId = shape.id;
    for (int i = 0; i < cells.size(); ++i) {
      S2ShapeIndex.Cell cell = cells.get(i);
      S2ClippedShape clipped = cell.findClipped(shapeId);
      if (clipped == null) {
        continue;
      }
      edges.ensureCapacity(edges.size() + clipped.numEdges());
      for (int j = 0; j < clipped.numEdges(); ++j) {
        edges.add(clipped.edge(j));
      }
    }
    
    if (cells.size() > 1) {
      // Sort and deduplicate the edge list.
      Collections.sort(edges);
      Iterator<Integer> it = edges.iterator();
      Integer last = null;
      while (it.hasNext()) {
        Integer edge = it.next();
        if (Objects.equals(edge, last)) {
          it.remove();
        } else {
          last = edge;
        }
      }
    }
    return !edges.isEmpty();
  }

  /**
   * Given a query edge AB, sets {@code edgeMap} to a map from the indexed shapes to a superset of
   * the edges for each shape that intersect AB. Returns false if no shapes intersect AB.
   */
  public boolean getCandidates(S2Point a, S2Point b, Map<S2Shape, List<Integer>> edgeMap) {
    // If there are only a few edges then it's faster to use brute force. We only bother with this
    // optimization when there is a single shape, since then we can also use some tricks to avoid
    // reallocating the edge map.
    if (index.numShapes() == 1) {
      // Typically this method is called many times, so it is worth checking whether the edge map is
      // empty or already consists of a single entry for this shape, and skip clearing 'edgeMap' in
      // that case.
      S2ShapeIndex.IdShape shape = index.shape(0);
      ArrayList<Integer> edges;
      if (edgeMap.size() != 1) {
        // 'edgeMap' must have been used to query some other S2ShapeIndex, so we need to clear its
        // current contents.
        edgeMap.clear();
        ArrayList<Integer> l = Lists.newArrayList();
        edgeMap.put(shape.shape, l);
        edges = l;
      } else if (edgeMap.get(shape.shape) == null) {
        ArrayList<Integer> l = Lists.newArrayList();
        edgeMap.put(shape.shape, l);
        edges = l;
      } else {
        edges = (ArrayList<Integer>) edgeMap.get(shape.shape);
      }
      // Note that we leave 'edgeMap' non-empty even if there are no candidates (i.e., there is a
      // single entry with an empty set of edges). This is an advantage for efficiency since it
      // avoids memory reallocation.
      return getCandidates(a, b, shape, edges);
    }
    
    // Compute the set of index cells intersected by the query edge.
    getCells(a, b);
    edgeMap.clear();
    if (cells.isEmpty()) {
      return false;
    }
    
    // Gather all the edges that intersect those cells.
    for (int i = 0; i < cells.size(); ++i) {
      S2ShapeIndex.Cell cell = cells.get(i);
      for (int s = 0; s < cell.numShapes(); ++s) {
        S2ClippedShape clipped = cell.clipped(s);
        S2Shape shape = index.shape(clipped.shapeId()).shape;
        List<Integer> edges = edgeMap.get(shape);
        if (edges == null) {
          List<Integer> edgeList = Lists.newArrayList();
          edgeMap.put(shape, edgeList);
          edges = edgeList;
        }
        for (int j = 0; j < clipped.numEdges(); ++j) {
          edges.add(clipped.edge(j));
        }
      }
    }
    
    // Sort the edges that were gathered. Note that if there is only one cell, then each clipped
    // shape's edges are added to edgeMap in sorted order already, because the clipped shape's edges
    // are stored in sorted order.
    if (cells.size() > 1) {
      for (List<Integer> edgeList : edgeMap.values()) {
        // Sort and de-duplicate 'edgeList'.
        Collections.sort(edgeList);
        Iterator<Integer> it = edgeList.iterator();
        Integer last = null;
        while (it.hasNext()) {
          Integer edge = it.next();
          if (Objects.equals(edge, last)) {
            it.remove();
          } else {
            last = edge;
          }
        }
      }
    }
    return !edgeMap.isEmpty();
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
   * Adds all cells to {@code cells} that might intersect the query edge from {@code a} to {@code b}
   * and the cell {@code root}. The {@code aVector} and {@code bVector} parameters are cached R2
   * versions of the [A, B] edge projected onto the same cube face as {@code root}.
   */
  @VisibleForTesting
  boolean getCells(S2PaddedCell root, S2Point a, R2Vector aVector, S2Point b, R2Vector bVector,
      List<S2ShapeIndex.Cell> cells) {
    this.cells.clear();
    if (S2EdgeUtil.clipToFace(a,  b,  root.id().face(), aVector, bVector)) {
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
    childBounds[0] = edgeBound;
    childBounds[1] = edgeBound;
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
}