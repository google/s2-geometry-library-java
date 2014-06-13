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
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2EdgeUtil.WedgeRelation;
import com.google.common.geometry.S2Error.Code;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.CellIterator;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;

import java.util.AbstractList;
import java.util.List;

/**
 * Utilities for working with S2Shape.
 */
@GwtCompatible
strictfp class S2ShapeUtil {
  /** Utility methods only. */
  private S2ShapeUtil() {}

  /**
   * S2EdgeVectorShape is an S2Shape representing a set of unrelated edges. It contains no area and
   * has no interior. Although it implements List<S2Edge>, only the {@link #add(S2Point, S2Point)}
   * method can mutate the list of edges.
   *
   * <p>It is mainly used for testing, but it can also be useful if you have, say, a collection of
   * polylines and don't care about memory efficiency (since this class would store most of the
   * vertices twice.) If the vertices are already stored somewhere else, you would be better off
   * writing your own subclass of S2Shape that points to the existing vertex data rather than
   * copying it.
   */
  static class S2EdgeVectorShape extends AbstractList<S2Edge> implements S2Shape {
    private final List<S2Edge> edges = Lists.newArrayList();

    /** Default constructor creates a vector with no edges. */
    public S2EdgeVectorShape() {}

    /** Convenience constructor for creating a vector of length 1. */
    public S2EdgeVectorShape(S2Point a, S2Point b) {
      add(a, b);
    }

    /** Add an edge to the vector. */
    public void add(S2Point a, S2Point b) {
      Preconditions.checkArgument(!a.equalsPoint(b));
      edges.add(new S2Edge(a, b));
    }

    @Override
    public void getEdge(int index, MutableEdge result) {
      S2Edge edge = edges.get(index);
      result.set(edge.getStart(), edge.getEnd());
    }

    @Override
    public boolean hasInterior() {
      return false;
    }

    @Override
    public boolean containsOrigin() {
      return false;
    }

    @Override
    public int numEdges() {
      return edges.size();
    }

    @Override
    public S2Edge get(int index) {
      return edges.get(index);
    }

    @Override
    public int size() {
      return edges.size();
    }
  }

  /** A simple implementation of the MutableEdge interface. */
  static class Edge implements MutableEdge {
    /**
     * Endpoints of this edge last set by passing this instance to
     * {@link S2Shape#getEdge(int, MutableEdge)}.
     */
    S2Point a, b;

    @Override
    public S2Point getStart() {
      return a;
    }

    @Override
    public S2Point getEnd() {
      return b;
    }

    @Override
    public void set(S2Point start, S2Point end) {
      this.a = start;
      this.b = end;
    }
  }

  /**
   * Given an S2ShapeIndex containing a single loop, return true if the loop has a self-intersection
   * (including duplicate vertices) and set "error" to a human-readable error message. Otherwise
   * return false and leave "error" unchanged.
   */
  static boolean findSelfIntersection(S2ShapeIndex index, S2Loop loop, S2Error error) {
    Preconditions.checkArgument(1 == index.shapes.size());
    for (CellIterator it = index.iterator(); !it.done(); it.next()) {
      if (findSelfIntersection(it.cell().clipped(0), loop, error)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Given an S2ShapeIndex containing a set of loops, return true if any loop has a
   * self-intersection (including duplicate vertices) or crosses any other loop (including vertex
   * crossings and duplicate edges) and set "error" to a human-readable error message. Otherwise
   * return false and leave "error" unchanged.
   */
  static boolean findAnyCrossing(S2ShapeIndex index, List<S2Loop> loops, S2Error error) {
    for (CellIterator it = index.iterator(); !it.done(); it.next()) {
      if (findSelfIntersection(loops, it.cell(), error)) {
        return true;
      }
      if (it.cell().numShapes() >= 2 && findLoopCrossing(loops, it.cell(), error)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Test for crossings between all edge pairs that do not share a vertex. This means that (a) the
   * loop edge indices must differ by 2 or more, and (b) the pair cannot consist of the first and
   * last loop edges. Part (b) is worthwhile in the case of very small loops; e.g. it reduces the
   * number of crossing tests in a loop with four edges from two to one.
   */
  static boolean findSelfIntersection(S2ClippedShape aClipped, S2Loop aLoop, S2Error error) {
    int aNumClipped = aClipped.numEdges();
    for (int i = 0; i < aNumClipped - 1; i++) {
      int ai = aClipped.edge(i);
      int j = i + 1;
      if (aClipped.edge(j) == ai + 1) {
        // Adjacent edges
        if (++j >= aNumClipped) {
          continue;
        }
      }

      EdgeCrosser crosser = new EdgeCrosser(aLoop.vertex(ai), aLoop.vertex(ai + 1));
      for (int ajPrev = -2; j < aNumClipped; j++) {
        int aj = aClipped.edge(j);
        if (aj - ai == aLoop.numVertices() - 1) {
          // First and last edges
          break;
        }
        if (aj != ajPrev + 1) {
          crosser.restartAt(aLoop.vertex(aj));
        }
        ajPrev = aj;

        // This test also catches duplicate vertices.
        int crossing = crosser.robustCrossing(aLoop.vertex(aj + 1));
        if (crossing < 0) {
          continue;
        }
        if (crossing == 0) {
          error.init(Code.DUPLICATE_VERTICES, "Edge %d has duplicate vertex with edge %d", ai, aj);
        } else {
          error.init(Code.LOOP_SELF_INTERSECTION, "Edge %d crosses edge %d", ai, aj);
        }

        return true;
      }
    }

    return false;
  }

  /**
   * Returns true if any of the given loops has a self-intersection (including a duplicate vertex),
   * and set "error" to a human-readable error message. Otherwise return false and leave "error"
   * unchanged. All tests are limited to edges that intersect the given cell.
   */
  static boolean findSelfIntersection(List<S2Loop> loops, Cell cell, S2Error error) {
    for (int a = 0; a < cell.numShapes(); a++) {
      S2ClippedShape aClipped = cell.clipped(a);
      if (findSelfIntersection(aClipped, loops.get(aClipped.shapeId()), error)) {
        error.init(error.code(), "Loop %d: %s", aClipped.shapeId(), error.text());
        return true;
      }
    }
    return false;
  }

  /**
   * Given two loop edges for which RobustCrossing returned a non-negative result "crossing",
   * returns true if there is a crossing and sets "error" to a human-readable error message,
   * otherwise returns false.
   */
  static boolean getCrossingError(
      S2Loop aLoop, int aShapeId, int ai,
      S2Loop bLoop, int bShapeId, int bj,
      int crossing, S2Error error) {
    if (crossing > 0) {
      error.init(Code.POLYGON_LOOPS_CROSS,
          "Loop %d edge %d crosses loop %d edge %d", aShapeId, ai, bShapeId, bj);
      return true;
    }

    // Loops are not allowed to share edges or cross at vertices.  We only need to check this once
    // per edge pair, so we also require that the two edges have the same end vertex.  (This is only
    // valid because we are iterating over all the cells in the index.)
    if (aLoop.vertex(ai + 1).equalsPoint(bLoop.vertex(bj + 1))) {
      if (aLoop.vertex(ai).equalsPoint(bLoop.vertex(bj))
          || aLoop.vertex(ai).equalsPoint(bLoop.vertex(bj + 2))) {
        // The second edge index is sometimes off by one, hence "near".
        error.init(Code.POLYGON_LOOPS_SHARE_EDGE,
            "Loop %d edge %d has duplicate near loop %d edge %d",
            aShapeId, ai, bShapeId, bj);
        return true;
      }

      // Note that we don't need to maintain any state regarding loop crossings because duplicate
      // edges are not allowed.
      if (WedgeRelation.WEDGE_PROPERLY_OVERLAPS == S2EdgeUtil.getWedgeRelation(
          aLoop.vertex(ai), aLoop.vertex(ai + 1), aLoop.vertex(ai + 2),
          bLoop.vertex(bj), bLoop.vertex(bj + 2))) {
        error.init(Code.POLYGON_LOOPS_CROSS,
            "Loop %d edge %d crosses loop %d edge %d", aShapeId, ai, bShapeId, bj);
        return true;
      }
    }

    return false;
  }

  /**
   * Returns true if any of the given loops crosses a different loop (including vertex crossings) or
   * two loops share a common edge, and sets "error" to a human-readable error message. Otherwise
   * return false and leaves "error" unchanged. All tests are limited to edges that intersect the
   * given cell.
   */
  static boolean findLoopCrossing(List<S2Loop> loops, Cell cell, S2Error error) {
    // Possible optimization:
    // Sort the ClippedShapes by edge count to reduce the number of calls to S2.robustCCW.  If n is
    // the total number of shapes in the cell, n_i is the number of edges in shape i, and c_i is the
    // number of continuous chains formed by these edges, the total number of calls is
    //
    //   sum(n_i * (1 + c_j + n_j), i=0..n-2, j=i+1..n-1)
    //
    // So for example if n=2, shape 0 has one chain of 1 edge, and shape 1 has one chain of 8 edges,
    // the number of calls to RobustCCW is 1*10=10 if the shapes are sorted by edge count, and
    // 8*3=24 otherwise.
    for (int a = 0; a < cell.numShapes() - 1; a++) {
      S2ClippedShape aClipped = cell.clipped(a);
      S2Loop aLoop = loops.get(aClipped.shapeId());
      int aNumClipped = aClipped.numEdges();
      for (int i = 0; i < aNumClipped; i++) {
        int ai = aClipped.edge(i);
        EdgeCrosser crosser = new EdgeCrosser(aLoop.vertex(ai), aLoop.vertex(ai + 1));
        for (int b = a + 1; b < cell.numShapes(); b++) {
          S2ClippedShape bClipped = cell.clipped(b);
          S2Loop bLoop = loops.get(bClipped.shapeId());
          int bjPrev = -2;
          int bNumClipped = bClipped.numEdges();
          for (int j = 0; j < bNumClipped; j++) {
            int bj = bClipped.edge(j);
            if (bj != bjPrev + 1) {
              crosser.restartAt(bLoop.vertex(bj));
            }
            bjPrev = bj;
            int crossing = crosser.robustCrossing(bLoop.vertex(bj + 1));
            if (crossing < 0) {
              // No crossing
              continue;
            }
            if (getCrossingError(
                aLoop, aClipped.shapeId(), ai,
                bLoop, bClipped.shapeId(), bj,
                crossing, error)) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }
}
