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
import com.google.common.base.Supplier;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Sets;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2EdgeUtil.WedgeRelation;
import com.google.common.geometry.S2Error.Code;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2Shape.ReferencePoint;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/** Utilities for working with S2Shape. */
@GwtCompatible
public strictfp class S2ShapeUtil {
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

    /** Adds an edge to the vector. */
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
    public int numChains() {
      return edges.size();
    }

    @Override
    public int getChainStart(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      return chainId;
    }

    @Override
    public int getChainLength(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      return 1;
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge result) {
      Preconditions.checkElementIndex(offset, getChainLength(chainId));
      getEdge(chainId, result);
    }

    @Override
    public S2Point getChainVertex(int chainId, int edgeOffset) {
      Preconditions.checkElementIndex(edgeOffset, getChainLength(chainId));
      S2Edge edge = edges.get(chainId);
      return edgeOffset == 0 ? edge.getStart() : edge.getEnd();
    }

    @Override
    public int dimension() {
      return 1;
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

  /**
   * Given an S2ShapeIndex containing a single loop, return true if the loop has a self-intersection
   * (including duplicate vertices) and set "error" to a human-readable error message. Otherwise
   * return false and leave "error" unchanged.
   */
  static boolean findSelfIntersection(S2ShapeIndex index, S2Loop loop, S2Error error) {
    Preconditions.checkArgument(1 == index.shapes.size());
    for (S2Iterator<S2ShapeIndex.Cell> it = index.iterator(); !it.done(); it.next()) {
      if (findSelfIntersection(it.entry().clipped(0), loop, error)) {
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
    for (S2Iterator<S2ShapeIndex.Cell> it = index.iterator(); !it.done(); it.next()) {
      if (findSelfIntersection(loops, it.entry(), error)) {
        return true;
      }
      if (it.entry().numShapes() >= 2 && findLoopCrossing(loops, it.entry(), error)) {
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
   * Returns the index of {@code shape} in {@code shapes} as {@link List#indexOf(Object)}, but using
   * identity instead of equality to honor the semantics of S2ShapeIndex (where adding two S2Loops
   * that are equal but not the same instance is treated as adding two separate shapes, with
   * distinct shape IDs.)
   */
  static int indexOf(List<? extends S2Shape> shapes, S2Shape shape) {
    for (int i = 0; i < shapes.size(); i++) {
      if (shapes.get(i) == shape) {
        return i;
      }
    }
    return -1;
  }

  /** A filter of indexes. */
  public interface IntPredicate {
    /** Returns whether the given value tests as true. */
    boolean test(int index);
  }

  /** Returns the lowest index in the range {@code [low, high)} not smaller than a target. */
  public static int lowerBound(int low, int high, IntPredicate targetIsGreater) {
    while (low < high) {
      int middle = low + (high - low) / 2;
      if (targetIsGreater.test(middle)) {
        low = middle + 1;
      } else {
        high = middle;
      }
    }
    return low;
  }

  /** Returns the lowest index in the range {@code [low, high)} greater than a target. */
  public static int upperBound(int low, int high, IntPredicate targetIsSmaller) {
    while (low < high) {
      int middle = low + (high - low) / 2;
      if (targetIsSmaller.test(middle)) {
        high = middle;
      } else {
        low = middle + 1;
      }
    }
    return low;
  }

  /**
   * Returns true if any of the given loops has a self-intersection (including a duplicate vertex),
   * and set "error" to a human-readable error message. Otherwise return false and leave "error"
   * unchanged. All tests are limited to edges that intersect the given cell.
   */
  static boolean findSelfIntersection(List<S2Loop> loops, Cell cell, S2Error error) {
    for (int a = 0; a < cell.numShapes(); a++) {
      S2ClippedShape aClipped = cell.clipped(a);
      S2Loop loop = (S2Loop) aClipped.shape();
      if (findSelfIntersection(aClipped, loop, error)) {
        error.init(error.code(), "Loop %d: %s", indexOf(loops, loop), error.text());
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
      List<S2Loop> loops, S2Loop aLoop, int ai, S2Loop bLoop, int bj, int crossing, S2Error error) {
    if (crossing > 0) {
      error.init(
          Code.POLYGON_LOOPS_CROSS,
          "Loop %d edge %d crosses loop %d edge %d",
          indexOf(loops, aLoop),
          ai,
          indexOf(loops, bLoop),
          bj);
      return true;
    }

    // Loops are not allowed to share edges or cross at vertices.  We only need to check this once
    // per edge pair, so we also require that the two edges have the same end vertex.  (This is only
    // valid because we are iterating over all the cells in the index.)
    if (aLoop.vertex(ai + 1).equalsPoint(bLoop.vertex(bj + 1))) {
      if (aLoop.vertex(ai).equalsPoint(bLoop.vertex(bj))
          || aLoop.vertex(ai).equalsPoint(bLoop.vertex(bj + 2))) {
        // The second edge index is sometimes off by one, hence "near".
        error.init(
            Code.POLYGON_LOOPS_SHARE_EDGE,
            "Loop %d edge %d has duplicate near loop %d edge %d",
            indexOf(loops, aLoop),
            ai,
            indexOf(loops, bLoop),
            bj);
        return true;
      }

      // Note that we don't need to maintain any state regarding loop crossings because duplicate
      // edges are not allowed.
      if (WedgeRelation.WEDGE_PROPERLY_OVERLAPS
          == S2EdgeUtil.getWedgeRelation(
              aLoop.vertex(ai),
              aLoop.vertex(ai + 1),
              aLoop.vertex(ai + 2),
              bLoop.vertex(bj),
              bLoop.vertex(bj + 2))) {
        error.init(
            Code.POLYGON_LOOPS_CROSS,
            "Loop %d edge %d crosses loop %d edge %d",
            indexOf(loops, aLoop),
            ai,
            indexOf(loops, bLoop),
            bj);
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
    // Sort the ClippedShapes by edge count to reduce the number of calls to S2Predicates.sign().
    // If n is the total number of shapes in the cell, n_i is the number of edges in shape i, and
    // c_i is the number of continuous chains formed by these edges, the total number of calls is
    //
    //   sum(n_i * (1 + c_j + n_j), i=0..n-2, j=i+1..n-1)
    //
    // So for example if n=2, shape 0 has one chain of 1 edge, and shape 1 has one chain of 8 edges,
    // the number of calls to S2Predicates.sign() is 1*10=10 if the shapes are sorted by edge count,
    // and 8*3=24 otherwise.
    for (int a = 0; a < cell.numShapes() - 1; a++) {
      S2ClippedShape aClipped = cell.clipped(a);
      S2Loop aLoop = (S2Loop) aClipped.shape();
      int aNumClipped = aClipped.numEdges();
      for (int i = 0; i < aNumClipped; i++) {
        int ai = aClipped.edge(i);
        EdgeCrosser crosser = new EdgeCrosser(aLoop.vertex(ai), aLoop.vertex(ai + 1));
        for (int b = a + 1; b < cell.numShapes(); b++) {
          S2ClippedShape bClipped = cell.clipped(b);
          S2Loop bLoop = (S2Loop) bClipped.shape();
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
            if (getCrossingError(loops, aLoop, ai, bLoop, bj, crossing, error)) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  /**
   * Returns true if all methods of the two S2Shapes return identical results, except for id() and
   * typeTag(). Also returns true if both instances are null.
   */
  public static boolean equals(S2Shape a, S2Shape b) {
    // Check null on either side.
    if (a == null) {
      return b == null;
    } else if (b == null) {
      return a == null;
    }

    // Check geometry type properties of the shapes.
    if (a.hasInterior() != b.hasInterior()) {
      return false;
    }
    if (a.hasInterior() && (a.containsOrigin() != b.containsOrigin())) {
      return false;
    }
    if (a.dimension() != b.dimension()) {
      return false;
    }

    // Check chain methods.
    if (a.numChains() != b.numChains()) {
      return false;
    }
    for (int i = 0; i < a.numChains(); i++) {
      if (a.getChainStart(i) != b.getChainStart(i)) {
        return false;
      }
      if (a.getChainLength(i) != b.getChainLength(i)) {
        return false;
      }
    }

    // Check edge methods.
    if (a.numEdges() != b.numEdges()) {
      return false;
    }
    MutableEdge edge = new MutableEdge();
    for (int i = 0; i < a.numEdges(); i++) {
      a.getEdge(i, edge);
      S2Point p = edge.a;
      S2Point q = edge.b;
      b.getEdge(i, edge);
      if (!p.equalsPoint(edge.a) || !q.equalsPoint(edge.b)) {
        return false;
      }
    }

    return true;
  }

  /**
   * Returns true if the lists 'a' and 'b' have identical shapes according to {@link
   * #equals(S2Shape, S2Shape)}.
   */
  public static boolean equals(List<S2Shape> a, List<S2Shape> b) {
    if (a.size() != b.size()) {
      return false;
    }
    for (int i = 0; i < a.size(); i++) {
      if (!equals(a.get(i), b.get(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if the clipped shapes 'a' and 'b' have identical edge offsets.
   *
   * <p>This method does not check that {@code a.shape()} and {@code b.shape()} are equal.
   */
  public static boolean equals(S2ClippedShape a, S2ClippedShape b) {
    if (a.containsCenter() != b.containsCenter()) {
      return false;
    }
    if (a.numEdges() != b.numEdges()) {
      return false;
    }
    for (int i = 0; i < a.numEdges(); i++) {
      if (a.edge(i) != b.edge(i)) {
        return false;
      }
    }
    return true;
  }

  /** Returns true if the index cells 'a' and 'b' contain identical contents. */
  public static boolean equals(S2ShapeIndex.Cell a, S2ShapeIndex.Cell b) {
    if (a.id() != b.id()) {
      return false;
    }
    if (a.numShapes() != b.numShapes()) {
      return false;
    }
    for (int i = 0; i < a.numShapes(); i++) {
      if (!equals(a.clipped(i), b.clipped(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if all methods of the two S2ShapeIndex values return identical results, including
   * all the S2Shapes in both indexes.
   */
  public static boolean equals(S2ShapeIndex a, S2ShapeIndex b) {
    // Check that both indexes have identical shapes.
    if (!equals(a.getShapes(), b.getShapes())) {
      return false;
    }

    // In order to test that the shapes referenced by all clipped shapes of all cells in each index
    // are equal, we build a Multimap<S2Shape,S2Shape> where each S2Shape is distinguished by its
    // identity hash code.
    Multimap<S2Shape, S2Shape> shapes =
        Multimaps.newSetMultimap(Maps.newIdentityHashMap(), Sets::newIdentityHashSet);
    for (int i = 0; i < a.getShapes().size(); i++) {
      shapes.put(a.getShapes().get(i), b.getShapes().get(i));
    }

    // Check that both indexes have identical cell contents.
    S2Iterator<S2ShapeIndex.Cell> aIt = a.iterator();
    S2Iterator<S2ShapeIndex.Cell> bIt = b.iterator();
    for (; !aIt.done(); aIt.next(), bIt.next()) {
      if (bIt.done()) {
        return false;
      }
      if (!equals(aIt.entry(), bIt.entry())) {
        return false;
      }
      // Check that each clipped shape references the same shape.
      for (int i = 0; i < aIt.entry().numShapes(); i++) {
        if (!shapes.containsEntry(aIt.entry().clipped(i).shape(), bIt.entry().clipped(i).shape())) {
          return false;
        }
      }
    }
    if (!bIt.done()) {
      return false;
    }

    return true;
  }

  /** Compares edges by start point, and then by end point. */
  private static final Comparator<S2Edge> EDGE_ORDER = (e1, e2) -> {
    int result = e1.getStart().compareTo(e2.getStart());
    if (result != 0) {
      return result;
    } else {
      return e1.getEnd().compareTo(e2.getEnd());
    }
  };

  /**
   * Returns true if the given shape contains the given point. Most clients should not use this
   * method, since its running time is linear in the number of shape edges. Instead clients should
   * create an S2ShapeIndex and use {@link S2ContainsPointQuery}, since that strategy is much more
   * efficient when many points need to be tested.
   *
   * <p>Polygon boundaries are treated as being semi-open. See {@link
   * S2ContainsPointQuery.S2VertexModel} for other options.
   */
  public static boolean containsBruteForce(S2Shape shape, S2Point point) {
    if (shape.dimension() < 2) {
      return false;
    }
    ReferencePoint refPoint = shape.getReferencePoint();
    if (refPoint.equalsPoint(point)) {
      return refPoint.contained();
    }
    EdgeCrosser crosser = new EdgeCrosser(refPoint, point);
    boolean inside = refPoint.contained();
    MutableEdge edge = new MutableEdge();
    for (int i = 0; i < shape.numEdges(); i++) {
      shape.getEdge(i, edge);
      inside ^= crosser.edgeOrVertexCrossing(edge.getStart(), edge.getEnd());
    }
    return inside;
  }

  /**
   * This is a helper function for implementing S2Shape.getReferencePoint().
   *
   * <p>Given a shape consisting of closed polygonal loops, the interior of the shape is defined as
   * the region to the left of all edges (which must be oriented consistently). This function then
   * chooses an arbitrary point and returns true if that point is contained by the shape.
   *
   * <p>Unlike S2Loop and S2Polygon, this method allows duplicate vertices and edges, which requires
   * some extra care with definitions. The rule that we apply is that an edge and its reverse edge
   * "cancel" each other: the result is the same as if that edge pair were not present. Therefore
   * shapes that consist only of degenerate loop(s) are either empty or full; by convention, the
   * shape is considered full if and only if it contains an empty loop (see S2LaxPolygonShape for
   * details).
   *
   * <p>Determining whether a loop on the sphere contains a point is harder than the corresponding
   * problem in 2D plane geometry. It cannot be implemented just by counting edge crossings because
   * there is no such thing as a "point at infinity" that is guaranteed to be outside the loop.
   */
  public static ReferencePoint getReferencePoint(S2Shape shape) {
    assert shape.dimension() == 2;

    if (shape.numEdges() == 0) {
      // A shape with no edges is defined to be full if and only if it contains at least one chain.
      return ReferencePoint.create(shape.numChains() > 0);
    }

    // Define a "matched" edge as one that can be paired with a corresponding reversed edge.
    // Define a vertex as "balanced" if all of its edges are matched. In order to determine
    // containment, we must find an unbalanced vertex.  Often every vertex is unbalanced, so we
    // start by trying an arbitrary vertex.
    MutableEdge edge = new MutableEdge();
    shape.getEdge(0, edge);
    Boolean result;
    if (null != (result = getReferencePointAtVertex(shape, edge.a))) {
      return ReferencePoint.create(edge.a, result);
    }

    // That didn't work, so now we do some extra work to find an unbalanced vertex (if any).
    // Essentially we gather a list of edges and a list of reversed edges, and then sort them.
    // The first edge that appears in one list but not the other is guaranteed to be unmatched.
    int n = shape.numEdges();
    List<S2Edge> fwdEdges = new ArrayList<>(n);
    List<S2Edge> revEdges = new ArrayList<>(n);
    for (int i = 0; i < n; ++i) {
      shape.getEdge(i, edge);
      fwdEdges.add(new S2Edge(edge.a, edge.b));
      revEdges.add(new S2Edge(edge.b, edge.a));
    }
    Collections.sort(fwdEdges, EDGE_ORDER);
    Collections.sort(revEdges, EDGE_ORDER);
    for (int i = 0; i < n; i++) {
      S2Edge fwd = fwdEdges.get(i);
      S2Edge rev = revEdges.get(i);
      S2Point v;
      int cmp = EDGE_ORDER.compare(fwd, rev);
      if (cmp < 0) {
        // fwd is unmatched
        v = fwd.getStart();
      } else if (cmp > 0) {
        // rev is unmatched
        v = rev.getStart();
      } else {
        // This edge is matched, so move on to the next one.
        continue;
      }
      // We have an unbalanced vertex, so reference it and return.
      return ReferencePoint.create(v, getReferencePointAtVertex(shape, v));
    }

    // All vertices are balanced, so this polygon is either empty or full except for degeneracies.
    // By convention it is defined to be full if it contains any chain with no edges.
    for (int i = 0; i < shape.numChains(); i++) {
      if (shape.getChainLength(i) == 0) {
        return ReferencePoint.create(true);
      }
    }
    return ReferencePoint.create(false);
  }

  /**
   * Returns null if 'vtest' is balanced (see definition above), otherwise 'vtest' is unbalanced and
   * the return value indicates whether it is contained by 'shape'.
   */
  private static Boolean getReferencePointAtVertex(S2Shape shape, S2Point vtest) {
    // Let P be an unbalanced vertex.  Vertex P is defined to be inside the region if the region
    // contains a particular direction vector starting from P, namely the direction of
    // S2.ortho(target). This can be calculated using S2ContainsVertexQuery.
    S2ContainsVertexQuery query = new S2ContainsVertexQuery(vtest);
    MutableEdge edge = new MutableEdge();
    int n = shape.numEdges();
    for (int e = 0; e < n; ++e) {
      shape.getEdge(e, edge);
      if (vtest.equalsPoint(edge.a)) {
        query.addOutgoing(edge.b);
      }
      if (vtest.equalsPoint(edge.b)) {
        query.addIncoming(edge.a);
      }
    }
    int containsSign = query.containsSign();
    if (containsSign == 0) {
      // There are no unmatched edges incident to this vertex.
      return null;
    } else {
      return containsSign > 0;
    }
  }

  /**
   * Returns a multimap of {@link S2Shape} from {@code index} to the shape's ID (i.e., its position
   * within {@code index.shapes}).
   */
  static Multimap<S2Shape, Integer> shapeToShapeId(S2ShapeIndex index) {
    Multimap<S2Shape, Integer> shapeToShapeId =
        Multimaps.newListMultimap(Maps.newIdentityHashMap(), Lists::newArrayList);
    for (S2Shape shape : index.shapes) {
      shapeToShapeId.put(shape, shapeToShapeId.size());
    }
    return shapeToShapeId;
  }

  /** A consumer of triangles. Implementations may sum area, turning angle, etc. */
  public interface TriangleConsumer {
    void accept(S2Point a, S2Point b, S2Point c);
  }

  // It is surprisingly difficult to compute the area of a loop robustly. The main issues are
  //
  // (1) whether degenerate loops are considered to be CCW or not
  //     i.e., whether their area is close to 0 or 4*Pi, and
  // (2) computing the areas of small loops with good relative accuracy.
  //
  // With respect to degeneracies, we would like getArea() to be consistent with
  // S2Loop.contains(S2Point), such that loops containing many points should have large areas, and
  // loops that contain few points should have small areas.
  //
  // For example, if a degenerate triangle is considered CCW according to S2Predicates.sign(), then
  // it will contain very few points and its area should be approximately zero. On the other hand if
  // it is considered clockwise, then it will contain virtually all points and so its area should
  // be approximately 4*Pi.
  //
  // More precisely, let U be the set of S2Points for which S2.isUnitLength() is true, let P(U) be
  // the projection of those points onto the mathematical unit sphere, and let V(P(U)) be the
  // Voronoi diagram of the projected points. Then for every loop x, we would like getArea() to
  // approximately equal the sum of the areas of the Voronoi regions of the points p for which
  // x.contains(p) is true.
  //
  // The second issue is that we want to compute the area of small loops accurately. This requires
  // having good relative precision rather than good absolute precision. For example, if the area
  // of a loop is 1e-12 and the error is 1e-15, then the area only has 3 digits of accuracy. (For
  // reference, 1e-12 is about 40 square meters on the surface of the earth.) We would like to
  // have good relative accuracy even for small loops.
  //
  // To achieve these goals, we combine two different methods of computing the area. This first
  // method is based on the Gauss-Bonnet theorem, which says that the area enclosed by the loop
  // equals 2*Pi minus the total geodesic curvature of the loop (i.e., the sum of the "turning
  // angles" at all the loop vertices). The big advantage of this method is that as long as we
  // use S2Predicates.sign() to compute the turning angle at each vertex, then degeneracies are
  // always handled correctly. In other words, if a degenerate loop is CCW according to the symbolic
  // perturbations used by S2Predicates.sign(), then its turning angle will be approximately 2*Pi.
  //
  // The disadvantage of the Gauss-Bonnet method is that its absolute error is about 2e-15 times
  // the number of vertices (see S2Loop.getTurningAngleMaxError). So, it cannot compute the area
  // of small loops accurately.
  //
  // The second method is based on splitting the loop into triangles and summing the area of each
  // triangle. To avoid the difficulty and expense of decomposing the loop into a union of non-
  // overlapping triangles, instead we compute a signed sum over triangles that may overlap (see
  // the comments for S2Loop.visitSurfaceIntegral). The advantage of this method is that the area
  // of each triangle can be computed with much better relative accuracy (using l'Huilier's
  // theorem). The disadvantage is that the result is a signed area: CCW loops may yield a small
  // positive value, while CW loops may yield a small negative value (which is converted to a
  // positive area by adding 4*Pi). This means that small errors in computing the signed area may
  // translate into a very large error in the result (if the sign of the sum is incorrect).
  //
  // So, our strategy is to combine these two methods as follows. First we compute the area using
  // the "signed sum over triangles" approach (since it is generally more accurate). We also
  // estimate the maximum error in this result. If the signed area is too close to zero (i.e.,
  // zero is within the error bounds), then we double-check the sign of the result using the
  // Gauss-Bonnet method. (In fact we just call isNormalized(), which is based on this method.)
  // If the two methods disagree, we return either 0 or 4*Pi based on the result of
  // isNormalized(). Otherwise we return the area that we computed originally.

  /** A collector of the steradian area. */
  public static final class AreaMeasure implements TriangleConsumer {
    private double area;

    @Override
    public void accept(S2Point a, S2Point b, S2Point c) {
      area += S2.signedArea(a, b, c);
    }

    /** Returns the area. Only call after all triangles have been consumed. */
    public double value(int numVertices, Supplier<Boolean> isNormalized) {
      // TODO(user): This error estimate is very approximate.  There are two issues:
      //
      // (1) signedArea needs some improvements to ensure that its error is actually never higher
      //     than girardArea, and
      // (2) although the number of triangles in the sum is typically N-2, in theory it could be as
      //     high as 2*N for pathological inputs.
      //
      // But in other respects this error bound is very conservative since it assumes that the
      // maximum error is achieved on every triangle.
      double maxError = S2.getTurningAngleMaxError(numVertices);

      // The signed area should be between approximately -4*Pi and 4*Pi.
      assert (Math.abs(area) <= 4 * S2.M_PI + maxError);
      if (area < 0) {
        // We have computed the negative of the area of the loop exterior.
        area += 4 * S2.M_PI;
      }
      area = Math.max(0.0, Math.min(4 * S2.M_PI, area));

      // If the area is close enough to zero or 4*Pi so that the loop orientation
      // is ambiguous, then we compute the loop orientation explicitly.
      if (area < maxError && !isNormalized.get()) {
        return 4 * S2.M_PI;
      } else if (area > (4 * S2.M_PI - maxError) && isNormalized.get()) {
        return 0.0;
      } else {
        return area;
      }
    }
  }

  /** A collector of the center of mass. */
  public static final class CentroidMeasure implements TriangleConsumer {
    private final double[] sum = new double[3];

    @Override
    public void accept(S2Point a, S2Point b, S2Point c) {
      S2Point centroid = S2.trueCentroid(a, b, c);
      sum[0] += centroid.x;
      sum[1] += centroid.y;
      sum[2] += centroid.z;
    }

    /** Returns the centroid. Only call after all triangles have been consumed. */
    public S2Point value() {
      return new S2Point(sum[0], sum[1], sum[2]);
    }
  }

  /** A collector of both combined area and centroid values. */
  public static final class AreaCentroidMeasure implements TriangleConsumer {
    private final AreaMeasure area = new AreaMeasure();
    private final CentroidMeasure centroid = new CentroidMeasure();

    @Override
    public void accept(S2Point a, S2Point b, S2Point c) {
      area.accept(a, b, c);
      centroid.accept(a, b, c);
    }

    /** Returns the area and centroid. Only call after all triangles have been consumed. */
    public S2AreaCentroid value(int numVertices, Supplier<Boolean> isNormalized) {
      return new S2AreaCentroid(area.value(numVertices, isNormalized), centroid.value());
    }
  }

  /**
   * Visits the surface integral of the vertices, that is, a collection of oriented triangles,
   * possibly overlapping.
   *
   * <p>Let the sign of a triangle be +1 if it is CCW and -1 otherwise, and let the sign of a point
   * "x" be the sum of the signs of the triangles containing "x". Then the collection of triangles T
   * is chosen such that either:
   *
   * <ol>
   *   <li>Each point in the loop interior has sign +1, and sign 0 otherwise; or
   *   <li>Each point in the loop exterior has sign -1, and sign 0 otherwise.
   * </ol>
   *
   * <p>The triangles basically consist of a "fan" from vertex 0 to every loop edge that does not
   * include vertex 0. These triangles will always satisfy either (1) or (2).
   *
   * <p>However, what makes this a bit tricky is that spherical edges become numerically unstable as
   * their length approaches 180 degrees. Of course there is not much we can do if the loop itself
   * contains such edges, but we would like to make sure that all the triangle edges under our
   * control (i.e., the non-loop edges) are stable. For example, consider a loop around the equator
   * consisting of four equally spaced points. This is a well-defined loop, but we cannot just split
   * it into two triangles by connecting vertex 0 to vertex 2.
   *
   * <p>We handle this type of situation by moving the origin of the triangle fan whenever we are
   * about to create an unstable edge. We choose a new location for the origin such that all
   * relevant edges are stable. We also create extra triangles with the appropriate orientation so
   * that the sum of the triangle signs is still correct at every point.
   */
  public static void visitSurfaceIntegral(List<S2Point> vertices, TriangleConsumer consumer) {
    if (vertices.size() < 3) {
      return;
    }
    // The maximum length of an edge for it to be considered numerically stable. The exact value is
    // fairly arbitrary since it depends on the stability of the consumer's processing. The value
    // below is quite conservative but could be reduced further if desired.
    final double maxLength = S2.M_PI - 1e-5;
    S2Point v0 = vertices.get(0);
    S2Point origin = v0;

    // Let V_i be vertices.get(i), let O be the current origin, and let length(A,B) be the length of
    // edge (A,B). At the start of each loop iteration, the "leading edge" of the triangle fan is
    // (O,V_i), and we want to extend the triangle fan so that the leading edge is (O,V_i+1).
    //
    // Invariants:
    //  1. length(O,V_i) < kMaxLength for all (i > 1).
    //  2. Either O == V_0, or O is approximately perpendicular to V_0.
    //  3. "sum" is the oriented integral of f over the area defined by
    //     (O, V_0, V_1, ..., V_i).
    int numVertices = vertices.size();
    for (int i = 1; i + 1 < numVertices; i++) {
      assert (i == 1 || origin.angle(vertices.get(i)) < maxLength);
      assert (origin.equals(v0) || Math.abs(origin.dotProd(v0)) < 1e-15);
      S2Point v1 = vertices.get(i);
      S2Point v2 = vertices.get(i + 1);
      if (v2.angle(origin) > maxLength) {
        // We are about to create an unstable edge, so choose a new origin O' for the triangle fan.
        S2Point oldOrigin = origin;
        if (origin.equalsPoint(v0)) {
          // The following point is well-separated from V_i and V_0 (and therefore V_i+1 as well).
          origin = S2Point.normalize(S2.robustCrossProd(v0, v1));
        } else if (v1.angle(v0) < maxLength) {
          // All edges of tri (O, V_0, V_i) are stable, so we can revert to using V_0 as the origin.
          origin = v0;
        } else {
          // (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are perpendicular.
          // Therefore V_0.crossProd(O) is approximately perpendicular to all of
          // {O, V_0, V_i, V_i+1}, and we can choose this point O' as the new origin.
          origin = S2Point.crossProd(v0, oldOrigin);

          // Advance the edge (V_0,O) to (V_0,O').
          consumer.accept(v0, oldOrigin, origin);
        }
        // Advance the edge (O,V_i) to (O',V_i).
        consumer.accept(oldOrigin, v1, origin);
      }
      // Advance the edge (O,V_i) to (O,V_i+1).
      consumer.accept(origin, v1, v2);
    }

    // If the origin is not V_0, we need to sum one more triangle.
    if (!origin.equalsPoint(v0)) {
      // Advance the edge (O,V_n-1) to (O,V_0).
      consumer.accept(origin, vertices.get(numVertices - 1), v0);
    }
  }
}
