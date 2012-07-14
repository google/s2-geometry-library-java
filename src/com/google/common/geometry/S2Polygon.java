/*
 * Copyright 2006 Google Inc.
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

import com.google.common.base.Preconditions;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultimap;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * An S2Polygon is an {@link S2Region} object that represents a polygon. A
 * polygon consists of zero or more {@link S2Loop}s representing "shells" and
 * "holes". All loops should be oriented CCW, i.e. the shell or hole is on the
 * left side of the loop. Loops may be specified in any order. A point is
 * defined to be inside the polygon if it is contained by an odd number of
 * loops.
 *
 * <p>Polygons have the following restrictions:
 * <ul>
 * <li>Loops may not cross, i.e. the boundary of a loop may not intersect both
 *     the interior and exterior of any other loop.
 * <li>Loops may not share edges, i.e. if a loop contains an edge AB, then no
 *     other loop may contain AB or BA.
 * <li>No loop may cover more than half the area of the sphere. This ensures
 *     that no loop properly contains the complement of any other loop, even if
 *     the loops are from different polygons. (Loops that represent exact
 *     hemispheres are allowed.)
 * </ul>
 * <p>Loops may share vertices, however no vertex may appear twice in a single
 * loop.
 *
 */
public final strictfp class S2Polygon implements S2Region, Comparable<S2Polygon> {
  private static final Logger log = Logger.getLogger(S2Polygon.class.getCanonicalName());

  private final List<S2Loop> loops = Lists.newArrayList();

  private S2LatLngRect bound = S2LatLngRect.empty();
  private boolean hasHoles = false;
  private int numVertices = 0;

  /**
   * Creates an empty polygon that should be initialized by calling
   * {@link #init}.
   */
  public S2Polygon() {}

  /**
   * Convenience constructor that calls {@link #init} with the given loops.
   * Clears the given list.
   */
  public S2Polygon(List<S2Loop> loops) {
    init(loops);
  }

  /**
   * Copy constructor.
   */
  public S2Polygon(S2Loop loop) {
    this.bound = loop.getRectBound();
    this.numVertices = loop.numVertices();

    loops.add(loop);
  }

  /**
   * Copy constructor.
   */
  public S2Polygon(S2Polygon src) {
    this.bound = src.getRectBound();
    this.hasHoles = src.hasHoles;
    this.numVertices = src.numVertices;

    for (int i = 0; i < src.numLoops(); ++i) {
      loops.add(new S2Loop(src.loop(i)));
    }
  }

  /**
   * Comparator (needed by Comparable interface). For two polygons to be
   * compared as equal:
   * <ul>
   * <li>They must have the same number of loops
   * <li>The loops must be ordered in the same way (this is guaranteed by the
   *     total ordering imposed by {@link #sortValueLoops})
   * <li>Loops must be logically equivalent (even if ordered with a different
   *     starting point, e.g. ABCD and BCDA).
   * </ul>
   */
  @Override
  public int compareTo(S2Polygon other) {
    // If number of loops differ, use that.
    if (this.numLoops() != other.numLoops()) {
      return this.numLoops() - other.numLoops();
    }
    for (int i = 0; i < this.numLoops(); ++i) {
      int compare = this.loops.get(i).compareTo(other.loops.get(i));
      if (compare != 0) {
        return compare;
      }
    }
    return 0;
  }

  /**
   * Initializes a polygon by taking ownership of the given loops and clearing
   * the given list. This method figures out the loop-nesting hierarchy and then
   * reorders the loops by following a preorder traversal. This implies that
   * each loop is immediately followed by its descendants in the nesting
   * hierarchy. (See also {@link #getParent} and {@link #getLastDescendant}.)
   */
  public void init(List<S2Loop> loops) {
    // assert isValid(loops);
    // assert (this.loops.isEmpty());

    Map<S2Loop, List<S2Loop>> loopMap = Maps.newHashMap();
    // Yes, a null key is valid. It is used here to refer to the root of the
    // loopMap
    loopMap.put(null, Lists.<S2Loop>newArrayList());

    for (S2Loop loop : loops) {
      insertLoop(loop, null, loopMap);
      this.numVertices += loop.numVertices();
    }
    loops.clear();

    // Sort all of the lists of loops; in this way we guarantee a total ordering
    // on loops in the polygon. Loops will be sorted by their natural ordering,
    // while also preserving the requirement that each loop is immediately
    // followed by its descendants in the nesting hierarchy.
    //
    // TODO(andriy): as per kirilll in CL 18750833 code review comments:
    // This should work for now, but I think it's possible to guarantee the
    // correct order inside insertLoop by searching for the correct position in
    // the children list before inserting.
    sortValueLoops(loopMap);

    // Reorder the loops in depth-first traversal order.
    // Starting at null == starting at the root
    initLoop(null, -1, loopMap);

    // TODO(dbeaumont): Add tests or preconditions for these asserts (here and elesewhere).
    // forall i != j : containsChild(loop(i), loop(j), loopMap) == loop(i).containsNested(loop(j)));

    // Compute the bounding rectangle of the entire polygon.
    hasHoles = false;
    bound = S2LatLngRect.empty();
    for (int i = 0; i < numLoops(); ++i) {
      if (loop(i).sign() < 0) {
        hasHoles = true;
      } else {
        bound = bound.union(loop(i).getRectBound());
      }
    }
  }

  /**
   * Initializes a polygon from a set of {@link S2Loop}s.
   *
   * <p>Unlike {@link #init} this method assumes the caller already knows the
   * nesting of loops within other loops.  The passed-in map maps from parents
   * to their immediate child loops, with {@code null} mapping to the list of
   * top-most shell loops.  Immediate child loops must be completely spatially
   * contained within their parent loop, but not contained in any other loop,
   * except for ancestors of the parent.  This method avoids the cost of
   * determining nesting internally, but if the passed in nesting is wrong,
   * future operations on the S2Polygon may be arbitrarily incorrect.
   *
   * <p>Note that unlike {@link #init}, the passed-in container of loops is not
   * cleared; however, the passed-in loops become owned by the S2Polygon and
   * should not be modified by the caller after calling this method.
   *
   * @param nestedLoops loops with nesting.
   */
  public void initWithNestedLoops(Map<S2Loop, List<S2Loop>> nestedLoops) {
    Preconditions.checkState(numLoops() == 0);
    initLoop(null, -1, nestedLoops);

    // Empty the map as an indication we have taken ownership of the loops.
    nestedLoops.clear();

    // Compute the bounding rectangle of the entire polygon.
    hasHoles = false;
    bound = S2LatLngRect.empty();
    numVertices = 0;
    for (int i = 0; i < loops.size(); ++i) {
      numVertices += loops.get(i).numVertices();
      if (loops.get(i).sign() < 0) {
        hasHoles = true;
      } else {
        bound = bound.union(loops.get(i).getRectBound());
      }
    }
  }

  /**
   * Releases ownership of the loops of this polygon by appending them to the
   * given list. Resets the polygon to be empty.
   */
  public void release(List<S2Loop> loops) {
    loops.addAll(this.loops);
    this.loops.clear();
    bound = S2LatLngRect.empty();
    hasHoles = false;
    numVertices = 0;
  }

  /**
   * Returns true if each loop on this polygon is valid, and if the
   * relationships between all loops are valid.
   *
   * <p>Specifically, this verifies that {@link S2Loop#isValid} is true for each
   * {@link S2Loop}, and that {@link S2Polygon#isValid(List)} is true for the
   * whole list of loops.
   */
  public boolean isValid() {
    for (S2Loop loop : loops) {
      if (!loop.isValid()) {
        return false;
      }
    }
    return S2Polygon.isValid(loops);
  }

  /**
   * Returns true if the given loops form a valid polygon. Assumes that all
   * of the given loops have already been validated.
   */
  public static boolean isValid(List<S2Loop> loops) {
    // If a loop contains an edge AB, then no other loop may contain AB or BA.
    // We only need this test if there are at least two loops, assuming that
    // each loop has already been validated.
    if (loops.size() > 1) {
      Map<UndirectedEdge, LoopVertexIndexPair> edges = Maps.newHashMap();
      for (int i = 0; i < loops.size(); ++i) {
        S2Loop lp = loops.get(i);
        for (int j = 0; j < lp.numVertices(); ++j) {
          UndirectedEdge key = new UndirectedEdge(lp.vertex(j), lp.vertex(j + 1));
          LoopVertexIndexPair value = new LoopVertexIndexPair(i, j);
          if (edges.containsKey(key)) {
            LoopVertexIndexPair other = edges.get(key);
            log.info(
                "Duplicate edge: loop " + i + ", edge " + j + " and loop " + other.getLoopIndex()
                    + ", edge " + other.getVertexIndex());
            return false;
          } else {
            edges.put(key, value);
          }
        }
      }
    }

    // Verify that no loop covers more than half of the sphere, and that no
    // two loops cross.
    for (int i = 0; i < loops.size(); ++i) {
      if (!loops.get(i).isNormalized()) {
        log.info("Loop " + i + " encloses more than half the sphere");
        return false;
      }
      for (int j = i + 1; j < loops.size(); ++j) {
        // This test not only checks for edge crossings, it also detects
        // cases where the two boundaries cross at a shared vertex.
        if (loops.get(i).containsOrCrosses(loops.get(j)) < 0) {
          log.info("Loop " + i + " crosses loop " + j);
          return false;
        }
      }
    }
    return true;
  }

  public int numLoops() {
    return loops.size();
  }

  public S2Loop loop(int k) {
    return loops.get(k);
  }

  /**
   * Returns the index of the parent of loop {@code k}, or -1 if it has no
   * parent.
   */
  public int getParent(int k) {
    int depth = loop(k).depth();
    if (depth == 0) {
      return -1; // Optimization.
    }
    while (--k >= 0 && loop(k).depth() >= depth) {
      // spin
    }
    return k;
  }

  /**
   * Returns the index of the last loop that is contained within loop {@code k}.
   * Returns {@code numLoops() - 1} if {@code k < 0}. Note that loops
   * are indexed according to a preorder traversal of the nesting hierarchy, so
   * the immediate children of loop {@code k} can be found by iterating over
   * loops {@code (k+1)..getLastDescendant(k)} and selecting those whose depth
   * is equal to {@code (loop(k).depth() + 1)}.
   */
  public int getLastDescendant(int k) {
    if (k < 0) {
      return numLoops() - 1;
    }
    int depth = loop(k).depth();
    while (++k < numLoops() && loop(k).depth() > depth) {
      // spin
    }
    return k - 1;
  }

  private S2AreaCentroid getAreaCentroid(boolean doCentroid) {
    double areaSum = 0;
    S2Point centroidSum = S2Point.ORIGIN;
    for (int i = 0; i < numLoops(); ++i) {
      S2AreaCentroid areaCentroid = doCentroid ? loop(i).getAreaAndCentroid() : null;
      double loopArea = doCentroid ? areaCentroid.getArea() : loop(i).getArea();

      int loopSign = loop(i).sign();
      areaSum += loopSign * loopArea;
      if (doCentroid) {
        S2Point currentCentroid = areaCentroid.getCentroid();
        centroidSum =
            new S2Point(centroidSum.x + loopSign * currentCentroid.x,
                centroidSum.y + loopSign * currentCentroid.y,
                centroidSum.z + loopSign * currentCentroid.z);
      }
    }

    return new S2AreaCentroid(areaSum, doCentroid ? centroidSum : null);
  }

  /**
   * Returns the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops (the area is between 0 and 4*Pi) and the true
   * centroid of the polygon, weighted by the area of the polygon
   * (see s2.h for details on centroids). Note that the centroid might not be
   * contained by the polygon.
   */
  public S2AreaCentroid getAreaAndCentroid() {
    return getAreaCentroid(true);
  }

  /**
   * Returns the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops. The return value is between 0 and 4*Pi.
   */
  public double getArea() {
    return getAreaCentroid(false).getArea();
  }

  /**
   * Returns the true centroid of the polygon, weighted by the area of the
   * polygon (see s2.h for details on centroids). Note that the centroid might
   * not be contained by the polygon.
   */
  public S2Point getCentroid() {
    return getAreaCentroid(true).getCentroid();
  }

  /**
   * Returns the shortest distance from a point P to this polygon, given as the
   * angle formed between P, the origin, and the nearest point on the polygon to
   * P. This angle in radians is equivalent to the arclength along the unit
   * sphere.
   *
   * <p>If the point is contained inside the polygon, the distance returned is
   * 0.
   */
  public S1Angle getDistance(S2Point p) {
    if (contains(p)) {
      return S1Angle.radians(0);
    }

    // The furthest point from p on the sphere is its antipode, which is an
    // angle of PI radians. This is an upper bound on the angle.
    S1Angle minDistance = S1Angle.radians(Math.PI);
    for (int i = 0; i < numLoops(); i++) {
      minDistance = S1Angle.min(minDistance, loop(i).getDistance(p));
    }

    return minDistance;
  }

  /**
   * Returns a point on the polygon that is closest to point P. The distance
   * between these two points should be the result of
   * {@link #getDistance(S2Point)}.
   *
   * <p>If point P is contained within the loop, it is returned.
   *
   * <p>The polygon must not be empty.
   */
  public S2Point project(S2Point p) {
    Preconditions.checkState(!loops.isEmpty());
    if (contains(p)) {
      return p;
    }
    S2Point normalized = S2Point.normalize(p);

    // The furthest point from p on the sphere is its antipode, which is an
    // angle of PI radians. This is an upper bound on the angle.
    S1Angle minDistance = S1Angle.radians(Math.PI);
    int minLoopIndex = 0;
    int minVertexIndex = 0;
    for (int loopIndex = 0; loopIndex < loops.size(); loopIndex++) {
      S2Loop loop = loops.get(loopIndex);
      for (int vertexIndex = 0; vertexIndex < loop.numVertices(); vertexIndex++) {
        S1Angle distanceToSegment = S2EdgeUtil.getDistance(
            normalized, loop.vertex(vertexIndex), loop.vertex(vertexIndex + 1));
        if (minDistance.greaterThan(distanceToSegment)) {
          minDistance = distanceToSegment;
          minLoopIndex = loopIndex;
          minVertexIndex = vertexIndex;
        }
      }
    }
    S2Loop minLoop = loop(minLoopIndex);
    S2Point closestPoint = S2EdgeUtil.getClosestPoint(
        p, minLoop.vertex(minVertexIndex), minLoop.vertex(minVertexIndex + 1));
    return closestPoint;
  }


  /**
   * Returns true if this polygon contains the given other polygon, i.e., if
   * polygon A contains all points contained by polygon B.
   */
  public boolean contains(S2Polygon b) {
    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.contains does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).contains(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more
    // efficient S2Loop::Contains method (rather than ContainsOrCrosses),
    // but it's worthwhile to do our own bounds check first.
    if (!bound.contains(b.getRectBound())) {
      // If the union of the bounding boxes spans the full longitude range,
      // it is still possible that polygon A contains B. (This is only
      // possible if at least one polygon has multiple shells.)
      if (!bound.lng().union(b.getRectBound().lng()).isFull()) {
        return false;
      }
    }
    if (!hasHoles && !b.hasHoles) {
      for (int j = 0; j < b.numLoops(); ++j) {
        if (!anyLoopContains(b.loop(j))) {
          return false;
        }
      }
      return true;
    }

    // This could be implemented more efficiently for polygons with lots of
    // holes by keeping a copy of the LoopMap computed during initialization.
    // However, in practice most polygons are one loop, and multiloop polygons
    // tend to consist of many shells rather than holes. In any case, the real
    // way to get more efficiency is to implement a sub-quadratic algorithm
    // such as building a trapezoidal map.

    // Every shell of B must be contained by an odd number of loops of A,
    // and every hole of A must be contained by an even number of loops of B.
    return containsAllShells(b) && b.excludesAllHoles(this);
  }

  /**
   * Returns true if this polygon intersects the given other polygon, i.e., if
   * there is a point that is contained by both polygons.
   */
  public boolean intersects(S2Polygon b) {
    // A.intersects(B) if and only if !complement(A).contains(B). However,
    // implementing a complement() operation is trickier than it sounds,
    // and in any case it's more efficient to test for intersection directly.

    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.intersects does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).intersects(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more
    // efficient S2Loop.intersects method. The polygons intersect if and
    // only if some pair of loop regions intersect.
    if (!bound.intersects(b.getRectBound())) {
      return false;
    }
    if (!hasHoles && !b.hasHoles) {
      for (int i = 0; i < numLoops(); ++i) {
        for (int j = 0; j < b.numLoops(); ++j) {
          if (loop(i).intersects(b.loop(j))) {
            return true;
          }
        }
      }
      return false;
    }

    // Otherwise if any shell of B is contained by an odd number of loops of A,
    // or any shell of A is contained by an odd number of loops of B, there is
    // an intersection.
    return intersectsAnyShell(b) || b.intersectsAnyShell(this);
  }

  /**
   * Indexing structure to efficiently {@link #clipEdge} of a polygon. This is
   * an abstract class because we need to use if for both polygons (for
   * {@link #initToIntersection} and friends) and for sets of lists of points
   * (for initToSimplified() future?).
   *
   * <p>Usage: In your subclass, create an array of vertex counts for each loop
   * in the loop sequence and pass it to this constructor. Overwrite
   * {@link #edgeFromTo}, calling {@link #decodeIndex} and use the resulting two
   * indices to access your vertices.
   */
  private abstract static class S2LoopSequenceIndex extends S2EdgeIndex {
    /**
     * Map from the uni-dimensional edge index to the loop this edge belongs to.
     */
    private final int[] indexToLoop;

    /**
     * Reverse of {@link #indexToLoop}: maps a loop index to the uni-dimensional
     * index of the first edge in the loop.
     */
    private final int[] loopToFirstIndex;

    /**
     * Must be called by each subclass with the array of vertices per loop. The
     * length of the array is the number of loops, and the {@code i}
     * <sup>th</sup> loop's vertex count is in the {@code i} <sup>th</sup> index
     * of the array.
     */
    public S2LoopSequenceIndex(int[] numVertices) {
      int totalEdges = 0;
      for (int edges : numVertices) {
        totalEdges += edges;
      }
      indexToLoop = new int[totalEdges];
      loopToFirstIndex = new int[numVertices.length];

      totalEdges = 0;
      for (int j = 0; j < numVertices.length; j++) {
        loopToFirstIndex[j] = totalEdges;
        for (int i = 0; i < numVertices[j]; i++) {
          indexToLoop[totalEdges] = j;
          totalEdges++;
        }
      }
    }

    public final LoopVertexIndexPair decodeIndex(int index) {
      int loopIndex = indexToLoop[index];
      int vertexInLoop = index - loopToFirstIndex[loopIndex];
      return new LoopVertexIndexPair(loopIndex, vertexInLoop);
    }

    @Override
    public final int getNumEdges() {
      return indexToLoop.length;
    }

    /**
     * Mark the {@link #edgeFromTo} method abstract again, so children of this
     * class <b>must</b> implement it without using {@link #edgeFrom(int)}
     * and {@link #edgeTo(int)}.
     */
    @Override
    public abstract S2Edge edgeFromTo(int index);

    @Override
    public S2Point edgeFrom(int index) {
      return edgeFromTo(index).getStart();
    }

    @Override
    protected S2Point edgeTo(int index) {
      return edgeFromTo(index).getEnd();
    }
  }

  /**
   * Indexing structure for an {@link S2Polygon}.
   */
  private static final class S2PolygonIndex extends S2LoopSequenceIndex {
    private final S2Polygon poly;
    private final boolean reverse;

    /**
     * Returns number of vertices per loop.
     */
    private static int[] getVertices(S2Polygon poly) {
      int[] vertices = new int[poly.numLoops()];
      for (int i = 0; i < vertices.length; i++) {
        vertices[i] = poly.loop(i).numVertices();
      }
      return vertices;
    }

    public S2PolygonIndex(S2Polygon poly, boolean reverse) {
      super(getVertices(poly));
      this.poly = poly;
      this.reverse = reverse;
    }

    @Override
    public S2Edge edgeFromTo(int index) {
      LoopVertexIndexPair indices = decodeIndex(index);
      int loopIndex = indices.getLoopIndex();
      int vertexInLoop = indices.getVertexIndex();
      S2Loop loop = poly.loop(loopIndex);
      int fromIndex;
      int toIndex;
      if (loop.isHole() ^ reverse) {
        fromIndex = loop.numVertices() - 1 - vertexInLoop;
        toIndex = 2 * loop.numVertices() - 2 - vertexInLoop;
      } else {
        fromIndex = vertexInLoop;
        toIndex = vertexInLoop + 1;
      }
      S2Point from = loop.vertex(fromIndex);
      S2Point to = loop.vertex(toIndex);
      return new S2Edge(from, to);
    }
  }

  /**
   * Clips the boundary of A to the interior of B, and adds the resulting edges
   * to "builder". Shells are directed CCW and holes are directed clockwise,
   * unless {@code reverseA} or {@code reverseB} is true, in which case these
   * directions in the corresponding polygon are reversed. If {@code invertB}
   * is true, the boundary of A is clipped to the exterior, rather than to the
   * interior of B. If {@code adSharedEdges} is true, then the output will
   * include any edges that are shared between A and B (both edges must be in
   * the same direction after any edge reversals are taken into account).
   */
  private static void clipBoundary(
      final S2Polygon a,
      boolean reverseA,
      final S2Polygon b,
      boolean reverseB,
      boolean invertB,
      boolean addSharedEdges,
      S2PolygonBuilder builder) {
    S2PolygonIndex bIndex = new S2PolygonIndex(b, reverseB);
    bIndex.predictAdditionalCalls(a.getNumVertices());

    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    for (S2Loop aLoop : a.loops) {
      int n = aLoop.numVertices();
      int dir = (aLoop.isHole() ^ reverseA) ? -1 : 1;
      boolean inside = b.contains(aLoop.vertex(0)) ^ invertB;
      for (int j = (dir > 0) ? 0 : n; n > 0; --n, j += dir) {
        S2Point a0 = aLoop.vertex(j);
        S2Point a1 = aLoop.vertex(j + dir);
        intersections.clear();
        bIndex.clipEdge(a0, a1, addSharedEdges, intersections);

        if (inside) {
          intersections.add(new ParametrizedS2Point(0.0, a0));
        }
        inside = ((intersections.size() & 0x1) == 0x1);
        // assert ((b.contains(a1) ^ invertB) == inside);
        if (inside) {
          intersections.add(new ParametrizedS2Point(1.0, a1));
        }

        // Remove duplicates and produce a list of unique intersections.
        Collections.sort(intersections);
        for (int size = intersections.size(), i = 1; i < size; i += 2) {
          builder.addEdge(intersections.get(i - 1).getPoint(), intersections.get(i).getPoint());
        }
      }
    }
  }

  /**
   * Returns the total number of vertices in all loops.
   */
  public int getNumVertices() {
    return this.numVertices;
  }

  /**
   * Initializes this polygon to the intersection, union, or difference (A - B)
   * of the given two polygons. The {@code vertexMergeRadius} determines how
   * close two vertices must be to be merged together and how close a vertex
   * must be to an edge in order to be spliced into it (see
   * {@link S2PolygonBuilder} for details). By default, the merge radius is just
   * large enough to compensate for errors that occur when computing
   * intersection points between edges
   * ({@link S2EdgeUtil#DEFAULT_INTERSECTION_TOLERANCE}).
   *
   * <p>If you are going to convert the resulting polygon to a lower-precision
   * format, it is necessary to increase the merge radius in order to get a
   * valid result after rounding (i.e., no duplicate vertices, etc). For
   * example, if you are going to convert them to {@code geostore.PolygonProto}
   * format, then {@code S1Angle.e7(1)} is a good value for
   * {@code vertexMergeRadius}.
   */
  public void initToIntersection(final S2Polygon a, final S2Polygon b) {
    initToIntersectionSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  public void initToIntersectionSloppy(
      final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    Preconditions.checkState(numLoops() == 0);
    if (!a.bound.intersects(b.bound)) {
      return;
    }

    // We want the boundary of A clipped to the interior of B,
    // plus the boundary of B clipped to the interior of A,
    // plus one copy of any directed edges that are in both boundaries.

    S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
        .setMergeDistance(vertexMergeRadius)
        .build();
    S2PolygonBuilder builder = new S2PolygonBuilder(options);
    clipBoundary(a, false, b, false, false, true, builder);
    clipBoundary(b, false, a, false, false, false, builder);
    if (!builder.assemblePolygon(this, null)) {
      // TODO(andriy): do something more meaningful here.
      log.severe("Bad directed edges");
    }
  }

  public void initToUnion(final S2Polygon a, final S2Polygon b) {
    initToUnionSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  public void initToUnionSloppy(final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    Preconditions.checkState(numLoops() == 0);

    // We want the boundary of A clipped to the exterior of B,
    // plus the boundary of B clipped to the exterior of A,
    // plus one copy of any directed edges that are in both boundaries.

    S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
        .setMergeDistance(vertexMergeRadius)
        .build();
    S2PolygonBuilder builder = new S2PolygonBuilder(options);
    clipBoundary(a, false, b, false, true, true, builder);
    clipBoundary(b, false, a, false, true, false, builder);
    if (!builder.assemblePolygon(this, null)) {
      // TODO(andriy): do something more meaningful here.
      log.severe("Bad directed edges");
    }
  }

  public void initToDifference(final S2Polygon a, final S2Polygon b) {
    initToDifferenceSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  public void initToDifferenceSloppy(
      final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    Preconditions.checkState(numLoops() == 0);

    // We want the boundary of A clipped to the exterior of B,
    // plus the reversed boundary of B clipped to the interior of A,
    // plus one copy of any edge in A that is also a reverse edge in B.

    S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
        .setMergeDistance(vertexMergeRadius)
        .build();
    S2PolygonBuilder builder = new S2PolygonBuilder(options);
    clipBoundary(a, false, b, true, true, true, builder);
    clipBoundary(b, true, a, false, false, false, builder);
    if (!builder.assemblePolygon(this, null)) {
      log.severe("Bad directed edges");
    }
  }

  /**
   * Returns a polygon that is the union of the given polygons. Note: clears the
   * List!
   */
  public static S2Polygon destructiveUnion(List<S2Polygon> polygons) {
    return destructiveUnionSloppy(polygons, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Returns a polygon that is the union of the given polygons; combines
   * vertices that form edges that are almost identical, as defined by
   * {@code vertexMergeRadius}. Note: clears the List!
   */
  public static S2Polygon destructiveUnionSloppy(
      List<S2Polygon> polygons, S1Angle vertexMergeRadius) {
    // Effectively create a priority queue of polygons in order of number of
    // vertices. Repeatedly union the two smallest polygons and add the result
    // to the queue until we have a single polygon to return.

    // map: # of vertices -> polygon
    TreeMultimap<Integer, S2Polygon> queue = TreeMultimap.create();

    for (S2Polygon polygon : polygons) {
      queue.put(polygon.getNumVertices(), polygon);
    }
    polygons.clear();

    Set<Map.Entry<Integer, S2Polygon>> queueSet = queue.entries();
    while (queueSet.size() > 1) {
      // Pop two simplest polygons from queue.
      queueSet = queue.entries();
      Iterator<Map.Entry<Integer, S2Polygon>> smallestIter = queueSet.iterator();

      Map.Entry<Integer, S2Polygon> smallest = smallestIter.next();
      int aSize = smallest.getKey().intValue();
      S2Polygon aPolygon = smallest.getValue();
      smallestIter.remove();

      smallest = smallestIter.next();
      int bSize = smallest.getKey().intValue();
      S2Polygon bPolygon = smallest.getValue();
      smallestIter.remove();

      // Union and add result back to queue.
      S2Polygon unionPolygon = new S2Polygon();
      unionPolygon.initToUnionSloppy(aPolygon, bPolygon, vertexMergeRadius);
      int unionSize = aSize + bSize;
      queue.put(unionSize, unionPolygon);
      // We assume that the number of vertices in the union polygon is the
      // sum of the number of vertices in the original polygons, which is not
      // always true, but will almost always be a decent approximation, and
      // faster than recomputing.
    }

    if (queue.isEmpty()) {
      return new S2Polygon();
    } else {
      return queue.get(queue.asMap().firstKey()).first();
    }
  }

  /**
   * Intersects this polygon with the {@link S2Polyline} {@code in} and returns
   * the resulting zero or more polylines. The polylines are ordered in the
   * order they would be encountered by traversing {@code in} from beginning to
   * end. Note that the output may include polylines with only one vertex, but
   * there will not be any zero-vertex polylines.
   *
   * <p>This is equivalent to calling {@link #intersectWithPolylineSloppy} with
   * the {@code vertexMergeRadius} set to
   * {@link S2EdgeUtil#DEFAULT_INTERSECTION_TOLERANCE}.
   */
  public List<S2Polyline> intersectWithPolyline(S2Polyline in) {
    return intersectWithPolylineSloppy(in, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Similar to {@link #intersectWithPolyline}, except that vertices will be
   * dropped as necessary to ensure that all adjacent vertices in the
   * sequence obtained by concatenating the output polylines will be
   * farther than {@code vertexMergeRadius} apart.  Note that this can change
   * the number of output polylines and/or yield single-vertex polylines.
   */
  public List<S2Polyline> intersectWithPolylineSloppy(
      S2Polyline in, S1Angle vertexMergeRadius) {
    return internalClipPolyline(false, in, vertexMergeRadius);
  }

  /**
   * Same as {@link #intersectWithPolyline}, but subtracts this polygon from
   * the given polyline.
   */
  public List<S2Polyline> subtractFromPolyline(S2Polyline in) {
    return subtractFromPolylineSloppy(in, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Same as {@link #intersectWithPolylineSloppy}, but subtracts this polygon
   * from the given polyline.
   */
  public List<S2Polyline> subtractFromPolylineSloppy(
      S2Polyline in, S1Angle vertexMergeRadius) {
    return internalClipPolyline(true, in, vertexMergeRadius);
  }

  /**
  * Clips the {@link S2Polyline} {@code a} to the interior of this polygon.
  * The resulting polyline(s) will be returned. If {@code invert} is
  * {@code true}, we clip {@code a} to the exterior of this polygon instead.
  * Vertices will be dropped such that adjacent vertices will not be closer
  * than {@code mergeRadius}.
  *
  * <p>We do the intersection/subtraction by walking the polyline edges.
  * For each edge, we compute all intersections with the polygon boundary
  * and sort them in increasing order of distance along that edge.
  * We then divide the intersection points into pairs, and output a
  * clipped polyline segment for each one. We keep track of whether we're
  * inside or outside of the polygon at all times to decide which segments
  * to output.
  */
  private List<S2Polyline> internalClipPolyline(
      boolean invert, S2Polyline a, S1Angle mergeRadius) {

    List<S2Polyline> out = Lists.newArrayList();
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    List<S2Point> vertices = Lists.newArrayList();
    S2PolygonIndex polyIndex = new S2PolygonIndex(this, false);
    int n = a.numVertices();
    boolean inside = contains(a.vertex(0)) ^ invert;

    for (int j = 0; j < n - 1; j++) {
      S2Point a0 = a.vertex(j);
      S2Point a1 = a.vertex(j + 1);
      polyIndex.clipEdge(a0, a1, true, intersections);
      if (inside) {
        intersections.add(new ParametrizedS2Point(0, a0));
      }
      inside = (intersections.size() & 1) != 0;
      // assert ((contains(a1) ^ invert) == inside);
      if (inside) {
        intersections.add(new ParametrizedS2Point(1, a1));
      }
      Collections.sort(intersections);
      // At this point we have a sorted array of vertex pairs representing
      // the edge(s) obtained after clipping (a0,a1) against the polygon.
      for (int k = 0; k < intersections.size(); k += 2) {
        if (intersections.get(k).equals(intersections.get(k + 1))) {
          continue;
        }
        S2Point v0 = intersections.get(k).getPoint();
        S2Point v1 = intersections.get(k + 1).getPoint();
        // If the gap from the previous vertex to this one is large
        // enough, start a new polyline.
        if (!vertices.isEmpty()) {
          S2Point back = vertices.get(vertices.size() - 1);
          if (back.angle(v0) > mergeRadius.radians()) {
            out.add(new S2Polyline(vertices));
            vertices.clear();
          }
        }
        // Append this segment to the current polyline, ignoring any
        // vertices that are too close to the previous vertex.
        if (vertices.isEmpty()) {
          vertices.add(v0);
        }
        S2Point back2 = vertices.get(vertices.size() - 1);
        if (back2.angle(v1) > mergeRadius.radians()) {
          vertices.add(v1);
        }
      }
      intersections.clear();
    }

    if (!vertices.isEmpty()) {
      out.add(new S2Polyline(vertices));
    }

    return out;
  }

  public boolean isNormalized() {
    Multiset<S2Point> vertices = HashMultiset.<S2Point>create();
    S2Loop lastParent = null;
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop child = loop(i);
      if (child.depth() == 0) {
        continue;
      }
      S2Loop parent = loop(getParent(i));
      if (parent != lastParent) {
        vertices.clear();
        for (int j = 0; j < parent.numVertices(); ++j) {
          vertices.add(parent.vertex(j));
        }
        lastParent = parent;
      }
      int count = 0;
      for (int j = 0; j < child.numVertices(); ++j) {
        if (vertices.count(child.vertex(j)) > 0) {
          ++count;
        }
      }
      if (count > 1) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if two polygons have the same boundary, except for vertex
   * perturbations. Both polygons must have loops with the same cyclic vertex
   * order and the same nesting hierarchy, but the vertex locations are allowed
   * to differ by up to {@code maxError}. Note: This method mostly useful only
   * for testing purposes.
   */
  boolean boundaryApproxEquals(S2Polygon b, double maxError) {
    if (numLoops() != b.numLoops()) {
      log.severe(
          "!= loops: " + Integer.toString(numLoops()) + " vs. " + Integer.toString(b.numLoops()));
      return false;
    }

    // For now, we assume that there is at most one candidate match for each
    // loop. (So far this method is just used for testing.)
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop aLoop = loop(i);
      boolean success = false;
      for (int j = 0; j < numLoops(); ++j) {
        S2Loop bLoop = b.loop(j);
        if (bLoop.depth() == aLoop.depth() && bLoop.boundaryApproxEquals(aLoop, maxError)) {
          success = true;
          break;
        }
      }
      if (!success) {
        return false;
      }
    }
    return true;
  }

  // S2Region interface (see S2Region.java for details):

  /** Returns a bounding spherical cap. */
  @Override
  public S2Cap getCapBound() {
    return bound.getCapBound();
  }


  /** Returns a bounding latitude-longitude rectangle. */
  @Override
  public S2LatLngRect getRectBound() {
    return bound;
  }

  /**
   * If this method returns true, the region completely contains the given cell.
   * Otherwise, either the region does not contain the cell or the containment
   * relationship could not be determined.
   */
  @Override
  public boolean contains(S2Cell cell) {
    if (numLoops() == 1) {
      return loop(0).contains(cell);
    }
    S2LatLngRect cellBound = cell.getRectBound();
    if (!bound.contains(cellBound)) {
      return false;
    }

    S2Loop cellLoop = new S2Loop(cell, cellBound);
    S2Polygon cellPoly = new S2Polygon(cellLoop);
    return contains(cellPoly);
  }

  /**
   * If this method returns false, the region does not intersect the given cell.
   * Otherwise, either region intersects the cell, or the intersection
   * relationship could not be determined.
   */
  @Override
  public boolean mayIntersect(S2Cell cell) {
    if (numLoops() == 1) {
      return loop(0).mayIntersect(cell);
    }
    S2LatLngRect cellBound = cell.getRectBound();
    if (!bound.intersects(cellBound)) {
      return false;
    }

    S2Loop cellLoop = new S2Loop(cell, cellBound);
    S2Polygon cellPoly = new S2Polygon(cellLoop);
    return intersects(cellPoly);
  }

  /**
   * The point {@code p} does not need to be normalized.
   */
  public boolean contains(S2Point p) {
    if (numLoops() == 1) {
      return loop(0).contains(p); // Optimization.
    }
    if (!bound.contains(p)) {
      return false;
    }
    boolean inside = false;
    for (int i = 0; i < numLoops(); ++i) {
      inside ^= loop(i).contains(p);
      if (inside && !hasHoles) {
        break; // Shells are disjoint.
      }
    }
    return inside;
  }

  /**
   * For each map entry, sorts the value list.
   */
  private static void sortValueLoops(Map<S2Loop, List<S2Loop>> loopMap) {
    for (S2Loop key : loopMap.keySet()) {
      Collections.sort(loopMap.get(key));
    }
  }

  private static void insertLoop(S2Loop newLoop, S2Loop parent, Map<S2Loop, List<S2Loop>> loopMap) {
    List<S2Loop> children = loopMap.get(parent);

    if (children == null) {
      children = Lists.newArrayList();
      loopMap.put(parent, children);
    }

    for (S2Loop child : children) {
      if (child.containsNested(newLoop)) {
        insertLoop(newLoop, child, loopMap);
        return;
      }
    }

    // No loop may contain the complement of another loop. (Handling this case
    // is significantly more complicated.)
    // assert (parent == null || !newLoop.containsNested(parent));

    // Some of the children of the parent loop may now be children of
    // the new loop.
    List<S2Loop> newChildren = loopMap.get(newLoop);
    for (int i = 0; i < children.size();) {
      S2Loop child = children.get(i);
      if (newLoop.containsNested(child)) {
        if (newChildren == null) {
          newChildren = Lists.newArrayList();
          loopMap.put(newLoop, newChildren);
        }
        newChildren.add(child);
        children.remove(i);
      } else {
        ++i;
      }
    }
    children.add(newLoop);
  }

  private void initLoop(S2Loop loop, int depth, Map<S2Loop, List<S2Loop>> loopMap) {
    if (loop != null) {
      loop.setDepth(depth);
      loops.add(loop);
    }
    List<S2Loop> children = loopMap.get(loop);
    if (children != null) {
      for (S2Loop child : children) {
        initLoop(child, depth + 1, loopMap);
      }
    }
  }

  private int containsOrCrosses(S2Loop b) {
    boolean inside = false;
    for (int i = 0; i < numLoops(); ++i) {
      int result = loop(i).containsOrCrosses(b);
      if (result < 0) {
        return -1; // The loop boundaries intersect.
      }
      if (result > 0) {
        inside ^= true;
      }
    }
    return inside ? 1 : 0; // True if loop B is contained by the polygon.
  }

  /** Return true if any loop contains the given loop. */
  private boolean anyLoopContains(S2Loop b) {
    for (int i = 0; i < numLoops(); ++i) {
      if (loop(i).contains(b)) {
        return true;
      }
    }
    return false;
  }

  /** Return true if this polygon (A) contains all the shells of B. */
  private boolean containsAllShells(S2Polygon b) {
    for (int j = 0; j < b.numLoops(); ++j) {
      if (b.loop(j).sign() < 0) {
        continue;
      }
      if (containsOrCrosses(b.loop(j)) <= 0) {
        // Shell of B is not contained by A, or the boundaries intersect.
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if this polygon (A) excludes (i.e., does not intersect) all
   * holes of B.
   */
  private boolean excludesAllHoles(S2Polygon b) {
    for (int j = 0; j < b.numLoops(); ++j) {
      if (b.loop(j).sign() > 0) {
        continue;
      }
      if (containsOrCrosses(b.loop(j)) != 0) {
        // Hole of B is contained by A, or the boundaries intersect.
        return false;
      }
    }
    return true;
  }

  /** Return true if this polygon (A) intersects any shell of B. */
  private boolean intersectsAnyShell(S2Polygon b) {
    for (int j = 0; j < b.numLoops(); ++j) {
      if (b.loop(j).sign() < 0) {
        continue;
      }
      if (containsOrCrosses(b.loop(j)) != 0) {
        // Shell of B is contained by A, or the boundaries intersect.
        return true;
      }
    }
    return false;
  }

  /**
   * Returns a human readable representation of the polygon.
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("Polygon: (").append(numLoops()).append(") loops:\n");
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop s2Loop = loop(i);
      sb.append("loop <\n");
      for (int v = 0; v < s2Loop.numVertices(); ++v) {
        S2Point s2Point = s2Loop.vertex(v);
        sb.append(s2Point.toDegreesString());
        sb.append("\n"); // end of vertex
      }
      sb.append(">\n"); // end of loop
    }
    return sb.toString();
  }

  private static final class UndirectedEdge {
    // Note: An UndirectedEdge and an S2Edge can never be considered equal (in
    // terms of the equals() method) and hence they re not be related types.
    // If you need to convert between the types then separate conversion
    // methods should be introduced.

    private final S2Point a;
    private final S2Point b;

    public UndirectedEdge(S2Point start, S2Point end) {
      this.a = start;
      this.b = end;
    }

    public S2Point getStart() {
      return a;
    }

    public S2Point getEnd() {
      return b;
    }

    @Override
    public String toString() {
      return String.format("Edge: (%s <-> %s)\n   or [%s <-> %s]",
          a.toDegreesString(), b.toDegreesString(), a, b);
    }

    @Override
    public boolean equals(Object o) {
      if (o == null || !(o instanceof UndirectedEdge)) {
        return false;
      }
      UndirectedEdge other = (UndirectedEdge) o;
      return ((getStart().equals(other.getStart()) && getEnd().equals(other.getEnd()))
          || (getStart().equals(other.getEnd()) && getEnd().equals(other.getStart())));
    }

    @Override
    public int hashCode() {
      return getStart().hashCode() + getEnd().hashCode();
    }
  }

  private static final class LoopVertexIndexPair {
    private final int loopIndex;
    private final int vertexIndex;

    public LoopVertexIndexPair(int loopIndex, int vertexIndex) {
      this.loopIndex = loopIndex;
      this.vertexIndex = vertexIndex;
    }

    public int getLoopIndex() {
      return loopIndex;
    }

    public int getVertexIndex() {
      return vertexIndex;
    }
  }
}
