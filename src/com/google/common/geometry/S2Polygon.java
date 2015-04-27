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

import static com.google.common.geometry.S2Projections.PROJ;

import com.google.common.annotations.GwtCompatible;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultimap;
import com.google.common.geometry.S2Projections.FaceSiTi;
import com.google.common.geometry.S2Shape.MutableEdge;

import java.io.Serializable;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Logger;

/**
 * <p>An S2Polygon is an S2Region object that represents a polygon. A polygon is defined by zero or
 * more loops; recall that the interior of a loop is defined to be its left-hand side (see
 * {@link S2Loop}.)
 *
 * <p>There are two different conventions for creating an S2Polygon:
 *
 * <p>First, {@link #initNested(List)} expects the input loops to be nested hierarchically. The
 * polygon interior then consists of the set of points contained by an odd number of loops. So for
 * example, a circular region with a hole in it would be defined as two CCW loops, with one loop
 * containing the other. The loops can be provided in any order.
 *
 * <p>When the orientation of the input loops is unknown, the nesting requirement is typically met
 * by calling {@link S2Loop#normalize()} on each loop (which inverts the loop if necessary so that
 * it encloses at most half the sphere). But in fact any set of loops can be used as long as (1)
 * there is no pair of loops that cross, and (2) there is no pair of loops whose union is the entire
 * sphere.
 *
 * <p>Second, {@link #initOriented(List)} expects the input loops to be oriented such that the
 * polygon interior is on the left-hand side of every loop. So for example, a circular region with a
 * hole in it would be defined using a CCW outer loop and a CW inner loop. The loop orientations
 * must all be consistent; for example, it is not valid to have one CCW loop nested inside another
 * CCW loop, because the region between the two loops is on the left-hand side of one loop and the
 * right-hand side of the other.
 *
 * <p>Most clients will not call these methods directly; instead they should use
 * {@link S2PolygonBuilder}, which has better support for dealing with imperfect data.
 *
 * <p>When the polygon is initialized, the given loops are automatically converted into a canonical
 * form consisting of "shells" and "holes". Shells and holes are both oriented CCW, and are nested
 * hierarchically. The loops are reordered to correspond to a preorder traversal of the nesting
 * hierarchy; initOriented may also invert some loops.
 *
 * <p>Polygons may represent any region of the sphere with a polygonal boundary, including the
 * entire sphere (known as the "full" polygon). The full polygon consists of a single full loop (see
 * {@link S2Loop#full()}), whereas the empty polygon has no loops at all.
 *
 * <p>Polygons have the following restrictions:
 *
 * <ul>
 * <li>Loops may not cross, i.e. the boundary of a loop may not intersect both the interior and
 * exterior of any other loop.
 * <li>Loops may not share edges, i.e. if a loop contains an edge AB, then no other loop may contain
 * AB or BA.
 * <li>Loops may share vertices, however no vertex may appear twice in a single loop (see S2Loop).
 * <li>No loop may be empty. The full loop may appear only in the full polygon.
 * </ul>
 *
 */
@GwtCompatible(serializable = true)
public final strictfp class S2Polygon implements S2Region, Comparable<S2Polygon>, Serializable {
  private static final Logger log = Platform.getLoggerForClass(S2Polygon.class);

  /**
   * The loops of this polygon. There is no total ordering of the loops, but a nested loop always
   * follows its containing loop, and all loops between parent and child are nested somewhere under
   * the parent.
   */
  private final List<S2Loop> loops = Lists.newArrayList();

  /**
   * {@code bound} is a conservative bound on all points contained by this polygon:
   * If A.contains(P), then A.bound.contains(new S2LatLng(P)).
   */
  private S2LatLngRect bound;

  /**
   * Since "bound" is not exact, it is possible that a polygon A contains another polygon B whose
   * bounds are slightly larger. "subregionBound" has been expanded sufficiently to account for this
   * error, i.e. if A.Contains(B), then A.subregionBound.contains(B.bound).
   */
  private S2LatLngRect subregionBound;

  /** The spatial index for this S2Polygon. */
  @VisibleForTesting
  transient S2ShapeIndex index;

  /**
   * In general we build the index the first time it is needed, but we make an exception for
   * contains(S2Point) because this method has a simple brute force implementation that is
   * relatively cheap. For this one method we keep track of the number of calls made and only build
   * the index once enough calls have been made that we think an index would be worthwhile.
   */
  private AtomicInteger unindexedContainsCalls = new AtomicInteger();

  /** True if this polygon has at least one hole. */
  private boolean hasHoles = false;

  /** Total number of vertices in all loops. */
  private int numVertices = 0;

  /** Creates an empty polygon. It can be made non-empty by calling {@link #init(List)}. */
  public S2Polygon() {
    bound = S2LatLngRect.empty();
    subregionBound = S2LatLngRect.empty();
    initIndex();
  }

  /** Creates an S2Polygon for a given cell. */
  public S2Polygon(S2Cell cell) {
    loops.add(new S2Loop(cell));
    initOneLoop();
  }

  /**
   * Creates an empty polygon and then calls {@link #initNested(List)} with the given loops.
   * Clears the given list.
   */
  public S2Polygon(List<S2Loop> loops) {
    initNested(loops);
  }

  /** Copy constructor. */
  public S2Polygon(S2Loop loop) {
    this.numVertices = loop.numVertices();
    this.bound = loop.getRectBound();
    this.subregionBound = loop.getSubregionBound();
    loops.add(loop);
    initIndex();
  }

  /** Copy constructor. */
  public S2Polygon(S2Polygon src) {
    copy(src);
  }

  /** Initializes this polygon to a copy of the given polygon. */
  void copy(S2Polygon src) {
    this.bound = src.bound;
    this.subregionBound = src.subregionBound;
    this.hasHoles = src.hasHoles;
    this.numVertices = src.numVertices;
    for (int i = 0; i < src.numLoops(); ++i) {
      loops.add(new S2Loop(src.loop(i)));
    }
    initIndex();
  }

  private void initIndex() {
    // See S2Loop for the details behind 'maxUnindexedContainsCalls'.
    int maxUnindexedContainsCalls;
    if (numVertices <= 8) {
      maxUnindexedContainsCalls = 10;
    } else if (numVertices <= 8192) {
      maxUnindexedContainsCalls = 50;
    } else if (numVertices <= 50000) {
      maxUnindexedContainsCalls = 10;
    } else {
      maxUnindexedContainsCalls = 2;
    }
    this.unindexedContainsCalls.set(maxUnindexedContainsCalls);

    // assert (index == null || index.numShapeIds() == 0);
    index = new S2ShapeIndex();
    for (int i = 0; i < numLoops(); ++i) {
      index.add(loop(i));
    }
    // assert (isValid());
  }

  /** Returns the same instance after initializing transient fields. */
  private Object readResolve() {
    initIndex();
    return this;
  }

  @Override
  public boolean equals(Object o) {
    if (o instanceof S2Polygon) {
      S2Polygon that = (S2Polygon) o;
      return this.numVertices == that.numVertices
          && this.bound.equals(that.bound)
          && this.loops.equals(that.loops);
    } else {
      return false;
    }
  }

  @Override
  public int hashCode() {
    return bound.hashCode();
  }

  /**
   * Comparator (needed by Comparable interface). For two polygons to be compared as equal:
   * <ul>
   * <li>They must have the same number of loops
   * <li>The loops must be ordered in the same way (this is guaranteed by the total ordering imposed
   *     by {@link #sortValueLoops})
   * <li>Loops must be logically equivalent (even if ordered with a different starting point, e.g.
   *     ABCD and BCDA).
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

  /** Initializes a polygon by calling {@link #initNested(List)}. */
  public void init(List<S2Loop> loops) {
    initNested(loops);
  }

  /**
   * <p>Initializes this polygon from a set of hierarchically nested loops. The polygon interior
   * consists of the points contained by an odd number of loops. (Recall that a loop contains the
   * set of points on its left-hand side.)
   *
   * <p>This method takes ownership of the given loops and clears the given list. It then figures
   * out the loop nesting hierarchy and assigns every loop a depth. Shells have even depths, and
   * holes have odd depths. Note that the loops are reordered so the hierarchy can be traversed more
   * easily (see {@link #getParent(int)}, {@link #getLastDescendant(int)}, and
   * {@link S2Loop#depth()}).
   *
   * <p>This method may be called more than once, in which case any existing loops are deleted
   * before being replaced by the input loops.
   */
  public void initNested(List<S2Loop> loops) {
    // assert isValid(loops);
    clearLoops();

    if (loops.size() == 1) {
      this.loops.clear();
      // Since we know size()==1, use remove(0) instead of get(0) followed by clear().
      this.loops.add(loops.remove(0));
      initOneLoop();
      return;
    }

    Map<S2Loop, List<S2Loop>> loopMap = Maps.newIdentityHashMap();
    // Yes, a null key is valid. It is used here to refer to the root of the loopMap.
    loopMap.put(null, Lists.<S2Loop>newArrayList());

    for (S2Loop loop : loops) {
      insertLoop(loop, null, loopMap);
    }
    loops.clear();

    // Sort all of the lists of loops; in this way we guarantee a total ordering on loops in the
    // polygon. Loops will be sorted by their natural ordering, while also preserving the
    // requirement that each loop is immediately followed by its descendants in the nesting
    // hierarchy.
    //
    // TODO(andriy): as per kirilll in CL 18750833 code review comments: This should work for now,
    // but I think it's possible to guarantee the correct order inside insertLoop by searching for
    // the correct position in the children list before inserting.
    sortValueLoops(loopMap);

    // Reorder the loops in depth-first traversal order. Starting at null == starting at the root.
    initLoop(null, -1, loopMap);

    // TODO(dbeaumont): Add tests or preconditions for these asserts (here and elsewhere).
    // forall i != j : containsChild(loop(i), loop(j), loopMap) == loop(i).containsNested(loop(j)));

    initLoopProperties();
  }

  /**
   * Like {@link #initNested(List)}, but expects loops to be oriented such that the polygon interior
   * is on the left-hand side of all loops. This implies that shells and holes should have opposite
   * orientations in the input to this method. (During initialization, loops representing holes will
   * automatically be inverted.)
   */
  public void initOriented(List<S2Loop> loops) {
    // Here is the algorithm:
    //
    // 1. Remember which of the given loops contain S2.origin().
    //
    // 2. Invert loops as necessary to ensure that they are nestable (i.e., no loop contains the
    //    complement of any other loop).  This may result in a set of loops corresponding to the
    //    complement of the given polygon, but we will fix that problem later.
    //
    //    We make the loops nestable by first normalizing all the loops (i.e., inverting any loops
    //    whose turning angle is negative).  This handles all loops except those whose turning angle
    //    is very close to zero (within the maximum error tolerance).  Any such loops are inverted
    //    if and only if they contain S2.origin().  (In theory this step is only necessary if there
    //    are at least two such loops.)  The resulting set of loops is guaranteed to be nestable.
    //
    // 3. Build the polygon.  This yields either the desired polygon or its complement.
    //
    // 4. If there is at least one loop, we find a loop L that is adjacent to S2.origin() (where
    //    "adjacent" means that there exists a path connecting S2.origin() to some vertex of L such
    //    that the path does not cross any loop).  There may be a single such adjacent loop, or
    //    there may be several (in which case they should all have the same contains_origin()
    //    value).  We choose L to be the loop containing the origin whose depth is greatest, or
    //    loop(0) (a top-level shell) if no such loop exists.
    //
    // 5. If (L originally contained origin) != (polygon contains origin), we invert the polygon.
    //    This is done by inverting a top-level shell whose turning angle is minimal and then fixing
    //    the nesting hierarchy.  Note that because we normalized all the loops initially, this step
    //    is only necessary if the polygon requires at least one non-normalized loop to represent
    //    it.
    Preconditions.checkState(this.loops.isEmpty());

    Set<S2Loop> containedOrigin = Sets.newIdentityHashSet();
    for (S2Loop loop : loops) {
      if (loop.containsOrigin()) {
        containedOrigin.add(loop);
      }
      double angle = loop.getTurningAngle();
      if (Math.abs(angle) > 1e-14) {
        // TODO(user): GetTurningAngleMaxError()
        // Normalize the loop.
        if (angle < 0) {
          loop.invert();
        }
      } else {
        // Ensure that the loop does not contain the origin.
        if (loop.containsOrigin()) {
          loop.invert();
        }
      }
    }
    initNested(loops);
    if (numLoops() > 0) {
      S2Loop originLoop = loop(0);
      boolean polygonContainsOrigin = false;
      for (int i = 0; i < numLoops(); ++i) {
        if (loop(i).containsOrigin()) {
          polygonContainsOrigin ^= true;
          originLoop = loop(i);
        }
      }
      if (containedOrigin.contains(originLoop) != polygonContainsOrigin) {
        invert();
      }
    }

    // Verify that the original loops had consistent shell/hole orientations. Each original loop L
    // should have been inverted if and only if it now represents a hole.
    for (S2Loop loop : loops) {
      assert (containedOrigin.contains(loop) != loop.containsOrigin()) == loop.isHole();
    }
  }

  /** Computes hasHoles, numVertices, bound, subregionBound, and the index.. */
  private void initLoopProperties() {
    hasHoles = false;
    numVertices = 0;
    S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();
    for (S2Loop loop : loops) {
      if (loop.isHole()) {
        hasHoles = true;
      } else {
        builder.union(loop.getRectBound());
      }
      numVertices += loop.numVertices();
    }
    bound = builder.build();
    subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(bound);
    initIndex();
  }

  /** Given that loops contains a single loop, initializes all other fields. */
  private void initOneLoop() {
    assert 1 == loops.size();
    S2Loop loop = loops.get(0);
    loop.setDepth(0);
    hasHoles = false;
    numVertices = loop.numVertices();
    bound = loop.getRectBound();
    subregionBound = loop.getSubregionBound();
    initIndex();
  }

  /**
   * Initializes a polygon from a set of {@link S2Loop}s.
   *
   * <p>Unlike {@link #init} this method assumes the caller already knows the nesting of loops
   * within other loops.  The passed-in map maps from parents to their immediate child loops, with
   * {@code null} mapping to the list of top-most shell loops.  Immediate child loops must be
   * completely spatially contained within their parent loop, but not contained in any other loop,
   * except for ancestors of the parent.  This method avoids the cost of determining nesting
   * internally, but if the passed in nesting is wrong, future operations on the S2Polygon may be
   * arbitrarily incorrect.
   *
   * <p>Note that unlike {@link #init}, the passed-in container of loops is not cleared; however,
   * the passed-in loops become owned by the S2Polygon and should not be modified by the caller
   * after calling this method.
   *
   * @param nestedLoops loops with nesting.
   */
  public void initWithNestedLoops(Map<S2Loop, List<S2Loop>> nestedLoops) {
    Preconditions.checkState(numLoops() == 0);
    initLoop(null, -1, nestedLoops);

    // Empty the map as an indication we have taken ownership of the loops.
    nestedLoops.clear();

    initLoopProperties();
  }

  /**
   * Appends the loops of this polygon to the given list and resets this polygon to be empty.
   */
  public void release(List<S2Loop> loops) {
    loops.addAll(this.loops);
    // Reset the polygon to be empty.
    this.loops.clear();
    bound = S2LatLngRect.empty();
    subregionBound = S2LatLngRect.empty();
    hasHoles = false;
    numVertices = 0;
    initIndex();
  }

  private void clearLoops() {
    loops.clear();
    initIndex();
  }

  /**
   * Returns true if the given loops form a valid polygon, including checking whether the loops
   * themselves are valid.
   */
  public static boolean isValid(List<S2Loop> loops) {
    return new S2Polygon(Lists.newArrayList(loops)).isValid();
  }

  /**
   * Returns true if each loop on this polygon is valid, and if the relationships between all loops
   * are valid.
   *
   * <p>Specifically, this verifies that {@link S2Loop#isValid} is true for each {@link S2Loop}, and
   * that {@link S2Polygon#isValid(List)} is true for the whole list of loops.
   */
  public boolean isValid() {
    S2Error error = new S2Error();
    return !findValidationError(error);
  }

  /**
   * Returns true if this is *not* a valid polygon and sets {@code error} appropriately. Otherwise,
   * returns false and leaves {@code error} unchanged.
   */
  public boolean findValidationError(S2Error error) {
    for (int i = 0; i < numLoops(); ++i) {
      // Check for loop errors that don't require building an S2ShapeIndex.
      if (loop(i).findValidationErrorNoIndex(error)) {
        error.init(error.code(), "Loop " + i + ": " + error.text());
        return true;
      }
      // Check that no loop is empty, and that the full loop only appears in the full polygon.
      if (loop(i).isEmpty()) {
        error.init(S2Error.Code.POLYGON_EMPTY_LOOP,
            "Loop " + i + ": empty loops are not allowed.");
        return true;
      }
      if (loop(i).isFull() && numLoops() > 1) {
        error.init(S2Error.Code.POLYGON_EXCESS_FULL_LOOP,
            "Loop " + i + ": full loop appears in non-full polygon");
        return true;
      }
    }
    // Finally, check for loop self-intersections and loop pairs that cross (including duplicate
    // edges and vertices).
    //
    // TODO(user): Also verify the nesting hierarchy.
    return S2ShapeUtil.findAnyCrossing(index, loops, error);
  }

  public boolean isEmpty() {
    return loops.isEmpty();
  }

  public boolean isFull() {
    return loops.size() == 1 && loops.get(0).isFull();
  }

  public int numLoops() {
    return loops.size();
  }

  /**
   * Returns the loop at the given index. Note that during initialization, the given loops are
   * reordered according to a preorder traversal of the loop nesting hierarchy. This implies that
   * every loop is immediately followed by its descendants. This hierarchy can be traversed using
   * the methods {@link #getParent(int)}, {@link #getLastDescendant(int)}, and
   * {@link S2Loop#depth()}.
   */
  public S2Loop loop(int k) {
    return loops.get(k);
  }

  /** Returns the index of the parent of loop {@code k}, or -1 if it has no parent. */
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
   * Returns the index of the last loop that is contained within loop {@code k}. Returns {@code
   * numLoops() - 1} if {@code k < 0}. Note that loops are indexed according to a preorder traversal
   * of the nesting hierarchy, so the immediate children of loop {@code k} can be found by iterating
   * over loops {@code (k+1)..getLastDescendant(k)} and selecting those whose depth is equal to
   * {@code (loop(k).depth() + 1)}.
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
   * Returns the area of the polygon interior, i.e. the region on the left side of an odd number of
   * loops (the area is between 0 and 4*Pi) and the true centroid of the polygon, weighted by the
   * area of the polygon (see s2.h for details on centroids). Note that the centroid might not be
   * contained by the polygon.
   */
  public S2AreaCentroid getAreaAndCentroid() {
    return getAreaCentroid(true);
  }

  /**
   * Returns the area of the polygon interior, i.e. the region on the left side of an odd number of
   * loops. The return value is between 0 and 4*Pi.
   */
  public double getArea() {
    return getAreaCentroid(false).getArea();
  }

  /**
   * Returns the true centroid of the polygon, weighted by the area of the polygon (see s2.h for
   * details on centroids). Note that the centroid might not be contained by the polygon.
   */
  public S2Point getCentroid() {
    return getAreaCentroid(true).getCentroid();
  }

  /**
   * If all of the polygon's vertices happen to be the centers of S2Cells at some level, then
   * returns that level, otherwise returns -1. See also {@link #initToSnapped(S2Polygon, int)} and
   * {@link S2PolygonBuilder.Options.Builder#setSnapToCellCenters(boolean)}. Returns -1 if the
   * polygon has no vertices.
   */
  public int getSnapLevel() {
    int snapLevel = -1;
    for (S2Loop loop : loops) {
      for (int j = 0; j < loop.numVertices(); j++) {
        S2Point p = loop.vertex(j);
        FaceSiTi faceSiTi = S2Projections.PROJ.xyzToFaceSiTi(p);
        int level = S2Projections.PROJ.levelIfCenter(faceSiTi, p);
        if (level < 0) {
          // Vertex is not a cell center.
          return level;
        }
        if (level != snapLevel) {
          if (snapLevel < 0) {
            // First vertex.
            snapLevel = level;
          } else {
            // Vertices at more than one cell level.
            return -1;
          }
        }
      }
    }
    return snapLevel;
  }

  /**
   * Returns the shortest distance from a point P to this polygon, given as the angle formed between
   * P, the origin, and the nearest point on the polygon to P. This angle in radians is equivalent
   * to the arclength along the unit sphere.
   *
   * <p>If the point is contained inside the polygon, the distance returned is 0.
   */
  public S1Angle getDistance(S2Point p) {
    if (contains(p)) {
      return S1Angle.radians(0);
    }

    // The farthest point from p on the sphere is its antipode, which is an angle of PI radians.
    // This is an upper bound on the angle.
    S1Angle minDistance = S1Angle.radians(Math.PI);
    for (int i = 0; i < numLoops(); i++) {
      minDistance = S1Angle.min(minDistance, loop(i).getDistance(p));
    }

    return minDistance;
  }

  /**
   * Returns a point on the polygon that is closest to point P. The distance between these two
   * points should be the result of {@link #getDistance(S2Point)}.
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

    // The farthest point from p on the sphere is its antipode, which is an angle of PI radians.
    // This is an upper bound on the angle.
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
   * Returns true if this polygon contains the given other polygon, i.e., if polygon A contains all
   * points contained by polygon B.
   */
  public boolean contains(S2Polygon b) {
    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.contains does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).contains(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more efficient S2Loop.contains()
    // method (rather than CompareBoundary), but it's worthwhile to do our own bounds check first.
    if (!subregionBound.contains(b.getRectBound())) {
      // Even though bound(A) does not contain bound(B), it is still possible that A contains B.
      // This can only happen when the union of the two bounds spans all longitudes.  For example,
      // suppose that B consists of two shells with a longitude gap between them, while A consists
      // of one shell that surrounds both shells of B but goes the other way around the sphere (so
      // that it does not intersect the longitude gap).
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

    // Polygon A contains B iff B does not intersect the complement of A.  From the intersection
    // algorithm below, this means that the complement of A must exclude the entire boundary of B,
    // and B must exclude all shell boundaries of the complement of A.  (It can be shown that B must
    // then exclude the entire boundary of the complement of A.)  The first call below returns false
    // if the boundaries cross, therefore the second call does not need to check for any crossing
    // edges (which makes it cheaper).
    return containsBoundary(b) && b.excludesNonCrossingComplementShells(this);
  }

  /**
   * Returns true if this polygon (A) approximately contains the given other polygon (B). This is
   * true if it is possible to move the vertices of B no further than "vertexMergeRadius" such
   * that A contains the modified B.
   *
   * <p>For example, the empty polygon will contain any polygon whose maximum width is no more than
   * vertexMergeRadius.
   */
  public boolean approxContains(S2Polygon b, S1Angle vertexMergeRadius) {
    S2Polygon difference = new S2Polygon();
    difference.initToDifferenceSloppy(b, this, vertexMergeRadius);
    return difference.numLoops() == 0;
  }

  /**
   * Returns true if this polygon intersects the given other polygon, i.e., if there is a point that
   * is contained by both polygons.
   */
  public boolean intersects(S2Polygon b) {
    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.intersects does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).intersects(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more efficient S2Loop.intersects
    // method. The polygons intersect if and only if some pair of loop regions intersect.
    if (!bound.intersects(b.getRectBound())) {
      return false;
    }
    if (!hasHoles && !b.hasHoles) {
      for (S2Loop loop : b.loops) {
        if (anyLoopIntersects(loop)) {
          return true;
        }
      }
      return false;
    }

    // Polygon A is disjoint from B if A excludes the entire boundary of B and B excludes all shell
    // boundaries of A.  (It can be shown that B must then exclude the entire boundary of A.)  The
    // first call below returns false if the boundaries cross, therefore the second call does not
    // need to check for crossing edges.
    return !excludesBoundary(b) || !b.excludesNonCrossingShells(this);
  }

  /**
   * Indexing structure to efficiently {@link #clipEdge} of a polygon. This is an abstract class
   * because we need to use if for both polygons (for {@link #initToIntersection} and friends) and
   * for sets of lists of points (for initToSimplified() future?).
   *
   * <p>Usage: In your subclass, create an array of vertex counts for each loop in the loop sequence
   * and pass it to this constructor. Overwrite {@link #edgeFromTo}, calling {@link #decodeIndex}
   * and use the resulting two indices to access your vertices.
   */
  private abstract static class S2LoopSequenceIndex extends S2EdgeIndex {
    /** Map from the uni-dimensional edge index to the loop this edge belongs to. */
    private final int[] indexToLoop;

    /**
     * Reverse of {@link #indexToLoop}: maps a loop index to the uni-dimensional index of the first
     * edge in the loop.
     */
    private final int[] loopToFirstIndex;

    /**
     * Must be called by each subclass with the array of vertices per loop. The length of the array
     * is the number of loops, and the {@code i} <sup>th</sup> loop's vertex count is in the
     * {@code i} <sup>th</sup> index of the array.
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
     * Mark the {@link #edgeFromTo} method abstract again, so children of this class <b>must</b>
     * implement it without using {@link #edgeFrom(int)} and {@link #edgeTo(int)}.
     */
    @Override
    public abstract S2Edge edgeFromTo(int index);

    @Override
    public S2Point edgeFrom(int index) {
      return edgeFromTo(index).getStart();
    }

    @Override
    public S2Point edgeTo(int index) {
      return edgeFromTo(index).getEnd();
    }
  }

  /** Indexing structure for an {@link S2Polygon}. */
  public static final class S2PolygonIndex extends S2LoopSequenceIndex {
    private final S2Polygon poly;
    private final boolean reverse;

    /** Returns number of vertices per loop. */
    private static int[] getVertices(S2Polygon poly) {
      int[] vertices = new int[poly.numLoops()];
      for (int i = 0; i < vertices.length; i++) {
        vertices[i] = poly.loop(i).numVertices();
      }
      return vertices;
    }

    public S2PolygonIndex(S2Polygon poly) {
      this(poly, false);
    }

    S2PolygonIndex(S2Polygon poly, boolean reverse) {
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
   * Clips the boundary of A to the interior of B, and adds the resulting edges to {@code builder}.
   * Shells are directed CCW and holes are directed clockwise. If {@code reverseA} is true, these
   * directions are reversed in polygon A. If {@code invertB} is true, the boundary of A is clipped
   * to the exterior rather than the interior of B. If {@code addSharedEdges} is true, then the
   * output will include any edges that are shared between A and B (both edges must be in the same
   * direction after any edge reversals are taken into account).
   */
  private static void clipBoundary(final S2Polygon a, boolean reverseA, final S2Polygon b,
      boolean invertB, boolean addSharedEdges, S2PolygonBuilder builder) {
    EdgeClipper clipper = new EdgeClipper(b.index, addSharedEdges, true);
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    for (S2Loop aLoop : a.loops) {
      int n = aLoop.numVertices();
      int dir = (aLoop.isHole() ^ reverseA) ? -1 : 1;
      boolean inside = b.contains(aLoop.vertex(0)) ^ invertB;
      for (int j = (dir > 0) ? 0 : n; n > 0; --n, j += dir) {
        S2Point a0 = aLoop.vertex(j);
        S2Point a1 = aLoop.vertex(j + dir);
        clipper.clipEdge(a0,  a1, intersections);

        if (inside) {
          intersections.add(new ParametrizedS2Point(0, a0));
        }
        inside = ((intersections.size() & 0x1) == 0x1);
        // assert (b.contains(a1) ^ invertB == inside);
        if (inside) {
          intersections.add(new ParametrizedS2Point(1, a1));
        }

        Collections.sort(intersections);
        for (int k = 0; k < intersections.size(); k += 2) {
          S2Point x = intersections.get(k).getPoint();
          S2Point y = intersections.get(k + 1).getPoint();
          if (x.equalsPoint(y)) {
            continue;
          }
          builder.addEdge(x, y);
        }
        intersections.clear();
      }
    }
  }

  /** Returns the total number of vertices in all loops. */
  public int getNumVertices() {
    return this.numVertices;
  }

  /** Initializes this polygon to the complement of the given polygon. */
  public void initToComplement(S2Polygon a) {
    Preconditions.checkState(numLoops() == 0);
    copy(a);
    invert();
  }

  /**
   * Use S2PolygonBuilder to build this polygon by assembling the edges of a given polygon after
   * snapping its vertices to the center of leaf cells.  This will simplify the polygon with a
   * tolerance of {@code S2Projections.maxDiag.getValue(S2CellId.MAX_LEVEL)}, or approximately 0.13
   * microdegrees, or 1.5cm on the surface of the Earth.  Such a polygon can be efficiently
   * compressed when serialized.  The snap level can be changed to a non-leaf level if needed.
   */
  public void initToSnapped(final S2Polygon a, int snapLevel) {
    // TODO(user): Remove (tolerance * 0.1) from initToSimplified and use that instead.
    // Ensure that there will be no two vertices within the max leaf cell diagonal of each other,
    // therefore no two vertices in the same leaf cell, and that no vertex will cross an edge after
    // the points have been snapped to the centers of leaf cells.  Add 1e-15 to the tolerance so we
    // don't set a tighter than leaf cell level because of numerical inaccuracy.n
    S2PolygonBuilder.Options options = S2PolygonBuilder.Options.builder()
        .setRobustnessRadius(
            S1Angle.radians(PROJ.maxDiag.getValue(snapLevel) / 2.0 + 1e-15))
        .setSnapToCellCenters(true)
        .build();

    S2PolygonBuilder polygonBuilder = new S2PolygonBuilder(options);
    polygonBuilder.addPolygon(a);

    if (!polygonBuilder.assemblePolygon(this, null)) {
      log.severe("assemblePolygon failed in initToSnapped");
    }

    // If there are no loops, check whether the result should be the full
    // polygon rather than the empty one. (See InitToIntersectionSloppy.)
    if (numLoops() == 0) {
      if (a.bound.area() > 2.0 * S2.M_PI && a.getArea() > 2.0 * S2.M_PI) {
        invert();
      }
    }
  }

  /** Inverts this polygon (replacing it by its complement.) */
  private void invert() {
    // Inverting any one loop will invert the polygon.  The best loop to invert is the one whose
    // area is largest, since this yields the smallest area after inversion.  The loop with the
    // largest area is always at depth 0. The descendants of this loop all have their depth reduced
    // by 1, while the former siblings of this loop all have their depth increased by 1.

    // The empty and full polygons are handled specially.
    if (isEmpty()) {
      loops.add(S2Loop.full());
    } else if (isFull()) {
      clearLoops();
    } else {
      // Find the loop whose area is largest (i.e., whose turning angle is smallest), minimizing
      // calls to getTurningAngle().  In particular, for polygons with a single shell at level 0
      // there is not need to call GetTurningAngle() at all.  (This method is relatively expensive.)
      int best = -1;
      double bestAngle = 0;
      for (int i = 1; i < numLoops(); ++i) {
        S2Loop loop = loop(i);
        if (loop.depth() == 0) {
          // We defer computing the turning angle of loop 0 until we discover that the polygon has
          // another top-level shell.
          if (best == -1) {
            best = 0;
            bestAngle = loop(best).getTurningAngle();
          }
          double angle = loop.getTurningAngle();
          if (angle < bestAngle) {
            best = i;
            bestAngle = angle;
          }
        }
      }

      if (best < 0) {
        best = 0;
      }

      // Build the new loops vector, starting with the inverted loop.
      loop(best).invert();
      List<S2Loop> newLoops = Lists.newArrayListWithCapacity(numLoops());
      newLoops.add(loop(best));

      // Add the former siblings of this loop as descendants.
      int lastBest = getLastDescendant(best);
      for (int i = 0; i < numLoops(); ++i) {
        if (i < best || i > lastBest) {
          S2Loop loop = loop(i);
          loop.setDepth(loop.depth() + 1);
          newLoops.add(loop);
        }
      }

      // Add the former children of this loop as siblings.
      for (int i = 0; i < numLoops(); ++i) {
        if (i > best && i <= lastBest) {
          S2Loop loop = loop(i);
          loop.setDepth(loop.depth() - 1);
          newLoops.add(loop);
        }
      }
      Preconditions.checkState(loops.size() == newLoops.size());
      loops.clear();
      loops.addAll(newLoops);
    }

    initLoopProperties();
  }

  /**
   * Initializes this polygon to the intersection, union, or difference (A - B) of the given two
   * polygons. The {@code vertexMergeRadius} determines how close two vertices must be to be merged
   * together and how close a vertex must be to an edge in order to be spliced into it (see
   * {@link S2PolygonBuilder} for details). By default, the merge radius is just large enough to
   * compensate for errors that occur when computing intersection points between edges
   * ({@link S2EdgeUtil#DEFAULT_INTERSECTION_TOLERANCE}).
   *
   * <p>If you are going to convert the resulting polygon to a lower-precision format, it is
   * necessary to increase the merge radius in order to get a valid result after rounding (i.e., no
   * duplicate vertices, etc). For example, if you are going to convert them to
   * {@code geostore.PolygonProto} format, then {@code S1Angle.e7(1)} is a good value for
   * {@code vertexMergeRadius}.
   */
  public void initToIntersection(final S2Polygon a, final S2Polygon b) {
    initToIntersectionSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  public void initToIntersectionSloppy(
      final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    Preconditions.checkState(numLoops() == 0);
    if (a.isEmpty() || b.isEmpty()) {
      // From the check above we are already empty.
    } else if (a.isFull()) {
      copy(b);
    } else if (b.isFull()) {
      copy(a);
    } else {
      if (!a.bound.intersects(b.bound)) {
        return;
      }

      // We want the boundary of A clipped to the interior of B, plus the boundary of B clipped to
      // the interior of A, plus one copy of any directed edges that are in both boundaries.

      S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
          .setMergeDistance(vertexMergeRadius).build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      clipBoundary(a, false, b, false, true, builder);
      clipBoundary(b, false, a, false, false, builder);
      if (!builder.assemblePolygon(this, null)) {
        // TODO(andriy): do something more meaningful here.
        log.severe("Bad directed edges");
      }

      // If the result had a non-empty boundary then we are done.  Unfortunately, if the boundary is
      // empty then there are two possible results: the empty polygon or the full polygon.  This
      // choice would be trivial to resolve except for the existence of "vertex_merge_radius" and
      // also numerical errors when computing edge intersection points.  In particular:
      //
      //  - The intersection of two non-full polygons may be full.  For example, one or both
      //    polygons may have tiny cracks that are eliminated due to vertex merging/edge splicing.
      //
      //  - The intersection of two polygons that both contain S2.origin() (or any other point) may
      //    be empty.  For example, both polygons may have tiny shells that surround the common
      //    point but that are eliminated.
      //
      //  - Even before any vertex merging/edge splicing, the computed boundary edges are not useful
      //    in distinguishing almost-full polygons from almost-empty due to numerical errors in
      //    computing edge intersections. Such errors can reverse the orientation of narrow cracks
      //    or slivers.
      //
      // So instead we fall back to heuristics.  Essentially we compute the minimum and maximum
      // intersection area based on the areas of the two input polygons.  If only one of {0, 4*Pi}
      // is possible then we return that result.  If neither is possible (before vertex merging,
      // etc) then we return the one that is closest to being possible.  (It never true that
      // both are possible.)
      if (numLoops() == 0) {
        // We know that both polygons are non-empty due to the initial bounds check.  By far the
        // most common case is that the intersection is empty, so we want to make that case fast.
        // The intersection area satisfies:
        //
        //   max(0, A + B - 4*Pi) <= Intersection(A, B) <= min(A, B)
        //
        // where A, B refer to a polygon and/or its area.  Note that if either A or B is at most
        // 2*Pi, the result must be "empty".  We can use the bounding rectangle areas as upper
        // bounds on the polygon areas.
        if (a.bound.area() <= 2 * S2.M_PI || b.bound.area() <= 2 * S2.M_PI) {
          return;
        }
        double aArea = a.getArea();
        double bArea = b.getArea();
        double minArea = Math.max(0.0, aArea + bArea - 4 * S2.M_PI);
        double maxArea = Math.min(aArea, bArea);
        if (minArea > 4 * S2.M_PI - maxArea) {
          invert();
        }
      }
    }
  }

  public void initToUnion(final S2Polygon a, final S2Polygon b) {
    initToUnionSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  public void initToUnionSloppy(final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    Preconditions.checkState(numLoops() == 0);
    if (a.isEmpty() || b.isFull()) {
      copy(b);
    } else if (b.isEmpty() || a.isFull()) {
      copy(a);
    } else {
      // We want the boundary of A clipped to the exterior of B, plus the boundary of B clipped to
      // the exterior of A, plus one copy of any directed edges that are in both boundaries.
      S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
          .setMergeDistance(vertexMergeRadius)
          .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      clipBoundary(a, false, b, true, true, builder);
      clipBoundary(b, false, a, true, false, builder);
      if (!builder.assemblePolygon(this, null)) {
        // TODO(andriy): do something more meaningful here.
        log.severe("Bad directed edges");
      }

      if (numLoops() == 0) {
        // See comments in InitToIntersectionSloppy().  In this case, the union area satisfies:
        //
        //   max(A, B) <= Union(A, B) <= min(4*Pi, A + B)
        //
        // The most common case is that neither input polygon is empty, but the union is empty due
        // to vertex merging/simplification.
        if (a.bound.area() + b.bound.area() <= 2 * S2.M_PI) {
          return;
        }
        double aArea = a.getArea();
        double bArea = b.getArea();
        double minArea = Math.max(aArea, bArea);
        double maxArea = Math.min(4 * S2.M_PI, aArea + bArea);
        if (minArea > 4 * S2.M_PI - maxArea) {
          invert();
        }
      }
    }
  }

  public void initToDifference(final S2Polygon a, final S2Polygon b) {
    initToDifferenceSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  public void initToDifferenceSloppy(
      final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    Preconditions.checkState(numLoops() == 0);
    if (a.isEmpty() || b.isFull()) {
      // From the check above, this polygon is already empty.
    } else if (a.isFull()) {
      initToComplement(b);
    } else {
      // Note that we cannot short circuit the b.isEmpty() case because even something and nothing
      // might have no difference if the difference falls within the merge distance.

      // We want the boundary of A clipped to the exterior of B, plus the reversed boundary of B
      // clipped to the interior of A, plus one copy of any edge in A that is also a reverse edge in
      // B.
      S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
          .setMergeDistance(vertexMergeRadius)
          .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      clipBoundary(a, false, b, true, true, builder);
      clipBoundary(b, true, a, false, false, builder);
      if (!builder.assemblePolygon(this, null)) {
        log.severe("Bad directed edges");
      }

      if (numLoops() == 0) {
        // See comments in InitToIntersectionSloppy().  In this case, the difference area satisfies:
        //
        //   max(0, A - B) <= Difference(A, B) <= min(A, 4*Pi - B)
        //
        // By far the most common case is that result is empty.
        if (a.bound.area() <= 2 * S2.M_PI || b.bound.area() >= 2 * S2.M_PI) {
          return;
        }
        double aArea = a.getArea();
        double bArea = b.getArea();
        double minArea = Math.max(0.0, aArea - bArea);
        double maxArea = Math.min(aArea, 4 * S2.M_PI - bArea);
        if (minArea > 4 * S2.M_PI - maxArea) {
          invert();
        }
      }
    }
  }

  /**
   * Initializes this polygon to a polygon that contains fewer vertices and is within tolerance of
   * the polygon a, with some caveats.
   * <p>If {@code snapToCellCenters} is true, the vertices of this polygon will be snapped to the
   * centers of cells at the smallest level that is guaranteed to result in a valid polygon given
   * the specified tolerance.
   * <ul>
   * <li>If there is a very small island in the original polygon, it may disappear completely. Thus
   * some parts of the original polygon may not be close to the simplified polygon. Those parts are
   * small, though, and arguably don't need to be kept.
   * <li>However, if there are dense islands, they may all disappear, instead of replacing them by a
   * big simplified island.
   * <li>For small tolerances (compared to the polygon size), it may happen that the simplified
   * polygon has more vertices than the original, if the first step of the simplification creates
   * too many self-intersections. One can construct unrealistic cases where that happens to an
   * extreme degree.
   * </ul>
   */
  public void initToSimplified(S2Polygon a, S1Angle tolerance, boolean snapToCellCenters) {
    initToSimplifiedInternal(a, tolerance, snapToCellCenters, null);
  }

  /**
   * Initializes this polygon to a polygon that contains fewer vertices and is within tolerance of
   * the polygon a, while ensuring that the vertices on the cell boundary are preserved.
   *
   * <p>Precondition: Polygon a is contained in the cell.
   */
  public void initToSimplifiedInCell(S2Polygon a, final S2Cell cell, S1Angle tolerance) {
    initToSimplifiedInternal(a, tolerance, false, new S2VertexFilter() {
      /** Edges of the cell to test against. */
      private final S2Point[] corners = {
          cell.getVertex(0),
          cell.getVertex(1),
          cell.getVertex(2),
          cell.getVertex(3)};
      /** Max distance to test against. */
      private final S1Angle d = S1Angle.radians(1e-15);
      /** Returns true if the vertex is close to any edge of 'cell'. */
      @Override public boolean shouldKeepVertex(S2Point vertex) {
        for (int i = 0; i < 4; i++) {
          if (S2EdgeUtil.getDistance(vertex, corners[i], corners[(i + 1) % 4]).lessThan(d)) {
            return true;
          }
        }
        return false;
      }
    });
  }

  /**
   * Simplifies the polygon. The algorithm is straightforward and naive:
   * <ol>
   * <li>Simplify each loop by removing points while staying in the tolerance zone. This uses
   * {@link S2Polyline#subsampleVertices(S1Angle)}, which is not guaranteed to be optimal in terms
   * of number of points.
   * <li>Break any edge in pieces such that no piece intersects any other.
   * <li>Use the polygon builder to regenerate the full polygon.
   * <li>If {@code vertexFilter} is not null, the vertices for which it returns true are kept in the
   * simplified polygon.
   * </ol>
   */
  private void initToSimplifiedInternal(S2Polygon a, S1Angle tolerance, boolean snapToCellCenters,
      S2VertexFilter vertexFilter) {
    S2PolygonBuilder.Options.Builder options = S2PolygonBuilder.Options.UNDIRECTED_XOR.toBuilder();
    options.setValidate(false);
    if (vertexFilter != null) {
      // If there is a vertex filter, then we want to do as little vertex merging as possible so
      // that the vertices we want to keep don't move. But on the other hand, when we break
      // intersecting edges into pieces there is some error in the intersection point.
      // S2PolygonBuilder needs to be able to move vertices by up to this amount in order to produce
      // valid output.
      options.setMergeDistance(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
    } else {
      // Ideally, we would want to set the vertex merge radius of the builder roughly to tolerance
      // (and in fact forego the edge splitting step).  Alas, if we do that, we are liable to the
      // 'chain effect', where vertices are merged with close by vertices and so on, so that a
      // vertex can move by an arbitrary distance. So we remain conservative:
      options.setMergeDistance(S1Angle.radians(tolerance.radians() * 0.10));
      options.setSnapToCellCenters(snapToCellCenters);
    }

    // Simplify each loop separately and add to the edge index.
    S2PolygonBuilder builder = new S2PolygonBuilder(options.build());
    S2ShapeIndex index = new S2ShapeIndex();
    for (int i = 0; i < a.numLoops(); i++) {
      S2Loop simpler = a.loop(i).simplify(tolerance, vertexFilter);
      if (simpler != null) {
        index.add(simpler);
      }
    }

    if (!index.shapes.isEmpty()) {
      breakEdgesAndAddToBuilder(index, builder);
      if (!builder.assemblePolygon(this, null)) {
        log.warning("Bad edges in InitToSimplified.");
      }
      // If there are no loops, check whether the result should be the full
      // polygon rather than the empty one.  (See InitToIntersectionSloppy.)
      if (numLoops() == 0) {
        if (a.bound.area() > 2 * S2.M_PI && a.getArea() > 2 * S2.M_PI) {
          invert();
        }
      }
    }
  }

  /**
   * Takes a set of possibly intersecting edges, stored in the S2ShapeIndex, and breaks the edges
   * into small pieces so that there is no intersection anymore, and adds all these edges to the
   * builder.
   */
  public static void breakEdgesAndAddToBuilder(S2ShapeIndex index, S2PolygonBuilder builder) {
    // If there are self intersections, we add the pieces separately.
    // add_shared_edges ("true" below) can be false or true: it makes no
    // difference due to the way we call ClipEdge.
    EdgeClipper clipper = new EdgeClipper(index, true, false);
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    MutableEdge edge = new MutableEdge();
    for (S2Shape shape : index.getShapes()) {
      int numEdges = shape.numEdges();
      for (int e = 0; e < numEdges; e++) {
        shape.getEdge(e, edge);
        intersections.add(new ParametrizedS2Point(0, edge.a));
        clipper.clipEdge(edge.a, edge.b, intersections);
        intersections.add(new ParametrizedS2Point(1, edge.b));
        Collections.sort(intersections);
        for (int k = 0; k + 1 < intersections.size(); ++k) {
          S2Point p1 = intersections.get(k).getPoint();
          S2Point p2 = intersections.get(k + 1).getPoint();
          if (!p1.equalsPoint(p2)) {
            builder.addEdge(p1, p2);
          }
        }
        intersections.clear();
      }
    }
  }

  /** Returns a polygon that is the union of the given polygons. Note: clears the List! */
  public static S2Polygon destructiveUnion(List<S2Polygon> polygons) {
    return destructiveUnionSloppy(polygons, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Returns a polygon that is the union of the given polygons; combines vertices that form edges
   * that are almost identical, as defined by {@code vertexMergeRadius}. Note: clears the List!
   */
  public static S2Polygon destructiveUnionSloppy(
      List<S2Polygon> polygons, S1Angle vertexMergeRadius) {
    // Effectively create a priority queue of polygons in order of number of vertices. Repeatedly
    // union the two smallest polygons and add the result to the queue until we have a single
    // polygon to return.

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
      // We assume that the number of vertices in the union polygon is the sum of the number of
      // vertices in the original polygons, which is not always true, but will almost always be a
      // decent approximation, and faster than recomputing.
    }

    if (queue.isEmpty()) {
      return new S2Polygon();
    } else {
      return queue.get(queue.asMap().firstKey()).first();
    }
  }

  /**
   * Intersects this polygon with the {@link S2Polyline} {@code in} and returns the resulting zero
   * or more polylines. The polylines are ordered in the order they would be encountered by
   * traversing {@code in} from beginning to end. Note that the output may include polylines with
   * only one vertex, but there will not be any zero-vertex polylines.
   *
   * <p>This is equivalent to calling {@link #intersectWithPolylineSloppy} with the {@code
   * vertexMergeRadius} set to {@link S2EdgeUtil#DEFAULT_INTERSECTION_TOLERANCE}.
   */
  public List<S2Polyline> intersectWithPolyline(S2Polyline in) {
    return intersectWithPolylineSloppy(in, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Similar to {@link #intersectWithPolyline}, except that vertices will be dropped as necessary to
   * ensure that all adjacent vertices in the sequence obtained by concatenating the output
   * polylines will be farther than {@code vertexMergeRadius} apart. Note that this can change the
   * number of output polylines and/or yield single-vertex polylines.
   */
  public List<S2Polyline> intersectWithPolylineSloppy(
      S2Polyline in, S1Angle vertexMergeRadius) {
    return internalClipPolyline(false, in, vertexMergeRadius);
  }

  /** Same as {@link #intersectWithPolyline}, but subtracts this polygon from the given polyline. */
  public List<S2Polyline> subtractFromPolyline(S2Polyline in) {
    return subtractFromPolylineSloppy(in, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Same as {@link #intersectWithPolylineSloppy}, but subtracts this polygon from the given
   * polyline.
   */
  public List<S2Polyline> subtractFromPolylineSloppy(
      S2Polyline in, S1Angle vertexMergeRadius) {
    return internalClipPolyline(true, in, vertexMergeRadius);
  }

  /**
  * Clips the {@link S2Polyline} {@code a} to the interior of this polygon. The resulting
  * polyline(s) will be returned. If {@code invert} is {@code true}, we clip {@code a} to the
  * exterior of this polygon instead. Vertices will be dropped such that adjacent vertices will not
  * be closer than {@code mergeRadius}.
  *
  * <p>We do the intersection/subtraction by walking the polyline edges. For each edge, we compute
  * all intersections with the polygon boundary and sort them in increasing order of distance along
  * that edge. We then divide the intersection points into pairs, and output a clipped polyline
  * segment for each one. We keep track of whether we're inside or outside of the polygon at all
  * times to decide which segments to output.
  */
  private List<S2Polyline> internalClipPolyline(
      boolean invert, S2Polyline a, S1Angle mergeRadius) {

    List<S2Polyline> out = Lists.newArrayList();
    EdgeClipper clipper = new EdgeClipper(index, true, false);
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    List<S2Point> vertices = Lists.newArrayList();
    int n = a.numVertices();
    boolean inside = contains(a.vertex(0)) ^ invert;
    for (int j = 0; j < n - 1; j++) {
      S2Point a0 = a.vertex(j);
      S2Point a1 = a.vertex(j + 1);
      clipper.clipEdge(a0, a1, intersections);
      if (inside) {
        intersections.add(new ParametrizedS2Point(0, a0));
      }
      inside = (intersections.size() & 1) != 0;
      // assert ((contains(a1) ^ invert) == inside);
      if (inside) {
        intersections.add(new ParametrizedS2Point(1, a1));
      }
      Collections.sort(intersections);
      // At this point we have a sorted array of vertex pairs representing the edge(s) obtained
      // after clipping (a0, a1) against the polygon.
      for (int k = 0; k < intersections.size(); k += 2) {
        S2Point v0 = intersections.get(k).getPoint();
        S2Point v1 = intersections.get(k + 1).getPoint();
        if (!v0.equalsPoint(v1)) {
          // If the gap from the previous vertex to this one is large enough, start a new polyline.
          if (!vertices.isEmpty()
              && vertices.get(vertices.size() - 1).angle(v0) > mergeRadius.radians()) {
            out.add(new S2Polyline(vertices));
            vertices.clear();
          }
          // Append this segment to the current polyline, ignoring any vertices that are too close
          // to the previous vertex.
          if (vertices.isEmpty()) {
            vertices.add(v0);
          }
          if (vertices.get(vertices.size() - 1).angle(v1) > mergeRadius.radians()) {
            vertices.add(v1);
          }
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
    // TODO(user): The condition tested here is insufficient.  The correct condition is that each
    // *connected component* of child loops can share at most one vertex with their parent loop.
    // Example: suppose loop A has children B, C, D, and the following pairs are connected: AB, BC,
    // CD, DA. Then the polygon is not normalized.
    Set<S2Point> vertices = Sets.newHashSet();
    S2Loop lastParent = null;
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop child = loop(i);
      if (child.depth() == 0) {
        continue;
      }
      S2Loop parent = loop(getParent(i));
      // Test if the loops are different; we can use identity since any two loops in a valid polygon
      // can only be equal by value if they're equal by reference.
      if (parent != lastParent) {
        vertices.clear();
        for (int j = 0; j < parent.numVertices(); ++j) {
          vertices.add(parent.vertex(j));
        }
        lastParent = parent;
      }
      int count = 0;
      for (int j = 0; j < child.numVertices(); ++j) {
        if (vertices.contains(child.vertex(j))) {
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
   * Returns true if two polygons have the same boundary, except for vertex perturbations. Both
   * polygons must have loops with the same cyclic vertex order and the same nesting hierarchy, but
   * the vertex locations are allowed to differ by up to {@code maxError}. Note: This method mostly
   * useful only for testing purposes.
   */
  boolean boundaryApproxEquals(S2Polygon b, double maxError) {
    if (numLoops() != b.numLoops()) {
      log.severe(
          "!= loops: " + Integer.toString(numLoops()) + " vs. " + Integer.toString(b.numLoops()));
      return false;
    }

    // For now, we assume that there is at most one candidate match for each loop. (So far this
    // method is just used for testing.)
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

  /**
   * Returns true if two polygons have boundaries that are within {@code maxError} of each other
   * along their entire lengths. More precisely, there must be a bijection between the two sets of
   * loops such that {@code aLoop.boundaryNear(bLoop)} is true for each aLoop in {@code this} and
   * each bLoop in {@code b}.
   */
  boolean boundaryNear(S2Polygon b, double maxError) {
    if (numLoops() != b.numLoops()) {
      return false;
    }

    // For now, we assume that there is at most one candidate match for each loop.
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop aLoop = loop(i);
      boolean success = false;
      for (int j = 0; j < numLoops(); ++j) {
        S2Loop bLoop = b.loop(j);
        if (bLoop.depth() == aLoop.depth() && bLoop.boundaryNear(aLoop, maxError)) {
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

  /**
   * Returns a spherical cap that bounds this loop. It may be expanded slightly such that if the
   * loop contains a point P, then the bound contains P also.
   */
  @Override
  public S2Cap getCapBound() {
    return bound.getCapBound();
  }

  /**
   * Returns a fairly tight bounding latitude-longitude rectangle. It is not guaranteed to be as
   * tight as possible, to ensure that if the loop contains a point P, then the bound contains P
   * also.
   */
  @Override
  public S2LatLngRect getRectBound() {
    return bound;
  }

  /**
   * If this method returns true, the region completely contains the given cell. Otherwise, either
   * the region does not contain the cell or the containment relationship could not be determined.
   */
  @Override
  public boolean contains(S2Cell cell) {
    S2ShapeIndex.CellIterator it = index.iterator();
    S2ShapeIndex.CellRelation relation = it.locate(cell.id());

    // If 'cell' is disjoint from all index cells, it is not contained. Similarly, if 'cell' is
    // subdivided into one or more index cells, then it is not contained, since index cells are
    // subdivided only if they (nearly) intersect a sufficient number of edges. (But note that if
    // 'cell' itself is an index cell, then it may be contained, since it could be a cell with no
    // indexed edges in the polygon interior.)
    if (relation != S2ShapeIndex.CellRelation.INDEXED) {
      return false;
    }

    // Otherwise check if any edges intersect 'cell'. At this point, the iterator is guaranteed to
    // to point to an index cell containing 'cell'.
    if (boundaryApproxIntersects(it, cell)) {
      return false;
    }

    // Otherwise check if the polygon contains the center of 'cell'.
    return contains(it, cell.getCenter());
  }

  /**
   * If this method returns false, the region does not intersect the given cell. Otherwise, either
   * region intersects the cell, or the intersection relationship could not be determined.
   */
  @Override
  public boolean mayIntersect(S2Cell target) {
    S2ShapeIndex.CellIterator it = index.iterator();
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If 'target' does not overlap any index cell, there is no intersection.
    if (relation == S2ShapeIndex.CellRelation.DISJOINT) {
      return false;
    }

    // If 'target' properly contains an index cell, then there is an intersection to within the
    // S2ShapeIndex error bound (see contains). But if 'target' itself is an index cell, then it may
    // not be contained, since it could be a cell with no indexed edges that is contained by both a
    // shell and a hole, and is therefore in the polygon exterior.
    if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) {
      return true;
    }

    // Otherwise check if any edges intersect 'target'.
    if (boundaryApproxIntersects(it, target)) {
      return true;
    }

    // Otherwise check if the polygon contains the center of 'target'.
    return contains(it, target.getCenter());
  }

  /**
   * Returns true if the polygon boundary intersects {@code target}. It may also return true when
   * the polygon boundary does not intersect {@code target} but some edge comes within the
   * worst-case error tolerance.
   *
   * Requires: {@code it.id().contains(target.id())}
   * (This condition is true whenever {@code it.locate(target)} returns INDEXED.
   */
  private boolean boundaryApproxIntersects(S2ShapeIndex.CellIterator it, S2Cell target) {
    // assert (it.id().contains(target.id()))
    S2ShapeIndex.Cell cell = it.cell();
    R2Rect bound = target.getBoundUV().expanded(S2EdgeUtil.MAX_CELL_EDGE_ERROR);
    R2Vector v0 = new R2Vector();
    R2Vector v1 = new R2Vector();
    for (int a = 0; a < cell.numShapes(); ++a) {
      S2ShapeIndex.S2ClippedShape aClipped = cell.clipped(a);
      int aNumClipped = aClipped.numEdges();
      if (aNumClipped == 0) {
        continue;
      }

      // We can save some work if 'target' is the index cell itself (given that there is at least
      // one indexed edge).
      if (it.id().equals(target.id())) {
        return true;
      }

      // Otherwise, check whether any of the edges intersect 'target'.
      S2Loop aLoop = (S2Loop) aClipped.shape();
      for (int i = 0; i < aNumClipped; ++i) {
        int ai = aClipped.edge(i);
        if (S2EdgeUtil.clipToPaddedFace(aLoop.vertex(ai), aLoop.vertex(ai + 1), target.face(),
            S2EdgeUtil.MAX_CELL_EDGE_ERROR, v0, v1) && S2EdgeUtil.intersectsRect(v0, v1, bound)) {
          return true;
        }
      }
    }
    return false;
  }

  /** The point {@code p} does not need to be normalized. */
  public boolean contains(S2Point p) {
    // Note(user): A bounds check slows down the C++ version of this method by about 50%. It is
    // worthwhile only when it might allow us to delay building the index.
    if (!index.isFresh() && !bound.contains(p)) {
      return false;
    }

    // For small polygons it is faster to just check all the crossings. Otherwise we keep track of
    // the number of calls to contains() and only build the index once enough calls have been made
    // so that we think it is worth the effort. See S2Loop.contains(S2Point) for detailed comments.
    int maxBruteForceVertices = 32;
    if (getNumVertices() <= maxBruteForceVertices || (!index.isFresh()
        && unindexedContainsCalls.decrementAndGet() > 0)) {
      boolean inside = false;
      for (int i = 0; i < numLoops(); ++i) {
        // Use brute force to avoid building the loops' S2ShapeIndex.
        boolean loopContainsP = loop(i).bruteForceContains(p);
        // Set 'inside' to the bitwise XOR of 'inside' and 'loopContainsP'.
        inside ^= loopContainsP;
      }
      return inside;
    }

    // Otherwise we look up the S2ShapeIndex cell containing this point.
    S2ShapeIndex.CellIterator it = index.iterator();
    if (!it.locate(p)) {
      return false;
    }
    return contains(it, p);
  }

  /**
   * Given an iterator that is already positioned at the S2ShapeIndex.Cell containing {@code p},
   * return {@code contains(p)}.
   */
  private boolean contains(S2ShapeIndex.CellIterator it, S2Point p) {
    // Test containment by drawing a line segment from the cell center to the given point and
    // counting edge crossings.
    S2ShapeIndex.Cell cell = it.cell();
    boolean inside = false;
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(it.center(), p);
    for (int a = 0; a < cell.numShapes(); ++a) {
      S2ShapeIndex.S2ClippedShape aClipped = cell.clipped(a);
      inside ^= aClipped.containsCenter();
      int aNumClipped = aClipped.numEdges();
      if (aNumClipped > 0) {
        int aiPrev = -2;
        S2Loop aLoop = (S2Loop) aClipped.shape();
        for (int i = 0; i < aNumClipped; ++i) {
          int ai = aClipped.edge(i);
          if (ai != aiPrev + 1) {
            crosser.restartAt(aLoop.vertex(ai));
          }
          aiPrev = ai;
          inside ^= crosser.edgeOrVertexCrossing(aLoop.vertex(ai + 1));
        }
      }
    }
    return inside;
  }

  /** For each map entry, sorts the value list. */
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

    // No loop may contain the complement of another loop. (Handling this case is significantly more
    // complicated.)
    // assert (parent == null || !newLoop.containsNested(parent));

    // Some of the children of the parent loop may now be children of the new loop.
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

  /**
   * Returns +1 if this polygon (A) contains the boundary of B, -1 if A excludes the boundary of B,
   * and 0 if the boundaries of A and B cross.
   */
  int compareBoundary(S2Loop b) {
    int result = -1;
    for (S2Loop loop : loops) {
      // If B crosses any loop of A, the result is 0. Otherwise the result changes sign each time B
      // is contained by a loop of A.
      result *= -loop.compareBoundary(b);
      if (result == 0) {
        break;
      }
    }
    return result;
  }

  /** Returns true if this polygon (A) contains the entire boundary of B. */
  private boolean containsBoundary(S2Polygon b) {
    for (S2Loop loop : b.loops) {
      if (compareBoundary(loop) <= 0) {
        return false;
      }
    }
    return true;
  }

  /** Returns true if this polygon (A) excludes the entire boundary of B. */
  private boolean excludesBoundary(S2Polygon b) {
    for (S2Loop loop : b.loops) {
      if (compareBoundary(loop) >= 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Given a polygon A and a loop B whose boundaries do not cross, returns true if A contains the
   * boundary of B. Shared edges are handled according to the rule described in
   * {@link S2Loop#containsNonCrossingBoundary(S2Loop, boolean)}.
   */
  private boolean containsNonCrossingBoundary(S2Loop b, boolean bReverse) {
    boolean inside = false;
    for (S2Loop loop : loops) {
      inside ^= loop.containsNonCrossingBoundary(b, bReverse);
    }
    return inside;
  }

  /**
   * Given two polygons A and B such that the boundary of A does not cross any loop of B, returns
   * true if A excludes all shell boundaries of B.
   */
  private boolean excludesNonCrossingShells(S2Polygon b) {
    for (S2Loop loop : b.loops) {
      if (!loop.isHole()) {
        if (containsNonCrossingBoundary(loop, false)) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Given two polygons A and B such that the boundary of A does not cross any loop of B, returns
   * true if A excludes all shell boundaries of the complement of B.
   */
  private boolean excludesNonCrossingComplementShells(S2Polygon b) {
    // Special case to handle the complement of the empty or full polygons.
    if (b.isEmpty()) {
      return !isFull();
    }

    if (b.isFull()) {
      return true;
    }

    // Otherwise the complement of B may be obtained by inverting loop(0) and then swapping the
    // shell/hole status of all other loops.  This implies that the shells of the complement consist
    // of loop 0 plus all the holes of the original polygon.
    for (int j = 0; j < b.numLoops(); ++j) {
      if (j > 0 && !b.loop(j).isHole()) {
        continue;
      }

      // The interior of the complement is to the right of loop 0, and to the left of the loops that
      // were originally holes.
      if (containsNonCrossingBoundary(b.loop(j), j == 0)) {
        return false;
      }
    }

    return true;
  }

  /** Returns true if any loop contains the given loop. */
  private boolean anyLoopContains(S2Loop b) {
    for (S2Loop loop : loops) {
      if (loop.contains(b)) {
        return true;
      }
    }
    return false;
  }

  /** Returns true if any loop intersects the given loop. */
  private boolean anyLoopIntersects(S2Loop b) {
    for (S2Loop loop : loops) {
      if (loop.intersects(b)) {
        return true;
      }
    }
    return false;
  }

  /** Returns a human readable representation of the polygon. */
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
        // end of vertex
        sb.append("\n");
      }
      // end of loop
      sb.append(">\n");
    }
    return sb.toString();
  }

  @GwtCompatible(serializable = false)
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

  /**
   * EdgeClipper finds all the intersections of a given edge with the edges contained in an
   * S2ShapeIndex. It is used to implement polygon operations such as intersection and union.
   */
  private static class EdgeClipper {
    private final S2EdgeQuery query;
    private final boolean addSharedEdges;
    private final boolean reverseEdges;

    public EdgeClipper(S2ShapeIndex index, boolean addSharedEdges, boolean reverseEdges) {
      query = new S2EdgeQuery(index);
      this.addSharedEdges = addSharedEdges;
      this.reverseEdges = reverseEdges;
    }

    /**
     * Finds all points where the polygon B intersects the edge (a0, a1), and add the corresponding
     * parameter values (in the range [0,1]) to {@code intersections}. The result is unsorted.
     */
    public void clipEdge(S2Point a0, S2Point a1, List<ParametrizedS2Point> intersections) {
      Map<S2Shape, S2EdgeQuery.Edges> edgeMap = query.getCandidates(a0, a1);
      if (edgeMap.isEmpty()) {
        return;
      }

      // Iterate through the candidate loops, and then the candidate edges within each loop.
      S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a0, a1);
      MutableEdge result = new MutableEdge();
      for (Entry<S2Shape, S2EdgeQuery.Edges> entry : edgeMap.entrySet()) {
        final S2Loop bLoop = (S2Loop) entry.getKey();
        S2EdgeQuery.Edges edges = entry.getValue();
        int b1Prev = -2;
        while (!edges.isEmpty()) {
          int edge = edges.nextEdge();
          bLoop.getEdge(edge, result);
          if (edge != b1Prev + 1) {
            crosser.restartAt(result.getStart());
          }
          int crossing = crosser.robustCrossing(result.getEnd());
          if (crossing >= 0) {
            addIntersection(a0, a1, result.getStart(), result.getEnd(), bLoop, crossing,
                intersections);
            b1Prev = edge;
          }
        }
      }
    }

    /**
     * Given two edges A and B such that robustCrossing(A, B) >= 0, determines if they intersect and
     * adds any intersection point to {@code intersections}. {@code bShape} is the S2Shape
     * containing edge B, and {@code crossing} is the result of {@code robustCrossing(A, B)}.
     */
    private void addIntersection(S2Point a0, S2Point a1, S2Point b0, S2Point b1,
        S2Loop bLoop, int crossing, List<ParametrizedS2Point> intersections) {
      // assert (crossing >= 0);
      if (crossing > 0) {
        // There is a proper edge crossing.
        S2Point x = S2EdgeUtil.getIntersection(a0, a1, b0, b1);
        double t = S2EdgeUtil.getDistanceFraction(x, a0, a1);
        intersections.add(new ParametrizedS2Point(t, x));
      } else if (S2EdgeUtil.vertexCrossing(a0, a1, b0, b1)) {
        // There is a crossing at one of the vertices. The basic rule is simple: if a0 equals one of
        // the 'b' vertices, the crossing occurs at t = 0; otherwise, it occurs at t = 1.
        //
        // This has the effect that when two reversed edges exist (i.e., a0 == b1 and a1 == b0) and
        // each edge is clipped against the other, neither one is included in the output (which is
        // correct). However, when two shared edges exist (i.e., a0 == b0 and a1 == b1), both are
        // included in the output (which is incorrect). The 'addSharedEdges' flag gives explicit
        // over whether these shared edges are included in the output; if it is false, shared edges
        // are excluded by changing their intersection parameter from 0 to 1. This allows exactly
        // one copy of shared edges to be preserved, by calling ClipBoundary() twice with opposite
        // values of the 'addSharedEdges' flag.
        double t = (a0.equalsPoint(b0) || a0.equalsPoint(b1)) ? 0 : 1;
        S2Point p = (reverseEdges && bLoop.isHole()) ? b0 : b1;
        if (!addSharedEdges && a1.equalsPoint(p)) {
          // Excludes (a0, a1) from the output.
          t = 1;
        }
        intersections.add(new ParametrizedS2Point(t, t == 0 ? a0 : a1));
      }
    }
  }
}
