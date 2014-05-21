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

import com.google.common.annotations.GwtCompatible;
import com.google.common.base.Objects;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.geometry.S2EdgeIndex.DataEdgeIterator;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * An S2Loop represents a simple spherical polygon. It consists of a single chain of vertices where
 * the first vertex is implicitly connected to the last. All loops are defined to have a CCW
 * orientation, i.e. the interior of the loop is on the left side of the edges. This implies that a
 * clockwise loop enclosing a small area is interpreted to be a CCW loop enclosing a very large
 * area.
 *
 * <p>Loops are not allowed to have any duplicate vertices (whether adjacent or not), and non-
 * adjacent edges are not allowed to intersect. Loops must have at least 3 vertices (except for the
 * "empty" and "full" loops discussed below). Although these restrictions are not enforced in
 * optimized code, you may get unexpected results if they are violated.
 *
 * <p>There are two special loops: the "empty" loop contains no points, while the "full" loop
 * contains all points. These loops do not have any edges, but to preserve the invariant that every
 * loop can be represented as a vertex chain, they are defined as having exactly one vertex each (
 * {@link #empty()} and {@link #full()}.)
 *
 * <p>Point containment of loops is defined such that if the sphere is subdivided into faces
 * (loops), every point is contained by exactly one face. This implies that loops do not necessarily
 * contain their vertices.
 *
 */
@GwtCompatible(serializable = true)
public final strictfp class S2Loop implements S2Region, Comparable<S2Loop>, Serializable, S2Shape {
  private static final Logger log = Platform.getLoggerForClass(S2Loop.class);

  /**
   * Max angle that intersections can be off by and yet still be considered
   * colinear.
   */
  public static final double MAX_INTERSECTION_ERROR = 1e-15;

  /**
   * The single vertex that defines a loop that contains no area.
   */
  static final S2Point EMPTY_VERTEX = S2Point.Z_POS;

  /**
   * The single vertex that defines a loop that contains the whole sphere.
   */
  static final S2Point FULL_VERTEX = S2Point.Z_NEG;

  /**
   * Edge index used for performance-critical operations. For example,
   * contains() can determine whether a point is inside a loop in nearly
   * constant time, whereas without an edge index it is forced to compare the
   * query point against every edge in the loop.
   */
  private transient S2EdgeIndex index;

  /** Maps each S2Point to its order in the loop, from 1 to numVertices. */
  private Map<S2Point, Integer> vertexToIndex;

  private final S2Point[] vertices;
  private final int numVertices;

  /**
   * The index (into "vertices") of the vertex that comes first in the total
   * ordering of all vertices in this loop.
   */
  private int firstLogicalVertex;

  /**
   * A conservative bound on all points contained by this loop: if A.contains(P), then
   * A.bound.contains(new S2LatLng(P)).
   */
  private S2LatLngRect bound;

  /**
   * Since "bound" is not exact, it is possible that a loop A contains another loop B whose bounds
   * are slightly larger. "subregionBound" has been expanded sufficiently to account for this error,
   * i.e. if A.contains(B), then A.subregionBound.contains(B.bound).
   */
  private S2LatLngRect subregionBound;

  private boolean originInside;
  private int depth;

  /**
   * Initializes a loop with the given vertices. The last vertex is implicitly connected to the
   * first. All points should be unit length. Loops must have at least 3 vertices (except for the
   * "empty" and "full" loops; see {@link #empty()} and {@link #full()}.
   *
   * @param vertices the vertices for this new loop
   */
  public S2Loop(final List<S2Point> vertices) {
    this.numVertices = vertices.size();
    this.vertices = new S2Point[numVertices];
    vertices.toArray(this.vertices);
    this.depth = 0;
    initOriginAndBound();
    initFirstLogicalVertex();
  }

  /**
   * Fast/unsafe loop initialization.
   *
   * <p>
   * This constructor provides known good values for bounds and the originInside value. This is
   * intended to be a "fast loop creation" when we already know a lot about the loop. It is
   * primarily used in combination with the fast S2Polygon initializer (
   * {@link S2Polygon#initWithNestedLoops(java.util.Map)}). The last vertex is implicitly connected
   * to the first. All points should be unit length. Loops must have at least 3 vertices, except for
   * the empty and full loops (see {@link #empty()} and {@link #full()}.)
   *
   * @param vertices loop vertices
   * @param originInside true if the S2::origin() is inside the loop
   * @param bound the lat/long bounds of the loop
   * @return new loop.
   */
  public static S2Loop newLoopWithTrustedDetails(
      List<S2Point> vertices, boolean originInside, S2LatLngRect bound) {
    // This is a static method to try and discourage its use.
    return new S2Loop(vertices, originInside, bound);
  }

  /**
   * Create a circle of points with a given center, radius, and number of vertices.
   */
  public static S2Loop makeRegularLoop(S2Point center, S1Angle radius, int numVertices) {
    Matrix3x3 m = Matrix3x3.fromCols(S2.getFrame(center));
    List<S2Point> vertices = Lists.newArrayList();
    double radianStep = 2 * Math.PI / numVertices;
    // We create the vertices on the plane tangent to 'center', so the radius on that plane is
    // larger.
    double planarRadius = Math.tan(radius.radians());
    for (int vi = 0; vi < numVertices; ++vi) {
      double angle = vi * radianStep;
      S2Point p = new S2Point(planarRadius * Math.cos(angle), planarRadius * Math.sin(angle), 1);
      vertices.add(S2Point.normalize(S2.rotate(p, m)));
    }
    return new S2Loop(vertices);
  }

  private S2Loop(List<S2Point> vertices, boolean originInside, S2LatLngRect bound) {
    this.numVertices = vertices.size();
    this.vertices = new S2Point[numVertices];
    this.bound = bound;
    this.subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(bound);
    this.depth = 0;
    this.originInside = originInside;

    vertices.toArray(this.vertices);

    initFirstLogicalVertex();
  }

  /**
   * Initialize a loop corresponding to the given cell.
   */
  public S2Loop(S2Cell cell) {
    numVertices = 4;
    vertices = new S2Point[numVertices];
    vertexToIndex = null;
    index = null;
    depth = 0;
    for (int i = 0; i < 4; ++i) {
      vertices[i] = cell.getVertex(i);
    }
    // We compute the bounding rectangle ourselves, since S2Cell uses a different method and we need
    // all the bounds to be consistent.
    initOriginAndBound();
    initFirstLogicalVertex();
  }

  /**
   * Copy constructor.
   */
  public S2Loop(S2Loop src) {
    this.numVertices = src.numVertices();
    this.vertices = new S2Point[numVertices];
    for (int i = 0; i < numVertices; i++) {
      this.vertices[i] = src.vertices[i];
    }
    this.vertexToIndex = src.vertexToIndex;
    this.index = src.index;
    this.firstLogicalVertex = src.firstLogicalVertex;
    this.bound = src.getRectBound();
    this.subregionBound = src.subregionBound;
    this.originInside = src.originInside;
    this.depth = src.depth();
  }

  /**
   * Returns a new loop with one vertex that defines an empty loop (i.e., a loop with no edges
   * that contains no points.)
   */
  public static final S2Loop empty() {
    return new S2Loop(Collections.singletonList(EMPTY_VERTEX));
  }

  /**
   * Returns a new loop with one vertex that creates a full loop (i.e., a loop with no edges
   * that contains all points).  See {@link #empty()} for further details.
   */
  public static final S2Loop full() {
    return new S2Loop(Collections.singletonList(FULL_VERTEX));
  }

  // Note that this doesn't do anything smart: it just compares a few fields for equality, but
  // doesn't check for equivalent loops that were initialized in different ways, etc.
  @Override
  public boolean equals(Object obj) {
    if (obj instanceof S2Loop) {
      S2Loop that = (S2Loop) obj;
      return Arrays.equals(this.vertices, that.vertices)
          && Objects.equal(this.originInside, that.originInside)
          && Objects.equal(this.bound, that.bound);
    }
    return false;
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(this.vertices, this.originInside, this.bound);
  }

  public int depth() {
    return depth;
  }

  /**
   * The depth of a loop is defined as its nesting level within its containing
   * polygon. "Outer shell" loops have depth 0, holes within those loops have
   * depth 1, shells within those holes have depth 2, etc. This field is only
   * used by the S2Polygon implementation.
   *
   * @param depth
   */
  public void setDepth(int depth) {
    this.depth = depth;
  }

  /**
   * Return true if this loop represents a hole in its containing polygon.
   */
  public boolean isHole() {
    return (depth & 1) != 0;
  }

  /**
   * The sign of a loop is -1 if the loop represents a hole in its containing
   * polygon, and +1 otherwise.
   */
  public int sign() {
    return isHole() ? -1 : 1;
  }

  public int numVertices() {
    return numVertices;
  }

  /**
   * For convenience, we make two entire copies of the vertex list available:
   * vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == numVertices().
   */
  public S2Point vertex(int i) {
    try {
      return vertices[i >= vertices.length ? i - vertices.length : i];
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new IllegalStateException("Invalid vertex index");
    }
  }

  /**
   * Returns true if this is the special "empty" loop that contains no points.
   */
  public boolean isEmpty() {
    return isEmptyOrFull() && !originInside;
  }

  /**
   * Returns true if this is the special "full" loop that contains all points.
   */
  public boolean isFull() {
    return isEmptyOrFull() && originInside;
  }

  /** Returns true if this loop is either "empty" or "full". */
  public boolean isEmptyOrFull() {
    return numVertices == 1;
  }

  @Override
  public int numEdges() {
    return isEmptyOrFull() ? 0 : numVertices;
  }

  @Override
  public void getEdge(int index, MutableEdge result) {
    result.set(vertex(index), vertex(index + 1));
  }

  @Override
  public boolean hasInterior() {
    return !isEmpty();
  }

  @Override
  public boolean containsOrigin() {
    return originInside;
  }

  /**
   * Comparator (needed by Comparable interface)
   */
  @Override
  public int compareTo(S2Loop other) {
    if (numVertices() != other.numVertices()) {
      return this.numVertices() - other.numVertices();
    }
    // Compare the two loops' vertices, starting with each loop's
    // firstLogicalVertex. This allows us to always catch cases where logically
    // identical loops have different vertex orderings (e.g. ABCD and BCDA).
    int maxVertices = numVertices();
    int iThis = firstLogicalVertex;
    int iOther = other.firstLogicalVertex;
    for (int i = 0; i < maxVertices; ++i, ++iThis, ++iOther) {
      int compare = vertex(iThis).compareTo(other.vertex(iOther));
      if (compare != 0) {
        return compare;
      }
    }
    return 0;
  }

  /**
   * Calculates firstLogicalVertex, the vertex in this loop that comes first in
   * a total ordering of all vertices (by way of S2Point's compareTo function).
   */
  private void initFirstLogicalVertex() {
    int first = 0;
    for (int i = 1; i < numVertices; ++i) {
      if (vertex(i).compareTo(vertex(first)) < 0) {
        first = i;
      }
    }
    firstLogicalVertex = first;
  }

  /**
   * Return true if the loop is generally a left-turning aka counter-clockwise
   * loop.
   */
  public boolean isNormalized() {
    // Optimization: if the longitude span is less than 180 degrees, then the
    // loop covers less than half the sphere and is therefore normalized.
    if (bound.lng().getLength() < S2.M_PI) {
      return true;
    }

    // We allow some error so that hemispheres are always considered normalized.
    // TODO(user): This is no longer required by the S2Polygon implementation,
    // so alternatively we could create the invariant that a loop is normalized
    // if and only if its complement is not normalized.
    return getTurningAngle() >= -1e-14;
  }

  /**
   * Invert the loop if necessary so that the area enclosed by the loop is at
   * most 2*Pi.
   */
  public void normalize() {
    if (!isNormalized()) {
      invert();
    }
  }

  /**
   * Reverse the order of the loop vertices, effectively complementing the
   * region represented by the loop.
   */
  public void invert() {
    int last = numVertices() - 1;
    if (isEmptyOrFull()) {
      vertices[0] = isFull() ? EMPTY_VERTEX : FULL_VERTEX;
    } else {
      for (int i = (last - 1) / 2; i >= 0; --i) {
        S2Point t = vertices[i];
        vertices[i] = vertices[last - i];
        vertices[last - i] = t;
      }
    }
    vertexToIndex = null;
    index = null;
    originInside ^= true;
    if (bound.lat().lo() > -S2.M_PI_2 && bound.lat().hi() < S2.M_PI_2) {
      // The complement of this loop contains both poles.
      bound = S2LatLngRect.full();
      subregionBound = bound;
    } else {
      initBound();
    }
    initFirstLogicalVertex();
  }

  /**
   * Helper method to get area and optionally centroid.
   */
  private S2AreaCentroid getAreaCentroid(boolean doCentroid) {
    S2Point centroid = null;
    if (isEmptyOrFull()) {
      // Return steradian area of the full or empty loop, based on whether the origin is inside.
      return new S2AreaCentroid(originInside ? (4 * S2.M_PI) : 0D, centroid);
    } else if (numVertices() < 3) {
      // Don't crash even if loop is not well-defined.
      return new S2AreaCentroid(0D, centroid);
    }

    // The triangle area calculation becomes numerically unstable as the length
    // of any edge approaches 180 degrees. However, a loop may contain vertices
    // that are 180 degrees apart and still be valid, e.g. a loop that defines
    // the northern hemisphere using four points. We handle this case by using
    // triangles centered around an origin that is slightly displaced from the
    // first vertex. The amount of displacement is enough to get plenty of
    // accuracy for antipodal points, but small enough so that we still get
    // accurate areas for very tiny triangles.
    //
    // Of course, if the loop contains a point that is exactly antipodal from
    // our slightly displaced vertex, the area will still be unstable, but we
    // expect this case to be very unlikely (i.e. a polygon with two vertices on
    // opposite sides of the Earth with one of them displaced by about 2mm in
    // exactly the right direction). Note that the approximate point resolution
    // using the E7 or S2CellId representation is only about 1cm.

    S2Point origin = vertex(0);
    int axis = (origin.largestAbsComponent() + 1) % 3;
    double slightlyDisplaced = origin.get(axis) + S2.M_E * 1e-10;
    origin =
        new S2Point((axis == 0) ? slightlyDisplaced : origin.x,
            (axis == 1) ? slightlyDisplaced : origin.y, (axis == 2) ? slightlyDisplaced : origin.z);
    origin = S2Point.normalize(origin);

    double areaSum = 0;
    S2Point centroidSum = S2Point.ORIGIN;
    for (int i = 1; i <= numVertices(); ++i) {
      areaSum += S2.signedArea(origin, vertex(i - 1), vertex(i));
      if (doCentroid) {
        // The true centroid is already premultiplied by the triangle area.
        S2Point trueCentroid = S2.trueCentroid(origin, vertex(i - 1), vertex(i));
        centroidSum = S2Point.add(centroidSum, trueCentroid);
      }
    }
    // The calculated area at this point should be between -4*Pi and 4*Pi,
    // although it may be slightly larger or smaller than this due to
    // numerical errors.
    // assert (Math.abs(areaSum) <= 4 * S2.M_PI + 1e-12);

    if (areaSum < 0) {
      // If the area is negative, we have computed the area to the right of the
      // loop. The area to the left is 4*Pi - (-area). Amazingly, the centroid
      // does not need to be changed, since it is the negative of the integral
      // of position over the region to the right of the loop. This is the same
      // as the integral of position over the region to the left of the loop,
      // since the integral of position over the entire sphere is (0, 0, 0).
      areaSum += 4 * S2.M_PI;
    }
    // The loop's sign() does not affect the return result and should be taken
    // into account by the caller.
    if (doCentroid) {
      centroid = centroidSum;
    }
    return new S2AreaCentroid(areaSum, centroid);
  }

  /**
   * Return the area of the loop interior, i.e. the region on the left side of
   * the loop. The return value is between 0 and 4*Pi and the true centroid of
   * the loop multiplied by the area of the loop (see S2.java for details on
   * centroids). Note that the centroid may not be contained by the loop.
   */
  public S2AreaCentroid getAreaAndCentroid() {
    return getAreaCentroid(true);
  }

  /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops. The return value is between 0 and 4*Pi.
   */
  public double getArea() {
    return getAreaCentroid(false).getArea();
  }

  /**
   * Return the true centroid of the polygon multiplied by the area of the
   * polygon (see {@link S2} for details on centroids). Note that the centroid
   * may not be contained by the polygon.
   */
  public S2Point getCentroid() {
    return getAreaCentroid(true).getCentroid();
  }

  public double getTurningAngle() {
    // For empty and full loops, we return the limit value as the loop area approaches 0 or 4*Pi
    // respectively.
    if (isEmptyOrFull()) {
      return originInside ? (-2 * S2.M_PI) : (2 * S2.M_PI);
    }

    // Don't crash even if the loop is not well-defined.
    if (numVertices < 3) {
      return 0;
    }

    // To ensure that we get the same result when the loop vertex order is
    // rotated, and that we get the same result with the opposite sign when the
    // vertices are reversed, we need to be careful to add up the individual
    // turn angles in a consistent order.  In general, adding up a set of
    // numbers in a different order can change the sum due to rounding errors.
    int i = firstLogicalVertex;
    int n = numVertices;
    int dir;
    if (vertex(i + 1).compareTo(vertex(i + n - 1)) < 0) {
      dir = 1;
      // 0 <= first <= n-1, so (first+n*dir) <= 2*n-1.
    } else {
      dir = -1;
      i += n;
      // n <= first <= 2*n-1, so (first+n*dir) >= 0.
    }
    double angle = S2.turnAngle(vertex((i + n - dir) % n), vertex(i), vertex((i + dir) % n));
    while (--n > 0) {
      i += dir;
      angle += S2.turnAngle(vertex(i - dir), vertex(i), vertex(i + dir));
    }
    return dir * angle;
  }

  // The following are the possible relationships between two loops A and B:
  //
  // (1) A and B do not intersect.
  // (2) A contains B.
  // (3) B contains A.
  // (4) The boundaries of A and B cross (i.e. the boundary of A
  // intersects the interior and exterior of B and vice versa).
  // (5) (A union B) is the entire sphere (i.e. A contains the
  // complement of B and vice versa).
  //
  // More than one of these may be true at the same time, for example if
  // A == B or A == Complement(B).

  /**
   * Return true if the region contained by this loop is a superset of the
   * region contained by the given other loop.
   */
  public boolean contains(S2Loop b) {
    // For this loop A to contains the given loop B, all of the following must
    // be true:
    //
    // (1) There are no edge crossings between A and B except at vertices.
    //
    // (2) At every vertex that is shared between A and B, the local edge
    // ordering implies that A contains B.
    //
    // (3) If there are no shared vertices, then A must contain a vertex of B
    // and B must not contain a vertex of A. (An arbitrary vertex may be
    // chosen in each case.)
    //
    // The second part of (3) is necessary to detect the case of two loops whose
    // union is the entire sphere, i.e. two loops that contains each other's
    // boundaries but not each other's interiors.

    if (!subregionBound.contains(b.getRectBound())) {
      return false;
    }

    // Special cases to handle either loop being empty or full.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return isFull() || b.isEmpty();
    }

    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B. Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick rejection test.
    if (!contains(b.vertex(0)) && findVertex(b.vertex(0)) < 0) {
      return false;
    }

    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    if (checkEdgeCrossings(b, new S2EdgeUtil.WedgeContains()) <= 0) {
      return false;
    }

    // At this point we know that the boundaries of A and B do not intersect,
    // and that A contains a vertex of B. However we still need to check for
    // the case mentioned above, where (A union B) is the entire sphere.
    // Normally this check is very cheap due to the bounding box precondition.
    if (bound.union(b.getRectBound()).isFull()) {
      if (b.contains(vertex(0)) && b.findVertex(vertex(0)) < 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if the region contained by this loop intersects the region
   * contained by the given other loop.
   */
  public boolean intersects(S2Loop b) {
    // a.intersects(b) if and only if !a.complement().contains(b).
    if (isFull() && !b.isEmpty()) {
      return true;
    }

    // This code is similar to contains(), but is optimized for the case
    // where both loops enclose less than half of the sphere.

    if (!bound.intersects(b.getRectBound())) {
      return false;
    }

    // Normalize the arguments so that B has a smaller longitude span than A.
    // This makes intersection tests much more efficient in the case where
    // longitude pruning is used (see CheckEdgeCrossings).
    if (b.getRectBound().lng().getLength() > bound.lng().getLength()) {
      return b.intersects(this);
    }

    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B. Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick acceptance test.
    if (!b.isEmpty() && contains(b.vertex(0)) && findVertex(b.vertex(0)) < 0) {
      return true;
    }

    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    if (checkEdgeCrossings(b, new S2EdgeUtil.WedgeIntersects()) < 0) {
      return true;
    }

    // We know that A does not contain a vertex of B, and that there are no edge
    // crossings. Therefore the only way that A can intersect B is if B
    // entirely contains A. We can check this by testing whether B contains an
    // arbitrary non-shared vertex of A. Note that this check is cheap because
    // of the bounding box precondition and the fact that we normalized the
    // arguments so that A's longitude span is at least as long as B's.
    if (b.subregionBound.contains(bound)) {
      if (!isEmpty() && b.contains(vertex(0)) && b.findVertex(vertex(0)) < 0) {
        return true;
      }
    }

    return false;
  }

  /**
   * Returns true if the wedge (a0, ab1, a2) contains the edge (ab1, b2), where [a0, ab1, a2] are
   * a subest of the vertices of loop A, and [ab1, ab2, b2] are a subset of the vertices of loop B.
   *
   * <p>Shared edges are handled as follows: If XY is a shared edge, define reversed(XY) to be true
   * if this edge appears in opposite directions in A and B. Then A contains XY if and only if
   * reversed(XY) == "bReverse".
   */
  private static boolean wedgeContainsEdge(
      S2Point a0, S2Point ab1, S2Point a2, S2Point b2, boolean bReverse) {
    if (b2.equalsPoint(a0) || b2.equalsPoint(a2)) {
      // We have a shared or reversed edge.
      return (b2.equalsPoint(a0)) == bReverse;
    } else {
      return S2.orderedCCW(a0, a2, b2, ab1);
    }
  }

  /**
   * A {@link S2EdgeUtil.WedgeProcessor} for {@link #compareBoundary(S2Loop)}.  If resultKnown()
   * returns true, then compareBoundary() returns +1 if loop A contains loop B, -1 if their
   * boundaries cross each other, and 0 otherwise.
   */
  private static final class CompareBoundaryProcessor implements S2EdgeUtil.WedgeProcessor {
    /** True if loop B should be reversed. */
    private final boolean bReverse;
    /** True if any wedge was processed. */
    private boolean foundWedge = false;
    /** True if any edge of B is contained by A. */
    private boolean containsEdge = false;
    /** True if any edge of B is excluded by A. */
    private boolean excludesEdge = false;

    public CompareBoundaryProcessor(boolean bReverse) {
      this.bReverse = bReverse;
    }

    /**
     * Returns true if wedge processing alone was able to determine the compareBoundary() result.
     */
    public boolean resultKnown() {
      return foundWedge;
    }

    /**
     * Returns +1 if A contains the boundary of B, -1 if A excludes the boundary of B, and 0 if the
     * boundaries of A and B cross.
     */
    public int compareBoundary() {
      return containsEdge ? (excludesEdge ? 0 : 1) : -1;
    }

    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // Because we don't care about the interior of B, only its boundary, it is sufficient to check
      // whether A contains the edge (ab1, b2).
      foundWedge = true;
      if (wedgeContainsEdge(a0, ab1, a2, b2, bReverse)) {
        containsEdge = true;
      } else {
        excludesEdge = true;
      }
      return containsEdge & excludesEdge ? -1 : 0;
    }
  }

/**
   * Given two loops of a polygon, return true if A contains B. This version of Contains() is cheap
   * because it does not test for edge intersections. The loops must meet all the S2Polygon
   * requirements; for example this implies that their boundaries may not cross or have any shared
   * edges (although they may have shared vertices).
   */
  public boolean containsNested(S2Loop b) {
    if (!subregionBound.contains(b.getRectBound())) {
      return false;
    }

    // Special cases to handle either loop being empty or full.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return isFull() || b.isEmpty();
    }

    // We are given that A and B do not share any edges, and that either one
    // loop contains the other or they do not intersect.
    int m = findVertex(b.vertex(1));
    if (m < 0) {
      // Since b->vertex(1) is not shared, we can check whether A contains it.
      return contains(b.vertex(1));
    }
    // Check whether the edge order around b->vertex(1) is compatible with
    // A containin B.
    return (new S2EdgeUtil.WedgeContains()).test(
        vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(0), b.vertex(2)) > 0;
  }

  /**
   * Returns +1 if A contains the boundary of B, -1 if A excludes the bound of B, and 0 if the
   * boundaries of A and B cross.
   *
   * <p>Shared edges are handled as follows: If XY is a shared edge, define reversed(XY) to be true
   * if XY appears in opposite directions in A and B. Then A contains XY if and only if
   * reversed(XY) == B->isHole(). Intuitively, this checks whether A contains a vanishingly small
   * region extending from the boundary of B toward the interior of the polygon to which loop B
   * belongs.
   *
   * <p>This method is used for testing containment and intersection of multi-loop polygons. Note
   * that this method is not symmetric, since the result depends on the direction of loop A but not
   * on the direction of loop B (in the absence of shared edges).
   *
   * @param b the loop to compare against this loop; neither loop may be empty, and if {@code b} is
   * full, then it must not be a hole.
   */
  public int compareBoundary(S2Loop b) {
    Preconditions.checkArgument(!isEmpty() && !b.isEmpty());
    Preconditions.checkArgument(!b.isFull() || !b.isHole());

    // The bounds must intersect for containment or crossing.
    if (!bound.intersects(b.bound)) {
      return -1;
    }

    // Full loops are handled as though the loop surrounded the entire sphere.
    if (isFull()) {
      return 1;
    }
    if (b.isFull()) {
      return -1;
    }

    // Check whether there are any edge crossings, and also check the loop relationship at any
    // shared vertices.
    CompareBoundaryProcessor wedge = new CompareBoundaryProcessor(b.isHole());
    if (checkEdgeCrossings(b, wedge) == -1) {
      return 0;
    }
    if (wedge.resultKnown()) {
      return wedge.compareBoundary();
    }

    // There are no edge intersections or shared vertices, so we can check whether A contains an
    // arbitrary vertex of B.
    return contains(b.vertex(0)) ? 1 : -1;
  }

  /**
   * Given two loops whose boundaries do not cross (see {@link #compareBoundary(S2Loop)}, returns
   * true if A contains the boundary of B.
   *
   * <p>This method is cheaper than compareBoundary() because it does not test for edge
   * intersections.
   *
   * @param b the loop to test for containment by this loop; neither loop may be empty or have an
   * edge crossing with the other, and if b is full then bReverse must be false.
   * @param bReverse  If true, the boundary of B is reversed first (which only affects the result
   * when there are shared edges).
   */
  public boolean containsNonCrossingBoundary(S2Loop b, boolean bReverse) {
    assert (!isEmpty() && !b.isEmpty());
    assert (!b.isFull() || !bReverse);

    // The bounds must intersect for containment.
    if (!bound.intersects(b.bound)) {
      return false;
    }

    // Full loops are handled as though the loop surrounded the entire sphere.
    if (isFull()) {
      return true;
    }
    if (b.isFull()) {
      return false;
    }

    int m = findVertex(b.vertex(0));
    if (m < 0) {
      // Since vertex b0 is not shared, we can check whether A contains it.
      return contains(b.vertex(0));
    }

    // Otherwise check whether the edge (b0, b1) is contained by A.
    return wedgeContainsEdge(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(1), bReverse);
  }

  /**
   * Return +1 if A contains B (i.e. the interior of B is a subset of the
   * interior of A), -1 if the boundaries of A and B cross, and 0 otherwise.
   * Requires that A does not properly contain the complement of B, i.e. A and B
   * do not contain each other's boundaries. This method is used for testing
   * whether multi-loop polygons contain each other.
   */
  public int containsOrCrosses(S2Loop b) {
    // There can be containment or crossing only if the bounds intersect.
    if (!bound.intersects(b.getRectBound())) {
      return 0;
    }

    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices. Note that unlike Contains() or
    // Intersects(), we can't do a point containment test as a shortcut because
    // we need to detect whether there are any edge crossings.
    int result = checkEdgeCrossings(b, new S2EdgeUtil.WedgeContainsOrCrosses());

    // If there was an edge crossing or a shared vertex, we know the result
    // already. (This is true even if the result is 1, but since we don't
    // bother keeping track of whether a shared vertex was seen, we handle this
    // case below.)
    if (result <= 0) {
      return result;
    }

    // At this point we know that the boundaries do not intersect, and we are
    // given that (A union B) is a proper subset of the sphere. Furthermore
    // either A contains B, or there are no shared vertices (due to the check
    // above). So now we just need to distinguish the case where A contains B
    // from the case where B contains A or the two loops are disjoint.
    if (!subregionBound.contains(b.getRectBound())) {
      return 0;
    }
    if (!contains(b.vertex(0)) && findVertex(b.vertex(0)) < 0) {
      return 0;
    }

    return 1;
  }

  /**
   * Return true if two loops have the same boundary. This is true if and only if the loops have the
   * same vertices in the same cyclic order. The empty and full loops are considered to have
   * different boundaries. (For testing purposes.)
   */
  boolean boundaryEquals(S2Loop b) {
    if (numVertices != b.numVertices) {
      return false;
    }

    // Special case to handle empty or full loops.  Since they have the same number of vertices, if
    // one loop is empty/full then so is the other.
    if (isEmptyOrFull()) {
      return isEmpty() == b.isEmpty();
    }

    for (int offset = 0; offset < numVertices; ++offset) {
      if (vertex(offset).equalsPoint(b.vertex(0))) {
        // There is at most one starting offset since loop vertices are unique.
        for (int i = 0; i < numVertices; ++i) {
          if (!vertex(i + offset).equalsPoint(b.vertex(i))) {
            return false;
          }
        }
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if two loops have the same boundary except for vertex
   * perturbations. More precisely, the vertices in the two loops must be in the
   * same cyclic order, and corresponding vertex pairs must be separated by no
   * more than maxError. Note: This method mostly useful only for testing
   * purposes.
   */
  boolean boundaryApproxEquals(S2Loop b, double maxError) {
    final S2Loop a = this;
    if (a.numVertices() != b.numVertices()) {
      return false;
    }

    // Special case to handle empty or full loops.  Since they have the same
    // number of vertices, if one loop is empty/full then so is the other.
    if (isEmptyOrFull()) {
      return isEmpty() == b.isEmpty();
    }

    for (int offset = 0; offset < a.numVertices(); ++offset) {
      if (S2.approxEquals(a.vertex(offset), b.vertex(0), maxError)) {
        boolean success = true;
        for (int i = 0; i < a.numVertices(); ++i) {
          if (!S2.approxEquals(a.vertex(i + offset), b.vertex(i), maxError)) {
            success = false;
            break;
          }
        }
        if (success) {
          return true;
        }
        // Otherwise continue looping.  There may be more than one candidate
        // starting offset since vertices are only matched approximately.
      }
    }
    return false;
  }

  boolean boundaryApproxEquals(S2Loop loop) {
    return boundaryApproxEquals(loop, 1e-15);
  }

  /**
   * Offsets into two loops at which a boundary distance comparison will start.
   *
   * <p>Used only by {@code matchBoundaries}.
   */
  private static final class LoopOffsets {
    /** The offset of the first loop. */
    public final int first;
    /** The offset of the second loop. */
    public final int second;
    public LoopOffsets(int first, int second) {
      this.first = first;
      this.second = second;
    }
    @Override
    public int hashCode() {
      return first * 517 + second;
    }
    @Override
    public boolean equals(Object o) {
      if (o instanceof LoopOffsets) {
        LoopOffsets that = (LoopOffsets) o;
        return this.first == that.first && this.second == that.second;
      } else {
        return false;
      }
    }
  }

  /**
   * Helper method called by {@code boundaryNear()} to determine if this loop
   * and loop {@code b} remain within {@code maxError} of each other, starting
   * the comparison with this loop at vertex {@code a_offset} and loop
   * {@code b} at vertex 0.
   */
  boolean matchBoundaries(S2Loop b, int a_offset, double maxError) {
    final S2Loop a = this;

    // The state consists of a pair (i,j).  A state transition consists of
    // incrementing either "i" or "j".  "i" can be incremented only if
    // a(i+1+a_offset) is near the edge from b(j) to b(j+1), and a similar rule
    // applies to "j".  The function returns true iff we can proceed all the way
    // around both loops in this way.
    //
    // Note that when "i" and "j" can both be incremented, sometimes only one
    // choice leads to a solution.  We handle this using a stack and
    // backtracking.  We also keep track of which states have already been
    // explored to avoid duplicating work.

    List<LoopOffsets> pending = Lists.newArrayList();
    Multiset<LoopOffsets> done = HashMultiset.create();
    pending.add(new LoopOffsets(0, 0));
    while (!pending.isEmpty()) {
      LoopOffsets last = pending.remove(pending.size() - 1);
      int i = last.first;
      int j = last.second;
      if (i == a.numVertices() && j == b.numVertices()) {
        return true;
      }
      done.add(new LoopOffsets(i, j));

      // If (i == na && offset == na-1) where na == a.numVertices(), then
      // (i+1+offset) overflows the [0, 2*na-1] range allowed by vertex().
      // So we reduce the range if necessary.
      int io = i + a_offset;
      if (io >= a.numVertices()) {
        io -= a.numVertices();
      }

      if (i < a.numVertices() && done.count(new LoopOffsets(i + 1, j)) == 0 &&
          S2EdgeUtil.getDistance(a.vertex(io + 1),
              b.vertex(j),
              b.vertex(j + 1)).radians() <= maxError) {
        pending.add(new LoopOffsets(i + 1, j));
      }
      if (j < b.numVertices() && done.count(new LoopOffsets(i, j + 1)) == 0 &&
          S2EdgeUtil.getDistance(b.vertex(j + 1),
              a.vertex(io),
              a.vertex(io + 1)).radians() <= maxError) {
        pending.add(new LoopOffsets(i, j + 1));
      }
    }
    return false;
  }

  /**
   * Return true if the two loop boundaries are within "max_error" of each other
   * along their entire lengths. The two loops may have different numbers of
   * vertices. More precisely, this method returns true if the two loops have
   * parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2 such that
   * {@code distance(a(t), b(t)) <= max_error} for all t.
   *
   * <p>You can think of this as testing whether it is possible to drive two
   * cars all the way around the two loops such that no car ever goes backward
   * and the cars are always within "max_error" of each other.
   *
   * <p>(Package private, only used for testing purposes.)
   */
  boolean boundaryNear(S2Loop b, double max_error) {
    // Special case to handle empty or full loops.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return (isEmpty() && b.isEmpty()) || (isFull() && b.isFull());
    }

    for (int a_offset = 0; a_offset < numVertices(); ++a_offset) {
      if (matchBoundaries(b, a_offset, max_error)) {
        return true;
      }
    }
    return false;
  }

  boolean boundaryNear(S2Loop loop) {
    return boundaryNear(loop, 1e-15);
  }

  // S2Region interface (see {@code S2Region} for details):

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
   * Returns a slightly looser bounding latitude-longitude rectangle than that returned by
   * {@link #getRectBound()}. It is not guaranteed that if this loop contains a loop X, then the
   * subregion bound will contain X.getRectBound().
   */
  public S2LatLngRect getSubregionBound() {
    return subregionBound;
  }

  /**
   * If this method returns true, the region completely contains the given cell.
   * Otherwise, either the region does not contain the cell or the containment
   * relationship could not be determined.
   */
  @Override
  public boolean contains(S2Cell cell) {
    // It is faster to construct a bounding rectangle for an S2Cell than for
    // a general polygon.
    S2LatLngRect cellBound = cell.getRectBound();
    // We can't check subregionBound.contains(cell.getRectBound()) because S2Cell bounds are not
    // calculated using S2EdgeUtil::RectBounder.
    if (!bound.contains(cellBound)) {
      return false;
    }
    return contains(new S2Loop(cell));
  }

  /**
   * If this method returns false, the region does not intersect the given cell.
   * Otherwise, either region intersects the cell, or the intersection
   * relationship could not be determined.
   */
  @Override
  public boolean mayIntersect(S2Cell cell) {
    // Check the cell bound first, since it is faster to construct a bounding rectangle for an
    // S2Cell than for a general polygon.
    S2LatLngRect cellBound = cell.getRectBound();
    if (!bound.intersects(cellBound)) {
      return false;
    }
    return new S2Loop(cell).intersects(this);
  }

  /**
   * The point 'p' does not need to be normalized.
   */
  public boolean contains(S2Point p) {
    if (!bound.contains(p)) {
      return false;
    }

    boolean inside = originInside;
    S2Point origin = S2.origin();
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(origin, p,
        vertices[numVertices - 1]);

    // The s2edgeindex library is not optimized yet for long edges,
    // so the tradeoff to using it comes with larger loops.
    if (numVertices < 2000) {
      // For the full or empty loop, we only call edgeOrVertexCrossing once,
      // with the same vertex used to initialize it above, so the result is
      // always 'originInside'.
      for (int i = 0; i < numVertices; i++) {
        inside ^= crosser.edgeOrVertexCrossing(vertices[i]);
      }
    } else {
      DataEdgeIterator it = getEdgeIterator(numVertices);
      int previousIndex = -2;
      for (it.getCandidates(origin, p); it.hasNext(); it.next()) {
        int ai = it.index();
        if (previousIndex != ai - 1) {
          crosser.restartAt(vertices[ai]);
        }
        previousIndex = ai;
        inside ^= crosser.edgeOrVertexCrossing(vertex(ai + 1));
      }
    }

    return inside;
  }

  /**
   * Returns the shortest distance from a point P to this loop, given as the
   * angle formed between P, the origin and the nearest point on the loop to P.
   * This angle in radians is equivalent to the arclength along the unit sphere.
   */
  public S1Angle getDistance(S2Point p) {
    S2Point normalized = S2Point.normalize(p);

    // The furthest point from p on the sphere is its antipode, which is an
    // angle of PI radians. This is an upper bound on the angle.
    S1Angle minDistance = S1Angle.radians(Math.PI);
    for (int i = 0; i < numVertices(); i++) {
      minDistance =
          S1Angle.min(minDistance, S2EdgeUtil.getDistance(normalized, vertex(i), vertex(i + 1)));
    }
    return minDistance;
  }

  /**
   * Creates an edge index over the vertices, which by itself takes no time.
   * Then the expected number of queries is used to determine whether brute
   * force lookups are likely to be slower than really creating an index, and if
   * so, we do so. Finally an iterator is returned that can be used to perform
   * edge lookups.
   */
  private final DataEdgeIterator getEdgeIterator(int expectedQueries) {
    if (index == null) {
      index = new S2EdgeIndex() {
        @Override
        public int getNumEdges() {
          return numVertices;
        }

        @Override
        public S2Point edgeFrom(int index) {
          return vertex(index);
        }

        @Override
        public S2Point edgeTo(int index) {
          return vertex(index + 1);
        }
      };
    }
    index.predictAdditionalCalls(expectedQueries);
    return new S2EdgeIndex.DataEdgeIterator(index);
  }

  /**
   * Return true if the S2:origin() is inside this loop.
   *
   * Primarily used to serialize internal details about a loop for later fast initialization.
   */
  public boolean isOriginInside() {
    return originInside;
  }

  /** Return true if this loop is valid. */
  public boolean isValid() {
    // The subregionBound_ must be at least as large as bound.
    if (!subregionBound.contains(bound)) {
      log.info("Subregion bound not initialized correctly");
      return false;
    }

    // Loops must have at least 3 vertices (except for "empty" and "full").
    if (numVertices < 3) {
      if (isEmptyOrFull()) {
        return true;
      }
      log.info("Non-empty, non-full loops must have at least 3 vertices.");
      return false;
    }

    // All vertices must be unit length.
    for (int i = 0; i < numVertices; ++i) {
      if (!S2.isUnitLength(vertex(i))) {
        log.info("Vertex " + i + " is not unit length");
        return false;
      }
    }

    // Loops are not allowed to have any duplicate vertices.
    HashMap<S2Point, Integer> vmap = Maps.newHashMap();
    for (int i = 0; i < numVertices; ++i) {
      Integer previousVertexIndex = vmap.put(vertex(i), i);
      if (previousVertexIndex != null) {
        log.info("Duplicate vertices: " + previousVertexIndex + " and " + i);
        return false;
      }
    }

    // Non-adjacent edges are not allowed to intersect.
    boolean crosses = false;
    DataEdgeIterator it = getEdgeIterator(numVertices);
    for (int a1 = 0; a1 < numVertices; a1++) {
      int a2 = (a1 + 1) % numVertices;
      EdgeCrosser crosser = new EdgeCrosser(vertex(a1), vertex(a2), vertex(0));
      int previousIndex = -2;
      for (it.getCandidates(vertex(a1), vertex(a2)); it.hasNext(); it.next()) {
        int b1 = it.index();
        int b2 = (b1 + 1) % numVertices;
        // If either 'a' index equals either 'b' index, then these two edges
        // share a vertex. If a1==b1 then it must be the case that a2==b2, e.g.
        // the two edges are the same. In that case, we skip the test, since we
        // don't want to test an edge against itself. If a1==b2 or b1==a2 then
        // we have one edge ending at the start of the other, or in other words,
        // the edges share a vertex -- and in S2 space, where edges are always
        // great circle segments on a sphere, edges can only intersect at most
        // once, so we don't need to do further checks in that case either.
        if (a1 != b2 && a2 != b1 && a1 != b1) {
          // WORKAROUND(shakusa, ericv): S2.robustCCW() currently
          // requires arbitrary-precision arithmetic to be truly robust. That
          // means it can give the wrong answers in cases where we are trying
          // to determine edge intersections. The workaround is to ignore
          // intersections between edge pairs where all four points are
          // nearly colinear.
          double abc = S2.angle(vertex(a1), vertex(a2), vertex(b1));
          boolean abcNearlyLinear = S2.approxEquals(abc, 0D, MAX_INTERSECTION_ERROR) ||
              S2.approxEquals(abc, S2.M_PI, MAX_INTERSECTION_ERROR);
          double abd = S2.angle(vertex(a1), vertex(a2), vertex(b2));
          boolean abdNearlyLinear = S2.approxEquals(abd, 0D, MAX_INTERSECTION_ERROR) ||
              S2.approxEquals(abd, S2.M_PI, MAX_INTERSECTION_ERROR);
          if (abcNearlyLinear && abdNearlyLinear) {
            continue;
          }

          if (previousIndex != b1) {
            crosser.restartAt(vertex(b1));
          }

          // Beware, this may return the loop is valid if there is a
          // "vertex crossing".
          // TODO(user): Fix that.
          crosses = crosser.robustCrossing(vertex(b2)) > 0;
          previousIndex = b2;
          if (crosses ) {
            log.info("Edges " + a1 + " and " + b1 + " cross");
            log.info("Edge locations in degrees: " +
                new S2LatLng(vertex(a1)).toStringDegrees() +
                "-" + new S2LatLng(vertex(a2)).toStringDegrees() +
                " and " + new S2LatLng(vertex(b1)).toStringDegrees() +
                "-" + new S2LatLng(vertex(b2)).toStringDegrees());
            return false;
          }
        }
      }
    }

    return true;
  }

  /**
   * Static version of isValid(), to be used only when an S2Loop instance is not
   * available, but validity of the points must be checked.
   *
   * @return true if the given loop is valid. Creates an instance of S2Loop and
   *         defers this call to {@link #isValid()}.
   */
  public static boolean isValid(List<S2Point> vertices) {
    return new S2Loop(vertices).isValid();
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder("S2Loop, ");

    builder.append(vertices.length).append(" points. [");

    for (S2Point v : vertices) {
      builder.append(v.toDegreesString()).append(" ");
    }
    builder.append("]");

    return builder.toString();
  }

  private void initOriginAndBound() {
    if (numVertices < 3) {
      // Check for the special "empty" and "full" loops (which have one vertex).
      if (isEmptyOrFull()) {
        // If the vertex is in the southern hemisphere then the loop is full,
        // otherwise it is empty.
        originInside = vertex(0).getZ() < 0;
      } else {
        originInside = false;
      }
      // Make sure we have non-null bounds whether it's valid or not.
      if (originInside) {
        bound = S2LatLngRect.full();
        subregionBound = S2LatLngRect.full();
      } else {
        bound = S2LatLngRect.empty();
        subregionBound = S2LatLngRect.empty();
      }
      return;
    }

    // To ensure that every point is contained in exactly one face of a
    // subdivision of the sphere, all containment tests are done by counting the
    // edge crossings starting at a fixed point on the sphere (S2::Origin()).
    // We need to know whether this point is inside or outside of the loop.
    // We do this by first guessing that it is outside, and then seeing whether
    // we get the correct containment result for vertex 1. If the result is
    // incorrect, the origin must be inside the loop.
    //
    // A loop with consecutive vertices A,B,C contains vertex B if and only if
    // the fixed vector R = S2.ortho(B) is contained by the wedge ABC.
    // The test below is written so that B is inside if C=R but not if A=R.
    // (Note that we can't use S2.origin() as the fixed vector because of the
    // possibility that B == S2.origin().)
    originInside = false; // Initialize before calling contains().
    boolean v1Inside = S2.orderedCCW(S2.ortho(vertex(1)), vertex(0), vertex(2), vertex(1));

    // Note that we need to initialize bound_ with a temporary value since contains() does a
    // bounding rectangle check before doing anything else.
    bound = S2LatLngRect.full();

    if (v1Inside != contains(vertex(1))) {
      originInside = true;
    }

    initBound();
  }

  private void initBound() {
    // Check for the special "empty" and "full" loops.
    if (isEmptyOrFull()) {
      if (isEmpty()) {
        subregionBound = bound = S2LatLngRect.empty();
      } else {
        subregionBound = bound = S2LatLngRect.full();
      }
      return;
    }

    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices.  First, the maximal latitude may be
    // attained along the interior of an edge.  Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe).  Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    S2EdgeUtil.RectBounder bounder = new S2EdgeUtil.RectBounder();
    for (int i = 0; i <= numVertices(); ++i) {
      bounder.addPoint(vertex(i));
    }
    S2LatLngRect b = bounder.getBound();

    // Note that we need to initialize bound with a temporary value since
    // contains() does a bounding rectangle check before doing anything else.
    bound = S2LatLngRect.full();
    if (contains(S2Point.Z_POS)) {
      b = new S2LatLngRect(new R1Interval(b.lat().lo(), S2.M_PI_2), S1Interval.full());
    }

    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.lng().isFull() due to the test above.
    if (b.lng().isFull() && contains(S2Point.Z_NEG)) {
      b = new S2LatLngRect(new R1Interval(-S2.M_PI_2, b.lat().hi()), b.lng());
    }

    bound = b;
    subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(bound);
  }

  /**
   * Return the index of a vertex at point "p", or -1 if not found. The return
   * value is in the range 1..num_vertices_ if found.
   */
  private int findVertex(S2Point p) {
    if (vertexToIndex == null) {
      vertexToIndex = new HashMap<S2Point, Integer>();
      for (int i = 1; i <= numVertices; i++) {
        vertexToIndex.put(vertex(i), i);
      }
    }
    Integer index = vertexToIndex.get(p);
    if (index == null) {
      return -1;
    } else {
      return index;
    }
  }

  /**
   * This method encapsulates the common code for loop containment and
   * intersection tests. It is used in three slightly different variations to
   * implement contains(), intersects(), and containsOrCrosses().
   *
   *  In a nutshell, this method checks all the edges of this loop (A) for
   * intersection with all the edges of B. It returns -1 immediately if any edge
   * intersections are found. Otherwise, if there are any shared vertices, it
   * returns the minimum value of the given WedgeRelation for all such vertices
   * (returning immediately if any wedge returns -1). Returns +1 if there are no
   * intersections and no shared vertices.
   */
  private int checkEdgeCrossings(S2Loop b, S2EdgeUtil.WedgeProcessor relation) {
    DataEdgeIterator it = getEdgeIterator(b.numVertices);
    int result = 1;
    // since 'this' usually has many more vertices than 'b', use the index on
    // 'this' and loop over 'b'
    for (int j = 0; j < b.numVertices(); ++j) {
      S2EdgeUtil.EdgeCrosser crosser =
        new S2EdgeUtil.EdgeCrosser(b.vertex(j), b.vertex(j + 1), vertex(0));
      int previousIndex = -2;
      for (it.getCandidates(b.vertex(j), b.vertex(j + 1)); it.hasNext(); it.next()) {
        int i = it.index();
        if (previousIndex != i - 1) {
          crosser.restartAt(vertex(i));
        }
        previousIndex = i;
        int crossing = crosser.robustCrossing(vertex(i + 1));
        if (crossing < 0) {
          continue;
        }
        if (crossing > 0) {
          return -1; // There is a proper edge crossing.
        }
        if (vertex(i + 1).equalsPoint(b.vertex(j + 1))) {
          result = Math.min(result, relation.test(
              vertex(i), vertex(i + 1), vertex(i + 2), b.vertex(j), b.vertex(j + 2)));
          if (result < 0) {
            return result;
          }
        }
      }
    }
    return result;
  }
}
