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

import static com.google.common.geometry.S2.M_PI_2;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.S2ShapeUtil.RangeIterator;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.IOException;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicInteger;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;
import org.jspecify.annotations.Nullable;

/**
 * An S2Loop represents a simple spherical polygon. It consists of a single chain of vertices where
 * the first vertex is implicitly connected to the last. All loops are defined to have a CCW
 * orientation, i.e. the interior of the loop is on the left side of the edges. This implies that a
 * clockwise loop enclosing a small area is interpreted to be a CCW loop enclosing a very large
 * area.
 *
 * <p>Loops are not allowed to have any duplicate vertices (whether adjacent or not), and non-
 * adjacent edges are not allowed to intersect. Loops must have at least 3 vertices, except that
 * every single-vertex loop is considered to be either "full" or "empty", as discussed below.
 * Loop vertices must be valid and unit length. Although these restrictions are not enforced in
 * optimized code, you may get unexpected results if they are violated. Assertion errors may be
 * thrown if assertions are enabled.
 *
 * <p>S2Loop can also represent two special cases: The "empty" loop contains no points, while the
 * "full" loop contains all points. These loops do not have any edges, but to preserve the invariant
 * that every loop can be represented as a vertex chain, they are defined as single vertex loops.
 * Any single vertex loop represents either the empty or full loop. See {@link #isEmpty()} and
 * {@link #isFull()}.)
 *
 * <p>Note that {@link S2Polygon} has stricter conventions for full and empty loops. It only
 * considers single-vertex loops to be full or empty if the single vertex in the loop is equal to
 * the ones in the "canonical" empty or full loops, as returned by {@link #empty()} and
 * {@link #full()}. All other single-vertex loops are considered invalid to S2Polygon. That is
 * safer, as "accidental" single-vertex loops may arise from snapping or other operations.
 *
 * <p>Point containment of loops is defined such that if the sphere is subdivided into faces
 * (loops), every point is contained by exactly one face. This implies that loops do not necessarily
 * contain their vertices, and merely "touching" loops, with one or more shared vertices or edges,
 * are not considered intersecting.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
@SuppressWarnings("Assertion")
public final class S2Loop implements S2Region, Comparable<S2Loop>, Serializable, S2Shape {

  @VisibleForTesting static final byte LOSSLESS_ENCODING_VERSION = 1;

  /** Max angle that intersections can be off by and yet still be considered collinear. */
  public static final double MAX_INTERSECTION_ERROR = 1e-15;

  /**
   * A canonical vertex for defining the empty loop that contains no area. However, note that *any*
   * loop with a single vertex having a positive Z coordinate is considered to be the empty loop.
   */
  static final S2Point EMPTY_VERTEX = S2Point.Z_POS;

  /**
   * A canonical vertex for defining the full loop that contains the whole sphere. However, note
   * that *any* loop with a single vertex having a negative Z coordinate is considered to be the
   * full loop.
   */
  static final S2Point FULL_VERTEX = S2Point.Z_NEG;

  /** Spatial index for this loop. */
  @VisibleForTesting transient S2ShapeIndex index;

  /**
   * In general we build the index the first time it is needed, but we make an exception for
   * contains(S2Point) because this method has a simple brute force implementation that is
   * relatively cheap. For this one method we keep track of the number of calls made and only build
   * the index once enough calls have been made that we think an index would be worthwhile.
   */
  private final AtomicInteger unindexedContainsCalls = new AtomicInteger();

  private final S2Point[] vertices;
  private final int numVertices;

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
   * Constructs a loop with the given vertices. The last vertex is implicitly connected to the
   * first. All points must be unit length. Loops must have at least 3 vertices, except that single
   * vertex loops are allowed, because any loop with a single vertex is considered to be either
   * "empty" or "full": see {@link #empty()} and {@link #full()}.
   *
   * @param vertices the vertices for this new loop
   */
  @JsIgnore
  public S2Loop(final List<S2Point> vertices) {
    initIndex();
    this.numVertices = vertices.size();
    this.vertices = new S2Point[numVertices];
    vertices.toArray(this.vertices);
    this.depth = 0;
    initOriginAndBound();
  }

  /**
   * Similar to above, but takes ownership of the provided vertices without copying. Package
   * private, as the caller is not allowed to modify the array after this.
   */
  @JsIgnore
  S2Loop(final S2Point[] vertices) {
    initIndex();
    this.numVertices = vertices.length;
    this.vertices = vertices;
    this.depth = 0;
    initOriginAndBound();
  }

  /**
   * Fast/unsafe loop initialization.
   *
   * <p>This constructor accepts known good values for bounds and the originInside value. This is
   * intended to be a "fast loop creation" when we already know a lot about the loop. It is
   * primarily used in combination with the fast S2Polygon initializer ( {@link
   * S2Polygon#initWithNestedLoops(java.util.Map)}). The last vertex is implicitly connected to the
   * first. All points should be unit length. Loops must have at least 3 vertices, except for the
   * empty and full loops (see {@link #empty()} and {@link #full()}.)
   *
   * @param vertices loop vertices
   * @param originInside true if the S2.origin() is inside the loop
   * @param bound the lat/long bounds of the loop
   * @return new loop.
   */
  public static S2Loop newLoopWithTrustedDetails(
      List<S2Point> vertices, boolean originInside, S2LatLngRect bound) {
    // This is a static method to try to discourage its use.
    return new S2Loop(vertices, originInside, bound);
  }

  /** Create a circle of points with a given center, radius, and number of vertices. */
  public static S2Loop makeRegularLoop(S2Point center, S1Angle radius, int numVertices) {
    assert S2.isUnitLength(center);
    assert numVertices >= 3;
    return new S2Loop(makeRegularVertices(center, radius, numVertices));
  }

  /**
   * Like the function above, creates a circle of points with a given radius, and number of
   * vertices, but this version constructs a loop centered around the z-axis of the given coordinate
   * frame, with the first vertex in the direction of the positive x-axis. (This allows the loop to
   * be rotated for testing purposes.)
   */
  @JsIgnore
  public static S2Loop makeRegularLoop(Matrix frame, S1Angle radius, int numVertices) {
    return new S2Loop(makeRegularVertices(frame, radius, numVertices));
  }

  /**
   * Like {@link #makeRegularLoop(S2Point, S1Angle, int)}, creates a circle of points with a given
   * center, radius, and number of vertices, but returns the vertices as a list.
   */
  public static List<S2Point> makeRegularVertices(S2Point center, S1Angle radius, int numVertices) {
    return makeRegularVertices(S2.getFrame(center), radius, numVertices);
  }

  /**
   * Like the function above, but constructs the list of vertices around the z-axis of the given
   * coordinate frame, with the first vertex in the direction of the positive x-axis.
   */
  @JsIgnore
  public static List<S2Point> makeRegularVertices(Matrix frame, S1Angle radius, int numVertices) {
    List<S2Point> vertices = Lists.newArrayList();
    // We construct the loop in the given frame coordinates, with the center at (0, 0, 1). For a
    // loop of radius "r", the loop vertices have the form (x, y, z) where x^2 + y^2 = sin(r) and
    // z = cos(r). The distance on the sphere (arc length) from each vertex to the center is
    // acos(cos(r)) = r.
    double z = cos(radius.radians());
    double r = sin(radius.radians());
    double radianStep = 2 * PI / numVertices;
    for (int i = 0; i < numVertices; ++i) {
      double angle = i * radianStep;
      S2Point p = new S2Point(r * cos(angle), r * sin(angle), z);
      vertices.add(S2.fromFrame(frame, p).normalize());
    }
    return vertices;
  }

  /** Constructs an S2Loop, making a copy of the provided vertices in a new array. */
  private S2Loop(List<S2Point> vertices, boolean originInside, S2LatLngRect bound) {
    initIndex();
    this.numVertices = vertices.size();
    this.vertices = new S2Point[numVertices];
    this.bound = bound;
    this.subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(bound);
    this.depth = 0;
    this.originInside = originInside;

    vertices.toArray(this.vertices);
  }

  /** Constructs a loop corresponding to the given cell. */
  @JsIgnore
  public S2Loop(S2Cell cell) {
    initIndex();
    numVertices = 4;
    vertices = new S2Point[numVertices];
    depth = 0;
    for (int i = 0; i < 4; ++i) {
      vertices[i] = cell.getVertex(i);
    }
    // We compute the bounding rectangle ourselves, since S2Cell uses a different method and we need
    // all the bounds to be consistent.
    initOriginAndBound();
  }

  /** Copy constructor. */
  @JsIgnore
  public S2Loop(S2Loop src) {
    initIndex();
    this.numVertices = src.numVertices;
    this.vertices = new S2Point[numVertices];
    for (int i = 0; i < numVertices; i++) {
      this.vertices[i] = src.vertices[i];
    }
    this.bound = src.getRectBound();
    this.subregionBound = src.subregionBound;
    this.originInside = src.originInside;
    this.depth = src.depth();
  }

  private void initIndex() {
    // 'maxUnindexedContainsCalls' was tuned using the benchmarks. We wait until the cumulative time
    // we would have saved with an index approximately equals the cost of building the index, and
    // then build it. (This gives the optimal competitive ratio of 2; look up 'competitive
    // algorithms' for details.)  We err on the side of building the index too early, because
    // building the index may be forced anyways by other API calls. We select
    // 'maxUnindexedContainsCalls' based on the number of vertices since there is great variation in
    // contains() efficiency.
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

    index = new S2ShapeIndex();
    index.add(this);
  }

  /**
   * Returns the same instance after initializing transient fields. Called by Java deserialization.
   */
  @CanIgnoreReturnValue
  private Object readResolve() {
    initIndex();
    return this;
  }

  /**
   * Returns a new canonical single-vertex loop that defines an empty loop (i.e., a loop with no
   * edges that contains no points.)
   */
  public static final S2Loop empty() {
    return new S2Loop(Collections.singletonList(EMPTY_VERTEX));
  }

  /**
   * Returns a new canonical single-vertex loop that defines a full loop (i.e., a loop with no edges
   * that contains all points).
   */
  public static final S2Loop full() {
    return new S2Loop(Collections.singletonList(FULL_VERTEX));
  }

  /**
   * Note that this doesn't do anything smart: it just compares a few fields for equality, but
   * doesn't check for equivalent loops that were initialized in different ways, etc. For instance,
   * two loops with the same vertices in the same order but starting at different positions are not
   * considered equal here.
   *
   * <p>Also, any single-vertex loop defines either the empty or full loop, but this method
   * does not consider two "full" loops to be equal unless they use the same vertex, and likewise
   * for "empty" loops. TODO(torrey): That should be fixed.
   */
  @Override
  public boolean equals(Object obj) {
    if (obj instanceof S2Loop) {
      S2Loop that = (S2Loop) obj;
      return Arrays.equals(this.vertices, that.vertices)
          && (this.originInside == that.originInside)
          && Objects.equals(this.bound, that.bound);
    }
    return false;
  }

  /**
   * For performance, this implementation of hashCode does not use the loop vertices, but only
   * depends on the loop bounds, number of vertices, and originInside values.
   */
  @Override
  public int hashCode() {
    return ((bound.hashCode() * 7) + numVertices * 11) + (originInside ? 1 : 0);
  }

  public int depth() {
    return depth;
  }

  /**
   * The depth of a loop is defined as its nesting level within its containing polygon. "Outer
   * shell" loops have depth 0, holes within those loops have depth 1, shells within those holes
   * have depth 2, etc. This field is only used by the S2Polygon implementation.
   */
  public void setDepth(int depth) {
    this.depth = depth;
  }

  /** Return true if this loop represents a hole in its containing polygon. */
  public boolean isHole() {
    return (depth & 1) != 0;
  }

  /**
   * The sign of a loop is -1 if the loop represents a hole in its containing polygon, and +1
   * otherwise.
   */
  public int sign() {
    return isHole() ? -1 : 1;
  }

  public int numVertices() {
    return numVertices;
  }

  /**
   * For convenience, we make two entire copies of the vertex list available: vertex(n..2*n-1) is
   * mapped to vertex(0..n-1), where n == numVertices().
   */
  public S2Point vertex(int i) {
    try {
      return vertices[i >= vertices.length ? i - vertices.length : i];
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new IllegalStateException(
          "Invalid vertex index " + i + " for loop of " + vertices.length + " vertices.", e);
    }
  }

  /**
   * Like vertex(), but this method returns vertices in reverse order if the loop represents a
   * polygon hole. For example, arguments 0, 1, 2 are mapped to vertices n-1, n-2, n-3, where n ==
   * numVertices(). This ensures that the interior of the polygon is always to the left of the
   * vertex chain.
   */
  public S2Point orientedVertex(int i) {
    if (isHole()) {
      i = 2 * numVertices() - 1 - i;
    }
    return vertex(i);
  }

  /** Returns an unmodifiable view of the vertices of this polyline. */
  public List<S2Point> vertices() {
    return Collections.unmodifiableList(Arrays.asList(vertices));
  }

  /** Returns the vertices oriented such that left is on the inside. */
  public List<S2Point> orientedVertices() {
    return new AbstractList<S2Point>() {
      @Override
      public int size() {
        return numVertices;
      }

      @Override
      public S2Point get(int index) {
        return orientedVertex(index);
      }
    };
  }

  /**
   * Returns true if this is an "empty" loop that contains no points. Note that any single-vertex
   * loop is considered to be either empty or full.
   */
  @Override
  public boolean isEmpty() {
    return isEmptyOrFull() && !originInside;
  }

  /**
   * Returns true if this is a "full" loop that contains all points.  Note that any single-vertex
   * loop is considered to be either empty or full.
   */
  @Override
  public boolean isFull() {
    return isEmptyOrFull() && originInside;
  }

  /**
   * Returns true if this loop has a single vertex, and thus is defined to be either the "empty" or
   * "full" loop.
   */
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
    return true;
  }

  @Override
  public boolean containsOrigin() {
    return originInside;
  }

  @Override
  public int numChains() {
    return isEmpty() ? 0 : 1;
  }

  @Override
  public int getChainStart(int chainId) {
    Preconditions.checkElementIndex(chainId, numChains());
    return 0;
  }

  @Override
  public int getChainLength(int chainId) {
    Preconditions.checkElementIndex(chainId, numChains());
    return numEdges();
  }

  @Override
  public void getChainEdge(int chainId, int offset, MutableEdge result) {
    Preconditions.checkElementIndex(chainId, numChains());
    getEdge(offset, result);
  }

  @Override
  public S2Point getChainVertex(int chainId, int edgeOffset) {
    Preconditions.checkElementIndex(chainId, numChains());
    return vertex(edgeOffset);
  }

  @Override
  public void getChainPosition(int edgeId, ChainPosition result) {
    Preconditions.checkElementIndex(edgeId, numEdges());
    // All the edges are in the single chain.
    result.set(0, edgeId);
  }

  @Override
  public int dimension() {
    return 2;
  }

  /** Comparator (needed by Comparable interface) */
  @Override
  public int compareTo(S2Loop other) {
    if (numVertices != other.numVertices) {
      return numVertices - other.numVertices;
    } else if (numVertices == 0) {
      return 0;
    }
    // Compare the two loops' vertices, starting with each loop's firstLogicalVertex. This allows us
    // to always catch cases where logically identical loops have different vertex orderings (e.g.
    // ABCD and BCDA).
    int maxVertices = numVertices;
    int iThis = getCanonicalFirstVertex() % numVertices;
    int iOther = other.getCanonicalFirstVertex() % other.numVertices;
    for (int i = 0; i < maxVertices; ++i, ++iThis, ++iOther) {
      int compare = vertex(iThis).compareTo(other.vertex(iOther));
      if (compare != 0) {
        return compare;
      }
    }
    return 0;
  }

  /**
   * Returns a canonical minimum vertex such that the vertex sequence starting at this vertex does
   * not change when the loop vertex order is rotated or inverted. This allows the loop vertices to
   * be traversed in a canonical order. If the initial value is less than {@link #numVertices()},
   * then stable iteration should be toward larger indices, otherwise smaller indices.
   */
  private int getCanonicalFirstVertex() {
    int first = 0;
    for (int i = 1; i < numVertices; ++i) {
      if (vertex(i).compareTo(vertex(first)) < 0) {
        first = i;
      }
    }
    if (numVertices > 0 && vertex(first + 1).compareTo(vertex(first + numVertices - 1)) >= 0) {
      first += numVertices;
    }
    return first;
  }

  /** Return true if the loop is generally a left-turning, aka counter-clockwise loop. */
  public boolean isNormalized() {
    // Optimization: if the longitude span is less than 180 degrees, then the loop covers less than
    // half the sphere and is therefore normalized.
    if (bound.lng().getLength() < PI) {
      return true;
    }

    // We allow some error so that hemispheres are always considered normalized.
    // TODO(user): This is no longer required by the S2Polygon implementation, so
    // alternatively we could create the invariant that a loop is normalized if and only if its
    // complement is not normalized.
    return S2ShapeMeasures.turningAngle(vertices()) >= -S2.getTurningAngleMaxError(numVertices);
  }

  /** Invert the loop if necessary so that the area enclosed by the loop is at most 2*Pi. */
  @CanIgnoreReturnValue
  public S2Loop normalize() {
    if (!isNormalized()) {
      invert();
    }
    return this;
  }

  /**
   * Reverse the order of the loop vertices, effectively complementing the region represented by the
   * loop.
   */
  public void invert() {
    initIndex();
    int last = numVertices - 1;
    if (isEmptyOrFull()) {
      vertices[0] = isFull() ? EMPTY_VERTEX : FULL_VERTEX;
    } else {
      for (int i = (last - 1) / 2; i >= 0; --i) {
        S2Point t = vertices[i];
        vertices[i] = vertices[last - i];
        vertices[last - i] = t;
      }
    }
    originInside ^= true;
    if (bound.lat().lo() > -M_PI_2 && bound.lat().hi() < M_PI_2) {
      // The complement of this loop contains both poles.
      bound = S2LatLngRect.full();
      subregionBound = bound;
    } else {
      // initBound() requires 'bound' to be null in order for its contains(S2Point) calls to work.
      bound = null;
      initBound();
    }
  }

  /**
   * Returns the area of the loop interior, i.e. the region on the left side of the loop regardless
   * of whether it is a shell or a hole. This is measured in steradians, and the value is between 0
   * and 4*Pi, or explicitly 0 if the loop has fewer than 3 vertices.
   */
  public double getArea() {
    if (isEmptyOrFull()) {
      return originInside ? (4 * PI) : 0;
    } else if (numVertices < 3) {
      return 0;
    }
    return S2ShapeMeasures.area(this);
  }

  /**
   * Returns the true centroid of the loop multiplied by the area of the loop, or null if this loop
   * is empty, full, or has fewer than 3 vertices.
   *
   * <p>The result is not unit length, so you may want to normalize it. Also note that in general,
   * the centroid may not be contained by the loop. See {@link S2} for additional centroid details.
   *
   * <p>We prescale by the loop area for two reasons:
   *
   * <ol>
   *   <li>It is cheaper to compute this way, and
   *   <li>It makes it easier to compute the centroid of more complicated shapes (by splitting them
   *       into disjoint regions and summing their centroids).
   * </ol>
   *
   * <p>Note that the return value is not affected by whether this loop is a "hole" or a "shell".
   */
  public @Nullable S2Point getCentroid() {
    if (numVertices < 3) {
      return null;
    }
    return S2ShapeMeasures.centroid(this);
  }

  /**
   * Returns a pair of {@link #getArea()} and {@link #getCentroid()}. Note that the centroid is not
   * unit length: it is scaled by the area of the loop.
   */
  public S2AreaCentroid getAreaAndCentroid() {
    return new S2AreaCentroid(getArea(), getCentroid());
  }

  /**
   * Returns the sum of the turning angles at each vertex. The return value is positive if the loop
   * is counter-clockwise, negative if the loop is clockwise, and zero if the loop is a great
   * circle.
   *
   * <p>Degenerate and nearly-degenerate loops are handled consistently with {@link
   * S2Predicates#sign(S2Point, S2Point, S2Point)}.
   *
   * <p>For example, if a loop has zero area (i.e., it is a very small CCW loop) then the turning
   * angle will always be negative.
   *
   * <p>This quantity is also called the "geodesic curvature" of the loop.
   */
  public double getTurningAngle() {
    // For empty and full loops, we return the limit value as the loop area approaches 0 or 4*Pi
    // respectively.
    if (isEmptyOrFull()) {
      return originInside ? (-2 * PI) : (2 * PI);
    }
    return S2ShapeMeasures.turningAngle(vertices());
  }

  // The following are the possible relationships between two loops A and B:
  //
  // (1) A and B do not intersect.
  // (2) A contains B.
  // (3) B contains A.
  // (4) The boundaries of A and B cross (i.e. the boundary of A intersects the interior and
  // exterior of B and vice versa).
  // (5) (A union B) is the entire sphere (i.e. A contains the complement of B and vice versa).
  //
  // More than one of these may be true at the same time, for example if A == B or
  // A == Complement(B).

  /**
   * Return true if the region contained by this loop is a superset of the region contained by the
   * given other loop.
   */
  public boolean contains(S2Loop b) {
    // For this loop A to contain the given loop B, all of the following must be true:
    //
    // (1) There are no edge crossings between A and B except at vertices.
    //
    // (2) At every vertex that is shared between A and B, the local edge ordering implies that A
    // contains B.
    //
    // (3) If there are no shared vertices, then A must contain a vertex of B and B must not contain
    // a vertex of A. (An arbitrary vertex may be chosen in each case.)
    //
    // The second part of (3) is necessary to detect the case of two loops whose union is the entire
    // sphere, i.e. two loops that contain each other's boundaries but not each other's interiors.

    if (!subregionBound.contains(b.bound)) {
      return false;
    }

    // Special cases to handle either loop being empty or full.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return isFull() || b.isEmpty();
    }

    // Check whether there are any edge crossings, and also check the loop relationship at any
    // shared vertices.
    ContainsRelation relation = new ContainsRelation();
    if (hasCrossingRelation(this, b, relation)) {
      return false;
    }

    // There are no crossings, and if there are any shared vertices, then A contains B locally at
    // each shared vertex.
    if (relation.foundSharedVertex()) {
      return true;
    }

    // Since there are no edge intersections or shared vertices, we just need to test condition (3)
    // above. We can skip this test if we discovered that A contains at least one point of B while
    // checking for edge crossings.
    if (!contains(b.vertex(0))) {
      return false;
    }

    // We still need to check whether (A union B) is the entire sphere. Normally this check is very
    // cheap due to the bounding box precondition.
    if ((b.subregionBound.contains(bound) || b.bound.union(bound).isFull())
        && b.contains(vertex(0))) {
      return false;
    }

    return true;
  }

  /**
   * Return true if the region contained by this loop "A" intersects the region contained by the
   * given other loop "B". Note that loops that only "touch" at a common point or along a shared
   * (but reversed) edge are not considered intersecting.
   */
  public boolean intersects(S2Loop b) {
    // This code is similar to contains(), but is optimized for the case where both loops enclose
    // less than half of the sphere.
    if (!bound.intersects(b.bound)) {
      return false;
    }

    // Now check whether there are any edge crossings, and also check the loop relationship at any
    // shared vertices.
    IntersectsRelation relation = new IntersectsRelation();
    if (hasCrossingRelation(this, b, relation)) {
      return true;
    }
    if (relation.foundSharedVertex()) {
      return false;
    }

    // Since there are no edge intersections or shared vertices, the loops intersect only if this
    // loop A contains B, B contains A, or the two loops contain each other's boundaries. These
    // checks are usually cheap because of the bounding box preconditions. Note that neither loop is
    // empty (because of the bounding box check above), so it is safe to access vertex(0).

    // Check whether this loop A contains B, or A and B contain each other's boundaries. (Note that
    // A contains all the vertices of B in either case.)
    if (subregionBound.contains(b.bound) || bound.union(b.bound).isFull()) {
      if (contains(b.vertex(0))) {
        return true;
      }
    }
    // Check whether B contains this loop A.
    if (b.subregionBound.contains(bound)) {
      if (b.contains(vertex(0))) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if the wedge (a0, ab1, a2) contains the edge (ab1, b2), where [a0, ab1, a2] are a
   * subset of the vertices of loop A, and [ab1, ab2, b2] are a subset of the vertices of loop B.
   *
   * <p>Shared edges are handled as follows: If XY is a shared edge, define reversed(XY) to be true
   * if this edge appears in opposite directions in A and B. Then A contains XY if and only if
   * {@code reversed(XY) == bReversed}.
   */
  static boolean wedgeContainsSemiwedge(
      S2Point a0, S2Point ab1, S2Point a2, S2Point b2, boolean bReversed) {
    boolean b2EqualsA0 = b2.equalsPoint(a0);
    if (b2EqualsA0 || b2.equalsPoint(a2)) {
      // We have a shared or reversed edge.
      return b2EqualsA0 == bReversed;
    } else {
      return S2Predicates.orderedCCW(a0, a2, b2, ab1);
    }
  }

  /**
   * Given another loop "B" whose boundaries do not cross the boundaries of this loop "A", (see
   * {@link #compareBoundary(S2Loop)}, returns true if this loop A contains the boundary of B. If
   * "bReverse" is true, the boundary of B is reversed first (which only affects the result when
   * there are shared edges). This method is cheaper than compareBoundary() because it does not test
   * for edge intersections, but it does trigger building the index if there are edges in both
   * polygons and the bounds intersect.
   *
   * <p>Requires that neither loop is empty or has an edge crossing with the other, and if B is
   * full, then bReverse must be false.
   */
  boolean containsNonCrossingBoundary(S2Loop b, boolean bReverse) {
    assert !isEmpty();
    assert !b.isEmpty();
    assert !b.isFull() || !bReverse;

    // The bounds must intersect for containment.
    if (!getRectBound().intersects(b.getRectBound())) {
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
      // Since b.vertex(0) is not shared, we can check whether A contains it.
      return contains(b.vertex(0));
    }
    // Otherwise check whether the edge (b0, b1) is contained by A.
    return wedgeContainsSemiwedge(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(1), bReverse);
  }

  /**
   * Given that this loop "A" and the given loop "B" are both loops of a common polygon, return true
   * if this loop A contains B. This version of contains() is cheap because it does not test for
   * edge intersections. The loops must meet all the S2Polygon requirements; for example this
   * implies that their boundaries may not cross or have any shared edges (although they may have
   * shared vertices).
   */
  public boolean containsNested(S2Loop b) {
    if (!subregionBound.contains(b.bound)) {
      return false;
    }

    // Special cases to handle either loop being empty or full. We check if b.numVertices() < 2
    // in order to also handle the case where b.numVertices() == 0.
    if (isEmptyOrFull() || b.numVertices() < 2) {
      return isFull() || b.isEmpty();
    }

    // We are given that this loop A and B do not share any edges, and that either one loop contains
    // the other or they do not intersect.
    int m = findVertex(b.vertex(1));
    if (m < 0) {
      // Since b.vertex(1) is not shared, we can check whether A contains it.
      return contains(b.vertex(1));
    }
    // Check whether the edge order around b.vertex(1) is compatible with this loop A containing B.
    return new S2EdgeUtil.WedgeContains()
            .test(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(0), b.vertex(2))
        > 0;
  }

  /**
   * Returns +1 if this loop A contains the boundary of B, -1 if A excludes the boundary of B, and 0
   * if the boundaries of A and B cross.
   *
   * <p>Shared edges are handled as follows: If XY is a shared edge, define reversed(XY) to be true
   * if XY appears in opposite directions in A and B. Then A contains XY if and only if {@code
   * reversed(XY) == B.isHole()}. Intuitively, this checks whether A contains a vanishingly small
   * region extending from the boundary of B toward the interior of the polygon to which loop B
   * belongs.
   *
   * <p>This method is used for testing containment and intersection of multi-loop polygons. Note
   * that this method is not symmetric, since the result depends on the direction of this loop A but
   * not on the direction of loop B (in the absence of shared edges).
   *
   * @param b the loop to compare against this loop; neither loop may be empty, and if {@code b} is
   *     full, then it must not be a hole.
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
    CompareBoundaryRelation relation = new CompareBoundaryRelation(b.isHole());
    if (hasCrossingRelation(this, b, relation)) {
      return 0;
    }
    if (relation.foundSharedVertex()) {
      return relation.containsEdge() ? 1 : -1;
    }

    // There are no edge intersections or shared vertices, so we can check whether this loop A
    // contains an arbitrary vertex of B.
    return contains(b.vertex(0)) ? 1 : -1;
  }

  /**
   * Returns true if this loop and the given loop 'b' have the same boundary. This is true if and
   * only if the loops have the same vertices in the same cyclic order. The empty and full loops are
   * considered to have different boundaries. (For testing purposes.)
   */
  boolean boundaryEquals(S2Loop b) {
    if (numVertices != b.numVertices) {
      return false;
    }

    // Special case to handle empty or full loops. Since they have the same number of vertices, if
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
   * Returns true if this loop and the given loop 'b' have the same boundary except for vertex
   * perturbations. More precisely, the vertices in the two loops must be in the same cyclic order,
   * and corresponding vertex pairs must be separated by no more than maxErrorRadians. Note: This
   * method is mostly useful for testing purposes.
   */
  boolean boundaryApproxEquals(S2Loop b, double maxErrorRadians) {
    final S2Loop a = this;
    if (a.numVertices != b.numVertices) {
      return false;
    }

    // Special case to handle empty or full loops. Since they have the same number of vertices, if
    // one loop is empty/full then so is the other.
    if (isEmptyOrFull()) {
      return isEmpty() == b.isEmpty();
    }

    for (int offset = 0; offset < a.numVertices; ++offset) {
      if (S2.approxEquals(a.vertex(offset), b.vertex(0), maxErrorRadians)) {
        boolean success = true;
        for (int i = 0; i < a.numVertices; ++i) {
          if (!S2.approxEquals(a.vertex(i + offset), b.vertex(i), maxErrorRadians)) {
            success = false;
            break;
          }
        }
        if (success) {
          return true;
        }
        // Otherwise continue looping. There may be more than one candidate starting offset since
        // vertices are only matched approximately.
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
   * Helper method called by {@code boundaryNear()} to determine if this loop and loop {@code b}
   * remain within {@code maxErrorRadians} of each other, starting the comparison with this loop at
   * vertex {@code aOffset} and loop {@code b} at vertex 0.
   */
  boolean matchBoundaries(S2Loop b, int aOffset, double maxErrorRadians) {
    final S2Loop a = this;

    // The state consists of a pair (i,j). A state transition consists of incrementing either "i"
    // or "j". "i" can be incremented only if a(i+1+aOffset) is near the edge from b(j) to b(j+1),
    // and a similar rule applies to "j". The function returns true iff we can proceed all the way
    // around both loops in this way.
    //
    // Note that when "i" and "j" can both be incremented, sometimes only one choice leads to a
    // solution. We handle this using a stack and backtracking. We also keep track of which states
    // have already been explored to avoid duplicating work.

    List<LoopOffsets> pending = Lists.newArrayList();
    // TODO(torrey): This doesn't need to be a multiset. The value is only ever checked against 0.
    Multiset<LoopOffsets> done = HashMultiset.create();
    pending.add(new LoopOffsets(0, 0));
    while (!pending.isEmpty()) {
      LoopOffsets last = pending.remove(pending.size() - 1);
      int i = last.first;
      int j = last.second;
      if (i == a.numVertices && j == b.numVertices) {
        return true;
      }
      done.add(new LoopOffsets(i, j));

      // If (i == na && offset == na-1) where na == a.numVertices, then (i+1+offset) overflows the
      // [0, 2*na-1] range allowed by vertex(). So we reduce the range if necessary.
      int io = i + aOffset;
      if (io >= a.numVertices) {
        io -= a.numVertices;
      }

      if (i < a.numVertices
          && done.count(new LoopOffsets(i + 1, j)) == 0
          && S2EdgeUtil.getDistance(a.vertex(io + 1), b.vertex(j), b.vertex(j + 1)).radians()
              <= maxErrorRadians) {
        pending.add(new LoopOffsets(i + 1, j));
      }
      if (j < b.numVertices
          && done.count(new LoopOffsets(i, j + 1)) == 0
          && S2EdgeUtil.getDistance(b.vertex(j + 1), a.vertex(io), a.vertex(io + 1)).radians()
              <= maxErrorRadians) {
        pending.add(new LoopOffsets(i, j + 1));
      }
    }
    return false;
  }

  /**
   * Returns true if the this loop and the given loop 'b' have boundaries are within {@code
   * maxError} of each other along their entire lengths. The two loops may have different numbers of
   * vertices. More precisely, this method returns true if the two loops have parameterizations
   * a:[0,1] -> S^2, b:[0,1] -> S^2 such that {@code distance(a(t), b(t)) <= maxError} for all t.
   *
   * <p>You can think of this as testing whether it is possible to drive two cars all the way around
   * the two loops such that no car ever goes backward and the cars are always within {@code
   * maxError} of each other.
   *
   * <p>(Package private, only used for testing purposes.)
   */
  boolean boundaryNear(S2Loop b, double maxError) {
    // Special case to handle empty or full loops.
    if (isEmptyOrFull() || b.isEmptyOrFull()) {
      return (isEmpty() && b.isEmpty()) || (isFull() && b.isFull());
    }

    for (int aOffset = 0; aOffset < numVertices; ++aOffset) {
      if (matchBoundaries(b, aOffset, maxError)) {
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
   * also. However, it is not guaranteed that if this loop contains a loop X, then this loops' bound
   * contains X.getRectBound(): see {@link #getSubregionBound()} if that is required.
   */
  @Override
  public S2LatLngRect getRectBound() {
    return bound;
  }

  /**
   * Returns a slightly looser bounding latitude-longitude rectangle than that returned by {@link
   * #getRectBound()}. It is guaranteed that if this loop contains a loop X, then this loops'
   * subregionBound contains X.getRectBound().
   */
  public S2LatLngRect getSubregionBound() {
    return subregionBound;
  }

  /**
   * If this method returns true, the region completely contains the given cell. Otherwise, either
   * the region does not contain the cell or the containment relationship could not be determined.
   */
  @Override
  @JsIgnore
  public boolean contains(S2Cell target) {
    S2Iterator<S2ShapeIndex.Cell> it = index.iterator();
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If 'target' is disjoint from all index cells, it is not contained. Similarly, if 'target' is
    // subdivided into one or more index cells, then it is not contained, since index cells are
    // subdivided only if they (nearly) intersect a sufficient number of edges. (But note that if
    // 'target' itself is an index cell then it may be contained, since it could be a cell with no
    // edges in the loop interior.
    if (relation != S2ShapeIndex.CellRelation.INDEXED) {
      return false;
    }

    // Otherwise check if any edges intersect 'target'.
    if (boundaryApproxIntersects(it, target)) {
      return false;
    }

    // Otherwise check if the loop contains the center of 'target'.
    return contains(it, target.getCenter());
  }

  /**
   * If this method returns false, this region does not intersect the given cell. Otherwise, either
   * this region intersects the cell, or the intersection relationship could not be determined.
   */
  @Override
  public boolean mayIntersect(S2Cell target) {
    S2Iterator<S2ShapeIndex.Cell> it = index.iterator();
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If 'target' does not overlap any index cell, there is no intersection.
    if (relation == S2ShapeIndex.CellRelation.DISJOINT) {
      return false;
    }

    // If 'target' is subdivided into one or more index cells, there is an intersection to within
    // the S2ShapeIndex error bound (see contains).
    if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) {
      return true;
    }

    // If 'target' is an index cell, there is an intersection because index cells are created only
    // if they have at least one edge or they are entirely contained by that loop.
    if (it.compareTo(target.id()) == 0) {
      return true;
    }

    // Otherwise check if any edges intersect 'target'.
    if (boundaryApproxIntersects(it, target)) {
      return true;
    }

    // Otherwise check if the loop contains the center of 'target'.
    return contains(it, target.getCenter());
  }

  /**
   * Returns true if the loop boundary intersects 'target'. It may also return true when the loop
   * boundary does not intersect 'target' but some edge comes within the worst-case error tolerance.
   *
   * <p>Requires: it.id().contains(target.id()). (This condition is true whenever it.locate(target)
   * returns INDEXED.)
   */
  private boolean boundaryApproxIntersects(S2Iterator<S2ShapeIndex.Cell> it, S2Cell target) {
    assert it.id().contains(target.id());
    S2ClippedShape aClipped = it.entry().clipped(0);
    int aNumClipped = aClipped.numEdges();

    // If there are no edges, there is no intersection.
    if (aNumClipped == 0) {
      return false;
    }

    // We can save some work if 'target' is the index cell itself.
    if (it.compareTo(target.id()) == 0) {
      return true;
    }

    // Otherwise check whether any of the edges intersect 'target'.
    R2Rect bound = target.getBoundUV().expanded(S2EdgeUtil.MAX_CELL_EDGE_ERROR);
    R2Vector v0 = new R2Vector();
    R2Vector v1 = new R2Vector();
    for (int i = 0; i < aNumClipped; ++i) {
      int ai = aClipped.edge(i);
      if (S2EdgeUtil.clipToPaddedFace(
              vertex(ai), vertex(ai + 1), target.face(), S2EdgeUtil.MAX_CELL_EDGE_ERROR, v0, v1)
          && S2EdgeUtil.intersectsRect(v0, v1, bound)) {
        return true;
      }
    }

    return false;
  }

  /**
   * Returns a simplified loop, which may be self-intersecting, or null if the entire loop was
   * within the tolerance.
   *
   * <p>If self-intersections could occur and a valid result is needed, instead use {@link
   * S2Polygon#initToSimplified(S2Polygon, S1Angle, boolean)}.
   *
   * <p>Always keeps the first vertex from the loop, and if {@code vertexFilter} is not null, also
   * keeps vertices for which {@code vertexFilter.shouldKeepVertex()} is true.
   */
  @JsIgnore // Predicate is not a JsType.
  public @Nullable S2Loop simplify(S1Angle tolerance, Predicate<S2Point> vertexFilter) {
    if (vertices.length < 2) {
      // Unable to simplify further, just return whatever this loop is.
      return this;
    }

    // Add the last point so we can simplify along the final edge as well, and then remove the last
    // point if it's still there.
    List<S2Point> points = Lists.newArrayListWithCapacity(vertices.length + 1);
    Collections.addAll(points, vertices);
    points.add(vertices[0]);
    S2Polyline line = new S2Polyline(points);
    List<S2Point> simplified = line.subsampleVertices(tolerance).vertices();
    if (simplified.get(0).equalsPoint(Iterables.getLast(simplified))) {
      simplified = simplified.subList(0, simplified.size() - 1);
    }

    // Merge in vertices that we should keep, if a filter was provided.
    if (vertexFilter != null) {
      List<S2Point> toKeep = Lists.newArrayList();
      for (int i = 0, j = 0; i < numVertices(); i++) {
        S2Point p = vertex(i);
        if (simplified.get(j).equalsPoint(p)) {
          j++;
          toKeep.add(p);
        } else if (vertexFilter.apply(p)) {
          toKeep.add(p);
        }
      }
      simplified = toKeep;
    }

    return simplified.size() <= 2 ? null : new S2Loop(simplified);
  }

  /**
   * Returns true if the point is contained by the loop. The containment test is exact, placing
   * {@code p} arbitrarily within or without the loop depending on orientation of the edges, such
   * that given two loops sharing an edge, and a point on that edge, only one of the loops will
   * contain it. The point does not need to be normalized, but does need to be valid.
   */
  @Override
  public boolean contains(S2Point p) {
    // NOTE(ericv): A bounds check slows down the C++ version of this function by about 50%. It is
    // worthwhile only when it might allow us to delay building the index.
    if (!index.isFresh() && bound != null && !bound.contains(p)) {
      return false;
    }

    // We must use the bruteForceContains() during loop initialization because initOriginAndBound()
    // calls contains() before creating the index. Also, for small loops it is faster to just check
    // all the crossings. Otherwise, we keep track of the number of calls to contains() and only
    // build the index when enough calls have been made so that we think it is worth the effort.
    // Note that the code below is structured so that if many calls are made in parallel, only one
    // thread builds the index, while the rest continue using brute force until the index is
    // actually available.
    int maxBruteForceVertices = 32;
    if (numVertices <= maxBruteForceVertices
        || (!index.isFresh() && unindexedContainsCalls.decrementAndGet() > 0)) {
      return bruteForceContains(p);
    }

    // Otherwise we look up the S2ShapeIndex cell containing this point. Note that the index is
    // built automatically the first time an iterator is created.
    S2Iterator<S2ShapeIndex.Cell> it = index.iterator();
    if (!it.locate(p)) {
      return false;
    }
    return contains(it, p);
  }

  // Package-private. Used only by S2Polygon and S2Loop itself.
  boolean bruteForceContains(S2Point p) {
    // Empty and full loops don't need a special case, but invalid loops with zero vertices do, so
    // we might as well handle them all at once.
    if (numVertices < 3) {
      return originInside;
    }

    S2Point origin = S2.origin();
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(origin, p, vertex(0));
    boolean inside = originInside;
    for (int i = 1; i <= numVertices; ++i) {
      inside ^= crosser.edgeOrVertexCrossing(vertex(i));
    }
    return inside;
  }

  /**
   * Given an iterator that is already positioned at the {@link S2ShapeIndex.Cell} containing {@code
   * p}, returns contains(p).
   */
  private boolean contains(S2Iterator<S2ShapeIndex.Cell> it, S2Point p) {
    // Test containment by drawing a line segment from the cell center to the given point and
    // counting edge crossings.
    S2ClippedShape aClipped = it.entry().clipped(0);
    boolean inside = aClipped.containsCenter();
    int aNumClipped = aClipped.numEdges();
    if (aNumClipped > 0) {
      S2Point center = it.center();
      S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(center, p);
      int aiPrev = -2;
      for (int i = 0; i < aNumClipped; ++i) {
        int ai = aClipped.edge(i);
        if (ai != aiPrev + 1) {
          crosser.restartAt(vertex(ai));
        }
        aiPrev = ai;
        inside ^= crosser.edgeOrVertexCrossing(vertex(ai + 1));
      }
    }
    return inside;
  }

  /**
   * Returns the shortest distance from a point P to this loop, given as the angle formed between P,
   * the origin and the nearest point on the loop to P. This angle in radians is equivalent to the
   * arclength along the unit sphere. The point is not required to be normalized.
   */
  public S1Angle getDistance(S2Point p) {
    S2Point normalized = p.normalize();
    // The farthest point from p on the sphere is its antipode, which is an angle of PI radians.
    // This is an upper bound on the angle.
    S1Angle minDistance = S1Angle.radians(PI);
    for (int i = 0; i < numVertices; i++) {
      minDistance =
          S1Angle.min(minDistance, S2EdgeUtil.getDistance(normalized, vertex(i), vertex(i + 1)));
    }
    return minDistance;
  }

  /**
   * Return true if the S2:origin() is inside this loop.
   *
   * <p>Primarily used to serialize internal details about a loop for later fast initialization.
   */
  public boolean isOriginInside() {
    return originInside;
  }

  /**
   * Returns true if this loop is valid. To check validity and also get problem details, use {@link
   * #findValidationError(S2Error)}.
   */
  public boolean isValid() {
    return !findValidationError(new S2Error());
  }

  /**
   * Static version of {@link #isValid()} for checking a set of points before attempting to create a
   * loop from them. To check validity and also get problem details, use {@link
   * #findValidationError(List, S2Error)}.
   *
   * @return true if the given set of vertices forms a valid loop.
   */
  @JsIgnore
  public static boolean isValid(List<S2Point> vertices) {
    return verticesAreUnitLength(vertices) && new S2Loop(vertices).isValid();
  }

  /**
   * Returns true if all vertices are unit length. If Java assertions are enabled, this must be
   * checked first, otherwise a loop cannot even be constructed without triggering an assertion.
   */
  private static boolean verticesAreUnitLength(List<S2Point> vertices) {
    for (int i = 0; i < vertices.size(); ++i) {
      if (!S2.isUnitLength(vertices.get(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if this is *not* a valid loop and sets {@code error} appropriately. Otherwise
   * returns false and leaves {@code error} unchanged. Requires that error != null.
   */
  @CanIgnoreReturnValue
  public boolean findValidationError(S2Error error) {
    return findValidationErrorNoIndex(error)
        || S2ShapeUtil.findSelfIntersection(index, this, error);
  }

  /**
   * Static version of {@link #findValidationError(S2Error)}, for checking if a set of points is
   * valid before attempting to create a loop from them.
   *
   * <p>Returns true if the given list of vertices do *not* form a valid loop and sets {@code error}
   * appropriately. Otherwise returns false and leaves {@code error} unchanged. Requires that error
   * != null.
   */
  @JsIgnore
  public static boolean findValidationError(List<S2Point> vertices, S2Error error) {
    if (!verticesAreUnitLength(vertices)) {
      error.init(S2Error.Code.NOT_UNIT_LENGTH, "Vertex is not unit length.");
      return true;
    }
    return new S2Loop(vertices).findValidationError(error);
  }

  /**
   * Like findValidationError(), but skips any checks that would require building the S2ShapeIndex
   * (i.e., self-intersection tests). This will be used by the S2Polygon implementation, which uses
   * its own index to check for loop self-intersection.
   */
  @CanIgnoreReturnValue
  public boolean findValidationErrorNoIndex(S2Error error) {
    if (!bound.isValid() || !subregionBound.isValid()) {
      error.init(
          S2Error.Code.NOT_UNIT_LENGTH,
          "Invalid loop bound, possibly due to vertices that are not unit length.");
      return true;
    }

    // If the bounds are valid, subregionBound must be at least as large as bound. (This is an
    // internal consistency check rather than a test of client data.)
    assert subregionBound.contains(bound);

    // All vertices must be unit length.
    for (int i = 0; i < numVertices; ++i) {
      if (!S2.isUnitLength(vertex(i))) {
        error.init(S2Error.Code.NOT_UNIT_LENGTH, "Vertex %d is not unit length.", i);
        return true;
      }
    }

    // Valid loops must have at least three vertices, or one vertex. Every single-vertex loop is
    // considered to be either empty or full.
    if (numVertices < 3) {
      if (isEmptyOrFull()) {
        // Skip the remaining tests.
        return false;
      }
      error.init(
          S2Error.Code.LOOP_NOT_ENOUGH_VERTICES,
          "Non-empty, non-full loops must have at least 3 vertices");
      return true;
    }

    // Loops are not allowed to have any degenerate edges (edge with identical vertices).
    for (int i = 0; i < numVertices; ++i) {
      if (vertex(i).equalsPoint(vertex(i + 1))) {
        error.init(
            S2Error.Code.DUPLICATE_VERTICES, "Edge %d is degenerate (duplicate vertex).", i);
        return true;
      }
    }

    return false;
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder("S2Loop, depth=");
    builder.append(depth).append(", ");
    builder.append(vertices.length).append(" points. [");

    for (S2Point v : vertices) {
      builder.append(v.toDegreesString()).append(" ");
    }
    builder.append("]");

    return builder.toString();
  }

  private void initOriginAndBound() {
    if (numVertices < 3) {
      // Check for 'empty' and 'full' loops. Any loop with a single vertex is either empty or full.
      if (isEmptyOrFull()) {
        // If the vertex is in the southern hemisphere then the loop is full, otherwise it is empty.
        originInside = (vertex(0).z < 0);
      } else {
        originInside = false;
      }
    } else {
      // The brute force point containment algorithm works by counting edge crossings starting at a
      // fixed reference point (chosen as S2.origin() for historical reasons). Loop initialization
      // would be more efficient if we used a loop vertex such as vertex(0) as the reference point
      // instead, however making this change would be a lot of work because originInside is
      // currently part of the encoded format.
      //
      // In any case, we initialize originInside by first guessing that it is outside, and then
      // seeing whether we get the correct containment result for vertex 1. If the result is
      // incorrect, the origin must be inside the loop instead. Note that "inside" may not be well
      // defined if the loop is invalid.
      boolean v1Inside = validAngleContainsVertex(vertex(0), vertex(1), vertex(2));

      originInside = false; // Initialize before calling contains().
      if (vertex(1).isValid() && v1Inside != contains(vertex(1))) {
        originInside = true;
      }
    }
    initBound();
  }

  /**
   * Returns true if the angle formed by the three points is valid and contains the vertex.
   *
   * <p>Note that the S2Loop is not necessarily valid and so we need to check the requirements of
   * angleContainsVertex() first.
   */
  private boolean validAngleContainsVertex(S2Point a, S2Point b, S2Point c) {
    return a.isValid()
        && b.isValid()
        && c.isValid()
        && !a.equalsPoint(b)
        && !c.equalsPoint(b)
        && S2Predicates.angleContainsVertex(a, b, c);
  }

  /** Initializes the bound. Requires {@code bound == null}. */
  private void initBound() {
    if (numVertices < 3) {
      if (isFull()) {
        subregionBound = bound = S2LatLngRect.full();
      } else {
        subregionBound = bound = S2LatLngRect.empty();
      }
      return;
    }

    // The bounding rectangle of a loop is not necessarily the same as the bounding rectangle of its
    // vertices. First, the maximal latitude may be attained along the interior of an edge.
    // Second, the loop may wrap entirely around the sphere (e.g. a loop that defines two
    // revolutions of a candy-cane stripe). Third, the loop may include one or both poles. Note
    // that a small clockwise loop near the equator contains both poles.
    S2EdgeUtil.RectBounder bounder = new S2EdgeUtil.RectBounder();
    for (int i = 0; i <= numVertices; ++i) {
      bounder.addPoint(vertex(i));
    }
    S2LatLngRect b = bounder.getBound();

    // Note that this and the following contains calls only work correctly if 'bound' is already
    // null. When constructing the loop, it will be, but invert() has to clear the bound first.
    if (contains(S2Point.Z_POS)) {
      b = new S2LatLngRect(new R1Interval(b.lat().lo(), M_PI_2), S1Interval.full());
    }

    // If a loop contains the south pole, then either it wraps entirely around the sphere (full
    // longitude range), or it also contains the north pole in which case b.lng().isFull() due to
    // the test above.
    if (b.lng().isFull() && contains(S2Point.Z_NEG)) {
      b = new S2LatLngRect(new R1Interval(-M_PI_2, b.lat().hi()), b.lng());
    }

    bound = b;
    subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(bound);
  }

  /**
   * Return the index of a vertex at point "p", or -1 if not found. The return value is in the range
   * 1..numVertices if found.
   */
  @VisibleForTesting
  int findVertex(S2Point p) {
    // 10 is the max number of edges in a cell. As an optimization, it is directly entered here
    // instead of creating an S2ShapeIndex.Options and then calling getMaxEdgesPerCell.
    // TODO(user): Create and use Options.DEFAULT_MAX_EDGES_PER_CELL instead.
    if (numVertices < 10) {
      // Exhaustive search. Return value must be in the range [1..N].
      for (int i = 1; i <= numVertices; ++i) {
        if (vertex(i).equalsPoint(p)) {
          return i;
        }
      }
      return -1;
    }
    S2Iterator<S2ShapeIndex.Cell> it = index.iterator();
    if (!it.locate(p)) {
      return -1;
    }

    S2ClippedShape aClipped = it.entry().clipped(0);
    for (int i = aClipped.numEdges() - 1; i >= 0; --i) {
      int ai = aClipped.edge(i);
      // Return value must be in the range [1..N].
      if (vertex(ai).equalsPoint(p)) {
        return (ai == 0) ? numVertices : ai;
      }
      if (vertex(ai + 1).equalsPoint(p)) {
        return ai + 1;
      }
    }
    return -1;
  }

  private static class CompressedEncodingProperties {
    // Recomputing the bound multiplies the decode time taken per vertex by a factor of about 3.5.
    // Without recomputing the bound, decode takes approximately 125 ns / vertex. A loop with 63
    // vertices encoded without the bound will take ~30us to decode, which is acceptable.
    // At ~3.5 bytes / vertex without the bound, adding the bound will increase the size by <15%,
    // which is also acceptable.
    private static final int MIN_LOOP_VERTICES_FOR_BOUND = 64;

    public enum Property {
      ORIGIN_INSIDE(1L),
      BOUND_ENCODED(1L << 1);

      public final long bitValue;

      private Property(long bitValue) {
        this.bitValue = bitValue;
      }
    }

    private long bits = 0L;

    public CompressedEncodingProperties(S2Loop loop) {
      if (loop.containsOrigin()) {
        setProperty(Property.ORIGIN_INSIDE);
      }

      // Write whether there is a bound so we can change the threshold later.
      if (loop.numVertices() >= MIN_LOOP_VERTICES_FOR_BOUND) {
        setProperty(Property.BOUND_ENCODED);
      }
    }

    public CompressedEncodingProperties(long bits) {
      this.bits = bits;
    }

    public void setProperty(Property property) {
      bits ^= property.bitValue;
    }

    public boolean hasProperty(Property property) {
      return (bits & property.bitValue) != 0;
    }

    public long asLong() {
      return bits;
    }
  }

  void encodeCompressed(int level, LittleEndianOutput encoder) throws IOException {
    // Encode the number of vertices.
    encoder.writeVarint32(numVertices());

    // Encode the individual vertices.
    S2PointCompression.encodePointsCompressed(Arrays.asList(vertices), level, encoder);

    // Encode the compression properties.
    CompressedEncodingProperties properties = new CompressedEncodingProperties(this);
    encoder.writeVarint64(properties.asLong());

    // Encode the depth.
    encoder.writeVarint32(depth());

    // Optionally encode the bounds, if the properties indicate it should be encoded.
    if (properties.hasProperty(CompressedEncodingProperties.Property.BOUND_ENCODED)) {
      getRectBound().encode(encoder);
    }
  }

  static S2Loop decodeCompressed(int level, LittleEndianInput decoder) throws IOException {
    // Decode the number of vertices.
    int numVertices = decoder.readVarint32();
    if (numVertices < 0) {
      throw new IOException(
          "Invalid numVertices: "
              + numVertices
              + ". Loops with more than 2^31 - 1 vertices are not supported.");
    }

    // Decode the individual vertices.
    List<S2Point> vertices = S2PointCompression.decodePointsCompressed(numVertices, level, decoder);

    // Decode the compression properties.
    CompressedEncodingProperties properties =
        new CompressedEncodingProperties(decoder.readVarint64());

    // Decode the depth.
    int depth = decoder.readVarint32();

    // If the bounds are encoded, decode them and instantiate the loop with all the available
    // information. Otherwise, just create the new loop by passing in the vertices, and let it do
    // the full instantiation.
    S2Loop loop = null;
    if (properties.hasProperty(CompressedEncodingProperties.Property.BOUND_ENCODED)) {
      // TODO(user): Many loops will be small enough that the bound isn't encoded, but we
      // still want to init with the known originInside value to avoid wasted work.
      boolean originInside =
          properties.hasProperty(CompressedEncodingProperties.Property.ORIGIN_INSIDE);

      S2LatLngRect bound = S2LatLngRect.decode(decoder);
      // Since we have all the ingredients, we can speed up construction here.
      // TODO(user): Potential optimization: Avoid copying the contents of vertices when
      // constructing the loop.
      loop = S2Loop.newLoopWithTrustedDetails(vertices, originInside, bound);
    } else {
      loop = new S2Loop(vertices);
    }

    // Set the depth explicitly, since it's not set during instantiation.
    loop.setDepth(depth);

    return loop;
  }

  private static S2Loop decodeInternal(LittleEndianInput decoder) throws IOException {
    int numVertices = decoder.readInt();
    if (numVertices < 0) {
      throw new IOException(
          "Invalid numVertices: "
              + numVertices
              + ". Loops with more than 2^31 - 1 vertices are not supported.");
    }

    ArrayList<S2Point> vertices = new ArrayList<>(numVertices);
    for (int i = 0; i < numVertices; i++) {
      vertices.add(S2Point.decode(decoder));
    }

    boolean originInside = decoder.readByte() != 0;
    int depth = decoder.readInt();
    S2LatLngRect bound = S2LatLngRect.decode(decoder);
    // TODO(user): Potential optimization: Avoid copying the contents of vertices when
    // constructing the loop.
    S2Loop loop = S2Loop.newLoopWithTrustedDetails(vertices, originInside, bound);
    loop.setDepth(depth);

    // An initialized loop will have some non-zero count of vertices. An uninitialized loop has zero
    // vertices. This code supports encoding and decoding of uninitialized loops, but we only want
    // to call InitIndex for initialized loops. Otherwise we defer InitIndex until the call to
    // Init().
    if (numVertices > 0) {
      loop.initIndex();
    }

    return loop;
  }

  private void encodeInternal(LittleEndianOutput encoder) throws IOException {
    encoder.writeInt(numVertices);
    for (int i = 0; i < numVertices; i++) {
      vertex(i).encode(encoder);
    }
    encoder.writeByte(isOriginInside() ? (byte) 1 : (byte) 0);
    encoder.writeInt(depth);
    bound.encode(encoder);
  }

  /** Encodes this S2Loop using the lossless encoding. */
  void encode(LittleEndianOutput encoder) throws IOException {
    // Only LOSSLESS encoding is supported.
    encoder.writeByte(LOSSLESS_ENCODING_VERSION);
    encodeInternal(encoder);
  }

  /**
   * Returns a loop decoded from the given stream. Note S2Loops are intended to be serialized as
   * part of an S2Polygon; see {@link S2Polygon#decode(java.io.InputStream)}.
   *
   * <p>Also note that the decoded loop is not checked for validity. Clients may call {@link
   * S2Loop#isValid()} if they are decoding from an untrusted source.
   */
  static S2Loop decode(LittleEndianInput decoder) throws IOException {
    byte version = decoder.readByte();
    switch (version) {
      case LOSSLESS_ENCODING_VERSION:
        return decodeInternal(decoder);

      default:
        throw new IOException(
            "Unknown S2Loop encoding version encountered during decoding: " + version);
    }
  }

  /**
   * This method checks all edges of loop A for intersection against all edges of loop B. If there
   * is any shared vertex, the wedges centered at this vertex are set to {@code relation}.
   */
  private static boolean hasCrossingRelation(S2Loop a, S2Loop b, LoopRelation relation) {
    // We look for S2CellId ranges where the indexes of A and B overlap, and then test those edges
    // for crossings.
    LoopRangeIterator ai = new LoopRangeIterator(a.index);
    LoopRangeIterator bi = new LoopRangeIterator(b.index);
    // Tests edges of A against B.
    LoopCrosser ab = new LoopCrosser(a, b, relation, false);
    // Tests edges of B against A.
    LoopCrosser ba = new LoopCrosser(b, a, relation, true);
    while (!ai.done() || !bi.done()) {
      if (ai.rangeMax().lessThan(bi.rangeMin())) {
        // The A and B cells don't overlap, and A precedes B.
        ai.seekTo(bi);
      } else if (bi.rangeMax().lessThan(ai.rangeMin())) {
        // The A and B cells don't overlap, and B precedes A.
        bi.seekTo(ai);
      } else {
        // One cell contains the other. Determine which cell is larger.
        long abRelation = UnsignedLongs.compare(ai.id().lowestOnBit(), bi.id().lowestOnBit());
        if (abRelation > 0) {
          // A's index cell is larger.
          if (ab.hasCrossingRelation(ai, bi)) {
            return true;
          }
        } else if (abRelation < 0) {
          // B's index cell is larger.
          if (ba.hasCrossingRelation(bi, ai)) {
            return true;
          }
        } else {
          // The A and B cells are the same. Since the two cells have the same center point P, check
          // whether P satisfies the crossing targets.
          if (ab.aCrossingTarget() == (ai.containsCenter() ? 1 : 0)
              && ab.bCrossingTarget() == (bi.containsCenter() ? 1 : 0)) {
            return true;
          }
          // Otherwise test all the edge crossings directly.
          if (ai.numEdges() > 0
              && bi.numEdges() > 0
              && ab.cellCrossesCell(ai.clipped(), bi.clipped())) {
            return true;
          }
          ai.next();
          bi.next();
        }
      }
    }
    return false;
  }

  private static class LoopRangeIterator extends RangeIterator<Cell> {
    private S2ClippedShape clipped;

    public LoopRangeIterator(S2ShapeIndex index) {
      super(index.iterator());
    }

    @Override
    protected void refresh() {
      super.refresh();
      clipped = it.done() ? null : it.entry().clipped(0);
    }

    /** Various other convenience methods for the current cell. */
    public S2ClippedShape clipped() {
      return clipped;
    }

    public int numEdges() {
      return clipped().numEdges();
    }

    public boolean containsCenter() {
      return clipped().containsCenter();
    }
  }

  /**
   * LoopCrosser is a helper class for determining whether two loops cross. It is instantiated twice
   * for each pair of loops to be tested, once for the pair (A, B) and once for the pair (B, A), in
   * order to be able to process edges in either loop nesting order.
   */
  private static final class LoopCrosser {
    private final S2Loop a;
    private final S2Loop b;
    private final LoopRelation relation;
    private final boolean swapped;
    private final int aCrossingTarget;
    private final int bCrossingTarget;

    // State maintained by startEdge() and edgeCrossesCell().
    private S2EdgeUtil.EdgeCrosser crosser;
    private int aj;
    private int bjPrev;

    // Temporary data declared here to avoid repeated memory allocations.
    private final S2CrossingEdgeQuery bQuery;
    private final List<S2ShapeIndex.Cell> bCells;

    /**
     * If {@code swapped} is true, the loops A and B have been swapped. This affects how arguments
     * are passed to the given loop relation, since for example A.contains(B) is not the same as
     * B.contains(A).
     */
    public LoopCrosser(S2Loop a, S2Loop b, LoopRelation relation, boolean swapped) {
      this.a = a;
      this.b = b;
      this.relation = relation;
      this.swapped = swapped;
      aCrossingTarget = swapped ? relation.bCrossingTarget() : relation.aCrossingTarget();
      bCrossingTarget = swapped ? relation.aCrossingTarget() : relation.bCrossingTarget();
      bQuery = new S2CrossingEdgeQuery(b.index);
      bCells = Lists.newArrayList();
    }

    /**
     * Returns the crossing targets for the loop relation, taking into account whether the loops
     * have been swapped.
     */
    public int aCrossingTarget() {
      return aCrossingTarget;
    }

    public int bCrossingTarget() {
      return bCrossingTarget;
    }

    /**
     * Given two iterators positioned such that {@code ai.id().contains(bi.id())}, returns true if
     * there is a crossing relationship anywhere within {@code ai.id()}. Specifically, this method
     * returns true if there is an edge crossing, a wedge crossing, or a point P that matches both
     * "crossing targets". Advances both iterators past {@code ai.id()}.
     */
    public boolean hasCrossingRelation(LoopRangeIterator ai, LoopRangeIterator bi) {
      assert ai.id().contains(bi.id());
      if (ai.numEdges() == 0) {
        if (aCrossingTarget == (ai.containsCenter() ? 1 : 0)) {
          // All points within ai.id() satisfy the crossing target for A, so it's worth iterating
          // through the cells of B to see whether any cell centers also satisfy the crossing target
          // for B.
          S2CellId maxRange = ai.rangeMax();
          do {
            if (bCrossingTarget == (bi.containsCenter() ? 1 : 0)) {
              return true;
            }
            bi.next();
          } while (bi.id().lessOrEquals(maxRange));
        } else {
          // The crossing target for A is not satisfied, so we skip over the cells of B using
          // binary search.
          bi.seekBeyond(ai);
        }
      } else {
        // The current cell of A has at least one edge, so check for crossings.
        if (hasCrossing(ai, bi)) {
          return true;
        }
      }
      ai.next();
      return false;
    }

    /**
     * Given two index cells, returns true if there are any edge crossings or wedge crossings within
     * those cells.
     */
    public boolean cellCrossesCell(S2ClippedShape aClipped, S2ClippedShape bClipped) {
      // Test all edges of 'aClipped' against all edges of 'bClipped'.
      int aNumClipped = aClipped.numEdges();
      for (int i = 0; i < aNumClipped; ++i) {
        startEdge(aClipped.edge(i));
        if (edgeCrossesCell(bClipped)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Given two iterators positioned such that {@code ai.id().contains(bi.id())}, returns true if
     * there is an edge crossing or a wedge crossing anywhere within {@code ai.id()}. Advances
     * {@code bi} (only) past {@code ai.id()}.
     */
    private boolean hasCrossing(LoopRangeIterator ai, LoopRangeIterator bi) {
      assert ai.id().contains(bi.id());
      // If ai.id() intersects many edges of B, then it is faster to use S2CrossingEdgeQuery to
      // narrow down the candidates. But if it intersects only a few edges, it is faster to check
      // all the crossings directly. We handle this by advancing 'bi' and keeping track of how many
      // edges we would need to test.

      // Tuned using Caliper benchmarking.
      final int edgeQueryMinEdges = 40;
      int totalEdges = 0;
      bCells.clear();
      S2CellId maxRange = ai.rangeMax();
      do {
        if (bi.numEdges() > 0) {
          totalEdges += bi.numEdges();
          if (totalEdges >= edgeQueryMinEdges) {
            // There are too many edges to test them directly, so use S2CrossingEdgeQuery.
            if (cellCrossesAnySubcell(ai.clipped(), ai.id())) {
              return true;
            }
            bi.seekBeyond(ai);
            return false;
          }
          bCells.add(bi.entry());
        }
        bi.next();
      } while (bi.id().lessOrEquals(maxRange));

      // Test all the edge crossings directly.
      for (int c = 0; c < bCells.size(); ++c) {
        if (cellCrossesCell(ai.clipped(), bCells.get(c).clipped(0))) {
          return true;
        }
      }

      return false;
    }

    /**
     * Given an index cell of A, returns true if there are any edge or wedge crossings with any
     * index cell of B contained within {@code bId}.
     */
    private boolean cellCrossesAnySubcell(S2ClippedShape aClipped, S2CellId bId) {
      // Test all edges of 'aClipped' against all edges of B. The relevant B edges are guaranteed to
      // be children of 'bId', which lets us find the correct index cells more efficiently.
      S2PaddedCell bRoot = new S2PaddedCell(bId, 0);
      int aNumClipped = aClipped.numEdges();
      for (int i = 0; i < aNumClipped; ++i) {
        int aj = aClipped.edge(i);
        // Use an S2CrossingEdgeQuery starting at 'bRoot' to find the index cells of B that might
        // contain crossing edges.
        if (!bQuery.getCells(a.vertex(aj), a.vertex(aj + 1), bRoot, bCells)) {
          continue;
        }
        startEdge(aj);
        for (int c = 0; c < bCells.size(); ++c) {
          if (edgeCrossesCell(bCells.get(c).clipped(0))) {
            return true;
          }
        }
      }
      return false;
    }

    /** Prepares to check the given edge of loop A for crossings. */
    private void startEdge(int aj) {
      // Start testing the given edge of A for crossings.
      crosser = new S2EdgeUtil.EdgeCrosser(a.vertex(aj), a.vertex(aj + 1));
      this.aj = aj;
      bjPrev = -2;
    }

    /**
     * Checks the current edge of loop A for crossings with all edges of the given index cell of
     * loop B.
     */
    private boolean edgeCrossesCell(S2ClippedShape bClipped) {
      // Test the current edge of A against all edges of 'bClipped'.
      int bNumClipped = bClipped.numEdges();
      for (int j = 0; j < bNumClipped; ++j) {
        int bj = bClipped.edge(j);
        if (bj != bjPrev + 1) {
          crosser.restartAt(b.vertex(bj));
        }
        bjPrev = bj;
        int crossing = crosser.robustCrossing(b.vertex(bj + 1));
        if (crossing < 0) {
          continue;
        }
        if (crossing > 0) {
          return true;
        }
        // We only need to check each shared vertex once, so we only consider the case where
        // a.vertex(aj + 1).equalsPoint(b.vertex(bj + 1)).
        if (a.vertex(aj + 1).equalsPoint(b.vertex(bj + 1))) {
          if (swapped) {
            if (relation.wedgesCross(
                b.vertex(bj), b.vertex(bj + 1), b.vertex(bj + 2), a.vertex(aj), a.vertex(aj + 2))) {
              return true;
            }
          } else {
            if (relation.wedgesCross(
                a.vertex(aj), a.vertex(aj + 1), a.vertex(aj + 2), b.vertex(bj), b.vertex(bj + 2))) {
              return true;
            }
          }
        }
      }
      return false;
    }
  }

  /** A relation between two loops (e.g. Contains, Intersects, or CompareBoundary.) */
  private interface LoopRelation {
    /**
     * Optionally, {@code aCrossingTarget} and {@code bCrossingTarget} can specify an early-exit
     * condition for the loop relation. If any point P is found such that
     *
     * <p>{@code aCrossingTarget == (a.contains(P) ? 1 : 0) && bCrossingTarget == (b.contains(P) ? 1
     * : 0) }
     *
     * <p>then the loop relation is assumed to be the same as if a pair of crossing edges were
     * found. For example, the contains() relation has
     *
     * <p>{@code aCrossingTarget() == 0 bCrossingTarget() == 1 }
     *
     * <p>because if {@code !a.contains(P)} and {@code b.contains(P)} for any point P, then it is
     * equivalent to finding an edge crossing (i.e., since contains() returns false in both cases).
     *
     * <p>Loop relations that do not have an early-exit condition of this form should return -1 for
     * both crossing targets.
     */
    int aCrossingTarget();

    int bCrossingTarget();

    /**
     * Given a vertex {@code ab1} that is shared between the two loops, returns true if the two
     * associated wedges (a0, ab1, a2) and (b0, ab1, b2) are equivalent to an edge crossing. The
     * loop relation is also allowed to maintain its own internal state, and can return true if it
     * observes any sequence of wedges that are equivalent to an edge crossing.
     */
    boolean wedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2);
  }

  /** Loop relation for contains(). */
  private static final class ContainsRelation implements LoopRelation {
    private boolean foundSharedVertex = false;

    public boolean foundSharedVertex() {
      return foundSharedVertex;
    }

    /**
     * If A.contains(P) == false && B.contains(P) == true, it is equivalent to having an edge
     * crossing (i.e., contains() returns false).
     */
    @Override
    public int aCrossingTarget() {
      // signifies false
      return 0;
    }

    @Override
    public int bCrossingTarget() {
      // signifies true
      return 1;
    }

    @Override
    public boolean wedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      foundSharedVertex = true; // ab1 is shared by definition.
      return new S2EdgeUtil.WedgeContains().test(a0, ab1, a2, b0, b2) != 1;
    }
  }

  /** Loop relation for intersects(). */
  private static final class IntersectsRelation implements LoopRelation {
    private boolean foundSharedVertex = false;

    public boolean foundSharedVertex() {
      return foundSharedVertex;
    }

    /**
     * If A.contains(P) == false && B.contains(P) == true, it is equivalent to having an edge
     * crossing (i.e., intersects() returns true).
     */
    @Override
    public int aCrossingTarget() {
      // signifies true
      return 1;
    }

    @Override
    public int bCrossingTarget() {
      // signifies true
      return 1;
    }

    @Override
    public boolean wedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      foundSharedVertex = true; // ab1 is shared by definition.
      // TODO(user): Make the test methods of each WedgeProcessor available as a static
      // method, so we can call them directly without having to instantiate an object.
      return new S2EdgeUtil.WedgeIntersects().test(a0, ab1, a2, b0, b2) == -1;
    }
  }

  /** Loop relation for compareBoundary(). */
  private static final class CompareBoundaryRelation implements LoopRelation {
    /** True if loop B should be reversed. */
    private final boolean bReversed;

    /** True if any wedge was processed. */
    private boolean foundSharedVertex = false;

    /** True if any edge of B is contained by A. */
    private boolean containsEdge = false;

    /** True if any edge of B is excluded by A. */
    private boolean excludesEdge = false;

    public CompareBoundaryRelation(boolean reverseB) {
      this.bReversed = reverseB;
    }

    public boolean foundSharedVertex() {
      return foundSharedVertex;
    }

    public boolean containsEdge() {
      return containsEdge;
    }

    /**
     * The CompareBoundaryRelation does not have a useful early-exit condition, so we return -1 for
     * both crossing targets.
     *
     * <p>Aside: A possible early exit condition could be based on the following:
     *
     * <ul>
     *   <li>If A contains a point of both B and ~B, then A intersects Boundary(B).
     *   <li>If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
     *   <li>So if the intersections of {A, ~A} with {B, ~B} are all non-empty, the return value is
     *       0, i.e., Boundary(A) intersects Boundary(B).
     * </ul>
     *
     * Unfortunately, it isn't worth detecting this situation because by the time we have seen a
     * point in all four intersection regions, we are also guaranteed to have seen at least one pair
     * of crossing edges.
     */
    @Override
    public int aCrossingTarget() {
      // Signifies no early-exit condition.
      return -1;
    }

    @Override
    public int bCrossingTarget() {
      // Signifies no early-exit condition.
      return -1;
    }

    @Override
    public boolean wedgesCross(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // Because we don't care about the interior of B, only its boundary, it is sufficient to check
      // whether A contains the semiwedge (ab1, b2).
      foundSharedVertex = true;
      if (wedgeContainsSemiwedge(a0, ab1, a2, b2, bReversed)) {
        containsEdge = true;
      } else {
        excludesEdge = true;
      }
      return containsEdge && excludesEdge;
    }
  }
}
