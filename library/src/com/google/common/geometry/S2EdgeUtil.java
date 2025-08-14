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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.DBL_ERROR;
import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2.M_SQRT1_2;
import static com.google.common.geometry.S2.M_SQRT2;
import static com.google.common.geometry.S2.M_SQRT3;
import static com.google.common.geometry.S2.ROBUST_CROSS_PROD_ERROR;
import static com.google.common.geometry.S2.isUnitLength;
import static com.google.common.geometry.S2Point.ZERO;
import static com.google.common.geometry.S2Point.Z_POS;
import static com.google.common.geometry.S2Point.scalarTripleProduct;
import static com.google.common.geometry.S2Predicates.compareEdgeDistance;
import static com.google.common.geometry.S2Predicates.orderedCCW;
import static com.google.common.geometry.S2Predicates.sign;
import static com.google.common.geometry.S2RobustCrossProd.robustCrossProd;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import com.google.common.base.Preconditions;
import com.google.errorprone.annotations.InlineMe;

/**
 * This class contains various utility functions related to edges. It collects together common code
 * that is needed to implement polygonal geometry such as polylines, loops, and general polygons.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@SuppressWarnings({"Assertion", "IdentifierName"})
public class S2EdgeUtil {
  /**
   * IEEE floating-point operations have a maximum error of 0.5 ULPS (units in the last place). For
   * double-precision numbers, this works out to 2**-53 (about 1.11e-16) times the magnitude of the
   * result. It is possible to analyze the calculation done by getIntersection() and work out the
   * worst-case rounding error. I have done a rough version of this, and my estimate is that the
   * worst case distance from the intersection point X to the great circle through (a0, a1) is about
   * 12 ULPS, or about 1.332e-15. This needs to be increased by a factor of (1/0.866) to account for
   * the edgeSpliceFraction() in S2PolygonBuilder.
   */
  public static final S1Angle DEFAULT_INTERSECTION_TOLERANCE =
      S1Angle.radians(12 * DBL_EPSILON / 2 / 0.866);

  /**
   * Threshold for small angles, that help lenientCrossing to determine whether two edges are likely
   * to intersect.
   */
  private static final double MAX_DET_ERROR = 1e-14;

  /**
   * The maximum angle between a returned vertex and the nearest point on the exact edge AB. It is
   * equal to the maximum directional error in {@link S2RobustCrossProd#robustCrossProd(S2Point,
   * S2Point)}, plus the error when projecting points onto a cube face.
   */
  public static final double FACE_CLIP_ERROR_RADIANS = 3 * DBL_EPSILON;

  /**
   * The same angle as {@link #FACE_CLIP_ERROR_RADIANS}, expressed as a maximum distance in
   * (u,v)-space. In other words, a returned vertex is at most this far from the exact edge AB
   * projected into (u,v)-space.
   */
  public static final double FACE_CLIP_ERROR_UV_DIST = 9 * DBL_EPSILON;

  /**
   * The same angle as {@link #FACE_CLIP_ERROR_RADIANS}, expressed as the maximum error in an
   * individual u- or v-coordinate. In other words, for each returned vertex there is a point on the
   * exact edge AB whose u- and v-coordinates differ from the vertex by at most this amount.
   */
  public static final double FACE_CLIP_ERROR_UV_COORD = 9 * M_SQRT1_2 * DBL_EPSILON;

  /**
   * The maximum error in IntersectRect. If some point of AB is inside the rectangle by at least
   * this distance, the result is guaranteed to be true; if all points of AB are outside the
   * rectangle by at least this distance, the result is guaranteed to be false. This bound assumes
   * that "rect" is a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it (e.g., by
   * 1e-10 or less).
   */
  public static final double INTERSECTS_RECT_ERROR_UV_DIST = 3 * M_SQRT2 * DBL_EPSILON;

  /**
   * The maximum error in a clipped point's u- or v-coordinate compared to the exact result,
   * assuming that the points A and B are in the rectangle [-1,1]x[1,1] or slightly outside it (by
   * 1e-10 or less).
   */
  public static final double EDGE_CLIP_ERROR_UV_COORD = 2.25 * DBL_EPSILON;

  /**
   * The maximum error between a clipped edge or boundary point and the corresponding exact result.
   * It is equal to the error in a single coordinate because at most one coordinate is subject to
   * error.
   */
  public static final double EDGE_CLIP_ERROR_UV_DIST = 2.25 * DBL_EPSILON;

  /** Max error allowed when checking if a loop boundary approximately intersects a target cell. */
  public static final double MAX_CELL_EDGE_ERROR =
      FACE_CLIP_ERROR_UV_COORD + INTERSECTS_RECT_ERROR_UV_DIST;

  // TODO(user): In C++, INTERSECTION_ERROR is 8 * DBL_ERR, where DBL_ERR is half
  // DBL_EPSILON. The value here is twice as large. However, changing it may have performance
  // implications.
  /**
   * An upper bound on the distance (in radians) from the intersection point returned by {@link
   * #getIntersection(S2Point, S2Point, S2Point, S2Point)} to the true intersection point.
   *
   * <p>{@code INTERSECTION_ERROR} can be set somewhat arbitrarily, because the algorithm uses more
   * precision than necessary in order to achieve the specified error. The only strict requirement
   * is that {@code INTERSECTION_ERROR >= 2 * DBL_EPSILON} radians. However, using a larger error
   * tolerance makes the algorithm more efficient because it reduces the number of cases where exact
   * arithmetic is needed.
   */
  public static final double INTERSECTION_ERROR = 8 * DBL_EPSILON;

  /**
   * GET_POINT_ON_LINE_ERROR is an upper bound on the distance between the point returned by
   * {@link #getPointOnLine(S2Point, S2Point, S1Angle)} and the corresponding true infinite-
   * precision result.
   */
  public static final S1Angle GET_POINT_ON_LINE_ERROR =
      S1Angle.radians((4 + (2 / M_SQRT3)) * DBL_ERROR).add(ROBUST_CROSS_PROD_ERROR);

  /**
   * PROJECT_PERPENDICULAR_ERROR is an upper bound on the distance from the point returned by
   * {@link #project(S2Point, S2Point, S2Point)} to the edge AB. Note that it only bounds the error
   * perpendicular to the edge, not the error parallel to it.
   */
  public static final S1Angle PROJECT_PERPENDICULAR_ERROR =
      S1Angle.radians((2 + (2 / M_SQRT3)) * DBL_ERROR).add(ROBUST_CROSS_PROD_ERROR);

  /**
   * GET_POINT_ON_RAY_PERPENDICULAR_ERROR is an upper bound on the distance from the point returned
   * by getPointOnRay() to the ray itself. Note that it only bounds the error perpendicular to the
   * ray, not the error parallel to it.
   */
  public static final S1Angle GET_POINT_ON_RAY_PERPENDICULAR_ERROR = S1Angle.radians(
      3 * DBL_ERROR);

  /**
   * This value can be used as the S2Builder snapRadius() to ensure that edges that have been
   * displaced by up to INTERSECTION_ERROR are merged back together again. For example, this can
   * happen when geometry is intersected with a set of tiles and then unioned. It is equal to twice
   * the intersection error because input edges might have been displaced in opposite directions.
   */
  public static final S1Angle INTERSECTION_MERGE_RADIUS = S1Angle.radians(2 * INTERSECTION_ERROR);

  /** Used to denote which point should be used when finding distances/points. */
  private enum ClosestPoint {
    A0,
    A1,
    B0,
    B1,
    NONE
  }

  /**
   * Used to efficiently test a fixed edge AB against an edge chain. To use it, {@link
   * #init(S2Point, S2Point) initialize} with the edge AB, and call
   * {@link #robustCrossing(S2Point, S2Point)} or {@link #edgeOrVertexCrossing(S2Point, S2Point)}
   * with each edge of the chain.
   *
   * <p>This class is <strong>not</strong> thread-safe.
   */
  public static final class EdgeCrosser {
    private S2Point a;
    private S2Point b;
    private S2Point aCrossB;

    /** Previous vertex in the vertex chain. */
    private S2Point c;

    /** The orientation of the triangle ACB, i.e. the orientation around the current vertex. */
    private int acb;

    /**
     * The orientation of triangle BDA. This is used to return an extra value from
     * robustCrossingInternal().
     */
    int bdaReturn;

    /**
     * True if the tangents have been computed. To reduce the number of calls to
     * {@link S2Predicates.Sign#expensive(S2Point, S2Point, S2Point, boolean)}, we compute an
     * outward-facing tangent at A and B if necessary. If the plane perpendicular to one of these
     * tangents separates AB from CD (i.e., one edge on each side) then there is no intersection.
     */
    private boolean haveTangents;

    /** Outward-facing tangent at A. */
    private S2Point aTangent;

    /** Outward-facing tangent at B. */
    private S2Point bTangent;

    /**
     * Constructs an uninitialized edge crosser. Invoke {@link #init(S2Point, S2Point)} before
     * calling the other methods.
     */
    public EdgeCrosser() {}

    /** Convenience constructor that calls init() with the given fixed edge AB. */
    public EdgeCrosser(S2Point a, S2Point b) {
      init(a, b);
    }

    /**
     * AB is the given fixed edge, and C is the first vertex of the vertex chain. Equivalent to
     * using the two-arg constructor and calling restartAt(c).
     */
    public EdgeCrosser(S2Point a, S2Point b, S2Point c) {
      this(a, b);
      restartAt(c);
    }

    /**
     * Initialize this edge crosser with the given endpoints, available via {@link #a()} and
     * {@link #b()}, and computes the approximate edge normal, available via {@link #normal()}.
     */
    public void init(S2Point a, S2Point b) {
      assert S2.skipAssertions || S2.isUnitLength(a);
      assert S2.skipAssertions || S2.isUnitLength(b);
      this.a = a;
      this.b = b;
      this.c = null;
      this.aCrossB = a.crossProd(b);
      this.haveTangents = false;
    }

    /** Returns the first point passed to {#link #init). */
    public S2Point a() {
      return a;
    }

    /** Returns the second point passed to {@link #init}. */
    public S2Point b() {
      return b;
    }

    /** Returns the last 'c' point checked. */
    public S2Point c() {
      return c;
    }

    /** Returns the approximate normal of the AB edge. */
    public S2Point normal() {
      return aCrossB;
    }

    /** Call this method when your chain 'jumps' to a new place. */
    public void restartAt(S2Point c) {
      assert S2.skipAssertions || S2.isUnitLength(c);
      this.c = c;
      acb = -S2Predicates.Sign.triage(aCrossB, c);
    }

    /**
     * This method is equivalent to calling the {@link #robustCrossing(S2Point, S2Point)} function
     * (defined below) on the edges AB and CD. It returns +1 if there is a crossing, -1 if there is
     * no crossing, and 0 if two points from different edges are the same. If an edge is degenerate
     * (A == B or C == D), the return value is 0 if two vertices from different edges are the same
     * and -1 otherwise. As a side effect, it saves vertex D to be used as the next vertex C.
     */
    public int robustCrossing(S2Point d) {
      assert S2.skipAssertions || S2.isUnitLength(d);
      // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must all be oriented the
      // same way (CW or CCW). We keep the orientation of ACB as part of our state. When each new
      // point D arrives, we compute the orientation of BDA and check whether it matches ACB. This
      // checks whether the points C and D are on opposite sides of the great circle through AB.

      // Recall that triageSign is invariant with respect to rotating its arguments, i.e. ABD has
      // the same orientation as BDA.
      int bda = S2Predicates.Sign.triage(aCrossB, d);
      if (this.acb == -bda && bda != 0) {
        // The most common case: triangles have opposite orientations. Save the current vertex D as
        // the next vertex C, and also save the orientation of the new triangle ACB (which is
        // opposite to the current triangle BDA).
        this.c = d;
        this.acb = -bda;
        return -1;
      }
      this.bdaReturn = bda;
      return robustCrossingInternal(d);
    }

    /**
     * As {@link #robustCrossing(S2Point)}, but restarts at {@code c} if that is not the previous
     * endpoint. Returns +1 if there is a crossing, -1 if there is no crossing, and 0 if two points
     * from different edges are the same. If an edge is degenerate (A == B or C == D), the return
     * value is 0 if two vertices from different edges are the same and -1 otherwise.
     */
    @SuppressWarnings("ReferenceEquality")
    public int robustCrossing(S2Point c, S2Point d) {
      if (this.c != c) {
        // Comparison by reference may sometimes cause us to do slightly extra work, but the vast
        // majority of the time if the points are equal by value, they are exactly the same
        // reference as well.
        restartAt(c);
      }
      return robustCrossing(d);
    }

    /**
     * Returns true if AB crosses CD, either within the edge or by a vertex crossing at a shared
     * vertex. This method is equivalent to the {@link #edgeOrVertexCrossing(S2Point, S2Point)}
     * method defined below. It is similar to {@link #robustCrossing(S2Point, S2Point)}, but handles
     * cases where two vertices are identical in a way that makes it easy to implement
     * point-in-polygon containment tests. Like {@link #robustCrossing(S2Point, S2Point)}, it saves
     * the vertex passed to the method as the next vertex C.
     */
    public boolean edgeOrVertexCrossing(S2Point d) {
      // Copy c, since the reference may be replaced by robustCrossing().
      S2Point c2 = c;

      int crossing = robustCrossing(d);
      if (crossing < 0) {
        return false;
      }
      if (crossing > 0) {
        return true;
      }

      return vertexCrossing(a, b, c2, d);
    }

    /**
     * Returns true if AB crosses CD, either within the edge or by a vertex crossing at a shared
     * vertex. Like {@link #edgeOrVertexCrossing(S2Point)}, but restarts at {@code c} if that is not
     * the previous endpoint.
     */
    @SuppressWarnings("ReferenceEquality")
    public boolean edgeOrVertexCrossing(S2Point c, S2Point d) {
      // Test by reference since the same value in different references is very rare.
      if (this.c != c) {
        restartAt(c);
      }
      return edgeOrVertexCrossing(d);
    }

    /**
     * Like edgeOrVertexCrossing(), but returns -1 if AB crosses CD from left to right, +1 if AB
     * crosses CD from right to left, and 0 otherwise. This implies that if CD bounds some region
     * according to the "interior is on the left" rule, this function returns -1 when AB exits the
     * region and +1 when AB enters.
     *
     * <p>This method allows computing the change in winding number from point A to point B by
     * summing the signed edge crossings of AB with the edges of the loop(s) used to define the
     * winding number.
     *
     * <p>
     */
    @SuppressWarnings("ReferenceEquality")
    public int signedEdgeOrVertexCrossing(S2Point c, S2Point d) {
      // Test by reference since the same value in different references is very rare.
      if (this.c != c) {
        restartAt(c);
      }
      return signedEdgeOrVertexCrossing(d);
    }

    /**
     * Like edgeOrVertexCrossing above, but uses the last vertex passed to one of the crossing
     * methods (or restartAt) as the first vertex of the current edge.
     */
    public int signedEdgeOrVertexCrossing(S2Point d) {
      // We need to copy c since it is changed by robustCrossing().
      S2Point c2 = c;

      int crossing = robustCrossing(d);
      if (crossing < 0) {
        return 0;
      }
      if (crossing > 0) {
        // When AB crosses CD, the crossing sign is sign(ABC). EdgeCrosser doesn't store this, but
        // it does store the sign of the *next* triangle ACB. These two values happen to be the
        // same.
        return acb;
      }
      return signedVertexCrossing(a, b, c2, d);
    }

    /**
     * If the preceding call to robustCrossing() returned +1 (indicating that the edge crosses the
     * edge CD), this method returns -1 if AB crossed CD from left to right and +1 if AB crossed CD
     * from right to left. Otherwise its return value is undefined.
     */
    public int lastInteriorCrossingSign() {
      // When AB crosses CD, the crossing sign is sign(ABC). EdgeCrosser doesn't store this, but
      // it does store the sign of the *next* triangle ACB. These two values happen to be the
      // same.
      return acb;
    }

    /**
     * Compute the actual result, and then save the current vertex D as the next vertex C, and save
     * the orientation of the next triangle ACB (which is opposite to the current triangle BDA).
     */
    private int robustCrossingInternal(S2Point d) {
      int result = robustCrossingInternal2(d);
      this.c = d;
      this.acb = -bdaReturn;
      return result;
    }

    private int robustCrossingInternal2(S2Point d) {
      // At this point it is still very likely that CD does not cross AB. Two common situations are
      // (1) CD crosses the great circle through AB but does not cross AB itself, or (2) A,B,C,D are
      // four points on a line such that AB does not overlap CD. For example, the latter happens
      // when a line or curve is sampled finely, or when geometry is constructed by computing the
      // union of S2CellIds.
      //
      // Most of the time, we can determine that AB and CD do not intersect by computing the two
      // outward-facing tangents at A and B (parallel to AB) and testing whether AB and CD are on
      // opposite sides of the plane perpendicular to one of these tangents. This is somewhat
      // expensive but still much cheaper S2Predicates.Sign.expensive().
      if (!haveTangents) {
        S2Point norm = robustCrossProd(a, b).normalize();
        aTangent = a.crossProd(norm);
        bTangent = norm.crossProd(b);
        haveTangents = true;
      }

      // The error in robustCrossProd() is insignificant. The maximum error in the call to
      // crossProd() (i.e. the maximum norm of the error vector) is (0.5 + 1/sqrt(3)) * DBL_EPSILON.
      // The maximum error in each call to dotProd() below is DBL_EPSILON. (There is also a small
      // relative error term that is insignificant because we are comparing the result against a
      // constant that is very close to zero.)
      // TODO(user): Test coverage to ensure that this is sufficiently large for correctness.
      final double kError = (1.5 + 1 / sqrt(3)) * DBL_EPSILON;
      if ((c.dotProd(aTangent) > kError && d.dotProd(aTangent) > kError)
          || (c.dotProd(bTangent) > kError && d.dotProd(bTangent) > kError)) {
        return -1;
      }

      // Otherwise, eliminate the cases where two vertices from different edges are equal. (These
      // cases could be handled in the code below, but we would rather avoid calling Sign.expensive
      // whenever possible.)
      if (a.equalsPoint(c) || a.equalsPoint(d) || b.equalsPoint(c) || b.equalsPoint(d)) {
        return 0;
      }

      // Eliminate cases where an input edge is degenerate. (Note that in most cases, if CD is
      // degenerate then this method is not even called because acb and bda have different signs.)
      if (a.equalsPoint(b) || c.equalsPoint(d)) {
        return -1;
      }

      // Otherwise it's time to break out the big guns.
      if (acb == 0) {
        acb = -S2Predicates.Sign.expensive(a, b, c, true);
        assert acb != 0;
      }
      if (bdaReturn == 0) {
        bdaReturn = S2Predicates.Sign.expensive(a, b, d, true);
        assert bdaReturn != 0;
      }
      if (bdaReturn != acb) {
        return -1;
      }

      S2Point cCrossD = c.crossProd(d);
      int cbd = -sign(c, d, b, cCrossD);
      assert cbd != 0;
      if (cbd != acb) {
        return -1;
      }

      int dac = sign(c, d, a, cCrossD);
      assert dac != 0;
      return (dac == acb) ? 1 : -1;
    }

    /** Helper that checks the sign of ABC, using a precomputed cross product for AxB. */
    private static int sign(S2Point a, S2Point b, S2Point c, S2Point aCrossB) {
      // Here we assume that 'a', 'b', and 'c' are unit length, as they were asserted to be so when
      // supplied to the public methods above.
      int ccw = S2Predicates.Sign.triage(aCrossB, c);
      if (ccw == 0) {
        ccw = S2Predicates.Sign.expensive(a, b, c, true);
      }
      return ccw;
    }
  }

  /**
   * This class computes a bounding rectangle that contains all edges defined by a vertex chain v0,
   * v1, v2, ... All vertices must be unit length. Note that the bounding rectangle of an edge can
   * be larger than the bounding rectangle of its endpoints, e.g. consider an edge that passes
   * through the north pole.
   *
   * <p>The bounds are calculated conservatively to account for numerical errors when S2Points are
   * converted to S2LatLngs. For example, this class guarantees that if L is a closed edge chain (a
   * loop) such that the interior of the loop does not contain either pole, and P is any point such
   * that L contains P, then the RectBounder of all edges in L will contain S2LatLng(P).
   */
  public static class RectBounder {
    /** The accumulated bounds, initially empty. */
    private S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();

    /** The previous vertex in the chain. */
    private S2Point a;

    /** The corresponding latitude-longitude. */
    private S2LatLng aLatLng;

    /** Temporary storage for the longitude range spanned by AB. */
    private final S1Interval lngAB = new S1Interval();

    /** Temporary storage for the latitude range spanned by AB. */
    private final R1Interval latAB = new R1Interval();

    public RectBounder() {}

    /**
     * This method is called to add each vertex to the chain. This method is much faster than
     * {@link #addPoint(S2Point)}, since converting S2LatLng to an S2Point is much faster than the
     * other way around.
     */
    public void addPoint(S2LatLng b) {
      addPoint(b.toPoint(), b);
    }

    /**
     * This method is called to add each vertex to the chain. Prefer calling
     * {@link #addPoint(S2LatLng)} if you have that type available. The point must be unit length.
     */
    public void addPoint(S2Point b) {
      addPoint(b, new S2LatLng(b));
    }

    /**
     * Internal implementation of addPoint that takes both the point and latLng representation, by
     * whichever path provided them, and expands the bounds accordingly.
     */
    private void addPoint(S2Point b, S2LatLng bLatLng) {
      assert S2.skipAssertions || isUnitLength(b);
      // TODO(user): Verify the analysis and bounds from the C++ implementation are valid
      // here.
      if (builder.isEmpty()) {
        builder.addPoint(bLatLng);
      } else {
        // First compute the cross product N = A x B robustly. This is the normal to the great
        // circle through A and B. We don't use S2RobustCrossProd.robustCrossProd() since that
        // method returns an arbitrary vector orthogonal to A if the two vectors are proportional,
        // and we want the zero vector in that case.
        // N = 2 * (A x B)
        S2Point n = a.sub(b).crossProd(a.add(b));

        // The relative error in N gets large as its norm gets very small (i.e. when the two points
        // are nearly identical or antipodal). We handle this by choosing a maximum allowable error,
        // and if the error is greater than this we fall back to a different technique. Since it
        // turns out that the other sources of error add up to at most 1.16 * DBL_EPSILON, and it is
        // desirable to have the total error be a multiple of DBL_EPSILON, we have chosen the
        // maximum error threshold here to be 3.84 * DBL_EPSILON. It is possible to show that the
        // error is less than this when
        //
        //   n.norm() >= 8 * sqrt(3) / (3.84 - 0.5 - sqrt(3)) * DBL_EPSILON
        //            = 1.91346e-15 (about 8.618 * DBL_EPSILON)
        double nNorm = n.norm();
        if (nNorm < 1.91346e-15) {
          // A and B are either nearly identical or nearly antipodal (to within 4.309 * DBL_EPSILON,
          // or about 6 nanometers on the earth's surface).
          if (a.dotProd(b) < 0) {
            // The two points are nearly antipodal. The easiest solution is to assume that the edge
            // between A and B could go in any direction around the sphere.
            builder.setFull();
          } else {
            // The two points are nearly identical (to within 4.309 * DBL_EPSILON). In this case we
            // can just use the bounding rectangle of the points, since after the expansion done by
            // getBound() this rectangle is guaranteed to include the (lat,lng) values of all points
            // along AB.
            builder.union(S2LatLngRect.fromPointPair(aLatLng, bLatLng));
          }
        } else {
          // Compute the longitude range spanned by AB.
          lngAB.initFromPointPair(aLatLng.lng().radians(), bLatLng.lng().radians());
          if (lngAB.getLength() >= PI - 2 * DBL_EPSILON) {
            // The points lie on nearly opposite lines of longitude to within the maximum error of
            // the calculation. (Note that this test relies on the fact that Math.PI is slightly
            // less than the true value of Pi, and that representable values near PI are 2 *
            // DBL_EPSILON apart.) The easiest solution is to assume that AB could go on either side
            // of the pole.
            lngAB.setFull();
          }

          // Next we compute the latitude range spanned by the edge AB. We start with the range
          // spanning the two endpoints of the edge:
          latAB.initFromPointPair(aLatLng.lat().radians(), bLatLng.lat().radians());

          // This is the desired range unless the edge AB crosses the plane through N and the Z-axis
          // (which is where the great circle through A and B attains its minimum and maximum
          // latitudes). To test whether AB crosses this plane, we compute a vector M perpendicular
          // to this plane and then project A and B onto it.
          S2Point m = n.crossProd(Z_POS);
          double mDotA = m.dotProd(a);
          double mDotB = m.dotProd(b);

          // We want to test the signs of "mDotA" and "mDotB", so we need to bound the error in
          // these calculations. It is possible to show that the total error is bounded by
          //
          //  (1 + sqrt(3)) * DBL_EPSILON * nNorm + 8 * sqrt(3) * (DBL_EPSILON**2)
          //    = 6.06638e-16 * nNorm + 6.83174e-31
          double mError = 6.06638e-16 * nNorm + 6.83174e-31;

          if (mDotA * mDotB < 0 || abs(mDotA) <= mError || abs(mDotB) <= mError) {
            // Minimum/maximum latitude *may* occur in the edge interior.
            //
            // The maximum latitude is 90 degrees minus the latitude of N. We compute this directly
            // using atan2 in order to get maximum accuracy near the poles.
            //
            // Our goal is compute a bound that contains the computed latitudes of all S2Points P
            // that pass the point-in-polygon containment test. There are three sources of error we
            // need to consider:
            //  - the directional error in N (at most 3.84 * DBL_EPSILON)
            //  - converting N to a maximum latitude
            //  - computing the latitude of the test point P
            // The latter two sources of error are at most 0.955 * DBL_EPSILON individually, but it
            // is possible to show by a more complex analysis that together they can add up to at
            // most 1.16 * DBL_EPSILON, for a total error of 5 * DBL_EPSILON.
            //
            // We add 3 * DBL_EPSILON to the bound here, and getBound() will pad the bound by
            // another 2 * DBL_EPSILON.
            double maxLat =
                min(
                    M_PI_2,
                    3 * DBL_EPSILON
                        + atan2(sqrt(n.getX() * n.getX() + n.getY() * n.getY()), abs(n.getZ())));

            // In order to get tight bounds when the two points are close together, we also bound
            // the min/max latitude relative to the latitudes of the endpoints A and B. First we
            // compute the distance between A and B, and then we compute the maximum change in
            // latitude between any two points along the great circle that are separated by this
            // distance. This gives us a latitude change "budget". Some of this budget must be
            // spent getting from A to B; the remainder bounds the round-trip distance (in latitude)
            // from A or B to the min or max latitude attained along the edge AB.
            double latBudget = 2 * asin(0.5 * a.sub(b).norm() * sin(maxLat));
            // TODO(user): Ensure that this maxDelta value is valid for Java.
            double maxDelta = 0.5 * (latBudget - latAB.getLength()) + DBL_EPSILON;

            // Test whether AB passes through the point of maximum latitude or minimum latitude. If
            // the dot product(s) are small enough then the result may be ambiguous.
            if (mDotA <= mError && mDotB >= -mError) {
              latAB.setHi(min(maxLat, latAB.hi() + maxDelta));
            }
            if (mDotB <= mError && mDotA >= -mError) {
              latAB.setLo(max(-maxLat, latAB.lo() - maxDelta));
            }
          }
          builder.union(new S2LatLngRect(latAB, lngAB));
        }
      }
      a = b;
      aLatLng = bLatLng;
    }

    /**
     * Returns the bounding rectangle of the edge chain that connects the vertices defined so far.
     * TODO(user): Verify the analysis and bounds from the C++ implementation are valid here.
     *
     * <p>This bound satisfies the guarantee made above, i.e. if the edge chain defines a loop, then
     * the bound contains the S2LatLng coordinates of all S2Points contained by the loop.
     */
    public S2LatLngRect getBound() {
      // To save time, we ignore numerical errors in the computed S2LatLngs while accumulating the
      // bounds and then account for them here.
      //
      // Two claims are made here for the accuracy of constructing S2LatLngs from S2Points:

      // 1. S2LatLng(S2Point) has a maximum error of 0.955 * DBL_EPSILON in latitude. In the worst
      // case, we might have rounded "inwards" when computing the bound and "outwards" when
      // computing the latitude of a contained point P, therefore we expand the latitude bounds by
      // 2 * DBL_EPSILON in each direction. (A more complex analysis shows that 1.5 * DBL_EPSILON
      // is enough, but the expansion amount should be a multiple of DBL_EPSILON in order to avoid
      // rounding errors during the expansion itself.)
      //
      // 2. S2LatLng(S2Point) has a maximum error of S2.DBL_EPSILON in longitude, which is simply
      // the maximum rounding error for results in the range [-Pi, Pi].
      //
      // In the C++ implementation these claims are true because the Gnu implementation of atan2()
      // comes from the IBM Accurate Mathematical Library, which implements correct rounding for
      // this intrinsic (i.e., it returns the infinite precision result rounded to the nearest
      // representable value, with ties rounded to even values). This implies that we don't need to
      // expand the longitude bounds at all, since we only guarantee that the bound contains the
      // *rounded* latitudes of contained points. The *true* latitudes of contained points may lie
      // up to S2.DBL_EPSILON outside of the returned bound.
      //
      // However, the Java Math specification for atan2() requires only that the computed result be
      // within 2 ulps of the exact result, which is less precise, so it is not clear if this
      // expansion is sufficient.
      S2LatLng expansion = S2LatLng.fromRadians(2 * DBL_EPSILON, 0);
      return builder.build().expanded(expansion).polarClosure();
    }

    /**
     * Returns the maximum error in getBound() provided that the result does not include either
     * pole. It is only to be used for testing purposes (e.g., by passing it to
     * {@link S2LatLngRect#approxEquals(S2LatLngRect)}).
     *
     * TODO(user): Verify the analysis and bounds from the C++ implementation are valid here.
     */
    static S2LatLng maxErrorForTests() {
      // The maximum error in the latitude calculation is
      //    3.84 * DBL_EPSILON   for the robustCrossProd calculation
      //    0.96 * DBL_EPSILON   for the latitude() calculation
      //    5    * DBL_EPSILON   added by AddPoint/GetBound to compensate for error
      //    ------------------
      //    9.80 * DBL_EPSILON   maximum error in result
      //
      // The maximum error in the longitude calculation is DBL_EPSILON. GetBound does not do any
      // expansion because this isn't necessary in order to bound the *rounded* longitudes of
      // contained points.
      return S2LatLng.fromRadians(10 * DBL_EPSILON, 1 * DBL_EPSILON);
    }

    /**
     * Expand a bound returned by getBound() so that it is guaranteed to contain the bounds of any
     * subregion whose bounds are computed using this class. For example, consider a loop L that
     * defines a square. GetBound() ensures that if a point P is contained by this square, then
     * S2LatLng(P) is contained by the bound. But now consider a diamond shaped loop S contained by
     * L. It is possible that GetBound() returns a larger* bound for S than it does for L, due to
     * rounding errors. This method expands the bound for L so that it is guaranteed to contain the
     * bounds of any subregion S.
     *
     * <p>More precisely, if L is a loop that does not contain either pole, and S is a loop such
     * that {@code L.contains(S)}, then {@code
     * expandForSubregions(RectBound(L)).contains(RectBound(S))}.
     *
     * <p>TODO(user): Verify the analysis and bounds from the C++ implementation are valid
     * here.
     */
    static S2LatLngRect expandForSubregions(S2LatLngRect bound) {
      // Empty bounds don't need expansion.
      if (bound.isEmpty()) {
        return bound;
      }

      // First we need to check whether the bound B contains any nearly-antipodal points (to within
      // 4.309 * DBL_EPSILON). If so then we need to return S2LatLngRect.full(), since the subregion
      // might have an edge between two such points, and addPoint() returns full() for such edges.
      // Note that this can happen even if B is not full(); for example, consider a loop that
      // defines a 10km strip straddling the equator extending from longitudes -100 to +100 degrees.
      //
      // It is easy to check whether B contains any antipodal points, but checking for
      // nearly-antipodal points is trickier. Essentially we consider the original bound B and its
      // reflection through the origin B', and then test whether the minimum distance between B and
      // B' is less than 4.309 * DBL_EPSILON.

      // "lngGap" is a lower bound on the longitudinal distance between B and its reflection B'.
      // (2.5 * DBL_EPSILON is the maximum combined error of the endpoint longitude calculations and
      // the getLength() call.)
      double lngGap = max(0.0, PI - bound.lng().getLength() - 2.5 * DBL_EPSILON);

      // "minAbsLat" is the minimum distance from B to the equator (if zero or negative, then B
      // straddles the equator).
      // TODO(user): Unit tests do not currently check that this is correct.
      double minAbsLat = max(bound.lat().lo(), -bound.lat().hi());

      // "latGap1" and "latGap2" measure the minimum distance from B to the south and north poles
      // respectively.
      double latGap1 = M_PI_2 + bound.lat().lo();
      double latGap2 = M_PI_2 - bound.lat().hi();

      if (minAbsLat >= 0) {
        // The bound B does not straddle the equator. In this case the minimum distance is between
        // one endpoint of the latitude edge in B closest to the equator and the other endpoint of
        // that edge in B'. The latitude distance between these two points is 2*minAbsLat, and the
        // longitude distance is lngGap. We could compute the distance exactly using the Haversine
        // formula, but then we would need to bound the errors in that calculation. Since we only
        // need accuracy when the distance is very small (close to 4.309 * DBL_EPSILON), we
        // substitute the Euclidean distance instead. This gives us a right triangle XYZ with two
        // edges of length x = 2*minAbsLat and y ~= lngGap. The desired distance is the
        // length of the third edge "z", and we have
        //
        //         z  ~=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
        //
        // Therefore the region may contain nearly antipodal points only if
        //
        //  2*minAbsLat + lngGap  <  sqrt(2) * 4.309 * DBL_EPSILON
        //                           ~= 1.354e-15
        //
        // Note that because the given bound B is conservative, "minAbsLat" and "lngGap" are both
        // lower bounds on their true values so we do not need to make any adjustments for their
        // errors.
        if (2 * minAbsLat + lngGap < 1.354e-15) {
          return S2LatLngRect.full();
        }
      } else if (lngGap >= M_PI_2) {
        // B spans at most Pi/2 in longitude. The minimum distance is always between one corner of B
        // and the diagonally opposite corner of B'. We use the same distance approximation that we
        // used above; in this case we have an obtuse triangle XYZ with two edges of length x =
        // latGap1 and y = latGap2, and angle Z >= Pi/2 between them. We then have
        //
        //         z  >=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
        //
        // Unlike the case above, "latGap1" and "latGap2" are not lower bounds (because of the extra
        // addition operation, and because M_PI_2 is not exactly equal to Pi/2); they can exceed
        // their true values by up to 0.75 * DBL_EPSILON. Putting this all together, the region may
        // contain nearly antipodal points only if
        //
        //   latGap1 + latGap2  <  (sqrt(2) * 4.309 + 1.5) * DBL_EPSILON
        //                        ~= 1.687e-15
        if (latGap1 + latGap2 < 1.687e-15) {
          return S2LatLngRect.full();
        }
      } else {
        // Otherwise we know that (1) the bound straddles the equator and (2) its width in longitude
        // is at least Pi/2. In this case the minimum distance can occur either between a corner of
        // B and the diagonally opposite corner of B' (as in the case above), or between a corner of
        // B and the opposite longitudinal edge reflected in B'. It is sufficient to only consider
        // the corner-edge case, since this distance is also a lower bound on the corner-corner
        // distance when that case applies.

        // Consider the spherical triangle XYZ where X is a corner of B with minimum absolute
        // latitude, Y is the closest pole to X, and Z is the point closest to X on the opposite
        // longitudinal edge of B'. This is a right triangle (Z = Pi/2), and from the spherical law
        // of sines we have
        //
        //     sin(z) / sin(Z)  =  sin(y) / sin(Y)
        //     sin(maxLatGap) / 1  =  sin(d_min) / sin(lngGap)
        //     sin(d_min)  =  sin(maxLatGap) * sin(lngGap)
        //
        // where "maxLatGap" = max(latGap1, latGap2) and "d_min" is the desired minimum distance.
        // Now using the facts that sin(t) >= (2/Pi)*t for 0 <= t <= Pi/2, that we only need an
        // accurate approximation when at least one of "maxLatGap" or "lngGap" is extremely small
        // (in which case sin(t) ~= t), and recalling that "maxLatGap" has an error of up to
        // 0.75 * DBL_EPSILON, we want to test whether
        //
        //   maxLatGap * lngGap  <  (4.309 + 0.75) * (Pi/2) * DBL_EPSILON
        //                          ~= 1.765e-15
        if (max(latGap1, latGap2) * lngGap < 1.765e-15) {
          return S2LatLngRect.full();
        }
      }

      // Next we need to check whether the subregion might contain any edges that span (M_PI - 2 *
      // DBL_EPSILON) radians or more in longitude, since AddPoint sets the longitude bound to
      // Full() in that case. This corresponds to testing whether (lngGap <= 0) in "lngExpansion"
      // below.

      // Otherwise, the maximum latitude error in AddPoint is 4.8 * DBL_EPSILON. In the worst case,
      // the errors when computing the latitude bound for a subregion could go in the opposite
      // direction as the errors when computing the bound for the original region, so we need to
      // double this value. (More analysis shows that it's okay to round down to a multiple of
      // DBL_EPSILON.)
      //
      // For longitude, we rely on the fact that atan2 is correctly rounded and therefore no
      // additional bounds expansion is necessary.

      double latExpansion = 9 * DBL_EPSILON;
      double lngExpansion = (lngGap <= 0) ? PI : 0;
      return bound.expanded(S2LatLng.fromRadians(latExpansion, lngExpansion)).polarClosure();
    }
  }

  /**
   * The purpose of this class is to find edges that intersect a given XYZ bounding box. It can be
   * used as an efficient rejection test when attempting to find edges that intersect a given
   * region. It accepts a vertex chain v0, v1, v2, ... and returns a boolean value indicating
   * whether each edge intersects the specified bounding box.
   *
   * <p>We use XYZ intervals instead of something like longitude intervals because it is cheap to
   * collect from S2Point lists and any slicing strategy should give essentially equivalent results.
   * See S2Loop for an example of use.
   */
  public static class XYZPruner {
    private S2Point lastVertex;

    // The region to be tested against.
    private boolean boundSet;
    private double xmin;
    private double ymin;
    private double zmin;
    private double xmax;
    private double ymax;
    private double zmax;
    private double maxDeformation;

    public XYZPruner() {
      boundSet = false;
    }

    /**
     * Accumulate a bounding rectangle from provided edges.
     *
     * @param from start of edge
     * @param to end of edge.
     */
    public void addEdgeToBounds(S2Point from, S2Point to) {
      if (!boundSet) {
        boundSet = true;
        xmin = xmax = from.x;
        ymin = ymax = from.y;
        zmin = zmax = from.z;
      }
      xmin = min(xmin, min(to.x, from.x));
      ymin = min(ymin, min(to.y, from.y));
      zmin = min(zmin, min(to.z, from.z));
      xmax = max(xmax, max(to.x, from.x));
      ymax = max(ymax, max(to.y, from.y));
      zmax = max(zmax, max(to.z, from.z));

      // Because our arcs are really geodesics on the surface of the earth an edge can have
      // intermediate points outside the xyz bounds implicit in the end points. Based on the length
      // of the arc we compute a generous bound for the maximum amount of deformation. For small
      // edges it will be very small but for some large arcs (ie. from (1N,90W) to (1N,90E) the path
      // can be wildly deformed. I did a bunch of experiments with geodesics to get safe bounds for
      // the deformation.
      double approxArcLen = abs(from.x - to.x) + abs(from.y - to.y) + abs(from.z - to.z);
      if (approxArcLen < 0.025) { // less than 2 degrees
        maxDeformation = max(maxDeformation, approxArcLen * 0.0025);
      } else if (approxArcLen < 1.0) { // less than 90 degrees
        maxDeformation = max(maxDeformation, approxArcLen * 0.11);
      } else {
        maxDeformation = approxArcLen * 0.5;
      }
    }

    public void setFirstIntersectPoint(S2Point v0) {
      xmin = xmin - maxDeformation;
      ymin = ymin - maxDeformation;
      zmin = zmin - maxDeformation;
      xmax = xmax + maxDeformation;
      ymax = ymax + maxDeformation;
      zmax = zmax + maxDeformation;
      this.lastVertex = v0;
    }

    /**
     * Returns true if the edge going from the last point to this point passes through the pruner
     * bounding box, otherwise returns false. So the method returns false if we are certain there is
     * no intersection, but it may return true when there turns out to be no intersection.
     */
    public boolean intersects(S2Point v1) {
      boolean result = true;

      if ((v1.x < xmin && lastVertex.x < xmin) || (v1.x > xmax && lastVertex.x > xmax)) {
        result = false;
      } else if ((v1.y < ymin && lastVertex.y < ymin) || (v1.y > ymax && lastVertex.y > ymax)) {
        result = false;
      } else if ((v1.z < zmin && lastVertex.z < zmin) || (v1.z > zmax && lastVertex.z > zmax)) {
        result = false;
      }

      lastVertex = v1;
      return result;
    }
  }

  /**
   * The purpose of this class is to find edges that intersect a given longitude interval. It can be
   * used as an efficient rejection test when attempting to find edges that intersect a given
   * region. It accepts a vertex chain v0, v1, v2, ... and returns a boolean value indicating
   * whether each edge intersects the specified longitude interval.
   *
   * <p>This class is not currently used as the XYZPruner is preferred for S2Loop, but this should
   * be usable in similar circumstances. Be wary of the cost of atan2() in conversions from S2Point
   * to longitude!
   */
  public static class LongitudePruner {
    // The interval to be tested against.
    private final S1Interval interval;

    // The longitude of the next v0.
    private double lng0;

    /**
     * 'interval' is the longitude interval to be tested against, and 'v0' is the first vertex of
     * edge chain.
     */
    public LongitudePruner(S1Interval interval, S2Point v0) {
      this.interval = interval;
      this.lng0 = S2LatLng.longitude(v0).radians();
    }

    /**
     * Returns true if the edge (v0, v1) intersects the given longitude interval, and then saves
     * 'v1' to be used as the next 'v0'.
     */
    public boolean intersects(S2Point v1) {
      double lng1 = S2LatLng.longitude(v1).radians();
      boolean result = interval.intersects(S1Interval.fromPointPair(lng0, lng1));
      lng0 = lng1;
      return result;
    }
  }

  /** Spatial containment relationships between a wedge A to another wedge B. */
  enum WedgeRelation {
    /** A and B are equal. */
    WEDGE_EQUALS,
    /** A is a strict superset of B. */
    WEDGE_PROPERLY_CONTAINS,
    /** A is a strict subset of B. */
    WEDGE_IS_PROPERLY_CONTAINED,
    /** A-B, B-A, and A intersect B are non-empty. */
    WEDGE_PROPERLY_OVERLAPS,
    /** A and B are disjoint. */
    WEDGE_IS_DISJOINT,
  }

  /** Returns the relation from wedge A to B. */
  public static WedgeRelation getWedgeRelation(
      S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
    // There are 6 possible edge orderings at a shared vertex (all of these orderings are circular,
    // i.e. abcd == bcda):
    //
    //  (1) a2 b2 b0 a0: A contains B
    //  (2) a2 a0 b0 b2: B contains A
    //  (3) a2 a0 b2 b0: A and B are disjoint
    //  (4) a2 b0 a0 b2: A and B intersect in one wedge
    //  (5) a2 b2 a0 b0: A and B intersect in one wedge
    //  (6) a2 b0 b2 a0: A and B intersect in two wedges
    //
    // We do not distinguish between 4, 5, and 6. We pay extra attention when some of the edges
    // overlap. When edges overlap, several of these orderings can be satisfied, and we take the
    // most specific.
    if (a0.equalsPoint(b0) && a2.equalsPoint(b2)) {
      return WedgeRelation.WEDGE_EQUALS;
    }

    if (orderedCCW(a0, a2, b2, ab1)) {
      // The cases with this vertex ordering are 1, 5, and 6, although case 2 is also possible if
      // a2 == b2.
      if (orderedCCW(b2, b0, a0, ab1)) {
        return WedgeRelation.WEDGE_PROPERLY_CONTAINS;
      }

      // We are in case 5 or 6, or case 2 if a2 == b2.
      if (a2.equalsPoint(b2)) {
        return WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED;
      } else {
        return WedgeRelation.WEDGE_PROPERLY_OVERLAPS;
      }
    }

    // We are in case 2, 3, or 4.
    if (orderedCCW(a0, b0, b2, ab1)) {
      return WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED;
    }
    if (orderedCCW(a0, b0, a2, ab1)) {
      return WedgeRelation.WEDGE_IS_DISJOINT;
    } else {
      return WedgeRelation.WEDGE_PROPERLY_OVERLAPS;
    }
  }

  /**
   * Wedge processors are used to determine the local relationship between two polygons that share a
   * common vertex.
   *
   * <p>Given an edge chain (x0, x1, x2), the wedge at x1 is the region to the left of the edges.
   * More precisely, it is the set of all rays from x1x0 (inclusive) to x1x2 (exclusive) in the
   * *clockwise* direction.
   *
   * <p>Implementations compare two *non-empty* wedges that share the same middle vertex: A=(a0,
   * ab1, a2) and B=(b0, ab1, b2).
   *
   * <p>All wedge processors require that a0 != a2 and b0 != b2. Other degenerate cases (such as a0
   * == b2) are handled as expected. The parameter "ab1" denotes the common vertex a1 == b1.
   */
  public interface WedgeProcessor {
    /**
     * A wedge processor's test method accepts two edge chains A=(a0,a1,a2) and B=(b0,b1,b2) where
     * a1==b1, and returns either -1, 0, or 1 to indicate the relationship between the region to the
     * left of A and the region to the left of B.
     */
    int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2);
  }

  /**
   * Returns true if wedge A contains wedge B. Equivalent to but faster than {@code
   * getWedgeRelation() == WEDGE_PROPERLY_CONTAINS || WEDGE_EQUALS}.
   */
  public static class WedgeContains implements WedgeProcessor {
    /**
     * Given two edge chains, this function returns +1 if the region to the left of A contains the
     * region to the left of B, and 0 otherwise.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // For A to contain B (where each wedge interior is defined to be its left side), the CCW edge
      // order around ab1 must be a2 b2 b0 a0. We split this test into two parts that test three
      // vertices each.
      return orderedCCW(a2, b2, b0, ab1) && orderedCCW(b0, a0, a2, ab1) ? 1 : 0;
    }
  }

  /**
   * Returns true if wedge A intersects wedge B. Equivalent to but faster than {@code
   * getWedgeRelation() != WEDGE_IS_DISJOINT}.
   */
  public static class WedgeIntersects implements WedgeProcessor {
    /**
     * Given two edge chains (see WedgeRelation above), this function returns -1 if the region to
     * the left of A intersects the region to the left of B, and 0 otherwise. Note that regions are
     * defined such that points along a boundary are contained by one side or the other, not both.
     * So for example, if A,B,C are distinct points ordered CCW around a vertex O, then the wedges
     * BOA, AOC, and COB do not intersect.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // For A not to intersect B, (where each wedge interior is defined to be its left side), the
      // CCW edge order around ab1 must be a0 b2 b0 a2. Note that it's important to check these
      // conditions as negatives, i.e. Not orderedCCW(a0, b2, b0, ab1) rather than
      // orderedCCW(b0, b2, a0, ab1) to get correct results when two vertices are the same.
      return (orderedCCW(a0, b2, b0, ab1) && orderedCCW(b0, a2, a0, ab1) ? 0 : -1);
    }
  }

  /** Determines if wedge A contains or intersects wedge B. */
  public static class WedgeContainsOrIntersects implements WedgeProcessor {
    /**
     * Given two edge chains (see WedgeRelation above), this function returns +1 if A contains B, 0
     * if A and B are disjoint, and -1 if A intersects but does not contain B.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // This is similar to WedgeContainsOrCrosses, except that we want to distinguish cases
      //  * (1) [A contains B],
      //  * (3) [A and B are disjoint], and
      //  * (2,4,5,6) [A intersects but does not contain B].

      if (orderedCCW(a0, a2, b2, ab1)) {
        // We are in case 1, 5, or 6, or case 2 if a2 == b2.
        return orderedCCW(b2, b0, a0, ab1) ? 1 : -1; // Case 1 vs. 2,5,6.
      }
      // We are in cases 2, 3, or 4.
      if (!orderedCCW(a2, b0, b2, ab1)) {
        return 0; // Case 3.
      }

      // We are in case 2 or 4, or case 3 if a2 == b0.
      return a2.equalsPoint(b0) ? 0 : -1; // Case 3 vs. 2,4.
    }
  }

  /** Determines if the crossing or containment relationship between two edge chains. */
  public static class WedgeContainsOrCrosses implements WedgeProcessor {
    /**
     * Given two edge chains (see WedgeRelation above), this function returns +1 if A contains B, 0
     * if B contains A or the two wedges do not intersect, and -1 if the edge chains A and B cross
     * each other (i.e. if A intersects both the interior and exterior of the region to the left of
     * B). In degenerate cases where more than one of these conditions is satisfied, the maximum
     * possible result is returned. For example, if A == B then the result is +1.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // There are 6 possible edge orderings at a shared vertex (all of these orderings are
      // circular, i.e. abcd == bcda):
      //
      // (1) a2 b2 b0 a0: A contains B
      // (2) a2 a0 b0 b2: B contains A
      // (3) a2 a0 b2 b0: A and B are disjoint
      // (4) a2 b0 a0 b2: A and B intersect in one wedge
      // (5) a2 b2 a0 b0: A and B intersect in one wedge
      // (6) a2 b0 b2 a0: A and B intersect in two wedges
      //
      // In cases (4-6), the boundaries of A and B cross (i.e. the boundary of A intersects the
      // interior and exterior of B and vice versa). Thus we want to distinguish cases (1), (2-3),
      // and (4-6).
      //
      // Note that the vertices may satisfy more than one of the edge orderings above if two or
      // more vertices are the same. The tests below are written so that we take the most favorable
      // interpretation, i.e. preferring (1) over (2-3) over (4-6). In particular note that if
      // orderedCCW(a,b,c,o) returns true, it may be possible that orderedCCW(c,b,a,o) is also true
      // (if a == b or b == c).

      if (orderedCCW(a0, a2, b2, ab1)) {
        // The cases with this vertex ordering are 1, 5, and 6, although case 2 is also possible if
        // a2 == b2.
        if (orderedCCW(b2, b0, a0, ab1)) {
          return 1; // Case 1 (A contains B)
        }

        // We are in case 5 or 6, or case 2 if a2 == b2.
        return a2.equalsPoint(b2) ? 0 : -1; // Case 2 vs. 5,6.
      }
      // We are in case 2, 3, or 4.
      return orderedCCW(a0, b0, a2, ab1) ? 0 : -1; // Case 2,3 vs. 4.
    }
  }

  /**
   * FaceSegment represents an edge AB clipped to an S2 cube face. It is represented by a face index
   * and a pair of (u,v) coordinates.
   */
  static final class FaceSegment {
    int face;
    final R2Vector a = new R2Vector();
    final R2Vector b = new R2Vector();

    /** Returns an array of newly created FaceSegments. */
    static FaceSegment[] allFaces() {
      FaceSegment[] faces = new FaceSegment[6];
      for (int i = 0; i < faces.length; i++) {
        faces[i] = new FaceSegment();
      }
      return faces;
    }
  }

  // The three functions below all compare a sum (u + v) to a third value w. They are implemented in
  // such a way that they produce an exact result even though all calculations are done with
  // ordinary floating-point operations. Here are the principles on which these functions are based:
  //
  // A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
  //
  // B. If u + v < w in exact arithmetic, then at least one of the following expressions is true in
  //    floating-point:
  //      u + v < w
  //      u < w - v
  //      v < w - u
  //
  //   Proof: By rearranging terms and substituting ">" for "<", we can assume that all values are
  //   non-negative. Now clearly "w" is not the smallest value, so assume that "u" is the smallest.
  //   We want to show that u < w - v in floating-point. If v >= w/2, the calculation of w - v is
  //   exact since the result is smaller in magnitude than either input value, so the result holds.
  //   Otherwise we have u <= v < w/2 and w - v >= w/2 (even in floating point), so the result also
  //   holds.

  /** Returns true if u + v == w exactly. */
  static boolean sumEquals(double u, double v, double w) {
    return (u + v == w) && (u == w - v) && (v == w - u);
  }

  /**
   * Returns true if a given directed line L intersects the cube face F. The line L is defined by
   * its normal N in the (u,v,w) coordinates of F.
   */
  static boolean intersectsFace(S2Point n) {
    // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot products of N with the
    // four corner vertices (-1,-1,1), (1,-1,1), (1,1,1), and (-1,1,1) do not all have the same
    // sign. This is true exactly when |Nu| + |Nv| >= |Nw|. The code below evaluates this expression
    // exactly (see comments above).
    double u = abs(n.x);
    double v = abs(n.y);
    double w = abs(n.z);
    // We only need to consider the cases where u or v is the smallest value, since if w is the
    // smallest then both expressions below will have a positive LHS and a negative RHS.
    return (v >= w - u) && (u >= w - v);
  }

  /**
   * Given a directed line L intersecting a cube face F, return true if L intersects two opposite
   * edges of F (including the case where L passes exactly through a corner vertex of F). The line L
   * is defined by its normal N in the (u,v,w) coordinates of F.
   */
  static boolean intersectsOppositeEdges(S2Point n) {
    // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if and only exactly
    // two of the corner vertices lie on each side of L. This is true exactly when
    // ||Nu| - |Nv|| >= |Nw|. The code below evaluates this expression exactly (see comments above).
    double u = abs(n.x);
    double v = abs(n.y);
    double w = abs(n.z);
    // If w is the smallest, the following line returns an exact result.
    if (abs(u - v) != w) {
      return abs(u - v) >= w;
    }
    // Otherwise u - v = w exactly, or w is not the smallest value. In either case the following
    // line returns the correct result.
    return (u >= v) ? (u - w >= v) : (v - w >= u);
  }

  /**
   * Given cube face F and a directed line L (represented by its CCW normal N in the (u,v,w)
   * coordinates of F), compute the axis of the cube face edge where L exits the face: return 0 if L
   * exits through the u=-1 or u=+1 edge, and 1 if L exits through the v=-1 or v=+1 edge. Either
   * result is acceptable if L exits exactly through a corner vertex of the cube face.
   */
  static int getExitAxis(S2Point n) {
    assert intersectsFace(n);
    if (intersectsOppositeEdges(n)) {
      // The line passes through opposite edges of the face. It exits through the v=+1 or v=-1 edge
      // if the u-component of N has a larger absolute magnitude than the v-component.
      return (abs(n.x) >= abs(n.y)) ? 1 : 0;
    } else {
      // The line passes through two adjacent edges of the face. It exits the v=+1 or v=-1 edge if
      // an even number of the components of N are negative. We test this using signbit() rather
      // than multiplication to avoid the possibility of underflow.
      assert n.x != 0 && n.y != 0 && n.z != 0;
      return ((n.x < 0) ^ (n.y < 0) ^ (n.z < 0)) ? 0 : 1;
    }
  }

  /**
   * Given a cube face F, a directed line L (represented by its CCW normal N in the (u,v,w)
   * coordinates of F), and result of {@link #getExitAxis(S2Point)}, set {@code result} to the (u,v)
   * coordinates of the point where L exits the cube face.
   */
  static void getExitPoint(S2Point n, int axis, R2Vector result) {
    if (axis == 0) {
      result.x = (n.y > 0) ? 1.0 : -1.0;
      result.y = (-result.x * n.x - n.z) / n.y;
    } else {
      result.y = (n.x < 0) ? 1.0 : -1.0;
      result.x = (-result.y * n.y - n.z) / n.x;
    }
  }

  /**
   * Given a line segment AB whose origin A has been projected onto a given cube face, determine
   * whether it is necessary to project A onto a different face instead. This can happen because the
   * normal of the line AB is not computed exactly, so that the line AB (defined as the set of
   * points perpendicular to the normal) may not intersect the cube face containing A. Even if it
   * does intersect the face, the "exit point" of the line from that face may be on the wrong side
   * of A (i.e., in the direction away from B). If this happens, we reproject A onto the adjacent
   * face where the line AB approaches A most closely. This moves the origin by a small amount, but
   * never more than the error tolerances documented in the header file.
   */
  static int moveOriginToValidFace(int face, S2Point a, S2Point ab, R2Vector aUv) {
    // Fast path: if the origin is sufficiently far inside the face, it is always safe to use it.
    final double kMaxSafeUVCoord = 1 - FACE_CLIP_ERROR_UV_COORD;
    double au = aUv.x;
    double av = aUv.y;
    if (max(abs(au), abs(av)) <= kMaxSafeUVCoord) {
      return face;
    }

    // Otherwise check whether the normal AB even intersects this face.
    S2Point n = S2Projections.faceXyzToUvw(face, ab);
    if (intersectsFace(n)) {
      // Check whether the point where the line AB exits this face is on the wrong side of A (by
      // more than the acceptable error tolerance).
      getExitPoint(n, getExitAxis(n), aUv);
      S2Point exit = S2Projections.faceUvToXyz(face, aUv);
      S2Point aTangent = ab.normalize().crossProd(a);
      if (exit.sub(a).dotProd(aTangent) >= -FACE_CLIP_ERROR_RADIANS) {
        // We can use the given face, but first put the original values back.
        aUv.x = au;
        aUv.y = av;
        return face;
      }
    }

    // Otherwise we reproject A to the nearest adjacent face. (If line AB does not pass through a
    // given face, it must pass through all adjacent faces.)
    if (abs(au) >= abs(av)) {
      face = S2Projections.getUVWFace(face, 0 /*U axis*/, au > 0 ? 1 : 0);
    } else {
      face = S2Projections.getUVWFace(face, 1 /*V axis*/, av > 0 ? 1 : 0);
    }
    assert intersectsFace(S2Projections.faceXyzToUvw(face, ab));
    S2Projections.validFaceXyzToUv(face, a, aUv);
    aUv.set(max(-1.0, min(1.0, aUv.x)), max(-1.0, min(1.0, aUv.y)));
    return face;
  }

  /**
   * Return the next face that should be visited by getFaceSegments, given that we have just visited
   * "face" and we are following the line AB (represented by its normal N in the (u,v,w) coordinates
   * of that face). The other arguments include the point where AB exits "face", the corresponding
   * exit axis, and the "target face" containing the destination point B.
   */
  static int getNextFace(int face, R2Vector exit, int axis, S2Point n, int targetFace) {
    // We return the face that is adjacent to the exit point along the given axis. If line AB exits
    // *exactly* through a corner of the face, there are two possible next faces. If one is the
    // "target face" containing B, then we guarantee that we advance to that face directly.
    //
    // The three conditions below check that (1) AB exits approximately through a corner, (2) the
    // adjacent face along the non-exit axis is the target face, and (3) AB exits *exactly* through
    // the corner. (The sumEquals() code checks whether the dot product of (u,v,1) and "n" is
    // exactly zero.)
    if (abs(exit.get(1 - axis)) == 1
        && S2Projections.getUVWFace(face, 1 - axis, exit.get(1 - axis) > 0 ? 1 : 0) == targetFace
        && sumEquals(exit.x * n.x, exit.y * n.y, -n.z)) {
      return targetFace;
    }

    // Otherwise return the face that is adjacent to the exit point in the direction of the exit
    // axis.
    return S2Projections.getUVWFace(face, axis, exit.get(axis) > 0 ? 1 : 0);
  }

  /**
   * Subdivide the given edge AB at every point where it crosses the boundary between two S2 cube
   * faces, returning the number of FaceSegments entries used (all entries must be prefilled). The
   * segments are returned in order from A toward B. The input points must be unit length.
   *
   * <p>This method guarantees that the returned segments form a continuous path from A to B, and
   * that all vertices are within FACE_CLIP_ERROR_UV_DIST of the line AB. All vertices lie within
   * the [-1,1]x[-1,1] cube face rectangles. The results are consistent with
   * {@link S2Predicates.Sign#expensive(S2Point, S2Point, S2Point, boolean)}, i.e. the edge is well-
   * defined even if its endpoints are antipodal.
   */
  static int getFaceSegments(S2Point a, S2Point b, FaceSegment[] segments) {
    assert isUnitLength(a);
    assert isUnitLength(b);

    // Fast path: both endpoints are on the same face.
    FaceSegment seg = segments[0];
    seg.face = S2Projections.xyzToFace(a);
    S2Projections.validFaceXyzToUv(seg.face, a, seg.a);
    int bFace = S2Projections.xyzToFace(b);
    S2Projections.validFaceXyzToUv(bFace, b, seg.b);
    if (seg.face == bFace) {
      return 1;
    } else {
      // Starting at A, we follow AB from face to face until we reach the face containing B. The
      // following code is designed to ensure that we always reach B, even in the presence of
      // numerical errors.
      //
      // First we compute the normal to the plane containing A and B. This normal becomes the
      // ultimate definition of the line AB; it is used to resolve all questions regarding where
      // exactly the line goes. Unfortunately due to numerical errors, the line may not quite
      // intersect the faces containing the original endpoints. We handle this by moving A and/or
      // B slightly if necessary so that they are on faces intersected by the line AB.
      S2Point ab = robustCrossProd(a, b);
      seg.face = moveOriginToValidFace(seg.face, a, ab, seg.a);
      bFace = moveOriginToValidFace(bFace, b, ab.neg(), seg.b);

      // Save b in the last possible segment.
      segments[5].b.set(seg.b);

      // Now we simply follow AB from face to face until we reach B.
      int size = 1;
      while (seg.face != bFace) {
        // Complete the current segment by finding the point where AB exits the current face.
        S2Point n = S2Projections.faceXyzToUvw(seg.face, ab);
        int exitAxis = getExitAxis(n);
        getExitPoint(n, exitAxis, seg.b);

        // Compute the next face intersected by AB, and translate the exit point of the current
        // segment into the (u,v) coordinates of the next face. This becomes the first point of the
        // next segment.
        int newFace = getNextFace(seg.face, seg.b, exitAxis, n, bFace);
        S2Point oldExitXyz = S2Projections.faceUvToXyz(seg.face, seg.b);
        S2Point newExitUvw = S2Projections.faceXyzToUvw(newFace, oldExitXyz);

        // Set up the first half of the next segment.
        seg = segments[size++];
        seg.face = newFace;
        seg.a.set(newExitUvw.x, newExitUvw.y);
      }

      // Finish the last segment.
      seg.b.set(segments[5].b);
      return size;
    }
  }

  /**
   * This helper function does two things. First, it clips the line segment AB to find the clipped
   * destination B' on a given face. (The face is specified implicitly by expressing *all arguments*
   * in the (u,v,w) coordinates of that face.) Second, it partially computes whether the segment AB
   * intersects this face at all. The actual condition is fairly complicated, but it turns out that
   * it can be expressed as a "score" that can be computed independently when clipping the two
   * endpoints A and B. This function returns the score for the given endpoint, which is an integer
   * ranging from 0 to 3. If the sum of the two scores is 3 or more, then AB does not intersect this
   * face. See the calling function for the meaning of the various parameters.
   */
  static int clipDestination(
      S2Point a,
      S2Point b,
      S2Point nScaled,
      S2Point aTangent,
      S2Point bTangent,
      double uvScale,
      R2Vector uv) {
    assert intersectsFace(nScaled);

    // Optimization: if B is within the safe region of the face, use it.
    final double kMaxSafeUVCoord = 1 - FACE_CLIP_ERROR_UV_COORD;
    if (b.z > 0) {
      uv.set(b.x / b.z, b.y / b.z);
      if (max(abs(uv.x), abs(uv.y)) <= kMaxSafeUVCoord) {
        return 0;
      }
    }

    // Otherwise find the point B' where the line AB exits the face.
    getExitPoint(nScaled, getExitAxis(nScaled), uv);
    uv.x *= uvScale;
    uv.y *= uvScale;
    S2Point p = new S2Point(uv.x, uv.y, 1);

    // Determine if the exit point B' is contained within the segment. We do this by computing the
    // dot products with two inward-facing tangent vectors at A and B. If either dot product is
    // negative, we say that B' is on the "wrong side" of that point. As the point B' moves around
    // the great circle AB past the segment endpoint B, it is initially on the wrong side of B only;
    // as it moves further it is on the wrong side of both endpoints; and then it is on the wrong
    // side of A only. If the exit point B' is on the wrong side of either endpoint, we can't use
    // it; instead the segment is clipped at the original endpoint B.
    //
    // We reject the segment if the sum of the scores of the two endpoints is 3 or more. Here is
    // what that rule encodes:
    //  - If B' is on the wrong side of A, then the other clipped endpoint A' must be in the
    //    interior of AB (otherwise AB' would go the wrong way around the circle). There is a
    //    similar rule for A'.
    //  - If B' is on the wrong side of either endpoint (and therefore we must use the original
    //    endpoint B instead), then it must be possible to project B onto this face (i.e., its
    //    w-coordinate must be positive).
    // This rule is only necessary to handle certain zero-length edges (A=B).
    int score = 0;
    if (p.sub(a).dotProd(aTangent) < 0) {
      score = 2; // B' is on wrong side of A.
    } else if (p.sub(b).dotProd(bTangent) < 0) {
      score = 1; // B' is on wrong side of B.
    }
    if (score > 0) { // B' is not in the interior of AB.
      if (b.z <= 0) {
        score = 3; // B cannot be projected onto this face.
      } else {
        uv.set(b.x / b.z, b.y / b.z);
      }
    }
    return score;
  }

  /**
   * As {@link #clipToFace(S2Point, S2Point, int, R2Vector, R2Vector)}, but rather than clipping to
   * the square [-1,1]x[-1,1] in (u,v) space, this method clips to [-R,R]x[-R,R] where
   * R=(1+padding).
   *
   * <p>If the given edge does not intersect the face, returns false and does NOT fill aUv and bUv.
   */
  public static boolean clipToPaddedFace(
      S2Point aXyz, S2Point bXyz, int face, double padding, R2Vector aUv, R2Vector bUv) {
    assert padding >= 0;
    // Fast path: both endpoints are on the given face.
    if (S2Projections.xyzToFace(aXyz) == face && S2Projections.xyzToFace(bXyz) == face) {
      S2Projections.validFaceXyzToUv(face, aXyz, aUv);
      S2Projections.validFaceXyzToUv(face, bXyz, bUv);
      return true;
    }

    // Convert everything into the (u,v,w) coordinates of the given face. Note that the cross
    // product *must* be computed in the original (x,y,z) coordinate system because RobustCrossProd
    // (unlike the mathematical cross product) can produce different results in different coordinate
    // systems when one argument is a linear multiple of the other, due to the use of symbolic
    // perturbations.
    S2Point n = S2Projections.faceXyzToUvw(face, robustCrossProd(aXyz, bXyz));
    S2Point a = S2Projections.faceXyzToUvw(face, aXyz);
    S2Point b = S2Projections.faceXyzToUvw(face, bXyz);

    // Padding is handled by scaling the u- and v-components of the normal. Letting R=1+padding,
    // this means that when we compute the dot product of the normal with a cube face vertex (such
    // as (-1,-1,1)), we will actually compute the dot product with the scaled vertex (-R,-R,1).
    // This allows methods such as intersectsFace(), getExitAxis(), etc, to handle padding with no
    // further modifications.
    final double uvScale = 1 + padding;
    S2Point nScaled = new S2Point(uvScale * n.x, uvScale * n.y, n.z);
    if (!intersectsFace(nScaled)) {
      return false;
    }

    n = n.normalize();
    S2Point aTangent = n.crossProd(a);
    S2Point bTangent = b.crossProd(n);
    // As described above, if the sum of the scores from clipping the two endpoints is 3 or more,
    // then the segment does not intersect this face.
    int aScore = clipDestination(b, a, nScaled.neg(), bTangent, aTangent, uvScale, aUv);
    int bScore = clipDestination(a, b, nScaled, aTangent, bTangent, uvScale, bUv);
    return aScore + bScore < 3;
  }

  /**
   * Returns true if the edge AB intersects the given (closed) rectangle to within the error bound
   * below.
   */
  static boolean intersectsRect(R2Vector a, R2Vector b, R2Rect rect) {
    // First check whether the bound of AB intersects "rect".
    R2Rect bound = R2Rect.fromPointPair(a, b);
    if (!rect.intersects(bound)) {
      return false;
    }

    // Otherwise AB intersects "rect" if and only if all four vertices of "rect" do not lie on the
    // same side of the extended line AB. We test this by finding the two vertices of "rect" with
    // minimum and maximum projections onto the normal of AB, and computing their dot products with
    // the edge normal.
    R2Vector n = R2Vector.sub(b, a).ortho();
    int i = (n.x >= 0) ? 1 : 0;
    int j = (n.y >= 0) ? 1 : 0;
    double max = n.dotProd(R2Vector.sub(rect.getVertex(i, j), a));
    double min = n.dotProd(R2Vector.sub(rect.getVertex(1 - i, 1 - j), a));
    return (max >= 0) && (min <= 0);
  }

  /** Moves an endpoint of the given bound to the given value. */
  static boolean updateEndpoint(R1Interval bound, boolean slopeNegative, double value) {
    if (!slopeNegative) {
      if (bound.hi() < value) {
        return false;
      }
      if (bound.lo() < value) {
        bound.setLo(value);
      }
    } else {
      if (bound.lo() > value) {
        return false;
      }
      if (bound.hi() > value) {
        bound.setHi(value);
      }
    }
    return true;
  }

  /**
   * Given a line segment from (a0,a1) to (b0,b1) and a bounding interval for each axis, clip the
   * segment further if necessary so that "bound0" does not extend outside the given interval
   * "clip". "diag" is a a precomputed helper variable that indicates which diagonal of the bounding
   * box is spanned by AB: it is 0 if AB has positive slope, and 1 if AB has negative slope.
   */
  static boolean clipBoundAxis(
      double a0,
      double b0,
      R1Interval bound0,
      double a1,
      double b1,
      R1Interval bound1,
      boolean slopeNegative,
      R1Interval clip0) {
    if (bound0.lo() < clip0.lo()) {
      if (bound0.hi() < clip0.lo()) {
        return false;
      }
      bound0.setLo(clip0.lo());
      if (!updateEndpoint(bound1, slopeNegative, interpolateDouble(clip0.lo(), a0, b0, a1, b1))) {
        return false;
      }
    }
    if (bound0.hi() > clip0.hi()) {
      if (bound0.lo() > clip0.hi()) {
        return false;
      }
      bound0.setHi(clip0.hi());
      if (!updateEndpoint(bound1, !slopeNegative, interpolateDouble(clip0.hi(), a0, b0, a1, b1))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Given an edge AB and a rectangle "clip", return the bounding rectangle of the portion of AB
   * intersected by "clip". The resulting bound may be empty. This is a convenience function built
   * on top of clipEdgeBound.
   */
  static R2Rect getClippedEdgeBound(R2Vector a, R2Vector b, R2Rect clip) {
    R2Rect bound = R2Rect.fromPointPair(a, b);
    if (clipEdgeBound(a, b, clip, bound)) {
      return bound;
    }
    return R2Rect.empty();
  }

  /**
   * This function can be used to clip an edge AB to sequence of rectangles efficiently. It
   * represents the clipped edges by their bounding boxes rather than as a pair of endpoints.
   * Specifically, let A'B' be some portion of an edge AB, and let "bound" be a tight bound of A'B'.
   * This function updates "bound" (in place) to be a tight bound of A'B' intersected with a given
   * rectangle "clip". If A'B' does not intersect "clip", returns false and does not necessarily
   * update "bound".
   *
   * <p>The given bound must be a tight bounding rectangle for some portion of AB. (This condition
   * is automatically satisfied if you start with the bounding box of AB and clip to a sequence of
   * rectangles, stopping when the method returns false.)
   */
  static boolean clipEdgeBound(R2Vector a, R2Vector b, R2Rect clip, R2Rect bound) {
    // "slopeNegative" indicates which diagonal of the bounding box is spanned by AB: it is false if
    // AB has positive slope, and true if AB has negative slope. This is used to determine which
    // interval endpoints need to be updated each time the edge is clipped.
    boolean slopeNegative = (a.x > b.x) != (a.y > b.y);
    return clipBoundAxis(a.x, b.x, bound.x(), a.y, b.y, bound.y(), slopeNegative, clip.x())
        && clipBoundAxis(a.y, b.y, bound.y(), a.x, b.x, bound.x(), slopeNegative, clip.y());
  }

  /**
   * Given an edge AB, assigns the portion of AB that is contained by the given rectangle "clip" to
   * the aClipped and bClipped arguments, and returns true if there is an intersection.
   */
  static boolean clipEdge(
      R2Vector a, R2Vector b, R2Rect clip, R2Vector aClipped, R2Vector bClipped) {
    // Compute the bounding rectangle of AB, clip it, and then extract the new endpoints from the
    // clipped bound.
    R2Rect bound = R2Rect.fromPointPair(a, b);
    if (clipEdgeBound(a, b, clip, bound)) {
      int iEnd = a.x > b.x ? 1 : 0;
      int jEnd = a.y > b.y ? 1 : 0;
      aClipped.set(bound.getVertex(iEnd, jEnd));
      bClipped.set(bound.getVertex(1 - iEnd, 1 - jEnd));
      return true;
    }
    return false;
  }

  /**
   * Given a value x that is some linear combination of a and b, return the value x1 that is the
   * same linear combination of a1 and b1. This function makes the following guarantees:
   *
   * <ol>
   *   <li>If x == a, then x1 = a1 (exactly).
   *   <li>If x == b, then x1 = b1 (exactly).
   *   <li>If If a <= x <= b and a1 <= b1, then a1 <= x1 <= b1 (even if a1 == b1).
   *   <li>More generally, if x is between a and b, then x1 is between a1 and b1.
   * </ol>
   *
   * <p>Results are undefined if a==b. When a <= x <= b or b <= x <= a we can prove the error
   * bound on the resulting value is 2.25*DBL_EPSILON. The error for extrapolating an x value
   * outside of a and b can be much worse. See the gappa proof in s2edge_clipping.h in the C++
   * implementation of S2.
   */
  static double interpolateDouble(double x, double a, double b, double a1, double b1) {
    assert a != b;
    // To get results that are accurate near both A and B, we interpolate starting from the closer
    // of the two points.
    if (abs(a - x) <= abs(b - x)) {
      return a1 + (b1 - a1) * ((x - a) / (b - a));
    } else {
      return b1 + (a1 - b1) * ((x - b) / (a - b));
    }
  }

  /**
   * Given an edge AB and a face, return the (u,v) coordinates in aUv and bUv for the portion of AB
   * that intersects that face. This method guarantees that the clipped vertices lie within the
   * [-1,1]x[-1,1] cube face rectangle and are within FACE_CLIP_ERROR_UV_DIST of the line AB, but
   * the results may differ from those produced by getFaceSegments.
   *
   * <p>Returns false and does not fill aUv and bUv if AB does not intersect the given face.
   */
  public static boolean clipToFace(S2Point a, S2Point b, int face, R2Vector aUv, R2Vector bUv) {
    return clipToPaddedFace(a, b, face, 0.0, aUv, bUv);
  }

  /**
   * Return true if edge AB crosses CD at a point that is interior to both edges. Properties:
   *
   * <ul>
   *   <li>simpleCrossing(b,a,c,d) == simpleCrossing(a,b,c,d)
   *   <li>simpleCrossing(c,d,a,b) == simpleCrossing(a,b,c,d)
   * </ul>
   */
  public static boolean simpleCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    // We compute simpleCCW() for triangles ACB, CBD, BDA, and DAC. All of these triangles need to
    // have the same orientation (CW or CCW) for an intersection to exist. Note that this is
    // slightly more restrictive than the corresponding definition for planar edges, since we need
    // to exclude pairs of line segments that would otherwise "intersect" by crossing two antipodal
    // points.

    S2Point ab = a.crossProd(b);
    double acb = -ab.dotProd(c);
    double bda = ab.dotProd(d);
    if (acb * bda <= 0) {
      return false;
    }

    S2Point cd = c.crossProd(d);
    double cbd = -cd.dotProd(b);
    if (acb * cbd <= 0) {
      return false;
    }

    double dac = cd.dotProd(a);
    return acb * dac > 0;
  }

  /**
   * This function determines whether the edge AB intersects the edge CD. Returns +1 if AB crosses
   * CD at a point that is interior to both edges. Returns 0 if any two vertices from different
   * edges are the same. Returns -1 otherwise.
   *
   * <p>Note that if an edge is degenerate (A == B or C == D), the return value is 0 if two vertices
   * from different edges are the same and -1 otherwise.
   *
   * <p>Properties of robustCrossing:
   *
   * <ul>
   *   <li>{@code robustCrossing(b,a,c,d) == robustCrossing(a,b,c,d)}
   *   <li>{@code robustCrossing(c,d,a,b) == robustCrossing(a,b,c,d)}
   *   <li>{@code robustCrossing(a,b,c,d) == 0} if {@code a==c, a==d, b==c, b==d}
   *   <li>{@code robustCrossing(a,b,c,d) <= 0} if {@code a==b or c==d} (see above)
   * </ul>
   *
   * <p>This function implements an exact, consistent perturbation model such that no three points
   * are ever considered to be collinear. This means that even if you have 4 points A, B, C, D that
   * lie exactly in a line (say, around the equator), C and D will be treated as being slightly to
   * one side or the other of AB. This is done in a way such that the results are always consistent
   * (see {@link S2Predicates.Sign}).
   *
   * <p>Note that if you want to check an edge against a *chain* of other edges, it is much more
   * efficient to use an EdgeCrosser (above).
   */
  public static int robustCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    EdgeCrosser crosser = new EdgeCrosser(a, b, c);
    return crosser.robustCrossing(d);
  }

  /**
   * Given two edges AB and CD where at least two vertices are identical (i.e.
   * robustCrossing(a,b,c,d) == 0), this function defines whether the two edges "cross" in such a
   * way that point-in-polygon containment tests can be implemented by counting the number of edge
   * crossings. The basic rule is that a "crossing" occurs if AB is encountered after CD during a
   * CCW sweep around the shared vertex starting from a fixed reference point.
   *
   * <p>Note that according to this rule, if AB crosses CD then in general CD does not cross AB.
   * However, this leads to the correct result when counting polygon edge crossings. For example,
   * suppose that A,B,C are three consecutive vertices of a CCW polygon. If we now consider the edge
   * crossings of a segment BP as P sweeps around B, the crossing number changes parity exactly when
   * BP crosses BA or BC.
   *
   * <p>Useful properties of vertexCrossing (VC):
   *
   * <ul>
   *   <li>VC(a,a,c,d) == VC(a,b,c,c) == false
   *   <li>VC(a,b,a,b) == VC(a,b,b,a) == true
   *   <li>VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
   *   <li>If exactly one of a,b equals one of c,d, then exactly one of VC(a,b,c,d) and VC(c,d,a,b)
   *       is true
   * </ul>
   *
   * <p>It is an error to call this method with 4 distinct vertices.
   */
  public static boolean vertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    // If A == B or C == D there is no intersection. We need to check this case first in case 3 or
    // more input points are identical.
    if (a.equalsPoint(b) || c.equalsPoint(d)) {
      return false;
    }

    // If any other pair of vertices is equal, there is a crossing if and only if orderedCCW()
    // indicates that the edge AB is further CCW around the shared vertex than the edge CD.
    if (a.equalsPoint(d)) {
      return orderedCCW(S2.refDir(a), c, b, a);
    }
    if (b.equalsPoint(c)) {
      return orderedCCW(S2.refDir(b), d, a, b);
    }
    if (a.equalsPoint(c)) {
      return orderedCCW(S2.refDir(a), d, b, a);
    }
    if (b.equalsPoint(d)) {
      return orderedCCW(S2.refDir(b), c, a, b);
    }

    assert false : "vertexCrossing called with 4 distinct vertices, which is not allowed";
    return false;
  }

  /**
   * Like {@link #vertexCrossing(S2Point, S2Point, S2Point, S2Point)} but returns -1 if AB crosses
   * CD from left to right, +1 if AB crosses CD from right to left, and 0 otherwise. This implies
   * that if CD bounds some region according to the "interior is on the left" rule, this function
   * returns -1 when AB exits the region and +1 when AB enters.
   *
   * <p>This is a helper method that allows computing the change in winding number from point A to
   * point B by summing the signed edge crossings of AB with the edges of the loop(s) used to define
   * the winding number.
   *
   * <p>Useful properties of signedVertexCrossing (SVC):
   *
   * <ul>
   *   <li>SVC(a,a,c,d) == SVC(a,b,c,c) == 0
   *   <li>SVC(a,b,a,b) == +1
   *   <li>SVC(a,b,b,a) == -1
   *   <li>SVC(a,b,c,d) == -SVC(a,b,d,c) == -SVC(b,a,c,d) == SVC(b,a,d,c)
   *   <li>If exactly one of a,b equals one of c,d, then exactly one of SVC(a,b,c,d) and
   *       SVC(c,d,a,b) is non-zero
   * </ul>
   *
   * <p>It is an error to call this method with 4 distinct vertices.
   */
  public static int signedVertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    if (a.equalsPoint(b) || c.equalsPoint(d)) {
      return 0;
    }

    // See vertexCrossing. The sign of the crossing is +1 if both edges are outgoing or both edges
    // are incoming with respect to the common vertex and -1 otherwise.
    if (a.equalsPoint(c)) {
      return (b.equalsPoint(d) || orderedCCW(S2.refDir(a), d, b, a)) ? 1 : 0;
    }
    if (b.equalsPoint(d)) {
      return orderedCCW(S2.refDir(b), c, a, b) ? 1 : 0;
    }
    if (a.equalsPoint(d)) {
      return (b.equalsPoint(c) || orderedCCW(S2.refDir(a), c, b, a)) ? -1 : 0;
    }
    if (b.equalsPoint(c)) {
      return orderedCCW(S2.refDir(b), d, a, b) ? -1 : 0;
    }

    assert false; // signedVertexCrossing called with 4 distinct vertices
    return 0;
  }

  /**
   * A convenience function that calls robustCrossing() to handle cases where all four vertices are
   * distinct, and VertexCrossing() to handle cases where two or more vertices are the same. This
   * defines a crossing function such that point-in-polygon containment tests can be implemented by
   * simply counting edge crossings.
   */
  public static boolean edgeOrVertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    int crossing = robustCrossing(a, b, c, d);
    if (crossing < 0) {
      return false;
    }
    if (crossing > 0) {
      return true;
    }
    return vertexCrossing(a, b, c, d);
  }

  /**
   * Finds the closest acceptable endpoint to a given point. An endpoint is acceptable if it lies
   * between the endpoints of the other line segment.
   */
  static S2Point closestAcceptableEndpoint(
      S2Point a0, S2Point a1, S2Point aNorm, S2Point b0, S2Point b1, S2Point bNorm, S2Point x) {
    CloserResult r = new CloserResult(Double.POSITIVE_INFINITY, x);
    if (orderedCCW(b0, a0, b1, bNorm)) {
      r.replaceIfCloser(x, a0);
    }
    if (orderedCCW(b0, a1, b1, bNorm)) {
      r.replaceIfCloser(x, a1);
    }
    if (orderedCCW(a0, b0, a1, aNorm)) {
      r.replaceIfCloser(x, b0);
    }
    if (orderedCCW(a0, b1, a1, aNorm)) {
      r.replaceIfCloser(x, b1);
    }
    return r.getVmin();
  }

  static class CloserResult {
    private double dmin2;
    private S2Point vmin;

    public double getDmin2() {
      return dmin2;
    }

    public S2Point getVmin() {
      return vmin;
    }

    public CloserResult(double dmin2, S2Point vmin) {
      this.dmin2 = dmin2;
      this.vmin = vmin;
    }

    public void replaceIfCloser(S2Point x, S2Point y) {
      // If the squared distance from x to y is less than dmin2, then replace vmin by y and update
      // dmin2 accordingly.
      double d2 = x.sub(y).norm2();
      if (d2 < dmin2 || (d2 == dmin2 && y.lessThan(vmin))) {
        dmin2 = d2;
        vmin = y;
      }
    }
  }

  /** Returns true if ab possibly crosses cd, by clipping tiny angles to zero. */
  public static final boolean lenientCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    assert isUnitLength(a);
    assert isUnitLength(b);
    assert isUnitLength(c);

    double acb = scalarTripleProduct(b, a, c);
    if (abs(acb) < MAX_DET_ERROR) {
      return true;
    }
    double bda = scalarTripleProduct(a, b, d);
    if (abs(bda) < MAX_DET_ERROR) {
      return true;
    }
    if (acb * bda < 0) {
      return false;
    }
    double cbd = scalarTripleProduct(d, c, b);
    if (abs(cbd) < MAX_DET_ERROR) {
      return true;
    }
    double dac = scalarTripleProduct(c, d, a);
    if (abs(dac) < MAX_DET_ERROR) {
      return true;
    }
    return (acb * cbd >= 0) && (acb * dac >= 0);
  }

  /**
   * Given two edges "ab" and "cd" such that robustCrossing() is true, return their intersection
   * point. Useful properties of getIntersection (GI):
   *
   * <ul>
   *   <li>GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
   *   <li>GI(c,d,a,b) == GI(a,b,c,d)
   * </ul>
   *
   * The returned intersection point X is guaranteed to be very close to the true intersection point
   * of "ab" and "cd", even if the edges intersect at a very small angle. See
   * {@link INTERSECTION_ERROR} above for details.
   */
  public static S2Point getIntersection(S2Point a, S2Point b, S2Point c, S2Point d) {
    return getIntersection(a, b, c, d, new ResultError());
  }

  /**
   * Helper for {@link #getIntersection(S2Point, S2Point, S2Point, S2Point)} with provided result
   * error parameter for testing and benchmarking purposes.
   */
  static S2Point getIntersection(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, ResultError resultError) {
    Preconditions.checkArgument(
        robustCrossing(a0, a1, b0, b1) > 0,
        "Input edges a0a1 and b0b1 must have a true robustCrossing.");

    // It is difficult to compute the intersection point of two edges accurately when the angle
    // between the edges is very small. Previously we handled this by only guaranteeing that the
    // returned intersection point is within INTERSECTION_ERROR of each edge. However, this means
    // that when the edges cross at a very small angle, the computed result may be very far from the
    // true intersection point.
    //
    // Instead this function now guarantees that the result is always within INTERSECTION_ERROR of
    // the true intersection. This requires using more sophisticated techniques and in some cases
    // extended precision.
    //
    // Two different techniques are used:
    //
    // <ul>
    //  <li>getIntersectionStable() computes the intersection point using projection and
    // interpolation, taking care to minimize cancellation error.
    //  <li>getIntersectionExact() computes the intersection point using exact arithmetic and
    // converts the final result back to an S2Point.
    // </ul>
    //
    // Our strategy is to first call getIntersectionStable(). If the result has an error bound
    // greater than INTERSECTION_ERROR, we fall back to exact arithmetic.
    S2Point result = getIntersectionApprox(a0, a1, b0, b1, resultError);
    if (resultError.error > INTERSECTION_ERROR) {
      result = getIntersectionExact(a0, a1, b0, b1);
    }
    return correctIntersectionSign(a0, a1, b0, b1, result);
  }

  /** Returns intersection result with sign corrected (if necessary). */
  static S2Point correctIntersectionSign(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, S2Point intersectionResult) {
    // Make sure the intersection point is on the correct side of the sphere. Since all vertices are
    // unit length, and both edge lengths are less than 180 degrees, (a0 + a1) and (b0 + b1) both
    // have positive dot product with the intersection point. We use the sum of all vertices to
    // make sure that the result is unchanged when the edges are swapped or reversed.
    if (intersectionResult.dotProd(a0.add(a1).add(b0.add(b1))) < 0) {
      intersectionResult = intersectionResult.neg();
    }
    return intersectionResult;
  }

  /**
   * Given a point X and an edge AB, return the distance ratio AX / (AX + BX). If X happens to be on
   * the line segment AB, this is the fraction "t" such that X == Interpolate(A, B, t). Requires
   * that A and B are distinct.
   */
  public static double getDistanceFraction(S2Point x, S2Point a, S2Point b) {
    Preconditions.checkArgument(!a.equalsPoint(b));
    double d0 = x.angle(a);
    double d1 = x.angle(b);
    return d0 / (d0 + d1);
  }

  /**
   * Return the minimum distance from X to any point on the edge AB. The result is very accurate for
   * small distances but may have some numerical error if the distance is large (approximately Pi/2
   * or greater). The case A == B is handled correctly.
   *
   * @throws IllegalArgumentException Thrown if the parameters are not all unit length.
   */
  public static S1Angle getDistance(S2Point x, S2Point a, S2Point b) {
    // TODO(user): Improve accuracy near Pi.
    Preconditions.checkArgument(isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(isUnitLength(b), "S2Point not normalized: %s", b);
    return S1Angle.radians(getDistanceRadians(x, a, b, robustCrossProd(a, b)));
  }

  /**
   * A slightly more efficient version of getDistance() where the cross product of the two endpoints
   * has been precomputed. The cross product does not need to be normalized, but should be computed
   * using S2RobustCrossProd.robustCrossProd() for the most accurate results.
   *
   * @throws IllegalArgumentException Thrown if the parameters are not all unit length.
   */
  public static S1Angle getDistance(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    Preconditions.checkArgument(isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(isUnitLength(b), "S2Point not normalized: %s", b);
    return S1Angle.radians(getDistanceRadians(x, a, b, aCrossB));
  }

  /** Gets the distance from {@code p} to {@code e}. */
  public static S1ChordAngle getDistance(S2Point p, S2Edge e) {
    return updateMinDistance(p, e, S1ChordAngle.INFINITY);
  }

  /**
   * Gets the minimum of the distance from {@code p} to {@code e} and {@code minDistance}.
   * TODO(user): Deprecate this and other methods when distance measurer is available.
   */
  public static S1ChordAngle updateMinDistance(S2Point p, S2Edge e, S1ChordAngle minDistance) {
    return updateMinDistance(p, e.getStart(), e.getEnd(), minDistance);
  }

  /**
   * Return the minimum of the distance from {@code x} to any point on edge ab and the given {@code
   * minDistance}. The case {@code a.equals(b)} is handled correctly. If the given minDistance is
   * the smaller, the returned S1ChordAngle is the same instance as the given minDistance. Thus
   * callers may safely compare the returned value to minDistance by reference to see if the minimum
   * distance was updated.
   *
   * <p>Note this is based on the AlwaysUpdateMinDistance method in the C++ implementation.
   * <p>TODO(user): Deprecate this and other methods when distance measurer is available.
   *
   * @throws IllegalArgumentException if the parameters are not all unit length.
   */
  @SuppressWarnings("ReferenceEquality")
  public static S1ChordAngle updateMinDistance(
      S2Point x, S2Point a, S2Point b, S1ChordAngle minDistance) {
    Preconditions.checkArgument(isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(isUnitLength(b), "S2Point not normalized: %s", b);

    double xa2 = x.sub(a).norm2();
    double xb2 = x.sub(b).norm2();

    S1ChordAngle dist = maybeUpdateMinInteriorDistance(x, a, b, xa2, xb2, minDistance);
    if (dist != minDistance) { // NOTE: deliberate use of reference equality
      // The minimum distance is less than minDistance, and is attained along the edge interior.
      return dist;
    }

    // Otherwise the minimum distance is to one of the endpoints.
    double dist2 = min(xa2, xb2);
    if (dist2 >= minDistance.getLength2()) {
      return minDistance;
    }
    return S1ChordAngle.fromLength2(dist2);
  }

  /**
   * Returns true if the distance from {@code x} to any point on edge {@code ab} is less than the
   * given {@code minDistance}.
   *
   * @throws IllegalArgumentException if the parameters are not all unit length.
   */
  @SuppressWarnings("ReferenceEquality")
  public static boolean isDistanceLess(S2Point x, S2Point a, S2Point b, S1ChordAngle minDistance) {
    // updateMinDistance is guaranteed to return the same instance if the distance is not less.
    return minDistance != updateMinDistance(x, a, b, minDistance);
  }

  /**
   * Returns true if the minimum distance from {@code x} to the edge {@code ab} is attained at an
   * interior point of {@code ab} (i.e., not an endpoint), and that distance is less than
   * "minDistance". (Specify minDistance.successor() for "less than or equal to".)
   *
   * @throws IllegalArgumentException if the parameters are not all unit length.
   */
  @SuppressWarnings("ReferenceEquality")
  public static boolean isInteriorDistanceLess(
      S2Point x, S2Point a, S2Point b, S1ChordAngle minDistance) {
    // updateMinInteriorDistance is guaranteed to return the same instance if the distance is not
    // less.
    return minDistance != updateMinInteriorDistance(x, a, b, minDistance);
  }

  /**
   * If the minimum distance from X to AB is attained at an interior point of AB (i.e., not an
   * endpoint), and that distance is less than the given "minDistance", then this method returns a new
   * S1ChordAngle containing that smaller distance. Otherwise, returns "minDistance" (the same
   * instance).
   *
   * @throws IllegalArgumentException if the parameters are not all unit length.
   */
  public static S1ChordAngle updateMinInteriorDistance(
      S2Point x, S2Point a, S2Point b, S1ChordAngle minDistance) {
    double xa2 = x.sub(a).norm2();
    double xb2 = x.sub(b).norm2();
    return maybeUpdateMinInteriorDistance(x, a, b, xa2, xb2, minDistance);
  }

  /**
   * If the minimum distance from X to AB is attained at an interior point of AB (i.e., not an
   * endpoint), and that distance is less than "minDistance", then returns a new S1ChordAngle with
   * that minimum distance. Otherwise, returns "minDistance", as the same instance.
   *
   * <p>Note that this is based on the AlwaysUpdateMinInteriorDistance method in the C++
   * implementation, with always_update==false.
   *
   * @throws IllegalArgumentException if the parameters are not all unit length.
   */
  private static S1ChordAngle maybeUpdateMinInteriorDistance(
      S2Point x, S2Point a, S2Point b, double xa2, double xb2, S1ChordAngle minDistance) {
    Preconditions.checkArgument(isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(isUnitLength(b), "S2Point not normalized: %s", b);
    assert xa2 == x.sub(a).norm2();
    assert xb2 == x.sub(b).norm2();

    // The closest point on AB could either be one of the two vertices (the "vertex case") or in the
    // interior (the "interior case"). Let C = A x B. If X is in the spherical wedge extending from
    // A to B around the axis through C, then we are in the interior case. Otherwise we are in the
    // vertex case.
    //
    // Check whether we might be in the interior case. For this to be true, XAB and XBA must both be
    // acute angles. Checking this condition exactly is expensive, so instead we consider the planar
    // triangle ABX (which passes through the sphere's interior). The planar angles XAB and XBA are
    // always less than the corresponding spherical angles, so if we are in the interior case then
    // both of these angles must be acute.
    //
    // We check this by computing the squared edge lengths of the planar triangle ABX, and testing
    // whether angles XAB and XBA are both acute using the law of cosines:
    //
    //            | XA^2 - XB^2 | < AB^2      (*)
    //
    // This test must be done conservatively (taking numerical errors into account) since otherwise
    // we might miss a situation where the true minimum distance is achieved by a point on the edge
    // interior.
    //
    // There are two sources of error in the expression above (*). The first is that points are not
    // normalized exactly; they are only guaranteed to be within 2 * DBL_EPSILON of unit length.
    // Under the assumption that the two sides of (*) are nearly equal, the total error due to
    // normalization errors can be shown to be at most
    //
    //        2 * DBL_EPSILON * (XA^2 + XB^2 + AB^2) + 8 * DBL_EPSILON ^ 2 .
    //
    // The other source of error is rounding of results in the calculation of (*). Each of XA^2,
    // XB^2, AB^2 has a maximum relative error of 2.5 * DBL_EPSILON, plus an additional relative
    // error of 0.5 * DBL_EPSILON in the final subtraction which we further bound as 0.25 *
    // DBL_EPSILON * (XA^2 + XB^2 + AB^2) for convenience. This yields a final error bound of
    //
    //        4.75 * DBL_EPSILON * (XA^2 + XB^2 + AB^2) + 8 * DBL_EPSILON ^ 2 .
    double ab2 = a.sub(b).norm2();
    double maxError = (4.75 * DBL_EPSILON * (xa2 + xb2 + ab2) + 8 * DBL_EPSILON * DBL_EPSILON);
    if (abs(xa2 - xb2) >= (ab2 + maxError)) {
      return minDistance;
    }

    // The minimum distance might be to a point on the edge interior. Let R be closest point to X
    // that lies on the great circle through AB. Rather than computing the geodesic distance along
    // the surface of the sphere, instead we compute the "chord length" through the sphere's
    // interior. If the squared chord length exceeds min_dist.length2() then we can return "null"
    // immediately.
    //
    // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q is the point X
    // projected onto the plane through the great circle AB. The distance XQ^2 can be written as
    // (X.C)^2 / |C|^2 where C = A x B. We ignore the QR^2 term and instead use XQ^2 as a lower
    // bound, since it is faster and the corresponding distance on the Earth's surface is accurate
    // to within 1% for distances up to about 1800km.
    S2Point c = robustCrossProd(a, b);
    double c2 = c.norm2();
    double xDotC = x.dotProd(c);
    double xDotC2 = xDotC * xDotC;
    if (xDotC2 > c2 * minDistance.getLength2()) {
      // The closest point on the great circle AB is too far away. We need to test this using ">"
      // rather than ">=" because the actual minimum bound on the distance is (xDotC2 / c2), which
      // can be rounded differently than the (more efficient) multiplicative test above.
      return minDistance;
    }

    // Otherwise we do the exact, more expensive test for the interior case. This test is very
    // likely to succeed because of the conservative planar test we did initially.
    //
    // TODO(ericv): Ensure that the errors in test are accurately reflected in
    // getUpdateMinInteriorDistanceMaxError().
    S2Point cx = c.crossProd(x);
    if (a.sub(x).dotProd(cx) >= 0 || b.sub(x).dotProd(cx) <= 0) {
      return minDistance;
    }

    // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above). This calculation has good
    // accuracy for all chord lengths since it is based on both the dot product and cross product
    // (rather than deriving one from the other). However, note that the chord length representation
    // itself loses accuracy as the angle approaches Pi.
    double qr = 1 - sqrt(cx.norm2() / c2);
    double dist2 = (xDotC2 / c2) + (qr * qr);
    if (dist2 >= minDistance.getLength2()) {
      return minDistance;
    }
    return S1ChordAngle.fromLength2(dist2);
  }

  /**
   * Returns the maximum of the distance from {@code x} to any point on edge AB and the given {@code
   * maxDistance}. The case {@code a.equals(b)} is handled correctly.
   */
  public static S1ChordAngle updateMaxDistance(
      S2Point x, S2Point a, S2Point b, S1ChordAngle maxDistance) {
    S1ChordAngle dist = S1ChordAngle.max(new S1ChordAngle(x, a), new S1ChordAngle(x, b));

    // If the x to AB distance is more than a quarter of the way around the sphere, recompute for
    // greater accuracy by computing the minimum distance from x.neg() to AB, and subtracting that
    // from STRAIGHT.
    if (dist.compareTo(S1ChordAngle.RIGHT) > 0) {
      // Get the distance from the antipodal
      dist = updateMinDistance(x.neg(), a, b, S1ChordAngle.INFINITY);
      dist = S1ChordAngle.sub(S1ChordAngle.STRAIGHT, dist);
    }
    if (maxDistance.compareTo(dist) >= 0) {
      return maxDistance;
    }
    return dist;
  }

  /**
   * Like {@link #updateMinDistance(S2Point, S2Point, S2Point, S1ChordAngle)}, but returns the
   * minimum of the given minDist and the distance between the given pair of edges. (If the two
   * edges cross, the distance is zero.) The cases {@code a0.equals(a1)} and {@code b0.equals(b1)}
   * are handled correctly.
   *
   * TODO(torrey): This method should be renamed "updateMinEdgePairDistance".
   */
  public static S1ChordAngle getEdgePairMinDistance(
      final S2Point a0,
      final S2Point a1,
      final S2Point b0,
      final S2Point b1,
      S1ChordAngle minDist) {
    if (minDist.equals(S1ChordAngle.ZERO)) {
      return minDist;
    }

    // If they cross, distance is 0 and no end point is closest.
    if (robustCrossing(a0, a1, b0, b1) >= 0) {
      return S1ChordAngle.ZERO;
    }

    // Otherwise, the minimum distance is achieved at an endpoint of at least one of the two edges.
    // The calculation below computes each of the six distances twice (this could be optimized).
    minDist = updateMinDistance(a0, b0, b1, minDist);
    minDist = updateMinDistance(a1, b0, b1, minDist);
    minDist = updateMinDistance(b0, a0, a1, minDist);
    minDist = updateMinDistance(b1, a0, a1, minDist);
    return minDist;
  }

  /** Gets distance between edges with no minimum distance. */
  public static S1ChordAngle getEdgePairDistance(
      final S2Point a0, final S2Point a1, final S2Point b0, final S2Point b1) {
    return getEdgePairMinDistance(a0, a1, b0, b1, S1ChordAngle.INFINITY);
  }

  // TODO(user): Create an EdgePairDistanceLess class for repeatedly testing a single left
  // edge against many right edges.
  /**
   * Returns true if the minimum distance between the two edges is less than the given distance. May
   * be significantly faster than computing the minimum distance and comparing it to the given
   * distance.
   */
  public static boolean isEdgePairDistanceLess(
      final S2Point a0,
      final S2Point a1,
      final S2Point b0,
      final S2Point b1,
      S1ChordAngle distance) {
     if (distance.isZero()) {
      return false;
    }

    // If they cross, distance is 0.
    if (robustCrossing(a0, a1, b0, b1) >= 0) {
      return distance.compareTo(S1ChordAngle.ZERO) != 0;
    }

    // Otherwise the minimum distance is achieved at an endpoint of at least one of the endpoints of
    // the two edges.
    double r2 = distance.getLength2();
    return compareEdgeDistance(a0, b0, b1, r2) < 0
        || compareEdgeDistance(a1, b0, b1, r2) < 0
        || compareEdgeDistance(b0, a0, a1, r2) < 0
        || compareEdgeDistance(b1, a0, a1, r2) < 0;
  }

  /**
   * Updates the {@code result} with points that achieve the minimum distance between edges a0a1 and
   * b0b1, where {@code a} is a point on a0a1 and {@code b} is a point on b0b1. If the two edges
   * intersect, {@code a} and {@code b} are both equal to the intersection point. Handles {@code
   * a0.equals(a1)} and {@code b0.equals(b1)} correctly.
   */
  @SuppressWarnings("ReferenceEquality")
  public static void getEdgePairClosestPoints(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, S2Shape.MutableEdge result) {
    // If they cross, distance is 0 and no end point is closest.
    if (robustCrossing(a0, a1, b0, b1) > 0) {
      S2Point intersection = getIntersection(a0, a1, b0, b1);
      result.a = intersection;
      result.b = intersection;
      return;
    }

    S1ChordAngle actualMin = S1ChordAngle.INFINITY;
    ClosestPoint closest = ClosestPoint.NONE;
    S1ChordAngle newMin = updateMinDistance(a0, b0, b1, actualMin);
    // Compare by reference is safe here, as newMin refers to the same S1ChordAngle instance as
    // actualMin if the minimum distance was not updated.
    if (newMin != actualMin) {
      closest = ClosestPoint.A0;
      actualMin = newMin;
    }
    newMin = updateMinDistance(a1, b0, b1, actualMin);
    // As above, compare by reference is safe here.
    if (newMin != actualMin) {
      closest = ClosestPoint.A1;
      actualMin = newMin;
    }
    newMin = updateMinDistance(b0, a0, a1, actualMin);
    // As above, compare by reference is safe here.
    if (newMin != actualMin) {
      closest = ClosestPoint.B0;
      actualMin = newMin;
    }
    newMin = updateMinDistance(b1, a0, a1, actualMin);
    // As above, compare by reference is safe here.
    if (newMin != actualMin) {
      closest = ClosestPoint.B1;
    }

    switch (closest) {
      case A0:
        result.a = a0;
        result.b = project(a0, b0, b1);
        return;
      case A1:
        result.a = a1;
        result.b = project(a1, b0, b1);
        return;
      case B0:
        result.a = project(b0, a0, a1);
        result.b = b0;
        return;
      case B1:
        result.a = project(b1, a0, a1);
        result.b = b1;
        return;
      default:
        throw new IllegalArgumentException(
            Platform.formatString(
                "Unknown ClosestPoint case when finding closest points of %s:%s and %s:%s",
                a0, a1, b0, b1));
    }
  }

  /**
   * Returns true if every point on edge B=b0b1 is no further than "tolerance" from some point on
   * edge A=a0a1. Equivalently, returns true if the directed Hausdorff distance from B to A is no
   * more than "tolerance". Requires that tolerance is less than 90 degrees.
   */
  public static boolean isEdgeBNearEdgeA(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, S1Angle tolerance) {
    // TODO(torrey): Optimize this function to use S1ChordAngle rather than S1Angle. Same in C++.
    assert tolerance.radians() < PI / 2;
    assert tolerance.radians() > 0;

    // The point on edge B=b0b1 furthest from edge A=a0a1 is either b0, b1, or some interior point
    // on B. If it is an interior point on B, then it must be one of the two points where the great
    // circle containing B (circ(B)) is furthest from the great circle containing A (circ(A)). At
    // these points, the distance between circ(B) and circ(A) is the angle between the planes
    // containing them.
    S2Point aOrtho = robustCrossProd(a0, a1).normalize();
    S2Point aNearestB0 = project(b0, a0, a1, aOrtho);
    S2Point aNearestB1 = project(b1, a0, a1, aOrtho);

    // If aNearestB0 and aNearestB1 have opposite orientation from a0 and a1, we invert aOrtho so
    // that it points in the same direction as aNearestB0 x aNearestB1. This helps us handle the
    // case where A and B are oppositely oriented but otherwise might be near each other. We check
    // orientation and invert rather than computing aNearestB0 x aNearestB1 because those two points
    // might be equal, and have an unhelpful cross product.
    if (sign(aOrtho, aNearestB0, aNearestB1) < 0) {
      aOrtho = aOrtho.neg();
    }

    // To check if all points on B are within tolerance of A, we first check to see if the endpoints
    // of B are near A. If they are not, B is not near A.
    S1Angle b0Distance = new S1Angle(b0, aNearestB0);
    S1Angle b1Distance = new S1Angle(b1, aNearestB1);
    if (b0Distance.greaterThan(tolerance) || b1Distance.greaterThan(tolerance)) {
      return false;
    }

    // If b0 and b1 are both within tolerance of A, we check to see if the angle between the planes
    // containing B and A is greater than tolerance. If it is not, no point on B can be further than
    // tolerance from A (recall that we already know that b0 and b1 are close to A, and S2Edges are
    // all shorter than 180 degrees). The angle between the planes containing circ(A) and circ(B) is
    // the angle between their normal vectors.
    S2Point bOrtho = robustCrossProd(b0, b1).normalize();
    S1Angle planarAngle = new S1Angle(aOrtho, bOrtho);
    if (planarAngle.lessOrEquals(tolerance)) {
      return true;
    }

    // When planarAngle >= PI/2, there are only two possible scenarios:
    //
    //  1.) b0 and b1 are closest to A at distinct endpoints of A, in which case the opposite
    //      orientation of aOrtho and bOrtho means that A and B are in opposite hemispheres and
    //      hence not close to each other.
    //
    //  2.) b0 and b1 are closest to A at the same endpoint of A, in which case the orientation of
    //      aOrtho was chosen arbitrarily to be that of a0 cross a1. B must be shorter than
    //      2*tolerance and all points in B are close to one endpoint of A, and hence to A.
    //
    // Note that this logic *must* be used when planar_angle >= Pi/2 because the code beyond does
    // not handle the case where the maximum distance is attained at the interior point of B that is
    // equidistant from the endpoints of A. This happens when B intersects the perpendicular
    // bisector of the endpoints of A in the hemisphere opposite A's midpoint.
    if (planarAngle.radians() >= M_PI_2) {
      boolean b0NearerA0 = new S1Angle(b0, a0).lessThan(new S1Angle(b0, a1));
      boolean b1NearerA0 = new S1Angle(b1, a0).lessThan(new S1Angle(b1, a1));
      // Are the endpoints of B both closest to the same endpoint of A?
      return b0NearerA0 == b1NearerA0;
    }

    // Otherwise, if either of the two points on circ(B) where circ(B) is furthest from circ(A) lie
    // on edge B, edge B is not near edge A.
    //
    // The normalized projection of aOrtho onto the plane of circ(B) is one of the two points along
    // circ(B) where it is furthest from circ(A). The other is -1 times the normalized projection.

    // Note that the formula (A - (A.B) * B) loses accuracy when |A.B| ~= 1, so instead we compute
    // it using two cross products. (The first product does not need RobustCrossProd since its
    // arguments are perpendicular.)
    S2Point furthest = bOrtho.crossProd(robustCrossProd(aOrtho, bOrtho)).normalize();
    S2Point furthestInv = furthest.neg();

    // A point p lies on B if you can proceed from bOrtho to b0 to p to b1 and back to bOrtho
    // without ever turning right. We test this for furthest and furthestInv, and return true if
    // neither point lies on B.
    return !((sign(bOrtho, b0, furthest) > 0 && sign(furthest, b1, bOrtho) > 0)
        || (sign(bOrtho, b0, furthestInv) > 0 && sign(furthestInv, b1, bOrtho) > 0));
  }

  /**
   * Like {@link #updateMaxDistance(S2Point, S2Point, S2Point, S1ChordAngle)}, but computes the
   * maximum of the given 'maxDist' and the maximum distance between the given pair of edges. The
   * maximum distance between edges is usually between the furthest separated pair of endpoints, but
   * due to the curve of the sphere it is possible for the maximum distance to be between one
   * edge's endpoint and the interior of the other edge, or even between the interiors of both edges
   * if one edge is antipodal to the other at some point. The cases {@code a0.equals(a1)} and
   * {@code b0.equals(b1)} are handled correctly.
   */
  public static S1ChordAngle getEdgePairMaxDistance(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, S1ChordAngle maxDist) {
    // If maxDist is already the maximum it can't be increased.
    if (maxDist.equals(S1ChordAngle.STRAIGHT)) {
      return maxDist;
    }

    // If one edge intersects the reflection of the other, they are antipodal at one point.
    if (S2EdgeUtil.robustCrossing(a0, a1, b0.neg(), b1.neg()) >= 0) {
      return S1ChordAngle.STRAIGHT;
    }

    // Otherwise, the maximum distance is achieved at an endpoint of at least one of the two edges.
    // The calculation below computes all six distances twice (this could be optimized).
    maxDist = updateMaxDistance(a0, b0, b1, maxDist);
    maxDist = updateMaxDistance(a1, b0, b1, maxDist);
    maxDist = updateMaxDistance(b0, a0, a1, maxDist);
    maxDist = updateMaxDistance(b1, a0, a1, maxDist);
    return maxDist;
  }

  /**
   * A more efficient version of getDistance() where the cross product of the endpoints has been
   * precomputed and the result is returned as a direct radian measure rather than wrapping it in an
   * S1Angle. This is the recommended method for making large numbers of back-to-back edge distance
   * tests, since it allocates no objects. The inputs are assumed to be unit length; results are
   * undefined if they are not.
   */
  public static double getDistanceRadians(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    // There are three cases. If X is located in the spherical wedge defined by A, B, and the axis
    // A x B, then the closest point is on the segment AB. Otherwise the closest point is either A
    // or B; the dividing line between these two cases is the great circle passing through (A x B)
    // and the midpoint of AB.
    if (scalarTripleProduct(a, x, aCrossB) > 0 && scalarTripleProduct(b, aCrossB, x) > 0) {
      // The closest point to X lies on the segment AB. We compute the distance to the corresponding
      // great circle. The result is accurate for small distances but not necessarily for large
      // distances (approaching Pi/2).
      double sinDist = abs(x.dotProd(aCrossB)) / aCrossB.norm();
      return asin(min(1.0, sinDist));
    }

    // Otherwise, the closest point is either A or B. The cheapest method is just to compute the
    // minimum of the two linear (as opposed to spherical) distances and convert the result to an
    // angle. Again, this method is accurate for small but not large distances (approaching Pi).
    double linearDist2 = min(distance2(x, a), distance2(x, b));
    return 2 * asin(min(1.0, 0.5 * sqrt(linearDist2)));
  }

  /** Returns the squared distance from {@code a} to {@code b}. */
  private static final double distance2(S2Point a, S2Point b) {
    double dx = a.getX() - b.getX();
    double dy = a.getY() - b.getY();
    double dz = a.getZ() - b.getZ();
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * As {@link #getClosestPoint(S2Point, S2Point, S2Point)}, returns the point on edge AB closest to
   * X, but faster if the cross product between 'a' and 'b' has already been computed. All points
   * must be unit length; results are undefined if that is not the case.
   *
   * @deprecated Use 'project' instead.
   */
  @Deprecated
  public static S2Point getClosestPoint(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    assert isUnitLength(a);
    assert isUnitLength(b);
    assert isUnitLength(x);

    // Find the closest point to X along the great circle through AB.
    S2Point p = x.sub(aCrossB.mul(x.dotProd(aCrossB) / aCrossB.norm2()));

    // If this point is on the edge AB, then it's the closest point.
    if (scalarTripleProduct(a, p, aCrossB) > 0 && scalarTripleProduct(b, aCrossB, p) > 0) {
      return p.normalize();
    }

    // Otherwise, the closest point is either A or B.
    return x.getDistance2(a) <= x.getDistance2(b) ? a : b;
  }

  /**
   * Returns the point on edge AB closest to X. All points must be unit length; results are
   * undefined if that is not the case.
   *
   * @deprecated Use 'project' instead.
   */
  @Deprecated
  public static S2Point getClosestPoint(S2Point x, S2Point a, S2Point b) {
    return getClosestPoint(x, a, b, robustCrossProd(a, b));
  }

  /**
   * Returns the point along the edge AB that is closest to the point X. Requires that A, B, and X
   * have unit length.
   *
   * <p>The fractional distance of this point along the edge AB can be obtained using {@link
   * #getDistanceFraction(S2Point, S2Point, S2Point)}.
   */
  public static S2Point project(S2Point x, S2Point a, S2Point b) {
    return project(x, a, b, robustCrossProd(a, b));
  }

  /**
   * Returns the point along the edge AB that is closest to the point X. The cross product aCrossB
   * does not need to be normalized, but should be computed using {@link
   * S2RobustCrossProd#robustCrossProd(S2Point, S2Point)} for the most accurate results. Requires
   * that A, B, and X have unit length. This version is slightly more efficient if the cross product
   * of the two endpoints has been precomputed.
   *
   * <p>The result is within {@link #PROJECT_PERPENDICULAR_ERROR} of the edge AB. Note that
   * this only bounds the error perpendicular to the edge, not the error parallel to it, and in
   * particular, if X is nearly perpendicular to the plane containing AB, the accuracy of the
   * returned result may be low.
   *
   * <p>The fractional distance of the returned point along the edge AB can be obtained using {@link
   * #getDistanceFraction(S2Point, S2Point, S2Point)}.
   */
  public static S2Point project(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    assert isUnitLength(a);
    assert isUnitLength(b);
    assert isUnitLength(x);

    // TODO(user): When X is nearly perpendicular to the plane containing AB, the result is
    // guaranteed to be close to the edge AB but may be far from the true projected result. This
    // could be fixed by computing the product (A x B) x X x (A x B) using methods similar to
    // S2RobustCrossProd.robustCrossProd() and S2.getIntersection(). However note that the error
    // tolerance would need to be significantly larger in order for this calculation to succeed in
    // double precision most of the time. For example to avoid higher precision when X is within 60
    // degrees of AB the minimum error would be 18 * DBL_ERR, and to avoid higher precision when X
    // is within 87 degrees of AB the minimum error would be 120 * DBL_ERR.

    // The following is not necessary to meet accuracy guarantees but helps to avoid unexpected
    // results in unit tests.
    if (x.equalsPoint(a) || x.equalsPoint(b)) {
      return x;
    }

    // Find the closest point to X along the great circle through AB.  Note that we use "n" rather
    // than aCrossB in the final cross product in order to avoid the possibility of underflow.
    // TODO(user): Require caller to normalize aCrossB instead of repeatedly doing it here.
    S2Point n = aCrossB.normalize();
    S2Point p = robustCrossProd(n, x).crossProd(n).normalize();

    // If this point is on the edge AB, then it's the closest point.
    S2Point pn = p.crossProd(n);
    if (S2Predicates.sign(p, n, a, pn) > 0 && S2Predicates.sign(p, n, b, pn) < 0) {
      return p;
    }

    // Otherwise, the closest point is either A or B.
    return (x.getDistance2(a) <= x.getDistance2(b)) ? a : b;
  }

  /**
   * Returns the normalized point at distance "r" along the ray with the given origin and direction.
   * "dir" is required to be perpendicular to "origin" (since this is how directions on the sphere
   * are represented).
   *
   * <p>This function is similar to getPointOnLine() except that (1) the first two arguments are
   * required to be perpendicular and (2) it is much faster. It can be used as an alternative to
   * repeatedly calling getPointOnLine() by computing "dir" as:
   *
   * {@snippet :
   * S2Point dir = S2RobustCrossProd.robustCrossProd(a, b).crossProd(a).normalize();
   * }
   *
   * REQUIRES: "origin" and "dir" are perpendicular to within the tolerance of the calculation
   * above.
   */
  public static S2Point getPointOnRay(S2Point origin, S2Point dir, S1ChordAngle r) {
    assert isUnitLength(origin);
    assert isUnitLength(dir);
    // The error bound below includes the error in computing the dot product.
    assert origin.dotProd(dir) <= ROBUST_CROSS_PROD_ERROR.radians() + 0.75 * DBL_EPSILON;

    // Mathematically the result should already be unit length, but we normalize it anyway to ensure
    // that the error is within acceptable bounds. (Otherwise errors can build up when the result of
    // one interpolation is fed into another interpolation.)
    //
    // Note that it is much cheaper to compute the sine and cosine of an S1ChordAngle than an
    // S1Angle.
    // (cos(r) * origin + sin(r) * dir).normalize();
    return origin.mul(S1ChordAngle.cos(r)).add(dir.mul(S1ChordAngle.sin(r))).normalize();
  }

  /**
   * Returns the normalized point at distance "r" along the ray with the given origin and direction.
   * "dir" is required to be perpendicular to "origin" (since this is how directions on the sphere
   * are represented).
   *
   * <p>This function is similar to getPointOnLine() except that (1) the first two arguments are
   * required to be perpendicular and (2) it is much faster. It can be used as an alternative to
   * repeatedly calling getPointOnLine() by computing "dir" as:
   *
   * <pre>
   *   S2Point dir = S2RobustCrossProd.robustCrossProd(a, b).crossProd(a).normalize();
   * </pre>
   *
   * REQUIRES: "origin" and "dir" are perpendicular to within the tolerance of the calculation
   * above.
   */
  public static S2Point getPointOnRay(S2Point origin, S2Point dir, S1Angle r) {
    // See comments above.
    assert isUnitLength(origin);
    assert isUnitLength(dir);
    assert origin.dotProd(dir) <= ROBUST_CROSS_PROD_ERROR.radians() + 0.75 * DBL_EPSILON;

    // (cos(r) * origin + sin(r) * dir).normalize();
    return origin.mul(r.cos()).add(dir.mul(r.sin())).normalize();
  }

  /**
   * Returns the normalized point at distance "r" from A along the line AB. Not quite as fast as the
   * variant taking S1ChordAngle to specify distance (i.e. 140 ns vs. 116 ns).
   *
   * <p>Note that the line AB has a well-defined direction even when A and B are antipodal or nearly
   * so. If A == B then an arbitrary direction is chosen.
   */
  public static S2Point getPointOnLine(S2Point a, S2Point b, S1Angle r) {
    // Use robustCrossProd() to compute the tangent vector at A towards B. This technique is robust
    // even when A and B are antipodal or nearly so.
    S2Point dir = robustCrossProd(a, b).crossProd(a).normalize();
    return getPointOnRay(a, dir, r);
  }

  /**
   * Returns the normalized point at distance "r" from A along the line AB. Slightly faster than the
   * variant taking an S1Angle to specify distance, but cannot accurately represent distances near
   * 180 degrees due to the limitations of S1ChordAngle.
   *
   * <p>Note that the line AB has a well-defined direction even when A and B are antipodal or nearly
   * so. If A == B then an arbitrary direction is chosen.
   */
  public static S2Point getPointOnLine(S2Point a, S2Point b, S1ChordAngle r) {
    // Use robustCrossProd() to compute the tangent vector at A towards B. This technique is robust
    // even when A and B are antipodal or nearly so.
    S2Point dir = robustCrossProd(a, b).crossProd(a).normalize();
    return getPointOnRay(a, dir, r);
  }

  /**
   * Returns the normalized S2Point to the left of the edge from "a" to "b" which is the distance
   * "r" away from "a", orthogonal to the specified edge:
   *
   * <pre>
   *     c (result)
   *     ^
   *     |  r
   *     |
   *     a --------> b
   * </pre>
   */
  public static S2Point getPointToLeft(S2Point a, S2Point b, S1Angle r) {
    return getPointOnRay(a, robustCrossProd(a, b).normalize(), r);
  }

  /**
   * Returns the normalized S2Point to the left of the edge from "a" to "b" which is the distance
   * "r" away from "a", orthogonal to the specified edge. This version is faster than the version
   * using S1Angle.
   *
   * <pre>
   *     c (result)
   *     ^
   *     |  r
   *     |
   *     a --------> b
   * </pre>
   */
  public static S2Point getPointToLeft(S2Point a, S2Point b, S1ChordAngle r) {
    return getPointOnRay(a, robustCrossProd(a, b).normalize(), r);
  }

  /**
   * Returns the normalized S2Point to the right of the edge from "a" to "b" which is the distance
   * "r" away from "a", orthogonal to the specified edge:
   *
   * <pre>
   *     a --------> b
   *     |
   *     |  r
   *     v
   *     c (result)
   * </pre>
   */
  public static S2Point getPointToRight(S2Point a, S2Point b, S1Angle r) {
    return getPointOnRay(a, robustCrossProd(b, a).normalize(), r);
  }

  /**
   * Returns the normalized S2Point to the right of the edge from "a" to "b" which is the distance
   * "r" away from "a", orthogonal to the specified edge. This version is faster than the version
   * using S1Angle.
   *
   * <pre>
   *     a --------> b
   *     |
   *     |  r
   *     v
   *     c (result)
   * </pre>
   */
  public static S2Point getPointToRight(S2Point a, S2Point b, S1ChordAngle r) {
    return getPointOnRay(a, robustCrossProd(b, a).normalize(), r);
  }

  /**
   * Returns a normalized point X along the geodesic through the points A and B. The distance from A
   * to X is specified by the parameter "ax". If 'ax' is not in the range 0 to the angle between a
   * and b, this will extrapolate along the great circle through 'a' and 'b', past 'b' for larger
   * angles, or past 'a' for negative angles. Requires that {@code a} and {@code b} are unit length.
   *
   * <p>This is a slightly more efficient version of {@link #interpolateAtDistance(S1Angle, S2Point,
   * S2Point)} that can be used when the distance AB is already known. See also {@link
   * #interpolate(double, S2Point, S2Point)} if you want to interpolate along AB to a point
   * expressed as a fraction of the distance AB.
   *
   * @deprecated This method is not as accurate as getPointOnLine(a, b, ax), and only slightly
   *     faster. If computing 'ab' is included, this method is slower. Just call getPointOnLine.
   */
  @Deprecated
  public static S2Point interpolateAtDistance(S1Angle ax, S2Point a, S2Point b, S1Angle ab) {
    assert isUnitLength(a);
    assert isUnitLength(b);

    double axRadians = ax.radians();
    double abRadians = ab.radians();

    // If the distance "ax" is zero, then produce "a" directly.
    // This avoids producing NaN in the case where both ax = 0 and ab = 0.
    if (axRadians == 0.0) {
      return a;
    }

    // The result X is some linear combination X = e*A + f*B of the input points. The fractions "e"
    // and "f" can be derived by looking at the components of this equation that are parallel and
    // perpendicular to A. Let E = e*A and F = f*B. Then OEXF is a parallelogram. You can obtain the
    // distance f = OF by considering the similar triangles produced by
    // dropping perpendiculars from the segments OF and OB to OA.
    double f = sin(axRadians) / sin(abRadians);

    // Form the dot product of the first equation with A to obtain A.X = e*A.A + f*A.B. Since A, B,
    // and X are all unit vectors, cos(ax) = e*1 + f*cos(ab), so
    double e = cos(axRadians) - f * cos(abRadians);

    // Mathematically speaking, if "a" and "b" are unit length then the result is unit length as
    // well. But we normalize it anyway to prevent points from drifting away from unit length when
    // multiple interpolations are done in succession (i.e. the result of one interpolation is fed
    // into another).
    return a.mul(e).add(b.mul(f)).normalize();
  }

  /**
   * Like {@link #interpolate(double, S2Point, S2Point)}, returns the normalized point X along the
   * line segment AB, except the distance from A is specified by the parameter "ax" which represents
   * the desired distance from A to the result X, rather than a fraction between 0 and 1. Requires
   * that {@code a} and {@code b} are unit length. Angles 'ax' that are not in the range 0 to the
   * angle between a and b will result in extrapolation along the great circle through 'a' and 'b',
   * either past 'b' for larger angles, or past 'a' for negative angles.
   *
   * @deprecated Call getPointOnLine(a, b, ax).
   */
  @Deprecated
  @InlineMe(
      replacement = "S2EdgeUtil.getPointOnLine(a, b, ax)",
      imports = "com.google.common.geometry.S2EdgeUtil")
  public static S2Point interpolateAtDistance(S1Angle ax, S2Point a, S2Point b) {
    return getPointOnLine(a, b, ax);
  }

  /**
   * Return the normalized point X along the line segment AB whose distance from A is the given
   * fraction "t" of the distance AB. Does NOT require that "t" be between 0 and 1. Note that all
   * distances are measured on the surface of the sphere, so this is more complicated than just
   * computing (1-t)*a + t*b and normalizing the result.
   */
  public static S2Point interpolate(S2Point a, S2Point b, double t) {
    if (t == 0) {
      return a;
    }
    if (t == 1) {
      return b;
    }
    S1Angle ab = new S1Angle(a, b);
    return getPointOnLine(a, b, S1Angle.radians(t * ab.radians()));
  }

  /**
   * Return the normalized point X along the line segment AB whose distance from A is the given
   * fraction "t" of the distance AB. Does NOT require that "t" be between 0 and 1. Note that all
   * distances are measured on the surface of the sphere, so this is more complicated than just
   * computing (1-t)*a + t*b and normalizing the result.
   *
   * @deprecated Call interpolate(S2Point, S2Point, double).
   */
  @Deprecated
  @InlineMe(
      replacement = "S2EdgeUtil.interpolate(a, b, t)",
      imports = {"com.google.common.geometry.S2EdgeUtil"})
  public static S2Point interpolate(double t, S2Point a, S2Point b) {
    return S2EdgeUtil.interpolate(a, b, t);
  }

  /**
   * Returns the maximum error in the result of {@link #updateMinInteriorDistance(S2Point, S2Point,
   * S2Point, S1ChordAngle)}, assuming that all input points are normalized to within the bounds
   * guaranteed by {@link S2Point#normalize()}. The error can be added or subtracted from an
   * S1ChordAngle "x" using {@code x.plusError(error)}.
   */
  static double getUpdateMinInteriorDistanceMaxError(S1ChordAngle distance) {
    // If a point is more than 90 degrees from an edge, then the minimum distance is always to one
    // of the endpoints, not to the edge interior.
    if (distance.compareTo(S1ChordAngle.RIGHT) >= 0) {
      return 0.0;
    }

    // This bound includes all source of error, assuming that the input points are normalized to
    // within the bounds guaranteed to S2Point::Normalize(). "a" and "b" are components of chord
    // length that are perpendicular and parallel to plane containing the edge respectively.
    double b = min(1.0, 0.5 * distance.getLength2());
    double a = sqrt(b * (2 - b));
    return ((2.5 + 2 * sqrt(3) + 8.5 * a) * a
            + (2 + 2 * sqrt(3) / 3 + 6.5 * (1 - b)) * b
            + (23 + 16 / sqrt(3)) * DBL_EPSILON)
        * DBL_EPSILON;
  }

  /**
   * Returns the maximum error in the result of {@link #updateMinDistance(S2Point, S2Point, S2Point,
   * S1ChordAngle)}, assuming that all input points are normalized to within the bounds guaranteed
   * by {@link S2Point#normalize()}.
   *
   * <p>Note that accuracy goes down as the distance approaches 0 degrees or 180 degrees (for
   * different reasons). Near 0 degrees the error is acceptable for all practical purposes (about
   * 1.2e-15 radians ~= 8 nanometers). For exactly antipodal points the maximum error is quite high
   * (0.5 meters), but this error drops rapidly as the points move away from antipodality
   * (approximately 1 millimeter for points that are 50 meters from antipodal, and 1 micrometer for
   * points that are 50km from antipodal).
   */
  static double getUpdateMinDistanceMaxError(S1ChordAngle dist) {
    // There are two cases for the maximum error in updateMinDistance(), depending on whether the
    // closest point is interior to the edge.
    return max(getUpdateMinInteriorDistanceMaxError(dist), dist.getS2PointConstructorMaxError());
  }

  /**
   * Compute the intersection point of (a0, a1) and (b0, b1) using exact arithmetic. Note that the
   * result is not exact because it is rounded to double precision. Also, the intersection point is
   * not guaranteed to have the correct sign (i.e., the return value may need to be negated).
   */
  static S2Point getIntersectionExact(S2Point a0, S2Point a1, S2Point b0, S2Point b1) {
    // Since we are using exact arithmetic, we don't need to worry about numerical stability.
    ExactIntersection exact = ExactIntersection.create(a0, a1, b0, b1);

    // The last two operations are done in double precision, which creates a directional error of
    // up to 2 * DBL_EPSILON. (BigPoint.toS2Point() and S2Point.normalize() each contribute up to
    // DBL_EPSILON of directional error.)
    S2Point x = exact.result.toS2Point().normalize();

    if (x.equals(ZERO)) {
      // The two edges are exactly collinear, but we still consider them to be "crossing" because of
      // simulation of simplicity. Out of the four endpoints, exactly two lie in the interior of
      // the other edge. Of those two we return the one that is lexicographically smallest.
      x = new S2Point(10, 10, 10); // Greater than any valid S2Point
      S2Point aNorm = exact.aNorm.toS2Point().normalize();
      S2Point bNorm = exact.bNorm.toS2Point().normalize();
      // Note: To support antipodal edges properly, we would need to add a crossProd() function that
      // computes the cross product using simulation of simplicity and rounds the result to the
      // nearest floating-point representation.
      Preconditions.checkArgument(
          !(aNorm.equals(ZERO) || bNorm.equals(ZERO)),
          "Exactly antipodal edges not supported by getIntersectionExact");
      x = closestAcceptableEndpoint(a0, a1, aNorm, b0, b1, bNorm, x);
    }
    assert isUnitLength(x);
    return x;
  }

  /** An exact intersection result and the normals. */
  static final class ExactIntersection {
    public final BigPoint result;
    public final BigPoint aNorm;
    public final BigPoint bNorm;

    ExactIntersection(BigPoint result, BigPoint aNorm, BigPoint bNorm) {
      this.result = result;
      this.aNorm = aNorm;
      this.bNorm = bNorm;
    }

    /** Computes the intersection and residual normals from the given points. */
    public static ExactIntersection create(S2Point a0, S2Point a1, S2Point b0, S2Point b1) {
      BigPoint aNormBp = new BigPoint(a0).crossProd(new BigPoint(a1));
      BigPoint bNormBp = new BigPoint(b0).crossProd(new BigPoint(b1));
      BigPoint xBp = aNormBp.crossProd(bNormBp);
      return new ExactIntersection(xBp, aNormBp, bNormBp);
    }
  }

  /**
   * Returns the approximate intersection point of the edges (a0,a1) and (b0,b1), and writes to
   * resultError a bound on its error.
   *
   * <p>The intersection point is not guaranteed to have the correct sign, i.e., it may need to be
   * negated.
   */
  static S2Point getIntersectionApprox(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, ResultError resultError) {
    // Sort the two edges so that (a0,a1) is longer, breaking ties in a deterministic way that does
    // not depend on the ordering of the endpoints. This is desirable for two reasons:
    //  - So that the result doesn't change when edges are swapped or reversed.
    //  - It reduces error, since the first edge is used to compute the edge normal (where a longer
    //    edge means less error), and the second edge is used for interpolation (where a shorter
    //    edge means less error).
    double aLen2 = a1.getDistance2(a0);
    double bLen2 = b1.getDistance2(b0);
    if ((aLen2 < bLen2) || ((aLen2 == bLen2) && compareEdges(a0, a1, b0, b1))) {
      return getIntersectionApproxSorted(b0, b1, a0, a1, resultError);
    } else {
      return getIntersectionApproxSorted(a0, a1, b0, b1, resultError);
    }
  }

  /**
   * Returns true if (a0,a1) is less than (b0,b1) with respect to a total ordering on edges that is
   * invariant under edge reversals.
   */
  private static boolean compareEdges(S2Point a0, S2Point a1, S2Point b0, S2Point b1) {
    if (a1.lessThan(a0)) {
      S2Point temp = a0;
      a0 = a1;
      a1 = temp;
    }
    if (b1.lessThan(b0)) {
      S2Point temp = b0;
      b0 = b1;
      b1 = temp;
    }
    return a0.lessThan(b0) || (a0.equalsPoint(b0) && b0.lessThan(b1));
  }

  /**
   * Returns the approximate intersection point of the edges (a0,a1) and (b0,b1), and writes to
   * resultError a bound on its error.
   *
   * <p>Requires that the edges (a0,a1) and (b0,b1) have been sorted so that the first edge is
   * longer.
   *
   * <p>The intersection point is not guaranteed to have the correct sign, i.e., it may need to be
   * negated.
   */
  private static S2Point getIntersectionApproxSorted(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, ResultError resultError) {
    assert a1.getDistance2(a0) >= b1.getDistance2(b0);

    // Compute the normal of the plane through (a0, a1) in a stable way.
    S2Point aNormal = robustCrossProd(a0, a1);
    double aNormalLen = aNormal.norm();
    double bLen = b1.getDistance(b0);

    // Compute the projection (i.e., signed distance) of b0 and b1 onto the plane through (a0, a1).
    // Distances are scaled by the length of aNormal.
    ResultError b0ResultError = new ResultError();
    ResultError b1ResultError = new ResultError();
    double b0Dist = getProjection(b0, aNormal, aNormalLen, a0, a1, b0ResultError);
    double b1Dist = getProjection(b1, aNormal, aNormalLen, a0, a1, b1ResultError);

    // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is |b0Dist - b1Dist|.
    // Note that b0Dist and b1Dist generally have opposite signs because b0 and b1 are on opposite
    // sides of (a0, a1). The code below finds the intersection point by interpolating along the
    // edge (b0, b1) to a fractional distance of b0Dist / (b0Dist - b1Dist).
    //
    // It can be shown that the maximum error in the interpolation fraction is
    //
    //     (b0Dist * b1ResultError.error - b1Dist * b0ResultError.error) /
    //        (distSum * (distSum - errorSum))
    //
    // We save ourselves some work by scaling the result and the error bound by "distSum", since the
    // result is normalized to be unit length anyway.
    double distSum = abs(b0Dist - b1Dist);
    double errorSum = b0ResultError.error + b1ResultError.error;
    if (distSum <= errorSum) {
      // Error is unbounded in this case. Return an all-zeros S2Point with infinite error.
      resultError.error = Double.POSITIVE_INFINITY;
      return ZERO;
    }
    S2Point x = b1.mul(b0Dist).sub(b0.mul(b1Dist));

    // Finally we normalize the result and compute the corresponding error.
    double xLen2 = x.norm2();
    if (xLen2 < Double.MIN_NORMAL) {
      // If x.norm2() is less than double's minimum norm value, xLen might lose precision and the
      // result might fail to satisfy isUnitLength(). Return an all-zeros S2Point with infinite
      // error.
      resultError.error = Double.POSITIVE_INFINITY;
      return ZERO;
    }
    double xLen = sqrt(xLen2);
    double scaledInterpFactor =
        abs(b0Dist * b1ResultError.error - b1Dist * b0ResultError.error) / (distSum - errorSum);
    resultError.error =
        ((bLen * scaledInterpFactor + 2 * DBL_EPSILON * distSum) / xLen) + DBL_EPSILON;
    return x.mul(1 / xLen);
  }

  /**
   * Returns 2x the dot product of x and aNormal, and writes to resultError a bound on the error
   * given that aNormal was calculated using
   * {@link S2RobustCrossProd#robustCrossProd(S2Point, S2Point)}.
   *
   * <p>The remaining parameters allow this calculation to be computed more accurately and
   * efficiently. They include the length of aNormal (aNormalLen) and the edge endpoints a0 and a1.
   *
   * <p>Note that the 2x factor mentioned above is the result of an error reducing step. Rescaling
   * the result would result in a loss of accuracy and efficiency, and thus is not performed.
   */
  static double getProjection(
      S2Point x,
      S2Point aNormal,
      double aNormalLen,
      S2Point a0,
      S2Point a1,
      ResultError resultError) {
    // The error in the dot product is proportional to the lengths of the input vectors, so rather
    // than using x itself (a unit-length vector) we use the vectors from x to the closer of the
    // two edge endpoints. This typically reduces the error by a huge factor.
    S2Point x0 = x.sub(a0);
    S2Point x1 = x.sub(a1);
    double x0Dist2 = x0.norm2();
    double x1Dist2 = x1.norm2();

    // If both distances are the same, we need to be careful to choose one endpoint
    // deterministically so that the result does not change if the order of the endpoints is
    // reversed.
    double dist;
    double result;
    if ((x0Dist2 < x1Dist2) || (x0Dist2 == x1Dist2 && x0.lessThan(x1))) {
      dist = sqrt(x0Dist2);
      result = x0.dotProd(aNormal);
    } else {
      dist = sqrt(x1Dist2);
      result = x1.dotProd(aNormal);
    }
    // This calculation bounds the error from all sources: the computation of the normal, the
    // subtraction of one endpoint, and the dot product itself.
    //
    // For reference, the bounds that went into this calculation are:
    // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * DBL_EPSILON) * DBL_EPSILON
    // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * DBL_EPSILON
    // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * DBL_EPSILON
    resultError.error =
        DBL_EPSILON
            * (dist * ((3.5 + 2 * sqrt(3)) * aNormalLen + 32 * sqrt(3) * DBL_EPSILON)
                + 1.5 * abs(result));
    return result;
  }

  /**
   * Encapsulation of a mutable error value.
   *
   * <p>Used as an output parameter for methods that calculate double error values for their return
   * values.
   */
  static final class ResultError {
    double error;
  }

  /** Constructor is private so that this class is never instantiated. */
  private S2EdgeUtil() {}
}
