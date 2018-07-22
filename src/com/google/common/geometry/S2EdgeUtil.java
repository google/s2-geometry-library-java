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

import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import com.google.common.annotations.GwtCompatible;
import com.google.common.base.Preconditions;

/**
 * This class contains various utility functions related to edges. It collects
 * together common code that is needed to implement polygonal geometry such as
 * polylines, loops, and general polygons.
 *
 */
@GwtCompatible(serializable = false)
public strictfp class S2EdgeUtil {
  /**
   * IEEE floating-point operations have a maximum error of 0.5 ULPS (units in
   * the last place). For double-precision numbers, this works out to 2**-53
   * (about 1.11e-16) times the magnitude of the result. It is possible to
   * analyze the calculation done by getIntersection() and work out the
   * worst-case rounding error. I have done a rough version of this, and my
   * estimate is that the worst case distance from the intersection point X to
   * the great circle through (a0, a1) is about 12 ULPS, or about 1.3e-15. This
   * needs to be increased by a factor of (1/0.866) to account for the
   * edgeSpliceFraction() in S2PolygonBuilder. Note that the maximum error
   * measured by the unittest in 1,000,000 trials is less than 3e-16.
   */
  public static final S1Angle DEFAULT_INTERSECTION_TOLERANCE = S1Angle.radians(1.5e-15);

  /**
   * Threshold for small angles, that help lenientCrossing to determine whether
   * two edges are likely to intersect.
   */
  private static final double MAX_DET_ERROR = 1e-14;

  /**
   * The maximum angle between a returned vertex and the nearest point on the exact edge AB.  It is
   * equal to the maximum directional error in S2.robustCrossProd, plus the error when projecting
   * points onto a cube face.
   */
  public static final double FACE_CLIP_ERROR_RADIANS = 3 * S2.DBL_EPSILON;

  /**
   * The same angle as {@link #FACE_CLIP_ERROR_RADIANS}, expressed as a maximum distance in
   * (u,v)-space.  In other words, a returned vertex is at most this far from the exact edge AB
   * projected into (u,v)-space.
   */
  public static final double FACE_CLIP_ERROR_UV_DIST = 9 * S2.DBL_EPSILON;

  /**
   * The same angle as {@link #FACE_CLIP_ERROR_RADIANS}, expressed as the maximum error in an
   * individual u- or v-coordinate.  In other words, for each returned vertex there is a point on
   * the exact edge AB whose u- and v-coordinates differ from the vertex by at most this amount.
   */
  public static final double FACE_CLIP_ERROR_UV_COORD = 9 * S2.M_SQRT1_2 * S2.DBL_EPSILON;

  /**
   * The maximum error in IntersectRect.  If some point of AB is inside the rectangle by at least
   * this distance, the result is guaranteed to be true; if all points of AB are outside the
   * rectangle by at least this distance, the result is guaranteed to be false.  This bound assumes
   * that "rect" is a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it (e.g., by
   * 1e-10 or less).
   */
  public static final double INTERSECTS_RECT_ERROR_UV_DIST = 3 * S2.M_SQRT2 * S2.DBL_EPSILON;

  /**
   * The maximum error in a clipped point's u- or v-coordinate compared to the exact result,
   * assuming that the points A and B are in the rectangle [-1,1]x[1,1] or slightly outside it (by
   * 1e-10 or less).
   */
  public static final double EDGE_CLIP_ERROR_UV_COORD = 2.25 * S2.DBL_EPSILON;

  /**
   * The maximum error between a clipped edge or boundary point and the corresponding exact result.
   * It is equal to the error in a single coordinate because at most one coordinate is subject to
   * error.
   */
  public static final double EDGE_CLIP_ERROR_UV_DIST = 2.25 * S2.DBL_EPSILON;

  /** Max error allowed when checking if a loop boundary approximately intersects a target cell */
  public static final double MAX_CELL_EDGE_ERROR =
      S2EdgeUtil.FACE_CLIP_ERROR_UV_COORD + S2EdgeUtil.INTERSECTS_RECT_ERROR_UV_DIST;

  /**
   * Used to efficiently test a fixed edge AB against an edge chain. To use it, initialize with the
   * edge AB, call {@link #restartAt(S2Point)} with the first vertex of the chain, followed by
   * {@link #robustCrossing(S2Point)} with each subsequent vertex of the chain.
   *
   * <p>This class is <strong>not</strong> thread-safe.
   */
  public static final class EdgeCrosser {
    private final S2Point a;
    private final S2Point b;
    private final S2Point aCrossB;

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
     * True if the tangents have been computed. To reduce the number of calls to S2.expensiveCCW(),
     * we compute an outward-facing tangent at A and B if necessary.  If the plane perpendicular to
     * one of these tangents separates AB from CD (i.e., one edge on each side) then there is no
     * intersection.
     */
    private boolean hasTangents;

    /** Outward-facing tangent at A. */
    private S2Point aTangent;

    /** Outward-facing tangent at B. */
    private S2Point bTangent;

    /** Convenience constructor that calls init() with the given fixed edge AB. */
    public EdgeCrosser(S2Point a, S2Point b) {
      this.a = a;
      this.b = b;
      this.aCrossB = S2Point.crossProd(a, b);
    }

    /**
     * AB is the given fixed edge, and C is the first vertex of the vertex chain.  Equivalent to
     * using the two-arg constructor and calling restartAt(c).
     */
    public EdgeCrosser(S2Point a, S2Point b, S2Point c) {
      this(a, b);
      restartAt(c);
    }

    /** Call this method when your chain 'jumps' to a new place. */
    public void restartAt(S2Point c) {
      this.c = c;
      acb = -S2.robustCCW(a, b, c, aCrossB);
    }

    /**
     * This method is equivalent to calling the S2EdgeUtil.robustCrossing()
     * function (defined below) on the edges AB and CD. It returns +1 if there
     * is a crossing, -1 if there is no crossing, and 0 if two points from
     * different edges are the same. Returns 0 or -1 if either edge is
     * degenerate. As a side effect, it saves vertex D to be used as the next
     * vertex C.
     */
    public int robustCrossing(S2Point d) {
      // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
      // all be oriented the same way (CW or CCW). We keep the orientation
      // of ACB as part of our state. When each new point D arrives, we
      // compute the orientation of BDA and check whether it matches ACB.
      // This checks whether the points C and D are on opposite sides of the
      // great circle through AB.

      // Recall that triageCCW is invariant with respect to rotating its
      // arguments, i.e. ABD has the same orientation as BDA.
      int bda = S2.triageCCW(aCrossB.dotProd(d));
      if (this.acb == -bda && bda != 0) {
        // The most common case -- triangles have opposite orientations.  Save the
        // current vertex D as the next vertex C, and also save the orientation of
        // the new triangle ACB (which is opposite to the current triangle BDA).
        this.c = d;
        this.acb = -bda;
        return -1;
      }
      this.bdaReturn = bda;
      return robustCrossingInternal(d);
    }

    /**
     * This method is equivalent to the S2EdgeUtil.edgeOrVertexCrossing() method
     * defined below. It is similar to robustCrossing, but handles cases where
     * two vertices are identical in a way that makes it easy to implement
     * point-in-polygon containment tests.
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
      // RobustCCW is very expensive, so we avoid calling it if at all possible.
      // First eliminate the cases where two vertices are equal.
      if (a.equalsPoint(c) || a.equalsPoint(d) || b.equalsPoint(c) || b.equalsPoint(d)) {
        return 0;
      }

      // At this point, a very common situation is that A,B,C,D are four points on a line such that
      // AB does not overlap CD.  (For example, this happens when a line or curve is sampled finely,
      // or when geometry is constructed by computing the union of S2CellIds.)  Most of the time, we
      // can determine that AB and CD do not intersect by computing the two outward-facing tangents
      // at A and B (parallel to AB) and testing whether AB and CD are on opposite sides of the
      // plane perpendicular to one of these tangents.  This is moderately expensive but still much
      // cheaper than S2.expensiveCCW().
      if (!hasTangents) {
        S2Point norm = S2Point.normalize(S2.robustCrossProd(a, b));
        aTangent = S2Point.crossProd(a, norm);
        bTangent = S2Point.crossProd(norm, b);
        hasTangents = true;
      }

      // The error in robustCrossProd() is insignificant.  The maximum error in the call to
      // crossProd() (i.e., the maximum norm of the error vector) is
      // (0.5 + 1/sqrt(3)) * S2.DBL_EPSILON.  The maximum error in each call to dotProd() below is
      // S2.DBL_EPSILON.  (There is also a small relative error term that is insignificant because
      // we are comparing the result against a constant that is very close to zero.)
      final double kError = (1.5 + 1 / Math.sqrt(3)) * S2.DBL_EPSILON;
      if ((c.dotProd(aTangent) > kError && d.dotProd(aTangent) > kError)
          || (c.dotProd(bTangent) > kError && d.dotProd(bTangent) > kError)) {
        return -1;
      }

      // Otherwise it's time to break out the big guns.
      if (acb == 0) {
        acb = -S2.expensiveCCW(a, b, c);
      }
      if (bdaReturn == 0) {
        bdaReturn = S2.expensiveCCW(a, b, d);
      }
      if (bdaReturn != acb) {
        return -1;
      }

      S2Point cCrossD = S2Point.crossProd(c, d);
      int cbd = -S2.robustCCW(c, d, b, cCrossD);
      if (cbd != acb) {
        return -1;
      }

      int dac = S2.robustCCW(c, d, a, cCrossD);
      return (dac == acb) ? 1 : -1;
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

    public RectBounder() {
    }

    /**
     * This method is called to add each vertex to the chain. This method is much faster than
     * {@link #addPoint(S2Point)}, since converting S2LatLng to an S2Point is much faster than the
     * other way around..
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
      // assert (S2.isUnitLength(b));
      if (builder.isEmpty()) {
        builder.addPoint(bLatLng);
      } else {
        // First compute the cross product N = A x B robustly.  This is the normal
        // to the great circle through A and B.  We don't use S2.RobustCrossProd()
        // since that method returns an arbitrary vector orthogonal to A if the two
        // vectors are proportional, and we want the zero vector in that case.
        // N = 2 * (A x B)
        S2Point n = S2Point.crossProd(S2Point.sub(a, b), S2Point.add(a, b));

        // The relative error in N gets large as its norm gets very small (i.e.,
        // when the two points are nearly identical or antipodal).  We handle this
        // by choosing a maximum allowable error, and if the error is greater than
        // this we fall back to a different technique.  Since it turns out that
        // the other sources of error add up to at most 1.16 * DBL_EPSILON, and it
        // is desirable to have the total error be a multiple of DBL_EPSILON, we
        // have chosen the maximum error threshold here to be 3.84 * DBL_EPSILON.
        // It is possible to show that the error is less than this when
        //
        //   n.norm() >= 8 * sqrt(3) / (3.84 - 0.5 - sqrt(3)) * DBL_EPSILON
        //            = 1.91346e-15 (about 8.618 * DBL_EPSILON)
        double nNorm = n.norm();
        if (nNorm < 1.91346e-15) {
          // A and B are either nearly identical or nearly antipodal (to within
          // 4.309 * DBL_EPSILON, or about 6 nanometers on the earth's surface).
          if (a.dotProd(b) < 0) {
            // The two points are nearly antipodal.  The easiest solution is to
            // assume that the edge between A and B could go in any direction
            // around the sphere.
            builder.setFull();
          } else {
            // The two points are nearly identical (to within 4.309 * DBL_EPSILON).
            // In this case we can just use the bounding rectangle of the points,
            // since after the expansion done by GetBound() this rectangle is
            // guaranteed to include the (lat,lng) values of all points along AB.
            builder.union(S2LatLngRect.fromPointPair(aLatLng, bLatLng));
          }
        } else {
          // Compute the longitude range spanned by AB.
          lngAB.initFromPointPair(aLatLng.lng().radians(), bLatLng.lng().radians());
          if (lngAB.getLength() >= S2.M_PI - 2 * S2.DBL_EPSILON) {
            // The points lie on nearly opposite lines of longitude to within the
            // maximum error of the calculation.  (Note that this test relies on
            // the fact that M_PI is slightly less than the true value of Pi, and
            // that representable values near M_PI are 2 * DBL_EPSILON apart.)
            // The easiest solution is to assume that AB could go on either side
            // of the pole.
            lngAB.setFull();
          }

          // Next we compute the latitude range spanned by the edge AB.  We start
          // with the range spanning the two endpoints of the edge:
          latAB.initFromPointPair(aLatLng.lat().radians(), bLatLng.lat().radians());

          // This is the desired range unless the edge AB crosses the plane
          // through N and the Z-axis (which is where the great circle through A
          // and B attains its minimum and maximum latitudes).  To test whether AB
          // crosses this plane, we compute a vector M perpendicular to this
          // plane and then project A and B onto it.
          S2Point m = S2Point.crossProd(n, S2Point.Z_POS);
          double mDotA = m.dotProd(a);
          double mDotB = m.dotProd(b);

          // We want to test the signs of "mDotA" and "mDotB", so we need to bound
          // the error in these calculations.  It is possible to show that the
          // total error is bounded by
          //
          //  (1 + sqrt(3)) * DBL_EPSILON * nNorm + 8 * sqrt(3) * (DBL_EPSILON**2)
          //    = 6.06638e-16 * nNorm + 6.83174e-31

          double mError = 6.06638e-16 * nNorm + 6.83174e-31;
          if (mDotA * mDotB < 0 || Math.abs(mDotA) <= mError || Math.abs(mDotB) <= mError) {
            // Minimum/maximum latitude *may* occur in the edge interior.
            //
            // The maximum latitude is 90 degrees minus the latitude of N.  We
            // compute this directly using atan2 in order to get maximum accuracy
            // near the poles.
            //
            // Our goal is compute a bound that contains the computed latitudes of
            // all S2Points P that pass the point-in-polygon containment test.
            // There are three sources of error we need to consider:
            //  - the directional error in N (at most 3.84 * DBL_EPSILON)
            //  - converting N to a maximum latitude
            //  - computing the latitude of the test point P
            // The latter two sources of error are at most 0.955 * DBL_EPSILON
            // individually, but it is possible to show by a more complex analysis
            // that together they can add up to at most 1.16 * DBL_EPSILON, for a
            // total error of 5 * DBL_EPSILON.
            //
            // We add 3 * DBL_EPSILON to the bound here, and getBound() will pad
            // the bound by another 2 * DBL_EPSILON.
            double maxLat = Math.min(S2.M_PI_2, 3 * S2.DBL_EPSILON + Math.atan2(
                Math.sqrt(n.getX() * n.getX() + n.getY() * n.getY()), Math.abs(n.getZ())));

            // In order to get tight bounds when the two points are close together,
            // we also bound the min/max latitude relative to the latitudes of the
            // endpoints A and B.  First we compute the distance between A and B,
            // and then we compute the maximum change in latitude between any two
            // points along the great circle that are separated by this distance.
            // This gives us a latitude change "budget".  Some of this budget must
            // be spent getting from A to B; the remainder bounds the round-trip
            // distance (in latitude) from A or B to the min or max latitude
            // attained along the edge AB.
            double latBudget = 2 * Math.asin(0.5 * S2Point.sub(a, b).norm() * Math.sin(maxLat));
            double maxDelta = 0.5 * (latBudget - latAB.getLength()) + S2.DBL_EPSILON;

            // Test whether AB passes through the point of maximum latitude or
            // minimum latitude.  If the dot product(s) are small enough then the
            // result may be ambiguous.
            if (mDotA <= mError && mDotB >= -mError) {
              latAB.setHi(Math.min(maxLat, latAB.hi() + maxDelta));
            }
            if (mDotB <= mError && mDotA >= -mError) {
              latAB.setLo(Math.max(-maxLat, latAB.lo() - maxDelta));
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
     * This bound satisfies the guarantee made above, i.e. if the edge chain defines a loop, then
     * the bound contains the S2LatLng coordinates of all S2Points contained by the loop.
     */
    public S2LatLngRect getBound() {
      // To save time, we ignore numerical errors in the computed S2LatLngs while
      // accumulating the bounds and then account for them here.
      //
      // S2LatLng(S2Point) has a maximum error of 0.955 * S2.DBL_EPSILON in latitude.
      // In the worst case, we might have rounded "inwards" when computing the
      // bound and "outwards" when computing the latitude of a contained point P,
      // therefore we expand the latitude bounds by 2 * S2.DBL_EPSILON in each
      // direction.  (A more complex analysis shows that 1.5 * S2.DBL_EPSILON is
      // enough, but the expansion amount should be a multiple of S2.DBL_EPSILON in
      // order to avoid rounding errors during the expansion itself.)
      //
      // S2LatLng(S2Point) has a maximum error of S2.DBL_EPSILON in longitude, which
      // is simply the maximum rounding error for results in the range [-Pi, Pi].
      // This is true because the Gnu implementation of atan2() comes from the IBM
      // Accurate Mathematical Library, which implements correct rounding for this
      // intrinsic (i.e., it returns the infinite precision result rounded to the
      // nearest representable value, with ties rounded to even values).  This
      // implies that we don't need to expand the longitude bounds at all, since
      // we only guarantee that the bound contains the *rounded* latitudes of
      // contained points.  The *true* latitudes of contained points may lie up to
      // S2.DBL_EPSILON outside of the returned bound.

      S2LatLng expansion = S2LatLng.fromRadians(2 * S2.DBL_EPSILON, 0);
      return builder.build().expanded(expansion).polarClosure();
    }

    /**
     * Returns the maximum error in getBound() provided that the result does not include either
     * pole. It is only to be used for testing purposes (e.g., by passing it to
     * {@link S2LatLngRect#approxEquals}).
     */
    static S2LatLng maxErrorForTests() {
      // The maximum error in the latitude calculation is
      //    3.84 * DBL_EPSILON   for the robustCrossProd calculation
      //    0.96 * DBL_EPSILON   for the latitude() calculation
      //    5    * DBL_EPSILON   added by AddPoint/GetBound to compensate for error
      //    ------------------
      //    9.80 * DBL_EPSILON   maximum error in result
      //
      // The maximum error in the longitude calculation is DBL_EPSILON.  GetBound
      // does not do any expansion because this isn't necessary in order to
      // bound the *rounded* longitudes of contained points.
      return S2LatLng.fromRadians(10 * S2.DBL_EPSILON, 1 * S2.DBL_EPSILON);
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
     * that {@code L.contains(S)}, then
     * {@code expandForSubregions(RectBound(L)).contains(RectBound(S))}.
     */
    static S2LatLngRect expandForSubregions(S2LatLngRect bound) {
      // Empty bounds don't need expansion.
      if (bound.isEmpty()) {
        return bound;
      }

      // First we need to check whether the bound B contains any nearly-antipodal
      // points (to within 4.309 * S2.DBL_EPSILON).  If so then we need to return
      // S2LatLngRect.full(), since the subregion might have an edge between two
      // such points, and addPoint() returns full() for such edges.  Note that
      // this can happen even if B is not full(); for example, consider a loop
      // that defines a 10km strip straddling the equator extending from
      // longitudes -100 to +100 degrees.
      //
      // It is easy to check whether B contains any antipodal points, but checking
      // for nearly-antipodal points is trickier.  Essentially we consider the
      // original bound B and its reflection through the origin B', and then test
      // whether the minimum distance between B and B' is less than 4.309 * DBL_EPSILON.

      // "lngGap" is a lower bound on the longitudinal distance between B and its
      // reflection B'.  (2.5 * S2.DBL_EPSILON is the maximum combined error of the
      // endpoint longitude calculations and the GetLength() call.)
      double lngGap = Math.max(0.0, S2.M_PI - bound.lng().getLength() - 2.5 * S2.DBL_EPSILON);

      // "minAbsLat" is the minimum distance from B to the equator (if zero or
      // negative, then B straddles the equator).
      double minAbsLat = Math.max(bound.lat().lo(), -bound.lat().hi());

      // "latGap1" and "latGap2" measure the minimum distance from B to the
      // south and north poles respectively.
      double latGap1 = S2.M_PI_2 + bound.lat().lo();
      double latGap2 = S2.M_PI_2 - bound.lat().hi();

      if (minAbsLat >= 0) {
        // The bound B does not straddle the equator.  In this case the minimum
        // distance is between one endpoint of the latitude edge in B closest to
        // the equator and the other endpoint of that edge in B'.  The latitude
        // distance between these two points is 2*minAbsLat, and the longitude
        // distance is lngGap.  We could compute the distance exactly using the
        // Haversine formula, but then we would need to bound the errors in that
        // calculation.  Since we only need accuracy when the distance is very
        // small (close to 4.309 * S2.DBL_EPSILON), we substitute the Euclidean
        // distance instead.  This gives us a right triangle XYZ with two edges of
        // length x = 2*minAbsLat and y ~= lngGap.  The desired distance is the
        // length of the third edge "z", and we have
        //
        //         z  ~=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
        //
        // Therefore the region may contain nearly antipodal points only if
        //
        //  2*minAbsLat + lngGap  <  sqrt(2) * 4.309 * S2.DBL_EPSILON
        //                           ~= 1.354e-15
        //
        // Note that because the given bound B is conservative, "minAbsLat" and
        // "lngGap" are both lower bounds on their true values so we do not need
        // to make any adjustments for their errors.
        if (2 * minAbsLat + lngGap < 1.354e-15) {
          return S2LatLngRect.full();
        }
      } else if (lngGap >= S2.M_PI_2) {
        // B spans at most Pi/2 in longitude.  The minimum distance is always
        // between one corner of B and the diagonally opposite corner of B'.  We
        // use the same distance approximation that we used above; in this case
        // we have an obtuse triangle XYZ with two edges of length x = latGap1
        // and y = latGap2, and angle Z >= Pi/2 between them.  We then have
        //
        //         z  >=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
        //
        // Unlike the case above, "latGap1" and "latGap2" are not lower bounds
        // (because of the extra addition operation, and because M_PI_2 is not
        // exactly equal to Pi/2); they can exceed their true values by up to
        // 0.75 * S2.DBL_EPSILON.  Putting this all together, the region may
        // contain nearly antipodal points only if
        //
        //   latGap1 + latGap2  <  (sqrt(2) * 4.309 + 1.5) * S2.DBL_EPSILON
        //                        ~= 1.687e-15
        if (latGap1 + latGap2 < 1.687e-15) {
          return S2LatLngRect.full();
        }
      } else {
        // Otherwise we know that (1) the bound straddles the equator and (2) its
        // width in longitude is at least Pi/2.  In this case the minimum
        // distance can occur either between a corner of B and the diagonally
        // opposite corner of B' (as in the case above), or between a corner of B
        // and the opposite longitudinal edge reflected in B'.  It is sufficient
        // to only consider the corner-edge case, since this distance is also a
        // lower bound on the corner-corner distance when that case applies.

        // Consider the spherical triangle XYZ where X is a corner of B with
        // minimum absolute latitude, Y is the closest pole to X, and Z is the
        // point closest to X on the opposite longitudinal edge of B'.  This is a
        // right triangle (Z = Pi/2), and from the spherical law of sines we have
        //
        //     sin(z) / sin(Z)  =  sin(y) / sin(Y)
        //     sin(maxLatGap) / 1  =  sin(d_min) / sin(lngGap)
        //     sin(d_min)  =  sin(maxLatGap) * sin(lngGap)
        //
        // where "maxLatGap" = max(latGap1, latGap2) and "d_min" is the
        // desired minimum distance.  Now using the facts that sin(t) >= (2/Pi)*t
        // for 0 <= t <= Pi/2, that we only need an accurate approximation when
        // at least one of "maxLatGap" or "lngGap" is extremely small (in
        // which case sin(t) ~= t), and recalling that "maxLatGap" has an error
        // of up to 0.75 * S2.DBL_EPSILON, we want to test whether
        //
        //   maxLatGap * lngGap  <  (4.309 + 0.75) * (Pi/2) * S2.DBL_EPSILON
        //                          ~= 1.765e-15
        if (Math.max(latGap1, latGap2) * lngGap < 1.765e-15) {
          return S2LatLngRect.full();
        }
      }

      // Next we need to check whether the subregion might contain any edges that
      // span (M_PI - 2 * S2.DBL_EPSILON) radians or more in longitude, since AddPoint
      // sets the longitude bound to Full() in that case.  This corresponds to
      // testing whether (lngGap <= 0) in "lng_expansion" below.

      // Otherwise, the maximum latitude error in AddPoint is 4.8 * S2.DBL_EPSILON.
      // In the worst case, the errors when computing the latitude bound for a
      // subregion could go in the opposite direction as the errors when computing
      // the bound for the original region, so we need to double this value.
      // (More analysis shows that it's okay to round down to a multiple of
      // S2.DBL_EPSILON.)
      //
      // For longitude, we rely on the fact that atan2 is correctly rounded and
      // therefore no additional bounds expansion is necessary.

      double latExpansion = 9 * S2.DBL_EPSILON;
      double lngExpansion = (lngGap <= 0) ? S2.M_PI : 0;
      return bound.expanded(S2LatLng.fromRadians(latExpansion, lngExpansion)).polarClosure();
    }
  }

  /**
   * The purpose of this class is to find edges that intersect a given XYZ
   * bounding box. It can be used as an efficient rejection test when attempting to
   * find edges that intersect a given region. It accepts a vertex chain v0, v1,
   * v2, ... and returns a boolean value indicating whether each edge intersects
   * the specified bounding box.
   *
   * <p>We use XYZ intervals instead of something like longitude intervals because
   * it is cheap to collect from S2Point lists and any slicing strategy should
   * give essentially equivalent results.  See S2Loop for an example of use.
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
      xmin = Math.min(xmin, Math.min(to.x, from.x));
      ymin = Math.min(ymin, Math.min(to.y, from.y));
      zmin = Math.min(zmin, Math.min(to.z, from.z));
      xmax = Math.max(xmax, Math.max(to.x, from.x));
      ymax = Math.max(ymax, Math.max(to.y, from.y));
      zmax = Math.max(zmax, Math.max(to.z, from.z));

      // Because our arcs are really geodesics on the surface of the earth
      // an edge can have intermediate points outside the xyz bounds implicit
      // in the end points.  Based on the length of the arc we compute a
      // generous bound for the maximum amount of deformation.  For small edges
      // it will be very small but for some large arcs (ie. from (1N,90W) to
      // (1N,90E) the path can be wildly deformed.  I did a bunch of
      // experiments with geodesics to get safe bounds for the deformation.
      double approxArcLen =
          Math.abs(from.x - to.x) + Math.abs(from.y - to.y) + Math.abs(from.z - to.z);
      if (approxArcLen < 0.025) { // less than 2 degrees
        maxDeformation = Math.max(maxDeformation, approxArcLen * 0.0025);
      } else if (approxArcLen < 1.0) { // less than 90 degrees
        maxDeformation = Math.max(maxDeformation, approxArcLen * 0.11);
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
     * Returns true if the edge going from the last point to this point passes
     * through the pruner bounding box, otherwise returns false.  So the
     * method returns false if we are certain there is no intersection, but it
     * may return true when there turns out to be no intersection.
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
   * The purpose of this class is to find edges that intersect a given longitude
   * interval. It can be used as an efficient rejection test when attempting to
   * find edges that intersect a given region. It accepts a vertex chain v0, v1,
   * v2, ... and returns a boolean value indicating whether each edge intersects
   * the specified longitude interval.
   *
   * <p>This class is not currently used as the XYZPruner is preferred for
   * S2Loop, but this should be usable in similar circumstances.  Be wary
   * of the cost of atan2() in conversions from S2Point to longitude!
   */
  public static class LongitudePruner {
    // The interval to be tested against.
    private S1Interval interval;

    // The longitude of the next v0.
    private double lng0;

    /**
     *'interval' is the longitude interval to be tested against, and 'v0' is
     * the first vertex of edge chain.
     */
    public LongitudePruner(S1Interval interval, S2Point v0) {
      this.interval = interval;
      this.lng0 = S2LatLng.longitude(v0).radians();
    }

    /**
     * Returns true if the edge (v0, v1) intersects the given longitude
     * interval, and then saves 'v1' to be used as the next 'v0'.
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
    // There are 6 possible edge orderings at a shared vertex (all
    // of these orderings are circular, i.e. abcd == bcda):
    //
    //  (1) a2 b2 b0 a0: A contains B
    //  (2) a2 a0 b0 b2: B contains A
    //  (3) a2 a0 b2 b0: A and B are disjoint
    //  (4) a2 b0 a0 b2: A and B intersect in one wedge
    //  (5) a2 b2 a0 b0: A and B intersect in one wedge
    //  (6) a2 b0 b2 a0: A and B intersect in two wedges
    //
    // We do not distinguish between 4, 5, and 6.
    // We pay extra attention when some of the edges overlap.  When edges
    // overlap, several of these orderings can be satisfied, and we take
    // the most specific.
    if (a0.equalsPoint(b0) && a2.equalsPoint(b2)) {
      return WedgeRelation.WEDGE_EQUALS;
    }

    if (S2.orderedCCW(a0, a2, b2, ab1)) {
      // The cases with this vertex ordering are 1, 5, and 6,
      // although case 2 is also possible if a2 == b2.
      if (S2.orderedCCW(b2, b0, a0, ab1)) {
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
    if (S2.orderedCCW(a0, b0, b2, ab1)) {
      return WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED;
    }
    if (S2.orderedCCW(a0, b0, a2, ab1)) {
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
   * <p>Implementations compare two *non-empty* wedges that share the same middle vertex:
   * A=(a0, ab1, a2) and B=(b0, ab1, b2).
   *
   * <p>All wedge processors require that a0 != a2 and b0 != b2. Other degenerate cases (such as
   * a0 == b2) are handled as expected. The parameter "ab1" denotes the common vertex a1 == b1.
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
   * Returns true if wedge A contains wedge B. Equivalent to but faster than
   * {@code getWedgeRelation() == WEDGE_PROPERLY_CONTAINS || WEDGE_EQUALS}.
   */
  public static class WedgeContains implements WedgeProcessor {
    /**
     * Given two edge chains, this function returns +1 if the region to the left of A contains the
     * region to the left of B, and 0 otherwise.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // For A to contain B (where each loop interior is defined to be its left
      // side), the CCW edge order around ab1 must be a2 b2 b0 a0. We split
      // this test into two parts that test three vertices each.
      return S2.orderedCCW(a2, b2, b0, ab1) && S2.orderedCCW(b0, a0, a2, ab1) ? 1 : 0;
    }
  }

  /**
   * Returns true if wedge A intersects wedge B.  Equivalent to but faster than
   * {@code getWedgeRelation() != WEDGE_IS_DISJOINT}.
   */
  public static class WedgeIntersects implements WedgeProcessor {
    /**
     * Given two edge chains (see WedgeRelation above), this function returns -1
     * if the region to the left of A intersects the region to the left of B,
     * and 0 otherwise. Note that regions are defined such that points along a
     * boundary are contained by one side or the other, not both. So for
     * example, if A,B,C are distinct points ordered CCW around a vertex O, then
     * the wedges BOA, AOC, and COB do not intersect.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // For A not to intersect B (where each loop interior is defined to be
      // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
      // Note that it's important to write these conditions as negatives
      // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
      // results when two vertices are the same.
      return (S2.orderedCCW(a0, b2, b0, ab1) && S2.orderedCCW(b0, a2, a0, ab1) ? 0 : -1);
    }
  }

  public static class WedgeContainsOrIntersects implements WedgeProcessor {
    /**
     * Given two edge chains (see WedgeRelation above), this function returns +1
     * if A contains B, 0 if A and B are disjoint, and -1 if A intersects but
     * does not contain B.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // This is similar to WedgeContainsOrCrosses, except that we want to
      // distinguish cases (1) [A contains B], (3) [A and B are disjoint],
      // and (2,4,5,6) [A intersects but does not contain B].

      if (S2.orderedCCW(a0, a2, b2, ab1)) {
        // We are in case 1, 5, or 6, or case 2 if a2 == b2.
        return S2.orderedCCW(b2, b0, a0, ab1) ? 1 : -1; // Case 1 vs. 2,5,6.
      }
      // We are in cases 2, 3, or 4.
      if (!S2.orderedCCW(a2, b0, b2, ab1)) {
        return 0; // Case 3.
      }

      // We are in case 2 or 4, or case 3 if a2 == b0.
      return (a2.equalsPoint(b0)) ? 0 : -1; // Case 3 vs. 2,4.
    }
  }

  public static class WedgeContainsOrCrosses implements WedgeProcessor {
    /**
     * Given two edge chains (see WedgeRelation above), this function returns +1
     * if A contains B, 0 if B contains A or the two wedges do not intersect,
     * and -1 if the edge chains A and B cross each other (i.e. if A intersects
     * both the interior and exterior of the region to the left of B). In
     * degenerate cases where more than one of these conditions is satisfied,
     * the maximum possible result is returned. For example, if A == B then the
     * result is +1.
     */
    @Override
    public int test(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2) {
      // There are 6 possible edge orderings at a shared vertex (all
      // of these orderings are circular, i.e. abcd == bcda):
      //
      // (1) a2 b2 b0 a0: A contains B
      // (2) a2 a0 b0 b2: B contains A
      // (3) a2 a0 b2 b0: A and B are disjoint
      // (4) a2 b0 a0 b2: A and B intersect in one wedge
      // (5) a2 b2 a0 b0: A and B intersect in one wedge
      // (6) a2 b0 b2 a0: A and B intersect in two wedges
      //
      // In cases (4-6), the boundaries of A and B cross (i.e. the boundary
      // of A intersects the interior and exterior of B and vice versa).
      // Thus we want to distinguish cases (1), (2-3), and (4-6).
      //
      // Note that the vertices may satisfy more than one of the edge
      // orderings above if two or more vertices are the same. The tests
      // below are written so that we take the most favorable
      // interpretation, i.e. preferring (1) over (2-3) over (4-6). In
      // particular note that if orderedCCW(a,b,c,o) returns true, it may be
      // possible that orderedCCW(c,b,a,o) is also true (if a == b or b == c).

      if (S2.orderedCCW(a0, a2, b2, ab1)) {
        // The cases with this vertex ordering are 1, 5, and 6,
        // although case 2 is also possible if a2 == b2.
        if (S2.orderedCCW(b2, b0, a0, ab1)) {
          return 1; // Case 1 (A contains B)
        }

        // We are in case 5 or 6, or case 2 if a2 == b2.
        return (a2.equalsPoint(b2)) ? 0 : -1; // Case 2 vs. 5,6.
      }
      // We are in case 2, 3, or 4.
      return S2.orderedCCW(a0, b0, a2, ab1) ? 0 : -1; // Case 2,3 vs. 4.
    }
  }

  /**
   * FaceSegment represents an edge AB clipped to an S2 cube face.  It is represented by a face
   * index and a pair of (u,v) coordinates.
   */
  static final class FaceSegment {
    final int face;
    final R2Vector a;
    final R2Vector b;
    public FaceSegment(int face, R2Vector a, R2Vector b) {
      this.face = face;
      this.a = a;
      this.b = b;
    }
  }

  // The three functions below all compare a sum (u + v) to a third value w.
  // They are implemented in such a way that they produce an exact result even
  // though all calculations are done with ordinary floating-point operations.
  // Here are the principles on which these functions are based:
  //
  // A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
  //
  // B. If u + v < w in exact arithmetic, then at least one of the following
  //    expressions is true in floating-point:
  //      u + v < w
  //      u < w - v
  //      v < w - u
  //
  //   Proof: By rearranging terms and substituting ">" for "<", we can assume
  //   that all values are non-negative.  Now clearly "w" is not the smallest
  //   value, so assume that "u" is the smallest.  We want to show that
  //   u < w - v in floating-point.  If v >= w/2, the calculation of w - v is
  //   exact since the result is smaller in magnitude than either input value,
  //   so the result holds.  Otherwise we have u <= v < w/2 and w - v >= w/2
  //   (even in floating point), so the result also holds.

  /** Returns true if u + v == w exactly. */
  static boolean sumEquals(double u, double v, double w) {
    return (u + v == w) && (u == w - v) && (v == w - u);
  }

  /**
   * Returns true if a given directed line L intersects the cube face F. The line L is defined by
   * its normal N in the (u,v,w) coordinates of F.
   */
  static boolean intersectsFace(S2Point n) {
    // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot
    // products of N with the four corner vertices (-1,-1,1), (1,-1,1), (1,1,1),
    // and (-1,1,1) do not all have the same sign.  This is true exactly when
    // |Nu| + |Nv| >= |Nw|.  The code below evaluates this expression exactly
    // (see comments above).
    double u = Math.abs(n.x);
    double v = Math.abs(n.y);
    double w = Math.abs(n.z);
    // We only need to consider the cases where u or v is the smallest value,
    // since if w is the smallest then both expressions below will have a
    // positive LHS and a negative RHS.
    return (v >= w - u) && (u >= w - v);
  }

  /**
   * Given a directed line L intersecting a cube face F, return true if L intersects two opposite
   * edges of F (including the case where L passes exactly through a corner vertex of F). The line L
   * is defined by its normal N in the (u,v,w) coordinates of F.
   */
  static boolean intersectsOppositeEdges(S2Point n) {
    // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if
    // and only exactly two of the corner vertices lie on each side of L.  This
    // is true exactly when ||Nu| - |Nv|| >= |Nw|.  The code below evaluates this
    // expression exactly (see comments above).
    double u = Math.abs(n.x);
    double v = Math.abs(n.y);
    double w = Math.abs(n.z);
    // If w is the smallest, the following line returns an exact result.
    if (Math.abs(u - v) != w) {
      return Math.abs(u - v) >= w;
    }
    // Otherwise u - v = w exactly, or w is not the smallest value.  In either
    // case the following line returns the correct result.
    return (u >= v) ? (u - w >= v) : (v - w >= u);
  }

  /**
   * Given cube face F and a directed line L (represented by its CCW normal N in the (u,v,w)
   * coordinates of F), compute the axis of the cube face edge where L exits the face: return 0 if L
   * exits through the u=-1 or u=+1 edge, and 1 if L exits through the v=-1 or v=+1 edge. Either
   * result is acceptable if L exits exactly through a corner vertex of the cube face.
   */
  static int getExitAxis(S2Point n) {
    // assert (intersectsFace(n));
    if (intersectsOppositeEdges(n)) {
      // The line passes through opposite edges of the face.
      // It exits through the v=+1 or v=-1 edge if the u-component of N has a
      // larger absolute magnitude than the v-component.
      return (Math.abs(n.x) >= Math.abs(n.y)) ? 1 : 0;
    } else {
      // The line passes through two adjacent edges of the face.
      // It exits the v=+1 or v=-1 edge if an even number of the components of N
      // are negative.  We test this using signbit() rather than multiplication
      // to avoid the possibility of underflow.
      // assert(n.x != 0 && n.y != 0  && n.z != 0);
      return ((n.x < 0) ^ (n.y < 0) ^ (n.z < 0)) ? 0 : 1;
    }
  }

  /**
   * Given a cube face F, a directed line L (represented by its CCW normal N in the (u,v,w)
   * coordinates of F), and result of getExitAxis(N), return the (u,v) coordinates of the point
   * where L exits the cube face.
   */
  static R2Vector getExitPoint(S2Point n, int axis) {
    if (axis == 0) {
      double u = (n.y > 0) ? 1.0 : -1.0;
      return new R2Vector(u, (-u * n.x - n.z) / n.y);
    } else {
      double v = (n.x < 0) ? 1.0 : -1.0;
      return new R2Vector((-v * n.y - n.z) / n.x, v);
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
    // Fast path: if the origin is sufficiently far inside the face, it is
    // always safe to use it.
    final double kMaxSafeUVCoord = 1 - FACE_CLIP_ERROR_UV_COORD;
    if (Math.max(Math.abs(aUv.x), Math.abs(aUv.y)) <= kMaxSafeUVCoord) {
      return face;
    }

    // Otherwise check whether the normal AB even intersects this face.
    S2Point n = S2Projections.faceXyzToUvw(face, ab);
    if (intersectsFace(n)) {
      // Check whether the point where the line AB exits this face is on the
      // wrong side of A (by more than the acceptable error tolerance).
      S2Point exit = S2Projections.faceUvToXyz(face, getExitPoint(n, getExitAxis(n)));
      S2Point aTangent = S2Point.crossProd(S2Point.normalize(ab), a);
      if (S2Point.sub(exit, a).dotProd(aTangent) >= -FACE_CLIP_ERROR_RADIANS) {
        return face;  // We can use the given face.
      }
    }

    // Otherwise we reproject A to the nearest adjacent face.  (If line AB does
    // not pass through a given face, it must pass through all adjacent faces.)
    if (Math.abs(aUv.x) >= Math.abs(aUv.y)) {
      face = S2Projections.getUVWFace(face, 0 /*U axis*/, aUv.x > 0 ? 1 : 0);
    } else {
      face = S2Projections.getUVWFace(face, 1 /*V axis*/, aUv.y > 0 ? 1 : 0);
    }
    // assert(intersectsFace(S2Projections.faceXyzToUvw(face, ab)));
    S2Projections.validFaceXyzToUv(face, a, aUv);
    aUv.set(
        Math.max(-1.0, Math.min(1.0, aUv.x)),
        Math.max(-1.0, Math.min(1.0, aUv.y)));
    return face;
  }

  /**
   * Return the next face that should be visited by getFaceSegments, given that we have just visited
   * "face" and we are following the line AB (represented by its normal N in the (u,v,w) coordinates
   * of that face). The other arguments include the point where AB exits "face", the corresponding
   * exit axis, and the "target face" containing the destination point B.
   */
  static int getNextFace(int face, R2Vector exit, int axis, S2Point n, int targetFace) {
    // We return the face that is adjacent to the exit point along the given
    // axis.  If line AB exits *exactly* through a corner of the face, there are
    // two possible next faces.  If one is the "target face" containing B, then
    // we guarantee that we advance to that face directly.
    //
    // The three conditions below check that (1) AB exits approximately through
    // a corner, (2) the adjacent face along the non-exit axis is the target
    // face, and (3) AB exits *exactly* through the corner.  (The sumEquals()
    // code checks whether the dot product of (u,v,1) and "n" is exactly zero.)
    if (Math.abs(exit.get(1 - axis)) == 1
        && S2Projections.getUVWFace(face, 1 - axis, exit.get(1 - axis) > 0 ? 1 : 0) == targetFace
        && sumEquals(exit.x * n.x, exit.y * n.y, -n.z)) {
      return targetFace;
    }

    // Otherwise return the face that is adjacent to the exit point in the
    // direction of the exit axis.
    return S2Projections.getUVWFace(face, axis, exit.get(axis) > 0 ? 1 : 0);
  }

  /**
   * Subdivide the given edge AB at every point where it crosses the boundary between two S2 cube
   * faces and return the corresponding FaceSegments. The segments are returned in order from A
   * toward B. The input points must be unit length.
   *
   * <p>This method guarantees that the returned segments form a continuous path from A to B, and
   * that all vertices are within kFaceClipErrorUVDist of the line AB. All vertices lie within the
   * [-1,1]x[-1,1] cube face rectangles. The results are consistent with S2.RobustCCW(), i.e. the
   * edge is well-defined even its endpoints are antipodal.
   */
  // TODO(user): Extend the implementation of S2.RobustCrossProd so that this statement is true.
  static int getFaceSegments(S2Point a, S2Point b, FaceSegment[] segments) {
    // assert(S2.IsUnitLength(a));
    // assert(S2.IsUnitLength(b));

    int size = 0;

    // Fast path: both endpoints are on the same face.
    int aFace = S2Projections.xyzToFace(a);
    R2Vector aSeg = S2Projections.validFaceXyzToUv(aFace, a);
    int bFace = S2Projections.xyzToFace(b);
    R2Vector bSeg = S2Projections.validFaceXyzToUv(bFace, b);
    if (aFace == bFace) {
      segments[size++] = new FaceSegment(aFace, aSeg, bSeg);
    } else {
      // Starting at A, we follow AB from face to face until we reach the face
      // containing B.  The following code is designed to ensure that we always
      // reach B, even in the presence of numerical errors.
      //
      // First we compute the normal to the plane containing A and B.  This normal
      // becomes the ultimate definition of the line AB; it is used to resolve all
      // questions regarding where exactly the line goes.  Unfortunately due to
      // numerical errors, the line may not quite intersect the faces containing
      // the original endpoints.  We handle this by moving A and/or B slightly if
      // necessary so that they are on faces intersected by the line AB.
      S2Point ab = S2.robustCrossProd(a, b);
      aFace = moveOriginToValidFace(aFace, a, ab, aSeg);
      bFace = moveOriginToValidFace(bFace, b, S2Point.neg(ab), bSeg);

      // Now we simply follow AB from face to face until we reach B.
      R2Vector bSaved = bSeg;
      for (int face = aFace; face != bFace; ) {
        // Complete the current segment by finding the point where AB exits the
        // current face.
        S2Point n = S2Projections.faceXyzToUvw(face, ab);
        int exitAxis = getExitAxis(n);
        bSeg = getExitPoint(n, exitAxis);
        segments[size++] = new FaceSegment(face, aSeg, bSeg);

        // Compute the next face intersected by AB, and translate the exit point
        // of the current segment into the (u,v) coordinates of the next face.
        // This becomes the first point of the next segment.
        S2Point exitXyz = S2Projections.faceUvToXyz(face, bSeg);
        face = getNextFace(face, bSeg, exitAxis, n, bFace);
        S2Point exitUvw = S2Projections.faceXyzToUvw(face, exitXyz);
        aSeg = new R2Vector(exitUvw.x, exitUvw.y);
      }

      // Finish the last segment.
      segments[size++] = new FaceSegment(bFace, aSeg, bSaved);
    }

    return size;
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
      S2Point a, S2Point b, S2Point nScaled,
      S2Point aTangent, S2Point bTangent, double uvScale,
      R2Vector uv) {
    // assert(intersectsFace(nScaled));

    // Optimization: if B is within the safe region of the face, use it.
    final double kMaxSafeUVCoord = 1 - FACE_CLIP_ERROR_UV_COORD;
    if (b.z > 0) {
      uv.set(b.x / b.z, b.y / b.z);
      if (Math.max(Math.abs(uv.x), Math.abs(uv.y)) <= kMaxSafeUVCoord) {
        return 0;
      }
    }

    // Otherwise find the point B' where the line AB exits the face.
    R2Vector exit = getExitPoint(nScaled, getExitAxis(nScaled));
    uv.set(uvScale * exit.x, uvScale * exit.y);
    S2Point p = new S2Point(uv.x, uv.y, 1);

    // Determine if the exit point B' is contained within the segment.  We do this
    // by computing the dot products with two inward-facing tangent vectors at A
    // and B.  If either dot product is negative, we say that B' is on the "wrong
    // side" of that point.  As the point B' moves around the great circle AB past
    // the segment endpoint B, it is initially on the wrong side of B only; as it
    // moves further it is on the wrong side of both endpoints; and then it is on
    // the wrong side of A only.  If the exit point B' is on the wrong side of
    // either endpoint, we can't use it; instead the segment is clipped at the
    // original endpoint B.
    //
    // We reject the segment if the sum of the scores of the two endpoints is 3
    // or more.  Here is what that rule encodes:
    //  - If B' is on the wrong side of A, then the other clipped endpoint A'
    //    must be in the interior of AB (otherwise AB' would go the wrong way
    //    around the circle).  There is a similar rule for A'.
    //  - If B' is on the wrong side of either endpoint (and therefore we must
    //    use the original endpoint B instead), then it must be possible to
    //    project B onto this face (i.e., its w-coordinate must be positive).
    //    This rule is only necessary to handle certain zero-length edges (A=B).
    int score = 0;
    if (S2Point.sub(p, a).dotProd(aTangent) < 0) {
      score = 2;  // B' is on wrong side of A.
    } else if (S2Point.sub(p, b).dotProd(bTangent) < 0) {
      score = 1;  // B' is on wrong side of B.
    }
    if (score > 0) {  // B' is not in the interior of AB.
      if (b.z <= 0) {
        score = 3;    // B cannot be projected onto this face.
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
   */
  static boolean clipToPaddedFace(
      S2Point aXyz, S2Point bXyz, int face, double padding, R2Vector aUv, R2Vector bUv) {
    // assert (padding >= 0);
    // Fast path: both endpoints are on the given face.
    if (S2Projections.xyzToFace(aXyz) == face && S2Projections.xyzToFace(bXyz) == face) {
      S2Projections.validFaceXyzToUv(face, aXyz, aUv);
      S2Projections.validFaceXyzToUv(face, bXyz, bUv);
      return true;
    }

    // Convert everything into the (u,v,w) coordinates of the given face.  Note
    // that the cross product *must* be computed in the original (x,y,z)
    // coordinate system because RobustCrossProd (unlike the mathematical cross
    // product) can produce different results in different coordinate systems
    // when one argument is a linear multiple of the other, due to the use of
    // symbolic perturbations.
    S2Point n = S2Projections.faceXyzToUvw(face, S2.robustCrossProd(aXyz, bXyz));
    S2Point a = S2Projections.faceXyzToUvw(face, aXyz);
    S2Point b = S2Projections.faceXyzToUvw(face, bXyz);

    // Padding is handled by scaling the u- and v-components of the normal.
    // Letting R=1+padding, this means that when we compute the dot product of
    // the normal with a cube face vertex (such as (-1,-1,1)), we will actually
    // compute the dot product with the scaled vertex (-R,-R,1).  This allows
    // methods such as intersectsFace(), getExitAxis(), etc, to handle padding
    // with no further modifications.
    final double uvScale = 1 + padding;
    S2Point nScaled = new S2Point(uvScale * n.x, uvScale * n.y, n.z);
    if (!intersectsFace(nScaled)) {
      return false;
    }

    // TODO(user): This is a temporary hack until I rewrite S2.RobustCrossProd;
    // it avoids loss of precision in normalize() when the vector is so small
    // that it underflows.
    if (Math.max(Math.abs(n.x), Math.max(Math.abs(n.y), Math.abs(n.z))) < Math.scalb(1d, -511)) {
      n = S2Point.mul(n, Math.scalb(1d, 563));
    }  // END OF HACK

    n = S2Point.normalize(n);
    S2Point aTangent = S2Point.crossProd(n, a);
    S2Point bTangent = S2Point.crossProd(b, n);
    // As described above, if the sum of the scores from clipping the two
    // endpoints is 3 or more, then the segment does not intersect this face.
    int aScore = clipDestination(b, a, S2Point.neg(nScaled), bTangent, aTangent, uvScale, aUv);
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

    // Otherwise AB intersects "rect" if and only if all four vertices of "rect"
    // do not lie on the same side of the extended line AB.  We test this by
    // finding the two vertices of "rect" with minimum and maximum projections
    // onto the normal of AB, and computing their dot products with the edge
    // normal.
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
      double a0, double b0, R1Interval bound0,
      double a1, double b1, R1Interval bound1,
      boolean slopeNegative, R1Interval clip0) {
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
    // "slopeNegative" indicates which diagonal of the bounding box is spanned by AB: it
    // is false if AB has positive slope, and true if AB has negative slope.  This is
    // used to determine which interval endpoints need to be updated each time
    // the edge is clipped.
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
    // Compute the bounding rectangle of AB, clip it, and then extract the new
    // endpoints from the clipped bound.
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
   * same linear combination of a1 and b1.  This function makes the following guarantees:
   *
   * <ol>
   * <li>If x == a, then x1 = a1 (exactly).
   * <li>If x == b, then x1 = b1 (exactly).
   * <li>If a <= x <= b, then a1 <= x1 <= b1 (even if a1 == b1).
   * </ol>
   *
   * <p>Results are undefined if a==b.
   */
  static double interpolateDouble(double x, double a, double b, double a1, double b1) {
    // assertTrue(a != b);
    // To get results that are accurate near both A and B, we interpolate starting from the closer
    // of the two points.
    if (Math.abs(a - x) <= Math.abs(b - x)) {
      return a1 + (b1 - a1) * (x - a) / (b - a);
    } else {
      return b1 + (a1 - b1) * (x - b) / (a - b);
    }
  }

  /**
   * Given an edge AB and a face, return the (u,v) coordinates for the portion of AB that intersects
   * that face. This method guarantees that the clipped vertices lie within the [-1,1]x[-1,1] cube
   * face rectangle and are within kFaceClipErrorUVDist of the line AB, but the results may differ
   * from those produced by getFaceSegments. Returns false if AB does not intersect the given face.
   */
  static boolean clipToFace(S2Point a, S2Point b, int face, R2Vector aUv, R2Vector bUv) {
    return clipToPaddedFace(a, b, face, 0.0, aUv, bUv);
  }

  /**
   * Return true if edge AB crosses CD at a point that is interior to both
   * edges. Properties:
   *
   * <ul>
   * <li>simpleCrossing(b,a,c,d) == simpleCrossing(a,b,c,d)
   * <li>simpleCrossing(c,d,a,b) == simpleCrossing(a,b,c,d)
   * </ul>
   */
  public static boolean simpleCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    // We compute simpleCCW() for triangles ACB, CBD, BDA, and DAC. All
    // of these triangles need to have the same orientation (CW or CCW)
    // for an intersection to exist. Note that this is slightly more
    // restrictive than the corresponding definition for planar edges,
    // since we need to exclude pairs of line segments that would
    // otherwise "intersect" by crossing two antipodal points.

    S2Point ab = S2Point.crossProd(a, b);
    double acb = -(ab.dotProd(c));
    double bda = ab.dotProd(d);
    if (acb * bda <= 0) {
      return false;
    }

    S2Point cd = S2Point.crossProd(c, d);
    double cbd = -(cd.dotProd(b));
    double dac = cd.dotProd(a);
    return (acb * cbd > 0) && (acb * dac > 0);
  }

  /**
   * Like SimpleCrossing, except that points that lie exactly on a line are
   * arbitrarily classified as being on one side or the other (according to the
   * rules of S2.robustCCW). It returns +1 if there is a crossing, -1 if there
   * is no crossing, and 0 if any two vertices from different edges are the
   * same. Returns 0 or -1 if either edge is degenerate. Properties of
   * robustCrossing:
   *
   * <ul>
   * <li>robustCrossing(b,a,c,d) == robustCrossing(a,b,c,d)
   * <li>robustCrossing(c,d,a,b) == robustCrossing(a,b,c,d)
   * <li>robustCrossing(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
   * <li>robustCrossing(a,b,c,d) <= 0 if a==b or c==d
   * </ul>
   *
   * <p>Note that if you want to check an edge against a *chain* of other edges,
   * it is much more efficient to use an EdgeCrosser (above).
   */
  public static int robustCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    EdgeCrosser crosser = new EdgeCrosser(a, b, c);
    return crosser.robustCrossing(d);
  }

  /**
   * Given two edges AB and CD where at least two vertices are identical (i.e.
   * robustCrossing(a,b,c,d) == 0), this function defines whether the two edges
   * "cross" in a such a way that point-in-polygon containment tests can be
   * implemented by counting the number of edge crossings. The basic rule is
   * that a "crossing" occurs if AB is encountered after CD during a CCW sweep
   * around the shared vertex starting from a fixed reference point.
   *
   * <p>Note that according to this rule, if AB crosses CD then in general CD does
   * not cross AB. However, this leads to the correct result when counting
   * polygon edge crossings. For example, suppose that A,B,C are three
   * consecutive vertices of a CCW polygon. If we now consider the edge
   * crossings of a segment BP as P sweeps around B, the crossing number changes
   * parity exactly when BP crosses BA or BC.
   *
   * <p>Useful properties of VertexCrossing (VC):
   *
   * <ul>
   * <li>VC(a,a,c,d) == VC(a,b,c,c) == false
   * <li> VC(a,b,a,b) == VC(a,b,b,a) == true
   * <li>VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
   * <li>If exactly one of a,b equals one of c,d, then exactly one of VC(a,b,c,d) and
   * VC(c,d,a,b) is true
   * </ul>
   *
   * <p>It is an error to call this method with 4 distinct vertices.
   */
  public static boolean vertexCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    // If A == B or C == D there is no intersection. We need to check this
    // case first in case 3 or more input points are identical.
    if (a.equalsPoint(b) || c.equalsPoint(d)) {
      return false;
    }

    // If any other pair of vertices is equal, there is a crossing if and only
    // if orderedCCW() indicates that the edge AB is further CCW around the
    // shared vertex than the edge CD.
    if (a.equalsPoint(d)) {
      return S2.orderedCCW(S2.ortho(a), c, b, a);
    }
    if (b.equalsPoint(c)) {
      return S2.orderedCCW(S2.ortho(b), d, a, b);
    }
    if (a.equalsPoint(c)) {
      return S2.orderedCCW(S2.ortho(a), d, b, a);
    }
    if (b.equalsPoint(d)) {
      return S2.orderedCCW(S2.ortho(b), c, a, b);
    }

    // assert (false);
    return false;
  }

  /**
   * A convenience function that calls robustCrossing() to handle cases where
   * all four vertices are distinct, and VertexCrossing() to handle cases where
   * two or more vertices are the same. This defines a crossing function such
   * that point-in-polygon containment tests can be implemented by simply
   * counting edge crossings.
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
      // If the squared distance from x to y is less than dmin2, then replace
      // vmin by y and update dmin2 accordingly.
      double d2 = S2Point.minus(x, y).norm2();
      if (d2 < dmin2 || (d2 == dmin2 && y.lessThan(vmin))) {
        dmin2 = d2;
        vmin = y;
      }
    }
  }

  /**
   * Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
   */
  public static final boolean lenientCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    // assert (S2.isUnitLength(a));
    // assert (S2.isUnitLength(b));
    // assert (S2.isUnitLength(c));

    double acb = S2Point.scalarTripleProduct(b, a, c);
    if (Math.abs(acb) < MAX_DET_ERROR) {
      return true;
    }
    double bda = S2Point.scalarTripleProduct(a, b, d);
    if (Math.abs(bda) < MAX_DET_ERROR) {
      return true;
    }
    if (acb * bda < 0) {
      return false;
    }
    double cbd = S2Point.scalarTripleProduct(d, c, b);
    if (Math.abs(cbd) < MAX_DET_ERROR) {
      return true;
    }
    double dac = S2Point.scalarTripleProduct(c, d, a);
    if (Math.abs(dac) < MAX_DET_ERROR) {
      return true;
    }
    return (acb * cbd >= 0) && (acb * dac >= 0);
  }

  /**
   * Given two edges AB and CD such that robustCrossing() is true, return their
   * intersection point. Useful properties of getIntersection (GI):
   *
   * <ul>
   * <li>GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
   * <li>GI(c,d,a,b) == GI(a,b,c,d)
   * </ul>
   *
   * <p>The returned intersection point X is guaranteed to be close to the edges AB
   * and CD, but if the edges intersect at a very small angle then X may not be
   * close to the true mathematical intersection point P. See the description of
   * "DEFAULT_INTERSECTION_TOLERANCE" below for details.
   */
  public static S2Point getIntersection(S2Point a0, S2Point a1, S2Point b0, S2Point b1) {
    Preconditions.checkArgument(robustCrossing(a0, a1, b0, b1) > 0,
        "Input edges a0a1 and b0b1 must have a true robustCrossing.");

    // We use robustCrossProd() to get accurate results even when two endpoints
    // are close together, or when the two line segments are nearly parallel.
    S2Point aNorm = S2Point.normalize(S2.robustCrossProd(a0, a1));
    S2Point bNorm = S2Point.normalize(S2.robustCrossProd(b0, b1));
    S2Point x = S2Point.normalize(S2.robustCrossProd(aNorm, bNorm));

    // Make sure the intersection point is on the correct side of the sphere.
    // Since all vertices are unit length, and edges are less than 180 degrees,
    // (a0 + a1) and (b0 + b1) both have positive dot product with the
    // intersection point. We use the sum of all vertices to make sure that the
    // result is unchanged when the edges are reversed or exchanged.
    if (x.dotProd(S2Point.add(S2Point.add(a0, a1), S2Point.add(b0, b1))) < 0) {
      x = S2Point.neg(x);
    }

    // The calculation above is sufficient to ensure that "x" is within
    // DEFAULT_INTERSECTION_TOLERANCE of the great circles through (a0,a1) and
    // (b0,b1).
    // However, if these two great circles are very close to parallel, it is
    // possible that "x" does not lie between the endpoints of the given line
    // segments. In other words, "x" might be on the great circle through
    // (a0,a1) but outside the range covered by (a0,a1). In this case we do
    // additional clipping to ensure that it does.

    if (S2.orderedCCW(a0, x, a1, aNorm) && S2.orderedCCW(b0, x, b1, bNorm)) {
      return x;
    }

    // Find the acceptable endpoint closest to x and return it. An endpoint is
    // acceptable if it lies between the endpoints of the other line segment.
    CloserResult r = new CloserResult(10, x);
    if (S2.orderedCCW(b0, a0, b1, bNorm)) {
      r.replaceIfCloser(x, a0);
    }
    if (S2.orderedCCW(b0, a1, b1, bNorm)) {
      r.replaceIfCloser(x, a1);
    }
    if (S2.orderedCCW(a0, b0, a1, aNorm)) {
      r.replaceIfCloser(x, b0);
    }
    if (S2.orderedCCW(a0, b1, a1, aNorm)) {
      r.replaceIfCloser(x, b1);
    }
    return r.getVmin();
  }

  /**
   * Given a point X and an edge AB, return the distance ratio AX / (AX + BX).
   * If X happens to be on the line segment AB, this is the fraction "t" such
   * that X == Interpolate(A, B, t). Requires that A and B are distinct.
   */
  public static double getDistanceFraction(S2Point x, S2Point a0, S2Point a1) {
    Preconditions.checkArgument(!a0.equalsPoint(a1));
    double d0 = x.angle(a0);
    double d1 = x.angle(a1);
    return d0 / (d0 + d1);
  }

  /**
   * Return the minimum distance from X to any point on the edge AB. The result
   * is very accurate for small distances but may have some numerical error if
   * the distance is large (approximately Pi/2 or greater). The case A == B is
   * handled correctly.
   *
   * @throws IllegalArgumentException Thrown if the parameters are not all unit
   * length.
   */
  public static S1Angle getDistance(S2Point x, S2Point a, S2Point b) {
    Preconditions.checkArgument(S2.isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(S2.isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(S2.isUnitLength(b), "S2Point not normalized: %s", b);
    return S1Angle.radians(getDistanceRadians(x, a, b, S2.robustCrossProd(a, b)));
  }

  /**
   * Return the minimum of the distance from X to any point on edge AB and the given
   * {@code minDistance}.
   *
   * @throws IllegalArgumentException Thrown if the parameters are not all unit
   * length.
   */
  public static S1ChordAngle updateMinDistance(S2Point x, S2Point a, S2Point b,
      S1ChordAngle minDistance) {
    Preconditions.checkArgument(S2.isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(S2.isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(S2.isUnitLength(b), "S2Point not normalized: %s", b);

    // We divide the problem into two cases, based on whether the closest point
    // on AB is one of the two vertices (the "vertex case") or in the interior
    // (the "interior case").  Let C = A x B.  If X is in the spherical wedge
    // extending from A to B around the axis through C, then we are in the
    // interior case.  Otherwise we are in the vertex case.

    // Check whether we might be in the interior case. For this to be true, XAB and XBA must both be
    // acute angles. Checking this condition exactly is expensive, so instead we consider the 3D
    // Euclidian triangle ABX (which passes through the sphere's interior). As can be observed from
    // the law of spherical excess, the planar angles XAB and XBA are always less than the
    // corresponding spherical angles, so if we are in the interior case then both of these angles
    // must be acute.
    //
    // We check this by computing the squared edge lengths of the 3D Euclidean triangle ABX, and
    // testing acuteness using the law of cosines:
    //
    //                      max(XA^2, XB^2) < AB^2 + min(XA^2, XB^2)
    // or equivalently:     XA^2 + XB^2 < AB^2 + 2 * min(XA^2, XB^2)
    //
    double xa2 = x.getDistance2(a), xb2 = x.getDistance2(b), ab2 = a.getDistance2(b);
    double dist2 = min(xa2, xb2);
    if (xa2 + xb2 < ab2 + 2 * dist2) {
      // The minimum distance might be to a point on the edge interior. Let R be the closest point
      // to X that lies on the great circle through AB. Rather than computing the geodesic distance
      // along the surface of the sphere, instead we compute the "chord length", the 3D Euclidian
      // length of the line passing through the sphere's interior. If the squared chord length
      // exceeds minDistance.getLength2() then we can return "false" immediately.
      //
      // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q is the point X
      // projected onto the plane through the great circle AB.
      // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
      // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it is faster and the
      // corresponding distance on the Earth's surface is accurate to within 1% for distances up to
      // about 1800km.
      S2Point c = S2.robustCrossProd(a, b);
      double c2 = c.norm2();
      double xDotC = x.dotProd(c);
      double xDotC2 = xDotC * xDotC;
      if (xDotC2 >= c2 * minDistance.getLength2()) {
        // The closest point on the great circle AB is too far away.
        return minDistance;
      }
      // Otherwise we do the exact, more expensive test for the interior case.
      // This test is very likely to succeed because of the conservative planar test we did
      // initially.
      S2Point cx = S2Point.crossProd(c, x);
      if (a.dotProd(cx) < 0 && b.dotProd(cx) > 0) {
        // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
        // This calculation has good accuracy for all chord lengths since it is based on both the
        // dot product and cross product (rather than deriving one from the other). However, note
        // that the chord length representation itself loses accuracy as the angle approaches Pi.
        double qr = 1 - sqrt(cx.norm2() / c2);
        dist2 = (xDotC2 / c2) + (qr * qr);
      }
    }
    if (dist2 >= minDistance.getLength2()) {
      return minDistance;
    }
    return S1ChordAngle.fromLength2(dist2);
  }

  /**
   * A slightly more efficient version of getDistance() where the cross product
   * of the two endpoints has been precomputed. The cross product does not need
   * to be normalized, but should be computed using S2.robustCrossProd() for the
   * most accurate results.
    *
   * @throws IllegalArgumentException Thrown if the parameters are not all unit
   * length.
  */
  public static S1Angle getDistance(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    Preconditions.checkArgument(S2.isUnitLength(x), "S2Point not normalized: %s", x);
    Preconditions.checkArgument(S2.isUnitLength(a), "S2Point not normalized: %s", a);
    Preconditions.checkArgument(S2.isUnitLength(b), "S2Point not normalized: %s", b);
    return S1Angle.radians(getDistanceRadians(x, a, b, aCrossB));
  }

  /**
   * A more efficient version of getDistance() where the cross product of the
   * endpoints has been precomputed and the result is returned as a direct
   * radian measure rather than wrapping it in an S1Angle. This is the
   * recommended method for making large numbers of back-to-back edge distance
   * tests, since it allocates no objects. The inputs are assumed to be unit
   * length; results are undefined if they are not.
   */
  public static double getDistanceRadians(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    // There are three cases. If X is located in the spherical wedge defined by
    // A, B, and the axis A x B, then the closest point is on the segment AB.
    // Otherwise the closest point is either A or B; the dividing line between
    // these two cases is the great circle passing through (A x B) and the
    // midpoint of AB.
    if (S2.simpleCCW(aCrossB, a, x) && S2.simpleCCW(x, b, aCrossB)) {
      // The closest point to X lies on the segment AB. We compute the distance
      // to the corresponding great circle. The result is accurate for small
      // distances but not necessarily for large distances (approaching Pi/2).
      double sinDist = Math.abs(x.dotProd(aCrossB)) / aCrossB.norm();
      return Math.asin(Math.min(1.0, sinDist));
    }

    // Otherwise, the closest point is either A or B. The cheapest method is
    // just to compute the minimum of the two linear (as opposed to spherical)
    // distances and convert the result to an angle. Again, this method is
    // accurate for small but not large distances (approaching Pi).
    double linearDist2 = Math.min(diffMag2(x, a), diffMag2(x, b));
    return 2 * Math.asin(Math.min(1.0, 0.5 * Math.sqrt(linearDist2)));
  }

  /** Returns the squared distance from {@code a} to {@code b}. */
  private static final double diffMag2(S2Point a, S2Point b) {
    double dx = a.getX() - b.getX();
    double dy = a.getY() - b.getY();
    double dz = a.getZ() - b.getZ();
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * As {@link #getClosestPoint(S2Point, S2Point, S2Point)}, but faster if the cross product between
   * a and b has already been computed. All points must be unit length; results are undefined if
   * that is not the case.
   */
  public static S2Point getClosestPoint(S2Point x, S2Point a, S2Point b, S2Point aCrossB) {
    // assert (S2.isUnitLength(a));
    // assert (S2.isUnitLength(b));
    // assert (S2.isUnitLength(x));

    // Find the closest point to X along the great circle through AB.
    S2Point p = S2Point.sub(x, S2Point.mul(aCrossB, x.dotProd(aCrossB) / aCrossB.norm2()));

    // If this point is on the edge AB, then it's the closest point.
    if (S2.simpleCCW(aCrossB, a, p) && S2.simpleCCW(p, b, aCrossB)) {
      return S2Point.normalize(p);
    }

    // Otherwise, the closest point is either A or B.
    return x.getDistance2(a) <= x.getDistance2(b) ? a : b;
  }

  /**
   * Returns the point on edge AB closest to X. All points must be unit length; results are
   * undefined if that is not the case.
   */
  public static S2Point getClosestPoint(S2Point x, S2Point a, S2Point b) {
    return getClosestPoint(x, a, b, S2.robustCrossProd(a, b));
  }

  /**
   * A slightly more efficient version of InterpolateAtDistance() that can be
   * used when the distance AB is already known.  Requires that all vectors
   * have unit length.
   */
  public static S2Point interpolateAtDistance(S1Angle ax, S2Point a, S2Point b, S1Angle ab) {
    // assert S2.isUnitLength(a);
    // assert S2.isUnitLength(b);

    double axRadians = ax.radians();
    double abRadians = ab.radians();

    // The result X is some linear combination X = e*A + f*B of the input
    // points.  The fractions "e" and "f" can be derived by looking at the
    // components of this equation that are parallel and perpendicular to A.
    // Let E = e*A and F = f*B.  Then OEXF is a parallelogram.  You can obtain
    // the distance f = OF by considering the similar triangles produced by
    // dropping perpendiculars from the segments OF and OB to OA.
    double f = Math.sin(axRadians) / Math.sin(abRadians);

    // Form the dot product of the first equation with A to obtain
    // A.X = e*A.A + f*A.B.  Since A, B, and X are all unit vectors,
    // cos(ax) = e*1 + f*cos(ab), so
    double e = Math.cos(axRadians) - f * Math.cos(abRadians);

    // Mathematically speaking, if "a" and "b" are unit length then the result
    // is unit length as well.  But we normalize it anyway to prevent points
    // from drifting away from unit length when multiple interpolations are done
    // in succession (i.e. the result of one interpolation is fed into another).
    return S2Point.normalize(S2Point.add(S2Point.mul(a, e), S2Point.mul(b, f)));
  }

  /**
   * Like Interpolate(), except that the parameter "ax" represents the desired
   * distance from A to the result X rather than a fraction between 0 and 1.
   * Requires that {@code a} and {@code b} are unit length.
   */
  public static S2Point interpolateAtDistance(S1Angle ax, S2Point a, S2Point b) {
    return interpolateAtDistance(ax, a, b, new S1Angle(a, b));
  }

  /**
   * Return the point X along the line segment AB whose distance from A is the
   * given fraction "t" of the distance AB. Does NOT require that "t" be between
   * 0 and 1. Note that all distances are measured on the surface of the sphere,
   * so this is more complicated than just computing (1-t)*a + t*b and
   * normalizing the result.
   */
  public static S2Point interpolate(double t, S2Point a, S2Point b) {
    if (t == 0) {
      return a;
    }
    if (t == 1) {
      return b;
    }
    S1Angle ab = new S1Angle(a, b);
    return interpolateAtDistance(S1Angle.radians(t * ab.radians()), a, b, ab);
  }

  /** Constructor is private so that this class is never instantiated. */
  private S2EdgeUtil() {
  }
}
