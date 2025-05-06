/*
 * Copyright 2005 Google Inc.
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

import static java.lang.Math.E;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.atan2;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.scalb;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;

import com.google.common.base.Preconditions;
import com.google.errorprone.annotations.Immutable;

/**
 * The S2 class is simply a namespace for constants and static utility functions related to
 * spherical geometry, such as area calculations and edge intersection tests. The name "S2" is
 * derived from the mathematical symbol for the two-dimensional unit sphere (note that the "2"
 * refers to the dimension of the surface, not the space it is embedded in).
 *
 * <p>This class also defines a framework for decomposing the unit sphere into a hierarchy of
 * "cells". Each cell is a quadrilateral bounded by four geodesics. The top level of the hierarchy
 * is obtained by projecting the six faces of a cube onto the unit sphere, and lower levels are
 * obtained by subdividing each cell into four children recursively.
 *
 * <p>This file also contains documentation of the various coordinate systems and conventions used.
 *
 * @author danieldanciu@google.com (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@SuppressWarnings({"Assertion", "NonFinalStaticField"})
public final class S2 {
  // Declare some frequently used constants
  public static final double M_PI = PI;
  public static final double M_1_PI = 1.0 / PI;
  public static final double M_PI_2 = PI / 2.0;
  public static final double M_PI_4 = PI / 4.0;

  /** Inverse of the root of 2. */
  public static final double M_SQRT1_2 = 1 / sqrt(2);

  public static final double M_SQRT2 = sqrt(2);
  public static final double M_SQRT3 = sqrt(3);
  public static final double M_E = E;

  /**
   * The smallest floating-point value {@code x} such that {@code (1 + x != 1)}, with a value of
   * 2.220446049250313E-16. The C++ value is the same, as you'd expect.
   */
  public static final double DBL_EPSILON = 2.220446049250313E-16;

  /**
   * DBL_ERROR is the maximum rounding error for arithmetic operations. Suppose the "true" result of
   * a calculation is some real number R. Suppose that there is no precise floating point
   * representation. Then there are two consecutive floating point number f1 and f2 such that
   * {@code f1 < R < f2}. Since {@code |f2 - f1| < DBL_EPSILON}, we see that the shortest distance
   * from R to either f1 or f2 must be less than or equal to {@code DBL_EPSILON / 2}. As such,
   * {@code DBL_ERROR = DBL_EPSILON / 2}.
   */
  public static final double DBL_ERROR = DBL_EPSILON / 2;

  /**
   * ROBUST_CROSS_PROD_ERROR is an upper bound on the angle between the vector returned by
   * robustCrossProd(a, b) and the true cross product of "a" and "b".
   *
   * <p>Note that cases where "a" and "b" are exactly proportional but not equal (e.g. a = -b or a =
   * (1 + epsilon) * b) are handled using symbolic perturbations in order to ensure that the result
   * is non-zero and consistent with S2.sign().
   */
  public static final S1Angle ROBUST_CROSS_PROD_ERROR = S1Angle.radians(8 * DBL_ERROR);

  public static final S1Angle EXACT_CROSS_PROD_ERROR = S1Angle.radians(DBL_ERROR);

  /**
   * MIN_NORM is the lower bound on the absolute error of the norm of cross product. If we compute
   * a cross product with a norm below MIN_NORM, then it's possible we have flipped signs and need
   * to fall back to more methods with more precision like calculating with Reals or BigDecimal.
   */
  public static final double MIN_NORM =
      (32 * M_SQRT3 * DBL_ERROR)
          / (ROBUST_CROSS_PROD_ERROR.radians() / DBL_ERROR - (1 + 2 * M_SQRT3));

  // This point is about 66km from the north pole towards the East Siberian Sea, at lat/lng
  // approximately (89.4081206, 165.4655449). See the unit test for more details. It is written here
  // using constant components to ensure the value is identical to other implementations of S2.
  // Warning: Not the same as S2Point.ORIGIN.
  @SuppressWarnings("FloatingPointLiteralPrecision") // Deliberate, same as other implementations.
  private static final S2Point ORIGIN =
      new S2Point(-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195);

  // Together these flags define a cell orientation. If SWAP_MASK is true, then canonical traversal
  // order is flipped around the diagonal (i.e. i and j are swapped with each other). If
  // INVERT_MASK is true, then the traversal order is rotated by 180 degrees (i.e. the bits of i and
  // j are inverted, or equivalently, the axis directions are reversed).
  public static final int SWAP_MASK = 0x01;
  public static final int INVERT_MASK = 0x02;

  /** Mapping Hilbert traversal order to orientation adjustment mask. */
  private static final int[] posToOrientation = {SWAP_MASK, 0, 0, INVERT_MASK + SWAP_MASK};

  /**
   * Returns an XOR bit mask indicating how the orientation of a child subcell is related to the
   * orientation of its parent cell. The returned value can be XOR'd with the parent cell's
   * orientation to give the orientation of the child cell.
   *
   * @param position the position of the subcell in the Hilbert traversal, in the range [0,3].
   * @return a bit mask containing some combination of {@link #SWAP_MASK} and {@link #INVERT_MASK}.
   * @throws IllegalArgumentException if position is out of bounds.
   */
  public static int posToOrientation(int position) {
    Preconditions.checkArgument(0 <= position && position < 4);
    return posToOrientation[position];
  }

  /** Mapping from cell orientation + Hilbert traversal to IJ-index. */
  private static final int[][] posToIj = {
    // 0 1 2 3
    {0, 1, 3, 2}, // canonical order: (0,0), (0,1), (1,1), (1,0)
    {0, 2, 3, 1}, // axes swapped: (0,0), (1,0), (1,1), (0,1)
    {3, 2, 0, 1}, // bits inverted: (1,1), (1,0), (0,0), (0,1)
    {3, 1, 0, 2}, // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
  };

  /**
   * Return the IJ-index of the subcell at the given position in the Hilbert curve traversal with
   * the given orientation. This is the inverse of {@link #ijToPos}.
   *
   * @param orientation the subcell orientation, in the range [0,3].
   * @param position the position of the subcell in the Hilbert traversal, in the range [0,3].
   * @return the IJ-index where {@code 0->(0,0), 1->(0,1), 2->(1,0), 3->(1,1)}.
   * @throws IllegalArgumentException if either parameter is out of bounds.
   */
  public static int posToIJ(int orientation, int position) {
    return posToIj[orientation][position];
  }

  /** Mapping from Hilbert traversal order + cell orientation to IJ-index. */
  private static final int[][] IJ_TO_POS = {
    // (0,0) (0,1) (1,0) (1,1)
    {0, 1, 3, 2}, // canonical order
    {0, 3, 1, 2}, // axes swapped
    {2, 3, 1, 0}, // bits inverted
    {2, 1, 3, 0}, // swapped & inverted
  };

  /**
   * Returns the order in which a specified subcell is visited by the Hilbert curve. This is the
   * inverse of {@link #posToIJ}.
   *
   * @param orientation the subcell orientation, in the range [0,3].
   * @param ijIndex the subcell index where {@code 0->(0,0), 1->(0,1), 2->(1,0), 3->(1,1)}.
   * @return the position of the subcell in the Hilbert traversal, in the range [0,3].
   * @throws IllegalArgumentException if either parameter is out of bounds.
   */
  public static final int ijToPos(int orientation, int ijIndex) {
    return IJ_TO_POS[orientation][ijIndex];
  }

  /**
   * True if some assertions are temporarily disabled to allow deliberate construction of invalid S2
   * objects.
   *
   * <p>Normally, when Java assertions are enabled, as they typically are when running unit tests,
   * creating S2 objects with invalid input parameters will cause AssertionErrors to be thrown. That
   * is usually the desired behavior for unit testing. However, it is necessary to write unit tests
   * that check code behavior when invalid S2 objects are created. Therefore, we need to be able to
   * temporarily disable these assertions.
   *
   * <p>The common case for this is directly in unit tests. For instance, in S2PolygonTest, invalid
   * polygons are created in order to test {@link S2Polygon#findValidationError(S2Error)}.
   *
   * <p>Unit tests that deliberately construct invalid S2 objects should use the public test-only
   * methods: {@link GeometryTestCase#unsafeCreate(Callable)} and {@link
   * GeometryTestCase#unsafeInitialize(Runnable)}, which set skipAssertions for the duration of
   * the provided Runnable or Callable.
   *
   * <p>Assertions within the S2 code that must be disabled to allow invalid S2 objects to be
   * constructed should be written like: {@code assert skipAssertions || isUnitLength(a); }
   */
  static boolean skipAssertions = false;


  /** Defines an area or a length cell metric. Immutable after construction. */
  @Immutable
  public static final class Metric {
    // NOTE: This isn't GWT serializable because writing custom field serializers for inner classes
    // is hard.

    private final double deriv;
    private final int dim;

    /** Defines a cell metric of the given dimension (1 == length, 2 == area). */
    public Metric(int dim, double deriv) {
      this.deriv = deriv;
      this.dim = dim;
    }

    /**
     * The "deriv" value of a metric is a derivative, and must be multiplied by a length or area in
     * (s,t)-space to get a useful value.
     */
    public double deriv() {
      return deriv;
    }

    /** Return the value of a metric for cells at the given level. */
    public double getValue(int level) {
      return scalb(deriv, -dim * level);
    }

    /**
     * Return the level at which the metric has approximately the given value. For example,
     * PROJ.avgEdge.getClosestLevel(0.1) returns the level at which the average cell edge length is
     * approximately 0.1. The return value is always a valid level.
     */
    public int getClosestLevel(double value) {
      return getMinLevel((dim == 1 ? M_SQRT2 : 2) * value);
    }

    /**
     * Return the minimum level such that the metric is at most the given value, or MAX_LEVEL if
     * there is no such level. For example, PROJ.maxDiag.getMinLevel(0.1) returns the minimum level
     * such that all cell diagonal lengths are 0.1 or smaller. The return value is always a valid
     * level.
     */
    public int getMinLevel(double value) {
      if (value <= 0) {
        return S2CellId.MAX_LEVEL;
      }

      // This code is equivalent to computing a floating-point "level" value and rounding up. The
      // getExponent() method returns the exponent corresponding to a fraction in the range [1,2).
      int exponent = Platform.getExponent(value / deriv);
      int level = max(0, min(S2CellId.MAX_LEVEL, -(exponent >> (dim - 1))));
      // assert (level == S2CellId.MAX_LEVEL || getValue(level) <= value);
      // assert (level == 0 || getValue(level - 1) > value);
      return level;
    }

    /**
     * Return the maximum level such that the metric is at least the given value, or zero if there
     * is no such level. For example, PROJ.minWidth.getMaxLevel(0.1) returns the maximum level such
     * that all cells have a minimum width of 0.1 or larger. The return value is always a valid
     * level.
     */
    public int getMaxLevel(double value) {
      if (value <= 0) {
        return S2CellId.MAX_LEVEL;
      }

      // This code is equivalent to computing a floating-point "level" value and rounding down.
      int exponent = Platform.getExponent(deriv / value);
      int level = max(0, min(S2CellId.MAX_LEVEL, exponent >> (dim - 1)));
      // assert (level == 0 || getValue(level) >= value);
      // assert (level == S2CellId.MAX_LEVEL || getValue(level + 1) < value);
      return level;
    }
  }

  /**
   * Return a unique "origin" on the sphere for operations that need a fixed reference point. It
   * should *not* be a point that is commonly used in edge tests in order to avoid triggering code
   * to handle degenerate cases. (This rules out the north and south poles.)
   *
   * <p>Warning: This is completely unrelated to S2Point.ORIGIN.
   */
  public static S2Point origin() {
    return ORIGIN;
  }

  /**
   * Return true if the given point is approximately unit length (this is mainly useful for
   * assertions).
   */
  public static boolean isUnitLength(S2Point p) {
    // normalize() is guaranteed to return a vector whose L2-norm differs from 1 by less than 2 *
    // DBL_EPSILON. Thus the squared L2-norm differs by less than 4 * DBL_EPSILON. The actual
    // calculated norm2() can have up to 1.5 * DBL_EPSILON of additional error. The total error of
    // 5.5 * DBL_EPSILON can then be rounded down since the result must be a representable
    // double-precision value.
    return abs(p.norm2() - 1) <= 5 * DBL_EPSILON; // About 1.11e-15
  }

  private static final S2Point[] ORTHO_BASES = {
    new S2Point(1, 0.0053, 0.00457), new S2Point(0.012, 1, 0.00457), new S2Point(0.012, 0.0053, 1)
  };

  /**
   * Returns a unit-length vector that is orthogonal to {@code a}. Satisfies {@code ortho(-a) =
   * -ortho(a)} for all {@code a}.
   */
  public static S2Point ortho(S2Point a) {
    int k = a.largestAbsComponent() - 1;
    if (k < 0) {
      k = 2;
    }
    return a.crossProd(ORTHO_BASES[k]).normalize();
  }

  /**
   * Returns a unit-length vector used as the reference direction for deciding whether a polygon
   * with semi-open boundaries contains the given vertex "a" (see S2ContainsVertexQuery). The result
   * is unit length and is guaranteed to be different from the given point "a".
   */
  public static S2Point refDir(S2Point a) {
    return ortho(a);
  }

  /**
   * Returns the area of triangle ABC. This method combines two different algorithms to get accurate
   * results for both large and small triangles. The maximum error is about 5e-15 (about 0.25 square
   * meters on the Earth's surface), the same as girardArea() below, but unlike that method it is
   * also accurate for small triangles. Example: when the true area is 100 square meters, area()
   * yields an error about 1 trillion times smaller than girardArea().
   *
   * <p>All points should be unit length, and no two points should be antipodal. The area is always
   * positive.
   */
  public static double area(S2Point a, S2Point b, S2Point c) {
    // assert isUnitLength(a) && isUnitLength(b) && isUnitLength(c);

    // This method is based on l'Huilier's theorem,
    //
    // tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
    //
    // where E is the spherical excess of the triangle (i.e. its area), a, b, c, are the side
    // lengths, and s is the semiperimeter (a + b + c) / 2 .
    //
    // The only significant source of error using l'Huilier's method is the cancellation error of
    // the terms (s-a), (s-b), (s-c). This leads to a *relative* error of about
    // 1e-16 * s / min(s-a, s-b, s-c). This compares to a relative error of about 1e-15 / E using
    // Girard's formula, where E is the true area of the triangle. Girard's formula can be even
    // worse than this for very small triangles, e.g. a triangle with a true area of 1e-30 might
    // evaluate to 1e-5.
    //
    // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
    // dmin = min(s-a, s-b, s-c). This basically includes all triangles except for extremely long
    // and skinny ones.
    //
    // Since we don't know E, we would like a conservative upper bound on the triangle area in terms
    // of s and dmin. It's possible to show that E <= k1 * s * sqrt(s * dmin), where
    // k1 = 2*sqrt(3)/Pi (about 1). Using this, it's easy to show that we should always use
    // l'Huilier's method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore, if
    // dmin < k2 * s^5, the triangle area is at most k3 * s^4, where k3 is about 0.1. Since the
    // best case error using Girard's formula is about 1e-15, this means that we shouldn't even
    // consider it unless s >= 3e-4 or so.
    //
    // The 'standard' formula for finding the angle between two vectors used by S2Point::angle is
    // known to have poor numerical properties when the vectors are nearly (anti-)parallel. Instead
    // we use an alternate formulation by Kahan which has much better properties in that case. See
    // {@link stableAngle()} for more information.
    final double sa = stableAngle(b, c);
    final double sb = stableAngle(c, a);
    final double sc = stableAngle(a, b);

    final double s = 0.5 * (sa + sb + sc);
    if (s >= 3e-4) {
      // Consider whether Girard's formula might be more accurate.
      double s2 = s * s;
      double dmin = s - max(sa, max(sb, sc));
      if (dmin < 1e-2 * s * s2 * s2) {
        // This triangle is skinny enough to consider using Girard's formula. We increase the area
        // by the approximate maximum error in the Girard calculation in order to ensure that this
        // test is conservative.
        double area = girardArea(a, b, c);
        if (dmin < s * (0.1 * (area + 5e-15))) {
          return area;
        }
      }
    }
    // Use l'Huilier's formula.
    return 4 * atan(sqrt(max(
        0.0,
        tan(0.5 * s)
            * tan(0.5 * (s - sa))
            * tan(0.5 * (s - sb))
            * tan(0.5 * (s - sc)))));
  }

  /**
   * Returns the area of the triangle computed using Girard's formula. All points should be unit
   * length, and no two points should be antipodal.
   *
   * <p>This method is about twice as fast as area() but has poor relative accuracy for small
   * triangles. The maximum error is about 5e-15 (about 0.25 square meters on the Earth's surface)
   * and the average error is about 1e-15. These bounds apply to triangles of any size, even as the
   * maximum edge length of the triangle approaches 180 degrees. But note that for such triangles,
   * tiny perturbations of the input points can change the true mathematical area dramatically.
   */
  public static double girardArea(S2Point a, S2Point b, S2Point c) {
    // This is equivalent to the usual Girard's formula but is slightly more accurate, faster to
    // compute, and handles a == b == c without a special case. RobustCrossProd() is necessary to
    // get good accuracy when two of the input points are very close together.
    S2Point ab = S2RobustCrossProd.robustCrossProd(a, b);
    S2Point bc = S2RobustCrossProd.robustCrossProd(b, c);
    S2Point ac = S2RobustCrossProd.robustCrossProd(a, c);
    return max(0.0, ab.angle(ac) - ab.angle(bc) + bc.angle(ac));
  }

  /**
   * Like area(), but returns a positive value for counterclockwise triangles and a negative value
   * otherwise.
   */
  public static double signedArea(S2Point a, S2Point b, S2Point c) {
    return S2Predicates.sign(a, b, c) * area(a, b, c);
  }

  // About centroids:
  // ----------------
  //
  // There are several notions of the "centroid" of a triangle. First, there is the planar
  // centroid, which is simply the centroid of the ordinary (non-spherical) triangle defined by the
  // three vertices. Second, there is the surface centroid, which is defined as the intersection of
  // the three medians of the spherical triangle. It is possible to show that this point is simply
  // the planar centroid projected to the surface of the sphere. Finally, there is the true
  // centroid (mass centroid), which is defined as the surface integral over the spherical triangle
  // of (x,y,z) divided by the triangle area. This is the point that the triangle would rotate
  // around if it was spinning in empty space.
  //
  // The best centroid for most purposes is the true centroid. Unlike the planar and surface
  // centroids, the true centroid behaves linearly as regions are added or subtracted. That is, if
  // you split a triangle into pieces and compute the average of their centroids (weighted by
  // triangle area), the result equals the centroid of the original triangle. This is not true of
  // the other centroids.
  //
  // Also note that the surface centroid may be nowhere near the intuitive "center" of a spherical
  // triangle. For example, consider the triangle with vertices A=(1,eps,0), B=(0,0,1),
  // C=(-1,eps,0) (a quarter-sphere). The surface centroid of this triangle is at S=(0, 2*eps, 1),
  // which is within a distance of 2*eps of the vertex B. Note that the median from A (the segment
  // connecting A to the midpoint of BC) passes through S, since this is the shortest path
  // connecting the two endpoints. On the other hand, the true centroid is at M=(0, 0.5, 0.5),
  // which when projected onto the surface is a much more reasonable interpretation of the "center"
  // of this triangle.

  /**
   * Return the centroid of the planar triangle ABC. This can be normalized to unit length to obtain
   * the "surface centroid" of the corresponding spherical triangle, i.e. the intersection of the
   * three medians. However, note that for large spherical triangles the surface centroid may be
   * nowhere near the intuitive "center" (see example above).
   */
  public static S2Point planarCentroid(S2Point a, S2Point b, S2Point c) {
    return new S2Point((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0);
  }

  /**
   * Returns the true centroid of the spherical geodesic edge AB multiplied by the length of the
   * edge AB. As with triangles, the true centroid of a collection of edges may be computed simply
   * by summing the result of this method for each edge.
   *
   * <p>Note that the planar centroid of a geodesic edge is simply 0.5 * (a + b), while the surface
   * centroid is (a + b).normalize(). However neither of these values is appropriate for computing
   * the centroid of a collection of edges (such as a polyline).
   *
   * <p>Also note that the result of this function is defined to be {@link S2Point#ORIGIN} if the
   * edge is degenerate (and that this is intended behavior).
   */
  public static S2Point trueCentroid(S2Point a, S2Point b) {
    // The centroid (multiplied by length) is a vector toward the midpoint of the edge, whose length
    // is twice the sine of half the angle between the two vertices.
    // Defining theta to be this angle, we have:
    S2Point vDiff = a.sub(b); // Length == 2*sin(theta)
    S2Point vSum = a.add(b); // Length == 2*cos(theta)
    double sin2 = vDiff.norm2();
    double cos2 = vSum.norm2();
    if (cos2 == 0) {
      return S2Point.ORIGIN; // Ignore antipodal edges.
    }
    return vSum.mul(sqrt(sin2 / cos2)); // Length == 2*sin(theta)
  }

  /**
   * Returns the true centroid of the spherical triangle ABC multiplied by the signed area of
   * spherical triangle ABC. The reasons for multiplying by the signed area are (1) this is the
   * quantity that needs to be summed to compute the centroid of a union or difference of triangles,
   * and (2) it's actually easier to calculate this way.
   */
  public static S2Point trueCentroid(S2Point a, S2Point b, S2Point c) {
    // I couldn't find any references for computing the true centroid of a spherical triangle... I
    // have a truly marvellous demonstration of this formula which this margin is too narrow to
    // contain :)

    // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));

    // Use angle() in order to get accurate results for small triangles.
    double aAngle = b.angle(c);
    double bAngle = c.angle(a);
    double cAngle = a.angle(b);
    double ra = (aAngle == 0) ? 1 : (aAngle / sin(aAngle));
    double rb = (bAngle == 0) ? 1 : (bAngle / sin(bAngle));
    double rc = (cAngle == 0) ? 1 : (cAngle / sin(cAngle));

    // Now compute a point M such that:
    //
    //  [Ax Ay Az] [Mx]                       [ra]
    //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
    //  [Cx Cy Cz] [Mz]                       [rc]
    //
    // To improve the numerical stability we subtract the first row (A) from the other two rows;
    // this reduces the cancellation error when A, B, and C are very close together. Then we solve
    // it using Cramer's rule.
    //
    // TODO(user): This code still isn't as numerically stable as it could be. The biggest
    // potential improvement is to compute B-A and C-A more accurately so that (B-A)x(C-A) is always
    // inside triangle ABC.
    S2Point x = new S2Point(a.x, b.x - a.x, c.x - a.x);
    S2Point y = new S2Point(a.y, b.y - a.y, c.y - a.y);
    S2Point z = new S2Point(a.z, b.z - a.z, c.z - a.z);
    S2Point r = new S2Point(ra, rb - ra, rc - ra);
    return new S2Point(
        0.5 * S2Point.scalarTripleProduct(r, y, z),
        0.5 * S2Point.scalarTripleProduct(r, z, x),
        0.5 * S2Point.scalarTripleProduct(r, x, y));
  }

  /**
   * Returns +1 if the edge AB is CCW around the origin, -1 if its clockwise, and 0 if the result is
   * indeterminate.
   */
  public static int planarCCW(R2Vector a, R2Vector b) {
    double sab = (a.dotProd(b) > 0) ? -1 : 1;
    R2Vector vab = R2Vector.add(a, R2Vector.mul(b, sab));
    double da = a.norm2();
    double db = b.norm2();
    double sign;
    if (da < db || (da == db && a.lessThan(b))) {
      sign = a.crossProd(vab) * sab;
    } else {
      sign = vab.crossProd(b);
    }
    if (sign > 0) {
      return 1;
    }
    if (sign < 0) {
      return -1;
    }
    return 0;
  }

  public static int planarOrderedCCW(R2Vector a, R2Vector b, R2Vector c) {
    int sum = 0;
    sum += planarCCW(a, b);
    sum += planarCCW(b, c);
    sum += planarCCW(c, a);
    if (sum > 0) {
      return 1;
    }
    if (sum < 0) {
      return -1;
    }
    return 0;
  }

  /**
   * Return the angle at the vertex B in the triangle ABC. The return value is always in the range
   * [0, Pi]. The points do not need to be normalized. Ensures that angle(a,b,c) == angle(c,b,a) for
   * all a,b,c.
   *
   * <p>The angle is undefined if A or C is diametrically opposite from B, and becomes numerically
   * unstable as the length of edge AB or BC approaches 180 degrees.
   */
  public static double angle(S2Point a, S2Point b, S2Point c) {
    // robustCrossProd() is necessary to get good accuracy when two of the input points are very
    // close together.
    return S2RobustCrossProd.robustCrossProd(a, b).angle(S2RobustCrossProd.robustCrossProd(c, b));
  }

  /**
   * Returns the exterior angle at the vertex B in the triangle ABC. The return value is positive if
   * ABC is counterclockwise and negative otherwise. If you imagine an ant walking from A to B to C,
   * this is the angle that the ant turns at vertex B (positive = left, negative = right).
   *
   * <p>Ensures that turnAngle(a,b,c) == -turnAngle(c,b,a) for all distinct a,b,c. The result is
   * undefined if (a == b || b == c), but is either -Pi or Pi if (a == c). All points should be
   * normalized.
   */
  public static double turnAngle(S2Point a, S2Point b, S2Point c) {
    // We use robustCrossProd() to get good accuracy when two points are very close together, and
    // S2Predicates.sign() to ensure that the sign is correct for turns that are close to 180
    // degrees.
    //
    // Unfortunately we can't save robustCrossProd(a, b) and pass it as the optional 4th argument to
    // S2Predicates.sign(), because it requires a.crossProd(b) exactly (the robust version differs
    // in magnitude).
    double angle =
        S2RobustCrossProd.robustCrossProd(a, b).angle(S2RobustCrossProd.robustCrossProd(b, c));

    // Don't return S2Predicates.sign() * angle because it is legal to have (a == c).
    return (S2Predicates.sign(a, b, c) > 0) ? angle : -angle;
  }

  /**
   * The {@link S2Point#angle(S2Point)} method uses atan2(|AxB|, A.B) to compute the angle between A
   * and B, which can lose about half its precision when A and B are nearly (anti-)parallel.
   *
   * <p>Kahan provides a much more stable form: 2 * atan2(| A*|B| - |A|*B |, | A*|B| + |A|*B |)
   *
   * <p>Since S2Points are unit magnitude by construction we can simplify further:
   * 2*atan2(|A-B|,|A+B|)
   *
   * <p>This likely can't replace S2Point.angle since it requires four magnitude calculations, each
   * of which takes 5 operations + a square root, plus 6 operations to find the sum and difference
   * of the vectors, for a total of 26 + 4 square roots. S2Point.angle requires 19 + 1 square root.
   *
   * <p>Since we always have unit vectors, we can elide two of those magnitude calculations for a
   * total of 16 + 2 square roots which is likely competitive with S2Point.angle performance.
   *
   * <p>Reference: Kahan, W. (2006, Jan 11). "How Futile are Mindless Assessments of Roundoff in
   * Floating-Point Computation?" (p. 47). https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf
   */
  public static double stableAngle(S2Point a, S2Point b) {
    // assert isUnitLength(a) && isUnitLength(b);
    return 2 * atan2(a.sub(b).norm(), a.add(b).norm());
  }

  /**
   * Returns the maximum error in {@link #turnAngle}. The returned value is proportional to the
   * number of vertices and the machine epsilon.
   */
  public static double getTurningAngleMaxError(int numVertices) {
    // The maximum error can be bounded as follows:
    //   2.24 * DBL_EPSILON    for robustCrossProd(b, a)
    //   2.24 * DBL_EPSILON    for robustCrossProd(c, b)
    //   3.25 * DBL_EPSILON    for angle()
    //   2.00 * DBL_EPSILON    for each addition in the Kahan summation
    //   ------------------
    //   9.73 * DBL_EPSILON
    // TODO(user): This error estimate is approximate. There are two issues: (1) SignedArea
    // needs some improvements to ensure that its error is actually never higher than GirardArea,
    // and (2) although the number of triangles in the sum is typically N-2, in theory it could be
    // as high as 2*N for pathological inputs. But in other respects this error bound is very
    // conservative since it assumes that the maximum error is achieved on every triangle.
    final double maxErrorPerVertex = 9.73 * DBL_EPSILON;
    return maxErrorPerVertex * numVertices;
  }

  /**
   * Returns a right-handed coordinate frame (three orthonormal vectors) based on a single point,
   * which will become the third axis.
   */
  public static Matrix getFrame(S2Point p0) {
    S2Point p1 = ortho(p0);
    S2Point p2 = p1.crossProd(p0).normalize();
    return Matrix.fromCols(p2, p1, p0);
  }

  /** Returns a normalized copy {@code p} after rotating it by the rotation matrix {@code r}. */
  static S2Point rotate(S2Point p, Matrix r) {
    Matrix rotated = r.mult(new Matrix(1, p.x, p.y, p.z));
    return new S2Point(rotated.get(0, 0), rotated.get(1, 0), rotated.get(2, 0)).normalize();
  }

  /**
   * Given an orthonormal basis "frame" of column vectors and a point "p", returns the coordinates
   * of "p" with respect to the basis "frame". The resulting point "q" satisfies the identity {@code
   * frame * q == p}.
   */
  static S2Point toFrame(Matrix frame, S2Point p) {
    // The inverse of an orthonormal matrix is its transpose.
    return frame.transpose().mult(Matrix.fromCols(p)).getCol(0);
  }

  /**
   * Given an orthonormal basis "frame" of column vectors and a point "q" with respect to that
   * basis, returns the equivalent point "p" with respect to the standard axis-aligned basis. The
   * result satisfies {@code p == frame * q}.
   */
  static S2Point fromFrame(Matrix frame, S2Point q) {
    return frame.mult(Matrix.fromCols(q)).getCol(0);
  }
  /**
   * Return true if two points are within the given distance in radians of each other. This is
   * mainly useful for testing.
   */
  public static boolean approxEquals(S2Point a, S2Point b, double maxErrorRadians) {
    return a.angle(b) <= maxErrorRadians;
  }

  /**
   * Return true if two points are within an angle of 1e-15 radians of each other. This is about 6.3
   * nanometers on the Earth's surface. This is mainly useful for testing.
   *
   * <p>Note: 1e-15 is a historical relic, related to floating point error that might be introduced
   * when computing intersection points between edges. However, it is insufficient for that purpose
   * and clients should not rely on it. See {@link S2EdgeUtil#DEFAULT_INTERSECTION_TOLERANCE}.
   */
  public static boolean approxEquals(S2Point a, S2Point b) {
    return approxEquals(a, b, 1e-15);
  }

  /**
   * Return true if the given values 'a' and 'b' are within the given 'maxError' of each other. This
   * is mainly useful for testing.
   */
  public static boolean approxEquals(double a, double b, double maxError) {
    return abs(a - b) <= maxError;
  }

  /**
   * Return true if the given values 'a' and 'b' are within 1e-15 of each other. This is mainly
   * useful for testing.
   */
  public static boolean approxEquals(double a, double b) {
    return approxEquals(a, b, 1e-15);
  }

  // Don't instantiate
  private S2() {}
}
