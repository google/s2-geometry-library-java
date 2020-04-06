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

import com.google.common.annotations.GwtCompatible;
import com.google.common.base.Preconditions;

@GwtCompatible
public final strictfp class S2 {
  // Declare some frequently used constants
  public static final double M_PI = Math.PI;
  public static final double M_1_PI = 1.0 / Math.PI;
  public static final double M_PI_2 = Math.PI / 2.0;
  public static final double M_PI_4 = Math.PI / 4.0;
  /** Inverse of the root of 2. */
  public static final double M_SQRT1_2 = 1 / Math.sqrt(2);

  public static final double M_SQRT2 = Math.sqrt(2);
  public static final double M_E = Math.E;

  /** The smallest floating-point value {@code x} such that {@code (1 + x != 1)}. */
  public static final double DBL_EPSILON;

  static {
    double machEps = 1.0d;
    do {
      machEps /= 2.0f;
    } while ((1.0 + (machEps / 2.0)) != 1.0);
    DBL_EPSILON = machEps;
  }

  // This point is about 66km from the north pole towards the East Siberian Sea.  See the unit test
  // for more details. It is written here using constant components to avoid computational errors
  // from producting a different value than other implementations of S2.
  private static final S2Point ORIGIN =
      new S2Point(-0.0099994664350250197, 0.0025924542609324121, 0.99994664350250195);

  // Together these flags define a cell orientation. If SWAP_MASK
  // is true, then canonical traversal order is flipped around the
  // diagonal (i.e. i and j are swapped with each other). If
  // INVERT_MASK is true, then the traversal order is rotated by 180
  // degrees (i.e. the bits of i and j are inverted, or equivalently,
  // the axis directions are reversed).
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

  /** Defines an area or a length cell metric. */
  @GwtCompatible(emulated = true, serializable = false)
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
      return Math.scalb(deriv, -dim * level);
    }

    /**
     * Return the level at which the metric has approximately the given value. For example,
     * S2::kAvgEdge.GetClosestLevel(0.1) returns the level at which the average cell edge length is
     * approximately 0.1. The return value is always a valid level.
     */
    public int getClosestLevel(double value) {
      return getMinLevel((dim == 1 ? S2.M_SQRT2 : 2) * value);
    }

    /**
     * Return the minimum level such that the metric is at most the given value, or
     * S2CellId::kMaxLevel if there is no such level. For example, S2::kMaxDiag.GetMinLevel(0.1)
     * returns the minimum level such that all cell diagonal lengths are 0.1 or smaller. The return
     * value is always a valid level.
     */
    public int getMinLevel(double value) {
      if (value <= 0) {
        return S2CellId.MAX_LEVEL;
      }

      // This code is equivalent to computing a floating-point "level"
      // value and rounding up.  The getExponent() method returns the
      // exponent corresponding to a fraction in the range [1,2).
      int exponent = Platform.getExponent(value / deriv);
      int level = Math.max(0, Math.min(S2CellId.MAX_LEVEL, -(exponent >> (dim - 1))));
      // assert (level == S2CellId.MAX_LEVEL || getValue(level) <= value);
      // assert (level == 0 || getValue(level - 1) > value);
      return level;
    }

    /**
     * Return the maximum level such that the metric is at least the given value, or zero if there
     * is no such level. For example, S2.kMinWidth.GetMaxLevel(0.1) returns the maximum level such
     * that all cells have a minimum width of 0.1 or larger. The return value is always a valid
     * level.
     */
    public int getMaxLevel(double value) {
      if (value <= 0) {
        return S2CellId.MAX_LEVEL;
      }

      // This code is equivalent to computing a floating-point "level"
      // value and rounding down.
      int exponent = Platform.getExponent(deriv / value);
      int level = Math.max(0, Math.min(S2CellId.MAX_LEVEL, exponent >> (dim - 1)));
      // assert (level == 0 || getValue(level) >= value);
      // assert (level == S2CellId.MAX_LEVEL || getValue(level + 1) < value);
      return level;
    }
  }

  /**
   * Return a unique "origin" on the sphere for operations that need a fixed reference point. It
   * should *not* be a point that is commonly used in edge tests in order to avoid triggering code
   * to handle degenerate cases. (This rules out the north and south poles.)
   */
  public static S2Point origin() {
    return ORIGIN;
  }

  /**
   * Return true if the given point is approximately unit length (this is mainly useful for
   * assertions).
   */
  public static boolean isUnitLength(S2Point p) {
    // Normalize() is guaranteed to return a vector whose L2-norm differs from 1
    // by less than 2 * DBL_EPSILON.  Thus the squared L2-norm differs by less
    // than 4 * DBL_EPSILON.  The actual calculated Norm2() can have up to 1.5 *
    // DBL_EPSILON of additional error.  The total error of 5.5 * DBL_EPSILON
    // can then be rounded down since the result must be a representable
    // double-precision value.
    return Math.abs(p.norm2() - 1) <= 5 * DBL_EPSILON; // About 1.11e-15
  }

  /**
   * Return true if edge AB crosses CD at a point that is interior to both edges. Properties:
   *
   * <p>(1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d) (2) SimpleCrossing(c,d,a,b) ==
   * SimpleCrossing(a,b,c,d)
   */
  public static boolean simpleCrossing(S2Point a, S2Point b, S2Point c, S2Point d) {
    // We compute SimpleCCW() for triangles ACB, CBD, BDA, and DAC. All
    // of these triangles need to have the same orientation (CW or CCW)
    // for an intersection to exist. Note that this is slightly more
    // restrictive than the corresponding definition for planar edges,
    // since we need to exclude pairs of line segments that would
    // otherwise "intersect" by crossing two antipodal points.

    S2Point ab = S2Point.crossProd(a, b);
    S2Point cd = S2Point.crossProd(c, d);
    double acb = -ab.dotProd(c);
    double cbd = -cd.dotProd(b);
    double bda = ab.dotProd(d);
    double dac = cd.dotProd(a);

    return (acb * cbd > 0) && (cbd * bda > 0) && (bda * dac > 0);
  }

  /**
   * Return a vector "c" that is orthogonal to the given unit-length vectors "a" and "b". This
   * function is similar to a.CrossProd(b) except that it does a better job of ensuring
   * orthogonality when "a" is nearly parallel to "b", and it returns a non-zero result even when a
   * == b or a == -b.
   *
   * <p>It satisfies the following properties (RCP == robustCrossProd):
   *
   * <p>(1) RCP(a,b) != 0 for all a, b (2) RCP(b,a) == -RCP(a,b) unless a == b or a == -b (3)
   * RCP(-a,b) == -RCP(a,b) unless a == b or a == -b (4) RCP(a,-b) == -RCP(a,b) unless a == b or a
   * == -b
   */
  public static S2Point robustCrossProd(S2Point a, S2Point b) {
    // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
    // approaches zero.  This leads to situations where a.CrossProd(b) is not
    // very orthogonal to "a" and/or "b".  We could fix this using Gram-Schmidt,
    // but we also want b.robustCrossProd(a) == -a.robustCrossProd(b).
    //
    // The easiest fix is to just compute the cross product of (b+a) and (b-a).
    // Mathematically, this cross product is exactly twice the cross product of
    // "a" and "b", but it has the numerical advantage that (b+a) and (b-a)
    // are always perpendicular (since "a" and "b" are unit length).  This
    // yields a result that is nearly orthogonal to both "a" and "b" even if
    // these two values differ only in the lowest bit of one component.
    // assert (isUnitLength(a) && isUnitLength(b));
    S2Point x = S2Point.crossProd(S2Point.add(b, a), S2Point.sub(b, a));
    if (!x.equalsPoint(S2Point.ORIGIN)) {
      return x;
    }

    // The only result that makes sense mathematically is to return zero, but
    // we find it more convenient to return an arbitrary orthogonal vector.
    return ortho(a);
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
    return S2Point.normalize(S2Point.crossProd(a, ORTHO_BASES[k]));
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
    // where E is the spherical excess of the triangle (i.e. its area),
    // a, b, c, are the side lengths, and
    // s is the semiperimeter (a + b + c) / 2 .
    //
    // The only significant source of error using l'Huilier's method is the
    // cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
    // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
    // to a relative error of about 1e-15 / E using Girard's formula, where E is
    // the true area of the triangle. Girard's formula can be even worse than
    // this for very small triangles, e.g. a triangle with a true area of 1e-30
    // might evaluate to 1e-5.
    //
    // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
    // dmin = min(s-a, s-b, s-c). This basically includes all triangles
    // except for extremely long and skinny ones.
    //
    // Since we don't know E, we would like a conservative upper bound on
    // the triangle area in terms of s and dmin. It's possible to show that
    // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
    // Using this, it's easy to show that we should always use l'Huilier's
    // method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
    // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
    // k3 is about 0.1. Since the best case error using Girard's formula
    // is about 1e-15, this means that we shouldn't even consider it unless
    // s >= 3e-4 or so.

    // We use volatile doubles to force the compiler to truncate all of these
    // quantities to 64 bits. Otherwise it may compute a value of dmin > 0
    // simply because it chose to spill one of the intermediate values to
    // memory but not one of the others.
    final double sa = b.angle(c);
    final double sb = c.angle(a);
    final double sc = a.angle(b);
    final double s = 0.5 * (sa + sb + sc);
    if (s >= 3e-4) {
      // Consider whether Girard's formula might be more accurate.
      double s2 = s * s;
      double dmin = s - Math.max(sa, Math.max(sb, sc));
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
    return 4
        * Math.atan(
            Math.sqrt(
                Math.max(
                    0.0,
                    Math.tan(0.5 * s)
                        * Math.tan(0.5 * (s - sa))
                        * Math.tan(0.5 * (s - sb))
                        * Math.tan(0.5 * (s - sc)))));
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
    // This is equivalent to the usual Girard's formula but is slightly more
    // accurate, faster to compute, and handles a == b == c without a special
    // case. RobustCrossProd() is necessary to get good accuracy when two of the
    // input points are very close together.
    S2Point ab = S2.robustCrossProd(a, b);
    S2Point bc = S2.robustCrossProd(b, c);
    S2Point ac = S2.robustCrossProd(a, c);
    return Math.max(0.0, ab.angle(ac) - ab.angle(bc) + bc.angle(ac));
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
  // There are several notions of the "centroid" of a triangle. First, there
  // is the planar centroid, which is simply the centroid of the ordinary
  // (non-spherical) triangle defined by the three vertices. Second, there is
  // the surface centroid, which is defined as the intersection of the three
  // medians of the spherical triangle. It is possible to show that this
  // point is simply the planar centroid projected to the surface of the
  // sphere. Finally, there is the true centroid (mass centroid), which is
  // defined as the surface integral over the spherical triangle of (x,y,z)
  // divided by the triangle area. This is the point that the triangle would
  // rotate around if it was spinning in empty space.
  //
  // The best centroid for most purposes is the true centroid. Unlike the
  // planar and surface centroids, the true centroid behaves linearly as
  // regions are added or subtracted. That is, if you split a triangle into
  // pieces and compute the average of their centroids (weighted by triangle
  // area), the result equals the centroid of the original triangle. This is
  // not true of the other centroids.
  //
  // Also note that the surface centroid may be nowhere near the intuitive
  // "center" of a spherical triangle. For example, consider the triangle
  // with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
  // The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
  // within a distance of 2*eps of the vertex B. Note that the median from A
  // (the segment connecting A to the midpoint of BC) passes through S, since
  // this is the shortest path connecting the two endpoints. On the other
  // hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
  // the surface is a much more reasonable interpretation of the "center" of
  // this triangle.

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
    return vSum.mul(Math.sqrt(sin2 / cos2)); // Length == 2*sin(theta)
  }

  /**
   * Returns the true centroid of the spherical triangle ABC multiplied by the signed area of
   * spherical triangle ABC. The reasons for multiplying by the signed area are (1) this is the
   * quantity that needs to be summed to compute the centroid of a union or difference of triangles,
   * and (2) it's actually easier to calculate this way.
   */
  public static S2Point trueCentroid(S2Point a, S2Point b, S2Point c) {
    // I couldn't find any references for computing the true centroid of a
    // spherical triangle... I have a truly marvellous demonstration of this
    // formula which this margin is too narrow to contain :)

    // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));

    // Use angle() in order to get accurate results for small triangles.
    double aAngle = b.angle(c);
    double bAngle = c.angle(a);
    double cAngle = a.angle(b);
    double ra = (aAngle == 0) ? 1 : (aAngle / Math.sin(aAngle));
    double rb = (bAngle == 0) ? 1 : (bAngle / Math.sin(bAngle));
    double rc = (cAngle == 0) ? 1 : (cAngle / Math.sin(cAngle));

    // Now compute a point M such that:
    //
    //  [Ax Ay Az] [Mx]                       [ra]
    //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
    //  [Cx Cy Cz] [Mz]                       [rc]
    //
    // To improve the numerical stability we subtract the first row (A) from the
    // other two rows; this reduces the cancellation error when A, B, and C are
    // very close together.  Then we solve it using Cramer's rule.
    //
    // TODO(user): This code still isn't as numerically stable as it could be.
    // The biggest potential improvement is to compute B-A and C-A more
    // accurately so that (B-A)x(C-A) is always inside triangle ABC.
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
   * [0, Pi]. The points do not need to be normalized. Ensures that Angle(a,b,c) == Angle(c,b,a) for
   * all a,b,c.
   *
   * <p>The angle is undefined if A or C is diametrically opposite from B, and becomes numerically
   * unstable as the length of edge AB or BC approaches 180 degrees.
   */
  public static double angle(S2Point a, S2Point b, S2Point c) {
    // robustCrossProd() is necessary to get good accuracy when two of the input
    // points are very close together.
    return robustCrossProd(a, b).angle(robustCrossProd(c, b));
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
    double angle = robustCrossProd(a, b).angle(robustCrossProd(b, c));

    // Don't return S2Predicates.sign() * angle because it is legal to have (a == c).
    return (S2Predicates.sign(a, b, c) > 0) ? angle : -angle;
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
    final double maxErrorPerVertex = 9.73 * S2.DBL_EPSILON;
    return maxErrorPerVertex * numVertices;
  }

  /**
   * Returns a right-handed coordinate frame (three orthonormal vectors) based on a single point,
   * which will become the third axis.
   */
  public static Matrix3x3 getFrame(S2Point p0) {
    S2Point p1 = ortho(p0);
    S2Point p2 = S2Point.normalize(S2Point.crossProd(p1, p0));
    return Matrix3x3.fromCols(p2, p1, p0);
  }

  /** Returns a normalized copy {@code p} after rotating it by the rotation matrix {@code r}. */
  static S2Point rotate(S2Point p, Matrix3x3 r) {
    Matrix3x3 rotated = r.mult(new Matrix3x3(1, p.x, p.y, p.z));
    return S2Point.normalize(new S2Point(rotated.get(0, 0), rotated.get(1, 0), rotated.get(2, 0)));
  }

  /** Converts 'p' to the basis given in 'frame'. */
  static S2Point toFrame(Matrix3x3 frame, S2Point p) {
    // The inverse of an orthonormal matrix is its transpose.
    return frame.transpose().mult(Matrix3x3.fromCols(p)).getCol(0);
  }

  /** Converts 'p' from the basis given in 'frame'. */
  static S2Point fromFrame(Matrix3x3 frame, S2Point q) {
    return frame.mult(Matrix3x3.fromCols(q)).getCol(0);
  }

  /**
   * Return true if two points are within the given distance of each other (mainly useful for
   * testing).
   */
  public static boolean approxEquals(S2Point a, S2Point b, double maxError) {
    return a.angle(b) <= maxError;
  }

  public static boolean approxEquals(S2Point a, S2Point b) {
    return approxEquals(a, b, 1e-15);
  }

  public static boolean approxEquals(double a, double b, double maxError) {
    return Math.abs(a - b) <= maxError;
  }

  public static boolean approxEquals(double a, double b) {
    return approxEquals(a, b, 1e-15);
  }

  // Don't instantiate
  private S2() {}
}
