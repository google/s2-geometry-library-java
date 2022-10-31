/*
 * Copyright 2018 Google Inc.
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
import static com.google.common.geometry.S2.M_SQRT1_2;
import static com.google.common.geometry.S2.M_SQRT2;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import java.math.BigDecimal;

/**
 * A collection of geometric predicates core to the robustness of the S2 library. In particular,
 *
 * <ul>
 *   <li>{@link #sign sign(A, B, C)}: Compute the orientation of a triple of points as clockwise
 *       (-1), colinear (0), or counter-clockwise (1).
 * </ul>
 */
public final strictfp class S2Predicates {
  /** Maximum rounding error of a 64 bit double. */
  private static final double DBL_ERR = DBL_EPSILON / 2;

  /**
   * Rounding error of numeric type used for computation. This is a distinct value from {@link
   * #DBL_ERR} to avoid recomputing the error bounds when we may later use higher precision (but
   * inexact) types for computations.
   */
  private static final double T_ERR = DBL_ERR;

  /** A predefined S1ChordAngle representing (approximately) 45 degrees. */
  private static final S1ChordAngle DEG_45 = S1ChordAngle.fromLength2(2 - M_SQRT2);

  // Preallocated constants.

  private static final BigDecimal QUARTER = new BigDecimal("0.25");
  private static final BigDecimal HALF = new BigDecimal("0.5");
  private static final BigDecimal TWO = new BigDecimal("2");
  private static final BigDecimal FOUR = new BigDecimal("4");

  /**
   * Returns +1 if the points A, B, C are counterclockwise, -1 if the points are clockwise, and 0 if
   * any two points are the same. This function is essentially like taking the sign of the
   * determinant of ABC, except that it has additional logic to make sure that the above properties
   * hold even when the three points are coplanar, and to deal with the limitations of
   * floating-point arithmetic.
   *
   * <p>Sign satisfies the following conditions:
   *
   * <ol>
   *   <li>{@code sign(a,b,c) == 0} iff {@code a.equals(b) || b.equals(c) || c.equals(a)}
   *   <li>{@code sign(b,c,a) == sign(a,b,c)}, for all a,b,c
   *   <li>{@code sign(c,b,a) == -sign(a,b,c)}, for all a,b,c
   * </ol>
   *
   * <p>In other words:
   *
   * <ol>
   *   <li>The result is zero if and only if two points are the same.
   *   <li>Rotating the order of the arguments does not affect the result.
   *   <li>Exchanging any two arguments inverts the result.
   * </ol>
   *
   * <p>On the other hand, note that it is not true in general that {@code sign(-a,b,c) ==
   * -sign(a,b,c)}, or any similar identities involving antipodal points.
   */
  public static int sign(S2Point a, S2Point b, S2Point c) {
    return Sign.sign(a, b, c, true);
  }

  /** Tests of whether three points represent a left turn (+1), right turn (-1), or neither (0). */
  public static class Sign {
    /** No instantiation. */
    private Sign() {}

    /**
     * This version of Sign returns +1 if the points are definitely CCW, -1 if they are definitely
     * CW, and 0 if two points are identical or the result is uncertain. Uncertain cases can be
     * resolved, if desired, by {@link #expensive}.
     *
     * <p>The purpose of this method is to allow additional cheap tests to be done, where possible,
     * in order to avoid calling more expensive operations unnecessarily.
     */
    public static int triage(S2Point a, S2Point b, S2Point c) {
      // There are 14 multiplications and additions to compute the determinant below. Since all
      // three points are normalized, it is possible to show that the average rounding error per
      // operation does not exceed 2**-54, the maximum rounding error for an operation whose result
      // magnitude is in the range [0.5,1). Therefore, if the absolute value of the determinant is
      // greater than 2*14*(2**-54), the determinant will have the same sign even if the arguments
      // are rotated (which produces a mathematically equivalent result but with potentially
      // different rounding errors).
      final double kMaxDetError = 1.6e-15; // 2 * 14 * 2**-54
      // assert S2.isUnitLength(a);
      // assert S2.isUnitLength(b);
      // assert S2.isUnitLength(c);

      double det = S2Point.scalarTripleProduct(c, a, b);

      // Double-check borderline cases in debug mode.
      // assert (abs(det) <= kMinAbsValue) || (abs(det) >= 100 * kMinAbsValue)
      //    || (det * expensiveCCW(a, b, c) > 0);

      if (det >= kMaxDetError) {
        return 1;
      }
      if (det <= -kMaxDetError) {
        return -1;
      }
      return 0;
    }

    /**
     * Returns the sign of the turn ABC. Exactly straight points will only result in 0 if 'perturb'
     * is false, otherwise the points are perturbed according to the rules in Simulation of
     * Simplicity, to provide a logically consistent non-zero result for all inputs.
     */
    public static int sign(S2Point a, S2Point b, S2Point c, boolean perturb) {
      // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));
      int ccw = Sign.triage(a, b, c);
      if (ccw == 0) {
        ccw = Sign.expensive(a, b, c, perturb);
      }
      return ccw;
    }

    /**
     * Returns the sign of the determinant using more expensive techniques. To be used when the
     * magnitude of the determinant is close enough to zero that its value is uncertain using faster
     * but less robust techniques.
     */
    public static int expensive(S2Point a, S2Point b, S2Point c, boolean perturb) {
      // Return zero if and only if two points are the same, the first property of sign().
      if (a.equalsPoint(b) || b.equalsPoint(c) || c.equalsPoint(a)) {
        return 0;
      }

      int sign = stable(a, b, c);
      if (sign != 0) {
        return sign;
      }

      return exact(a, b, c, perturb);
    }

    /**
     * Compute the determinant in a numerically stable way. Unlike {@link #triage}, this method can
     * usually compute the correct determinant sign even when all three points are as collinear as
     * possible. For example if three points are spaced 1km apart along a random line on the Earth's
     * surface using the nearest representable points, there is only a 0.4% chance that this method
     * will not be able to find the determinant sign. The probability of failure decreases as the
     * points get closer together; if the collinear points are 1 meter apart, the failure rate drops
     * to 0.0004%.
     *
     * <p>This method could be extended to also handle nearly-antipodal points (and in fact an
     * earlier version of this code did exactly that), but antipodal points are rare in practice so
     * it seems better to simply fall back to exact arithmetic in that case.
     */
    public static int stable(S2Point a, S2Point b, S2Point c) {
      S2Point ab = b.sub(a);
      S2Point bc = c.sub(b);
      S2Point ca = a.sub(c);
      double ab2 = ab.norm2();
      double bc2 = bc.norm2();
      double ca2 = ca.norm2();

      // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been cyclically
      // permuted if necessary so that AB is the longest edge. (This minimizes the magnitude of
      // cross product.) At the same time we also compute the maximum error in the determinant.
      // Using a similar technique to the one used for kMaxDetError, the error is at most
      //
      //   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
      //
      // where e = 0.5 * DBL_EPSILON. If the determinant magnitude is larger than this value then
      // we know its sign with certainty.
      final double detErrorMultiplier = 3.2321 * DBL_EPSILON; // see above
      double det;
      double maxError;
      if (ab2 >= bc2 && ab2 >= ca2) {
        // AB is the longest edge, so compute (A-C)x(B-C).C.
        det = -S2Point.scalarTripleProduct(c, ca, bc);
        maxError = detErrorMultiplier * sqrt(ca2 * bc2);
      } else if (bc2 >= ca2) {
        // BC is the longest edge, so compute (B-A)x(C-A).A.
        det = -S2Point.scalarTripleProduct(a, ab, ca);
        maxError = detErrorMultiplier * sqrt(ab2 * ca2);
      } else {
        // CA is the longest edge, so compute (C-B)x(A-B).B.
        det = -S2Point.scalarTripleProduct(b, bc, ab);
        maxError = detErrorMultiplier * sqrt(bc2 * ab2);
      }
      return signum(det, maxError);
    }

    /**
     * Computes the determinant using exact arithmetic and/or symbolic permutations. Requires that
     * the three points are distinct.
     */
    public static int exact(S2Point a, S2Point b, S2Point c, boolean perturb) {
      // assert !a.equalsPoint(b) && !b.equalsPoint(c) && !c.equalsPoint(a);

      // Check the determinant using any available Platform optimizations, which may return zero
      // but are faster than BigDecimal techniques. The current Platform implementation computes the
      // determinant using exact arithmetic with Reals.
      int sign = Platform.sign(a, b, c);
      if (sign != 0) {
        return sign;
      }

      // Sort the 3 points in lexicographic order, keeping track of the sign of the permutation.
      // (Each exchange inverts the sign of the determinant.)
      int permSign = 1;
      if (a.compareTo(b) > 0) {
        S2Point t = a;
        a = b;
        b = t;
        permSign = -permSign;
      }
      if (b.compareTo(c) > 0) {
        S2Point t = b;
        b = c;
        c = t;
        permSign = -permSign;
      }
      if (a.compareTo(b) > 0) {
        S2Point t = a;
        a = b;
        b = t;
        permSign = -permSign;
      }
      // assert a.compareTo(b) < 0 && b.compareTo(c) < 0;

      // The BigDecimal approach used below can't handle NaNs. Prior to use of this type, NaN always
      // resulted in -1, regardless of other coordinate values, so maintain that behavior here.
      if (!a.equalsPoint(a) || !b.equalsPoint(b) || !c.equalsPoint(c)) {
        return -permSign;
      }

      // Check the determinant using BigDecimal, which is exact but very slow. Save the cross
      // product for use in symbolic perturbation if necessary.
      BigPoint xa = big(a);
      BigPoint xb = big(b);
      BigPoint xc = big(c);
      BigPoint xbc = xb.crossProd(xc);
      sign = xbc.dotProd(xa).signum();
      if (sign != 0) {
        return permSign * sign;
      }

      // If !perturb, then the caller has requested no SoS tie breaking, so simply return 0 now.
      if (!perturb) {
        return 0;
      }

      // Resort to symbolic perturbations to resolve a stable non-zero result.
      sign = sos(xa, xb, xc, xbc);
      // assert 0 != sign;
      return permSign * sign;
    }

    /**
     * Returns the sign of the determinant of three column vectors A, B, C under a model where every
     * possible S2Point is slightly perturbed by a unique infinitesimal amount such that no three
     * perturbed points are collinear and no four points are coplanar. The perturbations are so
     * small that they do not change the sign of any determinant that was non-zero before the
     * perturbations, and therefore can be safely ignored unless the determinant of three points is
     * exactly zero (using multiple-precision arithmetic).
     *
     * <p>Since the symbolic perturbation of a given point is fixed (i.e., the perturbation is the
     * same for all calls to this method and does not depend on the other two arguments), the
     * results of this method are always self-consistent. It will never return results that would
     * correspond to an "impossible" configuration of non-degenerate points.
     *
     * <p>Requirements:
     *
     * <ul>
     *   <li>The 3x3 determinant of A, B, C must be exactly zero.
     *   <li>The points must be distinct, with {@code A < B < C} in lexicographic order.
     * </ul>
     *
     * <p>Returns +1 or -1 according to the sign of the determinant after the symbolic perturbations
     * are taken into account. This method is suitable as a tie breaker when all other tests fail to
     * choose a distinct turning direction.
     *
     * <p>Reference:
     *
     * <pre>
     * "Simulation of Simplicity"
     * Edelsbrunner and Muecke, ACM Transactions on Graphics, 1990
     * </pre>
     */
    public static int sos(BigPoint a, BigPoint b, BigPoint c, BigPoint bc) {
      // This method requires that the points are sorted in lexicographically increasing order. This
      // is because every possible S2Point has its own symbolic perturbation such that if A < B then
      // the symbolic perturbation for A is much larger than the perturbation for B.
      //
      // Alternatively, we could sort the points in this method and keep track of the sign of the
      // permutation, but it is more efficient to do this before converting the inputs to the
      // multi-precision representation, and this also lets us re-use the result of the cross
      // product B x C.
      //
      // assert a.compareTo(b) < 0 && b.compareTo(c) < 0;

      // Every input coordinate x[i] is assigned a symbolic perturbation dx[i]. We then compute the
      // sign of the determinant of the perturbed points,
      // i.e.
      //               | a[0]+da[0]  a[1]+da[1]  a[2]+da[2] |
      //               | b[0]+db[0]  b[1]+db[1]  b[2]+db[2] |
      //               | c[0]+dc[0]  c[1]+dc[1]  c[2]+dc[2] |
      //
      // The perturbations are chosen such that
      //
      //   da[2] > da[1] > da[0] > db[2] > db[1] > db[0] > dc[2] > dc[1] > dc[0]
      //
      // where each perturbation is so much smaller than the previous one that we don't even need to
      // consider it unless the coefficients of all previous perturbations are zero. In fact, it is
      // so small that we don't need to consider it unless the coefficient of all products of the
      // previous perturbations are zero. For example, we don't need to consider the coefficient of
      // db[1] unless the coefficient of db[2]*da[0] is zero.
      //
      // The following code simply enumerates the coefficients of the perturbations (and products of
      // perturbations) that appear in the determinant above, in order of decreasing perturbation
      // magnitude. The first non-zero coefficient determines the sign of the result. The easiest
      // way to enumerate the coefficients in the correct order is to pretend that each perturbation
      // is some tiny value "eps" raised to a power of two:
      //
      // eps**    1      2      4      8     16     32     64     128    256
      //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]  dc[2]  dc[1]  dc[0]
      //
      // Essentially we can then just count in binary and test the corresponding subset of
      // perturbations at each step. So for example, we must test the coefficient of db[2]*da[0]
      // before db[1] because eps**12 > eps**16.
      //
      // Of course, not all products of these perturbations appear in the determinant above, since
      // the determinant only contains the products of elements in distinct rows and columns. Thus
      // we don't need to consider da[2]*da[1], db[1]*da[1], etc. Furthermore, sometimes different
      // pairs of perturbations have the same coefficient in the determinant; for example,
      // da[1]*db[0] and db[1]*da[0] have the same coefficient (c[2]). Therefore we only need to
      // test this coefficient the first time we encounter it in the binary order above (which will
      // be db[1]*da[0]).
      //
      // The sequence of tests below also appears in Table 4-ii of the paper referenced above, if
      // you just want to look it up, with the following translations: [a,b,c] -> [i,j,k] and
      // [0,1,2] -> [1,2,3]. Also note that some of the signs are different because the opposite
      // cross product is used (e.g., B x C rather than C x B).

      int sign = bc.z.signum(); // da[2]
      if (sign != 0) {
        return sign;
      }
      sign = bc.y.signum(); // da[1]
      if (sign != 0) {
        return sign;
      }
      sign = bc.x.signum(); // da[0]
      if (sign != 0) {
        return sign;
      }

      sign = c.x.multiply(a.y).subtract(c.y.multiply(a.x)).signum(); // db[2]
      if (sign != 0) {
        return sign;
      }
      sign = c.x.signum(); // db[2] * da[1]
      if (sign != 0) {
        return sign;
      }
      sign = -c.y.signum(); // db[2] * da[0]
      if (sign != 0) {
        return sign;
      }
      sign = c.z.multiply(a.x).subtract(c.x.multiply(a.z)).signum(); // db[1]
      if (sign != 0) {
        return sign;
      }
      sign = c.z.signum(); // db[1] * da[0]
      if (sign != 0) {
        return sign;
      }

      // The following test is listed in the paper, but it is redundant because the previous tests
      // guarantee that C == (0, 0, 0).
      // assert 0 == cy.multiply(az).subtract(cz.multiply(ay)).signum(); // db[0]

      sign = a.x.multiply(b.y).subtract(a.y.multiply(b.x)).signum(); // dc[2]
      if (sign != 0) {
        return sign;
      }
      sign = -b.x.signum(); // dc[2] * da[1]
      if (sign != 0) {
        return sign;
      }
      sign = b.y.signum(); // dc[2] * da[0]
      if (sign != 0) {
        return sign;
      }
      sign = a.x.signum(); // dc[2] * db[1]
      if (sign != 0) {
        return sign;
      }
      return 1; // dc[2] * db[1] * da[0]
    }
  }

  /**
   * Return true if the edges OA, OB, and OC are encountered in that order while sweeping CCW around
   * the point O. You can think of this as testing whether {@code A <= B <= C} with respect to a
   * continuous CCW ordering around O.
   *
   * <p>Properties:
   *
   * <ol>
   *   <li>If {@code orderedCCW(a,b,c,o) && orderedCCW(b,a,c,o)}, then {@code a == b}
   *   <li>If {@code orderedCCW(a,b,c,o) && orderedCCW(a,c,b,o)}, then {@code b == c}
   *   <li>If {@code orderedCCW(a,b,c,o) && orderedCCW(c,b,a,o)}, then {@code a == b == c}
   *   <li>If {@code a == b || b == c}, then {@code orderedCCW(a,b,c,o)} is true
   *   <li>Otherwise if {@code a == c}, then {@code orderedCCW(a,b,c,o)} is false
   * </ol>
   */
  public static boolean orderedCCW(S2Point a, S2Point b, S2Point c, S2Point o) {
    // The last inequality below is ">" rather than ">=" so that we return true if A == B or B == C,
    // and otherwise false if A == C. Recall that sign(x,y,z) == -sign(z,y,x) for all x,y,z.

    int sum = 0;
    if (sign(b, o, a) >= 0) {
      ++sum;
    }
    if (sign(c, o, b) >= 0) {
      ++sum;
    }
    if (sign(a, o, c) > 0) {
      ++sum;
    }
    return sum >= 2;
  }

  /**
   * Returns -1, 0, or +1 according to whether {@code AX < BX, A == B}, or {@code AX > BX}
   * respectively. Distances are measured with respect to the positions of X, A, and B as though
   * they were reprojected to lie exactly on the surface of the unit sphere. Furthermore, this
   * method uses symbolic perturbations to ensure that the result is non-zero whenever
   * {@code A != B}, even when {@code AX == BX} exactly, or even when A and B project to the same
   * point on the sphere. Such results are guaranteed to be self-consistent, i.e. if {@code AB < BC}
   * and {@code BC < AC}, then {@code AB < AC}.
   */
  public static int compareDistances(S2Point x, S2Point a, S2Point b) {
    // We start by comparing distances using dot products (i.e., cosine of the angle), because (1)
    // this is the cheapest technique, and (2) it is valid over the entire range of possible angles.
    // (We can only use the sin^2 technique if both angles are less than 90 degrees or both angles
    // are greater than 90 degrees.)
    int sign = CompareDistances.triageCos(x, a, b);
    if (sign != 0) {
      return sign;
    }

    // Optimization for (a == b) to avoid falling back to exact arithmetic.
    if (a.equalsPoint(b)) {
      return 0;
    }

    // It is much better numerically to compare distances using cos(angle) if the distances are near
    // 90 degrees and sin^2(angle) if the distances are near 0 or 180 degrees. We only need to check
    // one of the two angles when making this decision because the fact that the test above failed
    // means that angles "a" and "b" are very close together.
    double cosAX = a.dotProd(x);
    if (cosAX > M_SQRT1_2) {
      // Angles < 45 degrees.
      sign = CompareDistances.triageSin2(x, a, b);
    } else if (cosAX < -M_SQRT1_2) {
      // Angles > 135 degrees. sin^2(angle) is decreasing in this range.
      sign = -CompareDistances.triageSin2(x, a, b);
    }
    if (sign != 0) {
      return sign;
    }
    // TODO(user): Use Real instead of BigDecimal.
    sign = CompareDistances.exact(x, a, b);
    if (sign != 0) {
      return sign;
    }
    return CompareDistances.sos(a, b);
  }

  /**
   * A set of tests to determine which of two points is closer to a reference point. Generally much
   * faster then computing even one exact distance, since most points are "obviously" ordered w.r.t.
   * the reference point. Returns -1, 0, or +1 according to whether AX < BX, A == B, or AX > BX
   * respectively, and may return 0 if the result is indeterminate.
   */
  static class CompareDistances {
    /** No instantiation. */
    private CompareDistances() {}

    /**
     * Returns a cosine-based test result. It handles all angles, is the fastest implementation, but
     * has a wide margin of uncertainty.
     */
    public static int triageCos(S2Point x, S2Point a, S2Point b) {
      double cosAX = cosDistance(a, x);
      double cosBX = cosDistance(b, x);
      return compare(cosBX, cosDistanceError(cosBX), cosAX, cosDistanceError(cosAX));
    }

    /**
     * Returns the test result using a more accurate sine strategy, which only allows angles either
     * both below -90 or both above +90 degrees.
     */
    public static int triageSin2(S2Point x, S2Point a, S2Point b) {
      double sin2AX = sin2Distance(a, x);
      double sin2BX = sin2Distance(b, x);
      return compare(sin2AX, sin2DistanceError(sin2AX), sin2BX, sin2DistanceError(sin2BX));
    }

    /** Returns a BigDecimal-based test result, which is slow but handles all input. */
    public static int exact(S2Point x, S2Point a, S2Point b) {
      return exact(big(x), big(a), big(b));
    }

    /** Returns a BigDecimal-based test result, which is slow but handles all input. */
    public static int exact(BigPoint x, BigPoint a, BigPoint b) {
      // This code produces the same result as though all points were reprojected to lie exactly on
      // the surface of the unit sphere. It is based on testing whether x.dotProd(a.normalize()) <
      // x.dotProd(b.normalize()), reformulated so that it can be evaluated using exact arithmetic.
      BigDecimal cosAX = x.dotProd(a);
      BigDecimal cosBX = x.dotProd(b);
      // If the values have different signs, handle that case now before squaring them below.
      int aSign = cosAX.signum();
      int bSign = cosBX.signum();
      if (aSign != bSign) {
        // If cos(AX) > cos(BX), then AX < BX.
        return Integer.compare(bSign, aSign);
      }
      int cmpSign =
          cosBX
              .multiply(cosBX)
              .multiply(a.norm2())
              .compareTo(cosAX.multiply(cosAX).multiply(b.norm2()));
      return aSign * cmpSign;
    }

    /**
     * Given that the exact test returned 0, returns a Simulation of Simplicity symbolic
     * perturbation-based test result to select a consistent non-zero result.
     */
    public static int sos(S2Point a, S2Point b) {
      // Our symbolic perturbation strategy is based on the following model. Similar to "simulation
      // of simplicity", we assign a perturbation to every point such that if A < B, then the
      // symbolic perturbation for A is much, much larger than the symbolic perturbation for B. We
      // imagine that rather than projecting every point to lie exactly on the unit sphere, instead
      // each point is positioned on its own tiny pedestal that raises it just off the surface of
      // the unit sphere. This means that the distance AX is actually the true distance AX plus the
      // (symbolic) heights of the pedestals for A and X. The pedestals are infinitesmally thin, so
      // they do not affect distance measurements except at the two endpoints. If several points
      // project to exactly the same point on the unit sphere, we imagine that they are placed on
      // separate pedestals placed close together, where the distance between pedestals is much,
      // much less than the height of any pedestal; as there are a finite number of S2Points, and
      // therefore a finite number of pedestals, this is possible.
      //
      // If A < B, then A is on a higher pedestal than B, and therefore AX > BX.
      return b.compareTo(a);
    }
  }

  /**
   * Returns -1, 0, or +1 according to whether the distance XY is less than, equal to, or greater
   * than the squared chord distance "r2" respectively. Distances are measured with respect the
   * positions of all points as though they are projected to lie exactly on the surface of the unit
   * sphere.
   */
  public static int compareDistance(S2Point x, S2Point y, double r2) {
    int sign = CompareDistance.triage(x, y, r2);
    if (sign != 0) {
      return sign;
    }
    // TODO(user): Use Real instead of BigDecimal.
    return CompareDistance.exact(x, y, r2);
  }

  /**
   * A set of tests to compare the distance XY and a previously computed distance. When doing many
   * distance tests, this saves a lot of work over computing exact distances only to throw away most
   * of them. Each test returns -1 if the distance XY is less than 'r2', +1 if the distance XY is
   * greater than 'r2', and 0 if the distances are exactly equal or the relation is indeterminate.
   */
  static class CompareDistance {
    /** No instantiation. */
    private CompareDistance() {}

    /**
     * Returns a cosine-based test result. This is the fastest test, it handle all angles, but has a
     * wide margin of uncertainty.
     */
    public static int triageCos(S2Point x, S2Point y, double r2) {
      double cosXY = cosDistance(x, y);
      double cosR = 1 - 0.5 * r2;
      return compare(cosR, 2 * T_ERR * cosR, cosXY, cosDistanceError(cosXY));
    }

    /**
     * Returns a sine-based test result, which has very good accuracy for small angles, although it
     * only handles angles below 90 degrees.
     */
    public static int triageSin2(S2Point x, S2Point y, double r2) {
      // assert r2 < 2.0; // Only valid for distance limits < 90 degrees.
      double xySin2 = sin2Distance(x, y);
      double rSin2 = r2 * (1 - 0.25 * r2);
      return compare(xySin2, sin2DistanceError(xySin2), rSin2, 3 * T_ERR * rSin2);
    }

    /**
     * Returns a test result based on first {@link #triageCos} and then {@link #triageSin2}, so it
     * is fast for any input where exact arithmetic isn't needed.
     */
    public static int triage(S2Point x, S2Point y, double r2) {
      // The Sin2 method is much more accurate for small distances, but it is only valid when the
      // actual distance and the distance limit are both less than 90 degrees. So we always start
      // with the Cos method.
      int sign = triageCos(x, y, r2);
      if (sign == 0 && r2 < DEG_45.getLength2()) {
        sign = triageSin2(x, y, r2);
      }
      return sign;
    }

    /** Calls {@link #exact(BigPoint, BigPoint, BigDecimal)} with its more precise types. */
    public static int exact(S2Point x, S2Point y, double r2) {
      return exact(big(x), big(y), big(r2));
    }

    /**
     * Returns a BigDecimal-based test result, which is exact for all inputs but very slow. In
     * particular, this code produces a result as though all points were reprojected to lie exactly
     * on the surface of the unit sphere. It is based on comparing the cosine of the angle XY (when
     * both points are projected to lie exactly on the sphere) to the given threshold.
     */
    public static int exact(BigPoint x, BigPoint y, BigDecimal r2) {
      BigDecimal cosXY = x.dotProd(y);
      BigDecimal cosR = BigDecimal.ONE.subtract(HALF.multiply(r2));
      // If the values have different signs, handle that case now before squaring them below.
      int xySign = cosXY.signum();
      int rSign = cosR.signum();
      if (xySign != rSign) {
        // If cos(XY) > cos(r), then XY < r.
        return Integer.compare(rSign, xySign);
      }
      int cmpSign = square(cosR).multiply(x.norm2().multiply(y.norm2())).compareTo(square(cosXY));
      return xySign * cmpSign;
    }
  }

  /**
   * Returns -1, 0, or +1 according to whether the distance from the point X to the edge AB is less
   * than, equal to, or greater than the squared chord distance "r2" respectively.
   *
   * <p>Distances are measured with respect to the positions of all points as though they were
   * projected to lie exactly on the surface of the unit sphere.
   *
   * <p>Requires that {@code A} and {@code B} do not project to antipodal points (e.g.,
   * {@code A != -B}). This can occur if for example {@code A == S * B}, for some constant
   * {@code S < 0}.
   *
   * <p>Note all of the predicates defined here could be extended to handle edges consisting of
   * antipodal points by implementing additional symbolic perturbation logic (similar to {@link
   * #sign}) in order to rigorously define the direction of such edges.
   */
  public static int compareEdgeDistance(S2Point x, S2Point a, S2Point b, double r2) {
    // Check that the edge does not consist of antipodal points. This catches the most common case.
    // The full test is in CompareEdgeDistance.exact().
    // assert !a.equals(b.neg());

    int sign = CompareEdgeDistance.triage(x, a, b, r2);
    if (sign != 0) {
      return sign;
    }

    // Optimization for the case where the edge is degenerate.
    if (a.equalsPoint(b)) {
      return compareDistance(x, a, r2);
    }
    // TODO(user): Use Real instead of BigDecimal.

    return CompareEdgeDistance.exact(x, a, b, r2);
  }

  /**
   * A test to compare the distance from point X to edge A with a previously computed distance. When
   * doing many edge distance tests, this saves a lot of work over computing exact distances only to
   * throw them away most of them.
   *
   * <p>Does not offer generally correct results for all inputs, so that multiple strategies may be
   * implemented for different classes of input.
   */
  static class CompareEdgeDistance {
    /** No instantiation. */
    private CompareEdgeDistance() {}

    /**
     * Returns -1, 0, or +1 according to whether the distance from the point X to the edge AB is
     * less than, equal to, or greater than "r2" respectively, and may return 0 if the relation is
     * indeterminate.
     *
     * <p>This test uses double arithmetic, which is reasonably precise but allocates a lot.
     */
    // TODO(user): x.angle(a)+a.angle(b) < r is usually false, so test that first?
    public static int triage(S2Point x, S2Point a, S2Point b, double r2) {
      // First we need to decide whether the closest point is an edge endpoint or somewhere in the
      // interior. To determine this we compute a plane perpendicular to (a, b) that passes through
      // X. Letting M be the normal to this plane, the closest point is in the edge interior if and
      // only if a.M < 0 and b.M > 0. Note that we can use "<" rather than "<=" because if a.M or
      // b.M is zero exactly then it doesn't matter which code path we follow (since the distance to
      // an endpoint and the distance to the edge interior are exactly the same in this case).
      S2Point n = ndCross(a, b);
      S2Point m = n.crossProd(x);
      // For better accuracy when the edge (a,b) is very short, we subtract "x" before computing the
      // dot products with M.
      S2Point aDir = a.sub(x);
      S2Point bDir = b.sub(x);
      double aSign = aDir.dotProd(m);
      double bSign = bDir.dotProd(m);
      double n2 = n.norm2();
      double n1 = sqrt(n2);
      double n1Error = ((3.5 + 8 / sqrt(3)) * n1 + 32 * sqrt(3) * DBL_ERR) * T_ERR;
      double aSignError = n1Error * aDir.norm();
      double bSignError = n1Error * bDir.norm();
      if (abs(aSign) < aSignError || abs(bSign) < bSignError) {
        // It is uncertain whether minimum distance is to an edge vertex or to the edge interior.
        // So compute both distances and check whether they yield the same result.
        int vertexSign = triageLineEndpoints(x, a, b, r2);
        int lineSign = triageLineInterior(x, a, b, r2, n, n1, n2);
        return (vertexSign == lineSign) ? lineSign : 0;
      }
      if (aSign >= 0 || bSign <= 0) {
        // The minimum distance is to an edge endpoint.
        return triageLineEndpoints(x, a, b, r2);
      } else {
        // The minimum distance is to the edge interior.
        return triageLineInterior(x, a, b, r2, n, n1, n2);
      }
    }

    /** Returns the min test result from XA and XB, assuming the projection is A or B. */
    static int triageLineEndpoints(S2Point x, S2Point a, S2Point b, double r2) {
      return min(CompareDistance.triage(x, a, r2), CompareDistance.triage(x, b, r2));
    }

    /** Returns the min test result from XA and XB, assuming the projection is between A and B. */
    static int triageLineInterior(
        S2Point x, S2Point a, S2Point b, double r2, S2Point n, double n1, double n2) {
      if (r2 < DEG_45.getLength2()) {
        return triageLineSin2(x, a, b, r2, n, n1, n2);
      } else {
        return triageLineCos2(x, r2, n, n1, n2);
      }
    }

    /**
     * Returns -1, 0, or +1 according to whether the distance from "x" to the great circle through
     * (a, b) is less than, equal to, or greater than the given squared chord length "r2". This
     * method computes the squared sines of the distances involved, which is more accurate when the
     * distances are small (less than 45 degrees).
     *
     * <p>The remaining parameters are functions of (a, b) and are passed in because they have
     * already been computed: n = (a - b) x (a + b), n1 = n.norm(), and n2 = n.norm2().
     */
    static int triageLineSin2(
        S2Point x, S2Point a, S2Point b, double r2, S2Point n, double n1, double n2) {
      // The minimum distance is to a point on the edge interior. Since the true distance to the
      // edge is always less than 90 degrees, we can return immediately if the limit is 90 degrees
      // or larger.
      if (r2 >= 2.0) {
        // distance < limit
        return -1;
      }

      // Otherwise we compute sin^2(distance to edge) to get the best accuracy when the distance
      // limit is small (e.g., S2EdgeUtil.INTERSECTION_ERROR).
      double n2Sin2R = n2 * r2 * (1 - 0.25 * r2);
      double n2Sin2RError = 6 * T_ERR * n2Sin2R;
      double[] ax2 = {0};
      double xDn = x.sub(closestVertex(x, a, b, ax2)).dotProd(n);
      double xDn2 = xDn * xDn;
      double c1 = ((3.5 + 2 * sqrt(3)) * n1 + 32 * sqrt(3) * DBL_ERR) * T_ERR * sqrt(ax2[0]);
      double xDn2Error = 4 * T_ERR * xDn2 + (2 * abs(xDn) + c1) * c1;

      // If we are using extended precision, then it is worthwhile to recompute the length of X more
      // accurately. Otherwise we use the fact that X is guaranteed to be unit length to with a
      // tolerance of 4 * DBL_ERR.
      if (T_ERR < DBL_ERR) {
        // Note this is effectively dead code in Java, pending a long double data type.
        n2Sin2R *= x.norm2();
        n2Sin2RError += 4 * T_ERR * n2Sin2R;
      } else {
        n2Sin2RError += 8 * DBL_ERR * n2Sin2R;
      }
      return compare(xDn2, xDn2Error, n2Sin2R, n2Sin2RError);
    }

    /**
     * Like triageLineSin2, but this method computes the squared cosines of the distances involved.
     * It is more accurate when the distances are large (greater than 45 degrees).
     */
    static int triageLineCos2(S2Point x, double r2, S2Point n, double n1, double n2) {
      // The minimum distance is to a point on the edge interior. Since the true distance to the
      // edge is always less than 90 degrees, return if the limit is 90 degrees or larger.
      if (r2 >= 2.0) {
        // distance < limit
        return -1;
      }

      // Otherwise we compute cos^2(distance to edge).
      double cosR = 1 - 0.5 * r2;
      double n2Cos2R = n2 * cosR * cosR;
      double n2Cos2RError = 7 * DBL_ERR * n2Cos2R;

      // The length of M = X.crossProd(N) is the cosine of the distance.
      double m2 = x.crossProd(n).norm2();
      double m1 = sqrt(m2);
      double m1Error = ((1 + 8 / sqrt(3)) * n1 + 32 * sqrt(3) * DBL_ERR) * DBL_ERR;
      double m2Error = 3 * DBL_ERR * m2 + (2 * m1 + m1Error) * m1Error;

      // If we are using extended precision, then it is worthwhile to recompute the length of X more
      // accurately. Otherwise we use the fact that X is guaranteed to be unit length to within a
      // tolerance of 4 * DBL_ERR.
      if (T_ERR < DBL_ERR) {
        // Note this is effectively dead code in Java, pending a long double data type.
        n2Cos2R *= x.norm2();
        n2Cos2RError += 4 * T_ERR * n2Cos2R;
      } else {
        n2Cos2RError += 8 * DBL_ERR * n2Cos2R;
      }
      return compare(n2Cos2R, n2Cos2RError, m2, m2Error);
    }

    /** Returns a BigDecimal-based test result, which is exact but very slow. */
    public static int exact(S2Point x, S2Point a, S2Point b, double r2) {
      // Even if previous calculations were uncertain, we might not need to do *all* the
      // calculations in exact arithmetic here. For example it may be easy to determine whether "x"
      // is closer to an endpoint or the edge interior. The only calculation where we always use
      // exact arithmetic is when measuring the distance to the extended line (great circle) through
      // "a" and "b", since it is virtually certain that the previous floating point calculations
      // failed in that case.

      // CompareEdgeDirections also checks that no edge has antipodal endpoints.
      if (compareEdgeDirections(a, b, a, x) > 0 && compareEdgeDirections(a, b, x, b) > 0) {
        return exactLineInterior(big(x), big(a), big(b), big(r2));
      } else {
        return exactLineEndpoints(x, a, b, r2);
      }
    }

    /** Returns a BigDecimal-based test result assuming the projection of X is onto A or B. */
    static int exactLineEndpoints(S2Point x, S2Point a, S2Point b, double r2) {
      return min(compareDistance(x, a, r2), compareDistance(x, b, r2));
    }

    /** Returns a BigDecimal-based test assuming the projection of "x" is between A and B. */
    static int exactLineInterior(BigPoint x, BigPoint a, BigPoint b, BigDecimal r2) {
      // Since we are given that the closest point is in the edge interior, the true distance is
      // always less than 90 degrees (which corresponds to a squared chord length of 2.0).
      if (r2.compareTo(TWO) >= 0) {
        // distance < limit
        return -1;
      }

      // Otherwise compute the edge normal
      BigPoint n = a.crossProd(b);
      BigDecimal sinD = x.dotProd(n);
      BigDecimal sin2R = r2.multiply(BigDecimal.ONE.subtract(QUARTER.multiply(r2)));
      return square(sinD).compareTo(sin2R.multiply(x.norm2()).multiply(n.norm2()));
    }
  }

  /**
   * Returns -1, 0, or +1 according to whether the normal of edge AB has negative, zero, or positive
   * dot product with the normal of edge CD. This essentially measures whether the edges AB and CD
   * are closer to proceeding in the same direction or in opposite directions around the sphere.
   *
   * <p>This method returns an exact result, i.e. the result is zero if and only if the two edges
   * are exactly perpendicular or at least one edge is degenerate. (i.e., both edge endpoints
   * project to the same point on the sphere).
   *
   * <p>However, this method does not use symbolic perturbations. Therefore it can return zero even
   * when A != B and C != D, e.g. if A == S * B exactly for some constant S > 0 (which is possible
   * even when both points are considered "normalized").
   *
   * <p>Edges may not consist of antipodal points (e.g., A != -B). See {@link #compareEdgeDistance}.
   */
  // TODO(user): Make this public once the library has moved onto it.
  static int compareEdgeDirections(S2Point a, S2Point b, S2Point c, S2Point d) {
    // Check that no edge consists of antipodal points. This catches the most common case; a full
    // test is in CompareEdgeDirections.exact.)
    // assert !a.equals(b.neg());
    // assert !c.equals(d.neg());

    int sign = CompareEdgeDirections.triage(a, b, c, d);
    if (sign != 0) {
      return sign;
    }

    // Optimization for the case where either edge is degenerate.
    if (a.equalsPoint(b) || c.equalsPoint(d)) {
      return 0;
    }

    // TODO(user): Use Real instead of BigDecimal.
    return CompareEdgeDirections.exact(a, b, c, d);
  }

  /**
   * A test to compare whether two edges are closer to proceeding in the same direction or in
   * opposite directions around the sphere, essentially signum((AxB)x(CxD)). Returns -1, 0, or +1
   * according to whether the normal of edge AB has negative, zero, or positive dot product with the
   * normal of edge CD, and may return 0 if the relation is indeterminate.
   */
  static class CompareEdgeDirections {
    /** No instantiation. */
    private CompareEdgeDirections() {}

    /** Returns a cosine-based test result. Fast but has a wide margin of uncertainty. */
    public static int triage(S2Point a, S2Point b, S2Point c, S2Point d) {
      S2Point n1 = ndCross(a, b);
      S2Point n2 = ndCross(c, d);
      double len1 = n1.norm();
      double len2 = n2.norm();
      double cos = n1.dotProd(n2);
      double cosError =
          ((5 + 4 * sqrt(3)) * len1 * len2 + 32 * sqrt(3) * DBL_ERR * (len1 + len2)) * T_ERR;
      return signum(cos, cosError);
    }

    /** Returns a BigDecimal-based test result. Exact but very slow. */
    public static int exact(S2Point a, S2Point b, S2Point c, S2Point d) {
      return exact(big(a), big(b), big(c), big(d));
    }

    /** Returns a BigDecimal-based test result. Exact but very slow. */
    public static int exact(BigPoint a, BigPoint b, BigPoint c, BigPoint d) {
      // assert !a.isAntipodal(b);
      // assert !c.isAntipodal(d);
      return a.crossProd(b).dotProd(c.crossProd(d)).signum();
    }
  }

  /**
   * Returns sign(P, Q, Z) where Z is the circumcenter of triangle ABC. The return value is -1 if Z
   * is to the left of edge PQ, and +1 if Z is to the right of edge PQ. The return value is zero if
   * A == B, B == C, or C == A (exactly), and also if P and Q project to identical points on the
   * sphere (e.g., P == Q).
   *
   * <p>The result is determined with respect to the positions of all points as though they were
   * projected to lie exactly on the surface of the unit sphere. Furthermore this method uses
   * symbolic perturbations to compute a consistent non-zero result even when Z lies exactly on edge
   * PQ.
   *
   * <p>Requires that P and Q do not project to antipodal points (e.g., P == -Q) (see comments in
   * compareEdgeDistance).
   */
  public static int edgeCircumcenterSign(S2Point p, S2Point q, S2Point a, S2Point b, S2Point c) {
    // Check that the edge does not consist of antipodal points. This catches the most common case,
    // see EdgeCircumcenterSign.exact for a full test.
    // assert !p.equals(q.neg());

    int abc = sign(a, b, c);
    int sign = EdgeCircumcenterSign.triage(p, q, a, b, c, abc);
    if (sign != 0) {
      return sign;
    }

    // Optimization for the cases that are going to return zero anyway, in order to avoid falling
    // back to exact arithmetic.
    if (p.equalsPoint(q) || a.equalsPoint(b) || b.equalsPoint(c) || c.equalsPoint(a)) {
      return 0;
    }

    // TODO(user): Use Real instead of BigDecimal.
    sign = EdgeCircumcenterSign.exact(p, q, a, b, c, abc);
    if (sign != 0) {
      return sign;
    }

    return EdgeCircumcenterSign.sos(p, q, a, b, c);
  }

  /**
   * A predicate for whether an edge PQ passes to the left, to the right, or through the center of
   * the circumcircle of triangle ABC. Useful to determine the orientation of an edge with respect
   * to the centers of a Voronoi diagram. Returns sign(P, Q, Z) where Z is the circumcenter of
   * triangle ABC. The return value is -1 if Z is to the left of edge PQ, and +1 if Z is to the
   * right of edge PQ. The return value is zero if the triangle has two or more exactly duplicate
   * vertices, or if the result is indeterminate.
   */
  static class EdgeCircumcenterSign {
    /** No instantiation. */
    private EdgeCircumcenterSign() {}

    /** Returns a double-based test result. Faster but has a larger margin of uncertainty. */
    public static int triage(S2Point p, S2Point q, S2Point a, S2Point b, S2Point c, int abc) {
      // Compute the circumcenter Z of triangle ABC, and then test which side of edge PQ it lies on.
      double[] zError = {0};
      S2Point z = circumcenter(a, b, c, zError);
      S2Point nx = ndCross(p, q);

      // When triangle ABC sign is negative, we have computed -Z and the result should be negated.
      double result = abc * nx.dotProd(z);
      double zLen = z.norm();
      double nxLen = nx.norm();
      double nxError = ((1 + 2 * sqrt(3)) * nxLen + 32 * sqrt(3) * DBL_ERR) * T_ERR;
      double resultError = ((3 * T_ERR * nxLen + nxError) * zLen + zError[0] * nxLen);
      return signum(result, resultError);
    }

    /** Returns a BigDecimal-based test result. Exact but very slow. */
    public static int exact(S2Point p, S2Point q, S2Point a, S2Point b, S2Point c, int abc) {
      return exact(big(p), big(q), big(a), big(b), big(c), abc);
    }

    /** Returns a BigDecimal-based test result. Exact but very slow. */
    public static int exact(BigPoint p, BigPoint q, BigPoint a, BigPoint b, BigPoint c, int abc) {
      // Return zero if the edge PQ is degenerate. (Also see the comments in sosTest.)
      if (p.isLinearlyDependent(q)) {
        // assert p.dotProd(q).signum() > 0; // Antipodal edges not allowed.
        return 0;
      }

      // The simplest predicate for testing whether the sign is positive is
      //
      // (1)  (P x Q) . (|C|(A x B) + |A|(B x C) + |B|(C x A)) > 0
      //
      // where |A| denotes A.norm() and the expression after the "." represents the circumcenter of
      // triangle ABC. (This predicate is terrible from a numerical accuracy point of view, but
      // that doesn't matter since we are going to use exact arithmetic.)  This predicate also
      // assumes that triangle ABC is CCW (positive sign); we correct for that below.
      //
      // The only problem with evaluating this inequality is that computing |A|, |B| and |C|
      // requires square roots. To avoid this problem we use the standard technique of rearranging
      // the inequality to isolate at least one square root and then squaring both sides. We need
      // to repeat this process twice in order to eliminate all the square roots, which leads to a
      // polynomial predicate of degree 20 in the input arguments.
      //
      // Rearranging (1) we get
      //
      //      (P x Q) . (|C|(A x B) + |A|(B x C)) > |B|(P x Q) . (A x C)
      //
      // Before squaring we need to check the sign of each side. If the signs are different then we
      // know the result without squaring, and if the signs are both negative then after squaring
      // both sides we need to invert the result. Define
      //
      //      dAB = (P x Q) . (A x B)
      //      dBC = (P x Q) . (B x C)
      //      dCA = (P x Q) . (C x A)
      //
      // Then we can now write the inequality above as
      //
      // (2)  |C| dAB + |A| dBC > -|B| dCA
      //
      // The RHS of (2) is positive if dCA < 0, and the LHS of (2) is positive if
      // (|C| dAB + |A| dBC) > 0. Since the LHS has square roots, we need to eliminate them using
      // the same process. Rewriting the LHS as
      //
      // (3)  |C| dAB > -|A| dBC
      //
      // we again the signs of both sides. Let's start with that. We also precompute the following
      // values because they are used repeatedly when squaring various expressions below:
      //
      //     abc2 = |A|^2 dBC^2
      //     bca2 = |B|^2 dCA^2
      //     cab2 = |C|^2 dAB^2
      BigPoint nx = p.crossProd(q);
      BigDecimal dab = nx.dotProd(a.crossProd(b));
      BigDecimal dbc = nx.dotProd(b.crossProd(c));
      BigDecimal dca = nx.dotProd(c.crossProd(a));
      BigDecimal abc2 = a.norm2().multiply(square(dbc));
      BigDecimal bca2 = b.norm2().multiply(square(dca));
      BigDecimal cab2 = c.norm2().multiply(square(dab));

      // If the two sides of (3) have different signs (including the case where one side is zero)
      // then we know the result. Also, if both sides are zero then we know the result. The
      // following logic encodes this.
      int lhsSign3 = dab.signum();
      int rhs3Sign = -dbc.signum();
      int lhsSign2 = max(-1, min(1, lhsSign3 - rhs3Sign));
      if (lhsSign2 == 0 && lhsSign3 != 0) {
        // Both sides of (3) have the same non-zero sign, so square both sides. If both sides were
        // negative then invert the result.
        lhsSign2 = cab2.compareTo(abc2) * lhsSign3;
        // cab2=4, abc2=2, cab2-abc2=2, compareTo(4,2)=1
      }
      // Now if the two sides of (2) have different signs then the result of this function is known.
      int rhsSign2 = -dca.signum();
      int result = max(-1, min(1, lhsSign2 - rhsSign2));
      if (result == 0 && lhsSign2 != 0) {
        // Both sides of (2) have the same non-zero sign, so square both sides. If both sides were
        // negative then we invert the result below. This gives
        //
        //        |C|^2 dAB^2 + |A|^2 dBC^2 + 2 |A| |C| dAB dBC > |B|^2 dCA^2
        //
        // This expression still has square roots (|A| and |C|), so we rewrite as
        //
        // (4)    2 |A| |C| dAB dBC > |B|^2 dCA^2 - |C|^2 dAB^2 - |A|^2 dBC^2 .
        //
        // Again, if the two sides have different signs then we know the result.
        int lhsSign4 = dab.signum() * dbc.signum();
        BigDecimal rhs4 = bca2.subtract(cab2).subtract(abc2);
        result = max(-1, min(1, lhsSign4 - rhs4.signum()));
        if (result == 0 && lhsSign4 != 0) {
          // Both sides of (4) have the same non-zero sign, so square both sides. If both sides were
          // negative then invert the result.
          result = FOUR.multiply(abc2).multiply(cab2).compareTo(square(rhs4)) * lhsSign4;
        }
        // Correct the sign if both sides of (2) were negative.
        result *= lhsSign2;
      }
      // If the sign of triangle ABC is negative, we have computed -Z so negate the result.
      return abc * result;
    }

    /**
     * Given the exact test resulted in 0, returns a Simulation of Simplicity-based test result,
     * that can only result in zero if P == Q, A == B, B == C, or C == A (the result will be nonzero
     * if these pairs are exactly proportional to each other but not equal.)
     */
    public static int sos(S2Point p, S2Point q, S2Point a, S2Point b, S2Point c) {
      // We use the same perturbation strategy as SymbolicCompareDistances. Note that pedestal
      // perturbations of P and Q do not affect the result, because Sign(P, Q, Z) does not change
      // when its arguments are scaled by a positive factor. Therefore we only need to consider
      // A, B, C. Suppose that A is the smallest lexicographically and therefore has the largest
      // perturbation. This has the effect of perturbing the circumcenter of ABC slightly towards
      // A, and since the circumcenter Z was previously exactly collinear with edge PQ, this implies
      // that after the perturbation sign(P, Q, Z) == unperturbedSign(P, Q, A). We want the result
      // to be zero if P, Q, and A are linearly dependent, rather than using symbolic perturbations,
      // because these perturbations are defined to be much, much smaller than the pedestal
      // perturbation of B and C that are considered below.)
      //
      // If A is also exactly collinear with edge PQ, then we move on to the next smallest point
      // lexicographically out of {B, C}. It is easy to see that as long as A, B, C are all
      // distinct, one of these three Sign calls will be nonzero, because if A, B, C are all
      // distinct and collinear with edge PQ then their circumcenter Z coincides with the normal of
      // PQ, and therefore Sign(P, Q, Z) is nonzero.
      //
      // This function could be extended to handle the case where P and Q are linearly dependent as
      // follows. First, suppose that every point has both a pedestal peturbation as described
      // above, and also the three axis-aligned perturbations described in the "Simulation of
      // Simplicity" paper, where all pedestal perturbations are defined to be much, much larger
      // than any axis-aligned perturbation. Note that since pedestal perturbations have no effect
      // on 'sign', we can use this model for *all* the S2 predicates, which ensures that all the
      // various predicates are fully consistent with each other.
      //
      // With this model, the strategy described above yields the correct result unless P and Q are
      // exactly linearly dependent. When that happens, then no perturbation (pedestal or
      // axis-aligned) of A,B,C affects the result, and no pedestal perturbation of P or Q affects
      // the result, therefore we need to consider the smallest axis-aligned perturbation of P or Q.
      // The first perturbation that makes P and Q linearly independent yields the result.
      // Supposing that P < Q, this is the perturbation of P[2] unless both points are multiples of
      // [0, 0, 1], in which case it is the perturbation of P[1]. The sign test can be implemented
      // by computing the perturbed cross product of P and Q and taking the dot product with the
      // exact value of Z. For example if P[2] is perturbed, the perturbed cross product is
      // proportional to (0, 0, 1) x Q = (-Q[1], Q[0], 0). Note that if the dot product with Z is
      // exactly zero, then it is still necessary to fall back to pedestal perturbations of A, B, C,
      // but one of these perturbations is now guaranteed to succeed.

      // If any two triangle vertices are equal, the result is zero.
      if (a.equalsPoint(b) || b.equalsPoint(c) || c.equalsPoint(a)) {
        return 0;
      }

      // Sort A, B, C in lexicographic order.
      if (b.compareTo(a) < 0) {
        S2Point temp = a;
        a = b;
        b = temp;
      }
      if (c.compareTo(b) < 0) {
        S2Point temp = b;
        b = c;
        c = temp;
      }
      if (b.compareTo(a) < 0) {
        S2Point temp = a;
        a = b;
        b = temp;
      }

      // Now consider the perturbations in decreasing order of size.
      int sign = Sign.sign(p, q, a, false);
      if (sign != 0) {
        return sign;
      }
      sign = Sign.sign(p, q, b, false);
      if (sign != 0) {
        return sign;
      }
      return Sign.sign(p, q, c, false);
    }
  }

  /**
   * Given two sites A and B that form the center of caps of radius 'r', this indicates which sites
   * are irrelevant to the Voronoi diagram relative to an edge PQ.
   */
  public enum Excluded {
    /** The first site is excluded, i.e. always further from the edge than the second site. */
    FIRST,
    /** The second site is excluded, i.e. always further from the edge than the first site. */
    SECOND,
    /** Neither site is excluded, i.e. both sites are closer to part of the edge than the other. */
    NEITHER,
    /** The algorithm could not robustly determine the exclusion, or A is exactly equal to B. */
    UNCERTAIN
  }

  /**
   * This is a specialized method that is used to compute the intersection of an edge PQ with the
   * Voronoi diagram of a set of points, where each Voronoi region is intersected with a disc of
   * fixed radius "r".
   *
   * <p>Given two sites A and B and an edge (P, Q) such that {@code d(A,P) < d(B,P)}, and both sites
   * are within the given distance "r" of edge PQ, this method intersects the Voronoi region of each
   * site with a disc of radius r and determines whether either region has an empty intersection
   * with edge PQ. It returns FIRST if site A has an empty intersection, SECOND if site B has an
   * empty intersection, NEITHER if neither site has an empty intersection, or UNCERTAIN if A == B
   * exactly. Note that it is not possible for both intersections to be empty because of the
   * requirement that both sites are within distance r of edge PQ. (For example, the only reason
   * that Voronoi region A can have an empty intersection with PQ is that site B is closer to all
   * points on PQ that are within radius r of site A.)
   *
   * <p>The result is determined with respect to the positions of all points as though they were
   * projected to lie exactly on the surface of the unit sphere. Furthermore this method uses
   * symbolic perturbations to compute a consistent non-zero result even when A and B lie on
   * opposite sides of PQ such that the Voronoi edge between them exactly coincides with edge PQ, or
   * when A and B are distinct but project to the same point on the sphere (i.e., they are linearly
   * dependent).
   *
   * <p>Requires that
   *
   * <ul>
   *   <li>{@code r < S1ChordAngle.RIGHT} (90 degrees)
   *   <li>{@code compareDistances(p, a, b) < 0}
   *   <li>{@code compareEdgeDistance(a, p, q, r) <= 0}
   *   <li>{@code compareEdgeDistance(b, p, q, r) <= 0}
   *   <li>P and Q do not project to antipodal points (e.g., P != -Q), see {@link
   *       CompareEdgeDistance} for details.
   * </ul>
   */
  public static Excluded getVoronoiSiteExclusion(
      S2Point a, S2Point b, S2Point p, S2Point q, double r2) {
    // assert r2 < 2.0; // Less than a right angle.
    // assert compareDistances(p, a, b) < 0;  // (implies a != b)
    // assert compareEdgeDistance(a, p, q, r2) <= 0;
    // assert compareEdgeDistance(b, p, q, r2) <= 0;

    // Check that the edge does not consist of antipodal points. This catches the most common case,
    // see VoronoiSiteExclusion.exact for the full test.
    // assert !p.equals(q.neg());

    // If one site is closer than the other to both endpoints of PQ, then it is closer to every
    // point on PQ. Note that this also handles the case where A and B are equidistant from every
    // point on PQ (i.e., PQ is the perpendicular bisector of AB), because CompareDistances uses
    // symbolic perturbations to ensure that A or B is considered closer, in a consistent way. This
    // also ensures that the choice of A or B does not depend on the direction of PQ.
    if (compareDistances(q, a, b) < 0) {
      // Site A is closer to every point on PQ.
      return Excluded.SECOND;
    } else {
      Excluded result = VoronoiSiteExclusion.triage(a, b, p, q, r2);

      // TODO(user): Use Real instead of BigDecimal.
      return result != Excluded.UNCERTAIN ? result : VoronoiSiteExclusion.exact(a, b, p, q, r2);
    }
  }

  /**
   * A test for which (if any) of two Voronoi sites within R of an edge PQ are covered by the other.
   *
   * <p>Does not offer generally correct results for all inputs, so that multiple strategies may be
   * implemented for different classes of input.
   */
  static class VoronoiSiteExclusion {
    /** No instantiation. */
    private VoronoiSiteExclusion() {}

    /** An exact representation of a right angle. */
    static final BigDecimal R90 = big(S1ChordAngle.RIGHT.getLength2());

    /** A site exclusion test using double arithmetic. */
    public static Excluded triage(S2Point a, S2Point b, S2Point p, S2Point q, double r2) {
      // Define the "coverage disc" of a site S to be the disc centered at S with radius r (i.e.
      // squared chord angle length r2). Similarly, define the "coverage interval" of S along an
      // edge PQ to be the intersection of PQ with the coverage disc of S. The coverage interval
      // can be represented as the point at the center of the interval and an angle that measures
      // the semi-width or "radius" of the interval.
      //
      // To test whether site A excludes site B along the input edge PQ, we test whether the
      // coverage interval of A contains the coverage interval of B. Let "ra" and "rb" be the radii
      // (semi-widths) of the two intervals, and let "d" be the angle between their center points.
      //  Then "a" properly contains "b" if (ra - rb > d), and "b" contains "a" if (rb - ra > d).
      // Note that only one of these conditions can be true. Therefore we can determine whether one
      // site excludes the other by checking whether
      //
      // (1)   |rb - ra| > d
      //
      // and use the sign of (rb - ra) to determine which site is excluded.
      //
      // The actual code is based on the following. Let A1 and B1 be the unit vectors A and B
      // scaled by cos(r) (these points are inside the sphere). The planes perpendicular to OA1 and
      // OA2 cut off two discs of radius r around A and B. Now consider the two lines (inside the
      // sphere) where these planes intersect the plane containing the input edge PQ, and let A2 and
      // B2 be the points on these lines that are closest to A and B. The coverage intervals of A
      // and B can be represented as an interval along each of these lines, centered at A2 and B2.
      // Let P1 and P2 be the endpoints of the coverage interval for A, and let Q1 and Q2 be the
      // endpoints of the coverage interval for B. We can view each coverage interval as either a
      // chord through the sphere's interior, or as a segment of the original edge PQ (by projecting
      // the chord onto the sphere's surface).
      //
      // To check whether B's interval is contained by A's interval, we test whether both endpoints
      // of B's interval (Q1 and Q2) are contained by A's interval. E.g., we could test whether
      // Qi.dotProd(A2) > A2.norm2().
      //
      // However rather than constructing the actual points A1, A2, and so on, it turns out to be
      // more efficient to compute the sines and cosines ("components") of the various angles and
      // then use trigonometric identities. Predicate (1) can be expressed as
      //
      //      |sin(rb - ra)| > sin(d)
      //
      // provided that |d| <= Pi/2 (which must be checked), and then expanded to
      //
      // (2)  |sin(rb) cos(ra) - sin(ra) cos(rb)| > sin(d) .
      //
      // The components of the various angles can be expressed using dot and cross products based on
      // the construction above:
      //
      //   sin(ra) = sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2) / |aXn|
      //   cos(ra) = cos(r) |a| |n| / |aXn|
      //   sin(rb) = sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) / |bXn|
      //   cos(rb) = cos(r) |b| |n| / |bXn|
      //   sin(d)  = (aXb).n |n| / (|aXn| |bXn|)
      //   cos(d)  = (aXn).(bXn) / (|aXn| |bXn|)
      //
      // Also, the squared chord length r2 is equal to 4 * sin^2(r / 2), which yields the following
      // relationships:
      //
      //   sin(r)  = sqrt(r2 (1 - r2 / 4))
      //   cos(r)  = 1 - r2 / 2
      //
      // We then scale both sides of (2) by |aXn| |bXn| / |n| (in order to minimize the number of
      // calculations and to avoid divisions), which gives:
      //
      //    cos(r) ||a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) -
      //            |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)| > (aXb).n
      //
      // Furthermore we can substitute |a| = |b| = 1 (as long as this is taken into account in the
      // error bounds), yielding
      //
      // (3)   cos(r) |sqrt(sin^2(r) |n|^2 - |b.n|^2) -
      //               sqrt(sin^2(r) |n|^2 - |a.n|^2)| > (aXb).n
      //
      // The code below is more complicated than this because many expressions have been modified
      // for better numerical stability. For example, dot products between unit vectors are
      // computed using (x - y).dotProd(x + y), and the dot product between a point P and the normal
      // N of an edge X is measured using (P - Xi).dotProd(N) where Xi is the endpoint of X that is
      // closer to P.

      S2Point n = ndCross(p, q); // 2 * p.crossProd(q)
      double n2 = n.norm2();
      double n1 = sqrt(n2);
      // This factor is used in the error terms of dot products with "n" below.
      double dnError = ((3.5 + 2 * sqrt(3)) * n1 + 32 * sqrt(3) * DBL_ERR) * T_ERR;

      double cosR = 1 - 0.5 * r2;
      double sin2r = r2 * (1 - 0.25 * r2);
      double n2sin2r = n2 * sin2r;

      // "ra" and "rb" denote sin(ra) and sin(rb) after the scaling above.
      double[] d = {0};
      double aDn = a.sub(closestVertex(a, p, q, d)).dotProd(n);
      double aDn2 = aDn * aDn;
      double aDnError = dnError * sqrt(d[0]);
      double ra2 = n2sin2r - aDn2;
      double ra2Error =
          (8 * DBL_ERR + 4 * T_ERR) * aDn2
              + (2 * abs(aDn) + aDnError) * aDnError
              + 6 * T_ERR * n2sin2r;
      // This is the minimum possible value of ra2, which is used to bound the derivative of
      // sqrt(ra2) in computing raError below.
      double minRa2 = ra2 - ra2Error;
      if (minRa2 < 0) {
        return Excluded.UNCERTAIN;
      }
      double ra = sqrt(ra2);
      // Includes the ra2 subtraction error above.
      double raError = 1.5 * T_ERR * ra + 0.5 * ra2Error / sqrt(minRa2);

      double[] bx2 = {0};
      double bDn = b.sub(closestVertex(b, p, q, bx2)).dotProd(n);
      double bDn2 = bDn * bDn;
      double bDnError = dnError * sqrt(bx2[0]);
      double rb2 = n2sin2r - bDn2;
      double rb2Error =
          (8 * DBL_ERR + 4 * T_ERR) * bDn2
              + (2 * abs(bDn) + bDnError) * bDnError
              + 6 * T_ERR * n2sin2r;
      double minRb2 = rb2 - rb2Error;
      if (minRb2 < 0) {
        return Excluded.UNCERTAIN;
      }
      double rb = sqrt(rb2);
      // Includes the rb2 subtraction error above.
      double rbError = 1.5 * T_ERR * rb + 0.5 * rb2Error / sqrt(minRb2);

      // The sign of LHS(3) determines which site may be excluded by the other.
      double lhs3 = cosR * (rb - ra);
      double absLhs3 = abs(lhs3);
      double lhs3Error = cosR * (raError + rbError) + 3 * T_ERR * absLhs3;

      // Now we evaluate the RHS of (3), which is proportional to sin(d).
      S2Point aXb = ndCross(a, b); // 2 * a.crossProd(b)
      double aXb1 = aXb.norm();
      double sinD = 0.5 * aXb.dotProd(n);
      double sinDError =
          (4 * DBL_ERR + (2.5 + 2 * sqrt(3)) * T_ERR) * aXb1 * n1
              + 16 * sqrt(3) * DBL_ERR * T_ERR * (aXb1 + n1);

      // If LHS(3) is definitely less than RHS(3), neither site excludes the other.
      double result = absLhs3 - sinD;
      double resultError = lhs3Error + sinDError;
      if (result < -resultError) {
        return Excluded.NEITHER;
      }

      // Otherwise, before proceeding further we need to check that |d| <= Pi/2. In fact, |d| < Pi/2
      // is enough because of the requirement that r < Pi/2. The following expression represents
      // cos(d) after scaling; it is equivalent to (aXn).(bXn) but has about 30% less error.
      double cosD = a.dotProd(b) * n2 - aDn * bDn;
      double cosDError =
          ((8 * DBL_ERR + 5 * T_ERR) * abs(aDn) + aDnError) * abs(bDn)
              + (abs(aDn) + aDnError) * bDnError
              + (8 * DBL_ERR + 8 * T_ERR) * n2;
      if (cosD <= -cosDError) {
        return Excluded.NEITHER;
      }

      // Potential optimization: if the sign of cos(d) is uncertain, then instead we could check
      // whether cos(d) >= cos(r). Unfortunately this is fairly expensive since it requires
      // computing denominator |aXn||bXn| of cos(d) and the associated error bounds. In any case
      // this case is relatively rare so it seems better to punt.
      if (cosD < cosDError) {
        return Excluded.UNCERTAIN;
      }

      // Normally we have d > 0 because the sites are sorted so that A is closer to P and B is
      // closer to Q. However if the edge PQ is longer than Pi/2, and the sites A and B are beyond
      // its endpoints, then AB can wrap around the sphere in the opposite direction from PQ. In
      // this situation d < 0 but each site is closest to one endpoint of PQ, so neither excludes
      // the other.
      //
      // It turns out that this can happen only when the site that is further away from edge PQ is
      // less than 90 degrees away from whichever endpoint of PQ it is closer to. It is provable
      // that if this distance is less than 90 degrees, then it is also less than r2, and therefore
      // the Voronoi regions of both sites intersect the edge.
      if (sinD < -sinDError) {
        double r90 = S1ChordAngle.RIGHT.getLength2();
        // "ca" is negative if Voronoi region A definitely intersects edge PQ.
        int ca = (lhs3 < -lhs3Error) ? -1 : CompareDistance.triageCos(a, p, r90);
        int cb = (lhs3 > lhs3Error) ? -1 : CompareDistance.triageCos(b, q, r90);
        if (ca < 0 && cb < 0) {
          return Excluded.NEITHER;
        }
        if (ca <= 0 && cb <= 0) {
          return Excluded.UNCERTAIN;
        }
        if (absLhs3 <= lhs3Error) {
          return Excluded.UNCERTAIN;
        }
      } else if (sinD <= sinDError) {
        return Excluded.UNCERTAIN;
      }

      // Now we can finish checking the results of predicate (3).
      if (result <= resultError) {
        return Excluded.UNCERTAIN;
      }
      // assert absLhs3 > lhs3Error;
      return (lhs3 > 0) ? Excluded.FIRST : Excluded.SECOND;
    }

    public static Excluded exact(S2Point a, S2Point b, S2Point p, S2Point q, double r2) {
      return exact(big(a), big(b), big(p), big(q), big(r2));
    }

    /** A site exclusion test using BigDecimal arithmetic. */
    public static Excluded exact(BigPoint a, BigPoint b, BigPoint p, BigPoint q, BigDecimal r2) {
      // assert !p.isAntipodal(q);

      // Recall that one site excludes the other if
      //
      // (1)  |sin(rb - ra)| > sin(d)
      //
      // and that the sign of (rb - ra) determines which site is excluded; see the comments in
      // triage(). To evaluate this using exact arithmetic, we expand this to the same predicate as
      // before:
      //
      // (2)    cos(r) ||a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) -
      //                |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)| > (aXb).n
      //
      // We also need to verify that d <= Pi/2, which is implemented by checking that sin(d) >= 0
      // and cos(d) >= 0.
      //
      // To eliminate the square roots we use the standard technique of rearranging the inequality
      // to isolate at least one square root and then squaring both sides. We need to repeat this
      // process twice in order to eliminate all the square roots, which leads to a polynomial
      // predicate of degree 20 in the input arguments (i.e. degree 4 in each of "a", "b", "p", "q",
      // and "r2").
      //
      // Before squaring we need to check the sign of each side. We also check the condition that
      // cos(d) >= 0. Given what else we need to compute, it is cheaper use the identity
      //
      //   (aXn).(bXn) = (a.b) |n|^2 - (a.n)(b.n)
      BigPoint n = p.crossProd(q);
      BigDecimal n2 = n.norm2();
      BigDecimal aDn = a.dotProd(n);
      BigDecimal bDn = b.dotProd(n);
      int cosDSign = a.dotProd(b).multiply(n2).compareTo(aDn.multiply(bDn));
      if (cosDSign < 0) {
        return Excluded.NEITHER;
      }

      // Otherwise we continue evaluating the LHS of (2), defining
      //    sa = |b| sqrt(sin^2(r) |a|^2 |n|^2 - |a.n|^2)
      //    sb = |a| sqrt(sin^2(r) |b|^2 |n|^2 - |b.n|^2) .
      // The sign of the LHS of (2) (before taking the absolute value) determines which coverage
      // interval is larger and therefore which site is potentially being excluded.
      BigDecimal a2 = a.norm2();
      BigDecimal b2 = b.norm2();
      BigDecimal n2Sin2R = r2.multiply(BigDecimal.ONE.subtract(QUARTER.multiply(r2))).multiply(n2);
      BigDecimal sa2 = b2.multiply(n2Sin2R.multiply(a2).subtract(square(aDn)));
      BigDecimal sb2 = a2.multiply(n2Sin2R.multiply(b2).subtract(square(bDn)));
      int lhsSign2 = sb2.compareTo(sa2);

      // If the RHS of (2) is negative (corresponding to sin(d) < 0), then we need to consider the
      // possibility that the edge AB wraps around the sphere in the opposite direction from edge
      // PQ, with the result that neither site excludes the other; see triage().
      BigDecimal rhs2 = a.crossProd(b).dotProd(n);
      int rhsSign2 = rhs2.signum();
      if (rhsSign2 < 0) {
        int ca = (lhsSign2 < 0) ? -1 : CompareDistance.exact(a, p, R90);
        int cb = (lhsSign2 > 0) ? -1 : CompareDistance.exact(b, q, R90);
        if (ca <= 0 && cb <= 0) {
          return Excluded.NEITHER;
        }
        // assert ca != 1 || cb != 1;
        return ca == 1 ? Excluded.FIRST : Excluded.SECOND;
      }
      if (lhsSign2 == 0) {
        // If the RHS of (2) is zero as well (i.e., d == 0) then both sites are equidistant from
        // every point on edge PQ. This case requires symbolic perturbations, but it should already
        // have been handled in getVoronoiSiteExclusion(); see the call to CompareDistances.
        // assert rhsSign2 > 0;
        return Excluded.NEITHER;
      }

      // Next we square both sides of (2), yielding
      //
      //      cos^2(r) (sb^2 + sa^2 - 2 sa sb) > (aXb.n)^2
      //
      // which can be rearranged to give
      //
      // (3)  cos^2(r) (sb^2 + sa^2) - (aXb.n)^2 > 2 cos^2(r) sa sb .
      //
      // The RHS of (3) is always non-negative, but we still need to check the sign of the LHS.
      BigDecimal cosR = BigDecimal.ONE.subtract(HALF.multiply(r2));
      BigDecimal cos2R = square(cosR);
      BigDecimal lhs3 = cos2R.multiply(sa2.add(sb2)).subtract(square(rhs2));
      if (lhs3.signum() < 0) {
        return Excluded.NEITHER;
      }

      // Otherwise we square both sides of (3) to obtain:
      //
      // (4)  LHS(3)^2  >  4 cos^4(r) sa^2 sb^2
      BigDecimal lhs4 = square(lhs3);
      BigDecimal rhs4 = FOUR.multiply(square(cos2R)).multiply(sa2).multiply(sb2);
      int result = lhs4.compareTo(rhs4);
      if (result < 0) {
        return Excluded.NEITHER;
      }
      if (result == 0) {
        // We have |rb - ra| = d and d > 0. This implies that one coverage interval contains the
        // other, but not properly: the two intervals share a common endpoint. The distance from
        // each site to that point is exactly "r", therefore we need to use symbolic perturbations.
        // Recall that site A is considered closer to an equidistant point if and only if A > B.
        // Therefore if (rb > ra && A > B) or (ra > rb && B > A) then each site is closer to at
        // least one point and neither site is excluded.
        //
        // Ideally this logic would be in a separate SymbolicVoronoiSiteExclusion method for better
        // testing, but this is not convenient because it needs lhsSign (which requires exact
        // arithmetic to compute).
        if ((lhsSign2 > 0) == (a.compareTo(b) > 0)) {
          return Excluded.NEITHER;
        }
      }

      // At this point we know that one of the two sites is excluded. The sign of the LHS of (2)
      // (before the absolute value) determines which one.
      return lhsSign2 > 0 ? Excluded.FIRST : Excluded.SECOND;
    }
  }

  /**
   * Returns (a-b).crossProd(a+b), which eliminates almost all of the error due to "x" and "y" being
   * not quite unit length. This method is extremely accurate for small distances; the *relative*
   * error in the result is O(DBL_ERR) for distances as small as DBL_ERR.
   */
  // TODO(user): Perhaps update callers to use S2.robustCrossProd(b,a).neg() instead?
  private static S2Point ndCross(S2Point a, S2Point b) {
    return a.sub(b).crossProd(a.add(b));
  }

  /** Returns cos(XY). Requires that "x" and "y" are {@link S2#isUnitLength}. */
  private static double cosDistance(S2Point x, S2Point y) {
    return x.dotProd(y);
  }

  /** Returns the error in a value returned from {@link #cosDistance}. */
  private static double cosDistanceError(double x) {
    return 9.5 * DBL_ERR * abs(x) + 1.5 * DBL_ERR;
  }

  /** Returns sin^2(XY), where XY=x.angle(y). Requires "x" and "y" to be {@link S2#isUnitLength}. */
  private static double sin2Distance(S2Point x, S2Point y) {
    return 0.25 * ndCross(x, y).norm2();
  }

  /** Returns the error in a value returned from {@link #sin2Distance}. */
  private static double sin2DistanceError(double x) {
    return (21 + 4 * sqrt(3)) * DBL_ERR * x
        + 32 * sqrt(3) * DBL_ERR * DBL_ERR * sqrt(x)
        + 768 * DBL_ERR * DBL_ERR * DBL_ERR * DBL_ERR;
  }

  /** Returns the same result as {@link Math#signum}, or 0 if 'value' is within 'error' of 0. */
  private static int signum(double value, double error) {
    return value > error ? 1 : value < -error ? -1 : 0;
  }

  /**
   * Returns the same result as {@link Double#compare}, or 0 if 'a' and 'b' are within their
   * measurement errors of each other.
   *
   * @param a the first value
   * @param aError the true value of 'a' must be at least a-aError and at most a+aError
   * @param b the second value
   * @param bError the true value of 'b' must be at least b-bError and at most b+bError
   */
  private static int compare(double a, double aError, double b, double bError) {
    double diff = a - b;
    double error = aError + bError;
    return (diff > error) ? 1 : (diff < -error) ? -1 : 0;
  }

  /**
   * Returns "a" or "b", whichever is closer to "x". Also returns the squared distance from the
   * returned point to "x" in "dx2[0]".
   */
  private static S2Point closestVertex(S2Point x, S2Point a, S2Point b, double[] dx2) {
    double ax2 = a.getDistance2(x);
    double bx2 = b.getDistance2(x);
    if (ax2 < bx2 || (ax2 == bx2 && a.lessThan(b))) {
      dx2[0] = ax2;
      return a;
    } else {
      dx2[0] = bx2;
      return b;
    }
  }

  /**
   * If triangle ABC has positive sign, returns its circumcenter. If ABC has negative sign, returns
   * the negated circumcenter.
   */
  private static S2Point circumcenter(S2Point a, S2Point b, S2Point c, double[] error) {
    // We compute the circumcenter using the intersection of the perpendicular bisectors of AB and
    // BC. The formula is essentially
    //
    //    Z = ((A x B) x (A + B)) x ((B x C) x (B + C)),
    //
    // except that we compute the cross product (A x B) as (A - B) x (A + B) (and similarly for
    // B x C) since this is much more stable when the inputs are unit vectors.
    S2Point abDiff = a.sub(b);
    S2Point abSum = a.add(b);
    S2Point bcDiff = b.sub(c);
    S2Point bcSum = b.add(c);
    S2Point nab = abDiff.crossProd(abSum);
    double nabLen = nab.norm();
    double abLen = abDiff.norm();
    S2Point nbc = bcDiff.crossProd(bcSum);
    double nbcLen = nbc.norm();
    double bcLen = bcDiff.norm();
    S2Point mab = nab.crossProd(abSum);
    S2Point mbc = nbc.crossProd(bcSum);
    error[0] =
        ((16 + 24 * sqrt(3)) * T_ERR + 8 * DBL_ERR * (abLen + bcLen)) * nabLen * nbcLen
            + 128 * sqrt(3) * DBL_ERR * T_ERR * (nabLen + nbcLen)
            + 3 * 4096 * DBL_ERR * DBL_ERR * T_ERR * T_ERR;
    return mab.crossProd(mbc);
  }

  /** Returns a BigDecimal-based representation of 'p'. */
  private static BigPoint big(S2Point p) {
    return new BigPoint(p);
  }

  /** Returns a BigDecimal-based representation of 'v'. */
  private static BigDecimal big(double v) {
    return Platform.newBigDecimal(v);
  }

  /** Returns v*v. */
  private static BigDecimal square(BigDecimal v) {
    return v.multiply(v);
  }

  // Don't instantiate
  private S2Predicates() {}
}
