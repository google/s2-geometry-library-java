/*
 * Copyright 2022 Google Inc.
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

import static com.google.common.geometry.S2.MIN_NORM;
import static com.google.common.primitives.Doubles.max;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;

import com.google.common.annotations.VisibleForTesting;
import java.math.BigDecimal;
import java.util.Optional;

/**
 * Class that implements RobustCrossProd which will attempt to calculate the cross product of two
 * S2Points with increasingly precise but expensive methods, as required to get a reliable result.
 * Namely:
 *
 * <ol>
 *   <li>Stable cross product: calculate (b+a) x (b-a) which should be more stable.
 *   <li>Check if a == b. If so return some orthogonal S2Point.
 *   <li>Calculate the cross product using Real.java which allows for more precision when dealing
 *       with values that greatly differ in magnitude.
 *   <li>Resort to calculating the cross product using BigDecimal.
 *   <li>Lastly, calculate a symbolic cross product.
 * </ol>
 */
public final class S2RobustCrossProd {
  /**
   * Return an S2Point "c" that is orthogonal to the given unit-length S2Points "a" and "b". This
   * function is similar to a.crossProd(b) except that it does a better job of ensuring
   * orthogonality when "a" is nearly parallel to "b", and it returns a non-zero result even when a
   * == b or a == -b.
   *
   * <p>Note- robustCrossProd makes no claims about the magnitude of the resulting S2Point. If this
   * is important, you must call normalize() on the result. The result should always be scaled such
   * that you can call normalize() without risking underflow.
   *
   * <p>It satisfies the following properties (RCP == robustCrossProd):
   *
   * <ol>
   *   <li>RCP(a,b) != 0 for all a, b
   *   <li>RCP(b,a) == -RCP(a,b) unless a == b or a == -b
   *   <li>RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
   *   <li>RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
   * </ol>
   */
  public static S2Point robustCrossProd(S2Point a, S2Point b) {
    // The direction of a.crossProd(b) becomes unstable as (a + b) or (a - b) approaches zero. This
    // leads to situations where a.crossProd(b) is not very orthogonal to "a" and/or "b". We could
    // fix this using Gram-Schmidt, but we also want b.robustCrossProd(a) == -a.robustCrossProd(b).
    //
    // The easiest fix is to just compute the cross product of (b+a) and (b-a). Mathematically, this
    // cross product is exactly twice the cross product of "a" and "b", but it has the numerical
    // advantage that (b+a) and (b-a) are always perpendicular (since "a" and "b" are unit length).
    // This yields a result that is nearly orthogonal to both "a" and "b" even if these two values
    // differ only in the lowest bit of one component.
    //
    // If a == b, we just return any S2Point orthogonal to a. Otherwise, we attempt to calculate
    // the cross product of (b+a) and (b-a) where the components are first converted to Reals. If
    // the result is still too small, we resort to using BigDecimal and then
    // assert isUnitLength(a) && isUnitLength(b);
    Optional<S2Point> result = stableCrossProd(a, b);
    if (result.isPresent()) {
      return result.get();
    }
    return exactCrossProd(a, b);
  }

  /**
   * Attempts to calculate the cross product in four steps. If the inputs are equal, then return
   * some orthogonal S2Point. Otherwise attempts to compute a valid result using Reals followed by
   * BigDecimal. If fail, falls back to doing a symbolically perturbed cross product.
   */
  @VisibleForTesting
  static S2Point exactCrossProd(S2Point a, S2Point b) {
    // Handle the (a == b) case now, before doing expensive arithmetic. The only result that makes
    // sense mathematically is to return zero, but it turns out to reduce the number of special
    // cases in client code if we instead return an arbitrary orthogonal S2Point.
    if (a.equalsPoint(b)) {
      return S2.ortho(a);
    }

    Optional<S2Point> realResult = realCrossProd(a, b);
    if (realResult.isPresent()) {
      return realResult.get();
    }

    Optional<S2Point> bigDecimalCrossProd = bigDecimalCrossProd(a, b);
    if (bigDecimalCrossProd.isPresent()) {
      return bigDecimalCrossProd.get();
    }

    return symbolicCrossProd(a, b);
  }

  /**
   * Evaluates the cross product of unit-length S2Points "a" and "b" in a numerically stable way,
   * returning a non-empty Optional if the error in the result is guaranteed to be at most
   * ROBUST_CROSS_PROD_ERROR.
   */
  @VisibleForTesting
  static Optional<S2Point> stableCrossProd(S2Point a, S2Point b) {
    // We compute the cross product (a - b) x (a + b). Mathematically this is exactly twice the
    // cross product of "a" and "b", but it has the numerical advantage that (a - b) and (a + b)
    // are nearly perpendicular (since "a" and "b" are unit length). This yields a result that is
    // nearly orthogonal to both "a" and "b" even if these two values differ only very slightly.
    //
    // The maximum directional error in radians when this calculation is done in precision T (where
    // T is a floating-point type) is:
    //
    //   (1 + 2 * sqrt(3) + 32 * sqrt(3) * DBL_ERR / ||N||) * T_ERR
    //
    // where ||N|| is the norm of the result. To keep this error to at most
    // ROBUST_CROSS_PROD_ERROR, assuming this value is much less than 1, we need
    //
    //   (1 + 2 * sqrt(3) + 32 * sqrt(3) * DBL_ERR / ||N||) * T_ERR <= kErr
    //
    //   ||N|| >= 32 * sqrt(3) * DBL_ERR / (kErr / T_ERR - (1 + 2 * sqrt(3)))
    //
    // From this you can see that in order for this calculation to ever succeed in double precision,
    // we must have kErr > (1 + 2 * sqrt(3)) * DBL_ERR, which is
    // ROBUST_CROSS_PROD_ERROR == 6 * DBL_ERR (== 3 * DBL_EPSILON) in order to minimize the number
    // of cases where higher precision is needed; in particular, higher precision is only necessary
    // when "a" and "b" are closer than about 18 * DBL_ERR == 9 * DBL_EPSILON. (80-bit precision can
    // handle inputs as close as 2.5 * LDBL_EPSILON.)
    S2Point result = b.add(a).crossProd(b.sub(a));
    if (result.norm2() < MIN_NORM * MIN_NORM) {
      return Optional.empty();
    } else {
      return Optional.of(result);
    }
  }

  /**
   * Calculate the cross product by using Reals which allow for extended precision compared to
   * simple double operations. Incurs a speed hit.
   *
   * <p>An absent Optional is returned if the either the calculated cross product results in the
   * zero S2Point S2Point.ORIGIN or if we detect Real error overflow when doing Real.strictMul.
   */
  @VisibleForTesting
  static Optional<S2Point> realCrossProd(S2Point a, S2Point b) {
    try {
      S2Point sum = a.add(b);
      S2Point difference = b.sub(a);

      Real cx =
          Real.strictMul(sum.getY(), difference.getZ())
              .sub(Real.strictMul(sum.getZ(), difference.getY()));
      Real cy =
          Real.strictMul(sum.getZ(), difference.getX())
              .sub(Real.strictMul(sum.getX(), difference.getZ()));
      Real cz =
          Real.strictMul(sum.getX(), difference.getY())
              .sub(Real.strictMul(sum.getY(), difference.getX()));
      S2Point realResult =
          normalizableFromReal(new S2Point(cx.doubleValue(), cy.doubleValue(), cz.doubleValue()));

      if (realResult.equalsPoint(S2Point.ORIGIN) || !isNormalizable(realResult)) {
        return Optional.empty();
      }
      return Optional.of(realResult);
    } catch (ArithmeticException e) {
      return Optional.empty();
    }
  }

  /**
   * Calculate the cross product by using BigDecimal which allows for a vast increase in precision
   * compared to simple double operations and Reals. Incurs a large speed penalty.
   */
  @VisibleForTesting
  static Optional<S2Point> bigDecimalCrossProd(S2Point a, S2Point b) {
    BigPoint ba = new BigPoint(a);
    BigPoint bb = new BigPoint(b);
    BigPoint axb = ba.crossProd(bb);
    S2Point normalizable = normalizeFromBigDecimal(axb.x, axb.y, axb.z);
    if (normalizable.equalsPoint(S2Point.ORIGIN)) {
      return Optional.empty();
    }
    return Optional.of(normalizable);
  }

  @VisibleForTesting
  static S2Point symbolicCrossProd(S2Point a, S2Point b) {
    if (a.lessThan(b)) {
      return ensureNormalizable(symbolicCrossProdSorted(a, b));
    }
    return ensureNormalizable(symbolicCrossProdSorted(b, a)).neg();
  }

  /**
   * Returns the cross product of "a" and "b" after symbolic perturbations. (These perturbations
   * only affect the result if "a" and "b" are exactly collinear, e.g. if a == -b or a == (1+eps) *
   * b.) The result may not be normalizable (i.e., ensureNormalizable() should be called on the
   * result).
   */
  @VisibleForTesting
  static S2Point symbolicCrossProdSorted(S2Point a, S2Point b) {
    // The following code uses the same symbolic perturbation model as S2Predicates.sign. The
    // particular sequence of tests below was obtained using Mathematica (although it would be easy
    // to do it by hand for this simple case).
    //
    // Just like the function S2Predicates.sos() every input coordinate x[i] is assigned a symbolic
    // perturbation dx[i]. We then compute the cross product
    //
    //     (a + da).crossProd(b + db) .
    //
    // The result is a polynomial in the perturbation symbols. For example if we did this in one
    // dimension, the result would be
    //
    //     a * b + b * da + a * db + da * db
    //
    // where "a" and "b" have numerical values and "da" and "db" are symbols. In 3 dimensions the
    // result is similar except that the coefficients are 3-Vectors rather than scalars.
    //
    // Every possible S2Point has its own symbolic perturbation in each coordinate (i.e., there are
    // about 3 * 2**192 symbols). The magnitudes of the perturbations are chosen such that if x < y
    // lexicographically, the perturbations for "y" are much smaller than the perturbations for "x".
    // Similarly, the perturbations for the coordinates of a given point x are chosen such that
    // dx[0] is much smaller than dx[1] which is much smaller than dx[2]. Putting this together
    // with the fact the inputs to this function have been sorted so that a < b lexicographically,
    // this tells us that
    //
    //     da[2] > da[1] > da[0] > db[2] > db[1] > db[0]
    //
    // where each perturbation is so much smaller than the previous one that we don't even need to
    // consider it unless the coefficients of all previous perturbations are zero. In fact, each
    // succeeding perturbation is so small that we don't need to consider it unless the coefficient
    // of all products of the previous perturbations are zero. For example, we don't need to
    // consider the coefficient of db[1] unless the coefficient of db[2]*da[0] is zero.
    //
    // The follow code simply enumerates the coefficients of the perturbations (and products of
    // perturbations) that appear in the cross product above, in order of decreasing perturbation
    // magnitude. The first non-zero coefficient determines the result. The easiest way to
    // enumerate the coefficients in the correct order is to pretend that each perturbation is
    // some tiny value "eps" raised to a power of two:
    //
    // eps**    1      2      4      8     16     32
    //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]
    //
    // Essentially we can then just count in binary and test the corresponding subset of
    // perturbations at each step. So for example, we must test the coefficient of db[2]*da[0]
    // before db[1] because eps**12 > eps**16.
    if (b.getX() != 0 || b.getY() != 0) { // da[2]
      return new S2Point(-b.getY(), b.getX(), 0);
    }
    if (b.getZ() != 0) { // da[1]
      return new S2Point(b.getZ(), 0, 0); // Note that b[0] == 0.
    }
    // None of the remaining cases can occur in practice, because we can only get to this point if b
    // = (0, 0, 0). Nevertheless, even (0, 0, 0) has a well-defined direction under the symbolic
    // perturbation model.
    if (a.getX() != 0 || a.getY() != 0) { // db[2]
      return new S2Point(a.getY(), -a.getX(), 0);
    }
    // The following coefficient is always non-zero, so we can stop here.
    //
    // It may seem strange that we are returning (1, 0, 0) as the cross product without even looking
    // at the sign of a[2]. (Wouldn't you expect (0, 0, -1) x (0, 0, 0) and (0, 0, 1) x (0, 0, 0)
    // to point in opposite directions?) It's worth pointing out that in this function there is *no
    // relationship whatsoever* between the S2Points "a" and "-a", because the perturbations applied
    // to these S2Points may be entirely different. This is why the identity "RCP(-a, b) == -RCP(a,
    // b)" does not hold whenever "a" and "b" are linearly dependent (i.e., proportional). [As it
    // happens the two cross products above actually do point in opposite directions, but for
    // example (1, 1, 1) x (2, 2, 2) = (-2, 2, 0) and (-1, -1, -1) x (2, 2, 2) = (-2, 2, 0) do not.]
    return new S2Point(1, 0, 0); // db[2] * da[1]
  }

  /**
   * Scales an S2Point as necessary to ensure that the result can be normalized without loss of
   * precision due to floating-point underflow.
   *
   * <p>REQUIRES: p != (0, 0, 0)
   */
  private static S2Point ensureNormalizable(S2Point point) {
    if (isNormalizable(point)) {
      return point;
    }
    // We can't just scale by a fixed factor because the smallest representable double is 2**-1074,
    // so if we multiplied by 2**(1074 - 242) then the result might be so large that we couldn't
    // square it without overflow.
    //
    // Note that we must scale by a power of two to avoid rounding errors. The code below scales
    // "p" such that the largest component is in the range [1, 2).
    double pmax = max(abs(point.getX()), abs(point.getY()), abs(point.getZ()));
    double factor = pow(2, -1 - Platform.getExponent(pmax));
    return point.mul(factor);
  }

  /**
   * Returns a normalizable S2Point from an S2Point obtained by real calculations. It scales the
   * result as necessary to ensure that the result can be normalized without loss of precision due
   * to floating-point underflow. (This method doesn't actually call normalize() since that would
   * create additional error in situations where normalization is not necessary.)
   */
  private static S2Point normalizableFromReal(S2Point point) {
    if (isNormalizable(point)) {
      return point;
    }

    double largestAbs = max(Math.abs(point.getX()), Math.abs(point.getY()), Math.abs(point.getZ()));

    return new S2Point(
        point.getX() / largestAbs, point.getY() / largestAbs, point.getZ() / largestAbs);
  }

  /**
   * Returns a normalized S2Point from a BigDecimal S2Point. It scales the result as necessary to
   * ensure that the result can be normalized without loss of precision due to floating-point
   * underflow. (This method doesn't actually call normalize() since that would create additional
   * error in situations where normalization is not necessary.)
   */
  private static S2Point normalizeFromBigDecimal(BigDecimal x, BigDecimal y, BigDecimal z) {
    S2Point noScaling = new S2Point(x.doubleValue(), y.doubleValue(), z.doubleValue());
    if (isNormalizable(noScaling)) {
      return noScaling;
    }

    int xExponent = -x.stripTrailingZeros().scale() + (x.stripTrailingZeros().precision() - 1);
    int yExponent = -y.stripTrailingZeros().scale() + (y.stripTrailingZeros().precision() - 1);
    int zExponent = -z.stripTrailingZeros().scale() + (z.stripTrailingZeros().precision() - 1);

    int maxExponent = Integer.MIN_VALUE;
    if (x.signum() != 0) {
      maxExponent = max(maxExponent, xExponent);
    }

    if (y.signum() != 0) {
      maxExponent = max(maxExponent, yExponent);
    }

    if (z.signum() != 0) {
      maxExponent = max(maxExponent, zExponent);
    }

    if (maxExponent == Integer.MIN_VALUE) {
      return S2Point.ORIGIN;
    }

    // Subtract 1 from the negative of the maxExponent to ensure the largest value is between
    // 0.1 and 1.
    double scaledX = x.scaleByPowerOfTen(-maxExponent - 1).doubleValue();
    double scaledY = y.scaleByPowerOfTen(-maxExponent - 1).doubleValue();
    double scaledZ = z.scaleByPowerOfTen(-maxExponent - 1).doubleValue();

    return new S2Point(scaledX, scaledY, scaledZ);
  }

  /**
   * Returns true if the given S2Point's magnitude is large enough such that the angle to another
   * S2Point of the same magnitude can be measured using angle() without loss of precision due to
   * floating-point underflow. (This requirement is also sufficient to ensure that normalize() can
   * be called without risk of precision loss.)
   */
  private static boolean isNormalizable(S2Point point) {
    // Let ab = robustCrossProd(a, b) and cd = robustCrossProd(cd). In order for ab.angle(cd) to
    // not lose precision, the squared magnitudes of ab and cd must each be at least 2**-484. This
    // ensures that the sum of the squared magnitudes of ab.crossProd(cd) and ab.dotProd(cd) is at
    // least 2**-968, which ensures that any denormalized terms in these two calculations do not
    // affect the accuracy of the result (since all denormalized numbers are smaller than 2**-1022,
    // which is less than DBL_ERR * 2**-968).
    //
    // The fastest way to ensure this is to test whether the largest component of the result has a
    // magnitude of at least 2**-242.
    return max(Math.abs(point.getX()), Math.abs(point.getY()), Math.abs(point.getZ()))
        >= Math.pow(2, -242);
  }

  private S2RobustCrossProd() {}
}
