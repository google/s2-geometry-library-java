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

import static com.google.common.geometry.S2RobustCrossProd.exactCrossProd;
import static com.google.common.geometry.S2RobustCrossProd.robustCrossProd;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Optional;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Verifies S2RobustCrossProd static methods. */
@RunWith(JUnit4.class)
@SuppressWarnings("UseCorrectAssertInTests")
public class S2RobustCrossProdTest extends GeometryTestCase {

  private static final S1ChordAngle MAX_ERROR =
      S1ChordAngle.fromS1Angle(
          S2.ROBUST_CROSS_PROD_ERROR
              .add(S2.EXACT_CROSS_PROD_ERROR)
              .add(S1Angle.radians(2 * S2.DBL_ERROR)));

  enum PrecisionLevel {
    DOUBLE,
    REAL,
    BIG_DECIMAL,
    SYMBOLIC
  };

  // Checks that the robustCrossProd result at one level of precision is consistent with the result
  // at the next higher level of precision. Also verifies that the results satisfy the identities
  // documented in the main file.
  private static PrecisionLevel checkRobustCrossProdError(S2Point a, S2Point b) {
    S2Point result = robustCrossProd(a, b).normalize();
    S2Point offset = result.mul(S2.ROBUST_CROSS_PROD_ERROR.radians());
    S2Point a90 = result.crossProd(a);
    assertEquals(1, S2Predicates.sign(a, b, result));
    assertGreaterThan(result.dotProd(a.add(offset)), 0.0);

    assertLessThan(result.dotProd(a.sub(offset)), 0.0);
    assertGreaterThan(result.dotProd(a90.add(offset)), 0.0);
    assertLessThan(result.dotProd(a90.sub(offset)), 0.0);

    // The following identities are true unless a, b are linearly dependent.
    Optional<S2Point> resultBigDecimal = S2RobustCrossProd.bigDecimalCrossProd(a, b);
    boolean haveBigDecimal = resultBigDecimal.isPresent();
    if (haveBigDecimal) {
      assertEquals(robustCrossProd(a.neg(), b).normalize(), result.neg());
      assertEquals(robustCrossProd(a, b.neg()).normalize(), result.neg());
    }

    S2Point resultExact;
    if (a.equalsPoint(b)) {
      // This case isn't handled by exactCrossProd() because it is supposed to be handled long
      // before we get to that point.
      resultExact = S2.ortho(a).normalize();
    } else {
      resultExact = exactCrossProd(a, b).normalize();
      assertEquals(robustCrossProd(b, a).normalize(), result.neg());
    }

    Optional<S2Point> stableResultOptional = S2RobustCrossProd.stableCrossProd(a, b);
    Optional<S2Point> realResultOptional = S2RobustCrossProd.realCrossProd(a, b);
    Optional<S2Point> bigDecimalResultOptional = S2RobustCrossProd.bigDecimalCrossProd(a, b);

    if (stableResultOptional.isPresent()) {
      S2Point stableResult = stableResultOptional.get().normalize();
      assertEquals(result, stableResult);
      assertLessThan(
          S2Predicates.compareDistance(stableResult, resultExact, MAX_ERROR.getLength2()), 0);
      return PrecisionLevel.DOUBLE;
    } else if (realResultOptional.isPresent()) {
      assertEquals(result, realResultOptional.get().normalize());
      assertLessThan(
          S2Predicates.compareDistance(
              realResultOptional.get().normalize(), resultExact, MAX_ERROR.getLength2()),
          0);
      return PrecisionLevel.REAL;
    } else if (bigDecimalResultOptional.isPresent()) {
      assertEquals(result, bigDecimalResultOptional.get().normalize());
      return PrecisionLevel.BIG_DECIMAL;
    } else {
      assertEquals(result, resultExact);
      return PrecisionLevel.SYMBOLIC;
    }
  }

  private static void checkRobustCrossProd(
      S2Point a, S2Point b, S2Point expectedResult, PrecisionLevel expectedPrecision) {
    assertEquals(expectedPrecision, checkRobustCrossProdError(a, b));
    assertEquals(1, S2Predicates.sign(a, b, expectedResult));
    assertEquals(expectedResult, robustCrossProd(a, b).normalize());
    assertEquals(expectedResult.neg(), robustCrossProd(b, a).normalize());
  }

  @Test
  public void testRobustCrossProdCoverage() {
    // Note that robustCrossProd() returns a non-unit length result. In "double" precision the
    // magnitude is twice the usual cross product. In exact precision it is equal to the usual cross
    // product unless the magnitude is so small that calling normalize() would lose precision, in
    // which case it is scaled up appropriately. When symbolic perturbations are used, the result
    // magnitude is arbitrary. For this reason we only check the result after calling normalize().
    checkRobustCrossProd(
        new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(0, 0, 1), PrecisionLevel.DOUBLE);
    checkRobustCrossProd(
        new S2Point(20 * S2.DBL_ERROR, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 0, 1),
        PrecisionLevel.DOUBLE);

    checkRobustCrossProd(
        new S2Point(16 * S2.DBL_ERROR, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 0, 1),
        PrecisionLevel.DOUBLE);

    checkRobustCrossProd(
        new S2Point(4.9E-324, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 0, 1),
        PrecisionLevel.REAL);

    checkRobustCrossProd(
        new S2Point(4.9E-324, 1, 0),
        new S2Point(4.9E-324, 1 - S2.DBL_ERROR, 0),
        new S2Point(0, 0, -1),
        PrecisionLevel.BIG_DECIMAL);

    checkRobustCrossProd(
        new S2Point(1, 0, 0),
        new S2Point(1 + S2.DBL_EPSILON, 0, 0),
        new S2Point(0, 1, 0),
        PrecisionLevel.SYMBOLIC);

    checkRobustCrossProd(
        new S2Point(0, 1 + S2.DBL_EPSILON, 0),
        new S2Point(0, 1, 0),
        new S2Point(1, 0, 0),
        PrecisionLevel.SYMBOLIC);
    checkRobustCrossProd(
        new S2Point(0, 0, 1),
        new S2Point(0, 0, -1),
        new S2Point(-1, 0, 0),
        PrecisionLevel.SYMBOLIC);

    assertEquals(
        S2RobustCrossProd.symbolicCrossProd(new S2Point(-1, 0, 0), new S2Point(0, 0, 0)),
        new S2Point(0, 1, 0));
    assertEquals(
        S2RobustCrossProd.symbolicCrossProd(new S2Point(0, 0, 0), new S2Point(0, -1, 0)),
        new S2Point(1, 0, 0));
    assertEquals(
        S2RobustCrossProd.symbolicCrossProd(new S2Point(0, 0, 0), new S2Point(0, 0, -1)),
        new S2Point(-1, 0, 0));
  }

  @Test
  public void testSymbolicCrossProdConsistentWithSign() {
    // Tests that robustCrossProd() is consistent with S2Predicates.sign() when symbolic
    // perturbations are necessary. This implies that the two vectors A and B must be linearly
    // dependent. We test all possible orderings of the components of vector A = (x, y, z) and all
    // possible orderings of the two vectors A and B.
    double[] ternary = {-1, 0, 1};
    double[] scales = {-1.0, 1 - S2.DBL_ERROR, 1 + 2 * S2.DBL_ERROR};
    for (double x : ternary) {
      for (double y : ternary) {
        for (double z : ternary) {
          S2Point a = new S2Point(x, y, z).normalize();
          if (a.equals(S2Point.ZERO)) {
            continue;
          }
          for (double scale : scales) {
            S2Point b = a.mul(scale);
            assertTrue(S2.isUnitLength(a));
            assertTrue(S2.isUnitLength(b));
            assertGreaterThan(
                S2Predicates.sign(a, b, S2RobustCrossProd.robustCrossProd(a, b).normalize()), 0);
          }
        }
      }
    }
  }

  @Test
  public void testRobustCrossProdMagnitude() {
    // Test that angles can be measured between vectors returned by robustCrossProd without loss of
    // precision due to underflow. The following magnitude (1e-100) is demonstrates that it is not
    // sufficient to simply ensure that the result of robustCrossProd() is normalizable, since in
    // this  case normalize() could be called on the result of both cross products without precision
    // loss but it is still not possible to measure the angle between them without scaling the
    // result.
    assertEquals(
        S2.M_PI_2,
        S2RobustCrossProd.robustCrossProd(new S2Point(1, 0, 0), new S2Point(1, 1e-100, 0))
            .angle(
                S2RobustCrossProd.robustCrossProd(new S2Point(1, 0, 0), new S2Point(1, 0, 1e-100))),
        0.0);

    assertEquals(
        S2.M_PI_2,
        S2RobustCrossProd.robustCrossProd(new S2Point(-1e-100, 0, 1), new S2Point(1e-100, 0, -1))
            .angle(
                S2RobustCrossProd.robustCrossProd(
                    new S2Point(0, -1e-100, 1), new S2Point(0, 1e-100, -1))),
        0.0);
  }

  // Chooses a random S2Point that is often near the intersection of one of the coordinates planes
  // or coordinate axes with the unit sphere. (It is possible to represent very small perturbations
  // near such points.)
  S2Point choosePoint() {
    S2Point initialPoint = data.getRandomPoint();
    double[] coords = {initialPoint.getX(), initialPoint.getY(), initialPoint.getZ()};
    for (int i = 0; i < 3; i++) {
      if (data.oneIn(4)) {
        coords[i] = coords[i] * Math.pow(2, -1022 - 53 * data.nextDouble());
      } else if (data.oneIn(3)) {
        coords[i] = coords[i] * Math.pow(2, -511 - 511 * data.nextDouble());
      } else if (data.oneIn(2)) {
        coords[i] = coords[i] * Math.pow(2, -100 * data.nextDouble());
      }
    }
    S2Point f = new S2Point(coords[0], coords[1], coords[2]);
    if (f.norm2() >= Math.pow(2, -968)) {
      return f.normalize();
    }
    return choosePoint();
  }

  // Perturbs the length of the given point slightly while ensuring that it still satisfies the
  // conditions of S2.isUnitLength.
  S2Point perturbLength(S2Point p) {
    S2Point q = p.mul(data.uniform(1 - 2 * S2.DBL_EPSILON, 1 + 2 * S2.DBL_EPSILON));
    BigDecimal bdX = BigDecimal.valueOf(q.getX());
    BigDecimal bdY = BigDecimal.valueOf(q.getY());
    BigDecimal bdZ = BigDecimal.valueOf(q.getZ());
    BigDecimal norm2 = bdX.pow(2).add(bdY.pow(2)).add(bdZ.pow(2));
    BigDecimal norm2Minus1 = norm2.subtract(BigDecimal.ONE);
    BigDecimal absNorm2Minus1 = norm2Minus1.abs();
    BigDecimal epsilonFour = BigDecimal.valueOf(S2.DBL_EPSILON * 4);
    // S2.isUnitLength() is not an exact test, since it tests the length using ordinary
    // double-precision arithmetic and allows for errors in its own calculations. This is fine for
    // most purposes, but since here we are explicitly trying to push the limits we instead use
    // exact arithmetic to test whether the point is within the error tolerance guaranteed by
    // normalize().
    if (absNorm2Minus1.compareTo(epsilonFour) < 0) {
      return q;
    }
    return p;
  }

  S2Point getPointOnRay(S2Point origin, S2Point dir, S1Angle r) {
    assertTrue(S2.isUnitLength(origin));
    assertTrue(S2.isUnitLength(dir));
    assertLessThan(
        origin.dotProd(dir), S2.ROBUST_CROSS_PROD_ERROR.radians() + 0.75 * S2.DBL_EPSILON);
    return origin.mul(r.cos()).add(dir.mul(r.sin())).normalize();
  }

  S2Point getPointOnLine(S2Point a, S2Point b, S1Angle r) {
    S2Point dir = S2RobustCrossProd.robustCrossProd(a, b).crossProd(a).normalize();
    return getPointOnRay(a, dir, r);
  }

  @Test
  public void testRobustCrossProdError() {
    // We repeatedly choose two points (usually linearly dependent, or almost so) and measure the
    // distance between the cross product returned by robustCrossProd and one returned by
    // exactCrossProd.
    //
    // NOTE: For this Java implementation, there is an issue at i = 10406 where
    // A = (-4.4E-323, -0.9997890404674511, 0.020539585223999478)
    // B = (1.7993361E-316, -0.9997890404674511, 0.020539585223999478)
    // Specifically, the test for `assertLessThan(result.dotProd(a.sub(offset)), 0.0)` in
    // checkRobustCrossProdError.
    HashMap<PrecisionLevel, Integer> counts = new HashMap<>();
    for (int i = 0; i < 10000; i++) {
      System.out.println(i);
      data.setSeed(i + 1);
      S2Point a;
      S2Point b;
      do {
        a = perturbLength(choosePoint());
        S2Point dir = choosePoint();
        S1Angle r = S1Angle.radians(S2.M_PI_2 * Math.pow(2, -53 * data.nextDouble()));
        if (data.oneIn(3)) {
          r = r.mul(Math.pow(2, -1022 * data.nextDouble()));
        }
        b = perturbLength(getPointOnLine(a, dir, r));
        if (data.oneIn(2)) {
          b = b.neg();
        }
      } while (b.equals(a));
      PrecisionLevel unused = checkRobustCrossProdError(a, b);
      counts.computeIfPresent(unused, (k, v) -> v + 1);
    }
  }
}
