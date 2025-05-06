/*
 * Copyright 2014 Google Inc.
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

import static java.lang.Math.PI;
import static java.lang.Math.scalb;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.math.BigDecimal;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link Real}. */
@RunWith(JUnit4.class)
public class RealTest extends GeometryTestCase {

  /** Check for long-running summation error. */
  @Test
  public void testSerialArithmetic() {
    BigDecimal bigSum = new BigDecimal(0);
    Real realSum = new Real(0);
    BigDecimal bigProduct = new BigDecimal(1);
    Real realProduct = new Real(1);
    for (int i = 0; i < 100; i++) {
      double d = data.nextDouble();
      bigSum = bigSum.add(new BigDecimal(d));
      realSum = realSum.add(new Real(d));
      bigProduct = bigProduct.multiply(new BigDecimal(d));
      realProduct = realProduct.mul(d);
    }
    assertExactly(bigSum.doubleValue(), realSum.doubleValue());
    assertExactly(bigProduct.doubleValue(), realProduct.doubleValue());
  }

  /** Check adding all ranges of values, without risking overflow. */
  @Test
  public void testSplitArithmetic() {
    double mantissa = PI / 4;
    for (int a = -400; a <= 400; a += 16) {
      for (int b = -400; b <= 400; b += 16) {
        double d1 = scalb(mantissa, a);
        double d2 = scalb(mantissa, b);
        BigDecimal b1 = new BigDecimal(d1);
        BigDecimal b2 = new BigDecimal(d2);
        assertEquals(b1.add(b2).stripTrailingZeros(), Real.add(d1, d2).bigValue());
        assertEquals(b1.multiply(b2).stripTrailingZeros(), Real.mul(d1, d2).bigValue());
      }
    }
  }

  /**
   * Verifies that naive and robust techniques produce nearly equal results all around the sphere.
   */
  @Test
  public void testRandomCCW() {
    for (int i = 0; i < 1000; i++) {
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      S2Point c = data.getRandomPoint();
      double naiveDeterminant = S2Point.crossProd(a, b).dotProd(c);
      double robustDeterminant = det(a, b, c).doubleValue();
      assertEquals(
          "Robust result not close enough to naive result",
          naiveDeterminant,
          robustDeterminant,
          1e-15);
    }
  }

  /**
   * Verifies that two random points and their midpoint are basically never truly collinear, using
   * inexact arithmetic, but that all combinations of arguments return logically consistent sign
   * results (e.g. if [a,b,c] forms a left turn at b, then [c,b,a] must form a right turn at b.)
   */
  @Test
  public void testNearlyCollinearCCW() {
    int collinearPoints = 0;
    for (int i = 0; i < 1000; i++) {
      S2Point a = data.getRandomPoint();
      S2Point c = data.getRandomPoint();
      S2Point b = S2Point.normalize(S2Point.add(a, c));
      int sign = sign(a, b, c);
      if (sign == 0) {
        collinearPoints++;
        assertTrue("Rounding error should prevent exactly collinear points", collinearPoints < 10);
      } else {
        // ABC == BCA == CAB
        assertEquals("ABC != BCA", sign, sign(b, c, a));
        assertEquals("ABC != CAB", sign, sign(c, a, b));
        // ABC == -ACB == -BAC == -CBA
        assertEquals("ABC != -ACB", sign, -sign(a, c, b));
        assertEquals("ABC != -BAC", sign, -sign(b, a, c));
        assertEquals("ABC != -CBA", sign, -sign(c, b, a));
      }
    }
  }

  /**
   * Verifies that exactly collinear points have sign()==0, and form left or right turns if any
   * nudge is applied.
   */
  @Test
  public void testExactlyCollinearCCW() {
    double spacing = scalb((double) 1, -3);
    // For a variety of power-of-2 coordinates.
    for (double x1 = -1; x1 <= 1; x1 += spacing) {
      for (double y1 = -1; x1 <= 1; x1 += spacing) {
        for (double z1 = -1; x1 <= 1; x1 += spacing) {
          // Construct p.
          S2Point p = new S2Point(x1, y1, z1);
          // For a variety of power-of-2 coordinates.
          for (double x2 = -1; x2 <= 1; x2 += spacing) {
            for (double y2 = -1; x2 <= 1; x2 += spacing) {
              for (double z2 = -1; x2 <= 1; x2 += spacing) {
                // Construct q as the average of p and q. Since p and q are powers of 2 on all
                // coordinates, they're exact, and so is q.
                S2Point q = new S2Point(x2, y2, z2);

                // Skip p/q combinations that are coincident.
                if (S2Point.normalize(p).angle(S2Point.normalize(q)) < 1e-10) {
                  continue;
                }

                // Ensure that sign(PQR) is 0 for all combinations of PQR.
                S2Point r = S2Point.add(p, q);
                assertEquals(0, sign(p, q, r));
                assertEquals(0, sign(q, r, p));
                assertEquals(0, sign(r, p, q));
                assertEquals(0, sign(r, q, p));
                assertEquals(0, sign(q, p, r));
                assertEquals(0, sign(p, r, q));

                // Start with the PQ normal, and shrink it to the smallest vector that perturbs r.
                for (S2Point n = S2Point.crossProd(p, q);
                    !S2Point.add(r, n).equals(r);
                    n = S2Point.mul(n, .5)) {
                  // Ensure that the change in r has the expected change in sign(PQR).
                  assertEquals(1, sign(p, q, S2Point.add(r, n)));
                  assertEquals(-1, sign(p, q, S2Point.sub(r, n)));
                }
              }
            }
          }
        }
      }
    }
  }

  @Test
  public void testStrictMul() {
    // Check first situation where error is 0 but the output is equal to one of the operands.
    double a = -2.594E-321;
    double b = 0.9991685425907498;

    // Check that non-strict Real multiplication does not error.
    Real nonStrict = Real.mul(a, b);
    assertAlmostEquals(nonStrict.doubleValue(), a);

    // Check that strict Real multiplication does error.
    try {
      Real unused = Real.strictMul(a, b);
      fail();
    } catch (ArithmeticException e) {
      assertEquals("twoProductError underflowed", e.getMessage());
    }

    // Check first situation where error is 0 but the output is equal to one of the operands.
    double eps = 2.594E-321;
    Real nonStrictZero = Real.mul(eps, eps);
    assertAlmostEquals(nonStrictZero.doubleValue(), 0.0);

    try {
      Real unused = Real.strictMul(eps, eps);
      fail();
    } catch (ArithmeticException e) {
      assertEquals("twoProductError underflowed", e.getMessage());
    }

    // Check edge cases for operands being 1.0 or 0.0
    Real isZero = Real.strictMul(a, 0);
    assertIdentical(0.0, isZero.doubleValue());

    Real isA = Real.strictMul(a, 1.0);
    assertIdentical(isA.doubleValue(), a);

    // Check commutativity.
    try {
      Real unused = Real.strictMul(b, a);
      fail();
    } catch (ArithmeticException e) {
      assertEquals("twoProductError underflowed", e.getMessage());
    }
  }

  @Test
  public void testStrictMulScale() {
    Real realA = new Real(-2.594E-321);
    double scaleA = 0.9991685425907498;
    try {
      Real unused = realA.strictMul(scaleA);
      fail();
    } catch (ArithmeticException e) {
      assertEquals("twoProductError underflowed", e.getMessage());
    }

    // Check edge cases for operands being 1.0 or 0.0
    Real isZero = realA.strictMul(0.0);
    assertIdentical(0.0, isZero.doubleValue());

    Real isA = realA.strictMul(1.0);
    assertIdentical(isA.doubleValue(), realA.doubleValue());

    // Check commutativity
    Real realB = new Real(0.9991685425907498);
    double scaleB = -2.594E-321;
    try {
      Real unused = realB.strictMul(scaleB);
      fail();
    } catch (ArithmeticException e) {
      assertEquals("twoProductError underflowed", e.getMessage());
    }
  }

  private static Real det(S2Point a, S2Point b, S2Point c) {
    Real bycz = Real.mul(b.y, c.z);
    Real bzcy = Real.mul(b.z, c.y);
    Real bzcx = Real.mul(b.z, c.x);
    Real bxcz = Real.mul(b.x, c.z);
    Real bxcy = Real.mul(b.x, c.y);
    Real bycx = Real.mul(b.y, c.x);
    Real bcx = bycz.sub(bzcy);
    Real bcy = bzcx.sub(bxcz);
    Real bcz = bxcy.sub(bycx);
    Real x = bcx.mul(a.x);
    Real y = bcy.mul(a.y);
    Real z = bcz.mul(a.z);
    return x.add(y).add(z);
  }

  private static int sign(S2Point a, S2Point b, S2Point c) {
    return det(a, b, c).signum();
  }

}
