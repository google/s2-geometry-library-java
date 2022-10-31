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

import com.google.common.collect.Lists;
import java.math.BigDecimal;
import java.util.Iterator;
import java.util.List;

/** Unit tests for {@link Real}. */
public strictfp class RealTest extends GeometryTestCase {
  public void benchmark() {
    // Naive: 0.07320955800000001 sec
    // Robust: 1.408399374 sec
    // Biggie D: 8.373768562 sec

    @SuppressWarnings("unused")
    int total = 0;
    int reps = 1000000;
    long start;

    List<S2Point> points = Lists.newArrayList();
    for (int i = 0; i < reps; i++) {
      S2Point a = data.getRandomPoint();
      S2Point c = data.getRandomPoint();
      S2Point b = S2Point.normalize(S2Point.add(a, c));
      points.add(a);
      points.add(b);
      points.add(c);
    }

    Iterator<S2Point> it = points.iterator();
    start = System.nanoTime();
    while (it.hasNext()) {
      S2Point a = it.next();
      S2Point b = it.next();
      S2Point c = it.next();
      total += (int) Math.signum(S2Point.crossProd(a, b).dotProd(c));
    }
    System.out.println("Naive: " + 1e-9 * (System.nanoTime() - start) + " sec");

    it = points.iterator();
    start = System.nanoTime();
    while (it.hasNext()) {
      S2Point a = it.next();
      S2Point b = it.next();
      S2Point c = it.next();
      total += sign(a, b, c);
    }
    System.out.println("Robust: " + 1e-9 * (System.nanoTime() - start) + " sec");

    it = points.iterator();
    start = System.nanoTime();
    while (it.hasNext()) {
      S2Point a = it.next();
      S2Point b = it.next();
      S2Point c = it.next();
      total += bigSign(a, b, c);
    }
    System.out.println("Biggie D: " + 1e-9 * (System.nanoTime() - start) + " sec");
  }

  /** Check for long-running summation error. */
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
  public void testExactlyCollinearCCW() {
    double spacing = scalb(1, -3);
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

  private static int bigSign(S2Point a, S2Point b, S2Point c) {
    BigDecimal ax = new BigDecimal(a.x);
    BigDecimal ay = new BigDecimal(a.y);
    BigDecimal az = new BigDecimal(a.z);
    BigDecimal bx = new BigDecimal(b.x);
    BigDecimal by = new BigDecimal(b.y);
    BigDecimal bz = new BigDecimal(b.z);
    BigDecimal cx = new BigDecimal(c.x);
    BigDecimal cy = new BigDecimal(c.y);
    BigDecimal cz = new BigDecimal(c.z);
    BigDecimal bycz = by.multiply(cz);
    BigDecimal bzcy = bz.multiply(cy);
    BigDecimal bzcx = bz.multiply(cx);
    BigDecimal bxcz = bx.multiply(cz);
    BigDecimal bxcy = bx.multiply(cy);
    BigDecimal bycx = by.multiply(cx);
    BigDecimal bcx = bycz.subtract(bzcy);
    BigDecimal bcy = bzcx.subtract(bxcz);
    BigDecimal bcz = bxcy.subtract(bycx);
    BigDecimal x = bcx.multiply(ax);
    BigDecimal y = bcy.multiply(ay);
    BigDecimal z = bcz.multiply(az);
    return x.add(y).add(z).signum();
  }
}
