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

import static com.google.common.geometry.GeometryTestCase.assertExactly;
import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.NaN;
import static java.lang.Double.POSITIVE_INFINITY;

import com.google.common.annotations.GwtIncompatible;
import java.math.BigDecimal;
import java.math.MathContext;
import junit.framework.TestCase;

/** Unit tests for {@link Platform} methods. */
public class PlatformTest extends TestCase {
  private static final Double NEGATIVE_ZERO = Double.parseDouble("-0");

  public void testExp() {
    assertEquals(-1023, Platform.getExponent(-0.000000));
    assertEquals(1, Platform.getExponent(-3.141593));
    assertEquals(3, Platform.getExponent(-12.566371));
    assertEquals(4, Platform.getExponent(-28.274334));
    assertEquals(5, Platform.getExponent(-50.265482));
    assertEquals(6, Platform.getExponent(-78.539816));
    assertEquals(6, Platform.getExponent(-113.097336));
    assertEquals(7, Platform.getExponent(-153.938040));
    assertEquals(7, Platform.getExponent(-201.061930));
    assertEquals(7, Platform.getExponent(-254.469005));
    assertEquals(1024, Platform.getExponent(POSITIVE_INFINITY));
    assertEquals(-1023, Platform.getExponent(0.000000));
    assertEquals(1, Platform.getExponent(3.141593));
    assertEquals(3, Platform.getExponent(12.566371));
    assertEquals(4, Platform.getExponent(28.274334));
    assertEquals(5, Platform.getExponent(50.265482));
    assertEquals(6, Platform.getExponent(78.539816));
    assertEquals(6, Platform.getExponent(113.097336));
    assertEquals(7, Platform.getExponent(153.938040));
    assertEquals(7, Platform.getExponent(201.061930));
    assertEquals(7, Platform.getExponent(254.469005));
    assertEquals(1024, Platform.getExponent(NEGATIVE_INFINITY));
    assertEquals(1024, Platform.getExponent(NaN));
  }

  public void testRemainder() {
    double[] numerators = {0, 1, -2, 7, Math.E, NaN, NEGATIVE_INFINITY};
    double[] denominators = {0, -3, 4, Math.PI, NaN, POSITIVE_INFINITY};
    // Expected results from the remainder of division of each [numerator, denominator] pair defined
    // above.
    double[] results = {
      NaN,
      0,
      0,
      0,
      NaN,
      0,
      NaN,
      1.0,
      1.0,
      1.0,
      NaN,
      1.0,
      NaN,
      1.0,
      -2.0,
      1.1415926535897931,
      NaN,
      -2.0,
      NaN,
      1.0,
      -1.0,
      0.7168146928204138,
      NaN,
      7.0,
      NaN,
      -0.2817181715409549,
      -1.281718171540955,
      -0.423310825130748,
      NaN,
      2.718281828459045,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN,
      NaN
    };
    int resultIndex = 0;
    for (double f1 : numerators) {
      for (double f2 : denominators) {
        double expected = results[resultIndex++];
        double actual = Platform.IEEEremainder(f1, f2);
        // Note that we can't just use assertEquals, since the GWT JUnit version returns false
        // for assertTrue(NaN, NaN).
        assertTrue(Double.isNaN(expected) ? Double.isNaN(actual) : expected == actual);
      }
    }
  }

  public void testUlp() {
    assertExactly(1.2689709186578246e-116, Platform.ulp(1e-100));
    assertExactly(1.2924697071141057E-26, Platform.ulp(1e-10));
    assertExactly(4.440892098500626e-16, Platform.ulp(Math.PI));
    assertExactly(1.9073486328125e-6, Platform.ulp(1e10));
    assertExactly(1.942668892225729e84, Platform.ulp(1e100));
  }

  public void testNextAfter() {
    // NaNs
    assertTrue(Double.isNaN(Platform.nextAfter(NaN, NaN)));
    assertTrue(Double.isNaN(Platform.nextAfter(NaN, NEGATIVE_INFINITY)));
    assertTrue(Double.isNaN(Platform.nextAfter(NaN, -MAX_VALUE)));
    assertTrue(Double.isNaN(Platform.nextAfter(NEGATIVE_INFINITY, NaN)));

    // -inf
    assertExactly(NEGATIVE_INFINITY, Platform.nextAfter(NEGATIVE_INFINITY, NEGATIVE_INFINITY));
    assertExactly(-MAX_VALUE, Platform.nextAfter(NEGATIVE_INFINITY, -MAX_VALUE));
    assertExactly(-MAX_VALUE, Platform.nextAfter(NEGATIVE_INFINITY, Math.PI));
    assertExactly(-MAX_VALUE, Platform.nextAfter(NEGATIVE_INFINITY, POSITIVE_INFINITY));
    assertExactly(NEGATIVE_INFINITY, Platform.nextAfter(-MAX_VALUE, NEGATIVE_INFINITY));

    // -max
    assertExactly(-MAX_VALUE, Platform.nextAfter(-MAX_VALUE, -MAX_VALUE));
    assertExactly(-1.7976931348623155E308, Platform.nextAfter(-MAX_VALUE, 0));
    assertExactly(-1.7976931348623155E308, Platform.nextAfter(-MAX_VALUE, POSITIVE_INFINITY));

    // -pi
    assertExactly(-3.1415926535897936, Platform.nextAfter(-Math.PI, -MAX_VALUE));
    assertExactly(-Math.PI, Platform.nextAfter(-Math.PI, -Math.PI));
    assertExactly(-3.1415926535897927, Platform.nextAfter(-Math.PI, 0));
    assertExactly(-3.1415926535897927, Platform.nextAfter(-Math.PI, POSITIVE_INFINITY));

    // -min
    assertExactly(-9.9E-324, Platform.nextAfter(-MIN_VALUE, NEGATIVE_INFINITY));
    assertExactly(-MIN_VALUE, Platform.nextAfter(-MIN_VALUE, -MIN_VALUE));
    assertExactly(NEGATIVE_ZERO, Platform.nextAfter(-MIN_VALUE, 0));

    // zero.
    assertExactly(-MIN_VALUE, Platform.nextAfter(0, NEGATIVE_INFINITY));
    assertExactly(-MIN_VALUE, Platform.nextAfter(-0, -MAX_VALUE));
    assertExactly(0.0, Platform.nextAfter(0, 0));
    assertExactly(MIN_VALUE, Platform.nextAfter(0, MIN_VALUE));

    // +min
    assertExactly(0.0, Platform.nextAfter(MIN_VALUE, -MAX_VALUE));
    assertExactly(0.0, Platform.nextAfter(MIN_VALUE, -MIN_VALUE));
    assertExactly(0.0, Platform.nextAfter(MIN_VALUE, 0));
    assertExactly(MIN_VALUE, Platform.nextAfter(MIN_VALUE, MIN_VALUE));
    assertExactly(9.9E-324, Platform.nextAfter(MIN_VALUE, Math.PI));

    // +pi
    assertExactly(3.1415926535897927, Platform.nextAfter(Math.PI, -Math.PI));
    assertExactly(Math.PI, Platform.nextAfter(Math.PI, Math.PI));
    assertExactly(3.1415926535897936, Platform.nextAfter(Math.PI, MAX_VALUE));

    // +max
    assertExactly(1.7976931348623155E308, Platform.nextAfter(MAX_VALUE, 0));
    assertExactly(MAX_VALUE, Platform.nextAfter(MAX_VALUE, MAX_VALUE));
    assertExactly(POSITIVE_INFINITY, Platform.nextAfter(MAX_VALUE, POSITIVE_INFINITY));

    // +inf
    assertExactly(MAX_VALUE, Platform.nextAfter(POSITIVE_INFINITY, 0));
    assertExactly(MAX_VALUE, Platform.nextAfter(POSITIVE_INFINITY, MAX_VALUE));
    assertExactly(POSITIVE_INFINITY, Platform.nextAfter(POSITIVE_INFINITY, POSITIVE_INFINITY));
  }

  public void testNewBigDecimal() {
    BigDecimal ideal = new BigDecimal("0.1000000000000000055511151231257827021181583404541015625");
    BigDecimal actual = Platform.newBigDecimal(0.1);

    // Assert precision is either 25 (GWT) or 55 (NonGWT)
    assertTrue((actual.precision() == 25) || (actual.precision() == 55));

    // Assert rounded value.
    assertEquals(ideal.round(new MathContext(actual.precision())), actual);
  }

  @GwtIncompatible("Javascript uses its own version of Platform#sign")
  public void testSign() {
    // Run some standard determinant checks.
    S2Point a = new S2Point(1.0, 0.0, 0.0);
    S2Point b = new S2Point(0.0, 2.0, 0.0);
    S2Point c = new S2Point(0.0, 0.0, 3.0);
    S2Point d = new S2Point(0.0, 0.0, 0.0);
    assertEquals(1, Platform.sign(a, b, c));
    assertEquals(-1, Platform.sign(a, c, b));
    assertEquals(0, Platform.sign(a, a, c));
    assertEquals(0, Platform.sign(a, b, d));

    // These should result in 0.
    S2Point q = new S2Point(-2.594E-321, 0.0, 1.0);
    S2Point r = new S2Point(1.0, 0.9991685425907498, 0.0);
    S2Point s = new S2Point(0.0, 0.9991685425907498, 1.0);
    int result = Platform.sign(q, r, s);
    assertEquals(0, result);
  }
}
