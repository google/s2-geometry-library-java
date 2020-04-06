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

import static java.lang.Double.MAX_VALUE;
import static java.lang.Double.MIN_VALUE;
import static java.lang.Double.NEGATIVE_INFINITY;
import static java.lang.Double.NaN;
import static java.lang.Double.POSITIVE_INFINITY;

import com.google.common.annotations.GwtCompatible;
import java.math.BigDecimal;
import java.math.MathContext;
import junit.framework.TestCase;

/** Verifies Platform test methods. */
@GwtCompatible
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
    // Expected results from the cross product of each [numerator, denominator] pair defined above.
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
    assertEquals(1.2689709186578246e-116, Platform.ulp(1e-100));
    assertEquals(1.2924697071141057E-26, Platform.ulp(1e-10));
    assertEquals(4.440892098500626e-16, Platform.ulp(Math.PI));
    assertEquals(1.9073486328125e-6, Platform.ulp(1e10));
    assertEquals(1.942668892225729e84, Platform.ulp(1e100));
  }

  public void testNextAfter() {
    // NaNs
    assertTrue(Double.isNaN(Platform.nextAfter(NaN, NaN)));
    assertTrue(Double.isNaN(Platform.nextAfter(NaN, NEGATIVE_INFINITY)));
    assertTrue(Double.isNaN(Platform.nextAfter(NaN, -MAX_VALUE)));
    assertTrue(Double.isNaN(Platform.nextAfter(NEGATIVE_INFINITY, NaN)));

    // -inf
    assertEquals(NEGATIVE_INFINITY, Platform.nextAfter(NEGATIVE_INFINITY, NEGATIVE_INFINITY));
    assertEquals(-MAX_VALUE, Platform.nextAfter(NEGATIVE_INFINITY, -MAX_VALUE));
    assertEquals(-MAX_VALUE, Platform.nextAfter(NEGATIVE_INFINITY, Math.PI));
    assertEquals(-MAX_VALUE, Platform.nextAfter(NEGATIVE_INFINITY, POSITIVE_INFINITY));
    assertEquals(NEGATIVE_INFINITY, Platform.nextAfter(-MAX_VALUE, NEGATIVE_INFINITY));

    // -max
    assertEquals(-MAX_VALUE, Platform.nextAfter(-MAX_VALUE, -MAX_VALUE));
    assertEquals(-1.7976931348623155E308, Platform.nextAfter(-MAX_VALUE, 0));
    assertEquals(-1.7976931348623155E308, Platform.nextAfter(-MAX_VALUE, POSITIVE_INFINITY));

    // -pi
    assertEquals(-3.1415926535897936, Platform.nextAfter(-Math.PI, -MAX_VALUE));
    assertEquals(-Math.PI, Platform.nextAfter(-Math.PI, -Math.PI));
    assertEquals(-3.1415926535897927, Platform.nextAfter(-Math.PI, 0));
    assertEquals(-3.1415926535897927, Platform.nextAfter(-Math.PI, POSITIVE_INFINITY));

    // -min
    assertEquals(-1.0E-323, Platform.nextAfter(-MIN_VALUE, NEGATIVE_INFINITY));
    assertEquals(-MIN_VALUE, Platform.nextAfter(-MIN_VALUE, -MIN_VALUE));
    assertEquals(NEGATIVE_ZERO, Platform.nextAfter(-MIN_VALUE, 0));

    // zero.
    assertEquals(-MIN_VALUE, Platform.nextAfter(0, NEGATIVE_INFINITY));
    assertEquals(-MIN_VALUE, Platform.nextAfter(-0, -MAX_VALUE));
    assertEquals(0.0, Platform.nextAfter(0, 0));
    assertEquals(MIN_VALUE, Platform.nextAfter(0, MIN_VALUE));

    // +min
    assertEquals(0.0, Platform.nextAfter(MIN_VALUE, -MAX_VALUE));
    assertEquals(0.0, Platform.nextAfter(MIN_VALUE, -MIN_VALUE));
    assertEquals(0.0, Platform.nextAfter(MIN_VALUE, 0));
    assertEquals(MIN_VALUE, Platform.nextAfter(MIN_VALUE, MIN_VALUE));
    assertEquals(1.0E-323, Platform.nextAfter(MIN_VALUE, Math.PI));

    // +pi
    assertEquals(3.1415926535897927, Platform.nextAfter(Math.PI, -Math.PI));
    assertEquals(Math.PI, Platform.nextAfter(Math.PI, Math.PI));
    assertEquals(3.1415926535897936, Platform.nextAfter(Math.PI, MAX_VALUE));

    // +max
    assertEquals(1.7976931348623155E308, Platform.nextAfter(MAX_VALUE, 0));
    assertEquals(MAX_VALUE, Platform.nextAfter(MAX_VALUE, MAX_VALUE));
    assertEquals(POSITIVE_INFINITY, Platform.nextAfter(MAX_VALUE, POSITIVE_INFINITY));

    // +inf
    assertEquals(MAX_VALUE, Platform.nextAfter(POSITIVE_INFINITY, 0));
    assertEquals(MAX_VALUE, Platform.nextAfter(POSITIVE_INFINITY, MAX_VALUE));
    assertEquals(POSITIVE_INFINITY, Platform.nextAfter(POSITIVE_INFINITY, POSITIVE_INFINITY));
  }

  public void testNewBigDecimal() {
    BigDecimal ideal = new BigDecimal("0.1000000000000000055511151231257827021181583404541015625");
    BigDecimal actual = Platform.newBigDecimal(0.1);

    // Assert precision is either 25 (GWT) or 55 (NonGWT)
    assertTrue((actual.precision() == 25) || (actual.precision() == 55));

    // Assert rounded value.
    assertEquals(ideal.round(new MathContext(actual.precision())), actual);
  }
}
