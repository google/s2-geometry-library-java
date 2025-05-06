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

import static java.lang.Math.PI;
import static org.junit.Assert.fail;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for methods in GeometryTestCase. */
@RunWith(JUnit4.class)
public final class GeometryTestCaseTest extends GeometryTestCase {

  private void assertNotAlmostEquals(double a, double b) {
    assertNotAlmostEquals("", a, b);
  }

  private void assertNotAlmostEquals(String message, double a, double b) {
    try {
      assertAlmostEquals(a, b);
    } catch (AssertionError expected) {
      return; // test passes.
    }
    fail("Expected assertAlmostEquals to fail, but it did not: " + message);
  }

  @Test
  public void testAssertAlmostEqualsFundamentals() {
    assertAlmostEquals(0, 0);
    // Negative zero is a special case.
    assertAlmostEquals(-0.0, 0);

    // Tiny numbers with the same sign are considered almost equal.
    assertAlmostEquals(0, nextUp(0));
    assertAlmostEquals(-0.0d, nextDown(0));

    // Tiny non-zero numbers with opposite signs are not considered almost equal, even if they're
    // within MAX_ULPS of each other.
    assertNotAlmostEquals(0, nextDown(0));
    assertNotAlmostEquals(-0.0d, nextUp(0));
    assertNotAlmostEquals(nextUp(0), nextDown(0));

    /* ULP-based comparison for some fixed epsilons.  */

    // 1e-12 is larger than 4 ULPs at 360, so ULP-based checking is stricter than a fixed epsilon of
    // 1e-12 for any reasonable value of degrees or radians.
    assertNotAlmostEquals(360, 360 + 1e-12);
    // But 1e-12 is less than 4 ULPs at 4096.
    assertAlmostEquals(4096, 4096 + 1e-12);

    // 1e-14 is larger than 4 ULPs at 4 * PI, so ULP-based checking is stricter than a fixed
    // epsilon of 1e-14 for any reasonable value of radians.
    assertNotAlmostEquals(4 * PI, 4 * PI + 1e-14);
    // However, for doubles representing a large angle in degrees, differences of 1e-14 could be
    // within 4 ULPs.
    assertAlmostEquals(360, 360 + 1e-14);
  }

  @Test
  public void testAssertAlmostEqualsEdgeCases() {
    assertNotAlmostEquals(Double.NaN, Double.NaN);
    assertAlmostEquals(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
    assertAlmostEquals(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
    assertNotAlmostEquals(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);

    // Almost Infinities
    double almostPositiveInfinity = nextDown(Double.POSITIVE_INFINITY);
    double almostNegativeInfinity = nextUp(Double.NEGATIVE_INFINITY);
    assertAlmostEquals(almostPositiveInfinity, Double.POSITIVE_INFINITY);
    assertAlmostEquals(almostNegativeInfinity, Double.NEGATIVE_INFINITY);
  }

  @Test
  public void testAssertAlmostEquals() {
    double[] testNumbers = {
          // Small positive numbers
          1.0, 1e-3d, 1e-5d, 1e-7, 1e-12d, 1e-20d, 1e-30d, 1e-40d,
          // Small negative numbers
          -1.0, -1e-3, -1e-5, -1e-7, -1e-12d, -1e-20d, -1e-30d, -1e-40,
          // Large positive numbers
          1e3, 1e5, 1e7, 1e10, 1e20d, 1e30d, 1e60d,
          // Large negative numbers
          -1e3, -1e5, -1e7, -1e10, -1e20d, -1e30d, -1e60d
          };

    for (double a : testNumbers) {
      double aUp = a;
      int ups = 0;
      for (; ups <= MAX_ULPS; ) {
        assertAlmostEquals(
            Platform.formatString("Expect doubleEquals(%s, %s), %s up", a, aUp, ups), a, aUp);
        aUp = nextUp(aUp);
        ups++;
      }
      assertNotAlmostEquals(
          Platform.formatString("Expect NOT doubleEquals(%s, %s) %s up", a, aUp, ups), a, aUp);

      double aDown = a;
      int downs = 0;
      for (; downs <= MAX_ULPS; ) {
        assertAlmostEquals(
            Platform.formatString("Expect doubleEquals(%s, %s), %s down", a, aDown, downs),
            a,
            aDown);
        aDown = nextDown(aDown);
        downs++;
      }
      assertNotAlmostEquals(
          Platform.formatString("Expect NOT doubleEquals(%s, %s) %s down", a, aDown, downs),
          a,
          aDown);
    }
  }
}