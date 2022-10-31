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


/**
 * Tests for methods in GeometryTestCase.
 */
public final class GeometryTestCaseTest extends GeometryTestCase {

  private void assertNotDoubleEquals(double a, double b) {
    assertNotDoubleEquals("", a, b);
  }

  private void assertNotDoubleEquals(String message, double a, double b) {
    try {
      assertDoubleEquals(a, b);
    } catch (AssertionError expected) {
      return; // test passes.
    }
    fail("Expected assertDoubleEquals to fail, but it did not: " + message);
  }

  public void testAssertDoubleEqualsFundamentals() {
    assertDoubleEquals(0, 0);
    // Negative zero is a special case.
    assertDoubleEquals(-0, 0);

    // Tiny numbers with the same sign are considered almost equal.
    assertDoubleEquals(0, nextUp(0));
    assertDoubleEquals(-0.0d, nextDown(0));

    // Tiny non-zero numbers with opposite signs are not considered almost equal, even if they're
    // within MAX_ULPS of each other.
    assertNotDoubleEquals(0, nextDown(0));
    assertNotDoubleEquals(-0.0d, nextUp(0));
    assertNotDoubleEquals(nextUp(0), nextDown(0));

    /* ULP-based comparison for some fixed epsilons.  */

    // 1e-12 is larger than 4 ULPs at 360, so ULP-based checking is stricter than a fixed epsilon of
    // 1e-12 for any reasonable value of degrees or radians.
    assertNotDoubleEquals(360, 360 + 1e-12);
    // But 1e-12 is less than 4 ULPs at 4096.
    assertDoubleEquals(4096, 4096 + 1e-12);

    // 1e-14 is larger than 4 ULPs at 4 * PI, so ULP-based checking is stricter than a fixed
    // epsilon of 1e-14 for any reasonable value of radians.
    assertNotDoubleEquals(4 * PI, 4 * PI + 1e-14);
    // However, for doubles representing a large angle in degrees, differences of 1e-14 could be
    // within 4 ULPs.
    assertDoubleEquals(360, 360 + 1e-14);
  }

  public void testAssertDoubleEqualsEdgeCases() {
    assertNotDoubleEquals(Double.NaN, Double.NaN);
    assertDoubleEquals(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
    assertDoubleEquals(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
    assertNotDoubleEquals(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);

    // Almost Infinities
    double almostPositiveInfinity = nextDown(Double.POSITIVE_INFINITY);
    double almostNegativeInfinity = nextUp(Double.NEGATIVE_INFINITY);
    assertDoubleEquals(almostPositiveInfinity, Double.POSITIVE_INFINITY);
    assertDoubleEquals(almostNegativeInfinity, Double.NEGATIVE_INFINITY);
  }

  public void testAssertDoubleEquals() {
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
        assertDoubleEquals(
            Platform.formatString("Expect doubleEquals(%s, %s), %s up", a, aUp, ups), a, aUp);
        aUp = nextUp(aUp);
        ups++;
      }
      assertNotDoubleEquals(
          Platform.formatString("Expect NOT doubleEquals(%s, %s) %s up", a, aUp, ups), a, aUp);

      double aDown = a;
      int downs = 0;
      for (; downs <= MAX_ULPS; ) {
        assertDoubleEquals(
            Platform.formatString("Expect doubleEquals(%s, %s), %s down", a, aDown, downs),
            a,
            aDown);
        aDown = nextDown(aDown);
        downs++;
      }
      assertNotDoubleEquals(
          Platform.formatString("Expect NOT doubleEquals(%s, %s) %s down", a, aDown, downs),
          a,
          aDown);
    }
  }
}