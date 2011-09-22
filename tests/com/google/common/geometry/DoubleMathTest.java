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

import com.google.common.geometry.DoubleMath.MantissaExponent;

import junit.framework.TestCase;

import java.util.Random;

/**
 * Tests the advanced floating point operations defined in DoubleMath.java.
 *
 */
public class DoubleMathTest extends TestCase {

  Random rnd = new Random(12345);

  public void testFrexp() {
    for (int i = 0; i < 10; ++i) {
      MantissaExponent me = DoubleMath.frexp(Math.pow(2, i));
      assertEquals(0.5, me.mantissa);
      assertEquals(i + 1, me.exp);
    }

    for (int i = 0; i < 10; ++i) {
      MantissaExponent me = DoubleMath.frexp(-Math.pow(2, i));
      assertEquals(-0.5, me.mantissa);
      assertEquals(i + 1, me.exp);
    }

    MantissaExponent me = DoubleMath.frexp(0);
    assertEquals(0.0, me.mantissa);
    assertEquals(0, me.exp);

    me = DoubleMath.frexp(3);
    assertEquals(0.75, me.mantissa);
    assertEquals(2, me.exp);

    me = DoubleMath.frexp(5);
    assertEquals(0.625, me.mantissa);
    assertEquals(3, me.exp);
  }

  public void testCompareTo() {
    for (int i = 0; i < 100; ++i) {
      double x = rnd.nextDouble() * 100 - 50;
      double y = rnd.nextDouble() * 100 - 50;
      MantissaExponent m1 = DoubleMath.frexp(x);
      MantissaExponent m2 = DoubleMath.frexp(y);

      assertEquals(new Double(x).compareTo(y), m1.compareTo(m2));
    }
  }
}
