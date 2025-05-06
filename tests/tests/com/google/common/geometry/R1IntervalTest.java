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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.geometry.R1Interval.Endpoint;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Verifies R1Interval. */
@RunWith(JUnit4.class)
public class R1IntervalTest extends GeometryTestCase {
  /**
   * Test all of the interval operations on the given pair of intervals.
   *
   * @param expected A sequence of "T" and "F" characters corresponding to the expected results of
   *     contains(), interiorContains(), intersects(), and interiorIntersects() respectively.
   */
  private static void testIntervalOps(R1Interval x, R1Interval y, String expected) {
    assertEquals(expected.charAt(0) == 'T', x.contains(y));
    assertEquals(expected.charAt(1) == 'T', x.interiorContains(y));
    assertEquals(expected.charAt(2) == 'T', x.intersects(y));
    assertEquals(expected.charAt(3) == 'T', x.interiorIntersects(y));
    assertEquals(x.contains(y), x.union(y).equals(x));
    assertEquals(x.intersects(y), !x.intersection(y).isEmpty());
  }

  @SuppressWarnings("SelfEquals")
  @Test
  public void testBasics() {
    // Constructors and accessors.
    R1Interval unit = new R1Interval(0, 1);
    R1Interval negunit = new R1Interval(-1, 0);
    assertExactly(0d, unit.lo());
    assertExactly(1d, unit.hi());
    assertExactly(-1d, negunit.getValue(Endpoint.LO));
    assertExactly(-1d, negunit.getValue(Endpoint.LO));
    assertExactly(0d, negunit.getValue(Endpoint.HI));
    R1Interval ten = new R1Interval(0, 0);
    ten.setValue(Endpoint.HI, 10);
    assertEquals(new R1Interval(0, 10), ten);
    ten.setLo(-10);
    assertEquals(new R1Interval(-10, 10), ten);
    ten.setHi(0);
    assertEquals(new R1Interval(-10, 0), ten);
    ten.set(0, 10);
    assertEquals(new R1Interval(0, 10), ten);

    // is_empty()
    R1Interval half = new R1Interval(0.5, 0.5);
    assertFalse(unit.isEmpty());
    assertFalse(half.isEmpty());
    R1Interval empty = R1Interval.empty();
    assertTrue(empty.isEmpty());

    // Equality.
    assertTrue(empty.equals(empty));
    assertTrue(unit.equals(unit));
    assertFalse(unit.equals(empty));
    assertFalse(new R1Interval(1, 2).equals(new R1Interval(1, 3)));

    // Check that the default R1Interval is identical to Empty().
    R1Interval defaultEmpty = new R1Interval();
    assertTrue(defaultEmpty.isEmpty());
    assertExactly(empty.lo(), defaultEmpty.lo());
    assertExactly(empty.hi(), defaultEmpty.hi());

    // getCenter(), getLength()
    assertExactly(0.5, unit.getCenter());
    assertExactly(0.5, half.getCenter());
    assertExactly(1.0, negunit.getLength());
    assertExactly(0.0, half.getLength());
    assertTrue(empty.getLength() < 0);

    // contains(double), interiorContains(double)
    assertTrue(unit.contains(0.5));
    assertTrue(unit.interiorContains(0.5));
    assertTrue(unit.contains(0));
    assertFalse(unit.interiorContains(0));
    assertTrue(unit.contains(1));
    assertFalse(unit.interiorContains(1));

    // contains(R1Interval), interiorContains(R1Interval)
    // intersects(R1Interval), interiorIntersects(R1Interval)
    testIntervalOps(empty, empty, "TTFF");
    testIntervalOps(empty, unit, "FFFF");
    testIntervalOps(unit, half, "TTTT");
    testIntervalOps(unit, unit, "TFTT");
    testIntervalOps(unit, empty, "TTFF");
    testIntervalOps(unit, negunit, "FFTF");
    testIntervalOps(unit, new R1Interval(0, 0.5), "TFTT");
    testIntervalOps(half, new R1Interval(0, 0.5), "FFTF");

    // addPoint()
    R1Interval r = empty;
    r = r.addPoint(5);
    assertExactly(5d, r.lo());
    assertExactly(5d, r.hi());
    r = r.addPoint(-1);
    assertExactly(-1d, r.lo());
    assertExactly(5d, r.hi());
    r = r.addPoint(0);
    assertExactly(-1d, r.lo());
    assertExactly(5d, r.hi());

    // unionInternal()
    r = R1Interval.empty();
    r.unionInternal(5);
    assertExactly(5d, r.lo());
    assertExactly(5d, r.hi());
    r.unionInternal(-1);
    assertExactly(-1d, r.lo());
    assertExactly(5d, r.hi());
    r.unionInternal(0);
    assertExactly(-1d, r.lo());
    assertExactly(5d, r.hi());

    // clampPoint()
    assertExactly(0.3, new R1Interval(0.1, 0.4).clampPoint(0.3));
    assertExactly(0.1, new R1Interval(0.1, 0.4).clampPoint(-7.0));
    assertExactly(0.4, new R1Interval(0.1, 0.4).clampPoint(0.6));

    // fromPointPair()
    assertEquals(new R1Interval(4, 4), R1Interval.fromPointPair(4, 4));
    assertEquals(new R1Interval(-2, -1), R1Interval.fromPointPair(-1, -2));
    assertEquals(new R1Interval(-5, 3), R1Interval.fromPointPair(-5, 3));

    // expanded()
    assertEquals(empty, empty.expanded(0.45));
    assertEquals(new R1Interval(-0.5, 1.5), unit.expanded(0.5));
    assertEquals(new R1Interval(0.5, 0.5), unit.expanded(-0.5));
    assertTrue(unit.expanded(-0.51).isEmpty());
    assertTrue(unit.expanded(-0.51).expanded(0.51).isEmpty());

    // union(), intersection()
    assertEquals(new R1Interval(99, 100), new R1Interval(99, 100).union(empty));
    assertEquals(new R1Interval(99, 100), empty.union(new R1Interval(99, 100)));
    assertTrue(new R1Interval(5, 3).union(new R1Interval(0, -2)).isEmpty());
    assertTrue(new R1Interval(0, -2).union(new R1Interval(5, 3)).isEmpty());
    assertEquals(unit, unit.union(unit));
    assertEquals(new R1Interval(-1, 1), unit.union(negunit));
    assertEquals(new R1Interval(-1, 1), negunit.union(unit));
    assertEquals(unit, half.union(unit));
    assertEquals(half, unit.intersection(half));
    assertEquals(new R1Interval(0, 0), unit.intersection(negunit));
    assertTrue(negunit.intersection(half).isEmpty());
    assertTrue(unit.intersection(empty).isEmpty());
    assertTrue(empty.intersection(unit).isEmpty());
  }

  @Test
  public void testApproxEquals() {
    // Choose two values kLo and kHi such that it's okay to shift an endpoint by kLo (i.e., the
    // resulting interval is equivalent) but not by kHi. The kLo bound is a bit closer to epsilon
    // in Java compared to C++. TODO(torrey): Investigate why, fix if possible.
    final double kLo = 2 * S2.DBL_EPSILON; // < max_error default
    final double kHi = 6 * S2.DBL_EPSILON; // > max_error default

    // Empty intervals.
    R1Interval empty = R1Interval.empty();
    assertTrue(empty.approxEquals(empty));
    assertTrue(new R1Interval(0, 0).approxEquals(empty));
    assertTrue(empty.approxEquals(new R1Interval(0, 0)));
    assertTrue(new R1Interval(1, 1).approxEquals(empty));
    assertTrue(empty.approxEquals(new R1Interval(1, 1)));
    assertFalse(empty.approxEquals(new R1Interval(0, 1)));
    assertTrue(empty.approxEquals(new R1Interval(1, 1 + 2 * kLo)));
    assertFalse(empty.approxEquals(new R1Interval(1, 1 + 2 * kHi)));

    // Singleton intervals.
    assertTrue(new R1Interval(1, 1).approxEquals(new R1Interval(1, 1)));
    assertTrue(new R1Interval(1, 1).approxEquals(new R1Interval(1 - kLo, 1 - kLo)));
    assertTrue(new R1Interval(1, 1).approxEquals(new R1Interval(1 + kLo, 1 + kLo)));
    assertFalse(new R1Interval(1, 1).approxEquals(new R1Interval(1 - kHi, 1)));
    assertFalse(new R1Interval(1, 1).approxEquals(new R1Interval(1, 1 + kHi)));
    assertTrue(new R1Interval(1, 1).approxEquals(new R1Interval(1 - kLo, 1 + kLo)));
    assertFalse(new R1Interval(0, 0).approxEquals(new R1Interval(1, 1)));

    // Other intervals.
    assertTrue(new R1Interval(1 - kLo, 2 + kLo).approxEquals(new R1Interval(1, 2)));
    assertTrue(new R1Interval(1 + kLo, 2 - kLo).approxEquals(new R1Interval(1, 2)));
    assertFalse(new R1Interval(1 - kHi, 2 + kLo).approxEquals(new R1Interval(1, 2)));
    assertFalse(new R1Interval(1 + kHi, 2 - kLo).approxEquals(new R1Interval(1, 2)));
    assertFalse(new R1Interval(1 - kLo, 2 + kHi).approxEquals(new R1Interval(1, 2)));
    assertFalse(new R1Interval(1 + kLo, 2 - kHi).approxEquals(new R1Interval(1, 2)));
  }

  @Test
  public void testOpposites() {
    assertEquals(Endpoint.LO, Endpoint.HI.opposite());
    assertEquals(Endpoint.HI, Endpoint.LO.opposite());
    assertEquals(Endpoint.LO, Endpoint.LO.opposite().opposite());
    assertEquals(Endpoint.HI, Endpoint.HI.opposite().opposite());
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testR1IntervalSerialization() {
    doSerializationTest(new R1Interval(1e5, 1e-6));
  }
}
