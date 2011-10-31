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

public strictfp class R1IntervalTest extends GeometryTestCase {


  /**
   * Test all of the interval operations on the given pair of intervals.
   * "expected_relation" is a sequence of "T" and "F" characters corresponding
   * to the expected results of contains(), interiorContains(), Intersects(),
   * and InteriorIntersects() respectively.
   */
  public void testIntervalOps(R1Interval x, R1Interval y, String expectedRelation) {
    assertEquals(x.contains(y), expectedRelation.charAt(0) == 'T');
    assertEquals(x.interiorContains(y), expectedRelation.charAt(1) == 'T');
    assertEquals(x.intersects(y), expectedRelation.charAt(2) == 'T');
    assertEquals(x.interiorIntersects(y), expectedRelation.charAt(3) == 'T');

    assertEquals(x.contains(y), x.union(y).equals(x));
    assertEquals(x.intersects(y), !x.intersection(y).isEmpty());
  }

  public void testBasic() {
    // Constructors and accessors.
    R1Interval unit = new R1Interval(0, 1);
    R1Interval negunit = new R1Interval(-1, 0);
    assertEquals(unit.lo(), 0.0);
    assertEquals(unit.hi(), 1.0);
    assertEquals(negunit.lo(), -1.0);
    assertEquals(negunit.hi(), 0.0);

    // is_empty()
    R1Interval half = new R1Interval(0.5, 0.5);
    assertTrue(!unit.isEmpty());
    assertTrue(!half.isEmpty());
    R1Interval empty = R1Interval.empty();
    assertTrue(empty.isEmpty());

    // GetCenter(), GetLength()
    assertEquals(unit.getCenter(), 0.5);
    assertEquals(half.getCenter(), 0.5);
    assertEquals(negunit.getLength(), 1.0);
    assertEquals(half.getLength(), 0.0);
    assertTrue(empty.getLength() < 0);

    // contains(double), interiorContains(double)
    assertTrue(unit.contains(0.5));
    assertTrue(unit.interiorContains(0.5));
    assertTrue(unit.contains(0));
    assertTrue(!unit.interiorContains(0));
    assertTrue(unit.contains(1));
    assertTrue(!unit.interiorContains(1));

    // contains(R1Interval), interiorContains(R1Interval)
    // Intersects(R1Interval), InteriorIntersects(R1Interval)
    testIntervalOps(empty, empty, "TTFF");
    testIntervalOps(empty, unit, "FFFF");
    testIntervalOps(unit, half, "TTTT");
    testIntervalOps(unit, unit, "TFTT");
    testIntervalOps(unit, empty, "TTFF");
    testIntervalOps(unit, negunit, "FFTF");
    testIntervalOps(unit, new R1Interval(0, 0.5), "TFTT");
    testIntervalOps(half, new R1Interval(0, 0.5), "FFTF");

    // addPoint()
    R1Interval r;
    r = empty.addPoint(5);
    assertTrue(r.lo() == 5.0 && r.hi() == 5.0);
    r = r.addPoint(-1);
    assertTrue(r.lo() == -1.0 && r.hi() == 5.0);
    r = r.addPoint(0);
    assertTrue(r.lo() == -1.0 && r.hi() == 5.0);

    // fromPointPair()
    assertEquals(R1Interval.fromPointPair(4, 4), new R1Interval(4, 4));
    assertEquals(R1Interval.fromPointPair(-1, -2), new R1Interval(-2, -1));
    assertEquals(R1Interval.fromPointPair(-5, 3), new R1Interval(-5, 3));

    // expanded()
    assertEquals(empty.expanded(0.45), empty);
    assertEquals(unit.expanded(0.5), new R1Interval(-0.5, 1.5));

    // union(), intersection()
    assertTrue(new R1Interval(99, 100).union(empty).equals(new R1Interval(99, 100)));
    assertTrue(empty.union(new R1Interval(99, 100)).equals(new R1Interval(99, 100)));
    assertTrue(new R1Interval(5, 3).union(new R1Interval(0, -2)).isEmpty());
    assertTrue(new R1Interval(0, -2).union(new R1Interval(5, 3)).isEmpty());
    assertTrue(unit.union(unit).equals(unit));
    assertTrue(unit.union(negunit).equals(new R1Interval(-1, 1)));
    assertTrue(negunit.union(unit).equals(new R1Interval(-1, 1)));
    assertTrue(half.union(unit).equals(unit));
    assertTrue(unit.intersection(half).equals(half));
    assertTrue(unit.intersection(negunit).equals(new R1Interval(0, 0)));
    assertTrue(negunit.intersection(half).isEmpty());
    assertTrue(unit.intersection(empty).isEmpty());
    assertTrue(empty.intersection(unit).isEmpty());
  }
}
