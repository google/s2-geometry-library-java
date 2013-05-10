/*
 * Copyright 2013 Google Inc.
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

import com.google.common.geometry.R1Interval.Endpoint;
import com.google.common.geometry.R2Rect.Axis;


/**
 * Verifies R2Rect. Most of the R2Rect methods have trivial implementations in terms of the
 * R1Interval class, so most of the testing is done in that unit test.
 *
 * <p>Package private because the underlying class is.
 */
public class R2RectTest extends GeometryTestCase {
  /**
   * Tests all of the interval operations on the given pair of intervals.
   *
   * @param x The first rectangle to test.
   * @param y The second rectangle to test.
   * @param expectedRexion a sequence of "T" and "F" characters corresponding to the expected
   *        results of contains, interior contains, intersects, and interiorIntersects predicates.
   * @param expectedUnion The expected union result.
   * @param expectedIntersection The expected intersection result.
   */
  private static void testIntervalOps(R2Rect x, R2Rect y, String expectedRexion,
      R2Rect expectedUnion, R2Rect expectedIntersection) {
    assertEquals(expectedRexion.charAt(0) == 'T', x.contains(y));
    assertEquals(expectedRexion.charAt(1) == 'T', x.interiorContains(y));
    assertEquals(expectedRexion.charAt(2) == 'T', x.intersects(y));
    assertEquals(expectedRexion.charAt(3) == 'T', x.interiorIntersects(y));

    assertEquals(x.union(y).equals(x), x.contains(y));
    assertEquals(!x.intersection(y).isEmpty(), x.intersects(y));

    assertEquals(expectedUnion, x.union(y));
    assertEquals(expectedIntersection, x.intersection(y));

    if (y.getSize().equals(new R2Vector(0, 0))) {
      R2Rect r = x;
      r.addPoint(y.lo());
      assertEquals(expectedUnion, r);
    }
  }

  public void testEmptyRectangles() {
    R2Rect empty = R2Rect.empty();
    assertTrue(empty.isValid());
    assertTrue(empty.isEmpty());
  }

  public void testConstructorsAndAccessors() {
    R2Rect r = new R2Rect(new R2Vector(0.1, 0), new R2Vector(0.25, 1));
    assertEquals(0.1, r.x().lo());
    assertEquals(0.25, r.x().hi());
    assertEquals(0.0, r.y().lo());
    assertEquals(1.0, r.y().hi());

    assertEquals(0.1, r.getInterval(Axis.X).getValue(Endpoint.LO));
    assertEquals(0.25, r.getInterval(Axis.X).getValue(Endpoint.HI));
    assertEquals(0.0, r.getInterval(Axis.Y).getValue(Endpoint.LO));
    assertEquals(1.0, r.getInterval(Axis.Y).getValue(Endpoint.HI));

    assertEquals(new R1Interval(0.1, 0.25), r.x());
    assertEquals(new R1Interval(0, 1), r.y());

    assertEquals(new R1Interval(0.1, 0.25), r.getInterval(Axis.X));
    assertEquals(new R1Interval(0, 1), r.getInterval(Axis.Y));

    r.getInterval(Axis.X).set(3, 4);
    assertEquals(new R1Interval(3, 4), r.getInterval(Axis.X));
    r.getInterval(Axis.Y).setValue(Endpoint.LO, 5);
    r.getInterval(Axis.Y).setValue(Endpoint.HI, 6);
    assertEquals(new R1Interval(5, 6), r.getInterval(Axis.Y));

    R2Rect r2 = new R2Rect();
    assertTrue(r2.isEmpty());
  }

  public void testFromCenterSize() {
    assertTrue(R2Rect.fromCenterSize(new R2Vector(0.3, 0.5), new R2Vector(0.2, 0.4)).
        approxEquals(new R2Rect(new R2Vector(0.2, 0.3), new R2Vector(0.4, 0.7))));
    assertTrue(R2Rect.fromCenterSize(new R2Vector(1, 0.1), new R2Vector(0, 2)).
        approxEquals(new R2Rect(new R2Vector(1, -0.9), new R2Vector(1, 1.1))));
  }

  public void testFromPoint() {
    R2Rect d1 = new R2Rect(new R2Vector(0.1, 0), new R2Vector(0.25, 1));
    assertEquals(new R2Rect(d1.lo(), d1.lo()), R2Rect.fromPoint(d1.lo()));
    assertEquals(new R2Rect(new R2Vector(0.15, 0.3), new R2Vector(0.35, 0.9)),
        R2Rect.fromPointPair(new R2Vector(0.15, 0.9), new R2Vector(0.35, 0.3)));
    assertEquals(new R2Rect(new R2Vector(0.12, 0), new R2Vector(0.83, 0.5)),
        R2Rect.fromPointPair(new R2Vector(0.83, 0), new R2Vector(0.12, 0.5)));
  }

  public void testSimplePredicates() {
    R2Vector sw1 = new R2Vector(0, 0.25);
    R2Vector ne1 = new R2Vector(0.5, 0.75);
    R2Rect r1 = new R2Rect(sw1, ne1);

    assertEquals(new R2Vector(0.25, 0.5), r1.getCenter());
    assertEquals(new R2Vector(0, 0.25), r1.getVertex(0));
    assertEquals(new R2Vector(0.5, 0.25), r1.getVertex(1));
    assertEquals(new R2Vector(0.5, 0.75), r1.getVertex(2));
    assertEquals(new R2Vector(0, 0.75), r1.getVertex(3));
    assertTrue(r1.contains(new R2Vector(0.2, 0.4)));
    assertFalse(r1.contains(new R2Vector(0.2, 0.8)));
    assertFalse(r1.contains(new R2Vector(-0.1, 0.4)));
    assertFalse(r1.contains(new R2Vector(0.6, 0.1)));
    assertTrue(r1.contains(sw1));
    assertTrue(r1.contains(ne1));
    assertFalse(r1.interiorContains(sw1));
    assertFalse(r1.interiorContains(ne1));

    // Make sure that GetVertex() returns vertices in CCW order.
    for (int k = 0; k < 4; ++k) {
      R2Vector a = r1.getVertex((k - 1) & 3);
      R2Vector b = r1.getVertex(k);
      R2Vector c = r1.getVertex((k + 1) & 3);
      assertTrue(R2Vector.sub(b, a).ortho().dotProd(R2Vector.sub(c, a)) > 0);
    }
  }

  public void testIntervalOperations() {
    R2Rect empty = R2Rect.empty();
    R2Vector sw1 = new R2Vector(0, 0.25);
    R2Vector ne1 = new R2Vector(0.5, 0.75);
    R2Rect r1 = new R2Rect(sw1, ne1);
    R2Rect r1Mid = new R2Rect(new R2Vector(0.25, 0.5), new R2Vector(0.25, 0.5));
    R2Rect rSw1 = new R2Rect(sw1, sw1);
    R2Rect rNe1 = new R2Rect(ne1, ne1);

    testIntervalOps(r1, r1Mid, "TTTT", r1, r1Mid);
    testIntervalOps(r1, rSw1, "TFTF", r1, rSw1);
    testIntervalOps(r1, rNe1, "TFTF", r1, rNe1);

    assertEquals(new R2Rect(new R2Vector(0, 0.25), new R2Vector(0.5, 0.75)), r1);
    testIntervalOps(r1, new R2Rect(new R2Vector(0.45, 0.1), new R2Vector(0.75, 0.3)), "FFTT",
        new R2Rect(new R2Vector(0, 0.1), new R2Vector(0.75, 0.75)),
        new R2Rect(new R2Vector(0.45, 0.25), new R2Vector(0.5, 0.3)));
    testIntervalOps(r1, new R2Rect(new R2Vector(0.5, 0.1), new R2Vector(0.7, 0.3)), "FFTF",
        new R2Rect(new R2Vector(0, 0.1), new R2Vector(0.7, 0.75)),
        new R2Rect(new R2Vector(0.5, 0.25), new R2Vector(0.5, 0.3)));
    testIntervalOps(r1, new R2Rect(new R2Vector(0.45, 0.1), new R2Vector(0.7, 0.25)), "FFTF",
        new R2Rect(new R2Vector(0, 0.1), new R2Vector(0.7, 0.75)),
        new R2Rect(new R2Vector(0.45, 0.25), new R2Vector(0.5, 0.25)));

    testIntervalOps(new R2Rect(new R2Vector(0.1, 0.2), new R2Vector(0.1, 0.3)),
        new R2Rect(new R2Vector(0.15, 0.7), new R2Vector(0.2, 0.8)), "FFFF",
        new R2Rect(new R2Vector(0.1, 0.2), new R2Vector(0.2, 0.8)),
        empty);

    // Check that the intersection of two rectangles that overlap in x but not y
    // is valid, and vice versa.
    testIntervalOps(new R2Rect(new R2Vector(0.1, 0.2), new R2Vector(0.4, 0.5)),
        new R2Rect(new R2Vector(0, 0), new R2Vector(0.2, 0.1)), "FFFF",
        new R2Rect(new R2Vector(0, 0), new R2Vector(0.4, 0.5)), empty);
    testIntervalOps(new R2Rect(new R2Vector(0, 0), new R2Vector(0.1, 0.3)),
        new R2Rect(new R2Vector(0.2, 0.1), new R2Vector(0.3, 0.4)), "FFFF",
        new R2Rect(new R2Vector(0, 0), new R2Vector(0.3, 0.4)), empty);
  }

  public void testAddPoint() {
    R2Vector sw1 = new R2Vector(0, 0.25);
    R2Vector ne1 = new R2Vector(0.5, 0.75);
    R2Rect r1 = new R2Rect(sw1, ne1);
    R2Rect r2 = R2Rect.empty();
    r2.addPoint(new R2Vector(0, 0.25));
    r2.addPoint(new R2Vector(0.5, 0.25));
    r2.addPoint(new R2Vector(0, 0.75));
    r2.addPoint(new R2Vector(0.1, 0.4));
    assertEquals(r1, r2);
  }

  public void testClampPoint() {
    R2Rect r1 = new R2Rect(new R1Interval(0, 0.5), new R1Interval(0.25, 0.75));
    assertEquals(new R2Vector(0, 0.25), r1.clampPoint(new R2Vector(-0.01, 0.24)));
    assertEquals(new R2Vector(0, 0.48), r1.clampPoint(new R2Vector(-5.0, 0.48)));
    assertEquals(new R2Vector(0, 0.75), r1.clampPoint(new R2Vector(-5.0, 2.48)));
    assertEquals(new R2Vector(0.19, 0.75), r1.clampPoint(new R2Vector(0.19, 2.48)));
    assertEquals(new R2Vector(0.5, 0.75), r1.clampPoint(new R2Vector(6.19, 2.48)));
    assertEquals(new R2Vector(0.5, 0.53), r1.clampPoint(new R2Vector(6.19, 0.53)));
    assertEquals(new R2Vector(0.5, 0.25), r1.clampPoint(new R2Vector(6.19, -2.53)));
    assertEquals(new R2Vector(0.33, 0.25), r1.clampPoint(new R2Vector(0.33, -2.53)));
    assertEquals(new R2Vector(0.33, 0.37), r1.clampPoint(new R2Vector(0.33, 0.37)));
  }

  public void testexpanded() {
    assertTrue(R2Rect.empty().expanded(new R2Vector(0.1, 0.3)).isEmpty());
    assertTrue(R2Rect.empty().expanded(new R2Vector(-0.1, -0.3)).isEmpty());
    assertTrue(new R2Rect(new R2Vector(0.2, 0.4), new R2Vector(0.3, 0.7))
        .expanded(new R2Vector(0.1, 0.3))
        .approxEquals(new R2Rect(new R2Vector(0.1, 0.1), new R2Vector(0.4, 1.0))));
    assertTrue(new R2Rect(new R2Vector(0.2, 0.4), new R2Vector(0.3, 0.7))
        .expanded(new R2Vector(-0.1, 0.3)).isEmpty());
    assertTrue(new R2Rect(new R2Vector(0.2, 0.4), new R2Vector(0.3, 0.7))
        .expanded(new R2Vector(0.1, -0.2)).isEmpty());
    assertTrue(new R2Rect(new R2Vector(0.2, 0.4), new R2Vector(0.3, 0.7))
        .expanded(new R2Vector(0.1, -0.1))
        .approxEquals(new R2Rect(new R2Vector(0.1, 0.5), new R2Vector(0.4, 0.6))));
    assertTrue(new R2Rect(new R2Vector(0.2, 0.4), new R2Vector(0.3, 0.7)).expanded(0.1)
        .approxEquals(new R2Rect(new R2Vector(0.1, 0.3), new R2Vector(0.4, 0.8))));
  }
}