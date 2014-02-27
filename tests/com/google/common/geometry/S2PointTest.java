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

import com.google.common.annotations.GwtCompatible;


/**
 * Verifies S2Point.
 */
@GwtCompatible
public class S2PointTest extends GeometryTestCase {
  public void testRotate() {
    // Check simple axial cases.
    S2Point p = makePoint("0:0");
    assertEquals(makePoint("0:90"), p.rotate(S2Point.Z_POS, Math.toRadians(90)), 1e-15);
    assertEquals(makePoint("-90:0"), p.rotate(S2Point.Y_POS, Math.toRadians(90)), 1e-15);
    assertEquals(makePoint("0:0"), p.rotate(S2Point.X_POS, Math.toRadians(90)), 1e-15);
    // Verify rotations in the plane containing the origin and two random points.
    for (int i = 0; i < 1000; i++) {
      S2Point a = randomPoint();
      S2Point b = randomPoint();
      S2Point axis = S2Point.crossProd(a, b);
      double angle = a.angle(b);
      // Rotate 'a' onto 'b'.
      assertEquals(b, a.rotate(axis, angle), 1e-15);
      // Rotate 'b' onto 'a'.
      assertEquals(a, b.rotate(axis, -angle), 1e-15);
      // Rotate 'a' to the midpoint of 'a' and 'b'.
      assertEquals(S2Point.normalize(S2Point.add(a, b)), a.rotate(axis, angle / 2), 1e-14);
      // Rotate 'b' to the antipodal point of 'a'.
      assertEquals(S2Point.neg(a), b.rotate(axis, Math.PI - angle), 1e-14);
    }
  }

  public void testS2Region() {
    S2Point point = new S2Point(1, 0, 0);

    S2Cap expectedCapBound = S2Cap.fromAxisHeight(point, 0);
    assertEquals(expectedCapBound, point.getCapBound());

    S2LatLng ll = new S2LatLng(point);
    S2LatLngRect expectedRect = new S2LatLngRect(ll, ll);
    assertEquals(expectedRect, point.getRectBound());

    // The leaf cell containing a point is still much larger than the point.
    S2Cell cell = new S2Cell(point);
    assertFalse(point.contains(cell));
    assertTrue(point.mayIntersect(cell));
  }
}
