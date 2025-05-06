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

import static org.junit.Assert.assertEquals;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class ProjectionTest extends GeometryTestCase {

  @Test
  public void testInterpolateArgumentsAreNotReversed() {
    assertEquals(
        new R2Vector(1.5, 6),
        Projection.interpolate(0.25, new R2Vector(1.0, 5.0), new R2Vector(3.0, 9.0)));
  }

  @Test
  public void testInterpolateExtrapolation() {
    assertEquals(
        new R2Vector(-3.0, 0),
        Projection.interpolate(-2, new R2Vector(1.0, 0.0), new R2Vector(3.0, 0.0)));
  }

  @Test
  public void testWrapDestination() {
    // Prefer traversing the Antemeridian rather than traversing the Prime Meridian.
    Projection proj = new Projection.PlateCarreeProjection(180);
    R2Vector easternPoint = new R2Vector(170, 0);
    R2Vector westernPoint = new R2Vector(-170, 0);
    R2Vector adjustedWesternPoint = proj.wrapDestination(easternPoint, westernPoint);
    assertEquals(new R2Vector(190, 0), adjustedWesternPoint);
  }

  @Test
  public void testInterpolateCheckSameLengthAtBothEndpoints() {
    // Check that interpolation is exact at both endpoints.
    R2Vector a = new R2Vector(1.234, -5.456e-20);
    R2Vector b = new R2Vector(2.1234e-20, 7.456);
    assertEquals(a, Projection.interpolate(0, a, b));
    assertEquals(b, Projection.interpolate(1, a, b));
  }

  @Test
  public void testMercatorUnproject() {
    Projection.MercatorProjection proj = new Projection.MercatorProjection(180);
    double inf = Double.POSITIVE_INFINITY;
    assertProjectUnproject(proj, new R2Vector(0, 0), new S2Point(1, 0, 0));
    assertProjectUnproject(proj, new R2Vector(180, 0), new S2Point(-1, 0, 0));
    assertProjectUnproject(proj, new R2Vector(90, 0), new S2Point(0, 1, 0));
    assertProjectUnproject(proj, new R2Vector(-90, 0), new S2Point(0, -1, 0));
    assertProjectUnproject(proj, new R2Vector(0, inf), new S2Point(0, 0, 1));
    assertProjectUnproject(proj, new R2Vector(0, -inf), new S2Point(0, 0, -1));
    assertProjectUnproject(
        proj, new R2Vector(0, 70.25557896783025), S2LatLng.fromRadians(1, 0).toPoint());
  }

  @Test
  public void testPlateCarreeUnproject() {
    Projection.PlateCarreeProjection proj = new Projection.PlateCarreeProjection(180);
    assertProjectUnproject(proj, new R2Vector(0, 0), new S2Point(1, 0, 0));
    assertProjectUnproject(proj, new R2Vector(180, 0), new S2Point(-1, 0, 0));
    assertProjectUnproject(proj, new R2Vector(90, 0), new S2Point(0, 1, 0));
    assertProjectUnproject(proj, new R2Vector(-90, 0), new S2Point(0, -1, 0));
    assertProjectUnproject(proj, new R2Vector(0, 90), new S2Point(0, 0, 1));
    assertProjectUnproject(proj, new R2Vector(0, -90), new S2Point(0, 0, -1));
  }

  /**
   * Asserts that roundtripping from projecting and unprojecting an S2 point yields the same value.
   */
  private static final void assertProjectUnproject(Projection projection, R2Vector px, S2Point x) {
    // The arguments are chosen such that projection is exact, but
    // unprojection may not be.
    assertEquals(px, projection.project(x));
    assertPointsWithinDistance(x, projection.unproject(px), 1e-15);
  }
}
