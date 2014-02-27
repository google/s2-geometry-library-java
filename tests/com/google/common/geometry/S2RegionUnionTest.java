/*
 * Copyright 2014 Google Inc.
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

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;

import com.google.common.annotations.GwtCompatible;

import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Test case for {@link S2RegionUnion}.
 */
@GwtCompatible
public class S2RegionUnionTest extends TestCase {
  public void testEmptyRegions() throws Exception {
    List<S2Region> regions = new ArrayList<>();
    S2RegionUnion regionUnion = new S2RegionUnion(regions);
    assertTrue(regionUnion.getCapBound().isEmpty());
    assertTrue(regionUnion.getRectBound().isEmpty());
  }

  public void testThreeRegions() throws Exception {
    S2LatLng latLng1 = S2LatLng.fromDegrees(35, -40);
    S2LatLng latLng2 = S2LatLng.fromDegrees(-35, 40);
    S2LatLng latLng3 = S2LatLng.fromDegrees(10, 10);
    List<S2Region> regions = new ArrayList<>();
    regions.add(latLng1.toPoint());
    regions.add(latLng2.toPoint());
    regions.add(latLng3.toPoint());

    S2RegionUnion regionUnion = new S2RegionUnion(regions);

    S2Cap expectedCap = S2Cap.empty();
    expectedCap = expectedCap.addPoint(latLng1.toPoint());
    expectedCap = expectedCap.addPoint(latLng2.toPoint());
    expectedCap = expectedCap.addPoint(latLng3.toPoint());
    assertTrue("Expected: " + expectedCap
        + "; actual: " + regionUnion.getCapBound(),
        regionUnion.getCapBound().approxEquals(expectedCap, 1e-3));

    S2LatLngRect expectedRect = new S2LatLngRect(
        S2LatLng.fromDegrees(-35, -40), S2LatLng.fromDegrees(35, 40));
    assertEquals(expectedRect, regionUnion.getRectBound());

    S2Cell face0 = S2Cell.fromFace(0);
    assertTrue(regionUnion.mayIntersect(face0));
    assertFalse(regionUnion.contains(face0));

    // The region intersects, but does not contain, the cells of the original
    // points.
    assertTrue(regionUnion.mayIntersect(new S2Cell(latLng1)));
    assertTrue(regionUnion.mayIntersect(new S2Cell(latLng2)));
    assertTrue(regionUnion.mayIntersect(new S2Cell(latLng3)));
    assertFalse(regionUnion.contains(new S2Cell(latLng1)));
    assertFalse(regionUnion.contains(new S2Cell(latLng2)));
    assertFalse(regionUnion.contains(new S2Cell(latLng3)));
  }
}
