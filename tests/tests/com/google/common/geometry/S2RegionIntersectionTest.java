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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Test case for {@link S2RegionIntersection}. */
@RunWith(JUnit4.class)
public class S2RegionIntersectionTest extends GeometryTestCase {
  @Test
  public void testEmptyRegions() throws Exception {
    List<S2Region> regions = new ArrayList<>();
    S2RegionIntersection regionIntersection = new S2RegionIntersection(regions);
    assertTrue(regionIntersection.getCapBound().isFull());
    assertTrue(regionIntersection.getRectBound().isFull());
  }

  private S2RegionIntersection threeRegions() throws Exception {
    S2LatLngRect rect1 =
        new S2LatLngRect(S2LatLng.fromDegrees(-35, -30), S2LatLng.fromDegrees(30, 15));
    S2LatLngRect rect2 =
        new S2LatLngRect(S2LatLng.fromDegrees(-45, -20), S2LatLng.fromDegrees(20, 25));
    S2LatLngRect rect3 =
        new S2LatLngRect(S2LatLng.fromDegrees(-25, -10), S2LatLng.fromDegrees(40, 35));
    List<S2Region> regions = new ArrayList<>();
    regions.add(rect1);
    regions.add(rect2);
    regions.add(rect3);

    return new S2RegionIntersection(regions);
  }

  @Test
  public void testThreeRegions() throws Exception {
    S2RegionIntersection regionIntersection = threeRegions();
    S2LatLngRect expectedRect =
        new S2LatLngRect(S2LatLng.fromDegrees(-25, -10), S2LatLng.fromDegrees(20, 15));
    assertEquals(expectedRect, regionIntersection.getRectBound());

    S2Cap expectedCap = expectedRect.getCapBound();
    assertEquals(expectedCap, regionIntersection.getCapBound());

    // The region intersects but does not contain face 0.
    S2Cell face0 = S2Cell.fromFace(0);
    assertTrue(regionIntersection.mayIntersect(face0));
    assertFalse(regionIntersection.contains(face0));

    // The region intersects and contains the cell of the center point of the intersection region.
    S2Cell center = new S2Cell(expectedRect.getCenter());
    assertTrue(regionIntersection.mayIntersect(center));
    assertTrue(regionIntersection.contains(center));
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2RegionIntersectionSerialization() throws Exception {
    S2RegionIntersection regionIntersection = threeRegions();
    doSerializationTest(regionIntersection);
  }
}
