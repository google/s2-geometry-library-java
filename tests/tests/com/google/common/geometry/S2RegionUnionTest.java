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
import com.google.common.collect.ImmutableList;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Test case for {@link S2RegionUnion}. */
@RunWith(JUnit4.class)
public class S2RegionUnionTest extends GeometryTestCase {
  @Test
  public void testEmptyRegions() throws Exception {
    List<S2Region> regions = new ArrayList<>();
    S2RegionUnion regionUnion = new S2RegionUnion(regions);
    assertTrue(regionUnion.getCapBound().isEmpty());
    assertTrue(regionUnion.getRectBound().isEmpty());
  }

  @Test
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
    assertTrue(
        "Expected: " + expectedCap + "; actual: " + regionUnion.getCapBound(),
        regionUnion.getCapBound().approxEquals(expectedCap, 1e-3));

    S2LatLngRect expectedRect =
        new S2LatLngRect(S2LatLng.fromDegrees(-35, -40), S2LatLng.fromDegrees(35, 40));
    assertEquals(expectedRect, regionUnion.getRectBound());

    S2Cell face0 = S2Cell.fromFace(0);
    assertTrue(regionUnion.mayIntersect(face0));
    assertFalse(regionUnion.contains(face0));

    // The region intersects, but does not contain, the cells of the original points.
    assertTrue(regionUnion.mayIntersect(new S2Cell(latLng1)));
    assertTrue(regionUnion.mayIntersect(new S2Cell(latLng2)));
    assertTrue(regionUnion.mayIntersect(new S2Cell(latLng3)));
    assertFalse(regionUnion.contains(new S2Cell(latLng1)));
    assertFalse(regionUnion.contains(new S2Cell(latLng2)));
    assertFalse(regionUnion.contains(new S2Cell(latLng3)));
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2RegionUnionSerialization() {
    ImmutableList<S2Point> vertices =
        ImmutableList.of(
            new S2Point(0, 1, 2).normalize(),
            new S2Point(.3, .4, 0).normalize(),
            new S2Point(5, 6, 7).normalize());
    S2Loop loop = new S2Loop(vertices);
    S2Polyline polyline = new S2Polyline(vertices);
    S2Cap cap = S2Cap.fromAxisHeight(new S2Point(2, 2, 2).normalize(), 0);
    S2Cell cell = S2Cell.fromFace(0);

    List<S2Loop> loops = new ArrayList<>();
    loops.add(loop);
    S2Polygon polygon = new S2Polygon(loops);

    ArrayList<Long> cellIds = new ArrayList<>();
    cellIds.add(S2CellId.fromFacePosLevel(1, 0, 1).id());
    cellIds.add(S2CellId.fromFacePosLevel(5, 0, 30).id());
    S2CellUnion cellUnion = new S2CellUnion();
    cellUnion.initFromIds(cellIds);

    S2Point point = new S2Point(3, 1, 0);
    S2LatLngRect rect = S2LatLngRect.full();
    List<S2Region> subRegions = new ArrayList<>();
    subRegions.add(point);
    subRegions.add(rect);
    S2RegionUnion regionUnion = new S2RegionUnion(subRegions);

    List<S2Region> regions = new ArrayList<>();
    regions.add(loop);
    regions.add(polyline);
    regions.add(regionUnion);
    regions.add(cap);
    regions.add(cell);
    regions.add(cellUnion);
    regions.add(point);
    regions.add(polygon);
    doSerializationTest(new S2RegionUnion(regions));
  }
}
