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

import static com.google.common.collect.ImmutableList.toImmutableList;
import static java.lang.Math.PI;
import static java.lang.Math.min;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/** Tests for {@link S2RegionCoverer}. */
public strictfp class S2RegionCovererTest extends GeometryTestCase {

  public void testRandomCells() {
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(1).build();

    // Test random cell ids at all levels.
    for (int i = 0; i < 10000; ++i) {
      S2CellId id = data.getRandomCellId();
      S2CellUnion covering = new S2CellUnion();
      coverer.getCovering(new S2Cell(id), covering.cellIds());
      assertEquals(1, covering.size());
      assertEquals(covering.cellId(0), id);
    }
  }

  public void checkCovering(
      S2RegionCoverer coverer, S2Region region, ArrayList<S2CellId> covering, boolean interior) {

    // Keep track of how many cells have the same coverer.min_level() ancestor.
    HashMap<S2CellId, Integer> minLevelCells = new HashMap<>();
    for (int i = 0; i < covering.size(); ++i) {
      int level = covering.get(i).level();
      assertTrue(level >= coverer.minLevel());
      assertTrue(level <= coverer.maxLevel());
      assertEquals(0, (level - coverer.minLevel()) % coverer.levelMod());
      S2CellId key = covering.get(i).parent(coverer.minLevel());
      if (!minLevelCells.containsKey(key)) {
        minLevelCells.put(key, 1);
      } else {
        minLevelCells.put(key, minLevelCells.get(key) + 1);
      }
    }
    if (covering.size() > coverer.maxCells()) {
      // If the covering has more than the requested number of cells, then check that the cell count
      // cannot be reduced by using the parent of some cell.
      for (Integer i : minLevelCells.values()) {
        assertEquals(1, i.intValue());
      }
    }

    if (interior) {
      for (int i = 0; i < covering.size(); ++i) {
        assertTrue(region.contains(new S2Cell(covering.get(i))));
      }
    } else {
      S2CellUnion cellUnion = new S2CellUnion();
      cellUnion.initFromCellIds(Lists.newArrayList(covering));
      checkCovering(region, cellUnion, true, new S2CellId());
    }
  }

  public void testRandomCaps() {
    final int kMaxLevel = S2CellId.MAX_LEVEL;
    for (int i = 0; i < 1000; ++i) {
      S2RegionCoverer.Builder covererBuilder =
          S2RegionCoverer.builder().setMaxCells(data.skewed(10)).setLevelMod(1 + data.random(3));
      do {
        covererBuilder.setMinLevel(data.random(kMaxLevel + 1));
        covererBuilder.setMaxLevel(data.random(kMaxLevel + 1));
      } while (covererBuilder.getMinLevel() > covererBuilder.getMaxLevel());
      S2RegionCoverer coverer = covererBuilder.build();
      double maxArea = min(
          4 * PI, (3 * coverer.maxCells() + 1) * S2Cell.averageArea(coverer.minLevel()));
      S2Cap cap = data.getRandomCap(0.1 * S2Cell.averageArea(kMaxLevel), maxArea);
      ArrayList<S2CellId> covering = new ArrayList<>();
      ArrayList<S2CellId> interior = new ArrayList<>();

      coverer.getCovering(cap, covering);
      checkCovering(coverer, cap, covering, false);

      coverer.getInteriorCovering(cap, interior);
      checkCovering(coverer, cap, interior, true);

      // Check that GetCovering is deterministic.
      ArrayList<S2CellId> covering2 = new ArrayList<>();
      coverer.getCovering(cap, covering2);
      assertEquals(covering, covering2);

      // Also check S2CellUnion.denormalize(). The denormalized covering may still be different and
      // smaller than "covering" because S2RegionCoverer does not guarantee that it will not output
      // all four children of the same parent.
      S2CellUnion cells = new S2CellUnion();
      cells.initFromCellIds(covering);
      ArrayList<S2CellId> denormalized = new ArrayList<>();
      cells.denormalize(coverer.minLevel(), coverer.levelMod(), denormalized);
      checkCovering(coverer, cap, denormalized, false);
    }
  }

  public void testSimpleCoverings() {
    final int kMaxLevel = S2CellId.MAX_LEVEL;
    for (int i = 0; i < 1000; ++i) {
      int level = data.random(kMaxLevel + 1);
      S2RegionCoverer coverer =
          S2RegionCoverer.builder()
              .setMaxCells(Integer.MAX_VALUE)
              .setMinLevel(level)
              .setMaxLevel(level)
              .build();
      double maxArea = min(4 * PI, 1000 * S2Cell.averageArea(level));
      S2Cap cap = data.getRandomCap(0.1 * S2Cell.averageArea(kMaxLevel), maxArea);
      ArrayList<S2CellId> covering = new ArrayList<>();
      S2RegionCoverer.getSimpleCovering(cap, cap.axis(), level, covering);
      checkCovering(coverer, cap, covering, false);
    }
  }

  public void testPolylineCovering() {
    S2RegionCoverer coverer = S2RegionCoverer.DEFAULT;
    S2Polyline line =
        new S2Polyline(
            Lists.newArrayList(
                S2LatLng.fromDegrees(0, 0).toPoint(), S2LatLng.fromDegrees(0, 1).toPoint()));
    List<String> tokens = Lists.newArrayList();
    for (S2CellId id : coverer.getCovering(line).cellIds()) {
      tokens.add(id.toToken());
    }
    assertEquals(
        tokens,
        Lists.newArrayList(
            "0555555555555555",
            "0fffffffffffffff",
            "10001",
            "10007",
            "10009",
            "1000a4",
            "1aaa09e01fffdf7f",
            "1aaaaaaaaaaaaaab"));
  }

  public void testPolylineCoveringJavaCcConsistency() {
    // This test ensures S2RegionCoverer implementation is consistent between Java and C++.
    S2Polyline polyline =
        new S2Polyline(
            ImmutableList.of(
                S2LatLng.fromDegrees(-33.8663457, 151.1960891).toPoint(),
                S2LatLng.fromDegrees(-33.866094000000004, 151.19517439999998).toPoint()));
    S2RegionCoverer coverer =
        S2RegionCoverer.builder()
            .setMinLevel(0)
            .setMaxLevel(22)
            .setMaxCells(Integer.MAX_VALUE)
            .build();
    ArrayList<S2CellId> ret = new ArrayList<>();
    coverer.getCovering(polyline, ret);
    List<String> expected =
        Arrays.asList(
            "6b12ae36313d",
            "6b12ae36313f",
            "6b12ae363141",
            "6b12ae363143",
            "6b12ae363145",
            "6b12ae363159",
            "6b12ae36315b",
            "6b12ae363343",
            "6b12ae363345",
            "6b12ae36334d",
            "6b12ae36334f",
            "6b12ae363369",
            "6b12ae36336f",
            "6b12ae363371",
            "6b12ae363377",
            "6b12ae363391",
            "6b12ae363393",
            "6b12ae36339b",
            "6b12ae36339d",
            "6b12ae3633e3",
            "6b12ae3633e5",
            "6b12ae3633ed",
            "6b12ae3633ef",
            "6b12ae37cc11",
            "6b12ae37cc13",
            "6b12ae37cc1b",
            "6b12ae37cc1d",
            "6b12ae37cc63",
            "6b12ae37cc65",
            "6b12ae37cc6d",
            "6b12ae37cc6f",
            "6b12ae37cc89",
            "6b12ae37cc8f",
            "6b12ae37cc91",
            "6b12ae37cc97",
            "6b12ae37ccb1",
            "6b12ae37ccb3",
            "6b12ae37ccbb",
            "6b12ae37ccbd",
            "6b12ae37cea5",
            "6b12ae37cea7",
            "6b12ae37cebb");
    assertEquals(expected, ret.stream().map(S2CellId::toToken).collect(toImmutableList()));
  }

  public void testInteriorCovering() {
    // We construct the region the following way. Start with S2 cell of level l. Remove from it
    // one of its grandchildren (level l+2). If we then set:
    //   minLevel = l
    //   maxLevel = l + 3
    //   maxCells = 3
    // the best interior covering should contain 3 children of the initial cell, that were not
    // affected by removal of a grandchild.
    final int level = 12;

    S2CellId smallCell = S2CellId.fromPoint(data.getRandomPoint()).parent(level + 2);
    S2CellId largeCell = smallCell.parent(level);

    S2CellUnion smallCellUnion = new S2CellUnion();
    smallCellUnion.initFromCellIds(Lists.newArrayList(smallCell));

    S2CellUnion largeCellUnion = new S2CellUnion();
    largeCellUnion.initFromCellIds(Lists.newArrayList(largeCell));

    // Because the Java S2CellUnion doesn't have getDifference(), we construct it manually by
    // taking all grandchildren except the missing one.
    ArrayList<S2CellId> diffList = Lists.newArrayList();
    for (S2CellId id = largeCell.childBegin(level + 2);
        !id.equals(largeCell.childEnd(level + 2));
        id = id.next()) {
      if (!id.equals(smallCell)) {
        diffList.add(id);
      }
    }

    S2CellUnion diff = new S2CellUnion();
    diff.initFromCellIds(diffList);

    S2RegionCoverer.Builder covererBuilder =
        S2RegionCoverer.builder().setMaxCells(3).setMaxLevel(level + 3).setMinLevel(level);

    S2RegionCoverer coverer = covererBuilder.setMaxCells(3).build();
    ArrayList<S2CellId> interior = Lists.newArrayList();
    coverer.getInteriorCovering(diff, interior);
    assertEquals(3, interior.size());

    for (int i = 0; i < 3; ++i) {
      assertEquals(level + 1, interior.get(i).level());
    }

    // For a full cover, we'd need 6 cells (3 at level+1, and 3 at level+2). If we allow 5, ensure
    // that all 3 large cells are there and we get 2 small ones.
    coverer = covererBuilder.setMaxCells(5).build();
    ArrayList<S2CellId> moreInterior = Lists.newArrayList();
    coverer.getInteriorCovering(diff, moreInterior);

    assertEquals(5, moreInterior.size());
    moreInterior.removeAll(interior);
    assertEquals(2, moreInterior.size());
    assertEquals(level + 2, moreInterior.get(0).level());
    assertEquals(level + 2, moreInterior.get(1).level());
  }
}
