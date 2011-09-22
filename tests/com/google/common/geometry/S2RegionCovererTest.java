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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

public strictfp class S2RegionCovererTest extends GeometryTestCase {
  private static Logger logger = Logger.getLogger(S2RegionCovererTest.class.getName());

  public void testRandomCells() {
    logger.info("TestRandomCells");

    S2RegionCoverer coverer = new S2RegionCoverer();
    coverer.setMaxCells(1);

    // Test random cell ids at all levels.
    for (int i = 0; i < 10000; ++i) {
      S2CellId id = getRandomCellId();
      S2CellUnion covering = new S2CellUnion();
      coverer.getCovering(new S2Cell(id), covering.cellIds());
      assertEquals(covering.size(), 1);
      assertEquals(covering.cellId(0), id);
    }
  }

  public void checkCovering(
      S2RegionCoverer coverer, S2Region region, ArrayList<S2CellId> covering, boolean interior) {

    // Keep track of how many cells have the same coverer.min_level() ancestor.
    HashMap<S2CellId, Integer> minLevelCells = new HashMap<S2CellId, Integer>();
    for (int i = 0; i < covering.size(); ++i) {
      int level = covering.get(i).level();
      assertTrue(level >= coverer.minLevel());
      assertTrue(level <= coverer.maxLevel());
      assertEquals((level - coverer.minLevel()) % coverer.levelMod(), 0);
      S2CellId key = covering.get(i).parent(coverer.minLevel());
      if (!minLevelCells.containsKey(key)) {
        minLevelCells.put(key, 1);
      } else {
        minLevelCells.put(key, minLevelCells.get(key) + 1);
      }
    }
    if (covering.size() > coverer.maxCells()) {
      // If the covering has more than the requested number of cells, then check
      // that the cell count cannot be reduced by using the parent of some cell.
      for (Integer i : minLevelCells.values()) {
        assertEquals(i.intValue(), 1);
      }
    }

    if (interior) {
      for (int i = 0; i < covering.size(); ++i) {
        assertTrue(region.contains(new S2Cell(covering.get(i))));
      }
    } else {
      S2CellUnion cellUnion = new S2CellUnion();
      cellUnion.initFromCellIds(covering);
      checkCovering(region, cellUnion, true, new S2CellId());
    }
  }

  public void testRandomCaps() {
    logger.info("TestRandomCaps");

    final int kMaxLevel = S2CellId.MAX_LEVEL;
    S2RegionCoverer coverer = new S2RegionCoverer();
    for (int i = 0; i < 1000; ++i) {
      do {
        coverer.setMinLevel(random(kMaxLevel + 1));
        coverer.setMaxLevel(random(kMaxLevel + 1));
      } while (coverer.minLevel() > coverer.maxLevel());
      coverer.setMaxCells(skewed(10));
      coverer.setLevelMod(1 + random(3));
      double maxArea = Math.min(
          4 * S2.M_PI, (3 * coverer.maxCells() + 1) * S2Cell.averageArea(coverer.minLevel()));
      S2Cap cap = getRandomCap(0.1 * S2Cell.averageArea(kMaxLevel), maxArea);
      ArrayList<S2CellId> covering = new ArrayList<S2CellId>();
      ArrayList<S2CellId> interior = new ArrayList<S2CellId>();

      coverer.getCovering(cap, covering);
      checkCovering(coverer, cap, covering, false);

      coverer.getInteriorCovering(cap, interior);
      checkCovering(coverer, cap, interior, true);


      // Check that GetCovering is deterministic.
      ArrayList<S2CellId> covering2 = new ArrayList<S2CellId>();
      coverer.getCovering(cap, covering2);
      assertTrue(covering.equals(covering2));

      // Also check S2CellUnion.denormalize(). The denormalized covering
      // may still be different and smaller than "covering" because
      // S2RegionCoverer does not guarantee that it will not output all four
      // children of the same parent.
      S2CellUnion cells = new S2CellUnion();
      cells.initFromCellIds(covering);
      ArrayList<S2CellId> denormalized = new ArrayList<S2CellId>();
      cells.denormalize(coverer.minLevel(), coverer.levelMod(), denormalized);
      checkCovering(coverer, cap, denormalized, false);
    }
  }

  public void testSimpleCoverings() {
    logger.info("TestSimpleCoverings");

    final int kMaxLevel = S2CellId.MAX_LEVEL;
    S2RegionCoverer coverer = new S2RegionCoverer();
    coverer.setMaxCells(Integer.MAX_VALUE);
    for (int i = 0; i < 1000; ++i) {
      int level = random(kMaxLevel + 1);
      coverer.setMinLevel(level);
      coverer.setMaxLevel(level);
      double maxArea = Math.min(4 * S2.M_PI, 1000 * S2Cell.averageArea(level));
      S2Cap cap = getRandomCap(0.1 * S2Cell.averageArea(kMaxLevel), maxArea);
      ArrayList<S2CellId> covering = new ArrayList<S2CellId>();
      S2RegionCoverer.getSimpleCovering(cap, cap.axis(), level, covering);
      checkCovering(coverer, cap, covering, false);
    }
  }
}
