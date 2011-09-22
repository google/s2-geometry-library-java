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

import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

public strictfp class S2CellUnionTest extends GeometryTestCase {
  public static Logger logger = Logger.getLogger(S2CellUnionTest.class.getName());

  public void testBasic() {
    logger.info("TestBasic");

    S2CellUnion empty = new S2CellUnion();
    ArrayList<S2CellId> ids = Lists.newArrayList();
    empty.initFromCellIds(ids);
    assertEquals(0, empty.size());

    S2CellId face1Id = S2CellId.fromFacePosLevel(1, 0, 0);
    S2CellUnion face1Union = new S2CellUnion();
    ids.add(face1Id);
    face1Union.initFromCellIds(ids);
    assertEquals(1, face1Union.size());
    assertEquals(face1Id, face1Union.cellId(0));

    S2CellId face2Id = S2CellId.fromFacePosLevel(2, 0, 0);
    S2CellUnion face2Union = new S2CellUnion();
    ArrayList<Long> cellids = Lists.newArrayList();
    cellids.add(face2Id.id());
    face2Union.initFromIds(cellids);
    assertEquals(1, face2Union.size());
    assertEquals(face2Id, face2Union.cellId(0));

    S2Cell face1Cell = new S2Cell(face1Id);
    S2Cell face2Cell = new S2Cell(face2Id);
    assertTrue(face1Union.contains(face1Cell));
    assertTrue(!face1Union.contains(face2Cell));
  }

  public void testContainsCellUnion() {
    logger.info("TestContainsCellUnion");

    Set<S2CellId> randomCells = new HashSet<S2CellId>();
    for (int i = 0; i < 100; i++) {
      randomCells.add(getRandomCellId(S2CellId.MAX_LEVEL));
    }

    S2CellUnion union = new S2CellUnion();
    union.initFromCellIds(Lists.newArrayList(randomCells));

    // Add one more
    while (!randomCells.add(getRandomCellId(S2CellId.MAX_LEVEL))) {      
    }

    S2CellUnion unionPlusOne = new S2CellUnion();
    unionPlusOne.initFromCellIds(Lists.newArrayList(randomCells));

    assertTrue(unionPlusOne.contains(union));
    assertFalse(union.contains(unionPlusOne));

    // Build the set of parent cells and check containment
    Set<S2CellId> parents = new HashSet<S2CellId>();
    for (S2CellId cellId : union) {
      parents.add(cellId.parent());
    }

    S2CellUnion parentUnion = new S2CellUnion();
    parentUnion.initFromCellIds(Lists.newArrayList(parents));

    assertTrue(parentUnion.contains(union));
    assertFalse(union.contains(parentUnion));
  }

  private void addCells(S2CellId id, boolean selected, List<S2CellId> input,
      ArrayList<S2CellId> expected) {
    // Decides whether to add "id" and/or some of its descendants to the
    // test case. If "selected" is true, then the region covered by "id"
    // *must* be added to the test case (either by adding "id" itself, or
    // some combination of its descendants, or both). If cell ids are to
    // the test case "input", then the corresponding expected result after
    // simplification is added to "expected".

    if (id.equals(S2CellId.none())) {
      // Initial call: decide whether to add cell(s) from each face.
      for (int face = 0; face < 6; ++face) {
        addCells(S2CellId.fromFacePosLevel(face, 0, 0), false, input, expected);
      }
      return;
    }
    if (id.isLeaf()) {
      // The rnd.OneIn() call below ensures that the parent of a leaf cell
      // will always be selected (if we make it that far down the hierarchy).
      assertTrue(selected);
      input.add(id);
      return;
    }
    // The following code ensures that the probability of selecting a cell
    // at each level is approximately the same, i.e. we test normalization
    // of cells at all levels.
    if (!selected && random(S2CellId.MAX_LEVEL - id.level()) != 0) {
      // Once a cell has been selected, the expected output is predetermined.
      // We then make sure that cells are selected that will normalize to
      // the desired output.
      expected.add(id);
      selected = true;
    }

    // With the rnd.OneIn() constants below, this function adds an average
    // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
    // level at which the cell was first selected (level 15 on average).
    // Therefore the average number of input cells in a test case is about
    // (5/6 * 15 * 6) = 75. The average number of output cells is about 6.

    // If a cell is selected, we add it to "input" with probability 5/6.
    boolean added = false;
    if (selected && random(6) != 0) {
      input.add(id);
      added = true;
    }
    int numChildren = 0;
    S2CellId child = id.childBegin();
    for (int pos = 0; pos < 4; ++pos, child = child.next()) {
      // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
      // This intentionally may result in a cell and some of its children
      // being included in the test case.
      //
      // If the cell is not selected, on average we recurse on one child.
      // We also make sure that we do not recurse on all 4 children, since
      // then we might include all 4 children in the input case by accident
      // (in which case the expected output would not be correct).
      if (random(selected ? 12 : 4) == 0 && numChildren < 3) {
        addCells(child, selected, input, expected);
        ++numChildren;
      }
      // If this cell was selected but the cell itself was not added, we
      // must ensure that all 4 children (or some combination of their
      // descendents) are added.
      if (selected && !added) {
        addCells(child, selected, input, expected);
      }
    }
  }

  public void testNormalize() {
    logger.info("TestNormalize");

    // Try a bunch of random test cases, and keep track of average
    // statistics for normalization (to see if they agree with the
    // analysis above).
    S2CellUnion cellunion = new S2CellUnion();
    double inSum = 0, outSum = 0;
    final int kIters = 2000;
    for (int i = 0; i < kIters; ++i) {
      ArrayList<S2CellId> input = Lists.newArrayList();
      ArrayList<S2CellId> expected = Lists.newArrayList();
      addCells(S2CellId.none(), false, input, expected);
      inSum += input.size();
      outSum += expected.size();
      cellunion.initFromCellIds(input);
      assertEquals(cellunion.size(), expected.size());

      assertEquals(expected, cellunion.cellIds());

      // Test GetCapBound().
      S2Cap cap = cellunion.getCapBound();
      for (int  k = 0; k < cellunion.size(); ++k) {
        assertTrue(cap.contains(new S2Cell(cellunion.cellId(k))));
      }

      // Test Contains(S2CellId) and Intersects(S2CellId).
      for (int j = 0; j < input.size(); ++j) {
        assertTrue(cellunion.contains(input.get(j)));
        assertTrue(cellunion.intersects(input.get(j)));
        if (!input.get(j).isFace()) {
          assertTrue(cellunion.intersects(input.get(j).parent()));
          if (input.get(j).level() > 1) {
            assertTrue(cellunion.intersects(input.get(j).parent().parent()));
            assertTrue(cellunion.intersects(input.get(j).parent(0)));
          }
        }
        if (!input.get(j).isLeaf()) {
          assertTrue(cellunion.contains(input.get(j).childBegin()));
          assertTrue(cellunion.intersects(input.get(j).childBegin()));
          assertTrue(cellunion.contains(input.get(j).childEnd().prev()));
          assertTrue(cellunion.intersects(input.get(j).childEnd().prev()));
          assertTrue(cellunion.contains(input.get(j).childBegin(S2CellId.MAX_LEVEL)));
          assertTrue(cellunion.intersects(input.get(j).childBegin(S2CellId.MAX_LEVEL)));
        }
      }
      for (int j = 0; j < expected.size(); ++j) {
        if (!expected.get(j).isFace()) {
          assertTrue(!cellunion.contains(expected.get(j).parent()));
          assertTrue(!cellunion.contains(expected.get(j).parent(0)));
        }
      }

      // Test contains(S2CellUnion) and intersects(S2CellUnion)
      ArrayList<S2CellId> x = Lists.newArrayList();
      ArrayList<S2CellId> y = Lists.newArrayList();
      ArrayList<S2CellId> xOrY = Lists.newArrayList();
      ArrayList<S2CellId> xAndY = Lists.newArrayList();
      for (int j = 0; j < input.size(); ++j) {
        boolean inX = random(2) == 0;
        boolean inY = random(2) == 0;
        if (inX) {
          x.add(input.get(j));
        }
        if (inY) {
          y.add(input.get(j));
        }
        if (inX || inY) {
          xOrY.add(input.get(j));
        }
      }
      S2CellUnion xCells = new S2CellUnion();
      S2CellUnion yCells = new S2CellUnion();
      S2CellUnion xOrYExpected = new S2CellUnion();
      S2CellUnion xAndYExpected = new S2CellUnion();
      xCells.initFromCellIds(x);
      yCells.initFromCellIds(y);
      xOrYExpected.initFromCellIds(xOrY);

      S2CellUnion xOrYCells = new S2CellUnion();
      xOrYCells.getUnion(xCells, yCells);
      assertEquals(xOrYExpected, xOrYCells);

      // Compute the intersection of "x" with each cell of "y",
      // check that this intersection is correct, and append the
      // results to xAndYExpected.
      for (int j = 0; j < yCells.size(); ++j) {
        S2CellId yId = yCells.cellId(j);
        S2CellUnion u = new S2CellUnion();
        u.getIntersection(xCells, yId);
        for (int k = 0; k < xCells.size(); ++k) {
          S2CellId xId = xCells.cellId(k);
          if (xId.contains(yId)) {
            assertEquals(1, u.size());
            assertEquals(yId, u.cellId(0));
          } else if (yId.contains(xId)) {
            if (!u.contains(xId)) {
              u.getIntersection(xCells, yId);
            }
            assertTrue(u.contains(xId));
          }
        }
        for (int k = 0; k < u.size(); ++k) {
          assertTrue(xCells.contains(u.cellId(k)));
          assertTrue(yId.contains(u.cellId(k)));
        }
        xAndY.addAll(u.cellIds());
      }
      xAndYExpected.initFromCellIds(xAndY);

      S2CellUnion xAndYCells = new S2CellUnion();
      xAndYCells.getIntersection(xCells, yCells);
      assertEquals(xAndYExpected, xAndYCells);

      ArrayList<S2CellId> test = Lists.newArrayList();
      ArrayList<S2CellId> dummy = Lists.newArrayList();

      addCells(S2CellId.none(), false, test, dummy);
      for (int j = 0; j < test.size(); ++j) {
        boolean contains = false, intersects = false;
        for (int k = 0; k < expected.size(); ++k) {
          if (expected.get(k).contains(test.get(j))) {
            contains = true;
          }
          if (expected.get(k).intersects(test.get(j))) {
            intersects = true;
          }
        }
        assertEquals(cellunion.contains(test.get(j)), contains);
        assertEquals(cellunion.intersects(test.get(j)), intersects);
      }

    }
  }

  double getMaxAngle(S2CellUnion covering, S2Point axis) {
    double maxAngle = 0;
    for (int i = 0; i < covering.size(); ++i) {
      S2Cell cell = new S2Cell(covering.cellId(i));
      S2Cap cellCap = cell.getCapBound();
      double angle = axis.angle(cellCap.axis()) + cellCap.angle().radians();
      maxAngle = Math.max(maxAngle, angle);
    }
    return maxAngle;
  }

  public void testExpand() {
    logger.info("TestExpand");

    // This test generates coverings for caps of random sizes, and expands
    // the coverings by a random radius, and then make sure that the new
    // covering covers the expanded cap. It also makes sure that the
    // new covering is not too much larger than expected.

    S2RegionCoverer coverer = new S2RegionCoverer();
    for (int i = 0; i < 1000; ++i) {
      S2Cap cap = getRandomCap(S2Cell.averageArea(S2CellId.MAX_LEVEL), 4 * S2.M_PI);

      // Expand the cap by a random factor whose log is uniformly distributed
      // between 0 and log(1e2).
      S2Cap expandedCap =
          S2Cap.fromAxisHeight(cap.axis(), Math.min(2.0, Math.pow(1e2, rand.nextDouble())
              * cap.height()));

      double radius = expandedCap.angle().radians() - cap.angle().radians();
      int maxLevelDiff = random(8);

      S2CellUnion covering = new S2CellUnion();
      coverer.setMaxCells(1 + skewed(10));
      coverer.getCovering(cap, covering);
      checkCovering(cap, covering, true, new S2CellId());

      double maxAngle = getMaxAngle(covering, cap.axis());
      int minLevel = S2CellId.MAX_LEVEL;
      for (int j = 0; j < covering.size(); ++j) {
        minLevel = Math.min(minLevel, covering.cellId(j).level());
      }
      covering.expand(S1Angle.radians(radius), maxLevelDiff);
      checkCovering(expandedCap, covering, false, new S2CellId());

      int expandLevel =
          Math.min(minLevel + maxLevelDiff, S2Projections.MIN_WIDTH.getMaxLevel(radius));
      double expandedMaxAngle = getMaxAngle(covering, cap.axis());

      // If the covering includes a tiny cell along the boundary, in theory the
      // maximum angle of the covering from the cap axis can increase by up to
      // twice the maximum length of a cell diagonal. We allow for an increase
      // of slightly more than this because cell bounding caps are not exact.
      assertTrue(expandedMaxAngle - maxAngle <= 2.01 * S2Projections.MAX_DIAG
          .getValue(expandLevel));
    }
  }

  public void testLeafCellsCovered() {
    S2CellUnion cellUnion = new S2CellUnion();

    // empty union
    assertEquals(0, cellUnion.leafCellsCovered());

    ArrayList<S2CellId> ids = Lists.newArrayList();
    ids.add(S2CellId.fromFacePosLevel(
        0, (1L << ((S2CellId.MAX_LEVEL << 1) - 1)), S2CellId.MAX_LEVEL));

    // One leaf on face 0.
    cellUnion.initFromCellIds(ids);
    assertEquals(1L, cellUnion.leafCellsCovered());

    // Face 0.
    ids.add(S2CellId.fromFacePosLevel(0, 0, 0));
    cellUnion.initFromCellIds(ids);
    assertEquals(1L << 60, cellUnion.leafCellsCovered());

    // Five faces.
    cellUnion.expand(0);
    assertEquals(5L << 60, cellUnion.leafCellsCovered());

    // Whole world.
    cellUnion.expand(0);
    assertEquals(6L << 60, cellUnion.leafCellsCovered());

    // Add some disjoint cells.
    ids.add(S2CellId.fromFacePosLevel(1, 0, 1));
    ids.add(S2CellId.fromFacePosLevel(2, 0, 2));
    ids.add(S2CellId.fromFacePosLevel(2, (1L << 60), 2));
    ids.add(S2CellId.fromFacePosLevel(3, 0, 14));
    ids.add(S2CellId.fromFacePosLevel(4, (1L << 60), 15));
    ids.add(S2CellId.fromFacePosLevel(4, 0, 27));
    ids.add(S2CellId.fromFacePosLevel(5, 0, 30));
    cellUnion.initFromCellIds(ids);
    long expected = 1L + (1L << 6) + (1L << 30) + (1L << 32) + (2L << 56) + (1L << 58) + (1L << 60);
    assertEquals(expected, cellUnion.leafCellsCovered());
  }


  public void testAverageBasedArea() {
    S2CellUnion cellUnion = new S2CellUnion();

    // empty union
    assertEquals(0.0, cellUnion.averageBasedArea());

    ArrayList<S2CellId> ids = Lists.newArrayList();
    ids.add(S2CellId.fromFacePosLevel(1, 0, 1));
    ids.add(S2CellId.fromFacePosLevel(5, 0, 30));
    cellUnion.initFromCellIds(ids);

    double expected = S2Cell.averageArea(S2CellId.MAX_LEVEL) * (1L + (1L << 58));
    assertEquals(expected, cellUnion.averageBasedArea());
  }

  public void testApproxArea() {
    S2CellUnion cellUnion = new S2CellUnion();

    // empty union
    assertEquals(0.0, cellUnion.approxArea());

    ArrayList<S2CellId> ids = Lists.newArrayList();
    ids.add(S2CellId.fromFacePosLevel(1, 0, 1));
    ids.add(S2CellId.fromFacePosLevel(5, 0, 30));
    cellUnion.initFromCellIds(ids);

    double expected = new S2Cell(ids.get(0)).approxArea() + new S2Cell(ids.get(1)).approxArea();
    assertEquals(expected, cellUnion.approxArea());
  }

  public void testExactArea() {
    S2CellUnion cellUnion = new S2CellUnion();

    // empty union
    assertEquals(0.0, cellUnion.exactArea());

    ArrayList<S2CellId> ids = Lists.newArrayList();
    ids.add(S2CellId.fromFacePosLevel(1, 0, 1));
    ids.add(S2CellId.fromFacePosLevel(5, 0, 30));
    cellUnion.initFromCellIds(ids);

    double expected = new S2Cell(ids.get(0)).exactArea() + new S2Cell(ids.get(1)).exactArea();
    assertEquals(expected, cellUnion.averageBasedArea());
  }
}
