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

import static com.google.common.geometry.S2Projections.PROJ;

import com.google.common.annotations.GwtCompatible;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/** Tests for S2CellUnion. */
@GwtCompatible
public strictfp class S2CellUnionTest extends GeometryTestCase {
  public void testBasic() {
    S2CellUnion empty = new S2CellUnion();
    ArrayList<S2CellId> ids = Lists.newArrayList();
    empty.initFromCellIds(ids);
    assertEquals(0, empty.size());

    S2CellId face1Id = S2CellId.fromFace(1);
    S2CellUnion face1Union = new S2CellUnion();
    ids.add(face1Id);
    face1Union.initFromCellIds(ids);
    assertEquals(1, face1Union.size());
    assertEquals(face1Id, face1Union.cellId(0));

    S2CellId face2Id = S2CellId.fromFace(2);
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
        addCells(S2CellId.fromFace(face), false, input, expected);
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
    // of 5/6 * (MAX_LEVEL - level) cells to "input" where "level" is the
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

  /** Returns the maximum geodesic distance from "axis" to any point of "covering". */
  private static double getRadius(S2CellUnion covering, S2Point axis) {
    double maxDist = 0;
    for (int i = 0; i < covering.size(); ++i) {
      S2Cell cell = new S2Cell(covering.cellId(i));
      for (int j = 0; j < 4; ++j) {
        S2Point a = cell.getVertex(j);
        S2Point b = cell.getVertex((j + 1) & 3);
        // The maximum distance is not always attained at a cell vertex: if at
        // least one vertex is in the opposite hemisphere from "axis" then the
        // maximum may be attained along an edge.  We solve this by computing
        // the minimum distance from the edge to (-axis) instead.  We can't
        // simply do this all the time because S2EdgeUtil.getDistance() has
        // poor accuracy when the result is close to Pi.
        //
        // TODO(user): Improve S2EdgeUtil.getDistance() accuracy near Pi.
        double dist;
        if (a.angle(axis) > S2.M_PI_2 || b.angle(axis) > S2.M_PI_2) {
          dist = S2.M_PI - S2EdgeUtil.getDistance(S2Point.neg(axis), a, b).radians();
        } else {
          dist = a.angle(axis);
        }
        maxDist = Math.max(maxDist, dist);
      }
    }
    return maxDist;
  }

  public void testExpand() {
    // This test generates coverings for caps of random sizes, expands
    // the coverings by a random radius, and then make sure that the new
    // covering covers the expanded cap.  It also makes sure that the
    // new covering is not too much larger than expected.

    S2RegionCoverer coverer = new S2RegionCoverer();
    for (int i = 0; i < 1000; i++) {
      S2Cap cap = getRandomCap(S2Cell.averageArea(S2CellId.MAX_LEVEL), 4 * S2.M_PI);

      // Expand the cap area by a random factor whose log is uniformly
      // distributed between 0 and log(1e2).
      S2Cap expandedCap = S2Cap.fromAxisHeight(
          cap.axis(), Math.min(2.0, Math.pow(1e2, rand.nextDouble()) * cap.height()));

      double radius = expandedCap.angle().radians() - cap.angle().radians();
      int maxLevelDiff = rand.nextInt(8);

      // Generate a covering for the original cap, and measure the maximum
      // distance from the cap center to any point in the covering.
      coverer.setMaxCells(1 + skewed(10));
      S2CellUnion covering = coverer.getCovering(cap);
      checkCovering(cap, covering, true);
      double coveringRadius = getRadius(covering, cap.axis());

      // This code duplicates the logic in Expand(minRadius, maxLevelDiff)
      // that figures out an appropriate cell level to use for the expansion.
      int minLevel = S2CellId.MAX_LEVEL;
      for (int j = 0; j < covering.size(); ++j) {
        minLevel = Math.min(minLevel, covering.cellId(j).level());
      }
      int expandLevel = Math.min(minLevel + maxLevelDiff, PROJ.minWidth.getMaxLevel(radius));

      // Generate a covering for the expanded cap, and measure the new maximum
      // distance from the cap center to any point in the covering.
      covering.expand(S1Angle.radians(radius), maxLevelDiff);
      checkCovering(expandedCap, covering, false);
      double expandedCoveringRadius = getRadius(covering, cap.axis());

      // If the covering includes a tiny cell along the boundary, in theory the
      // maximum angle of the covering from the cap axis can increase by up to
      // twice the maximum length of a cell diagonal.
      assertTrue(expandedCoveringRadius - coveringRadius
          <= 2 * PROJ.maxDiag.getValue(expandLevel));
    }
  }

  private static void checkInitFromMinMax(S2CellId minId, S2CellId maxId) {
    S2CellUnion cellUnion = new S2CellUnion();
    cellUnion.initFromMinMax(minId, maxId);
    List<S2CellId> cellIds = cellUnion.cellIds();
    assertTrue(cellIds.size() > 0);
    assertEquals(minId, cellIds.get(0).rangeMin());
    assertEquals(maxId, cellIds.get(cellIds.size() - 1).rangeMax());
    for (int i = 1; i < cellIds.size(); i++) {
      assertEquals(cellIds.get(i).rangeMin(), cellIds.get(i - 1).rangeMax().next());
    }
    assertFalse(cellUnion.normalize());
  }

  public void testInitFromRange() {
    // Check the very first leaf cell and face cell.
    S2CellId face1Id = S2CellId.fromFace(0);
    checkInitFromMinMax(face1Id.rangeMin(), face1Id.rangeMin());
    checkInitFromMinMax(face1Id.rangeMin(), face1Id.rangeMax());

    // Check the very last leaf cell and face cell.
    S2CellId idFace5 = S2CellId.fromFace(5);
    checkInitFromMinMax(idFace5.rangeMin(), idFace5.rangeMax());
    checkInitFromMinMax(idFace5.rangeMax(), idFace5.rangeMax());

    // Check random ranges of leaf cells.
    for (int i = 0; i < 100; ++i) {
      S2CellId x = getRandomCellId(S2CellId.MAX_LEVEL);
      S2CellId y = getRandomCellId(S2CellId.MAX_LEVEL);
      if (x.compareTo(y) > 0) {
        S2CellId temp = x;
        x = y;
        y = temp;
      }
      checkInitFromMinMax(x, y);
    }
  }

  public void testInitFromBeginEnd() {
    // Since InitFromRange() is implemented in terms of InitFromBeginEnd(), we
    // focus on test cases that generate an empty range.
    ArrayList<S2CellId> initialIds = Lists.newArrayList(S2CellId.fromFace(3));
    S2CellUnion cellUnion = new S2CellUnion();

    // Test an empty range before the minimum S2CellId.
    S2CellId idBegin = S2CellId.begin(S2CellId.MAX_LEVEL);
    cellUnion.initFromCellIds(initialIds);
    cellUnion.initFromBeginEnd(idBegin, idBegin);
    assertEquals(0, cellUnion.size());

    // Test an empty range after the maximum S2CellId.
    S2CellId idEnd = S2CellId.end(S2CellId.MAX_LEVEL);
    cellUnion.initFromCellIds(initialIds);
    cellUnion.initFromBeginEnd(idEnd, idEnd);
    assertEquals(0, cellUnion.size());

    // Test the full sphere.
    cellUnion.initFromBeginEnd(idBegin, idEnd);
    assertEquals(6, cellUnion.size());
    for (int i = 0; i < cellUnion.size(); ++i) {
      assertTrue(cellUnion.cellId(i).isFace());
    }
  }

  public void testLeafCellsCovered() {
    S2CellUnion cellUnion = new S2CellUnion();
    assertEquals(0, cellUnion.leafCellsCovered());

    ArrayList<S2CellId> ids = Lists.newArrayList();
    // One leaf cell on face 0.
    ids.add(S2CellId.fromFace(0).childBegin(S2CellId.MAX_LEVEL));
    cellUnion.initFromCellIds(ids);
    assertEquals(1L, cellUnion.leafCellsCovered());

    // Face 0 itself (which includes the previous leaf cell).
    ids.add(S2CellId.fromFace(0));
    cellUnion.initFromCellIds(ids);
    assertEquals(1L << 60, cellUnion.leafCellsCovered());

    // Five faces.
    cellUnion.expand(0);
    assertEquals(5L << 60, cellUnion.leafCellsCovered());

    // Whole world.
    cellUnion.expand(0);
    assertEquals(6L << 60, cellUnion.leafCellsCovered());

    // Add some disjoint cells.
    ids.add(S2CellId.fromFace(1).childBegin(1));
    ids.add(S2CellId.fromFace(2).childBegin(2));
    ids.add(S2CellId.fromFace(2).childEnd(2).prev());
    ids.add(S2CellId.fromFace(3).childBegin(14));
    ids.add(S2CellId.fromFace(4).childBegin(27));
    ids.add(S2CellId.fromFace(4).childEnd(15).prev());
    ids.add(S2CellId.fromFace(5).childBegin(30));
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
