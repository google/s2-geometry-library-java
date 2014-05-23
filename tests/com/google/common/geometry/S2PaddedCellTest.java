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

import com.google.common.annotations.GwtCompatible;
import com.google.common.geometry.R2Rect.Axis;

/**
 * Verifies S2PaddedCell.
 */
@GwtCompatible
public class S2PaddedCellTest extends GeometryTestCase {
  public void testS2CellMethods() {
    // Test the S2PaddedCell methods that have approximate S2Cell equivalents.
    final int iters = 1000;
    for (int iter = 0; iter < iters; ++iter) {
      S2CellId id = getRandomCellId();
      double padding = Math.pow(1e-15, rand.nextDouble());
      S2Cell cell = new S2Cell(id);
      S2PaddedCell pcell = new S2PaddedCell(id, padding);
      compareS2CellToPadded(cell, pcell, padding);
      if (!id.isLeaf()) {
        S2Cell[] children = new S2Cell[] {new S2Cell(), new S2Cell(), new S2Cell(), new S2Cell()};
        assertTrue(cell.subdivide(children));
        for (int pos = 0; pos < 4; ++pos) {
          int ij = S2.posToIJ(pcell.orientation(), pos);
          compareS2CellToPadded(children[pos], new S2PaddedCell(pcell, ij >> 1, ij & 1), padding);
        }
      }
    }
  }

  public void testGetEntryExitVertices() {
    final int iters = 1000;
    for (int iter = 0; iter < iters; ++iter) {
      S2CellId id = getRandomCellId();
      // Check that entry/exit vertices do not depend on padding.
      assertEquals(new S2PaddedCell(id, 0).getEntryVertex(),
          new S2PaddedCell(id, 0.5).getEntryVertex());
      assertEquals(new S2PaddedCell(id, 0).getExitVertex(),
          new S2PaddedCell(id, 0.5).getExitVertex());

      // Check that the exit vertex of one cell is the same as the entry vertex of the immediately
      // following cell.  (This also tests wrapping from the end to the start of the S2CellId curve
      // with high probability.)
      assertEquals(new S2PaddedCell(id, 0).getExitVertex(),
          new S2PaddedCell(id.nextWrap(), 0).getEntryVertex());

      // Check that the entry vertex of a cell is the same as the entry vertex of its first child,
      // and similarly for the exit vertex.
      if (!id.isLeaf()) {
        assertEquals(new S2PaddedCell(id, 0).getEntryVertex(),
            new S2PaddedCell(id.child(0), 0).getEntryVertex());
        assertEquals(new S2PaddedCell(id, 0).getExitVertex(),
            new S2PaddedCell(id.child(3), 0).getExitVertex());
      }
    }
  }

  public void testShrinkToFit() {
    final int iters = 1000;
    for (int iter = 0; iter < iters; ++iter) {
      // Start with the desired result and work backwards.
      S2CellId result = getRandomCellId();
      R2Rect uvResult = result.getBoundUV();
      R2Vector uvSize = uvResult.getSize();

      // Find the biggest rectangle that fits in "result" after padding.
      // (These calculations ignore numerical errors.)
      double maxPadding = 0.5 * Math.min(uvSize.x, uvSize.y);
      double padding = maxPadding * rand.nextDouble();
      R2Rect maxRect = uvResult.expanded(-padding);

      // Start with a random subset of the maximum rectangle.
      R2Rect rect = R2Rect.fromPointPair(
          new R2Vector(sampleInterval(maxRect.x()), sampleInterval(maxRect.y())),
          new R2Vector(sampleInterval(maxRect.x()), sampleInterval(maxRect.y())));
      if (!result.isLeaf()) {
        // If the result is not a leaf cell, we must ensure that no child of "result" also satisfies
        // the conditions of ShrinkToFit().  We do this by ensuring that "rect" intersects at least
        // two children of "result" (after padding).
        Axis axis = rand.nextInt(2) == 0 ? Axis.X : Axis.Y;
        double center = result.getCenterUV().get(axis.ordinal());

        // Find the range of coordinates that are shared between child cells along that axis.
        R1Interval shared = new R1Interval(center - padding, center + padding);
        double mid = sampleInterval(shared.intersection(maxRect.getInterval(axis)));
        rect.getInterval(axis).set(
            sampleInterval(new R1Interval(maxRect.getInterval(axis).lo(), mid)),
            sampleInterval(new R1Interval(mid, maxRect.getInterval(axis).hi())));
      }

      // Choose an arbitrary ancestor as the S2PaddedCell.
      S2CellId initialId = result.parent(rand.nextInt(result.level() + 1));
      assertEquals(result, new S2PaddedCell(initialId, padding).shrinkToFit(rect));
    }
  }

  private double sampleInterval(R1Interval x) {
    return uniform(x.lo(), x.hi());
  }

  private static void compareS2CellToPadded(S2Cell cell, S2PaddedCell pcell, double padding) {
    assertEquals(cell.id(), pcell.id());
    assertEquals(cell.level(), pcell.level());
    assertEquals(padding, pcell.padding());
    assertEquals(cell.getBoundUV().expanded(padding), pcell.bound());
    R2Vector uvCenter = cell.id().getCenterUV();
    assertEquals(R2Rect.fromPoint(uvCenter).expanded(padding), pcell.middle());
    assertEquals(cell.getCenter(), pcell.getCenter());
  }
}