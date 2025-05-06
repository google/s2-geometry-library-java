/*
 * Copyright 2024 Google Inc.
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

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2IndexCellData.EdgeAndIdChain;
import com.google.common.geometry.S2ShapeIndex.Cell;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2IndexCellData}. */
@RunWith(JUnit4.class)
public class S2IndexCellDataTest extends GeometryTestCase {
  S2ShapeIndex index =
      S2TextFormat.makeIndexOrDie(
          "0:0" // One point
              + "#1:1, 2:2" // One line
              + "#1:0, 0:1, -1:0, 0:-1" // One triangle.
          );

  /**
   * A trivial implementation of EdgeAndIdChain for testing the default methods on the interface.
   */
  private static final class EdgeAndIdChainImpl implements EdgeAndIdChain {
    final S2Point start;
    final S2Point end;
    final int edgeId;
    final int chainId;
    final int offset;

    public EdgeAndIdChainImpl(S2Point start, S2Point end, int edgeId, int chainId, int offset) {
      this.start = start;
      this.end = end;
      this.edgeId = edgeId;
      this.chainId = chainId;
      this.offset = offset;
    }

    @Override
    public S2Point start() {
      return start;
    }
    @Override
    public S2Point end() {
      return end;
    }
    @Override
    public int edgeId() {
      return edgeId;
    }
    @Override
    public int chainId() {
      return chainId;
    }
    @Override
    public int offset() {
      return offset;
    }
  }

  @Test
  public void testEdgeAndIdChainInterfaceDefaultMethods() {
    // e1 and e2 are equal, because only points are used for equality. e3 has a different endpoint
    // as e1 amd e2. e4 has a different start point as e1 and e2, with the reverse endpoints as e3.
    EdgeAndIdChainImpl e1 =
        new EdgeAndIdChainImpl(new S2Point(1, 0, 0), new S2Point(0, 1, 0), 1, 2, 3);
    EdgeAndIdChainImpl e2 =
        new EdgeAndIdChainImpl(new S2Point(1, 0, 0), new S2Point(0, 1, 0), 4, 5, 6);
    EdgeAndIdChainImpl e3 =
        new EdgeAndIdChainImpl(new S2Point(1, 0, 0), new S2Point(0, 0, 1), 7, 8, 9);
    EdgeAndIdChainImpl e4 =
        new EdgeAndIdChainImpl(new S2Point(0, 0, 1), new S2Point(1, 0, 0), 10, 20, 30);

    assertTrue(e1.isEqualTo(e2));
    assertFalse(e1.isEqualTo(e3));
    assertTrue(EdgeAndIdChain.equals(e1, e2));
    assertFalse(EdgeAndIdChain.equals(e1, e3));
    assertFalse(EdgeAndIdChain.equals(e2, e4));

    assertTrue(e3.isReverseOf(e4));
    assertTrue(e4.isReverseOf(e3));
    assertFalse(e1.isReverseOf(e2));
    assertFalse(e2.isReverseOf(e1));

    // Ordering is by starting vertex, breaking ties by ending vertex. S2Point.lessThan() is
    // lexicographic, comparing X, Y, Z. So e4 is the smallest, followed by e3, then e1 and e2 are
    // tied.
    assertTrue(EdgeAndIdChain.lessThan(e4, e3));
    assertTrue(EdgeAndIdChain.lessThan(e3, e1));
    assertFalse(EdgeAndIdChain.lessThan(e1, e2));
    assertFalse(EdgeAndIdChain.lessThan(e2, e1));

    // Verify hasEndpoint() works.
    assertTrue(e1.hasEndpoint(new S2Point(1, 0, 0)));
    assertTrue(e1.hasEndpoint(new S2Point(0, 1, 0)));
    assertFalse(e1.hasEndpoint(new S2Point(0, 0, 1)));
  }

  // Check that we get all dimensions by default.
  @Test
  public void testDimensionFilteringAllDimensionsDefault() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.loadCell(index, iter.id(), iter.entry());
    assertFalse(data.dimEdges(0).isEmpty());
    assertFalse(data.dimEdges(1).isEmpty());
    assertFalse(data.dimEdges(2).isEmpty());
  }

  // No dimensions should work too, we just don't decode edges.
  @Test
  public void testDimensionFilteringNoDimensions() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.setDimWanted(0, false);
    data.setDimWanted(1, false);
    data.setDimWanted(2, false);
    data.loadCell(index, iter.id(), iter.entry());
    assertTrue(data.dimEdges(0).isEmpty());
    assertTrue(data.dimEdges(1).isEmpty());
    assertTrue(data.dimEdges(2).isEmpty());
  }

  // Should be able to get an empty range for a dimension that is turned off.
  @Test
  public void testDimensionFilteringRangeForOffDimension() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.setDimWanted(0, false);
    data.setDimWanted(1, true);
    data.setDimWanted(2, true);
    data.loadCell(index, iter.id(), iter.entry());
    assertTrue(data.dimRangeEdges(0, 0).isEmpty());
    assertFalse(data.dimRangeEdges(0, 2).isEmpty());
  }

  @Test
  public void testDimensionFilteringRangeForLines() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.setDimWanted(0, false);
    data.setDimWanted(1, true);
    data.setDimWanted(2, false);
    data.loadCell(index, iter.id(), iter.entry());
    assertTrue(data.dimEdges(0).isEmpty());
    assertFalse(data.dimEdges(1).isEmpty());
    assertTrue(data.dimEdges(2).isEmpty());
  }

  @Test
  public void testDimensionFilteringRangeForPointsAndPolygons() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.setDimWanted(0, true);
    data.setDimWanted(1, false);
    data.setDimWanted(2, true);
    data.loadCell(index, iter.id(), iter.entry());
    assertFalse(data.dimEdges(0).isEmpty());
    assertTrue(data.dimEdges(1).isEmpty());
    assertFalse(data.dimEdges(2).isEmpty());
  }

  @Test
  public void testDimensionFilteringRangeForPoints() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.setDimWanted(0, true);
    data.setDimWanted(1, false);
    data.setDimWanted(2, false);
    data.loadCell(index, iter.id(), iter.entry());
    assertFalse(data.dimEdges(0).isEmpty());
    assertTrue(data.dimEdges(1).isEmpty());
    assertTrue(data.dimEdges(2).isEmpty());
  }

  @Test
  public void testDimensionFilteringRangeForPolygons() {
    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = index.iterator();
    data.setDimWanted(0, false);
    data.setDimWanted(1, false);
    data.setDimWanted(2, true);
    data.loadCell(index, iter.id(), iter.entry());
    assertTrue(data.dimEdges(0).isEmpty());
    assertTrue(data.dimEdges(1).isEmpty());
    assertFalse(data.dimEdges(2).isEmpty());
  }

  @Test
  public void testCellAndCenterRecomputed() {
    // A line between two faces will guarantee we get at least two cells.
    S2ShapeIndex twoFacesIndex = S2TextFormat.makeIndexOrDie("# 0:0, 0:-90 #");

    S2IndexCellData data = new S2IndexCellData();
    S2Iterator<Cell> iter = twoFacesIndex.iterator();

    data.loadCell(index, iter.id(), iter.entry());
    S2Point center0 = data.center();
    S2Cell cell0 = data.cell();

    iter.next();
    assertFalse(iter.done());

    data.loadCell(index, iter.id(), iter.entry());
    S2Point center1 = data.center();
    S2Cell cell1 = data.cell();

    assertFalse(cell0.equals(cell1));
    assertFalse(center0.equalsPoint(center1));

    // Load the same cell again, nothing should change.
    data.loadCell(index, iter.id(), iter.entry());
    S2Point center2 = data.center();
    S2Cell cell2 = data.cell();

    // Test with reference equality: the cell and point should be the same objects as before.
    assertSame(cell1, cell2);
    assertSame(center1, center2);
  }
}
