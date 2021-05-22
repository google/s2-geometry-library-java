/*
 * Copyright 2018 Google Inc.
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
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@GwtCompatible
public class S2ShapeIndexRegionTest extends GeometryTestCase {
  /**
   * When building the S2ShapeIndex, each UV rectangle around an index cell is padded slightly to
   * ensure edges are definitely included in each index cell they intersect, despite having
   * imprecise math doing the tests. This constant is slightly larger than that, so expanding a cell
   * rect outward by this amount will certainly insert it into the surrounding index cells, and
   * contracting by this amount will certainly exclude it from the surrounding index cells.
   */
  private static final double CELL_PADDING =
      2 * (S2EdgeUtil.FACE_CLIP_ERROR_UV_COORD + S2EdgeUtil.INTERSECTS_RECT_ERROR_UV_DIST);

  private S2ShapeIndex index;

  @Override
  protected void setUp() {
    index = new S2ShapeIndex();
  }

  public void testGetCapBound() {
    S2CellId id = S2CellId.fromDebugString("3/0123012301230123012301230123");

    // Add a polygon that is slightly smaller than the cell being tested.
    index.add(paddedCell(id, -CELL_PADDING));
    S2Cap cellBound = new S2Cell(id).getCapBound();
    S2Cap indexBound = new S2ShapeIndexRegion(index).getCapBound();
    assertTrue(indexBound.contains(cellBound));

    // Note that S2CellUnion.getCapBound returns a slightly larger bound than S2Cell.getBound even
    // when the cell union consists of a single S2CellId.
    assertTrue(
        indexBound.radius().toAngle().radians()
            <= 1.00001 * cellBound.radius().toAngle().radians());
  }

  public void testGetRectBound() {
    S2CellId id = S2CellId.fromDebugString("3/0123012301230123012301230123");

    // Add a polygon that is slightly smaller than the cell being tested.
    index.add(paddedCell(id, -CELL_PADDING));
    S2LatLngRect cellBound = new S2Cell(id).getRectBound();
    S2LatLngRect indexBound = new S2ShapeIndexRegion(index).getRectBound();
    assertEquals(indexBound, cellBound);
  }

  public void testGetCellUnionBoundMultipleFaces() {
    S2CellId[] ids = {S2CellId.fromDebugString("3/00123"), S2CellId.fromDebugString("2/11200013")};
    for (S2CellId id : ids) {
      index.add(paddedCell(id, -CELL_PADDING));
    }
    List<S2CellId> covering = new ArrayList<>();
    new S2ShapeIndexRegion(index).getCellUnionBound(covering);
    Arrays.sort(ids);
    assertEquals(Arrays.asList(ids), covering);
  }

  public void testGetCellUnionBoundOneFace() {
    // This tests consists of 3 pairs of S2CellIds.  Each pair is located within
    // one of the children of face 5, namely the cells 5/0, 5/1, and 5/3.
    // We expect GetCellUnionBound to compute the smallest cell that bounds the
    // pair on each face.
    S2CellId[] input = {
      S2CellId.fromDebugString("5/010"),
      S2CellId.fromDebugString("5/0211030"),
      S2CellId.fromDebugString("5/110230123"),
      S2CellId.fromDebugString("5/11023021133"),
      S2CellId.fromDebugString("5/311020003003030303"),
      S2CellId.fromDebugString("5/311020023"),
    };
    S2CellId[] expected = {
      S2CellId.fromDebugString("5/0"),
      S2CellId.fromDebugString("5/110230"),
      S2CellId.fromDebugString("5/3110200")
    };
    for (S2CellId id : input) {
      // Add each shape 3 times to ensure that the S2ShapeIndex subdivides.
      for (int copy = 0; copy < 3; ++copy) {
        index.add(paddedCell(id, -CELL_PADDING));
      }
    }
    List<S2CellId> actual = new ArrayList<>();
    new S2ShapeIndexRegion(index).getCellUnionBound(actual);
    assertEquals(Arrays.asList(expected), actual);
  }

  public void testContainsCellMultipleShapes() {
    // Create a debug cell, and shapes slightly smaller and larger than that cell.
    S2CellId id = S2CellId.fromDebugString("3/0123012301230123012301230123");
    S2Shape slightlySmaller = paddedCell(id, -CELL_PADDING);
    S2Shape slightlyLarger = paddedCell(id, CELL_PADDING);

    // Check that the index region of the smaller shape doesn't contain the cell.
    index.add(slightlySmaller);
    assertFalse(new S2ShapeIndexRegion(index).contains(new S2Cell(id)));

    // Check that the index region of the larger and smaller shapes does contain the cell.
    // Note that contains() should return true if *any* shape contains the cell.
    index = new S2ShapeIndex();
    index.add(slightlySmaller);
    index.add(slightlyLarger);
    assertTrue(new S2ShapeIndexRegion(index).contains(new S2Cell(id)));

    // Verify that all children of the cell are also contained.
    for (S2CellId child : id.childrenAtLevel(id.level() + 1)) {
      assertTrue(new S2ShapeIndexRegion(index).contains(new S2Cell(child)));
    }
  }

  public void testIntersectsShrunkenCell() {
    S2CellId target = S2CellId.fromDebugString("3/0123012301230123012301230123");

    // Add a polygon that is slightly smaller than the cell being tested.
    index.add(paddedCell(target, -CELL_PADDING));
    S2Region region = new S2ShapeIndexRegion(index);

    // Check that the index intersects the cell itself, but not any of the
    // neighboring cells.
    assertTrue(region.mayIntersect(new S2Cell(target)));
    List<S2CellId> nbrs = new ArrayList<>();
    target.getAllNeighbors(target.level(), nbrs);
    for (S2CellId id : nbrs) {
      assertFalse(region.mayIntersect(new S2Cell(id)));
    }
  }

  public void testIntersectsExactCell() {
    S2CellId target = S2CellId.fromDebugString("3/0123012301230123012301230123");

    // Adds a polygon that exactly follows a cell boundary.
    index.add(paddedCell(target, 0.0));
    S2Region region = new S2ShapeIndexRegion(index);

    // Check that the index intersects the cell and all of its neighbors.
    List<S2CellId> ids = new ArrayList<>();
    ids.add(target);
    target.getAllNeighbors(target.level(), ids);
    for (S2CellId id : ids) {
      assertTrue(region.mayIntersect(new S2Cell(id)));
    }
  }

  private static S2Shape paddedCell(S2CellId id, double paddingUv) {
    R2Rect uv = id.getBoundUV().expanded(paddingUv);
    return new S2Loop(
        new AbstractList<S2Point>() {
          @Override
          public int size() {
            return 4;
          }

          @Override
          public S2Point get(int offset) {
            return S2Point.normalize(S2Projections.faceUvToXyz(id.face(), uv.getVertex(offset)));
          }
        });
  }
}
