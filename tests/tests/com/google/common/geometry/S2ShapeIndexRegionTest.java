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

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.List;

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
    super.setUp();
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
    // This tests consists of 3 pairs of S2CellIds. Each pair is located within one of the children
    // of face 5, namely the cells 5/0, 5/1, and 5/3. We expect GetCellUnionBound to compute the
    // smallest cell that bounds the pair on each face.
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

    // Check that the index intersects the cell itself, but not any of the neighboring cells.
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
            return S2Projections.faceUvToXyz(id.face(), uv.getVertex(offset)).normalize();
          }
        });
  }

  /**
   * Tests that VisitIntersectingShapes() produces results that are consistent with mayIntersect()
   * and contains() for the given S2ShapeIndex.  It tests all cells in the given index, all
   * ancestors of those cells, and a randomly chosen subset of descendants of those cells.
   */
  class VisitIntersectingShapesTest {
    private final S2ShapeIndex index;
    private final S2Iterator<S2ShapeIndex.Cell> iter;
    private final S2ShapeIndexRegion region;
    private final IdentityHashMap<S2Shape, S2ShapeIndex> shapeIndexes;

    public VisitIntersectingShapesTest(S2ShapeIndex index) {
      this.index = index;
      this.iter = index.iterator();
      this.region = new S2ShapeIndexRegion(index);

      // Create an S2ShapeIndex for each shape in the original index, so that we can use
      // mayIntersect() and contains() to determine the status of individual shapes.
      shapeIndexes = new IdentityHashMap<>();
      for (S2Shape shape : index.getShapes()) {
        S2ShapeIndex shapeIndex = new S2ShapeIndex();
        shapeIndex.add(shape);
        shapeIndexes.put(shape, shapeIndex);
      }
    }

    public void run() {
      for (S2CellId id = S2CellId.begin(0); !id.equals(S2CellId.end(0)); id = id.next()) {
        testCell(new S2Cell(id));
      }
    }

    private void testCell(S2Cell target) {
      // Indicates whether each shape that intersects "target" also contains it.
      IdentityHashMap<S2Shape, Boolean> shapeContains = new IdentityHashMap<>();
      assertTrue(region.visitIntersectingShapes(target, new S2ShapeIndexRegion.ShapeVisitor() {
        @Override
        public boolean test(S2Shape shape, boolean containsTarget) {
          // Verify that each shape is visited at most once.
          assertFalse(shapeContains.containsKey(shape));
          shapeContains.put(shape, containsTarget);
          return true;
        }
      }));

      for (S2Shape shape : index.getShapes()) {
        S2ShapeIndexRegion shapeRegion = new S2ShapeIndexRegion(shapeIndexes.get(shape));
        if (!shapeRegion.mayIntersect(target)) {
          assertFalse(shapeContains.containsKey(shape));
        } else {
          assertEquals(shapeRegion.contains(target), shapeContains.get(shape).booleanValue());
        }
      }

      switch (iter.locate(target.id())) {
        case DISJOINT:
          return;

        case SUBDIVIDED:
          S2Cell[] children = new S2Cell[4];
          for (int i = 0; i < 4; ++i) {
            children[i] = new S2Cell();
          }
          assertTrue(target.subdivide(children));
          for (S2Cell child : children) {
            testCell(child);
          }
          return;

        case INDEXED:
          // We check a few random descendant cells by continuing randomly down one branch of the
          // tree for a few levels.
          if (target.isLeaf() || data.oneIn(3)) {
            return;
          }
          testCell(new S2Cell(target.id().child(data.uniform(4))));
          return;
      }
    }
  }

  /** Test visitIntersectingShapes with a single shape of 100 points. */
  public void testVisitIntersectingShapes_points() {
    List<S2Point> vertices = new ArrayList<>();
    for (int i = 0; i < 100; ++i) {
      vertices.add(data.getRandomPoint());
    }
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(S2Point.Shape.fromList(vertices));
    new VisitIntersectingShapesTest(index).run();
  }

  /** Test visitIntersectingShapes with 50 randomly located regular polyline shapes. */
  public void testVisitIntersectingShapes_polylines() {
    S2ShapeIndex index = new S2ShapeIndex();
    S2Cap centerCap = S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.radians(0.5));
    for (int i = 0; i < 50; ++i) {
      S2Point center = data.samplePoint(centerCap);
      List<S2Point> vertices = new ArrayList<>();
      if (data.oneIn(10)) {
        // Try a few degenerate polylines, having just two identical vertices.
        vertices.add(center);
        vertices.add(center);
      } else {
        vertices =
            S2Loop.makeRegularVertices(
                center, S1Angle.radians(data.nextDouble()), data.uniform(20) + 3);
      }
      index.add(new S2Polyline(vertices));
    }
    new VisitIntersectingShapesTest(index).run();
  }

  /**
   * Test visitIntersectingShapes with 10 fractal polygon loops plus one big polygon that contains
   * most of the others.
   */
  public void testVisitIntersectingShapes_polygons() {
    S2ShapeIndex index = new S2ShapeIndex();
    S2Cap centerCap = S2Cap.fromAxisAngle(S2Point.X_POS, S1Angle.radians(0.5));
    S2FractalBuilder fractalBuilder = new S2FractalBuilder(data.rand);
    for (int i = 0; i < 10; ++i) {
      fractalBuilder.setLevelForApproxMaxEdges(3 * 64);
      S2Point center = data.samplePoint(centerCap);
      index.add(
          fractalBuilder.makeLoop(
              data.getRandomFrameAt(center), S1Angle.radians(data.nextDouble())));
    }
    // Also add a big polygon containing most of the polygons above to ensure that we test
    // containment of cells that are ancestors of index cells.
    index.add(paddedCell(S2CellId.fromFace(0), 0));
    new VisitIntersectingShapesTest(index).run();
  }
}