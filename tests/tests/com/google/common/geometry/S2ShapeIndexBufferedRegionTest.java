/*
 * Copyright 2022 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;

public final class S2ShapeIndexBufferedRegionTest extends GeometryTestCase {

  @Override
  protected void setUp() {
    super.setUp();
  }

  public void testEmptyIndex() {
    // Test buffering an empty S2ShapeIndex.
    S2ShapeIndex index = new S2ShapeIndex();
    S1ChordAngle radius = S1ChordAngle.fromDegrees(2);
    S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, radius);
    S2RegionCoverer coverer = S2RegionCoverer.builder().build();
    S2CellUnion covering = coverer.getCovering(region);
    assertTrue(covering.isEmpty());
  }

  public void testFullPolygon() {
    // Test buffering an S2ShapeIndex that contains a full polygon.
    S2ShapeIndex index = makeIndexOrDie("# # full");
    S1ChordAngle radius = S1ChordAngle.fromDegrees(2);
    S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, radius);
    S2RegionCoverer coverer = S2RegionCoverer.builder().build();
    S2CellUnion covering = coverer.getCovering(region);
    assertEquals(6, covering.size());
    for (S2CellId id : covering.cellIds()) {
      assertTrue(id.isFace());
    }
  }

  public void testFullAfterBuffering() {
    // Test a region that becomes the full polygon after buffering.
    S2ShapeIndex index = makeIndexOrDie("0:0 | 0:90 | 0:180 | 0:-90 | 90:0 | -90:0 # #");
    S1ChordAngle radius = S1ChordAngle.fromDegrees(60);
    S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, radius);
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(1000).build();
    S2CellUnion covering = coverer.getCovering(region);
    assertEquals(6, covering.size());
    for (S2CellId id : covering.cellIds()) {
      assertTrue(id.isFace());
    }
  }

  public void testPointZeroRadius() {
    // Test that buffering a point using a zero radius produces a non-empty covering. (This
    // requires using "less than or equal to" distance tests.)
    S2ShapeIndex index = makeIndexOrDie("34:25 # #");
    S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, S1ChordAngle.ZERO);
    S2RegionCoverer coverer = S2RegionCoverer.builder().build();
    S2CellUnion covering = coverer.getCovering(region);
    assertEquals(1, covering.size());
    for (S2CellId id : covering.cellIds()) {
      assertTrue(id.isLeaf());
    }
  }

  public void testBufferedPointVsCap() {
    // Compute an S2Cell covering of a buffered S2Point, then make sure that the covering is
    // equivalent to the corresponding S2Cap.
    S2ShapeIndex index = makeIndexOrDie("3:5 # #");
    S2Point point = makePointOrDie("3:5");
    S1ChordAngle radius = S1ChordAngle.fromDegrees(2);
    S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, radius);
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(50).build();
    S2CellUnion covering = coverer.getCovering(region);
    S2Cap equivalentCap = S2Cap.fromAxisChord(point, radius);
    checkCovering(equivalentCap, covering, true);
  }

  // Puts the given shape in an index and creates a buffered covering. Returns a polygon constructed
  // from the resulting cells.
  S2Polygon getBufferedIndexPolygon(S2Shape shape, S1ChordAngle radius) {
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(shape);
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(100).build();
    S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, radius);
    S2CellUnion covering = coverer.getCovering(region);
    System.err.println(
        "Covering uses " + covering.size() + " cells vs. max of " + coverer.maxCells());

    // Compute an S2Polygon representing the union of the cells in the covering.
    return S2Polygon.fromCellUnionBorder(covering);
  }

  // Checks that the distance between the given "shape" and the boundary of the S2Polygon
  // "coveringPolygon" is at least "distance".
  void testCoveringDistanceToGeometry(
      S2Polygon coveringPolygon, S1ChordAngle distance, S2Shape shape) {
    S2ShapeIndex coveringIndex = new S2ShapeIndex();
    coveringIndex.add(S2LaxPolygonShape.create(coveringPolygon));
    S2ClosestEdgeQuery.Query query =
        S2ClosestEdgeQuery.builder().setIncludeInteriors(false).build(coveringIndex);
    S2ShapeIndex targetIndex = new S2ShapeIndex();
    targetIndex.add(shape);
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2ClosestEdgeQuery.createShapeIndexTarget(targetIndex);
    assertFalse(query.isDistanceLess(target, distance));
  }


  // Each of the following tests verifies that a S2ShapeIndex containing some kind of shape is
  // buffered correctly, by first converting the covering to an S2Polygon and then checking that
  // (a) the S2Polygon contains the original geometry and (b) the distance between the original
  // geometry and the boundary of the S2Polygon is at least "radius".

  // Test buffering a set of points.
  public void testPointSet() {
    S2Point.Shape pointShape =
        S2Point.Shape.fromList(S2TextFormat.parsePointsOrDie("10:20, 10:23, 10:26"));
    S1ChordAngle distance = S1ChordAngle.fromDegrees(5);
    S2Polygon coveringPolygon = getBufferedIndexPolygon(pointShape, distance);

    // (a) Check that the covering contains the original points.
    for (S2Point p : pointShape) {
      assertTrue(coveringPolygon.contains(p));
    }

    // (b) Check that the distance between the boundary of the covering and the original points
    // is at least "distance".
    testCoveringDistanceToGeometry(coveringPolygon, distance, pointShape);
  }

  // Test buffering a polyline.
  public void testPolyline() {
    String polylineString = "10:5, 20:30, -10:60, -60:100";
    S2LaxPolylineShape lineShape = S2TextFormat.makeLaxPolylineOrDie(polylineString);
    S1ChordAngle distance = S1ChordAngle.fromDegrees(2);
    S2Polygon coveringPolygon = getBufferedIndexPolygon(lineShape, distance);

    // (a) Check that the covering contains the original polyline.
    assertTrue(
        coveringPolygon
            .subtractFromPolyline(S2TextFormat.makePolylineOrDie(polylineString))
            .isEmpty());

    // (b) Check that the distance between the boundary of the covering and the original polyline
    // is at least "distance".
    testCoveringDistanceToGeometry(coveringPolygon, distance, lineShape);
  }

  // Test buffering a polygon with a hole.
  public void testPolygonWithHole() {
    String polygonString = "10:10, 10:100, 70:0; 11:11, 69:0, 11:99";
    S2LaxPolygonShape laxPolygonShape = S2TextFormat.makeLaxPolygonOrDie(polygonString);
    S1ChordAngle distance = S1ChordAngle.fromDegrees(2);
    S2Polygon coveringPolygon = getBufferedIndexPolygon(laxPolygonShape, distance);

    // (a) Check that the covering contains the original polygon.
    assertTrue(coveringPolygon.contains(S2TextFormat.makePolygonOrDie(polygonString)));

    // (b) Check that the distance between the boundary of the covering and the original polygon
    // is at least "distance".
    testCoveringDistanceToGeometry(coveringPolygon, distance, laxPolygonShape);
  }

  // Test buffering a single point by 200 degrees.
  public void testHugeBufferRadius() {
    S2Point.Shape pointShape = S2Point.Shape.fromList(S2TextFormat.parsePointsOrDie("10:20"));
    S1ChordAngle distance = S1ChordAngle.fromDegrees(200);
    S2Polygon coveringPolygon = getBufferedIndexPolygon(pointShape, distance);

    // (a) Check that the covering contains the original point.
    for (S2Point p : pointShape) {
      assertTrue(coveringPolygon.contains(p));
    }

    // (b) Check that the distance between the boundary of the covering and the original point
    // is at least "distance". In this case, the covering must be the whole sphere and will have no
    // boundary. Therefore, the distance will be S1ChordAngle.INFINITY, which is not less than
    // S1ChordAngle.STRAIGHT, which is what 200 degrees will be mapped to as an S1Angle.
    assertTrue(coveringPolygon.isFull());
    testCoveringDistanceToGeometry(coveringPolygon, distance, pointShape);
  }
}
