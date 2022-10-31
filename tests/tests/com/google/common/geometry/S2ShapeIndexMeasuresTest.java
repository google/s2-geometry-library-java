/*
 * Copyright 2019 Google Inc.
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

import static com.google.common.geometry.TestDataGenerator.makeIndex;
import static java.lang.Math.PI;

import com.google.common.collect.Lists;

public class S2ShapeIndexMeasuresTest extends GeometryTestCase {

  public void testDimensionEmpty() {
    assertEquals(-1, S2ShapeIndexMeasures.dimension(makeIndex("# #")));
  }

  public void testDimensionPoints() {
    assertEquals(0, S2ShapeIndexMeasures.dimension(makeIndex("0:0 # #")));

    // Create an index with an empty point set.
    S2ShapeIndex shapeIndex = new S2ShapeIndex();
    shapeIndex.add(S2Point.Shape.fromList(Lists.newArrayList()));
    assertEquals(0, S2ShapeIndexMeasures.dimension(shapeIndex));
  }

  public void testDimensionPointsAndLines() {
    assertEquals(1, S2ShapeIndexMeasures.dimension(makeIndex("0:0 # 1:1, 1:2 #")));

    // Note that a polyline with one vertex has no edges, so it is effectively empty for the purpose
    // of testing dimension().
    assertEquals(1, S2ShapeIndexMeasures.dimension(makeIndex("0:0 # 1:1 #")));
  }

  public void testDimensionPointsLinesAndPolygons() {
    assertEquals(2, S2ShapeIndexMeasures.dimension(makeIndex("0:0 # 1:1, 2:2 # 3:3, 3:4, 4:3")));

    assertEquals(2, S2ShapeIndexMeasures.dimension(makeIndex("# # empty")));
  }

  public void testLengthEmpty() {
    assertEquals(S1Angle.ZERO, S2ShapeIndexMeasures.length(makeIndex("# #")));
  }

  public void testLengthTwoLines() {
    assertEquals(
        S1Angle.degrees(2),
        S2ShapeIndexMeasures.length(makeIndex("4:4 # 0:0, 1:0 | 1:0, 2:0 # 5:5, 5:6, 6:5")));
  }

  public void testPerimeterEmpty() {
    assertEquals(S1Angle.ZERO, S2ShapeIndexMeasures.perimeter(makeIndex("# #")));
  }

  public void testPerimeterDegeneratePolygon() {
    assertDoubleNear(
        4,
        S2ShapeIndexMeasures.perimeter(makeIndex("4:4 # 0:0, 1:0 | 2:0, 3:0 # 0:1, 0:2, 0:3"))
            .degrees(),
        1e-12);
  }

  public void testAreaEmpty() {
    assertEquals(0.0, S2ShapeIndexMeasures.area(makeIndex("# #")), 0.0);
  }

  public void testAreaTwoFullPolygons() {
    assertDoubleEquals(8 * PI, S2ShapeIndexMeasures.area(makeIndex("# # full | full")));
  }

  public void testCentroidPoints() {
    S2Point actual = S2ShapeIndexMeasures.centroid(makeIndex("0:0 | 0:90 # #"));
    assertEquals(new S2Point(1, 1, 0), actual);
  }

  public void testCentroidPolyline() {
    // Checks that points are ignored when computing the centroid.
    S2Point actual = S2ShapeIndexMeasures.centroid(makeIndex("5:5 | 6:6 # 0:0, 0:90 #"));
    assertTrue(S2.approxEquals(new S2Point(1, 1, 0), actual));
  }

  public void testCentroidPolygon() {
    // Checks that points and polylines are ignored when computing the centroid.
    S2Point actual = S2ShapeIndexMeasures.centroid(makeIndex("5:5 # 6:6, 7:7 # 0:0, 0:90, 90:0"));
    assertTrue(S2.approxEquals(new S2Point(PI / 4, PI / 4, PI / 4), actual));
  }
}
