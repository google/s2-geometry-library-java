/*
 * Copyright 2020 Google Inc.
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
import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.List;

/** Tests for S2TextFormat. */
@GwtCompatible
public strictfp class S2TextFormatTest extends GeometryTestCase {
  private static final int ITERATIONS = 10000;

  /**
   * Verify that S2TextFormat.ToString() formats the given lat/lng with at most "maxDigits" after
   * the decimal point and has no trailing zeros.
   */
  private static void expectMaxDigits(S2LatLng ll, int maxDigits) {
    String result = S2TextFormat.toString(ll);
    String[] values = result.split(":", 0);
    assertEquals(result, 2, values.length);
    for (String value : values) {
      int numDigits = 0;
      if (value.contains(".")) {
        numDigits = value.length() - value.indexOf(".") - 1;
        assertFalse(value, value.endsWith("0"));
      }
      assertTrue(
          "Expected at most "
              + maxDigits
              + " but got "
              + numDigits
              + ", with result "
              + value
              + " for input "
              + ll.toStringDegrees(),
          numDigits <= maxDigits);
    }
  }

  private static void expectString(String expected, S2LatLng ll) {
    assertEquals(expected, S2TextFormat.toString(ll));
  }

  private static void expectString(String expected, double d) {
    assertEquals(expected, Platform.formatDouble(d));
  }

  @GwtIncompatible("Platform.formatDouble() behavior is inconsistent in Javascript.")
  public void testDoubleToString() {
    expectString("1.01234567890123e-05", 1.012345678901234e-5);
    expectString("0.000101234567890123", 1.012345678901234e-4);
    expectString("0.00101234567890123", 1.012345678901234e-3);
    expectString("0.0101234567890123", 1.012345678901234e-2);
    expectString("0.101234567890123", 1.012345678901234e-1);
    expectString("1.01234567890123", 1.012345678901234);

    expectString("1e-05", 1.0e-5);
    expectString("1.01e-05", 1.01e-5);
    expectString("1.01234567890123e-05", 1.012345678901234e-5);
  }

  public void testToStringSpecialCases() {
    expectString("0:0", S2LatLng.fromDegrees(0, 0));
    expectString("90:0", new S2LatLng(new S2Point(0, 0, 1)));
    expectString("1e-20:1e-30", S2LatLng.fromDegrees(1e-20, 1e-30));
  }

  @GwtIncompatible("Fails on GWT, unclear why")
  public void testToStringNegativeZeros() {
    // Verify that negative zero coordinates in S2Points are formatted identically to positive
    // zeros.  This ensure that whenever two S2Points compare equal to each other, their string
    // representations do as well.
    //
    // <p>Note that we do not require that negative zero coordinates in S2LatLngs are formatted
    // identically to positive zeros, since this can result from legitimate differences between
    // S2Points.
    assertEquals("0:0", S2TextFormat.toString(new S2Point(1., -0., 0.)));
    assertEquals("0:0", S2TextFormat.toString(new S2Point(1., 0, -0.)));
    assertEquals("0:0", S2TextFormat.toString(new S2Point(1., -0., -0.)));
    assertEquals("0:180", S2TextFormat.toString(new S2Point(-1., -0., 0.)));
    assertEquals("0:180", S2TextFormat.toString(new S2Point(-1., 0., -0.)));
    assertEquals("0:180", S2TextFormat.toString(new S2Point(-1., -0., -0.)));
    assertEquals("90:0", S2TextFormat.toString(new S2Point(-0., 0., 1.)));
    assertEquals("90:0", S2TextFormat.toString(new S2Point(0., -0., 1.)));
    assertEquals("90:0", S2TextFormat.toString(new S2Point(-0., -0., 1.)));
  }

  @GwtIncompatible("Platform.formatDouble() behavior")
  public void testToStringMinimalDigitsE5() {
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      S2LatLng ll = new S2LatLng(randomPoint());
      S2LatLng llE5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5());
      expectMaxDigits(llE5, 5);
    }
  }

  @GwtIncompatible("Platform.formatDouble() behavior")
  public void testToStringMinimalDigitsE6() {
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      S2LatLng ll = new S2LatLng(randomPoint());
      S2LatLng llE6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6());
      expectMaxDigits(llE6, 6);
    }
  }

  @GwtIncompatible("Platform.formatDouble() behavior")
  public void testToStringMinimalDigitsE7() {
    expectMaxDigits(S2LatLng.fromDegrees(0, 0), 7);
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      S2LatLng ll = new S2LatLng(randomPoint());
      S2LatLng llE7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7());
      expectMaxDigits(llE7, 7);
    }
  }

  @GwtIncompatible("Platform.formatDouble() behavior")
  public void testToStringMinimalDigitsDoubleConstants() {
    // Verify that points specified as floating-point literals in degrees using
    // up to 10 digits after the decimal point are formatted with the minimal
    // number of digits.
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      int maxDigits = uniform(11);
      long scale = Math.round(Math.pow(10, maxDigits));
      long lat = Math.round(uniform(-90 * scale, 90 * scale));
      long lng = Math.round(uniform(-180 * scale, 180 * scale));
      S2LatLng ll = S2LatLng.fromDegrees(lat / (double) scale, lng / (double) scale);
      expectMaxDigits(ll, maxDigits);
    }
  }

  public void testToStringEmptyLoop() {
    assertEquals("empty", S2TextFormat.toString(S2Loop.empty()));
  }

  public void testToStringFullLoop() {
    assertEquals("full", S2TextFormat.toString(S2Loop.full()));
  }

  public void testToStringEmptyPolyline() {
    S2Polyline polyline = new S2Polyline(new ArrayList<S2Point>());
    assertEquals("", S2TextFormat.toString(polyline));
  }

  public void testToStringEmptyPointList() {
    List<S2Point> points = new ArrayList<>();
    assertEquals("", S2TextFormat.s2PointsToString(points));
  }

  public void testToStringEmptyPolygon() {
    S2Polygon empty = new S2Polygon();
    assertEquals("empty", S2TextFormat.toString(empty));
  }

  public void testToStringFullPolygon() {
    S2Polygon full = new S2Polygon(S2Loop.full());
    assertEquals("full", S2TextFormat.toString(full));
  }

  @GwtIncompatible("Platform.formatDouble() behavior is inconsistent in Javascript.")
  public void testToStringS2PolygonLoopSeparator() {
    String kLoop1 = "0:0, 0:5, 5:0";
    String kLoop2 = "1:1, 1:4, 4:1"; // Shells and holes same direction.
    S2Polygon polygon = S2TextFormat.makePolygonOrDie(kLoop1 + "; " + kLoop2);

    assertEquals(kLoop1 + ";\n" + kLoop2, S2TextFormat.toString(polygon));
    assertEquals(kLoop1 + "; " + kLoop2, S2TextFormat.toString(polygon, "; "));
  }

  @GwtIncompatible("S2LaxPolygonShape, Platform.formatDouble")
  public void testToStringLaxPolygonLoopSeparator() {
    String kLoop1 = "0:0, 0:5, 5:0";
    String kLoop2 = "1:1, 4:1, 1:4"; // Interior on left of all loops.
    S2LaxPolygonShape polygon = S2TextFormat.makeLaxPolygonOrDie(kLoop1 + "; " + kLoop2);
    assertEquals(kLoop1 + ";\n" + kLoop2, S2TextFormat.toString(polygon));
    assertEquals(kLoop1 + "; " + kLoop2, S2TextFormat.toString(polygon, "; "));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testMakeLaxPolygonEmpty() {
    // Verify that "" and "empty" both create empty polygons.
    S2LaxPolygonShape shape = S2TextFormat.makeLaxPolygonOrDie("");
    assertEquals(0, shape.numChains());
    shape = S2TextFormat.makeLaxPolygonOrDie("empty");
    assertEquals(0, shape.numChains());
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testMakeLaxPolygonFull() {
    S2LaxPolygonShape shape = S2TextFormat.makeLaxPolygonOrDie("full");
    assertEquals(1, shape.numChains());
    assertEquals(0, shape.getChainLength(0));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testMakeLaxPolygonFullWithHole() {
    S2LaxPolygonShape shape = S2TextFormat.makeLaxPolygonOrDie("full; 0:0");
    assertEquals(2, shape.numChains());
    assertEquals(0, shape.getChainLength(0));
    assertEquals(1, shape.getChainLength(1));
    assertEquals(1, shape.numEdges());
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  void testS2ShapeIndex(String str) {
    assertEquals(str, S2TextFormat.toString(S2TextFormat.makeIndexOrDie(str)));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testToStringS2ShapeIndex() {
    testS2ShapeIndex("# #");
    testS2ShapeIndex("0:0 # #");
    testS2ShapeIndex("0:0 | 1:0 # #");
    testS2ShapeIndex("0:0 | 1:0 # #");
    testS2ShapeIndex("# 0:0, 0:0 #");
    testS2ShapeIndex("# 0:0, 0:0 | 1:0, 2:0 #");
    testS2ShapeIndex("# # 0:0");
    testS2ShapeIndex("# # 0:0, 0:1");
    testS2ShapeIndex("# # 0:0, 0:1, 1:0");
    testS2ShapeIndex("# # 0:0, 0:1, 1:0; 2:2");
    testS2ShapeIndex("# # full");
  }

  public void testMakePointValidInput() {
    S2Point point = S2TextFormat.makePoint("-20:150");
    assertNotNull(point);
    assertEquals(S2LatLng.fromDegrees(-20, 150).toPoint(), point);
  }

  public void testMakePointInvalidInput() {
    assertNull(S2TextFormat.makePoint("blah"));
  }

  public void testSafeParseLatLngsValidInput() {
    List<S2LatLng> latlngs = S2TextFormat.parseLatLngs("-20:150, -20:151, -19:150");
    assertNotNull(latlngs);
    assertEquals(3, latlngs.size());
    assertEquals(latlngs.get(0), S2LatLng.fromDegrees(-20, 150));
    assertEquals(latlngs.get(1), S2LatLng.fromDegrees(-20, 151));
    assertEquals(latlngs.get(2), S2LatLng.fromDegrees(-19, 150));
  }

  public void testSafeParseLatLngsInvalidInput() {
    assertNull(S2TextFormat.parseLatLngs("blah"));
  }

  public void testSafeParsePointsValidInput() {
    List<S2Point> vertices = S2TextFormat.parsePoints("-20:150, -20:151, -19:150");
    assertNotNull(vertices);
    assertEquals(3, vertices.size());
    assertEquals(vertices.get(0), S2LatLng.fromDegrees(-20, 150).toPoint());
    assertEquals(vertices.get(1), S2LatLng.fromDegrees(-20, 151).toPoint());
    assertEquals(vertices.get(2), S2LatLng.fromDegrees(-19, 150).toPoint());
  }

  public void testSafeParsePointsInvalidInput() {
    List<S2Point> vertices = S2TextFormat.parsePoints("blah");
    assertNull(vertices);
  }

  public void testSafeMakeLatLngRectValidInput() {
    S2LatLngRect rect = S2TextFormat.makeLatLngRect("-10:-10, 10:10");
    assertNotNull(rect);
    assertEquals(
        rect, new S2LatLngRect(S2LatLng.fromDegrees(-10, -10), S2LatLng.fromDegrees(10, 10)));
  }

  public void testSafeMakeLatLngRectInvalidInput() {
    assertNull(S2TextFormat.makeLatLngRect("blah"));
  }

  public void testSafeMakeLatLngValidInput() {
    S2LatLng latlng = S2TextFormat.makeLatLng("-12.3:45.6");
    assertNotNull(latlng);
    assertEquals(latlng, S2LatLng.fromDegrees(-12.3, 45.6));
  }

  public void testSafeMakeLatLngInvalidInput() {
    assertNull(S2TextFormat.makeLatLng("blah"));
  }

  public void testSafeMakeCellIdValidInput() {
    S2CellId cellId = S2TextFormat.makeCellId("3/");
    assertNotNull(cellId);
    assertEquals(cellId, S2CellId.fromFace(3));
  }

  public void testSafeMakeCellIdInvalidInput() {
    assertNull(S2TextFormat.makeCellId("blah"));
    assertNull(S2TextFormat.makeCellId("6/0"));
    assertNull(S2TextFormat.makeCellId("3/04"));
  }

  public void testSafeMakeCellUnionValidInput() {
    S2CellUnion cellUnion = S2TextFormat.makeCellUnion("1/3, 4/");
    assertNotNull(cellUnion);
    S2CellUnion expected = new S2CellUnion();
    expected.initFromCellIds(
        new ArrayList<>(ImmutableList.of(S2CellId.fromFace(1).child(3), S2CellId.fromFace(4))));
    assertEquals(expected, cellUnion);
  }

  public void testSafeMakeCellUnionInvalidInput() {
    assertNull(S2TextFormat.makeCellUnion("abc"));
    assertNull(S2TextFormat.makeCellUnion("3/1 4/1"));
  }

  public void testSafeMakeLoopValidInput() {
    S2Loop loop = S2TextFormat.makeLoop("-20:150, -20:151, -19:150");
    assertNotNull(loop);
    assertTrue(
        loop.boundaryApproxEquals(
            new S2Loop(
                ImmutableList.of(
                    S2LatLng.fromDegrees(-20, 150).toPoint(),
                    S2LatLng.fromDegrees(-20, 151).toPoint(),
                    S2LatLng.fromDegrees(-19, 150).toPoint()))));
  }

  public void testSafeMakeLoopInvalidInput() {
    assertNull(S2TextFormat.makeLoop("blah"));
  }

  public void testSafeMakeLoopEmpty() {
    // Verify that "empty" creates an empty loop.
    S2Loop loop = S2TextFormat.makeLoop("empty");
    assertNotNull(loop);
    assertTrue(loop.isEmpty());
  }

  public void testSafeMakeLoopFull() {
    // Verify that "full" creates a full loop.
    S2Loop loop = S2TextFormat.makeLoop("full");
    assertNotNull(loop);
    assertTrue(loop.isFull());
  }

  public void testSafeMakePolylineValidInput() {
    S2Polyline polyline = S2TextFormat.makePolyline("-20:150, -20:151, -19:150");
    assertNotNull(polyline);
    S2Polyline expected =
        new S2Polyline(
            ImmutableList.of(
                S2LatLng.fromDegrees(-20, 150).toPoint(),
                S2LatLng.fromDegrees(-20, 151).toPoint(),
                S2LatLng.fromDegrees(-19, 150).toPoint()));
    assertEquals(expected, polyline);
  }

  public void testSafeMakePolylineInvalidInput() {
    assertNull(S2TextFormat.makePolyline("blah"));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testSafeMakeLaxPolylineValidInput() {
    S2LaxPolylineShape laxPolyline = S2TextFormat.makeLaxPolyline("-20:150, -20:151, -19:150");
    assertNotNull(laxPolyline);
    // No easy equality check for LaxPolylines; check vertices instead.
    assertEquals(3, laxPolyline.numVertices());
    assertTrue(new S2LatLng(laxPolyline.vertex(0)).approxEquals(S2LatLng.fromDegrees(-20, 150)));
    assertTrue(new S2LatLng(laxPolyline.vertex(1)).approxEquals(S2LatLng.fromDegrees(-20, 151)));
    assertTrue(new S2LatLng(laxPolyline.vertex(2)).approxEquals(S2LatLng.fromDegrees(-19, 150)));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testSafeMakeLaxPolylineInvalidInput() {
    assertNull(S2TextFormat.makeLaxPolyline("blah"));
  }

  public void testSafeMakePolygonValidInput() {
    S2Polygon polygon = S2TextFormat.makePolygon("-20:150, -20:151, -19:150");
    assertNotNull(polygon);
    List<S2Point> vertices =
        ImmutableList.of(
            S2LatLng.fromDegrees(-20, 150).toPoint(),
            S2LatLng.fromDegrees(-20, 151).toPoint(),
            S2LatLng.fromDegrees(-19, 150).toPoint());
    S2Polygon expected = new S2Polygon(new S2Loop(vertices));
    assertEquals(expected, polygon);
  }

  public void testSafeMakePolygonInvalidInput() {
    assertNull(S2TextFormat.makePolygon("blah"));
  }

  public void testSafeMakePolygonEmpty() {
    // Verify that "" and "empty" both create empty polygons.
    S2Polygon polygon = S2TextFormat.makePolygon("");
    assertNotNull(polygon);
    assertTrue(polygon.isEmpty());
    polygon = S2TextFormat.makePolygon("empty");
    assertNotNull(polygon);
    assertTrue(polygon.isEmpty());
  }

  public void testSafeMakePolygonFull() {
    // Verify that "full" creates the full polygon.
    S2Polygon polygon = S2TextFormat.makePolygon("full");
    assertNotNull(polygon);
    assertTrue(polygon.isFull());
  }

  public void testSafeMakeVerbatimPolygonValidInput() {
    S2Polygon polygon = S2TextFormat.makeVerbatimPolygon("-20:150, -20:151, -19:150");
    List<S2Point> vertices =
        ImmutableList.of(
            S2LatLng.fromDegrees(-20, 150).toPoint(),
            S2LatLng.fromDegrees(-20, 151).toPoint(),
            S2LatLng.fromDegrees(-19, 150).toPoint());
    S2Polygon expected = new S2Polygon(new S2Loop(vertices));
    assertEquals(expected, polygon);
  }

  public void testSafeMakeVerbatimPolygonInvalidInput() {
    assertNull(S2TextFormat.makeVerbatimPolygon("blah"));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testSafeMakeLaxPolygonValidInput() {
    S2LaxPolygonShape laxPolygon = S2TextFormat.makeLaxPolygon("-20:150, -20:151, -19:150");
    assertNotNull(laxPolygon);

    // No easy equality check for LaxPolygons; check vertices, edges, and chains instead.
    // Expect three edges joining three vertices in one chain.
    assertEquals(1, laxPolygon.numChains());
    assertEquals(3, laxPolygon.numVertices());
    assertEquals(3, laxPolygon.numEdges());

    // Check each chain vertex.
    S2Point[] expectedVertices =
        new S2Point[] {
          S2LatLng.fromDegrees(-20, 150).toPoint(),
          S2LatLng.fromDegrees(-20, 151).toPoint(),
          S2LatLng.fromDegrees(-19, 150).toPoint()
        };

    assertTrue(
        new S2LatLng(laxPolygon.getChainVertex(0, 0))
            .approxEquals(new S2LatLng(expectedVertices[0])));
    assertTrue(
        new S2LatLng(laxPolygon.getChainVertex(0, 1))
            .approxEquals(new S2LatLng(expectedVertices[1])));
    assertTrue(
        new S2LatLng(laxPolygon.getChainVertex(0, 2))
            .approxEquals(new S2LatLng(expectedVertices[2])));

    // Ensure the chain has three edges, and starts at offset 0.
    assertEquals(3, laxPolygon.getChainLength(0));
    assertEquals(0, laxPolygon.getChainStart(0));

    // The chain edges connect the vertices in a closed loop.
    S2Point[][] expectedEdges =
        new S2Point[][] {
          {expectedVertices[0], expectedVertices[1]},
          {expectedVertices[1], expectedVertices[2]},
          {expectedVertices[2], expectedVertices[0]}
        };

    MutableEdge chainEdge = new MutableEdge();
    laxPolygon.getChainEdge(0, 0, chainEdge);
    assertEquals(expectedEdges[0][0], chainEdge.getStart());
    assertEquals(expectedEdges[0][1], chainEdge.getEnd());
    laxPolygon.getChainEdge(0, 1, chainEdge);
    assertEquals(expectedEdges[1][0], chainEdge.getStart());
    assertEquals(expectedEdges[1][1], chainEdge.getEnd());
    laxPolygon.getChainEdge(0, 2, chainEdge);
    assertEquals(expectedEdges[2][0], chainEdge.getStart());
    assertEquals(expectedEdges[2][1], chainEdge.getEnd());
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testSafeMakeLaxPolygonInvalidInput() {
    assertNull(S2TextFormat.makeLaxPolygon("blah"));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testSafeMakeIndexValidInput() {
    S2ShapeIndex index = S2TextFormat.makeIndex("# 0:0, 0:0 | 1:0, 2:0 #");
    assertNotNull(index);
    assertEquals("# 0:0, 0:0 | 1:0, 2:0 #", S2TextFormat.toString(index));
  }

  @GwtIncompatible("S2LaxPolylineShape, S2LaxPolygonShape")
  public void testSafeMakeIndexInvalidInput() {
    assertNull(S2TextFormat.makeIndex("# blah #"));
  }
}
