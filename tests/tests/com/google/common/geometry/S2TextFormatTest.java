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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for S2TextFormat. */
@RunWith(JUnit4.class)
public class S2TextFormatTest extends GeometryTestCase {
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

  // Under J2CL, this fails with expected:<[1.01234567890123e-05]> but was:<[0.0000101234567]>
  @GwtIncompatible("Platform.formatDouble() behavior is inconsistent in Javascript.")
  @Test
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

  @Test
  public void testToStringSpecialCases() {
    expectString("0:0", S2LatLng.fromDegrees(0, 0));
    expectString("90:0", new S2LatLng(new S2Point(0, 0, 1)));
    expectString("1e-20:1e-30", S2LatLng.fromDegrees(1e-20, 1e-30));
  }

  @Test
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

  /** Simple tests that lat/lng values are printed using a minimal number of digits. */
  @Test
  public void testToStringMinimalDigits() {
    assertEquals(
        "12.3456789:-98.7654321",
        S2TextFormat.toString(S2LatLng.fromDegrees(12.3456789, -98.7654321)));
    assertEquals("1.23456:0.98765", S2TextFormat.toString(S2LatLng.fromE5(123456, 98765)));
    assertEquals("1:2", S2TextFormat.toString(S2LatLng.fromDegrees(1, 2)));
    assertEquals("3:4", S2TextFormat.toString(S2LatLng.fromDegrees(3, 4)));
    assertEquals("0:1.5", S2TextFormat.toString(S2LatLng.fromDegrees(0, 1.5)));
    assertEquals("88:160", S2TextFormat.toString(S2LatLng.fromDegrees(88, 160)));
    assertEquals("13.45:-12.34", S2TextFormat.toString(S2LatLng.fromDegrees(13.45, -12.34)));
  }

  @Test
  public void testToStringMinimalDigitsE5() {
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      S2LatLng ll = new S2LatLng(data.getRandomPoint());
      S2LatLng llE5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5());
      expectMaxDigits(llE5, 5);
    }
  }

  @Test
  public void testToStringMinimalDigitsE6() {
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      S2LatLng ll = new S2LatLng(data.getRandomPoint());
      S2LatLng llE6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6());
      expectMaxDigits(llE6, 6);
    }
  }

  @Test
  public void testToStringMinimalDigitsE7() {
    expectMaxDigits(S2LatLng.fromDegrees(0, 0), 7);
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      S2LatLng ll = new S2LatLng(data.getRandomPoint());
      S2LatLng llE7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7());
      expectMaxDigits(llE7, 7);
    }
  }

  @Test
  public void testToStringMinimalDigitsDoubleConstants() {
    // Verify that points specified as floating-point literals in degrees using up to 10 digits
    // after the decimal point are formatted with the minimal number of digits.
    for (int iter = 0; iter < ITERATIONS; ++iter) {
      int maxDigits = data.uniform(11);
      double scale = Math.round(Math.pow(10, maxDigits));
      long lat = Math.round(data.uniform(-90 * scale, 90 * scale));
      long lng = Math.round(data.uniform(-180 * scale, 180 * scale));
      S2LatLng ll = S2LatLng.fromDegrees(lat / scale, lng / scale);
      expectMaxDigits(ll, maxDigits);
    }
  }

  @Test
  public void testToStringEmptyLoop() {
    assertEquals("empty", S2TextFormat.toString(S2Loop.empty()));
  }

  @Test
  public void testToStringFullLoop() {
    assertEquals("full", S2TextFormat.toString(S2Loop.full()));
  }

  @Test
  public void testToStringEmptyPolyline() {
    S2Polyline polyline = new S2Polyline(new ArrayList<S2Point>());
    assertEquals("", S2TextFormat.toString(polyline));
  }

  @Test
  public void testToStringEmptyPointList() {
    List<S2Point> points = new ArrayList<>();
    assertEquals("", S2TextFormat.s2PointsToString(points));
  }

  @Test
  public void testToStringEmptyPolygon() {
    S2Polygon empty = new S2Polygon();
    assertEquals("empty", S2TextFormat.toString(empty));
  }

  @Test
  public void testToStringFullPolygon() {
    S2Polygon full = new S2Polygon(S2Loop.full());
    assertEquals("full", S2TextFormat.toString(full));
  }

  @Test
  public void testToStringS2PolygonLoopSeparator() {
    String kLoop1 = "0:0, 0:5, 5:0";
    String kLoop2 = "1:1, 1:4, 4:1"; // Shells and holes same direction.
    S2Polygon polygon = S2TextFormat.makePolygonOrDie(kLoop1 + "; " + kLoop2);

    assertEquals(kLoop1 + ";\n" + kLoop2, S2TextFormat.toString(polygon));
    assertEquals(kLoop1 + "; " + kLoop2, S2TextFormat.toString(polygon, "; "));
  }

  @Test
  public void testToStringLaxPolygonLoopSeparator() {
    String kLoop1 = "0:0, 0:5, 5:0";
    String kLoop2 = "1:1, 4:1, 1:4"; // Interior on left of all loops.
    S2LaxPolygonShape polygon = S2TextFormat.makeLaxPolygonOrDie(kLoop1 + "; " + kLoop2);
    assertEquals(kLoop1 + ";\n" + kLoop2, S2TextFormat.toString(polygon));
    assertEquals(kLoop1 + "; " + kLoop2, S2TextFormat.toString(polygon, "; "));
  }

  @Test
  public void testMakeLaxPolygonEmpty() {
    // Verify that "" and "empty" both create empty polygons.
    S2LaxPolygonShape shape = S2TextFormat.makeLaxPolygonOrDie("");
    assertEquals(0, shape.numChains());
    shape = S2TextFormat.makeLaxPolygonOrDie("empty");
    assertEquals(0, shape.numChains());
  }

  @Test
  public void testMakeLaxPolygonFull() {
    S2LaxPolygonShape shape = S2TextFormat.makeLaxPolygonOrDie("full");
    assertEquals(1, shape.numChains());
    assertEquals(0, shape.getChainLength(0));
  }

  @Test
  public void testMakeLaxPolygonFullWithHole() {
    S2LaxPolygonShape shape = S2TextFormat.makeLaxPolygonOrDie("full; 0:0");
    assertEquals(2, shape.numChains());
    assertEquals(0, shape.getChainLength(0));
    assertEquals(1, shape.getChainLength(1));
    assertEquals(1, shape.numEdges());
  }

  void testS2ShapeIndex(String str) {
    assertEquals(str, S2TextFormat.toString(S2TextFormat.makeIndexOrDie(str)));
  }

  @Test
  public void testToStringS2Point() {
    assertEquals("0:0", S2TextFormat.toString(S2Point.X_POS));
    assertEquals("0:90", S2TextFormat.toString(S2Point.Y_POS));
    assertEquals("90:0", S2TextFormat.toString(S2Point.Z_POS));
    assertEquals("null", S2TextFormat.toString((S2Point) null));
  }

  @Test
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

  @Test
  public void testMakePointValidInput() {
    S2Point point = S2TextFormat.makePoint("-20:150");
    assertNotNull(point);
    assertEquals(S2LatLng.fromDegrees(-20, 150).toPoint(), point);
  }

  @Test
  public void testMakePointInvalidInput() {
    assertNull(S2TextFormat.makePoint("blah"));
  }

  @Test
  public void testSafeParseLatLngsValidInput() {
    List<S2LatLng> latlngs = S2TextFormat.parseLatLngs("-20:150, -20:151, -19:150");
    assertNotNull(latlngs);
    assertEquals(3, latlngs.size());
    assertEquals(latlngs.get(0), S2LatLng.fromDegrees(-20, 150));
    assertEquals(latlngs.get(1), S2LatLng.fromDegrees(-20, 151));
    assertEquals(latlngs.get(2), S2LatLng.fromDegrees(-19, 150));
  }

  @Test
  public void testSafeParseLatLngsInvalidInput() {
    assertNull(S2TextFormat.parseLatLngs("blah"));
  }

  @Test
  public void testSafeParsePointsValidInput() {
    List<S2Point> vertices = S2TextFormat.parsePoints("-20:150, -20:151, -19:150");
    assertNotNull(vertices);
    assertEquals(3, vertices.size());
    assertEquals(vertices.get(0), S2LatLng.fromDegrees(-20, 150).toPoint());
    assertEquals(vertices.get(1), S2LatLng.fromDegrees(-20, 151).toPoint());
    assertEquals(vertices.get(2), S2LatLng.fromDegrees(-19, 150).toPoint());
  }

  @Test
  public void testSafeParsePointsInvalidInput() {
    List<S2Point> vertices = S2TextFormat.parsePoints("blah");
    assertNull(vertices);
  }

  @Test
  public void testSafeMakeLatLngRectValidInput() {
    S2LatLngRect rect = S2TextFormat.makeLatLngRect("-10:-10, 10:10");
    assertNotNull(rect);
    assertEquals(
        rect, new S2LatLngRect(S2LatLng.fromDegrees(-10, -10), S2LatLng.fromDegrees(10, 10)));
  }

  @Test
  public void testSafeMakeLatLngRectInvalidInput() {
    assertNull(S2TextFormat.makeLatLngRect("blah"));
  }

  @Test
  public void testSafeMakeLatLngValidInput() {
    S2LatLng latlng = S2TextFormat.makeLatLng("-12.3:45.6");
    assertNotNull(latlng);
    assertEquals(latlng, S2LatLng.fromDegrees(-12.3, 45.6));
  }

  @Test
  public void testSafeMakeLatLngInvalidInput() {
    assertNull(S2TextFormat.makeLatLng("blah"));
  }

  @Test
  public void testSafeMakeCellIdValidInput() {
    S2CellId cellId = S2TextFormat.makeCellId("3/");
    assertNotNull(cellId);
    assertEquals(cellId, S2CellId.fromFace(3));
  }

  @Test
  public void testSafeMakeCellIdInvalidInput() {
    assertNull(S2TextFormat.makeCellId("blah"));
    assertNull(S2TextFormat.makeCellId("6/0"));
    assertNull(S2TextFormat.makeCellId("3/04"));
  }

  @Test
  public void testSafeMakeCellUnionValidInput() {
    S2CellUnion cellUnion = S2TextFormat.makeCellUnion("1/3, 4/");
    assertNotNull(cellUnion);
    S2CellUnion expected = new S2CellUnion();
    expected.initFromCellIds(
        new ArrayList<>(ImmutableList.of(S2CellId.fromFace(1).child(3), S2CellId.fromFace(4))));
    assertEquals(expected, cellUnion);
  }

  @Test
  public void testSafeMakeCellUnionInvalidInput() {
    assertNull(S2TextFormat.makeCellUnion("abc"));
    assertNull(S2TextFormat.makeCellUnion("3/1 4/1"));
  }

  @Test
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

  @Test
  public void testSafeMakeLoopInvalidInput() {
    assertNull(S2TextFormat.makeLoop("blah"));
  }

  @Test
  public void testSafeMakeLoopEmpty() {
    // Verify that "empty" creates an empty loop.
    S2Loop loop = S2TextFormat.makeLoop("empty");
    assertNotNull(loop);
    assertTrue(loop.isEmpty());
  }

  @Test
  public void testSafeMakeLoopFull() {
    // Verify that "full" creates a full loop.
    S2Loop loop = S2TextFormat.makeLoop("full");
    assertNotNull(loop);
    assertTrue(loop.isFull());
  }

  @Test
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

  @Test
  public void testSafeMakePolylineInvalidInput() {
    assertNull(S2TextFormat.makePolyline("blah"));
  }

  @Test
  public void testSafeMakeLaxPolylineValidInput() {
    S2LaxPolylineShape laxPolyline = S2TextFormat.makeLaxPolyline("-20:150, -20:151, -19:150");
    assertNotNull(laxPolyline);
    // No easy equality check for LaxPolylines; check vertices instead.
    assertEquals(3, laxPolyline.numVertices());
    assertTrue(new S2LatLng(laxPolyline.vertex(0)).approxEquals(S2LatLng.fromDegrees(-20, 150)));
    assertTrue(new S2LatLng(laxPolyline.vertex(1)).approxEquals(S2LatLng.fromDegrees(-20, 151)));
    assertTrue(new S2LatLng(laxPolyline.vertex(2)).approxEquals(S2LatLng.fromDegrees(-19, 150)));
  }

  @Test
  public void testSafeMakeLaxPolylineInvalidInput() {
    assertNull(S2TextFormat.makeLaxPolyline("blah"));
  }

  @Test
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

  @Test
  public void testSafeMakePolygonInvalidInput() {
    assertNull(S2TextFormat.makePolygon("blah"));
  }

  @Test
  public void testSafeMakePolygonEmpty() {
    // Verify that "" and "empty" both create empty polygons.
    S2Polygon polygon = S2TextFormat.makePolygon("");
    assertNotNull(polygon);
    assertTrue(polygon.isEmpty());
    polygon = S2TextFormat.makePolygon("empty");
    assertNotNull(polygon);
    assertTrue(polygon.isEmpty());
  }

  @Test
  public void testSafeMakePolygonFull() {
    // Verify that "full" creates the full polygon.
    S2Polygon polygon = S2TextFormat.makePolygon("full");
    assertNotNull(polygon);
    assertTrue(polygon.isFull());
  }

  @Test
  public void testSafeMakeVerbatimPolygonValidInput() {
    S2Polygon polygon = makeVerbatimPolygon("-20:150, -20:151, -19:150");
    List<S2Point> vertices =
        ImmutableList.of(
            S2LatLng.fromDegrees(-20, 150).toPoint(),
            S2LatLng.fromDegrees(-20, 151).toPoint(),
            S2LatLng.fromDegrees(-19, 150).toPoint());
    S2Polygon expected = new S2Polygon(new S2Loop(vertices));
    assertEquals(expected, polygon);
  }

  @Test
  public void testSafeMakeVerbatimPolygonInvalidInput() {
    assertNull(makeVerbatimPolygon("blah"));
  }

  @Test
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

  @Test
  public void testSafeMakeLaxPolygonInvalidInput() {
    assertNull(S2TextFormat.makeLaxPolygon("blah"));
  }

  @Test
  public void testSafeMakeIndexValidInput() {
    S2ShapeIndex index = S2TextFormat.makeIndex("# 0:0, 0:0 | 1:0, 2:0 #");
    assertNotNull(index);
    assertEquals("# 0:0, 0:0 | 1:0, 2:0 #", S2TextFormat.toString(index));
  }

  @Test
  public void testSafeMakeIndexInvalidInput() {
    assertNull(S2TextFormat.makeIndex("# blah #"));
  }
}
