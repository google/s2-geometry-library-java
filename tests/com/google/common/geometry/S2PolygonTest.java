/*
 * Copyright 2006 Google Inc.
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
import com.google.common.collect.Maps;

import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * Tests for {@link S2Polygon}.
 *
 */
@GwtCompatible
public strictfp class S2PolygonTest extends GeometryTestCase {
  private static Logger logger = Logger.getLogger(S2PolygonTest.class.getName());

  /**
   * The error margin used when comparing the actual and expected results of
   * constructive geometry operations.
   * The intersections in the expected data were computed in lat-lng
   * space, while the actual intersections are computed using geodesics.
   * The error due to this depends on the length and direction of the line
   * segment being intersected, and how close the intersection is to the
   * endpoints of the segment. The worst case is for a line segment between
   * two points at the same latitude, where the intersection point is in the
   * middle of the segment.  In this case the error is approximately
   * {@code (p * t^2) / 8}, where {@code p} is the absolute latitude in
   * radians, {@code t} is the longitude difference in radians, and both
   * {@code p} and {@code t} are small. The test cases all have small latitude
   * and longitude differences. If {@code p} and {@code t} are converted to
   * degrees, the following error bound is valid as long as
   * {@code (p * t^2 < 150)}.
   */
  private static final double OPERATIONS_MAX_ERROR = 1e-4;

  // A set of nested loops around the point 0:0 (lat:lng).
  // Every vertex of NEAR0 is a vertex of NEAR1.
  private static final String NEAR_POINT = "0:0";
  private static final String NEAR0 = "-1:0, 0:1, 1:0, 0:-1;";
  private static final String NEAR1 = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;";
  private static final String NEAR2 = "-1:-2, -2:5, 5:-2;";
  private static final String NEAR3 = "-2:-2, -3:6, 6:-3;";
  private static final String NEAR_HEMI = "0:-90, -90:0, 0:90, 90:0;";

  // A set of nested loops around the point 0:180 (lat:lng).
  // Every vertex of FAR0 and FAR2 belongs to FAR1, and all
  // the loops except FAR2 are non-convex.
  private static final String FAR0 = "0:179, 1:180, 0:-179, 2:-180;";
  private static final String FAR1 =
      "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;";
  private static final String FAR2 = "3:-178, 3:178, -1:179, -1:-179;";
  private static final String FAR3 = "-3:-178, 4:-177, 4:177, -3:178, -2:179;";
  private static final String FAR_HEMI = "0:-90, 60:90, -60:90;";

  // A set of nested loops around the point -90:0 (lat:lng).
  private static final String SOUTH_POINT = "-89.9999:0.001";
  private static final String SOUTH0A = "-90:0, -89.99:0.01, -89.99:0;";
  private static final String SOUTH0B = "-90:0, -89.99:0.03, -89.99:0.02;";
  private static final String SOUTH0C = "-90:0, -89.99:0.05, -89.99:0.04;";
  private static final String SOUTH1 = "-90:0, -89.9:0.1, -89.9:-0.1;";
  private static final String SOUTH2 = "-90:0, -89.8:0.2, -89.8:-0.2;";
  private static final String SOUTH_HEMI = "0:-180, 0:60, 0:-60;";

  // Two different loops that surround all the NEAR and FAR loops except
  // for the hemispheres.
  private static final String NEAR_FAR1 = "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, "
      + "1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;";
  private static final String NEAR_FAR2 =
      "-2:15, -2:170, -8:-175, 8:-175, 2:170, 2:15, 8:-4, -8:-4;";

  // Loops that result from intersection of other loops.
  private static final String FAR_SOUTH_H = "0:-180, 0:90, -60:90, 0:-90;";

  // Rectangles that form a cross, with only shared vertices, no crossing edges.
  // Optional holes outside the intersecting region.
  private static final String CROSS1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1;";
  private static final String CROSS1_SIDE_HOLE = "-1.5:0.5, -1.2:0.5, -1.2:-0.5, -1.5:-0.5;";
  private static final String CROSS2 = "1:-2, 1:-1, 1:1, 1:2, -1:2, -1:1, -1:-1, -1:-2;";
  private static final String CROSS2_SIDE_HOLE = "0.5:-1.5, 0.5:-1.2, -0.5:-1.2, -0.5:-1.5;";
  private static final String CROSS_CENTER_HOLE = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

  //Two rectangles that intersect, but no edges cross and there's always
  //local containment (rather than crossing) at each shared vertex.
  //In this ugly ASCII art, 1 is A+B, 2 is B+C:
  //     +---+---+---+
  //     | A | B | C |
  //     +---+---+---+
  private static final String OVERLAP1 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
  private static final String OVERLAP1_SIDE_HOLE = "0.2:0.8, 0.8:0.8, 0.8:0.2, 0.2:0.2;";
  private static final String OVERLAP2 = "1:1, 2:1, 3:1, 3:0, 2:0, 1:0;";
  private static final String OVERLAP2_SIDE_HOLE = "2.2:0.8, 2.8:0.8, 2.8:0.2, 2.2:0.2;";
  private static final String OVERLAP_CENTER_HOLE = "1.2:0.8, 1.8:0.8, 1.8:0.2, 1.2:0.2;";

  // Two rectangles that are "adjacent", but rather than having common edges,
  // those edges are slighly off. A third rectangle that is not adjacent to
  // either of the first two.
  private static final String ADJACENT0 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
  private static final String ADJACENT1 = "0:2, 1:2, 2:2, 2:1.01, 1:0.99, 0:1.01;";
  private static final String UN_ADJACENT = "10:10, 11:10, 12:10, 12:9, 11:9, 10:9;";

  // Shapes used to test comparison functions for polygons.
  private static final String RECTANGLE1 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
  private static final String RECTANGLE2 = "5:1, 6:1, 7:1, 7:0, 6:0, 5:0;";
  private static final String TRIANGLE = "15:0, 17:0, 16:2;";
  private static final String TRIANGLE_ROT = "17:0, 16:2, 15:0;";

  private final S2Polygon near10 = makeVerbatimPolygon(NEAR0 + NEAR1);
  private final S2Polygon near30 = makeVerbatimPolygon(NEAR3 + NEAR0);
  private final S2Polygon near32 = makeVerbatimPolygon(NEAR2 + NEAR3);
  private final S2Polygon near3210 = makeVerbatimPolygon(NEAR0 + NEAR2 + NEAR3 + NEAR1);
  private final S2Polygon nearH3210 = makeVerbatimPolygon(
      NEAR0 + NEAR2 + NEAR3 + NEAR_HEMI + NEAR1);

  private final S2Polygon far10 = makeVerbatimPolygon(FAR0 + FAR1);
  private final S2Polygon far21 = makeVerbatimPolygon(FAR2 + FAR1);
  private final S2Polygon far321 = makeVerbatimPolygon(FAR2 + FAR3 + FAR1);
  private final S2Polygon farH20 = makeVerbatimPolygon(FAR2 + FAR_HEMI + FAR0);
  private final S2Polygon farH3210 = makeVerbatimPolygon(FAR2 + FAR_HEMI + FAR0 + FAR1 + FAR3);

  private final S2Polygon south0ab = makeVerbatimPolygon(SOUTH0A + SOUTH0B);
  private final S2Polygon south2 = makeVerbatimPolygon(SOUTH2);
  private final S2Polygon south210b = makeVerbatimPolygon(SOUTH2 + SOUTH0B + SOUTH1);
  private final S2Polygon southH21 = makeVerbatimPolygon(SOUTH2 + SOUTH_HEMI + SOUTH1);
  private final S2Polygon southH20abc = makeVerbatimPolygon(
      SOUTH2 + SOUTH0B + SOUTH_HEMI + SOUTH0A + SOUTH0C);

  private final S2Polygon nf1n10f2s10abc = makeVerbatimPolygon(
      SOUTH0C + FAR2 + NEAR1 + NEAR_FAR1 + NEAR0 + SOUTH1 + SOUTH0B + SOUTH0A);

  private final S2Polygon nf2n2f210s210ab = makeVerbatimPolygon(
      FAR2 + SOUTH0A + FAR1 + SOUTH1 + FAR0 + SOUTH0B + NEAR_FAR2 + SOUTH2 + NEAR2);

  private final S2Polygon f32n0 = makeVerbatimPolygon(FAR2 + NEAR0 + FAR3);
  private final S2Polygon n32s0b = makeVerbatimPolygon(NEAR3 + SOUTH0B + NEAR2);

  private final S2Polygon adj0 = makeVerbatimPolygon(ADJACENT0);
  private final S2Polygon adj1 = makeVerbatimPolygon(ADJACENT1);
  private final S2Polygon unAdj = makeVerbatimPolygon(UN_ADJACENT);

  private final S2Polygon farH = makeVerbatimPolygon(FAR_HEMI);
  private final S2Polygon southH = makeVerbatimPolygon(SOUTH_HEMI);
  private final S2Polygon farHSouthH = makeVerbatimPolygon(FAR_SOUTH_H);

  private final S2Polygon cross1 = makeVerbatimPolygon(CROSS1);
  private final S2Polygon cross1SideHole = makeVerbatimPolygon(CROSS1 + CROSS1_SIDE_HOLE);
  private final S2Polygon cross1CenterHole = makeVerbatimPolygon(CROSS1 + CROSS_CENTER_HOLE);
  private final S2Polygon cross2 = makeVerbatimPolygon(CROSS2);
  private final S2Polygon cross2SideHole = makeVerbatimPolygon(CROSS2 + CROSS2_SIDE_HOLE);
  private final S2Polygon cross2CenterHole = makeVerbatimPolygon(CROSS2 + CROSS_CENTER_HOLE);

  private final S2Polygon overlap1 = makeVerbatimPolygon(OVERLAP1);
  private final S2Polygon overlap1SideHole = makeVerbatimPolygon(OVERLAP1 + OVERLAP1_SIDE_HOLE);
  private final S2Polygon overlap2 = makeVerbatimPolygon(OVERLAP2);
  private final S2Polygon overlap2SideHole = makeVerbatimPolygon(OVERLAP2 + OVERLAP2_SIDE_HOLE);
  private final S2Polygon overlap1CenterHole = makeVerbatimPolygon(OVERLAP1 + OVERLAP_CENTER_HOLE);
  private final S2Polygon overlap2CenterHole = makeVerbatimPolygon(OVERLAP2 + OVERLAP_CENTER_HOLE);

  private final S2Polygon empty = new S2Polygon();
  private final S2Polygon full = new S2Polygon(S2Loop.full());

  private static void checkContains(String aText, String bText) {
    S2Polygon a = makeVerbatimPolygon(aText);
    S2Polygon b = makeVerbatimPolygon(bText);
    assertTrue(a.contains(b));
    assertTrue(a.approxContains(b, S1Angle.radians(1e-15)));
  }

  private static void checkContainsPoint(String aText, String bText) {
    S2Polygon a = makePolygon(aText);
    assertTrue(a.contains(makePoint(bText)));
  }

  // Make sure we've set things up correctly.
  public void testInit() {
    checkContains(NEAR1, NEAR0);
    checkContains(NEAR2, NEAR1);
    checkContains(NEAR3, NEAR2);
    checkContains(NEAR_HEMI, NEAR3);
    checkContains(FAR1, FAR0);
    checkContains(FAR2, FAR1);
    checkContains(FAR3, FAR2);
    checkContains(FAR_HEMI, FAR3);
    checkContains(SOUTH1, SOUTH0A);
    checkContains(SOUTH1, SOUTH0B);
    checkContains(SOUTH1, SOUTH0C);
    checkContains(SOUTH_HEMI, SOUTH2);
    checkContains(NEAR_FAR1, NEAR3);
    checkContains(NEAR_FAR1, FAR3);
    checkContains(NEAR_FAR2, NEAR3);
    checkContains(NEAR_FAR2, FAR3);

    checkContainsPoint(NEAR0, NEAR_POINT);
    checkContainsPoint(NEAR1, NEAR_POINT);
    checkContainsPoint(NEAR2, NEAR_POINT);
    checkContainsPoint(NEAR3, NEAR_POINT);
    checkContainsPoint(NEAR_HEMI, NEAR_POINT);
    checkContainsPoint(SOUTH0A, SOUTH_POINT);
    checkContainsPoint(SOUTH1, SOUTH_POINT);
    checkContainsPoint(SOUTH2, SOUTH_POINT);
    checkContainsPoint(SOUTH_HEMI, SOUTH_POINT);
  }

  public void testOriginNearPole() {
    // S2Polygon operations are more efficient if S2.origin() is near a pole.
    // (Loops that contain a pole tend to have very loose bounding boxes because
    // they span the full longitude range.  S2Polygon canonicalizes all loops so
    // that they don't contain S2.origin(), thus by placing S2.origin() near a
    // pole we minimize the number of canonical loops which contain that pole.)
    assertTrue(S2LatLng.latitude(S2.origin()).degrees() >= 80);
  }

  private static class TestCase {
    final String a;
    final String b;
    final String a_and_b;
    final String a_or_b;
    final String a_minus_b;

    TestCase(
        String a, String b, String a_and_b, String a_or_b, String a_minus_b) {
      this.a = a;
      this.b = b;
      this.a_and_b = a_and_b;
      this.a_or_b = a_or_b;
      this.a_minus_b = a_minus_b;
    }
  }

  private static TestCase[] testCases = {
    // Two triangles that share an edge.
    new TestCase(
        "4:2, 3:1, 3:3;",
        "3:1, 2:2, 3:3;",
        "",
        "4:2, 3:1, 2:2, 3:3;",
        "4:2, 3:1, 3:3;"),
    // Two vertical bars and a horizontal bar connecting them.
    new TestCase(
        "0:0, 0:2, 3:2, 3:0;   0:3, 0:5, 3:5, 3:3;",
        "1:1, 1:4, 2:4, 2:1;",
        "1:1, 1:2, 2:2, 2:1;   1:3, 1:4, 2:4, 2:3;",
        "0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0;",
        "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;   " +
            "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;"),
    // Two vertical bars and two horizontal bars centered around S2.origin().
    new TestCase(
        "1:88, 1:93, 2:93, 2:88;   -1:88, -1:93, 0:93, 0:88;",
        "-2:89, -2:90, 3:90, 3:89;   -2:91, -2:92, 3:92, 3:91;",
        "1:89, 1:90, 2:90, 2:89;   1:91, 1:92, 2:92, 2:91;   " +
            "-1:89, -1:90, 0:90, 0:89;   -1:91, -1:92, 0:92, 0:91;",
        "-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92," +
            "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91," +
            "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88;   " +
            "0:90, 0:91, 1:91, 1:90;",
        "1:88, 1:89, 2:89, 2:88;   1:90, 1:91, 2:91, 2:90;   " +
            "1:92, 1:93, 2:93, 2:92;   -1:88, -1:89, 0:89, 0:88;   " +
            "-1:90, -1:91, 0:91, 0:90;   -1:92, -1:93, 0:93, 0:92;"),
    // Two interlocking square doughnuts centered around -S2.origin().
    new TestCase(
        "-1:-93, -1:-89, 3:-89, 3:-93;   0:-92, 0:-90, 2:-90, 2:-92;",
        "-3:-91, -3:-87, 1:-87, 1:-91;   -2:-90, -2:-88, 0:-88, 0:-90;",
        "-1:-91, -1:-90, 0:-90, 0:-91;   0:-90, 0:-89, 1:-89, 1:-90;",
        "-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93;   " +
            "0:-92, 0:-91, 1:-91, 1:-90, 2:-90, 2:-92;   " +
            "-2:-90, -2:-88, 0:-88, 0:-89, -1:-89, -1:-90;",
        "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89," +
            "3:-93;   -1:-90, -1:-89, 0:-89, 0:-90;"),
    // An incredibly thin triangle intersecting a square, such that the two
    // intersection points of the triangle with the square are identical.
    // This results in a degenerate loop that needs to be handled correctly.
    new TestCase(
        "10:44, 10:46, 12:46, 12:44;",
        "11:45, 89:45.00000000000001, 90:45;",
        "",  // Empty intersection!
        // Original square with extra vertex, and triangle disappears (due to
        // default vertex_merge_radius of
        // S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE).
        "10:44, 10:46, 12:46, 12:45, 12:44;",
        "10:44, 10:46, 12:46, 12:45, 12:44;")
  };

  public void testOperations() {
    S2Polygon farSouth = new S2Polygon();
    farSouth.initToIntersection(farH, southH);
    checkEqual(farSouth, farHSouthH, 1e-15);

    for (int testNumber = 0; testNumber < testCases.length; testNumber++) {
      TestCase test = testCases[testNumber];
      logger.info("Polygon operation test case " + testNumber);
      S2Polygon a = makeVerbatimPolygon(test.a);
      S2Polygon b = makeVerbatimPolygon(test.b);
      S2Polygon expectedAAndB = makeVerbatimPolygon(test.a_and_b);
      S2Polygon expectedAOrB = makeVerbatimPolygon(test.a_or_b);
      S2Polygon expectedAMinusB = makeVerbatimPolygon(test.a_minus_b);

      S2Polygon aAndB = new S2Polygon();
      S2Polygon aOrB = new S2Polygon();
      S2Polygon aMinusB = new S2Polygon();
      aAndB.initToIntersection(a, b);
      checkEqual(aAndB, expectedAAndB, OPERATIONS_MAX_ERROR);
      aOrB.initToUnion(a, b);
      checkEqual(aOrB, expectedAOrB, OPERATIONS_MAX_ERROR);
      checkDestructiveUnion(a, b);
      aMinusB.initToDifference(a, b);
      checkEqual(aMinusB, expectedAMinusB, OPERATIONS_MAX_ERROR);
    }
  }

  private void polylineIntersectionSharedEdgeTest(
      S2Polygon p, int startVertex, int direction) {
    logger.info("Polyline intersection shared edge test start=" +
                startVertex + " direction=" + direction);
    List<S2Point> points = Lists.newArrayList();
    points.add(p.loop(0).vertex(startVertex));
    points.add(p.loop(0).vertex(startVertex + direction));
    S2Polyline polyline = new S2Polyline(points);
    List<S2Polyline> polylines;
    if (direction < 0) {
      polylines = p.intersectWithPolyline(polyline);
      assertEquals(0, polylines.size());
      polylines = p.subtractFromPolyline(polyline);
      assertEquals(1, polylines.size());
      assertEquals(2, polylines.get(0).numVertices());
      assertEquals(points.get(0), polylines.get(0).vertex(0));
      assertEquals(points.get(1), polylines.get(0).vertex(1));
    } else {
      polylines = p.intersectWithPolyline(polyline);
      assertEquals(1, polylines.size());
      assertEquals(2, polylines.get(0).numVertices());
      assertEquals(points.get(0), polylines.get(0).vertex(0));
      assertEquals(points.get(1), polylines.get(0).vertex(1));
      polylines = p.subtractFromPolyline(polyline);
      assertEquals(0, polylines.size());
    }
  }

  /**
   * This tests polyline-polyline intersections.
   * It covers the same edge cases as {@code testOperations} and also adds some
   * extra tests for shared edges.
   */
  public void testPolylineIntersection() {
    for (int v = 0; v < 3; ++v) {
      polylineIntersectionSharedEdgeTest(cross1, v, 1);
      polylineIntersectionSharedEdgeTest(cross1, v + 1, -1);
      polylineIntersectionSharedEdgeTest(cross1SideHole, v, 1);
      polylineIntersectionSharedEdgeTest(cross1SideHole, v + 1, -1);
    }

    // This duplicates some of the tests in testOperations by
    // converting the outline of polygon A to a polyline then intersecting
    // it with the polygon B. It then converts B to a polyline and intersects
    // it with A. It then feeds all of the results into a polygon builder and
    // tests that the output is equal to doing an intersection between A and B.

    for (int testNumber = 0; testNumber < testCases.length; testNumber++) {
      TestCase test = testCases[testNumber];
      logger.info("Polyline intersection test case " + testNumber);
      S2Polygon a = makeVerbatimPolygon(test.a);
      S2Polygon b = makeVerbatimPolygon(test.b);
      S2Polygon expectedAAndB = makeVerbatimPolygon(test.a_and_b);

      List<S2Point> points = Lists.newArrayList();
      List<S2Polyline> polylines = Lists.newArrayList();
      for (int ab = 0; ab < 2; ab++) {
        S2Polygon tmp = (ab == 1) ? a : b;
        S2Polygon tmp2 = (ab == 1) ? b : a;
        for (int l = 0; l < tmp.numLoops(); l++) {
          points.clear();
          if (tmp.loop(l).isHole()) {
            for (int v = tmp.loop(l).numVertices(); v >= 0; v--) {
              points.add(tmp.loop(l).vertex(v));
            }
          } else {
            for (int v = 0; v <= tmp.loop(l).numVertices(); v++) {
              points.add(tmp.loop(l).vertex(v));
            }
          }
          polylines.addAll(tmp2.intersectWithPolyline(new S2Polyline(points)));
        }
      }

      S2PolygonBuilder builder = new S2PolygonBuilder(
          S2PolygonBuilder.Options.DIRECTED_XOR);
      for (int i = 0; i < polylines.size(); i++) {
        for (int j = 0; j < polylines.get(i).numVertices() - 1; j++) {
          builder.addEdge(polylines.get(i).vertex(j),
                          polylines.get(i).vertex(j + 1));
          logger.info(" ... Adding edge: " + polylines.get(i).vertex(j) +
                      " - " + polylines.get(i).vertex(j + 1));
        }
      }

      S2Polygon a_and_b = new S2Polygon();
      assertTrue(builder.assemblePolygon(a_and_b, null));
      checkEqual(a_and_b, expectedAAndB, OPERATIONS_MAX_ERROR);
    }
  }

  private static void checkEqual(S2Polygon a, S2Polygon b) {
    checkEqual(a, b, 0);
  }

  private static void checkEqual(S2Polygon a, S2Polygon b, final double maxError) {
    if (a.boundaryApproxEquals(b, maxError)) {
      return;
    }

    S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.DIRECTED_XOR);
    builder.addPolygon(a);
    S2Polygon a2 = new S2Polygon();
    assertTrue(builder.assemblePolygon(a2, null));
    builder.addPolygon(b);
    S2Polygon b2 = new S2Polygon();
    assertTrue(builder.assemblePolygon(b2, null));
    assertTrue(a2.boundaryApproxEquals(b2, maxError));
  }

  private static void checkComplementary(S2Polygon a, S2Polygon b) {
    S2Polygon b1 = new S2Polygon();
    b1.initToComplement(b);
    checkEqual(a, b1);
  }

  public void testApproxContains() {
    // Get a random S2Cell as a polygon.
    S2CellId id = S2CellId.fromLatLng(S2LatLng.fromE6(69852241, 6751108));
    S2Cell cell = new S2Cell(id.parent(10));
    S2Polygon cellAsPolygon = new S2Polygon(cell);

    // We want to roughly bisect the polygon, so we make a rectangle that is the
    // top half of the current polygon's bounding rectangle.
    S2LatLngRect bounds = cellAsPolygon.getRectBound();
    S2LatLngRect upperHalf = new S2LatLngRect(
        new R1Interval(
            bounds.lat().getCenter(),
            bounds.lat().hi()),
        bounds.lng());

    // Turn the S2LatLngRect into an S2Polygon
    List<S2Point> points = Lists.newArrayList();
    for (int i = 0; i < 4; i++) {
      points.add(upperHalf.getVertex(i).toPoint());
    }
    List<S2Loop> loops = Lists.newArrayList();
    loops.add(new S2Loop(points));
    S2Polygon upperHalfPolygon = new S2Polygon(loops);

    // Get the intersection. There is no guarantee that the intersection will be
    // contained by A or B.
    S2Polygon intersection = new S2Polygon();
    intersection.initToIntersection(cellAsPolygon, upperHalfPolygon);
    assertFalse(cellAsPolygon.contains(intersection));
    assertTrue(cellAsPolygon.approxContains(
        intersection, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
  }

  public void tryUnion(S2Polygon a, S2Polygon b) {
    S2Polygon union = new S2Polygon();
    union.initToUnion(a, b);

    List<S2Polygon> polygons = Lists.newArrayList();
    polygons.add(new S2Polygon(a));
    polygons.add(new S2Polygon(b));
    S2Polygon destructiveUnion = S2Polygon.destructiveUnion(polygons);

    checkEqual(union, destructiveUnion);
  }

  /**
   * Given a pair of polygons where A contains B, check that various identities involving union,
   * intersection, and difference operations hold true.
   */
  private static void checkOneNestedPair(S2Polygon a, S2Polygon b) {
    assertTrue(a.contains(b));
    assertEquals(!b.isEmpty(), a.intersects(b));
    assertEquals(!b.isEmpty(), b.intersects(a));

    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    checkEqual(c, a);

    S2Polygon d = new S2Polygon();
    d.initToIntersection(a, b);
    checkEqual(d, b);

    S2Polygon e = new S2Polygon();
    e.initToDifference(b, a);
    assertTrue(e.isEmpty());
  }

  /**
   * Given a pair of disjoint polygons A and B, check that various identities involving union,
   * intersection, and difference operations hold true.
   */
  private static void checkOneDisjointPair(S2Polygon a, S2Polygon b) {
    assertFalse(a.intersects(b));
    assertFalse(b.intersects(a));
    assertEquals(b.isEmpty(), a.contains(b));
    assertEquals(a.isEmpty(), b.contains(a));

    S2Polygon ab = new S2Polygon();
    S2Polygon c = new S2Polygon();
    S2Polygon d = new S2Polygon();
    S2Polygon e = new S2Polygon();
    S2Polygon f = new S2Polygon();
    S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.DIRECTED_XOR);
    builder.addPolygon(a);
    builder.addPolygon(b);
    assertTrue(builder.assemblePolygon(ab, null));

    c.initToUnion(a, b);
    checkEqual(c, ab);

    d.initToIntersection(a, b);
    assertTrue(d.isEmpty());

    e.initToDifference(a, b);
    checkEqual(e, a);

    f.initToDifference(b, a);
    checkEqual(f, b);
  }

  /**
   * Given polygons A and B whose union covers the sphere, check that various identities involving
   * union, intersection, and difference hold true.
   */
  private static void checkOneCoveringPair(S2Polygon a, S2Polygon b) {
    assertEquals(a.isFull(), a.contains(b));
    assertEquals(b.isFull(), b.contains(a));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    assertTrue(c.isFull());
  }

  /**
   * Given polygons A and B such that both A and its complement intersect both B and its complement,
   * check that various identities involving union, intersection, and difference hold true.
   */
  private static void checkOneOverlappingPair(S2Polygon a, S2Polygon b) {
    assertFalse(a.contains(b));
    assertFalse(b.contains(a));
    assertTrue(a.intersects(b));

    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    assertFalse(c.isFull());

    S2Polygon d = new S2Polygon();
    d.initToIntersection(a, b);
    assertFalse(d.isEmpty());

    S2Polygon e = new S2Polygon();
    e.initToDifference(b, a);
    assertFalse(e.isEmpty());
  }

  /**
   * Given a pair of polygons where A contains B, test various identities involving A, B, and their
   * complements.
   */
  private static void checkNestedPair(S2Polygon a, S2Polygon b) {
    S2Polygon a1 = new S2Polygon();
    S2Polygon b1 = new S2Polygon();
    a1.initToComplement(a);
    b1.initToComplement(b);
    checkOneNestedPair(a, b);
    checkOneNestedPair(b1, a1);
    checkOneDisjointPair(a1, b);
    checkOneCoveringPair(a, b1);
  }

  /**
   * Given a pair of disjoint polygons A and B, test various identities involving A, B, and their
   * complements.
   */
  private static void checkDisjointPair(S2Polygon a, S2Polygon b) {
    S2Polygon a1 = new S2Polygon();
    a1.initToComplement(a);
    S2Polygon b1 = new S2Polygon();
    b1.initToComplement(b);
    checkOneDisjointPair(a, b);
    checkOneCoveringPair(a1, b1);
    checkOneNestedPair(a1, b);
    checkOneNestedPair(b1, a);
  }

  /**
   * Given polygons A and B such that both A and its complement intersect both B and its complement,
   * test various identities involving these four polygons.
   */
  private static void checkOverlappingPair(S2Polygon a, S2Polygon b) {
    S2Polygon a1 = new S2Polygon();
    a1.initToComplement(a);
    S2Polygon b1 = new S2Polygon();
    b1.initToComplement(b);
    checkOneOverlappingPair(a, b);
    checkOneOverlappingPair(a1, b1);
    checkOneOverlappingPair(a1, b);
    checkOneOverlappingPair(a, b1);
  }

  /** "a1" is the complement of "a", and "b1" is the complement of "b". */
  private static void checkOneComplementPair(S2Polygon a, S2Polygon a1, S2Polygon b, S2Polygon b1) {
    // Check DeMorgan's Law and that subtraction is the same as intersection
    // with the complement.  This function is called multiple times in order to
    // test the various combinations of complements.
    S2Polygon aOrB = new S2Polygon();
    S2Polygon aAndB1 = new S2Polygon();
    S2Polygon aMinusB = new S2Polygon();
    aAndB1.initToIntersection(a, b1);
    aOrB.initToUnion(a1, b);
    aMinusB.initToDifference(a, b);
    checkComplementary(aOrB, aAndB1);
    checkEqual(aMinusB, aAndB1);
  }

  /**
   * Test identities that should hold for any pair of polygons A, B and their complements.
   */
  private static void checkComplements(S2Polygon a, S2Polygon b) {
    S2Polygon a1 = new S2Polygon();
    a1.initToComplement(a);
    S2Polygon b1 = new S2Polygon();
    b1.initToComplement(b);
    checkOneComplementPair(a, a1, b, b1);
    checkOneComplementPair(a1, a, b, b1);
    checkOneComplementPair(a, a1, b1, b);
    checkOneComplementPair(a1, a, b1, b);
  }

  private static void checkDestructiveUnion(S2Polygon a, S2Polygon b) {
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    List<S2Polygon> polygons = Lists.newArrayList();
    polygons.add(new S2Polygon(a));
    polygons.add(new S2Polygon(b));
    S2Polygon cDestructive = S2Polygon.destructiveUnion(polygons);
    checkEqual(c, cDestructive);
  }

  private static void checkRelationImpl(S2Polygon a, S2Polygon b,
      boolean contains, boolean contained, boolean intersects) {
    assertEquals(contains, a.contains(b));
    assertEquals(contained, b.contains(a));
    assertEquals(intersects, a.intersects(b));
    if (contains) {
      checkNestedPair(a, b);
    }
    if (contained) {
      checkNestedPair(b, a);
    }
    if (!intersects) {
      checkDisjointPair(a, b);
    }
    if (intersects && !(contains | contained)) {
      checkOverlappingPair(a, b);
    }
    checkDestructiveUnion(a, b);
    checkComplements(a, b);
  }

  private static void checkRelation(S2Polygon a, S2Polygon b,
      boolean contains, boolean contained, boolean intersects) {
    try {
      checkRelationImpl(a, b, contains, contained, intersects);
    } catch (AssertionError e) {
      System.err.println("args " + a  + ", " + b);
      throw e;
    }
  }

  public void testRelations() {
    checkRelation(near10, empty, true, false, false);
    checkRelation(near10, near10, true, true, true);
    checkRelation(full, near10, true, false, true);
    checkRelation(near10, near30, false, true, true);
    checkRelation(near10, near32, false, false, false);
    checkRelation(near10, near3210, false, true, true);
    checkRelation(near10, nearH3210, false, false, false);
    checkRelation(near30, near32, true, false, true);
    checkRelation(near30, near3210, true, false, true);
    checkRelation(near30, nearH3210, false, false, true);
    checkRelation(near32, near3210, false, true, true);
    checkRelation(near32, nearH3210, false, false, false);
    checkRelation(near3210, nearH3210, false, false, false);

    checkRelation(far10, far21, false, false, false);
    checkRelation(far10, far321, false, true, true);
    checkRelation(far10, farH20, false, false, false);
    checkRelation(far10, farH3210, false, false, false);
    checkRelation(far21, far321, false, false, false);
    checkRelation(far21, farH20, false, false, false);
    checkRelation(far21, farH3210, false, true, true);
    checkRelation(far321, farH20, false, false, true);
    checkRelation(far321, farH3210, false, false, true);
    checkRelation(farH20, farH3210, false, false, true);

    checkRelation(south0ab, south2, false, true, true);
    checkRelation(south0ab, south210b, false, false, true);
    checkRelation(south0ab, southH21, false, true, true);
    checkRelation(south0ab, southH20abc, false, true, true);
    checkRelation(south2, south210b, true, false, true);
    checkRelation(south2, southH21, false, false, true);
    checkRelation(south2, southH20abc, false, false, true);
    checkRelation(south210b, southH21, false, false, true);
    checkRelation(south210b, southH20abc, false, false, true);
    checkRelation(southH21, southH20abc, true, false, true);

    checkRelation(nf1n10f2s10abc, nf2n2f210s210ab, false, false, true);
    checkRelation(nf1n10f2s10abc, near32, true, false, true);
    checkRelation(nf1n10f2s10abc, far21, false, false, false);
    checkRelation(nf1n10f2s10abc, south0ab, false, false, false);
    checkRelation(nf1n10f2s10abc, f32n0, true, false, true);

    checkRelation(nf2n2f210s210ab, near10, false, false, false);
    checkRelation(nf2n2f210s210ab, far10, true, false, true);
    checkRelation(nf2n2f210s210ab, south210b, true, false, true);
    checkRelation(nf2n2f210s210ab, south0ab, true, false, true);
    checkRelation(nf2n2f210s210ab, n32s0b, true, false, true);

    checkRelation(cross1, cross2, false, false, true);
    checkRelation(cross1SideHole, cross2, false, false, true);
    checkRelation(cross1CenterHole, cross2, false, false, true);
    checkRelation(cross1, cross2SideHole, false, false, true);
    checkRelation(cross1, cross2CenterHole, false, false, true);
    checkRelation(cross1SideHole, cross2SideHole, false, false, true);
    checkRelation(cross1CenterHole, cross2SideHole, false, false, true);
    checkRelation(cross1SideHole, cross2CenterHole, false, false, true);
    checkRelation(cross1CenterHole, cross2CenterHole, false, false, true);

    // These cases, when either polygon has a hole, test a different code path
    // from the other cases.
    checkRelation(overlap1, overlap2, false, false, true);
    checkRelation(overlap1SideHole, overlap2, false, false, true);
    checkRelation(overlap1CenterHole, overlap2, false, false, true);
    checkRelation(overlap1, overlap2SideHole, false, false, true);
    checkRelation(overlap1, overlap2CenterHole, false, false, true);
    checkRelation(overlap1SideHole, overlap2SideHole, false, false, true);
    checkRelation(overlap1CenterHole, overlap2SideHole, false, false, true);
    checkRelation(overlap1SideHole, overlap2CenterHole, false, false, true);
    checkRelation(overlap1CenterHole, overlap2CenterHole, false, false, true);
  }

  public void testEmptyAndFull() {
    assertTrue(empty.isEmpty());
    assertFalse(full.isEmpty());
    assertFalse(empty.isFull());
    assertTrue(full.isFull());
    checkNestedPair(empty, empty);
    checkNestedPair(full, empty);
    checkNestedPair(full, full);
  }

  public void testUnionSloppySuccess() {
    List<S2Polygon> polygons = Lists.newArrayList();
    polygons.add(adj0);
    polygons.add(adj1);
    S2Polygon union = S2Polygon.destructiveUnionSloppy(polygons, S1Angle.degrees(0.1));

    assertEquals(1, union.numLoops());
    if (union.numLoops() != 1) {
      return;
    }
    S2Loop s2Loop = union.loop(0);
    assertEquals(8, s2Loop.numVertices());
    if (s2Loop.numVertices() != 8) {
      return;
    }
    S2Loop expected = makeLoop("2:0, 1:0, 0:0, 0:1, 0:2, 1:2, 2:2, 2:1");
    double maxError = S1Angle.degrees(0.01000001).radians();
    assertTrue(expected.boundaryApproxEquals(s2Loop, maxError));
  }

  public void testUnionSloppyFailure() {
    List<S2Polygon> polygons = Lists.newArrayList();
    polygons.add(adj0);
    polygons.add(unAdj);
    // The polygons are sufficiently far apart that this angle will not
    // bring them together:
    S2Polygon union = S2Polygon.destructiveUnionSloppy(polygons, S1Angle.degrees(0.1));

    assertEquals(2, union.numLoops());
  }

  public void testCompareTo() {
    // Polygons with same loops, but in different order:
    S2Polygon p1 = makePolygon(RECTANGLE1 + RECTANGLE2);
    S2Polygon p2 = makePolygon(RECTANGLE2 + RECTANGLE1);
    assertEquals(0, p1.compareTo(p2));

    // Polygons with same loops, but in different order and containins a
    // different number of points.
    S2Polygon p3 = makePolygon(RECTANGLE1 + TRIANGLE);
    S2Polygon p4 = makePolygon(TRIANGLE + RECTANGLE1);
    assertEquals(0, p3.compareTo(p4));

    // Polygons with same logical loop (but loop is reordered).
    S2Polygon p5 = makePolygon(TRIANGLE);
    S2Polygon p6 = makePolygon(TRIANGLE_ROT);
    assertEquals(0, p5.compareTo(p6));

    // Polygons with a differing number of loops
    S2Polygon p7 = makePolygon(RECTANGLE1 + RECTANGLE2);
    S2Polygon p8 = makePolygon(TRIANGLE);
    assertTrue(0 > p8.compareTo(p7));
    assertTrue(0 < p7.compareTo(p8));

    // Polygons with a differing number of loops (one a subset of the other)
    S2Polygon p9 = makePolygon(RECTANGLE1 + RECTANGLE2 + TRIANGLE);
    S2Polygon p10 = makePolygon(RECTANGLE1 + RECTANGLE2);
    assertTrue(0 < p9.compareTo(p10));
    assertTrue(0 > p10.compareTo(p9));
  }

  public void testGetDistance() {
    // Error margin since we're doing numerical computations
    double epsilon = 1e-15;

    // A rectangle with (lat,lng) vertices (3,1), (3,-1), (-3,-1) and (-3,1)
    String inner = "3:1, 3:-1, -3:-1, -3:1;";
    // A larger rectangle with (lat,lng) vertices (4,2), (4,-2), (-4,-2) and
    // (-4,s)
    String outer = "4:2, 4:-2, -4:-2, -4:2;";


    S2Polygon rect = makePolygon(inner);
    S2Polygon shell = makePolygon(inner + outer);

    // All of the vertices of a polygon should be distance 0
    for (int i = 0; i < shell.numLoops(); i++) {
      for (int j = 0; j < shell.loop(i).numVertices(); j++) {
        assertEquals(0d, shell.getDistance(shell.loop(i).vertex(j)).radians(), epsilon);
      }
    }

    // A non-vertex point on an edge should be distance 0
    assertEquals(0d, rect.getDistance(
        S2Point.normalize(S2Point.add(rect.loop(0).vertex(0), rect.loop(0).vertex(1)))).radians(),
        epsilon);

    S2Point origin = S2LatLng.fromDegrees(0, 0).toPoint();
    // rect contains the origin
    assertEquals(0d, rect.getDistance(origin).radians(), epsilon);

    // shell does NOT contain the origin, since it has a hole. The shortest
    // distance is to (1,0) or (-1,0), and should be 1 degree
    assertEquals(1d, shell.getDistance(origin).degrees(), epsilon);
  }

  public void testMultipleInit() {
    S2Polygon polygon = makePolygon("0:0, 0:2, 2:0");
    assertEquals(1, polygon.numLoops());
    assertEquals(3, polygon.getNumVertices());
    S2LatLngRect bound1 = polygon.getRectBound();

    List<S2Loop> loops = Lists.newArrayList();
    loops.add(makeLoop("10:0, -10:-20, -10:20"));
    loops.add(makeLoop("40:30, 20:10, 20:50"));
    polygon.initNested(loops);
    assertTrue(polygon.isValid());
    assertTrue(loops.isEmpty());
    assertEquals(2, polygon.numLoops());
    assertEquals(6, polygon.getNumVertices());
    assertTrue(!bound1.equals(polygon.getRectBound()));
  }

  public void testProject() {
    S2Polygon polygon = makeVerbatimPolygon(NEAR0 + NEAR2);
    double epsilon = 1e-15;
    S2Point point;
    S2Point projected;

    // The point inside the polygon should be projected into itself.
    point = makePoint("1.1:0");
    projected = polygon.project(point);
    assertTrue(point.aequal(projected, epsilon));

    // The point is on the outside of the polygon.
    point = makePoint("5.1:-2");
    projected  = polygon.project(point);
    assertTrue(makePoint("5:-2").aequal(projected, epsilon));

    // The point is inside the hole in the polygon. Note the expected value
    // is based on a plane, so it's not that accurate; thus, tolerance is
    // reduced to 1e-6.
    point = makePoint("-0.49:-0.49");
    projected = polygon.project(point);
    assertTrue(makePoint("-0.5:-0.5").aequal(projected, 1e-6));

    point = makePoint("0:-3");
    projected = polygon.project(point);
    assertTrue(makePoint("0:-2").aequal(projected, epsilon));
  }

  public void testProjectMatchesDistance() {
    S2Polygon polygon = makePolygon(NEAR0 + NEAR2);
    double epsilon = 1e-15;
    S2Point point;
    S2Point projected;

    // In the hole
    point = makePoint("-0.49:-0.49");
    projected = polygon.project(point);
    assertEquals(polygon.getDistance(point).radians(), point.angle(projected), epsilon);

    // Outside the polygon
    point = makePoint("10:15");
    projected = polygon.project(point);
    assertEquals(polygon.getDistance(point).radians(), point.angle(projected), epsilon);
  }

  public void testFastInit() {
    S2LatLngRect bound = null;
    Map<S2Loop, List<S2Loop>> nestedLoops = Maps.newHashMap();
    nestedLoops.put(null, Lists.<S2Loop>newArrayList());

    List<S2Point> vertices = Lists.newArrayList();

    bound = parseVertices("-2:-2, -3:6, 6:-3", vertices);
    S2Loop loop1 = S2Loop.newLoopWithTrustedDetails(vertices, false, bound);
    nestedLoops.get(null).add(loop1);
    nestedLoops.put(loop1, Lists.<S2Loop>newArrayList());

    vertices = Lists.newArrayList();
    bound = parseVertices("-1:-2, -2:5, 5:-2", vertices);
    S2Loop loop2 = S2Loop.newLoopWithTrustedDetails(vertices, false, bound);
    nestedLoops.get(loop1).add(loop2);
    nestedLoops.put(loop2, Lists.<S2Loop>newArrayList());

    vertices = Lists.newArrayList();
    bound = parseVertices("-1:0, 0:1, 1:0, 0:-1", vertices);
    S2Loop loop3 = S2Loop.newLoopWithTrustedDetails(vertices, false, bound);
    nestedLoops.get(loop2).add(loop3);
    nestedLoops.put(loop3, Lists.<S2Loop>newArrayList());

    S2Polygon polygon = new S2Polygon();
    polygon.initWithNestedLoops(nestedLoops);

    assertEquals(0, polygon.compareTo(makePolygon(NEAR0 + NEAR2 + NEAR3)));

    assertTrue(polygon.isValid());
    assertEquals(0, polygon.loop(0).depth());
    assertEquals(1, polygon.loop(1).depth());
    assertEquals(2, polygon.loop(2).depth());
    assertDoubleNear(0.003821967440517272, polygon.getArea());
  }

  public void testInitToSnappedWithSnapLevel() {
    S2Polygon polygon = makePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");

    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      S2Polygon snappedPolygon = new S2Polygon();
      snappedPolygon.initToSnapped(polygon, level);
      assertTrue(snappedPolygon.isValid());
      double cellAngle = PROJ.maxDiag.getValue(level);
      S1Angle mergeRadius = S1Angle.radians(cellAngle);
      assertTrue("snapped polygon should approx contain original polygon for"
          + "\nsnap level = " + level + ", mergeRadius = " + mergeRadius
          + "\noriginal polygon: " + polygon
          + "\nsnapped polygon: " + snappedPolygon,
          snappedPolygon.approxContains(polygon, mergeRadius));
    }
  }

  /**
   * Verifies that clipBoundary can succeed with duplicate adjacent vertices. Although such a case
   * means the polygon is invalid, it is common to fix invalidity issues by doign a self-
   * intersection to node crossings and drop duplicates.
   */
  public void testDuplicatePointClipping() {
    S2Polygon p = makePolygon("0:0, 0:0, 0:4, 4:4, 4:0");
    S2Polygon fixed = new S2Polygon();
    fixed.initToIntersection(p, p);
    assertEquals(makePolygon("0:0, 0:4, 4:4, 4:0"), fixed);
  }
}
