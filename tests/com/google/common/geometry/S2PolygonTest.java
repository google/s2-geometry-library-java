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

import com.google.common.collect.Lists;

import java.util.List;
import java.util.logging.Logger;

/**
 * Tests for {@link S2Polygon}.
 *
 */
public strictfp class S2PolygonTest extends GeometryTestCase {
  private static Logger logger = Logger.getLogger(S2PolygonTest.class.getName());

  // A set of nested loops around the point 0:0 (lat:lng).
  // Every vertex of NEAR0 is a vertex of NEAR1.
  private static final String NEAR0 = "-1:0, 0:1, 1:0, 0:-1;";
  private static final String NEAR1 = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;";
  private static final String NEAR2 = "5:-2, -2:5, -1:-2;";
  private static final String NEAR3 = "6:-3, -3:6, -2:-2;";
  private static final String NEAR_HEMI = "0:-90, -90:0, 0:90, 90:0;";

  // A set of nested loops around the point 0:180 (lat:lng).
  // Every vertex of FAR0 and FAR2 belongs to FAR1, and all
  // the loops except FAR2 are non-convex.
  private static final String FAR0 = "0:179, 1:180, 0:-179, 2:-180;";
  private static final String FAR1 =
      "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;";
  private static final String FAR2 = "-1:-179, -1:179, 3:178, 3:-178;"; // opposite
                                                                        // direction
  private static final String FAR3 = "-3:-178, -2:179, -3:178, 4:177, 4:-177;";
  private static final String FAR_HEMI = "0:-90, 60:90, -60:90;";

  // A set of nested loops around the point -90:0 (lat:lng).
  private static final String SOUTH0a = "-90:0, -89.99:0, -89.99:0.01;";
  private static final String SOUTH0b = "-90:0, -89.99:0.02, -89.99:0.03;";
  private static final String SOUTH0c = "-90:0, -89.99:0.04, -89.99:0.05;";
  private static final String SOUTH1 = "-90:0, -89.9:-0.1, -89.9:0.1;";
  private static final String SOUTH2 = "-90:0, -89.8:-0.2, -89.8:0.2;";
  private static final String SOUTH_HEMI = "0:-180, 0:60, 0:-60;";

  // Two different loops that surround all the Near and Far loops except
  // for the hemispheres.
  private static final String NEAR_FAR1 =
      "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, " + "1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;";
  private static final String NEAR_FAR2 =
      "-8:-4, 8:-4, 2:15, 2:170, 8:-175, -8:-175, -2:170, -2:15;";

  // Loops that result from intersection of other loops.
  private static final String FAR_SOUTH_H = "0:-180, 0:90, -60:90, 0:-90;";

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

  private void assertContains(String aStr, String bStr) {
    S2Polygon a = makePolygon(aStr);
    S2Polygon b = makePolygon(bStr);
    assertTrue(a.contains(b));
  }

  // Make sure we've set things up correctly.
  public void testInit() {
    assertContains(NEAR1, NEAR0);
    assertContains(NEAR2, NEAR1);
    assertContains(NEAR3, NEAR2);
    assertContains(NEAR_HEMI, NEAR3);
    assertContains(FAR1, FAR0);
    assertContains(FAR2, FAR1);
    assertContains(FAR3, FAR2);
    assertContains(FAR_HEMI, FAR3);
    assertContains(SOUTH1, SOUTH0a);
    assertContains(SOUTH1, SOUTH0b);
    assertContains(SOUTH1, SOUTH0c);
    assertContains(SOUTH_HEMI, SOUTH2);
    assertContains(NEAR_FAR1, NEAR3);
    assertContains(NEAR_FAR1, FAR3);
    assertContains(NEAR_FAR2, NEAR3);
    assertContains(NEAR_FAR2, FAR3);
  }

  S2Polygon near10 = makePolygon(NEAR0 + NEAR1);
  S2Polygon near30 = makePolygon(NEAR3 + NEAR0);
  S2Polygon near32 = makePolygon(NEAR2 + NEAR3);
  S2Polygon near3210 = makePolygon(NEAR0 + NEAR2 + NEAR3 + NEAR1);
  S2Polygon nearH3210 = makePolygon(NEAR0 + NEAR2 + NEAR3 + NEAR_HEMI + NEAR1);

  S2Polygon far10 = makePolygon(FAR0 + FAR1);
  S2Polygon far21 = makePolygon(FAR2 + FAR1);
  S2Polygon far321 = makePolygon(FAR2 + FAR3 + FAR1);
  S2Polygon farH20 = makePolygon(FAR2 + FAR_HEMI + FAR0);
  S2Polygon farH3210 = makePolygon(FAR2 + FAR_HEMI + FAR0 + FAR1 + FAR3);

  S2Polygon south0ab = makePolygon(SOUTH0a + SOUTH0b);
  S2Polygon south2 = makePolygon(SOUTH2);
  S2Polygon south210b = makePolygon(SOUTH2 + SOUTH0b + SOUTH1);
  S2Polygon southH21 = makePolygon(SOUTH2 + SOUTH_HEMI + SOUTH1);
  S2Polygon southH20abc = makePolygon(SOUTH2 + SOUTH0b + SOUTH_HEMI + SOUTH0a + SOUTH0c);

  S2Polygon nf1n10f2s10abc =
      makePolygon(SOUTH0c + FAR2 + NEAR1 + NEAR_FAR1 + NEAR0 + SOUTH1 + SOUTH0b + SOUTH0a);

  S2Polygon nf2n2f210s210ab =
      makePolygon(FAR2 + SOUTH0a + FAR1 + SOUTH1 + FAR0 + SOUTH0b + NEAR_FAR2 + SOUTH2 + NEAR2);

  S2Polygon f32n0 = makePolygon(FAR2 + NEAR0 + FAR3);
  S2Polygon n32s0b = makePolygon(NEAR3 + SOUTH0b + NEAR2);

  S2Polygon adj0 = makePolygon(ADJACENT0);
  S2Polygon adj1 = makePolygon(ADJACENT1);
  S2Polygon unAdj = makePolygon(UN_ADJACENT);

  S2Polygon farH = makePolygon(FAR_HEMI);
  S2Polygon southH = makePolygon(SOUTH_HEMI);
  S2Polygon farHSouthH = makePolygon(FAR_SOUTH_H);

  private void assertRelation(S2Polygon a, S2Polygon b, int contains, boolean intersects) {
    assertEquals(a.contains(b), contains > 0);
    assertEquals(b.contains(a), contains < 0);
    assertEquals(a.intersects(b), intersects);
  }

  public void testRelations() {
    assertRelation(near10, near30, -1, true);
    assertRelation(near10, near32, 0, false);
    assertRelation(near10, near3210, -1, true);
    assertRelation(near10, nearH3210, 0, false);
    assertRelation(near30, near32, 1, true);
    assertRelation(near30, near3210, 1, true);
    assertRelation(near30, nearH3210, 0, true);
    assertRelation(near32, near3210, -1, true);
    assertRelation(near32, nearH3210, 0, false);
    assertRelation(near3210, nearH3210, 0, false);

    assertRelation(far10, far21, 0, false);
    assertRelation(far10, far321, -1, true);
    assertRelation(far10, farH20, 0, false);
    assertRelation(far10, farH3210, 0, false);
    assertRelation(far21, far321, 0, false);
    assertRelation(far21, farH20, 0, false);
    assertRelation(far21, farH3210, -1, true);
    assertRelation(far321, farH20, 0, true);
    assertRelation(far321, farH3210, 0, true);
    assertRelation(farH20, farH3210, 0, true);

    assertRelation(south0ab, south2, -1, true);
    assertRelation(south0ab, south210b, 0, true);
    assertRelation(south0ab, southH21, -1, true);
    assertRelation(south0ab, southH20abc, -1, true);
    assertRelation(south2, south210b, 1, true);
    assertRelation(south2, southH21, 0, true);
    assertRelation(south2, southH20abc, 0, true);
    assertRelation(south210b, southH21, 0, true);
    assertRelation(south210b, southH20abc, 0, true);
    assertRelation(southH21, southH20abc, 1, true);

    assertRelation(nf1n10f2s10abc, nf2n2f210s210ab, 0, true);
    assertRelation(nf1n10f2s10abc, near32, 1, true);
    assertRelation(nf1n10f2s10abc, far21, 0, false);
    assertRelation(nf1n10f2s10abc, south0ab, 0, false);
    assertRelation(nf1n10f2s10abc, f32n0, 1, true);

    assertRelation(nf2n2f210s210ab, near10, 0, false);
    assertRelation(nf2n2f210s210ab, far10, 1, true);
    assertRelation(nf2n2f210s210ab, south210b, 1, true);
    assertRelation(nf2n2f210s210ab, south0ab, 1, true);
    assertRelation(nf2n2f210s210ab, n32s0b, 1, true);
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
  };

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
    checkEqual(farSouth, farHSouthH);

    int testNumber = 0;
    for (TestCase test : testCases) {
      logger.info("Polygon operation test case " + testNumber++);
      S2Polygon a = makePolygon(test.a);
      S2Polygon b = makePolygon(test.b);
      S2Polygon expected_a_and_b = makePolygon(test.a_and_b);
      S2Polygon expected_a_or_b = makePolygon(test.a_or_b);
      S2Polygon expected_a_minus_b = makePolygon(test.a_minus_b);

      // The intersections in the "expected" data were computed in lat-lng
      // space, while the actual intersections are computed using geodesics.
      // The error due to this depends on the length and direction of the line
      // segment being intersected, and how close the intersection is to the
      // endpoints of the segment.  The worst case is for a line segment between
      // two points at the same latitude, where the intersection point is in the
      // middle of the segment.  In this case the error is approximately
      // (p * t^2) / 8, where "p" is the absolute latitude in radians, "t" is
      // the longitude difference in radians, and both "p" and "t" are small.
      // The test cases all have small latitude and longitude differences.
      // If "p" and "t" are converted to degrees, the following error bound is
      // valid as long as (p * t^2 < 150).

      final double maxError = 1e-4;

      S2Polygon a_and_b = new S2Polygon();
      S2Polygon a_or_b = new S2Polygon();
      S2Polygon a_minus_b = new S2Polygon();
      a_and_b.initToIntersection(a, b);
      checkEqual(a_and_b, expected_a_and_b, maxError);
      a_or_b.initToUnion(a, b);
      tryUnion(a, b);
      checkEqual(a_or_b, expected_a_or_b, maxError);
      a_minus_b.initToDifference(a, b);
      checkEqual(a_minus_b, expected_a_minus_b, maxError);
    }
  }

  private void checkEqual(S2Polygon a, S2Polygon b) {
    checkEqual(a, b, 1e-31);
  }

  private void checkEqual(S2Polygon a, S2Polygon b, final double maxError) {
    if (a.isNormalized() && b.isNormalized()) {
      boolean r = a.boundaryApproxEquals(b, maxError);
      assertTrue(r);
    } else {
      S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.UNDIRECTED_XOR);
      S2Polygon a2 = new S2Polygon();
      S2Polygon b2 = new S2Polygon();
      builder.addPolygon(a);
      assertTrue(builder.assemblePolygon(a2, null));
      builder.addPolygon(b);
      assertTrue(builder.assemblePolygon(b2, null));
      assertTrue(a2.boundaryApproxEquals(b2, maxError));
    }
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

  public void testDisjoint() {
    S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.UNDIRECTED_XOR);
    builder.addPolygon(adj0);
    builder.addPolygon(unAdj);
    S2Polygon ab = new S2Polygon();
    assertTrue(builder.assemblePolygon(ab, null));

    S2Polygon union = new S2Polygon();
    union.initToUnion(adj0, unAdj);
    assertEquals(2, union.numLoops());

    checkEqual(ab, union);
    tryUnion(adj0, unAdj);
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
}
