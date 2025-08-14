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

import static com.google.common.collect.ImmutableSet.toImmutableSet;
import static com.google.common.geometry.S2CellId.FACE_CELLS;
import static com.google.common.geometry.S2Projections.MAX_DIAG;
import static com.google.common.geometry.S2Projections.MAX_WIDTH;
import static com.google.common.geometry.S2TextFormat.makeLoop;
import static com.google.common.geometry.S2TextFormat.makePoint;
import static com.google.common.geometry.S2TextFormat.makePolygon;
import static com.google.common.geometry.S2TextFormat.makePolygonOrDie;
import static com.google.common.geometry.S2TextFormat.makePolygonVerbatimOrDie;
import static com.google.common.geometry.S2TextFormat.parseVertices;
import static com.google.common.geometry.TestDataGenerator.concentricLoopsPolygon;
import static java.lang.Math.PI;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.S2CellIdSnapFunction;
import com.google.common.geometry.S2Shape.ChainPosition;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.io.BaseEncoding;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Tests for {@link S2Polygon}.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@RunWith(JUnit4.class)
@SuppressWarnings("IdentifierName")
public class S2PolygonTest extends GeometryTestCase {
  private static final Logger logger = Logger.getLogger(S2PolygonTest.class.getName());

  /**
   * The error margin used when comparing the actual and expected results of constructive geometry
   * operations. The intersections in the expected data were computed in lat-lng space, while the
   * actual intersections are computed using geodesics. The error due to this depends on the length
   * and direction of the line segment being intersected, and how close the intersection is to the
   * endpoints of the segment. The worst case is for a line segment between two points at the same
   * latitude, where the intersection point is in the middle of the segment. In this case the error
   * is approximately {@code (p * t^2) / 8}, where {@code p} is the absolute latitude in radians,
   * {@code t} is the longitude difference in radians, and both {@code p} and {@code t} are small.
   * The test cases all have small latitude and longitude differences. If {@code p} and {@code t}
   * are converted to degrees, the following error bound is valid as long as {@code (p * t^2 <
   * 150)}.
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
  // Every vertex of FAR0 and FAR2 belongs to FAR1, and all the loops except FAR2 are non-convex.
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

  // Two different loops that surround all the NEAR and FAR loops except for the hemispheres.
  private static final String NEAR_FAR1 =
      "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, " + "1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;";
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

  // Two rectangles that intersect, but no edges cross and there's always local containment (rather
  // than crossing) at each shared vertex. In this ugly ASCII art, 1 is A+B, 2 is B+C:
  //     +---+---+---+
  //     | A | B | C |
  //     +---+---+---+
  private static final String OVERLAP1 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
  private static final String OVERLAP1_SIDE_HOLE = "0.2:0.8, 0.8:0.8, 0.8:0.2, 0.2:0.2;";
  private static final String OVERLAP2 = "1:1, 2:1, 3:1, 3:0, 2:0, 1:0;";
  private static final String OVERLAP2_SIDE_HOLE = "2.2:0.8, 2.8:0.8, 2.8:0.2, 2.2:0.2;";
  private static final String OVERLAP_CENTER_HOLE = "1.2:0.8, 1.8:0.8, 1.8:0.2, 1.2:0.2;";

  // By symmetry, the intersection of the two polygons has almost half the area of either polygon.
  private static final String OVERLAP_3 = "-10:10, 0:10, 0:-10, -10:-10, -10:0";
  private static final String OVERLAP_4 = "-10:0, 10:0, 10:-10, -10:-10";

  // The intersection of polygon B on A should be ~.25 while the intersection of A on B should
  // be ~.5.
  private static final String OVERLAP_A = "-10:10, 0:10, 0:-10, -10:-10, -10:0";
  private static final String OVERLAP_B = "-10:0, 10:0, 10:-5, -10:-5";

  // Two rectangles that are "adjacent", but rather than having common edges, those edges are
  // slightly off. A third rectangle that is not adjacent to either of the first two.
  private static final String ADJACENT0 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
  private static final String ADJACENT1 = "0:2, 1:2, 2:2, 2:1.01, 1:0.99, 0:1.01;";
  private static final String UN_ADJACENT = "10:10, 11:10, 12:10, 12:9, 11:9, 10:9;";

  // Shapes used to test comparison functions for polygons.
  private static final String RECTANGLE1 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
  private static final String RECTANGLE2 = "5:1, 6:1, 7:1, 7:0, 6:0, 5:0;";
  private static final String TRIANGLE = "15:0, 17:0, 16:2;";
  private static final String TRIANGLE_ROT = "17:0, 16:2, 15:0;";

  private final S2Polygon empty = new S2Polygon();
  private final S2Polygon full = new S2Polygon(S2Loop.full());
  private final S2Polygon near0 = makeVerbatimPolygon(NEAR0);
  private final S2Polygon near10 = makeVerbatimPolygon(NEAR0 + NEAR1);
  private final S2Polygon near30 = makeVerbatimPolygon(NEAR3 + NEAR0);
  private final S2Polygon near32 = makeVerbatimPolygon(NEAR2 + NEAR3);
  private final S2Polygon near3210 = makeVerbatimPolygon(NEAR0 + NEAR2 + NEAR3 + NEAR1);
  private final S2Polygon nearH3210 =
      makeVerbatimPolygon(NEAR0 + NEAR2 + NEAR3 + NEAR_HEMI + NEAR1);

  private final S2Polygon far10 = makeVerbatimPolygon(FAR0 + FAR1);
  private final S2Polygon far21 = makeVerbatimPolygon(FAR2 + FAR1);
  private final S2Polygon far321 = makeVerbatimPolygon(FAR2 + FAR3 + FAR1);
  private final S2Polygon farH20 = makeVerbatimPolygon(FAR2 + FAR_HEMI + FAR0);
  private final S2Polygon farH3210 = makeVerbatimPolygon(FAR2 + FAR_HEMI + FAR0 + FAR1 + FAR3);

  private final S2Polygon south0ab = makeVerbatimPolygon(SOUTH0A + SOUTH0B);
  private final S2Polygon south2 = makeVerbatimPolygon(SOUTH2);
  private final S2Polygon south210b = makeVerbatimPolygon(SOUTH2 + SOUTH0B + SOUTH1);
  private final S2Polygon southH21 = makeVerbatimPolygon(SOUTH2 + SOUTH_HEMI + SOUTH1);
  private final S2Polygon southH20abc =
      makeVerbatimPolygon(SOUTH2 + SOUTH0B + SOUTH_HEMI + SOUTH0A + SOUTH0C);

  private final S2Polygon nf1n10f2s10abc =
      makeVerbatimPolygon(SOUTH0C + FAR2 + NEAR1 + NEAR_FAR1 + NEAR0 + SOUTH1 + SOUTH0B + SOUTH0A);

  private final S2Polygon nf2n2f210s210ab =
      makeVerbatimPolygon(
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

  private final S2Polygon overlapA = makeVerbatimPolygon(OVERLAP_A);
  private final S2Polygon overlapB = makeVerbatimPolygon(OVERLAP_B);

  static final S2Polygon POLYGON = makePolygon("0:0, 0:1, 1:0");
  static final S2Polygon SNAPPED_POLYGON =
      new S2Polygon(
          new S2Loop(
              ImmutableList.of(
                  FACE_CELLS[0].toPoint(), FACE_CELLS[1].toPoint(), FACE_CELLS[2].toPoint())));

  private static final int VALIDITY_ITERS = 100;

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

  private static byte[] encode(S2Polygon polygon) throws IOException {
    ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
    polygon.encode(outputStream);
    return outputStream.toByteArray();
  }

  private static S2Polygon decode(byte[] encoded) throws IOException {
    return S2Polygon.decode(new ByteArrayInputStream(encoded));
  }

  private static void encodeDecode(S2Polygon polygon) throws IOException {
    byte[] encoded = encode(polygon);
    S2Polygon decodedPolygon = decode(encoded);
    assertEquals(polygon, decodedPolygon);
    assertTrue(polygon.boundaryEquals(decodedPolygon));
    assertTrue(polygon.getRectBound().equals(decodedPolygon.getRectBound()));
  }

  @Test
  public void testCoderFast() {
    assertEquals(POLYGON, roundtrip(S2Polygon.FAST_CODER, POLYGON));
  }

  @Test
  public void testCoderCompact() {
    assertEquals(POLYGON, roundtrip(S2Polygon.COMPACT_CODER, POLYGON));
  }

  @Test
  public void testCoderShapeFast() {
    assertShapesEqual(POLYGON.shape(), roundtrip(S2Polygon.Shape.FAST_CODER, POLYGON.shape()));
  }

  @Test
  public void testCoderShapeCompact() {
    assertShapesEqual(
        SNAPPED_POLYGON.shape(), roundtrip(S2Polygon.Shape.COMPACT_CODER, SNAPPED_POLYGON.shape()));
  }

  // Make sure we've set things up correctly.
  @Test
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

  @Test
  public void testOverlapFractions() {
    double bOverlapA = S2Polygon.getOverlapFraction(overlapA, overlapB);
    double aOverlapB = S2Polygon.getOverlapFraction(overlapB, overlapA);
    assertEquals(0.25, bOverlapA, 0.01);
    assertEquals(0.5, aOverlapB, 0.01);

    // Empty polygons are defined as overlapping.
    S2Polygon a = makePolygon("");
    S2Polygon b = makePolygon("");
    aOverlapB = S2Polygon.getOverlapFraction(a, b);
    bOverlapA = S2Polygon.getOverlapFraction(b, a);
    assertEquals(1.0, aOverlapB, 0.01);
    assertEquals(1.0, bOverlapA, 0.01);

    b = makePolygon(OVERLAP_3);
    assertAlmostEquals(1.0, S2Polygon.getOverlapFraction(a, b));
    assertAlmostEquals(0.0, S2Polygon.getOverlapFraction(b, a));

    a = makePolygon(OVERLAP_4);
    assertDoubleNear(0.5, S2Polygon.getOverlapFraction(a, b), 1e-14);
    assertDoubleNear(0.5, S2Polygon.getOverlapFraction(b, a), 1e-14);

    // Two nearly identical polygons except that one is inverted.
    S2Polygon n =
        makeVerbatimPolygon(
            "33.596254:-112.1686194, 33.6106261:-112.1683897, 33.6106433:-112.1597676,"
                + " 33.6106034:-112.1510097, 33.5963451:-112.1511747");
    S2Polygon x =
        makeVerbatimPolygon(
            "33.596254:-112.1686194, 33.5963451:-112.1511747, 33.6106034:-112.1510097,"
                + " 33.6106433:-112.1597676, 33.6106261:-112.1683897");
    assertDoubleNear(0.0, S2Polygon.getOverlapFraction(n, x), 1e-9);
    assertDoubleNear(0.0, S2Polygon.getOverlapFraction(x, n), 1e-9);
  }

  @Test
  public void testEmptyAndFull() {
    assertTrue(empty.isEmpty());
    assertFalse(full.isEmpty());
    assertFalse(empty.isFull());
    assertTrue(full.isFull());

    checkNestedPair(empty, empty);
    checkNestedPair(full, empty);
    checkNestedPair(full, full);
  }

  @Test
  public void testIntersectionOverUnion() {
    double iouAb = S2Polygon.getIntersectionOverUnion(overlapA, overlapB);
    assertEquals(0.2, iouAb, 0.01);

    assertAlmostEquals(0.0, S2Polygon.getIntersectionOverUnion(empty, empty));

    S2Polygon b = makeVerbatimPolygon(OVERLAP_3);
    assertAlmostEquals(0.0, S2Polygon.getIntersectionOverUnion(b, empty));
  }

  @Test
  public void testOriginNearPole() {
    // S2Polygon operations are more efficient if S2.origin() is near a pole.
    // (Loops that contain a pole tend to have very loose bounding boxes because they span the full
    // longitude range.  S2Polygon canonicalizes all loops so that they don't contain S2.origin(),
    // thus by placing S2.origin() near a pole we minimize the number of canonical loops which
    // contain that pole.)
    assertTrue(S2LatLng.latitude(S2.origin()).degrees() >= 80);
  }

  private static class TestCase {
    final String a;
    final String b;
    final String aAndB;
    final String aOrB;
    final String aMinusB;
    final String aXorB;

    TestCase(String a, String b, String aAndB, String aOrB, String aMinusB, String aXorB) {
      this.a = a;
      this.b = b;
      this.aAndB = aAndB;
      this.aOrB = aOrB;
      this.aMinusB = aMinusB;
      this.aXorB = aXorB;
    }
  }

  private static final TestCase[] testCases = {
    // Two triangles that share an edge.
    new TestCase(
        /* a= */ "4:2, 3:1, 3:3;",
        /* b= */ "3:1, 2:2, 3:3;",
        /* aAndB= */ "",
        /* aOrB= */ "4:2, 3:1, 2:2, 3:3;",
        /* aMinusB= */ "4:2, 3:1, 3:3;",
        /* aXorB= */ "4:2, 3:1, 2:2, 3:3;"),
    // Two vertical bars and a horizontal bar connecting them.
    new TestCase(
        "0:0, 0:2, 3:2, 3:0;   0:3, 0:5, 3:5, 3:3;", // a
        "1:1, 1:4, 2:4, 2:1;", // b
        "1:1, 1:2, 2:2, 2:1;   1:3, 1:4, 2:4, 2:3;", // aAndB
        "0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0;", // aOrB
        "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;   " // aMinusB
            + "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;", //
        "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;   " // aXorB
            + "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;   " //
            + "1:2, 1:3, 2:3, 2:2"),
    // Two vertical bars and two horizontal bars centered around S2.origin().
    new TestCase(
        "1:88, 1:93, 2:93, 2:88;   -1:88, -1:93, 0:93, 0:88;", // a
        "-2:89, -2:90, 3:90, 3:89;   -2:91, -2:92, 3:92, 3:91;", // b
        "1:89, 1:90, 2:90, 2:89;   1:91, 1:92, 2:92, 2:91;   " // aAndB
            + "-1:89, -1:90, 0:90, 0:89;   -1:91, -1:92, 0:92, 0:91;",
        "-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92," // aOrB
            + "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91,"
            + "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88;   "
            + "0:90, 0:91, 1:91, 1:90;",
        "1:88, 1:89, 2:89, 2:88;   1:90, 1:91, 2:91, 2:90;   " // aMinusB
            + "1:92, 1:93, 2:93, 2:92;   -1:88, -1:89, 0:89, 0:88;   " //
            + "-1:90, -1:91, 0:91, 0:90;   -1:92, -1:93, 0:93, 0:92;", //
        "1:88, 1:89, 2:89, 2:88;   -1:88, -1:89, 0:89, 0:88;   " // aXorB
            + "1:90, 1:91, 2:91, 2:90;   -1:90, -1:91, 0:91, 0:90;   " //
            + "1:92, 1:93, 2:93, 2:92;   -1:92, -1:93, 0:93, 0:92;   " //
            + "-2:89, -2:90, -1:90, -1:89;   -2:91, -2:92, -1:92, -1:91;   " //
            + "0:89, 0:90, 1:90, 1:89;   0:91, 0:92, 1:92, 1:91;   " //
            + "2:89, 2:90, 3:90, 3:89;   2:91, 2:92, 3:92, 3:91;" //
        ),
    // Two interlocking square doughnuts centered around -S2.origin().
    new TestCase(
        "-1:-93, -1:-89, 3:-89, 3:-93;   0:-92, 0:-90, 2:-90, 2:-92;", // a
        "-3:-91, -3:-87, 1:-87, 1:-91;   -2:-90, -2:-88, 0:-88, 0:-90;", // b
        "-1:-91, -1:-90, 0:-90, 0:-91;   0:-90, 0:-89, 1:-89, 1:-90;", // aAndB
        "-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93;   " // aOrB
            + "0:-92, 0:-91, 1:-91, 1:-90, 2:-90, 2:-92;   "
            + "-2:-90, -2:-88, 0:-88, 0:-89, -1:-89, -1:-90;",
        "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, " // aMinusB
            + "1:-90, 1:-89, 3:-89, 3:-93;   " //
            + "-1:-90, -1:-89, 0:-89, 0:-90;", //
        "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, " // aXorB
            + "1:-90, 1:-89, 3:-89, 3:-93;   "
            + "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88, "
            + "-2:-88, -2:-90, -1:-90, -1:-91;   "
            + "-1:-90, -1:-89, 0:-89, 0:-90;   "
            + "1:-91, 0:-91, 0:-90, 1:-90;"),
    // An incredibly thin triangle intersecting a square, such that the two intersection points of
    // the triangle with the square are identical. This results in a degenerate loop that needs to
    // be handled correctly.
    new TestCase(
        "10:44, 10:46, 12:46, 12:44;", // a
        "11:45, 89:45.00000000000001, 90:45;", // b
        "", // aAndB - Empty intersection!
        // Original square with extra vertex, and triangle disappears (due to default
        // vertexMergeRadius of S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE).
        "10:44, 10:46, 12:46, 12:45, 12:44;", // aOrB
        "10:44, 10:46, 12:46, 12:45, 12:44;", // aMinusB
        "10:44, 10:46, 12:46, 12:45.001774937, 12:44;" // aXorB
        )
  };

  @Test
  public void testOperations() {
    // The error margin used when comparing the actual and expected results of constructive geometry
    // operations. The intersections in the expected data were computed in lat-lng space, while the
    // actual intersections are computed using geodesics. The error due to this depends on the
    // length and direction of the line segment being intersected, and how close the intersection is
    // to the endpoints of the segment. The worst case is for a line segment between two points at
    // the same latitude, where the intersection point is in the middle of the segment. In this case
    // the error is approximately {@code (p * t^2) / 8}, where {@code p} is the absolute latitude in
    // radians, {@code t} is the longitude difference in radians, and both {@code p} and {@code t}
    // are small. The test cases all have small latitude and longitude differences. If {@code p} and
    // {@code t} are converted to degrees, the following error bound is valid as long as {@code (p *
    // t^2 < 150)}.
    final double operationsMaxError = 1e-4;

    S2Polygon farSouth = new S2Polygon();
    farSouth.initToIntersection(farH, southH);
    checkEqual(farSouth, farHSouthH, 1e-15);

    for (int testNumber = 0; testNumber < testCases.length; testNumber++) {
      TestCase test = testCases[testNumber];
      logger.fine("Polygon operation test case " + testNumber);
      S2Polygon a = makeVerbatimPolygon(test.a);
      S2Polygon b = makeVerbatimPolygon(test.b);
      S2Polygon expectedAAndB = makeVerbatimPolygon(test.aAndB);
      S2Polygon expectedAOrB = makeVerbatimPolygon(test.aOrB);
      S2Polygon expectedAMinusB = makeVerbatimPolygon(test.aMinusB);
      S2Polygon expectedAXorB = makePolygon(test.aXorB);

      S2Polygon aAndB = new S2Polygon();
      S2Polygon aOrB = new S2Polygon();
      S2Polygon aMinusB = new S2Polygon();
      S2Polygon aXorB = new S2Polygon();

      aAndB.initToIntersection(a, b);
      checkEqual(aAndB, expectedAAndB, operationsMaxError);
      aOrB.initToUnion(a, b);
      checkEqual(aOrB, expectedAOrB, operationsMaxError);
      checkUnion(a, b);
      aMinusB.initToDifference(a, b);
      checkEqual(aMinusB, expectedAMinusB, operationsMaxError);
      aXorB.initToSymmetricDifference(a, b);
      checkEqual(aXorB, expectedAXorB, operationsMaxError);
    }
  }

  @Test
  public void testIntersectionSnapFunction() {
    // This tests that an intersection point is rounded to the nearest allowable vertex position
    // (using E0 coordinates, i.e. integer lat/lng values).
    S2Polygon a = makePolygon("0:0, 0:10, 1:10, 1:0");
    S2Polygon b = makePolygon("0:0, 0:10, 3:0");
    S2Polygon expected = makePolygon("0:0, 0:10, 1:7, 1:0");
    S2Polygon actual = new S2Polygon();
    actual.initToIntersection(a, b, new IntLatLngSnapFunction(0)); // E0 coords
    checkEqual(expected, actual);
  }

  @Test
  public void testIntersectionPreservesLoopOrder() {
    S2Polygon a = makePolygon("0:0, 0:10, 10:10, 10:0");
    S2Polygon b = makePolygon("1:1, 1:9, 9:5; 2:2, 2:8, 8:5");
    S2Polygon actual = new S2Polygon();
    actual.initToIntersection(a, b);
    assertEquals(S2TextFormat.toString(b), S2TextFormat.toString(actual));
  }

  /**
   * Verifies that the bounding rectangle optimization in initToIntersection() resets the result
   * polygon to be empty.
   */
  @Test
  public void testEmptyIntersectionClearsResult() {
    // The bounding rectangles of these two polygons do not intersect.
    S2Polygon a = makePolygon("0:0, 0:1, 1:0");
    S2Polygon b = makePolygon("3:3, 3:4, 4:3");

    // Initialize the result polygon to be non-empty, then verify that computing the intersection
    // clears the result.
    S2Polygon result = makePolygon("0:0, 0:1, 1:0");
    result.initToIntersection(a, b);
    assertTrue(result.isEmpty());

    // Repeat with the version of InitToIntersection that allows error reporting.
    S2Error error = new S2Error();
    result = makePolygon("0:0, 0:1, 1:0");
    assertTrue(result.initToIntersection(a, b, new IdentitySnapFunction(S1Angle.ZERO), error));
    assertTrue(result.isEmpty());
  }

  @Test
  public void testFullAndEmptyLoops() {
    S2Error error = new S2Error();

    // Canonical full and empty loops result in valid full and empty polygons.
    S2Polygon fullPoly = S2Polygon.fromLoops(ImmutableList.of(S2Loop.full()));
    assertTrue(fullPoly.isFull());
    assertTrue(fullPoly.isValid());
    assertEquals(1, fullPoly.numLoops());  // A full polygon has one full loop.

    S2Polygon emptyPoly = S2Polygon.fromLoops(ImmutableList.of(S2Loop.empty()));
    assertTrue(emptyPoly.isEmpty());
    assertTrue(emptyPoly.isValid());
    assertEquals(0, emptyPoly.numLoops());  // An empty polygon has no loops.

    // Points in the northern and southern hemispheres, not the canonical vertices for empty/full.
    S2Point pNorth = S2LatLng.fromDegrees(49, -112).toPoint();
    S2Point pSouth = S2LatLng.fromDegrees(-49, -112).toPoint();
    // Loops with those single points.
    S2Loop oneVertexLoopNorth = new S2Loop(ImmutableList.of(pNorth));
    S2Loop oneVertexLoopSouth = new S2Loop(ImmutableList.of(pSouth));

    // The one in the north is considered empty, the one in the south is considered full by S2Loop.
    assertTrue(oneVertexLoopNorth.isEmpty());
    assertTrue(oneVertexLoopSouth.isFull());

    // This "empty" loop isn't canonical, so it is added to an S2Polygon but considered an error.
    ArrayList<S2Loop> loops = new ArrayList<>();
    loops.add(oneVertexLoopNorth);
    S2Polygon polygonA = new S2Polygon(loops);
    assertEquals(1, polygonA.getLoops().size());
    assertTrue(polygonA.findValidationError(error));
    assertEquals(S2Error.Code.LOOP_NOT_ENOUGH_VERTICES, error.code());
    // Although the non-canonical single-vertex loop is considered empty by S2Loop, the S2Polygon
    // containing it is not empty.
    assertFalse(polygonA.isEmpty());

    error.clear();
    loops.clear();

    // Likewise, a non-canonical single-vertex loop in the southern hemisphere is added but
    // considered to be an error.
    loops.add(oneVertexLoopSouth);
    S2Polygon polygonB = new S2Polygon(loops);
    loops.add(oneVertexLoopSouth);
    assertEquals(1, polygonB.getLoops().size());
    assertTrue(polygonB.findValidationError(error));
    assertEquals(S2Error.Code.LOOP_NOT_ENOUGH_VERTICES, error.code());
    // Although the non-canonical single-vertex loop is considered full by S2Loop, the S2Polygon
    // containing it is not full.
    assertFalse(polygonB.isFull());
  }

  // Verifies that S2Polygon does not destroy or replace S2Loops, so callers can rely on using
  // references for loop identity.
  @Test
  public void testLoopPointers() {
    List<S2Loop> loops = new ArrayList<>();
    loops.add(S2TextFormat.makeLoopOrDie("4:4, 4:6, 6:6, 6:4"));
    loops.add(S2TextFormat.makeLoopOrDie("3:3, 3:7, 7:7, 7:3"));
    loops.add(S2TextFormat.makeLoopOrDie("2:2, 2:8, 8:8, 8:2"));
    loops.add(S2TextFormat.makeLoopOrDie("1:1, 1:9, 9:9, 9:1"));
    loops.add(S2TextFormat.makeLoopOrDie("10:10, 15:15, 20:10"));
    loops.add(S2TextFormat.makeLoopOrDie("-1:-1, -9:-1, -9:-9, -1:-9"));
    loops.add(S2TextFormat.makeLoopOrDie("-5:-5, -6:-5, -6:-6, -5:-6"));

    Set<S2Loop> loopReferences = Sets.newIdentityHashSet();
    loopReferences.addAll(loops);
    S2Polygon polygon = new S2Polygon(loops);

    // Check that loop pointers didn't change (but could've gotten reordered).
    assertEquals(loopReferences.size(), polygon.numLoops());
    for (int i = 0; i < polygon.numLoops(); i++) {
      assertTrue(loopReferences.contains(polygon.loop(i)));
    }
  }

  // The "Bug" tests are regression tests from previous versions of the algorithm, ported from C++.

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug1() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10531193335759943, -0.80522214810955617, 0.58354664670985534),
                p(-0.10531194840431297, -0.80522215192439039, 0.58354663873039425),
                p(-0.10531192794033867, -0.80522217497559767, 0.58354661061568747),
                p(-0.10531191284235047, -0.80522217121852058, 0.58354661852470402)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10531174240075937, -0.80522236320875284, 0.58354638436119843),
                p(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
                p(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508), // A
                p(-0.10531191288915481, -0.80522217116640804, 0.5835466185881667), // B
                p(-0.10531191288915592, -0.8052221711664066, 0.58354661858816803), // B
                p(-0.10531192039151964, -0.80522217309710431, 0.58354661457010204), // A
                p(-0.10531192794033779, -0.80522217497559878, 0.58354661061568636),
                p(-0.1053117575499668, -0.80522236690813498, 0.58354637652254981)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Given edges do not form loops (indegree != outdegree)
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug2() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10618951389689163, -0.80546461394606728, 0.58305277875939732),
                p(-0.10618904764039243, -0.8054645437464607, 0.58305296065497536),
                p(-0.10618862643748632, -0.80546451917975415, 0.58305307130470341),
                p(-0.10617606798507535, -0.80544758470051458, 0.58307875187433833)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10618668131028208, -0.80544613076731553, 0.58307882755616247),
                p(-0.10618910658843225, -0.80546454998744921, 0.58305294129732887),
                p(-0.10618904764039225, -0.80546454374646081, 0.58305296065497536),
                p(-0.10618898834264634, -0.80546453817003949, 0.58305297915823251)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Given edges do not form loops (indegree != outdegree)
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug3() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10703494861068318, -0.80542232562508131, 0.58295659972299307),
                p(-0.10703494998722708, -0.80542232255642865, 0.58295660370995028),
                p(-0.10703495367938694, -0.80542232008675829, 0.58295660644418046),
                p(-0.10703495869785147, -0.80542231887781635, 0.58295660719304865),
                p(-0.10703496369792719, -0.80542231925353791, 0.58295660575589636),
                p(-0.10703496733984781, -0.80542232111324863, 0.58295660251780734),
                p(-0.10703496864776367, -0.80542232395864055, 0.58295659834642488),
                p(-0.10703496727121976, -0.80542232702729322, 0.58295659435946767),
                p(-0.10703496357905991, -0.80542232949696357, 0.5829565916252375),
                p(-0.10703495856059538, -0.80542233070590552, 0.58295659087636931),
                p(-0.10703495356051966, -0.80542233033018396, 0.58295659231352159),
                p(-0.10703494991859903, -0.80542232847047324, 0.58295659555161061)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10703494861068762, -0.80542232562508098, 0.58295659972299274),
                p(-0.10703494998723152, -0.80542232255642832, 0.58295660370994995),
                p(-0.10703495367939138, -0.80542232008675796, 0.58295660644418013),
                p(-0.10703495869785591, -0.80542231887781601, 0.58295660719304832),
                p(-0.10703496369793163, -0.80542231925353758, 0.58295660575589603),
                p(-0.10703496733985225, -0.8054223211132483, 0.58295660251780701),
                p(-0.10703496864776811, -0.80542232395864022, 0.58295659834642455),
                p(-0.1070349672712242, -0.80542232702729288, 0.58295659435946734),
                p(-0.10703496357906438, -0.80542232949696346, 0.58295659162523727),
                p(-0.10703495856059982, -0.80542233070590519, 0.58295659087636897),
                p(-0.1070349535605241, -0.80542233033018362, 0.58295659231352126),
                p(-0.10703494991860348, -0.8054223284704729, 0.58295659555161028)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Given edges do not form loops (indegree != outdegree)
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug4() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10667065556339718, -0.80657502337947207, 0.58142764201754193),
                p(-0.10667064691895933, -0.80657502457251051, 0.58142764194845853),
                p(-0.10667064691930939, -0.80657502457246333, 0.58142764194845975),
                p(-0.10667065556339746, -0.80657502337947395, 0.5814276420175396),
                p(-0.10667077559567185, -0.80657589269604968, 0.58142641405029793),
                p(-0.10667077059539463, -0.80657589232162286, 0.58142641548708696),
                p(-0.10667063827452879, -0.80657502576554818, 0.58142764187937435),
                p(-0.10667063169531328, -0.80657498170361974, 0.58142770421053058),
                p(-0.10667064898418178, -0.8065749793175444, 0.58142770434869739)),
            ImmutableList.of(
                p(-0.10667064691897719, -0.80657502457250896, 0.58142764194845697),
                p(-0.10667063827452879, -0.80657502576554818, 0.58142764187937435),
                p(-0.10667064691861985, -0.80657502457255736, 0.58142764194845586)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10667064691896312, -0.80657502457251107, 0.58142764194845697),
                p(-0.10667064691896297, -0.80657502457251007, 0.58142764194845853),
                p(-0.10667064033974753, -0.80657498051058207, 0.58142770427961399),
                p(-0.10667064076268165, -0.80657498045444342, 0.58142770427989865),
                p(-0.10667051785242875, -0.80657409963649807, 0.58142894872603923),
                p(-0.1066707756642685, -0.80657588679775971, 0.58142642222003538)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Loop 1: Edge 1 crosses edge 3
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug5() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10574444273627338, -0.80816264611829447, 0.57938868667714882),
                p(-0.10574444845633162, -0.80816268110163325, 0.57938863683652475),
                p(-0.10574444825833453, -0.80816268112970524, 0.57938863683350494),
                p(-0.10574444253827629, -0.80816264614636646, 0.57938868667412902),
                p(-0.10574408792844124, -0.80816047738475361, 0.57939177648757634),
                p(-0.10574408812643833, -0.80816047735668162, 0.57939177649059592)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.1057440881264381, -0.80816047735668017, 0.57939177649059825),
                p(-0.10574408802743954, -0.80816047737071606, 0.57939177648908835),
                p(-0.10574408812649677, -0.8081604773570521, 0.57939177649006868),
                p(-0.10574408812649701, -0.80816047735705354, 0.57939177649006646),
                p(-0.10574408802703171, -0.80816047737077379, 0.57939177648908202),
                p(-0.10574408792844098, -0.80816047738475194, 0.57939177648757834),
                p(-0.10574408792838257, -0.80816047738438168, 0.5793917764881058),
                p(-0.1057440879283823, -0.80816047738438002, 0.57939177648810791),
                p(-0.10574407993470979, -0.80816042849578984, 0.57939184613891748),
                p(-0.10574408013270691, -0.80816042846771807, 0.57939184614193739)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Loop 0 edge 8 crosses loop 1 edge 0
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug6() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10618849949725141, -0.80552159562437586, 0.58297423747304822),
                p(-0.10618849959636036, -0.80552159561106063, 0.58297423747339361),
                p(-0.10618849949722192, -0.80552159562415893, 0.5829742374733532),
                p(-0.10618834540082922, -0.80552043435619214, 0.58297587011440333),
                p(-0.10618834559910612, -0.80552043432999554, 0.58297587011448437),
                p(-0.10618849969546933, -0.80552159559774539, 0.58297423747373922),
                p(-0.10618849969546955, -0.80552159559774716, 0.582974237473737),
                p(-0.10618849969549882, -0.80552159559796233, 0.58297423747343424),
                p(-0.10618849959710704, -0.80552159561096182, 0.58297423747339394),
                p(-0.10618849949725161, -0.80552159562437742, 0.58297423747304589)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10618856154870562, -0.80552206324314812, 0.58297358004005528),
                p(-0.10618849949722212, -0.80552159562416048, 0.58297423747335086),
                p(-0.10618849969549901, -0.80552159559796388, 0.58297423747343191),
                p(-0.10618856174698249, -0.8055220632169513, 0.58297358004013622),
                p(-0.10618857104277038, -0.80552213326985989, 0.58297348155149287),
                p(-0.10618857084449349, -0.80552213329605649, 0.58297348155141182)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Loop 0 edge 0 crosses loop 1 edge 4
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug7() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10651728339354898, -0.80806023027835039, 0.57938996589599123),
                p(-0.10651728368541774, -0.80806023024121265, 0.57938996589412783),
                p(-0.10651743884289547, -0.80806147782022508, 0.5793881973990701),
                p(-0.1065172793067945, -0.80806153133252501, 0.5793881520963412),
                p(-0.10651707335497011, -0.80806158532388361, 0.57938811465868356),
                p(-0.10651593657771009, -0.80806167503227055, 0.57938819853274059),
                p(-0.10651567693742285, -0.80806182530835402, 0.57938803667826444),
                p(-0.10651496089498214, -0.80806213485510237, 0.57938773659696563),
                p(-0.10651453461919227, -0.80806229235522298, 0.57938759530083062),
                p(-0.10651448583749658, -0.80806230280784852, 0.57938758969074455),
                p(-0.10651428153471061, -0.80806061225022852, 0.57938998503506256),
                p(-0.10651428161845182, -0.8080606122395747, 0.57938998503452654),
                p(-0.10651427761078044, -0.80806057978063328, 0.57939003104095654),
                p(-0.10651427761077951, -0.80806057978062562, 0.57939003104096709),
                p(-0.10651387099203104, -0.8080572864940091, 0.5793946988282096),
                p(-0.10651387099202798, -0.80805728649398445, 0.57939469882824468),
                p(-0.10651386444607201, -0.80805723347699177, 0.57939477397218053),
                p(-0.10651386444607169, -0.8080572334769891, 0.57939477397218409),
                p(-0.106513765993723, -0.80805643609199118, 0.57939590414857456),
                p(-0.10651376671438624, -0.8080564359989727, 0.57939590414581921),
                p(-0.10651368187839319, -0.80805575808078389, 0.57939686520139033),
                p(-0.10651465698432123, -0.80805552598235797, 0.57939700963750851),
                p(-0.1065149024434091, -0.80805548225095913, 0.57939702550292815),
                p(-0.10651504788182964, -0.80805555533715756, 0.5793968968362615),
                p(-0.10651511658091152, -0.80805559604710031, 0.57939682743066534),
                p(-0.10651517919248171, -0.80805562751022852, 0.57939677204023521),
                p(-0.10651528575974038, -0.80805561374213786, 0.57939677165077275),
                p(-0.10651648823358072, -0.80805539171529139, 0.57939686023850034),
                p(-0.10651666406737116, -0.80805537863686483, 0.57939684615295572),
                p(-0.10651674780673852, -0.80805605121551227, 0.57939589274577097),
                p(-0.10651674667750256, -0.80805605136137271, 0.57939589274994641),
                p(-0.10651678418140036, -0.80805634336988752, 0.57939547860450136),
                p(-0.10651680240261223, -0.80805648524178364, 0.57939527739240138),
                p(-0.10651680240261237, -0.80805648524178486, 0.57939527739239993)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10651727337444802, -0.80806023111043901, 0.57938996657744879),
                p(-0.10651727440799089, -0.80806022882029649, 0.57938996958144073),
                p(-0.10651679374955145, -0.80805648637258243, 0.57939527740611751),
                p(-0.10651677552833975, -0.80805634450068775, 0.57939547861821594),
                p(-0.10651673802444192, -0.80805605249217261, 0.57939589276366099),
                p(-0.10651674651102909, -0.80805605138312775, 0.5793958927502102),
                p(-0.10651673915225639, -0.80805605233507238, 0.57939589277542292),
                p(-0.10651665541288889, -0.80805537975642383, 0.57939684618260878),
                p(-0.10651667272185343, -0.80805537751730583, 0.57939684612330267),
                p(-0.1065167564612207, -0.8080560500959526, 0.57939589271611924),
                p(-0.1065167553320342, -0.80805605024202609, 0.57939589271998793),
                p(-0.10651679283446101, -0.80805634223908773, 0.57939547859078699),
                p(-0.10651681105567287, -0.80805648411098374, 0.57939527737868723),
                p(-0.10651680240318392, -0.80805648524170914, 0.5793952773924006),
                p(-0.10651680240261234, -0.80805648524178475, 0.57939527739239982),
                p(-0.1065168110556733, -0.80805648411098718, 0.57939527737868224),
                p(-0.10651729169518892, -0.80806022641135866, 0.57938996976297907),
                p(-0.10651729210462238, -0.80806022661896348, 0.579389969398166),
                p(-0.1065172934126499, -0.80806022944626155, 0.57938996521453356),
                p(-0.10651729203606744, -0.80806023249651726, 0.57938996121349717),
                p(-0.1065172883437291, -0.80806023495241674, 0.57938995846713126),
                p(-0.10651728332499401, -0.80806023615590394, 0.5793899577113224),
                p(-0.10651727832462815, -0.80806023578450537, 0.57938995914858893),
                p(-0.10651727468247554, -0.80806023393773707, 0.57938996239381635)),
            ImmutableList.of(
                p(-0.10651680240204828, -0.80805648524185858, 0.57939527739240082),
                p(-0.10651679861449742, -0.80805648573682254, 0.57939527739840524),
                p(-0.10651680240261419, -0.80805648524178353, 0.57939527739240138)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Loop 0: Edge 33 crosses edge 35
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug8() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10703872198218529, -0.80846112144645677, 0.57873424566545062),
                p(-0.10703872122182066, -0.80846111957630917, 0.57873424841857957),
                p(-0.10703873813385757, -0.80846111582010538, 0.57873425053786276),
                p(-0.1070387388942222, -0.80846111769025297, 0.57873424778473381),
                p(-0.10703873050793056, -0.80846111955286837, 0.57873424673382978),
                p(-0.1070387388942227, -0.80846111769025419, 0.57873424778473193),
                p(-0.10703919382477994, -0.80846223660916783, 0.57873260056976505),
                p(-0.10703917691274406, -0.80846224036537406, 0.57873259845047831)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10703917691274355, -0.80846224036537273, 0.57873259845047997),
                p(-0.1070391853685064, -0.8084622384873289, 0.57873259951008804),
                p(-0.10703919381027188, -0.80846223657409677, 0.57873260062144094),
                p(-0.10703919381027233, -0.80846223657409788, 0.57873260062143939),
                p(-0.10703918536876245, -0.80846223848727206, 0.57873259951012024),
                p(-0.10703919382478132, -0.80846223660917116, 0.57873260056976017),
                p(-0.10703957146434441, -0.80846316542623331, 0.57873123320737097),
                p(-0.10703955455230836, -0.8084631691824391, 0.57873123108808489)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(a));
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(b));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    //  Loop 1: Edge 1 crosses edge 3
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(c));
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  @Ignore("I suspect this test fails because of b/400510172")
  public void testBug9() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10639937100501309, -0.80810205676564995, 0.57935329437301375),
                p(-0.10639937101137514, -0.80810205688156922, 0.57935329421015713),
                p(-0.10639937101137305, -0.80810205688156944, 0.57935329421015713),
                p(-0.106399371005011, -0.80810205676565017, 0.57935329437301375)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10639937099530022, -0.8081020567669569, 0.57935329437297489),
                p(-0.10639937102108385, -0.80810205688026293, 0.5793532942101961),
                p(-0.10639937102108181, -0.80810205688026326, 0.5793532942101961),
                p(-0.10639937099529816, -0.80810205676695701, 0.57935329437297478)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Given edges do not form loops (indegree != outdegree)
    assertFalse(
        "\nS2Polygon A: "
            + S2TextFormat.toString(a)
            + "\nS2Polygon B: "
            + S2TextFormat.toString(b)
            + "\nS2Polygon C: "
            + S2TextFormat.toString(c),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug10() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10592889932808099, -0.80701394501854917, 0.58095400922339757),
                p(-0.10592787800899696, -0.8070140771413753, 0.58095401191158469),
                p(-0.1059270044681431, -0.80701419014619669, 0.58095401421031945),
                p(-0.10592685562894633, -0.80701420940058122, 0.58095401460194696),
                p(-0.10592685502239066, -0.80701420947920588, 0.58095401460332308),
                p(-0.10592681668594067, -0.80701421444855337, 0.5809540146902914),
                p(-0.10592586497682262, -0.8070143378130904, 0.58095401684902004),
                p(-0.10592586434121586, -0.80701433789547994, 0.58095401685046155),
                p(-0.10592585898876766, -0.80701428569270217, 0.58095409034224832),
                p(-0.10592585898876755, -0.80701428569270128, 0.58095409034224987),
                p(-0.10592571912106936, -0.8070129215545373, 0.58095601078971082),
                p(-0.10592571912106795, -0.80701292155452331, 0.58095601078973025),
                p(-0.10592546626664477, -0.80701045545315664, 0.58095948256783148),
                p(-0.10592546630689463, -0.80701045544795602, 0.58095948256771723),
                p(-0.10592538513536764, -0.80700975616910509, 0.58096046873415197),
                p(-0.10592564439344856, -0.80700971612782446, 0.58096047708524956),
                p(-0.1059267844512099, -0.80700966174311928, 0.58096034476466896),
                p(-0.10592686088387009, -0.80700965393230761, 0.58096034167862642),
                p(-0.10592691331665709, -0.80700961093727019, 0.58096039184274961),
                p(-0.10592705773734933, -0.80700947507458121, 0.58096055423665138),
                p(-0.10592721940752658, -0.80700934249808198, 0.58096070892049412),
                p(-0.10592756003095027, -0.80700933299293154, 0.58096066001769275),
                p(-0.10592832507751106, -0.80700935762745474, 0.58096048630521868),
                p(-0.1059284165295875, -0.80701007424011018, 0.58095947418602778),
                p(-0.10592841614913188, -0.80701007428931704, 0.58095947418704452),
                p(-0.10592864947042728, -0.8070119434176124, 0.58095683523192998),
                p(-0.1059286884898481, -0.80701225600079662, 0.58095639390519271),
                p(-0.10592868927069989, -0.80701225581371527, 0.58095639402269295),
                p(-0.10592869427137827, -0.80701225619024619, 0.58095639258785126),
                p(-0.10592869791375134, -0.80701225804491505, 0.58095638934738025),
                p(-0.10592869922184817, -0.80701226088076483, 0.5809563851695615),
                p(-0.10592869922184843, -0.80701226088076705, 0.58095638516955805),
                p(-0.10592869784516552, -0.80701226393793402, 0.58095638117383475),
                p(-0.10592869415258396, -0.80701226639725276, 0.58095637843085768),
                p(-0.10592868991437976, -0.80701226741266929, 0.58095637779310561)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10592564460843924, -0.80700972122716552, 0.58096046996257766),
                p(-0.10592539435053176, -0.80700975987840939, 0.58096046190138972),
                p(-0.10592547496472972, -0.80701045435596641, 0.58095948250602925),
                p(-0.10592546630689462, -0.80701045544795591, 0.58095948256771723),
                p(-0.10592546630693271, -0.80701045544826022, 0.58095948256728758),
                p(-0.1059254749287661, -0.80701045440038255, 0.5809594824508878),
                p(-0.10592572778318898, -0.80701292050174633, 0.58095601067279068),
                p(-0.1059257191207934, -0.80701292155455673, 0.58095601078973391),
                p(-0.1059257194541381, -0.80701292151405679, 0.58095601078521419),
                p(-0.10592572778319062, -0.80701292050176254, 0.58095601067276803),
                p(-0.10592586765088864, -0.80701428463992497, 0.58095409022530931),
                p(-0.10592585899855227, -0.80701428569151201, 0.58095409034211776),
                p(-0.10592585898857355, -0.80701428569272593, 0.58095409034225098),
                p(-0.10592586765088888, -0.80701428463992686, 0.58095409022530675),
                p(-0.10592587247896063, -0.80701433172842685, 0.58095402393347073),
                p(-0.10592681605007616, -0.80701420941876889, 0.58095402179319922),
                p(-0.10592685438651758, -0.80701420444942229, 0.58095402170623067),
                p(-0.10592685499307326, -0.80701420437079774, 0.58095402170485466),
                p(-0.10592685562894634, -0.80701420940058122, 0.58095401460194696),
                p(-0.10592685499689927, -0.80701420437030225, 0.58095402170484534),
                p(-0.10592700383609792, -0.80701418511591771, 0.58095402131321794),
                p(-0.10592787737695626, -0.80701407211109533, 0.58095401901448296),
                p(-0.10592889869604118, -0.80701393998826909, 0.58095401632629584),
                p(-0.10592889996012077, -0.80701395004882903, 0.58095400212049919),
                p(-0.10592787864104941, -0.80701408217165349, 0.58095400480868631),
                p(-0.10592787800903029, -0.80701407714164064, 0.58095401191120999),
                p(-0.10592787864103763, -0.80701408217165482, 0.5809540048086862),
                p(-0.10592700510019466, -0.80701419517647521, 0.58095400710742118),
                p(-0.1059270044681431, -0.80701419014619669, 0.58095401421031934),
                p(-0.10592700510018833, -0.8070141951764761, 0.58095400710742118),
                p(-0.10592685626275877, -0.80701421443063182, 0.58095400749904391),
                p(-0.10592685565826369, -0.80701421450898914, 0.58095400750041526),
                p(-0.10592685502239063, -0.80701420947920566, 0.58095401460332308),
                p(-0.10592685565826078, -0.80701421450898947, 0.58095400750041526),
                p(-0.10592681732181129, -0.80701421947833718, 0.58095400758738369),
                p(-0.10592681668594069, -0.80701421444855348, 0.58095401469029151),
                p(-0.10592681732180521, -0.80701421947833796, 0.58095400758738369),
                p(-0.10592586561269894, -0.80701434284287321, 0.58095400974611222),
                p(-0.10592586497746249, -0.80701433781815202, 0.58095401684187198),
                p(-0.10592586561268771, -0.80701434284287465, 0.58095400974611222),
                p(-0.10592586497708102, -0.80701434292526464, 0.58095400974755396),
                p(-0.10592586434121586, -0.80701433789548005, 0.58095401685046166),
                p(-0.10592585567909471, -0.80701433894825569, 0.58095401696740323),
                p(-0.1059258503266465, -0.80701428674547793, 0.58095409045919011),
                p(-0.10592571045894811, -0.80701292260731206, 0.58095601090665361),
                p(-0.10592571912060067, -0.80701292155459425, 0.58095601078971715),
                p(-0.10592571878923682, -0.80701292159485349, 0.58095601079421),
                p(-0.10592571045894694, -0.80701292260730051, 0.58095601090666993),
                p(-0.10592545760452345, -0.80701045650593073, 0.58095948268477515),
                p(-0.10592545764454649, -0.80701045650106651, 0.58095948268423492),
                p(-0.10592537647753246, -0.80700975726109381, 0.58096046879584118),
                p(-0.10592538513536764, -0.80700975616910509, 0.58096046873415197),
                p(-0.10592538413784101, -0.80700975119062324, 0.58096047583161736),
                p(-0.10592564339592514, -0.80700971114934217, 0.58096048418271495),
                p(-0.10592564439344856, -0.80700971612782446, 0.58096047708524956),
                p(-0.10592564496449927, -0.80700971099098684, 0.58096048411668999),
                p(-0.10592678502227458, -0.80700965660628099, 0.58096035179610783),
                p(-0.10592678388014524, -0.80700966687995779, 0.58096033773323019)),
            ImmutableList.of(
                p(-0.10592585898876757, -0.80701428569270128, 0.58095409034224987),
                p(-0.10592585897888845, -0.80701428569390288, 0.58095409034238166),
                p(-0.1059258503266465, -0.80701428674547793, 0.58095409045919011)),
            ImmutableList.of(
                p(-0.10592546626664477, -0.80701045545315664, 0.58095948256783148),
                p(-0.10592546623958927, -0.8070104554564449, 0.58095948256819674),
                p(-0.10592546626662946, -0.80701045545303429, 0.580959482568004)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(a));
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(b));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Inconsistent loop orientations detected
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(c));
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug11() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10727349803435572, -0.80875763107088172, 0.57827631008375979),
                p(-0.10727349807040805, -0.80875763112192245, 0.57827631000568813),
                p(-0.10727349807040625, -0.80875763112192278, 0.57827631000568813)),
            ImmutableList.of(
                p(-0.1072729603486537, -0.80875606054879057, 0.57827860629945249),
                p(-0.10727299870478688, -0.80875633377729705, 0.57827821705818028),
                p(-0.10727299875560981, -0.80875633413933223, 0.57827821654242495),
                p(-0.10727309272230967, -0.80875700360375646, 0.57827726282438607),
                p(-0.10727318660000487, -0.80875767243400742, 0.57827631000742785),
                p(-0.10727349802669105, -0.80875763101356435, 0.57827631016534387),
                p(-0.10727349803435525, -0.80875763107087817, 0.57827631008376468),
                p(-0.10727349803435572, -0.80875763107088172, 0.57827631008375979),
                p(-0.1072734980420204, -0.80875763112819909, 0.57827631000217561),
                p(-0.10727318657570066, -0.80875767255391384, 0.57827630984423972),
                p(-0.10727318651657966, -0.80875767256177711, 0.57827630984420975),
                p(-0.10727318650891528, -0.80875767250445951, 0.57827630992579371),
                p(-0.10727318640981781, -0.80875767251785957, 0.57827630992543622),
                p(-0.10727309252411468, -0.80875700363055636, 0.57827726282367087),
                p(-0.10727299855741491, -0.8087563341661328, 0.57827821654170874),
                p(-0.10727299850659211, -0.8087563338040985, 0.57827821705746318),
                p(-0.10727296014242577, -0.80875606051836801, 0.57827860638025652),
                p(-0.10727296024152315, -0.80875606050496729, 0.57827860638061501),
                p(-0.10727296023340849, -0.8087560604477102, 0.57827860646219797),
                p(-0.10727348576547496, -0.80875598914629976, 0.57827860869282954),
                p(-0.1072734857817042, -0.80875598926081438, 0.57827860852966395)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.1072734857735896, -0.80875598920355718, 0.5782786086112468),
                p(-0.10727348576547457, -0.80875598914629976, 0.57827860869282954),
                p(-0.10727839137361543, -0.80875532356817348, 0.57827862950694298),
                p(-0.10727839137881608, -0.80875532356471602, 0.57827862951081388),
                p(-0.10727839143632178, -0.80875532355090063, 0.5782786295194674),
                p(-0.10727839149361706, -0.80875532355509905, 0.57827862950296649),
                p(-0.1072783915353497, -0.80875532357618651, 0.57827862946573261),
                p(-0.10727839154773799, -0.80875532360290581, 0.57827862942606567),
                p(-0.10727848921795155, -0.80875531035110082, 0.57827862984032907),
                p(-0.1072784892332832, -0.80875531046514559, 0.57827862967798682),
                p(-0.10727971608197531, -0.8087551454635169, 0.57827863284376713),
                p(-0.10727986275126807, -0.80875539440654376, 0.57827825747332484),
                p(-0.10727959167812619, -0.80875599171505064, 0.57827747239052929),
                p(-0.10727974196569352, -0.80875625444235633, 0.57827707706958686),
                p(-0.10727993501555312, -0.80875677560355186, 0.57827631237878363),
                p(-0.10727870858143702, -0.80875693828645479, 0.57827631237896882),
                p(-0.1072787085493927, -0.80875693804871851, 0.5782763127174031),
                p(-0.10727615977928232, -0.80875727704955946, 0.57827631143112901),
                p(-0.10727615977915911, -0.80875727704957578, 0.57827631143112901),
                p(-0.10727349803435751, -0.80875763107088128, 0.57827631008375968),
                p(-0.10727349803435574, -0.80875763107088183, 0.57827631008375979),
                p(-0.10727318656803594, -0.80875767249659658, 0.57827630992582391),
                p(-0.10727318650891531, -0.80875767250445962, 0.57827630992579382),
                p(-0.10727309262321218, -0.80875700361715641, 0.57827726282402847),
                p(-0.10727299865651231, -0.80875633415273218, 0.57827821654206735),
                p(-0.10727299860568951, -0.80875633379069789, 0.57827821705782179),
                p(-0.10727296024152314, -0.80875606050496718, 0.57827860638061501)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Given edges do not form loops (indegree != outdegree)
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testBug12() {
    ImmutableList<List<S2Point>> aVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10772916872905106, -0.80699542608967267, 0.58064861015531188),
                p(-0.10772916892726483, -0.80699542606300401, 0.58064861015560143),
                p(-0.10772916892726613, -0.80699542606301333, 0.58064861015558844),
                p(-0.10772916872905235, -0.806995426089682, 0.58064861015529889)));
    ImmutableList<List<S2Point>> bVertices =
        ImmutableList.of(
            ImmutableList.of(
                p(-0.10772916872905348, -0.80699542608969022, 0.58064861015528724),
                p(-0.10772916892726496, -0.80699542606300489, 0.58064861015559999),
                p(-0.10772930108168739, -0.80699639165138115, 0.58064724364290399),
                p(-0.10772930088347589, -0.80699639167806647, 0.58064724364259113)));
    S2Polygon a = new S2Polygon(makeLoops(aVertices));
    S2Polygon b = new S2Polygon(makeLoops(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    // Given edges do not form loops (indegree != outdegree)
    assertFalse(
        "\nS2Polygon: " + S2TextFormat.toString(a) + "\nS2Polygon: " + S2TextFormat.toString(b),
        c.isEmpty());
  }

  /**
   * This tests polygon-polyline intersections where the polyline is a single edge of the polygon.
   */
  @Test
  public void testPolylineIntersectionForwardEdge() {
    for (int v = 0; v < 3; ++v) {
      polylineIntersectionSharedEdgeTest(cross1, v, 1);
      polylineIntersectionSharedEdgeTest(cross1SideHole, v, 1);
    }
  }

  /**
   * This tests polygon-polyline intersections where the polyline is a single reversed edge of the
   * polygon.
   */
  @Test
  public void testPolylineIntersectionReverseEdge() {
    for (int v = 0; v < 3; ++v) {
      polylineIntersectionSharedEdgeTest(cross1, v + 1, -1);
      polylineIntersectionSharedEdgeTest(cross1SideHole, v + 1, -1);
    }
  }

  /**
   * This tests polygon-polyline intersections. It covers the same edge cases as {@code
   * testOperations} and also adds some extra tests for shared edges.
   */
  @Test
  public void testMorePolylineIntersections() {
    // See comments in testOperations about the value of this constant.
    final double operationsMaxError = 1e-4;

    // This duplicates some of the tests in testOperations by converting the outline of polygon A to
    // a polyline then intersecting it with the polygon B. It then converts B to a polyline and
    // intersects it with A. It then feeds all of the results into a polygon builder and tests that
    // the output is equal to doing an intersection between A and B.

    for (int testNumber = 0; testNumber < testCases.length; testNumber++) {
      TestCase test = testCases[testNumber];
      logger.fine("Polyline intersection test case " + testNumber);
      S2Polygon a = makeVerbatimPolygon(test.a);
      S2Polygon b = makeVerbatimPolygon(test.b);
      S2Polygon expectedAAndB = makeVerbatimPolygon(test.aAndB);

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

      S2Builder builder = new S2Builder.Builder().build();
      S2PolygonLayer layer = new S2PolygonLayer();
      builder.startLayer(layer);

      for (S2Polyline polyLine : polylines) {
        builder.addPolyline(polyLine);
      }

      S2Error error = new S2Error();
      boolean built = builder.build(error);
      assertTrue(error.text(), built);
      S2Polygon aAndB = layer.getPolygon();
      checkEqual(aAndB, expectedAAndB, operationsMaxError);
    }
  }

  @Test
  public void testIsValidUnitLength() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = getConcentricLoops(1 + data.random(6), 3);
      List<S2Point> loop = loops.get(data.random(loops.size()));
      int index = data.random(loop.size());
      S2Point p = loop.get(index);
      switch (data.random(3)) {
        case 0:
          p = new S2Point(0, 0, 0);
          break;
        case 1:
          p = p.mul(1e-30 * pow(1e60, data.nextDouble()));
          break;
        default:
          p = new S2Point(Double.NaN, Double.NaN, Double.NaN);
          break;
      }
      loop.set(index, p);
      checkInvalid(loops, "unit length");
    }
  }

  @Test
  public void testIsValidVertexCount() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = Lists.newArrayList();
      List<S2Point> loopPoints = Lists.newArrayList();
      loopPoints.add(data.getRandomPoint());
      loopPoints.add(data.getRandomPoint());
      loops.add(loopPoints);
      checkInvalid(loops, S2Error.Code.LOOP_NOT_ENOUGH_VERTICES);
    }
  }

  @Test
  public void testIsValidDuplicateVertex() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = getConcentricLoops(1, 3);
      List<S2Point> loop = loops.get(0);
      int n = loop.size();
      int i = data.random(n);
      int j = data.random(n - 1);
      loop.set(i, loop.get(j + (j >= i ? 1 : 0)));
      checkInvalid(
          loops,
          ImmutableSet.of(
              // Duplicate caused a degenerate edges.
              S2Error.Code.DUPLICATE_VERTICES,
              // Duplicate caused the interior to flip.
              S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
              // Duplicate caused a duplicate polygon edge.
              S2Error.Code.OVERLAPPING_GEOMETRY));
    }
  }

  @Test
  public void testIsValidSelfIntersection() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      // Use multiple loops so that we can test both holes and shells. We need at least 5 vertices
      // so that the modified edges don't intersect any nested loops.
      List<List<S2Point>> loops = getConcentricLoops(1 + data.random(6), 5);
      List<S2Point> loop = loops.get(data.random(loops.size()));
      int n = loop.size();
      int i = data.random(n);
      Collections.swap(loop, i, (i + 1) % n);
      checkInvalid(
          loops,
          ImmutableSet.of(
              S2Error.Code.LOOP_SELF_INTERSECTION,
              S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
              S2Error.Code.OVERLAPPING_GEOMETRY));
    }
  }
  @Test
  public void testIsValidEmptyLoop() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = getConcentricLoops(data.random(5), 3);
      List<S2Point> emptyLoop = S2Loop.empty().vertices();
      loops.add(emptyLoop);
      // Empty loops should be ignored, leaving only the non-empty loops.
      S2Polygon polygon = checkValid(loops);
      assertEquals(loops.size() - 1, polygon.numLoops());
    }
  }

  @Test
  public void testIsValidFullLoop() {
    S2Loop fullLoop = S2Loop.full();
    List<S2Point> fullLoopPoints = Lists.newArrayList();
    for (int i = 0; i < fullLoop.numVertices(); ++i) {
      fullLoopPoints.add(fullLoop.vertex(i));
    }

    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      // This is only an error if there is at least one other loop.
      List<List<S2Point>> loops = getConcentricLoops(data.uniformInt(1, 6), 3);
      loops.add(fullLoopPoints);

      // We can't distinguish full/empty loops through the S2Shape API (only whether the shape as a
      // whole is the full or empty polygon). So when we have an extra full loop, we'll see it as an
      // empty loop via S2Shape.
      checkInvalid(
          loops,
          ImmutableSet.of(S2Error.Code.POLYGON_EXCESS_FULL_LOOP, S2Error.Code.POLYGON_EMPTY_LOOP));
    }
  }

  @Test
  public void testIsValidLoopsCrossing() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = getConcentricLoops(2, 4);
      // Both loops have the same number of vertices, and vertices at the same index position are
      // collinear with the center point, so we can create a crossing by simply exchanging two
      // vertices at the same index position.
      int n = loops.get(0).size();
      int i = data.random(n);
      S2Point temp = loops.get(0).get(i);
      loops.get(0).set(i, loops.get(1).get(i));
      loops.get(1).set(i, temp);
      if (data.oneIn(2)) {
        // By copying the two adjacent vertices from one loop to the other, we can ensure that the
        // crossings happen at vertices rather than edges.
        loops.get(0).set((i + 1) % n, loops.get(1).get((i + 1) % n));
        loops.get(0).set((i + n - 1) % n, loops.get(1).get((i + n - 1) % n));
      }
      checkInvalid(
          loops,
          ImmutableSet.of(
              S2Error.Code.OVERLAPPING_GEOMETRY,
              S2Error.Code.POLYGON_LOOPS_CROSS,
              S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS));
    }
  }

  @Test
  public void testIsValidDuplicateEdge() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = getConcentricLoops(2, 4);
      int n = loops.get(0).size();
      List<S2Point> loopA = loops.get(0);
      List<S2Point> loopB = loops.get(1);
      if (data.oneIn(2)) {
        // Create a shared edge (same direction in both loops).
        int i = data.random(n);
        loopA.set(i, loopB.get(i));
        loopA.set((i + 1) % n, loopB.get((i + 1) % n));
      } else {
        // Create a reversed edge (opposite direction in either loop) by cutting loop 0 into two
        // halves along one of its diagonals and replacing both loops with the result.
        int split = 2 + data.random(n - 3);
        loopB.clear();
        loopB.add(loopA.get(0));
        for (int s = split; s < n; ++s) {
          loopB.add(loopA.get(s));
        }
      }
      checkInvalid(
          loops,
          ImmutableSet.of(
              S2Error.Code.DUPLICATE_VERTICES,
              S2Error.Code.OVERLAPPING_GEOMETRY,
              S2Error.Code.POLYGON_LOOPS_SHARE_EDGE,
              S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS));
    }
  }

  @Test
  public void testIsValidInconsistentOrientations() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      List<List<S2Point>> loops = getConcentricLoops(data.uniformInt(2, 7), 3 /*min_vertices*/);
      S2Polygon polygon = new S2Polygon();
      uncheckedInitialize(() -> polygon.initOriented(loops(loops)));
      checkInvalid(polygon, S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
    }
  }

  @Test
  public void testIsValidLoopDepthNegative() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      S2Polygon poly =
          new S2Polygon(loops(getConcentricLoops(1 + data.random(4), 3 /*min_vertices*/)));
      int i = data.random(poly.numLoops());
      if (i == 0 || data.oneIn(3)) {
        poly.loop(i).setDepth(-1);
      } else {
        poly.loop(i).setDepth(poly.loop(i - 1).depth() + 2);
      }
      checkInvalid(poly, S2Error.Code.POLYGON_INVALID_LOOP_DEPTH);
    }
  }

  @Test
  public void testIsValidLoopNestingInvalid() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      int numLoops = data.uniformInt(2, 6);
      List<List<S2Point>> loops = getConcentricLoops(numLoops, /* minVertices= */ 3);

      // Randomly invert all the loops in order to generate cases where the outer loop encompasses
      // almost the entire sphere. This tests different code paths because bounding box checks are
      // not as useful.
      if (data.oneIn(2)) {
        for (List<S2Point> loop : loops) {
          Collections.reverse(loop);
        }
      }

      S2Polygon polygon = new S2Polygon(loops(loops));

      // Randomly invert one loop of the polygon to not match the others.
      int i = data.uniformInt(0, polygon.numLoops());
      polygon.loop(i).invert();

      checkInvalid(polygon, S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
    }
  }

  /**
   * Checks that the S2Loop / S2Polygon constructors, and {@code isValid()} methods don't crash when
   * they receive arbitrary invalid input. (We don't test large inputs; it is assumed that the
   * client enforces their own size limits before even attempting to construct geometric objects.)
   */
  @Test
  public void testIsValidFuzz() {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      int numLoops = 1 + data.random(10);
      List<List<S2Point>> loops = Lists.newArrayList();
      for (int i = 0; i < numLoops; ++i) {
        int numVertices = data.random(10);
        List<S2Point> loop = Lists.newArrayList();
        while (loop.size() < numVertices) {
          // Since the number of vertices is random, we automatically test empty loops, full loops,
          // and invalid vertex counts. Also, since most vertices are random, we automatically get
          // self-intersections and loop crossings. That leaves zero and NaN vertices, duplicate
          // vertices, and duplicate edges to be created explicitly.
          if (data.oneIn(10)) {
            // Zero vertex.
            loop.add(new S2Point(0, 0, 0));
          } else if (data.oneIn(10)) {
            // NaN vertex.
            loop.add(new S2Point(Double.NaN, Double.NaN, Double.NaN));
          } else if (data.oneIn(10) && !loop.isEmpty()) {
            // Duplicate vertex.
            loop.add(loop.get(data.random(loop.size())));
          } else if (data.oneIn(10) && loop.size() + 2 <= numVertices && !loops.isEmpty()) {
            // Try to copy an edge from a random loop.
            List<S2Point> other = loops.get(data.random(loops.size()));
            int n = other.size();
            if (n >= 2) {
              int k0 = data.random(n);
              int k1 = (k0 + 1) % n;
              if (data.oneIn(2)) {
                // Copy the reverse of this edge.
                int temp = k0;
                k0 = k1;
                k1 = temp;
              }
              loop.add(other.get(k0));
              loop.add(other.get(k1));
            }
          } else {
            // Add a random non-unit-length point.
            S2Point p = data.getRandomPoint();
            loop.add(p.mul(1e-30 * pow(1e60, data.nextDouble())));
          }
        }
        loops.add(loop);
      }
      // We could get any error message.
      checkInvalid(loops, "");
    }
  }

  private static S2Polygon simplify(String poly, double toleranceInDegrees) {
    return simplify(makePolygonOrDie(poly), toleranceInDegrees);
  }

  private static S2Polygon simplify(S2Polygon original, double toleranceInDegrees) {
    S2Polygon simplified = new S2Polygon();
    simplified.initToSimplified(
        original, new IdentitySnapFunction(S1Angle.degrees(toleranceInDegrees)));
    return simplified;
  }

  @Test
  public void testSimplifierNoSimplification() {
    S2Polygon original = makePolygonOrDie("0:0, 0:20, 20:20, 20:0");
    S2Polygon simplified = simplify(original, 1.0);
    assertEquals(4, simplified.numVertices());

    assertExactly(0, maximumDistanceInDegrees(simplified, original, 0));
    assertExactly(0, maximumDistanceInDegrees(original, simplified, 0));
  }

  // Here, 10:-2 will be removed and  0:0-20:0 will intersect two edges.
  // (The resulting polygon will in fact probably have more edges.)
  @Test
  public void testSimplifierSimplifiedLoopSelfIntersects() {
    S2Polygon original = makePolygonOrDie("0:0, 0:20, 10:-0.1, 20:20, 20:0, 10:-0.2");
    S2Polygon simplified = simplify(original, 0.22);

    // The simplified polygon has the same number of vertices but it should now consist of two loops
    // rather than one.
    assertEquals(2, simplified.numLoops());
    assertTrue(maximumDistanceInDegrees(simplified, original, 0) <= 0.22);
    assertTrue(maximumDistanceInDegrees(original, simplified, 0.22) <= 0.22);
  }

  @Test
  public void testSimplifierNoSimplificationManyLoops() {
    S2Polygon original =
        makePolygonOrDie( //
            "0:0,    0:1,   1:0;   0:20, 0:21, 1:20; " + "20:20, 20:21, 21:20; 20:0, 20:1, 21:0");
    S2Polygon simplified = simplify(original, 0.01);
    assertExactly(0, maximumDistanceInDegrees(simplified, original, 0));
    assertExactly(0, maximumDistanceInDegrees(original, simplified, 0));
  }

  @Test
  public void testSimplifierTinyLoopDisappears() {
    S2Polygon simplified = simplify("0:0, 0:1, 1:1, 1:0", 1.1);
    assertTrue(simplified.isEmpty());
  }

  @Test
  public void testSimplifierStraightLinesAreSimplified() {
    S2Polygon simplified =
        simplify("0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0," + "6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1", 0.01);
    assertEquals(4, simplified.numVertices());
  }

  @Test
  public void testSimplifierEdgeSplitInManyPieces() {
    // nearSquare's right four-point side will be simplified to a vertical line at lng=7.9, that
    // will cut the 9 teeth of the saw (the edge will therefore be broken into 19 pieces).
    String saw =
        "1:1, 1:8, 2:2, 2:8, 3:2, 3:8, 4:2, 4:8, 5:2, 5:8,"
            + "6:2, 6:8, 7:2, 7:8, 8:2, 8:8, 9:2, 9:8, 10:1";
    String nearSquare = "0:0, 0:7.9, 1:8.1, 10:8.1, 11:7.9, 11:0";
    S2Polygon original = makePolygonOrDie(saw + ";" + nearSquare);
    S2Polygon simplified = simplify(original, 0.21);

    assertTrue(simplified.isValid());
    assertTrue(maximumDistanceInDegrees(simplified, original, 0) <= 0.11);
    assertTrue(maximumDistanceInDegrees(original, simplified, 0) <= 0.11);
    // The resulting polygon's 9 little teeth are very small and disappear due to the
    // vertex_merge_radius of the polygon builder. There remains nine loops.
    assertEquals(9, simplified.numLoops());
  }

  @Test
  public void testSimplifierEdgesOverlap() {
    // Two loops, One edge of the second one ([0:1 - 0:2]) is part of an edge of the first one..
    S2Polygon simplified = simplify("0:0, 0:3, 1:0; 0:1, -1:1, 0:2", 0.01);
    S2Polygon truePoly = makePolygonOrDie("0:3, 1:0, 0:0, 0:1, -1:1, 0:2");
    assertTrue(simplified.boundaryApproxEquals(truePoly, S1Angle.radians(1e-15).radians()));
  }

  /**
   * Constructs a polygon (presumed to be invalid) from loops constructed from the given points.
   * Verifies that the S2Polygon constructor throws an AssertionError on the given loops, but does
   * not when wrapped in {@link S2TestSupport#uncheckedCreate(Supplier)}.
   *
   * <p>Checks that the resulting polygon is invalid, with an S2Error having a code that matches one
   * of the given 'codes'.
   */
  private static void checkInvalid(List<List<S2Point>> invalidLoops, Set<S2Error.Code> codes) {
    // First create a list of loops from the lists of points, and shuffle their order.
    List<S2Loop> invalidS2Loops = uncheckedLoops(invalidLoops);
    // Verify that S2Polygon construction from invalid loops throws an AssertionError.
    // TODO(user): After enabling assertions that polygons are valid, check that the
    // assertions work by enabling this "assertThrows".
    // assertThrows(AssertionError.class, () -> new S2Polygon(new ArrayList<>(invalidS2Loops)));

    // The same S2Polygon construction succeeds when wrapped in uncheckedCreate().
    S2Polygon invalidPolygon =
        uncheckedCreate(() -> new S2Polygon(new ArrayList<>(invalidS2Loops)));
    // Verify that findValidationError() finds the expected problem.
    checkInvalid(invalidPolygon, codes);
  }

  /**
   * Constructs a polygon (presumed to be invalid) from loops constructed from the given points.
   * Verifies that the S2Polygon constructor throws an AssertionError on the given loops, but does
   * not when wrapped in {@link S2TestSupport#uncheckedCreate(Supplier)}.
   *
   * <p>Checks that the resulting polygon is invalid, with an S2Error having a code that matches the
   * given 'code'.
   */
  private static void checkInvalid(List<List<S2Point>> invalidLoops, S2Error.Code code) {
    // First create a list of loops from the lists of points, and shuffle their order.
    List<S2Loop> invalidS2Loops = uncheckedLoops(invalidLoops);
    // Verify that S2Polygon construction from invalid loops throws an AssertionError.
    // TODO(user): After enabling assertions that polygons are valid, check that the
    // assertions work by enabling this "assertThrows".
    // assertThrows(AssertionError.class, () -> new S2Polygon(new ArrayList<>(invalidS2Loops)));

    // The same S2Polygon construction succeeds when wrapped in uncheckedCreate().
    S2Polygon invalidPolygon =
        uncheckedCreate(() -> new S2Polygon(new ArrayList<>(invalidS2Loops)));
    // Verify that findValidationError() finds the expected problem.
    checkInvalid(invalidPolygon, code);
  }

  /**
   * Checks that the given polygon is invalid, with an S2Error having one of the expected 'codes'.
   */
  private static void checkInvalid(S2Polygon polygon, Set<S2Error.Code> codes) {
    ImmutableSet<String> strCodes = codes.stream().map(Enum::toString).collect(toImmutableSet());
    S2Error error = new S2Error();
    boolean hasError = polygon.findValidationError(error);
    assertTrue(
        "No error found, but expected validation error with one of the codes: '"
            + strCodes
            + "'\nfor polygon "
            + S2TextFormat.toString(polygon),
        hasError);
    assertTrue(
        "Actual Error: '"
            + error.text()
            + "' with code "
            + error.code()
            + "\n but expected one of the codes: "
            + strCodes,
        codes.contains(error.code()));
  }

  /** Checks that the given polygon is invalid, with an S2Error having the expected 'code'. */
  private static void checkInvalid(S2Polygon polygon, S2Error.Code code) {
    S2Error error = new S2Error();
    boolean hasError = polygon.findValidationError(error);
    assertTrue(
        "No error found, but expected validation error with the code: '"
            + code
            + "'\nfor polygon "
            + S2TextFormat.toString(polygon),
        hasError);
    assertTrue(
        "Actual Error: '"
            + error.text()
            + "' with code "
            + error.code()
            + "\n but expected the code: "
            + code,
        code == error.code());
  }

  /**
   * Constructs a polygon (presumed to be invalid) from loops constructed from the given points.
   * Verifies that the S2Polygon constructor throws an AssertionError on the given loops, but does
   * not when wrapped in {@link S2TestSupport#uncheckedCreate(Supplier)}.
   *
   * <p>Checks that the resulting polygon is invalid, with an S2Error that matches the given
   * snippet, i.e. either error.text() contains the snippet, or error.code().name() equals the
   * snippet.
   */
  private static void checkInvalid(List<List<S2Point>> invalidLoops, String snippet) {
    // First create a list of loops from the lists of points, and shuffle their order.
    List<S2Loop> invalidS2Loops = uncheckedLoops(invalidLoops);
    // Verify that S2Polygon construction from invalid loops throws an AssertionError.
    // TODO(user): After enabling assertions that polygons are valid, check that the
    // assertions work by enabling this "assertThrows".
    // assertThrows(AssertionError.class, () -> new S2Polygon(new ArrayList<>(invalidS2Loops)));

    // The same S2Polygon construction succeeds when wrapped in uncheckedCreate().
    S2Polygon invalidPolygon =
        uncheckedCreate(() -> new S2Polygon(new ArrayList<>(invalidS2Loops)));
    // Verify that findValidationError() finds the expected problem.
    checkInvalid(invalidPolygon, snippet);
  }

  /**
   * Checks that the given polygon is invalid, with an S2Error that has text containing the given
   * snippet, or S2Error.Code.name() equal to the given snippet, like "DUPLICATE_VERTICES".
   */
  private static void checkInvalid(S2Polygon polygon, String snippet) {
    S2Error error = new S2Error();
    boolean hasError = polygon.findValidationError(error);
    assertTrue(
        "No error found, but expected validation error with substring or code: '"
            + snippet
            + "'\nfor polygon "
            + S2TextFormat.toString(polygon),
        hasError);
    assertTrue(
        "Actual Error: '" + error.text() + "' with code " + error.code()
            + "\n but expected substring: " + snippet,
        error.text().contains(snippet) || error.code().name().equals(snippet));
  }

  /**
   * Constructs S2Loops from the given lists of S2Points, shuffles those loops, and then constructs
   * a polygon (presumed to be valid) from those S2Loops. Verifies that the S2Polygon is valid, and
   * returns it.
   */
  private static S2Polygon checkValid(List<List<S2Point>> validLoops) {
    List<S2Loop> loops = Lists.newArrayList();
    for (int i = 0; i < validLoops.size(); ++i) {
      loops.add(new S2Loop(validLoops.get(i)));
    }
    Collections.shuffle(loops);

    S2Polygon polygon = new S2Polygon(loops);
    S2Error error = new S2Error();
    boolean hasError = polygon.findValidationError(error);
    assertFalse(error.text(), hasError);
    return polygon;
  }

  /**
   * Constructs a List of possibly invalid S2Loops from the given lists of points. The vertices are
   * allowed to not be unit length. Shuffles the order of the resulting loops.
   */
  private static List<S2Loop> uncheckedLoops(List<List<S2Point>> vertices) {
    List<S2Loop> loops = Lists.newArrayList();
    for (int i = 0; i < vertices.size(); ++i) {
      loops.add(uncheckedLoop(vertices.get(i)));
    }
    Collections.shuffle(loops);
    return loops;
  }

  private static S2Loop uncheckedLoop(List<S2Point> vertices) {
    return uncheckedCreate(() -> new S2Loop(vertices));
  }

  /** Constructs S2Loops and shuffles their order. All vertices must at least be unit length. */
  private static List<S2Loop> loops(List<List<S2Point>> vertices) {
    List<S2Loop> loops = Lists.newArrayList();
    for (int i = 0; i < vertices.size(); ++i) {
      loops.add(new S2Loop(vertices.get(i)));
    }
    Collections.shuffle(loops);
    return loops;
  }

  /** Constructs S2Loops and checks they are valid. */
  private static List<S2Loop> makeLoops(List<List<S2Point>> loopVertices) {
    List<S2Loop> result = new ArrayList<>();
    for (List<S2Point> vertices : loopVertices) {
      S2Loop loop = new S2Loop(vertices);
      result.add(loop);
      S2Error error = new S2Error();
      boolean invalid = loop.findValidationError(error);
      assertFalse("Loop " + (result.size() - 1) + ": " + error, invalid);
    }
    return result;
  }

  @Test
  public void testEncodeDecode() throws IOException {
    // Empty polygon.
    encodeDecode(empty);

    // Full polygon
    encodeDecode(full);

    for (TestCase testCase : testCases) {
      encodeDecode(makePolygon(testCase.a, S2CellId.MAX_LEVEL));
      encodeDecode(makePolygon(testCase.b, S2CellId.MAX_LEVEL));
      encodeDecode(makePolygon(testCase.aAndB, S2CellId.MAX_LEVEL));
      encodeDecode(makePolygon(testCase.aOrB, S2CellId.MAX_LEVEL));
      encodeDecode(makePolygon(testCase.aMinusB, S2CellId.MAX_LEVEL));
    }
  }

  @Test
  public void testEncodingSize_emptyPolygon() throws IOException {
    byte[] encoded = encode(new S2Polygon());
    // 1 byte for version, 1 for the level, 1 for the length.
    assertEquals(3, encoded.length);
  }

  @Test
  public void testEncodingSizeSnappedPolygon() throws IOException {
    S2Polygon snapped = makePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0", S2CellId.MAX_LEVEL);
    byte[] encoded = encode(snapped);

    // 2 loops, one with 3 vertices, one with 4.
    // Polygon:
    //   1 byte for version
    //   1 byte for level
    //   1 byte for numLoops
    // Loops:
    //   5 bytes overhead
    //   8 bytes per vertex
    assertEquals(1 + 1 + 1 + 2 * 5 + 7 * 8, encoded.length);
  }

  private void polylineIntersectionSharedEdgeTest(S2Polygon p, int startVertex, int direction) {
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
   * This tests polyline-polyline intersections with shared edges.
   */
  @Test
  public void testPolylineIntersectionSharedEdges() {
    for (int v = 0; v < 3; ++v) {
      polylineIntersectionSharedEdgeTest(cross1, v, 1);
      polylineIntersectionSharedEdgeTest(cross1, v + 1, -1);
      polylineIntersectionSharedEdgeTest(cross1SideHole, v, 1);
      polylineIntersectionSharedEdgeTest(cross1SideHole, v + 1, -1);
    }
  }

  /**
   * This tests polyline-polyline intersections. It covers the same edge cases as {@code
   * testOperations} and also adds some extra tests for shared edges.
   */
  @Test
  public void testPolylineIntersection() {
    // This duplicates some of the tests in testOperations by converting the outline of polygon A to
    // a polyline then intersecting it with the polygon B. It then converts B to a polyline and
    // intersects it with A. It then feeds all of the results into a polygon builder and tests that
    // the output is equal to doing an intersection between A and B.

    for (int testNumber = 0; testNumber < testCases.length; testNumber++) {
      TestCase test = testCases[testNumber];
      logger.fine("Polyline intersection test case " + testNumber);
      S2Polygon a = makeVerbatimPolygon(test.a);
      S2Polygon b = makeVerbatimPolygon(test.b);
      S2Polygon expectedAAndB = makeVerbatimPolygon(test.aAndB);

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
          S2Polyline polyline = new S2Polyline(points);
          try {
            polylines.addAll(tmp2.intersectWithPolyline(polyline));
          } catch (AssertionError e) {
            System.err.println("AssertionError " + e.getMessage());
            System.err.println(" during S2Polygon.intersectWithPolyline(S2Polyline) with");
            System.err.println(" S2Polygon = " + S2TextFormat.toString(tmp2.index()));
            System.err.println(" S2Polyline = " + S2TextFormat.toString(polyline));
            throw e;
          }
        }
      }

      S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.DIRECTED_XOR);
      for (S2Polyline polyLine : polylines) {
        for (int j = 0; j < polyLine.numVertices() - 1; j++) {
          builder.addEdge(polyLine.vertex(j), polyLine.vertex(j + 1));
          logger.fine(" ... Adding edge: " + polyLine.vertex(j) + " - " + polyLine.vertex(j + 1));
        }
      }

      S2Polygon aAndB = new S2Polygon();
      assertTrue(builder.assemblePolygon(aAndB, null));
      checkEqual(aAndB, expectedAAndB, OPERATIONS_MAX_ERROR);
    }
  }

  private static void checkEqual(S2Polygon a, S2Polygon b) {
    checkEqual(a, b, 0);
  }

  private static void checkEqual(S2Polygon a, S2Polygon b, final double maxError) {
    checkEqual("", a, b, maxError);
  }

  private static void checkEqual(String message, S2Polygon a, S2Polygon b, final double maxError) {
    if (a.boundaryApproxEquals(b, maxError)) {
      return;
    }
    message +=
        "Expected polygons to be equal: \nA:"
            + S2TextFormat.toString(a)
            + "\nB:"
            + S2TextFormat.toString(b);

    message += "\nBoundaries are not approximately equal. Trying reconstruction...\n";

    S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.DIRECTED_XOR);
    builder.addPolygon(a);
    S2Polygon a2 = new S2Polygon();
    assertTrue(
        message + "\nFailed to assemble the loops of A into a new polygon.",
        builder.assemblePolygon(a2, null));
    builder.addPolygon(b);
    S2Polygon b2 = new S2Polygon();
    assertTrue(
        message + "\nFailed to assemble the loops of B into a new polygon.",
        builder.assemblePolygon(b2, null));

    message +=
        "\nReconstructed polygons:\nA2:"
            + S2TextFormat.toString(a2)
            + "\nB2:"
            + S2TextFormat.toString(b2);

    assertTrue(
        message + "\nBoundaries of reconstructed polygons are not approximately equal.",
        a2.boundaryApproxEquals(b2, maxError));
  }

  private static void checkComplementary(S2Polygon a, S2Polygon b) {
    S2Polygon b1 = new S2Polygon();
    b1.initToComplement(b);
    checkEqual(a, b1);
  }

  @Test
  public void testApproxContains() {
    // Get a random S2Cell as a polygon.
    S2CellId id = S2CellId.fromLatLng(S2LatLng.fromE6(69852241, 6751108));
    S2Cell cell = new S2Cell(id.parent(10));
    S2Polygon cellAsPolygon = new S2Polygon(cell);

    // We want to roughly bisect the polygon, so we make a rectangle that is the top half of the
    // current polygon's bounding rectangle.
    S2LatLngRect bounds = cellAsPolygon.getRectBound();
    S2LatLngRect upperHalf =
        new S2LatLngRect(new R1Interval(bounds.lat().getCenter(), bounds.lat().hi()), bounds.lng());

    // Turn the S2LatLngRect into an S2Polygon
    List<S2Point> points = Lists.newArrayList();
    for (int i = 0; i < 4; i++) {
      points.add(upperHalf.getVertex(i).toPoint());
    }
    List<S2Loop> loops = Lists.newArrayList();
    loops.add(new S2Loop(points));
    S2Polygon upperHalfPolygon = new S2Polygon(loops);

    // Get the intersection. There is no guarantee that the intersection will be contained by A or
    // B.
    S2Polygon intersection = new S2Polygon();
    intersection.initToIntersection(cellAsPolygon, upperHalfPolygon);
    assertFalse(cellAsPolygon.contains(intersection));
    assertTrue(
        cellAsPolygon.approxContains(intersection, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
  }

  public void tryUnion(S2Polygon a, S2Polygon b) {
    S2Polygon union1 = new S2Polygon();
    union1.initToUnion(a, b);

    S2Polygon union2 = S2Polygon.union(Arrays.asList(new S2Polygon(a), new S2Polygon(b)));

    checkEqual(union1, union2);
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
    assertFalse("Expect disjoint polygons to not intersect", a.intersects(b));
    assertFalse("Expect disjoint polygons to not intersect", b.intersects(a));
    assertEquals("Expect disjoint polygon to not contain the other", b.isEmpty(), a.contains(b));
    assertEquals("Expect disjoint polygon to not contain the other", a.isEmpty(), b.contains(a));

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
    // Check DeMorgan's Law and that subtraction is the same as intersection with the complement.
    // This function is called multiple times in order to test the various combinations of
    // complements.
    S2Polygon aOrB = new S2Polygon();
    S2Polygon aAndB1 = new S2Polygon();
    S2Polygon aMinusB = new S2Polygon();
    aAndB1.initToIntersection(a, b1);
    aOrB.initToUnion(a1, b);
    aMinusB.initToDifference(a, b);
    checkComplementary(aOrB, aAndB1);
    checkEqual(aMinusB, aAndB1);
  }

  /** Test identities that should hold for any pair of polygons A, B and their complements. */
  private static void checkComplements(S2Polygon a, S2Polygon b) {
    S2Polygon a1 = new S2Polygon();
    a1.initToComplement(a);
    S2Polygon b1 = new S2Polygon();
    b1.initToComplement(b);
    checkOneComplementPair(a, a1, b, b1);
    checkOneComplementPair(a1, a, b, b1);
    checkOneComplementPair(a, a1, b1, b);
    checkOneComplementPair(a1, a, b1, b);

    // There is a lot of redundancy if we do this test for each complementary pair, so we just do it
    // once instead.
    S2Polygon aXorB1 = new S2Polygon();
    S2Polygon a1XorB = new S2Polygon();
    aXorB1.initToSymmetricDifference(a, b1);
    a1XorB.initToSymmetricDifference(a1, b);
    checkEqual(aXorB1, a1XorB);
  }

  private static void checkUnion(S2Polygon a, S2Polygon b) {
    S2Polygon aCopy = new S2Polygon(a);
    S2Polygon bCopy = new S2Polygon(b);
    // Check the copy constructor.
    checkEqual(a, aCopy);
    checkEqual(b, bCopy);

    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);

    List<S2Polygon> polygons = Lists.newArrayList();
    polygons.add(aCopy);
    polygons.add(bCopy);
    S2Polygon cUnion = S2Polygon.union(polygons);
    checkEqual(c, cUnion);
  }

  private static void checkRelationImpl(
      S2Polygon a, S2Polygon b, boolean contains, boolean contained, boolean intersects) {
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
    if (intersects && !(contains || contained)) {
      checkOverlappingPair(a, b);
    }
    checkUnion(a, b);
    checkComplements(a, b);
  }

  private static void checkRelation(
      S2Polygon a, S2Polygon b, boolean contains, boolean contained, boolean intersects) {
    try {
      checkRelationImpl(a, b, contains, contained, intersects);
    } catch (AssertionError e) {
      System.err.println("S2Polygon a: " + S2TextFormat.toString(a.index));
      System.err.println("S2Polygon b: " + S2TextFormat.toString(b.index));
      throw e;
    }
  }

  @Test
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

  @Test
  public void testUnionSloppySuccess() {
    S2Polygon union = S2Polygon.unionSloppy(Arrays.asList(adj0, adj1), S1Angle.degrees(0.1));

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

  @Test
  public void testUnionSloppyFailure() {
    // The polygons are sufficiently far apart that this angle will not
    // bring them together:
    S2Polygon union = S2Polygon.unionSloppy(Arrays.asList(adj0, unAdj), S1Angle.degrees(0.1));

    assertEquals(2, union.numLoops());
  }

  /**
   * Test for b/373414882. The two nearly-equal points shown in screen/6VXEja3nCcnvzVj generate
   * intersecting edges which, due to floating point errors, cause S2PolygonBuilder.spliceEdges to
   * loop infinitely if the default mergeRadius is too small.
   */
  @Test
  public void testUnionErrorTolerance() {
    var aPoints =
        new S2Point[] {
          new S2Point(-0.4233762535751846, -0.6818147440520357, 0.5965577949385107),
          new S2Point(-0.4233760956593306, -0.6818147800268495, 0.5965578658950028),
          new S2Point(-0.4233761112005162, -0.6818147082926758, 0.5965579368514861),
          new S2Point(-0.42337613281264846, -0.6818146817526661, 0.5965579518463546),
          new S2Point(-0.4233762907285011, -0.6818146457778499, 0.5965578808898605),
          new S2Point(-0.42337627518731696, -0.6818147175120258, 0.5965578099333793),
        };
    var bPoints =
        new S2Point[] {
          new S2Point(-0.42337629072850064, -0.6818146457778491, 0.5965578808898617),
          new S2Point(-0.4233761884587999, -0.6818147713665217, 0.5965578099333851),
          new S2Point(-0.4233761172714648, -0.6818147534868428, 0.5965578808898722),
          new S2Point(-0.4233761149633853, -0.6818147497698579, 0.5965578867761089),
          new S2Point(-0.423376217233087, -0.6818146241811852, 0.5965579577325839),
          new S2Point(-0.4233762884204217, -0.6818146420608652, 0.5965578867760968),
        };
    S2Polygon a = new S2Polygon(new S2Loop(aPoints));
    S2Polygon b = new S2Polygon(new S2Loop(bPoints));

    S2Polygon union = new S2Polygon();
    union.initToUnionSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);

    assertEquals(12, union.getNumVertices());
  }

  @Test
  public void testCompareTo() {
    // Polygons with same logical loop (but loop is reordered).
    S2Polygon p5 = makePolygon(TRIANGLE);
    S2Polygon p6 = makePolygon(TRIANGLE_ROT);
    assertEquals(0, p5.compareTo(p6));

    // Polygons with a differing number of loops
    S2Polygon p7 = makePolygon(RECTANGLE1 + RECTANGLE2);
    S2Polygon p8 = makePolygon(TRIANGLE);
    assertTrue(p8.compareTo(p7) < 0);
    assertTrue(p7.compareTo(p8) > 0);

    // Polygons with a differing number of loops (one a subset of the other)
    S2Polygon p9 = makePolygon(RECTANGLE1 + RECTANGLE2 + TRIANGLE);
    S2Polygon p10 = makePolygon(RECTANGLE1 + RECTANGLE2);
    assertTrue(p9.compareTo(p10) > 0);
    assertTrue(p10.compareTo(p9) < 0);
  }

  @Test
  public void testGetDistance() {
    // Error margin since we're doing numerical computations.
    double epsilon = 1e-15;

    // A rectangle with (lat,lng) vertices (3,1), (3,-1), (-3,-1) and (-3,1).
    String inner = "3:1, 3:-1, -3:-1, -3:1;";
    // A larger rectangle with (lat,lng) vertices (4,2), (4,-2), (-4,-2) and (-4,s).
    String outer = "4:2, 4:-2, -4:-2, -4:2;";

    S2Polygon rect = makePolygon(inner);
    S2Polygon shell = makePolygon(inner + outer);

    // All of the vertices of a polygon should be distance 0.
    for (int i = 0; i < shell.numLoops(); i++) {
      for (int j = 0; j < shell.loop(i).numVertices(); j++) {
        assertEquals(0d, shell.getDistance(shell.loop(i).vertex(j)).radians(), epsilon);
      }
    }

    // A non-vertex point on an edge should be distance 0.
    assertEquals(
        0d,
        rect.getDistance(rect.loop(0).vertex(0).add(rect.loop(0).vertex(1)).normalize()).radians(),
        epsilon);

    S2Point origin = S2LatLng.fromDegrees(0, 0).toPoint();
    // rect contains the origin.
    assertEquals(0d, rect.getDistance(origin).radians(), epsilon);

    // shell does NOT contain the origin, since it has a hole. The shortest distance is to (1,0) or
    // (-1,0), and should be 1 degree.
    assertEquals(1d, shell.getDistance(origin).degrees(), epsilon);
  }

  @Test
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

  @Test
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
    projected = polygon.project(point);
    assertTrue(makePoint("5:-2").aequal(projected, epsilon));

    // The point is inside the hole in the polygon. Note the expected value is based on a plane, so
    // it's not that accurate; thus, tolerance is reduced to 1e-6.
    point = makePoint("-0.49:-0.49");
    projected = polygon.project(point);
    assertTrue(makePoint("-0.5:-0.5").aequal(projected, 1e-6));

    point = makePoint("0:-3");
    projected = polygon.project(point);
    assertTrue(makePoint("0:-2").aequal(projected, epsilon));

    // Project for full or empty polygons, with no edges, are expected to return the given point.
    assertTrue(empty.project(point).equalsPoint(point));
    assertTrue(full.project(point).equalsPoint(point));
  }

  @Test
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

  @Test
  public void testFastInit() {
    S2LatLngRect bound;
    IdentityHashMap<S2Loop, List<S2Loop>> nestedLoops = new IdentityHashMap<>();
    nestedLoops.put(S2Polygon.ROOT, Lists.<S2Loop>newArrayList());

    List<S2Point> vertices = Lists.newArrayList();

    bound = parseVertices("-2:-2, -3:6, 6:-3", vertices);
    S2Loop loop1 = S2Loop.newLoopWithTrustedDetails(vertices, false, bound);
    nestedLoops.get(S2Polygon.ROOT).add(loop1);
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

    S2Error error = new S2Error();
    boolean hasError = polygon.findValidationError(error);
    assertFalse(error.text(), hasError);
    assertEquals(0, polygon.loop(0).depth());
    assertEquals(1, polygon.loop(1).depth());
    assertEquals(2, polygon.loop(2).depth());
    assertDoubleNear(0.003821967440517272, polygon.getArea());
  }

  @Test
  public void testInitToSnappedWithSnapLevel() {
    S2Polygon polygon = makePolygonOrDie("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");

    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      S2Polygon snappedPolygon = new S2Polygon();
      snappedPolygon.initToSnapped(polygon, level);
      assertTrue(snappedPolygon.isValid());
      double cellAngle = MAX_DIAG.getValue(level);
      S1Angle mergeRadius = S1Angle.min(
          S2BuilderSnapFunctions.maxSnapRadius(),
          S1Angle.radians(cellAngle));
      assertTrue(
          "snapped polygon should approx contain original polygon for"
              + "\nsnap level = "
              + level
              + ", mergeRadius = "
              + mergeRadius
              + "\noriginal polygon: "
              + polygon
              + "\nsnapped polygon: "
              + snappedPolygon,
          snappedPolygon.approxContains(polygon, mergeRadius));
    }
  }

  private static class S2PolygonSimplifierTest {
    final S2Polygon original;
    final S2Polygon simplified;

    public S2PolygonSimplifierTest(S2Polygon original, double toleranceInDegrees) {
      this.original = original;
      this.simplified = new S2Polygon();
      this.simplified.initToSimplified(
          original, new IdentitySnapFunction(S1Angle.degrees(toleranceInDegrees)));
    }

    public S2PolygonSimplifierTest(String poly, double toleranceInDegrees) {
      this(makeVerbatimPolygonOrDie(poly), toleranceInDegrees);
    }
  }

  @Test
  public void testInitToSimplifiedNoSimplification() {
    S2PolygonSimplifierTest s = new S2PolygonSimplifierTest("0:0, 0:20, 20:20, 20:0", 1.0);
    assertEquals(4, s.simplified.getNumVertices());

    assertExactly(0.0, maximumDistanceInDegrees(s.simplified, s.original, 0));
    assertExactly(0.0, maximumDistanceInDegrees(s.original, s.simplified, 0));
  }

  // Here, 10:-2 will be removed and  0:0-20:0 will intersect two edges. (The resulting polygon will
  // in fact probably have more edges.)
  @Test
  public void testInitToSimplifiedSimplifiedLoopSelfIntersects() {
    S2PolygonSimplifierTest s =
        new S2PolygonSimplifierTest("0:0, 0:20, 10:-0.1, 20:20, 20:0, 10:-0.2", 0.22);

    // The simplified polygon has the same number of vertices but it should now consists of two
    // loops rather than one.
    assertEquals(2, s.simplified.numLoops());
    assertGreaterOrEqual(0.22, maximumDistanceInDegrees(s.simplified, s.original, 0));
    assertGreaterOrEqual(0.22, maximumDistanceInDegrees(s.original, s.simplified, 0.22));
  }

  @Test
  public void testInitToSimplifiedNoSimplificationManyLoops() {
    S2PolygonSimplifierTest s =
        new S2PolygonSimplifierTest(
            "0:0,    0:1,   1:0;   0:20, 0:21, 1:20; 20:20, 20:21, 21:20; 20:0, 20:1, 21:0", 0.01);
    assertExactly(0.0, maximumDistanceInDegrees(s.simplified, s.original, 0));
    assertExactly(0.0, maximumDistanceInDegrees(s.original, s.simplified, 0));
  }

  @Test
  public void testInitToSimplifiedTinyLoopDisappears() {
    S2PolygonSimplifierTest s = new S2PolygonSimplifierTest("0:0, 0:1, 1:1, 1:0", 1.1);
    assertTrue(s.simplified.isEmpty());
  }

  // Check that initFromBuilder correctly chooses between an empty or a full polygon when
  // simplification produces a polygon with no loops.
  @Test
  public void initToSimplifiedFullPolygonPredicateTest() {
    // A narrow rectangle polygon that includes the north pole and extends to the equator. This has
    // very large bounds but a small area. Simplifying with 10 degree tolerance will collapse
    // the rectangle to nothing.
    S2Polygon bigBoundsSmallArea =
        makePolygonVerbatimOrDie("0:1, 85:1, 85:179, 85:-179, 85:-1, 0:-1");
    assertGreaterThan(bigBoundsSmallArea.getRectBound().area(), 2 * PI);
    assertLessThan(bigBoundsSmallArea.getArea(), 1.0);

    S2PolygonSimplifierTest s1 = new S2PolygonSimplifierTest(bigBoundsSmallArea, 10);
    assertTrue(s1.simplified.isEmpty());

    // A polygon that has both bounds and area > 2* pi because it is a small clockwise loop, i.e.
    // hole in the sphere. Simplifying with 10 degree tolerance will collapse it to nothing.
    S2Polygon hole = makePolygonVerbatimOrDie("10:10, 10:8, 12:8, 12:10");

    assertGreaterThan(hole.getRectBound().area(), 2 * PI);
    assertGreaterThan(hole.getArea(), 2 * PI);

    S2PolygonSimplifierTest s2 = new S2PolygonSimplifierTest(hole, 10);
    assertTrue(s2.simplified.isFull());
  }

  @Test
  public void testInitToSimplifiedStraightLinesAreSimplified() {
    S2PolygonSimplifierTest s =
        new S2PolygonSimplifierTest(
            "0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1", 0.01);
    assertEquals(4, s.simplified.getNumVertices());
  }

  @Test
  public void testInitToSimplifiedEdgeSplitInManyPieces() {
    // nearSquare's right four-point side will be simplified to a vertical line at lng=7.9, that
    // will cut the 9 teeth of the saw (the edge will therefore be broken into 19 pieces).
    String saw =
        "1:1, 1:8, 2:2, 2:8, 3:2, 3:8, 4:2, 4:8, 5:2, 5:8,"
            + "6:2, 6:8, 7:2, 7:8, 8:2, 8:8, 9:2, 9:8, 10:1";
    String nearSquare = "0:0, 0:7.9, 1:8.1, 10:8.1, 11:7.9, 11:0";
    S2PolygonSimplifierTest s = new S2PolygonSimplifierTest(saw + ";" + nearSquare, 0.21);

    assertTrue(s.simplified.isValid());
    assertGreaterOrEqual(0.11, maximumDistanceInDegrees(s.simplified, s.original, 0));
    assertGreaterOrEqual(0.11, maximumDistanceInDegrees(s.original, s.simplified, 0));
    // The resulting polygon's 9 little teeth are very small and disappear due to the
    // vertexMergeRadius of the builder. Nine loops remain.
    assertEquals(9, s.simplified.numLoops());
  }

  @Test
  public void testInitToSimplifiedEdgesOverlap() {
    // Two loops, One edge of the second one ([0:1 - 0:2]) is part of an edge of the first one.
    S2PolygonSimplifierTest s = new S2PolygonSimplifierTest("0:0, 0:3, 1:0; 0:1, -1:1, 0:2", 0.01);
    S2Polygon truePoly = makePolygonOrDie("0:3, 1:0, 0:0, 0:1, -1:1, 0:2");
    assertTrue(s.simplified.boundaryApproxEquals(truePoly, 1e-15));
  }

  /**
   * Returns the maximum distance from any vertex of polyA to polyB, that is, the directed Haussdorf
   * distance of the set of vertices of polyA to the boundary of polyB. Doesn't consider loops from
   * polyA that have diameter less than minDiameterInDegrees.
   */
  private static double maximumDistanceInDegrees(
      S2Polygon polyA, S2Polygon polyB, double minDiameterInDegrees) {
    double minDistance = 360;
    boolean hasBigLoops = false;
    for (int l = 0; l < polyA.numLoops(); ++l) {
      S2Loop aLoop = polyA.loop(l);
      if (loopDiameter(aLoop).degrees() <= minDiameterInDegrees) {
        continue;
      }
      hasBigLoops = true;
      for (int v = 0; v < aLoop.numVertices(); ++v) {
        double distance = polyB.getDistance(aLoop.vertex(v)).degrees();
        if (distance < minDistance) {
          minDistance = distance;
        }
      }
    }
    if (hasBigLoops) {
      return minDistance;
    } else {
      return 0.; // As if the first polygon were empty.
    }
  }

  /** Returns the diameter of a loop (maximum distance between any two points in the loop). */
  private static S1Angle loopDiameter(S2Loop loop) {
    S1Angle diameter = new S1Angle();
    for (int i = 0; i < loop.numVertices(); ++i) {
      S2Point testPoint = loop.vertex(i);
      for (int j = i + 1; j < loop.numVertices(); ++j) {
        diameter =
            S1Angle.max(
                diameter, S2EdgeUtil.getDistance(testPoint, loop.vertex(j), loop.vertex(j + 1)));
      }
    }
    return diameter;
  }

  @Test
  public void testGetBestSnapLevel() {
    // Unsnapped polygon.
    assertEquals(-1, makePolygon("10:10, 30:20, 20:30").getBestSnapLevel());
    // Polygons snapped at particular levels.
    assertEquals(20, makePolygon("10:10, 30:20, 20:30", 20).getBestSnapLevel());
    assertEquals(
        S2CellId.MAX_LEVEL,
        makePolygon("10:10, 30:20, 20:30", S2CellId.MAX_LEVEL).getBestSnapLevel());
    // Mixed snap levels for different loops.
    List<S2Loop> loops = Lists.newArrayList();
    loops.add(makeLoop("10:10, 20:30, 30:20", 10));
    loops.add(makeLoop("60:60, 70:80, 80:70", 12));
    assertEquals(10, new S2Polygon(loops).getBestSnapLevel());
  }

  @Test
  public void testInitToSimplifiedInCellPointsOnCellBoundaryKept() {
    S2CellId cellId = S2CellId.fromToken("89c25c");
    S2Cell cell = new S2Cell(cellId);
    S2Polygon polygon = makeCellPolygon(cell, "0.1:0, 0.2:0, 0.2:0.5");
    double angle = polygon.loop(0).vertex(0).angle(polygon.loop(0).vertex(1));
    S1Angle tolerance = S1Angle.radians(angle * 1.1);
    S2Polygon simplified = new S2Polygon();
    simplified.initToSimplified(polygon, tolerance, false);
    assertTrue(simplified.isEmpty());

    S2Polygon simplifiedInCell = new S2Polygon();
    simplifiedInCell.initToSimplifiedInCell(polygon, cell, tolerance);

    assertTrue(simplifiedInCell.boundaryEquals(polygon));
    assertEquals(3, simplifiedInCell.getNumVertices());
    assertEquals(-1, simplifiedInCell.getSnapLevel());
  }

  @Test
  public void testInitToSimplifiedInCellPointsInsideCellSimplified() {
    S2CellId cellId = S2CellId.fromToken("89c25c");
    S2Cell cell = new S2Cell(cellId);
    S2Polygon loop = makeCellPolygon(cell, "0.3:0, 0.4:0, 0.4:0.5, 0.4:0.8, 0.2:0.8");
    S1Angle tolerance = S1Angle.radians(loop.loop(0).vertex(0).angle(loop.loop(0).vertex(1)) * 1.1);
    S2Polygon simplifiedLoop = new S2Polygon();
    simplifiedLoop.initToSimplifiedInCell(loop, cell, tolerance);
    assertTrue(loop.boundaryNear(simplifiedLoop, 1e-15));
    assertEquals(4, simplifiedLoop.getNumVertices());
    assertEquals(-1, simplifiedLoop.getSnapLevel());
  }

  @Test
  public void testInitToSimplifiedInCellCellCornerKept() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("00001"));
    S2Polygon input = makeCellPolygon(cell, "1:0, 1:0.05, 0.99:0");
    S1Angle tolerance = new S1Angle(cell.getVertex(0), cell.getVertex(1)).mul(0.02);
    S2Polygon simplified = new S2Polygon();
    simplified.initToSimplifiedInCell(input, cell, tolerance);
    assertTrue(simplified.boundaryNear(input, 1e-15));
  }

  @Test
  public void testInitToSimplifiedInCellNarrowStripRemoved() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("00001"));
    S2Polygon input = makeCellPolygon(cell, "0.9:0, 0.91:0, 0.91:1, 0.9:1");
    S1Angle tolerance = new S1Angle(cell.getVertex(0), cell.getVertex(1)).mul(0.02);
    S2Polygon simplified = new S2Polygon();
    simplified.initToSimplifiedInCell(input, cell, tolerance);
    assertTrue(simplified.isEmpty());
  }

  @Test
  public void testInitToSimplifiedInCellNarrowGapRemoved() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("00001"));
    S2Polygon input =
        makeCellPolygon(cell, "0.7:0, 0.75:0, 0.75:1, 0.7:1", "0.76:0, 0.8:0, 0.8:1, 0.76:1");
    S2Polygon expected = makeCellPolygon(cell, "0.7:0, 0.8:0, 0.8:1, 0.7:1");
    S1Angle tolerance = new S1Angle(cell.getVertex(0), cell.getVertex(1)).mul(0.02);
    S2Polygon simplified = new S2Polygon();
    simplified.initToSimplifiedInCell(input, cell, tolerance);
    assertTrue(simplified.boundaryNear(expected, 1e-15));
  }

  @Test
  public void testInitToSimplifiedInCellCloselySpacedEdgeVerticesKept() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("00001"));
    S2Polygon input = makeCellPolygon(cell, "0:0.303, 0:0.302, 0:0.301, 0:0.3, 0.1:0.3, 0.1:0.4");
    S1Angle tolerance = new S1Angle(cell.getVertex(0), cell.getVertex(1)).mul(0.02);
    S2Polygon simplified = new S2Polygon();
    simplified.initToSimplifiedInCell(input, cell, tolerance);
    assertTrue(simplified.boundaryApproxEquals(input, 1e-15));
  }

  @Test
  public void testInitToSimplifiedInCellPolylineAssemblyBug() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("5701"));
    S2Polygon polygon =
        makePolygonOrDie(
            "55.8699252:-163.9412145, " // South-west corner of 5701
                + "54.7672352:-166.7579678, " // North-east corner of 5701
                /* Offending part: a tiny triangle near south-east corner */
                + "54.7109214:-164.6376338, " // forced vertex, on edge 4
                + "54.7140193:-164.6398404, "
                + "54.7113202:-164.6374015"); // forced vertex, on edge 4
    S1Angle tolerance = S1Angle.radians(2.138358e-05); // 136.235m
    S1Angle maxDist = S1Angle.radians(2.821947e-09); // 18mm
    S2Polygon simplifiedInCell = new S2Polygon();
    simplifiedInCell.initToSimplifiedInCell(polygon, cell, tolerance, maxDist);
    assertFalse(simplifiedInCell.isEmpty());
  }

  @Test
  public void testInitToSimplifiedInCellInteriorEdgesSnappedToBoundary() {
    S2Polygon polygon =
        makePolygonOrDie(
            "37.8011672:-122.3247322, 37.8011648:-122.3247399, "
                + "37.8011647:-122.3247403, 37.8011646:-122.3247408, "
                + "37.8011645:-122.3247411, 37.8011633:-122.3247449, "
                + "37.8011621:-122.3247334");
    S2Cell cell = new S2Cell(S2CellId.fromDebugString("4/001013300"));
    S1Angle snapRadius = metersToAngle(1.0);
    S1Angle boundaryTolerance =
        S1Angle.radians(0.5 * S2Projections.MAX_WIDTH.getValue(S2CellId.MAX_LEVEL - 1))
            .add(IntLatLngSnapFunction.minSnapRadiusForExponent(7));
    S2Polygon simplifiedPolygon = new S2Polygon();
    // simplifiedPolygon.set_s2debug_override(S2Debug::DISABLE);
    simplifiedPolygon.initToSimplifiedInCell(polygon, cell, snapRadius, boundaryTolerance);
    S2Error error = new S2Error();
    boolean invalid = simplifiedPolygon.findValidationError(error);
    assertFalse(error.text(), invalid);
  }

  @Test
  public void testInitToSimplifiedMaxVertices() {
    for (int n = 4; n <= 1024; n *= 2) {
      S2Polygon p = new S2Polygon(S2Loop.makeRegularLoop(S2Point.X_POS, S1Angle.degrees(1), n));
      S2Polygon simplified = p.simplify(24);
      assertLessOrEqual(simplified.getNumVertices(), 24);
      int level = simplified.getSnapLevel();
      assertEquals(n <= 24, level < 0);
      if (level >= 0) {
        S1Angle r = level < 0 ? S1Angle.ZERO : new S2CellIdSnapFunction(level).snapRadius();
        assertTrue(simplified.boundaryNear(p, r.radians()));
      }
    }
  }

  /**
   * Verifies fromCellUnionBorder. The main thing to check is that adjacent cells of different sizes
   * get merged correctly. To do this we generate two random adjacent cells, convert to polygon, and
   * make sure the polygon only has a single loop.
   */
  @Test
  public void testFromCellUnionBorder() {
    for (int iter = 0; iter < 200; ++iter) {
      // Choose a random non-leaf cell.
      S2CellId bigCell = data.getRandomCellId(data.uniform(S2CellId.MAX_LEVEL));
      // Get all neighbors at some smaller level.
      int smallLevel =
          data.uniformInt(bigCell.level(), min(bigCell.level() + 16, S2CellId.MAX_LEVEL));
      List<S2CellId> neighbors = new ArrayList<>();
      bigCell.getAllNeighbors(smallLevel, neighbors);
      // Pick one at random.
      S2CellId smallCell = neighbors.get(data.uniform(neighbors.size()));
      // If it's diagonally adjacent, bail out.
      S2CellId[] edgeNeighbors = new S2CellId[4];
      bigCell.getEdgeNeighbors(edgeNeighbors);
      boolean diagonal = true;
      for (int i = 0; i < 4; i++) {
        if (edgeNeighbors[i].contains(smallCell)) {
          diagonal = false;
        }
      }
      if (diagonal) {
        continue;
      }

      S2CellUnion cellUnion = new S2CellUnion();
      cellUnion.initFromIds(ImmutableList.of(bigCell.id(), smallCell.id()));
      assertEquals(2, cellUnion.size());
      S2Polygon poly = S2Polygon.fromCellUnionBorder(cellUnion);
      assertEquals(1, poly.numLoops());

      // If the conversion were perfect we could test containment, but due to rounding the polygon
      // won't always exactly contain both cells. We can at least test intersection.
      assertTrue(poly.mayIntersect(new S2Cell(bigCell)));
      assertTrue(poly.mayIntersect(new S2Cell(smallCell)));
    }
  }

  @Test
  public void testGetLoops() {
    // These loops share an edge (0:36, 0:39), so the polygon is actually invalid.
    S2Loop loopA = makeLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30");
    S2Loop loopB = makeLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34");
    ImmutableList<S2Loop> expected = ImmutableList.of(loopA, loopB);
    // Create a new mutable list because the constructor clears its input.
    List<S2Loop> listCopy = new ArrayList<>(expected);
    List<S2Loop> actual = uncheckedCreate(() -> new S2Polygon(listCopy)).getLoops();
    assertEquals(2, actual.size());
    assertTrue(actual.contains(loopA));
    assertTrue(actual.contains(loopB));
  }

  /**
   * Creates a polygon from multiple loops specified as comma separated lists of u:v coordinates
   * relative to a cell.
   */
  private static S2Polygon makeCellPolygon(S2Cell cell, String... strs) {
    List<S2Loop> loops = new ArrayList<>();
    R2Rect uv = cell.getBoundUV();

    for (String str : strs) {
      List<S2LatLng> points = S2TextFormat.parseLatLngsOrDie(str); // Actually u/vs, not lat/lngs
      List<S2Point> loopVertices = Lists.newArrayList();
      for (S2LatLng p : points) {
        double u = p.latDegrees();
        double v = p.lngDegrees();
        loopVertices.add(
            S2Point.normalize(
                S2Projections.faceUvToXyz(
                    cell.face(),
                    uv.x().lo() * (1 - u) + uv.x().hi() * u,
                    uv.y().lo() * (1 - v) + uv.y().hi() * v)));
      }
      loops.add(new S2Loop(loopVertices));
    }

    return new S2Polygon(loops);
  }

  /**
   * Verifies that clipBoundary can succeed with duplicate adjacent vertices. Although such a case
   * means the polygon is invalid, it is common to fix invalidity issues by doing a self-
   * intersection to node crossings and drop duplicates.
   */
  @Test
  public void testDuplicatePointClipping() {
    S2Polygon invalid = uncheckedCreate(() -> makePolygonOrDie("0:0, 0:0, 0:4, 4:4, 4:0"));
    S2Polygon fixed = new S2Polygon();
    fixed.initToIntersectionOld(invalid, invalid);
    S2Polygon valid = makePolygon("0:0, 0:4, 4:4, 4:0");
    assertEquals(valid, fixed);
  }

  /**
   * Create 'numLoops' nested regular loops around a common center point. All loops have the given
   * number of vertices 'numVertices'. Furthermore, the vertices at the same index position are
   * collinear with the common center point of all the loops. The loop radii decrease exponentially
   * in order to prevent accidental loop crossings when one of the loops is modified.
   */
  private List<List<S2Point>> getConcentricLoops(int numLoops, int minVertices) {
    List<List<S2Point>> loops = Lists.newArrayList();
    // Radii decrease exponentially.
    assertTrue(numLoops <= 10);
    S2Point center = data.getRandomPoint();
    int numVertices = minVertices + data.random(10);
    for (int i = 0; i < numLoops; ++i) {
      S1Angle radius = S1Angle.degrees(80 * pow(0.1, i));
      loops.add(S2Loop.makeRegularVertices(center, radius, numVertices));
    }
    return loops;
  }

  @Test
  public void testSplitting() {
    // It takes too long to test all the polygons in debug mode, so we just pick out some of the
    // more interesting ones.
    splitAndAssemble(near10);
    splitAndAssemble(nearH3210);
    splitAndAssemble(south0ab);
    splitAndAssemble(south210b);
    splitAndAssemble(nf1n10f2s10abc);
    splitAndAssemble(nf2n2f210s210ab);
    splitAndAssemble(farHSouthH);
    splitAndAssemble(southH);
    splitAndAssemble(farH);
    splitAndAssemble(southH20abc);
    splitAndAssemble(farH3210);
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++
  @Test
  public void testUnionWithAmbgiuousCrossings() {
    ImmutableList<S2Point> aVertices =
        ImmutableList.of(
            p(0.044856812877680216, -0.80679210859571904, 0.5891301722422051),
            p(0.044851868273159699, -0.80679240802900054, 0.5891301386444033),
            p(0.044854246527738666, -0.80679240292188514, 0.58912996457145106));
    ImmutableList<S2Point> bVertices =
        ImmutableList.of(
            p(0.044849715793028468, -0.80679253837178111, 0.58913012401412856),
            p(0.044855344598821352, -0.80679219751320641, 0.589130162266992),
            p(0.044854017712818696, -0.80679210327223405, 0.58913039235179754));
    S2Polygon a = new S2Polygon(new S2Loop(aVertices));
    S2Polygon b = new S2Polygon(new S2Loop(bVertices));
    S2Polygon c = new S2Polygon();
    c.initToUnion(a, b);
    assertFalse(c.isEmpty());
  }

  private void splitAndAssemble(S2Polygon polygon) {
    S2PolygonBuilder builder = new S2PolygonBuilder(S2PolygonBuilder.Options.DIRECTED_XOR);
    S2Polygon expected = new S2Polygon();
    builder.addPolygon(polygon);
    assertTrue(builder.assemblePolygon(expected, null));
    checkEqual(polygon, expected);

    for (int iter = 0; iter < 10; ++iter) {
      // Compute the minimum level such that the polygon's bounding cap is guaranteed to be cut.
      double diameter = 2 * polygon.getCapBound().angle().radians();
      int minLevel = MAX_WIDTH.getMinLevel(diameter);

      // Now choose a level that has up to 500 cells in the covering.
      int level = minLevel + data.random(6);
      S2RegionCoverer coverer =
          S2RegionCoverer.builder()
              .setMinLevel(minLevel)
              .setMaxLevel(level)
              .setMaxCells(500)
              .build();

      ArrayList<S2CellId> cells = Lists.newArrayList();
      coverer.getCovering(polygon, cells);
      S2CellUnion covering = new S2CellUnion().initFromCellIds(cells);
      checkCovering(polygon, covering, false);
      checkCoveringIsConservative(polygon, cells);
      logger.fine(cells.size() + " cells in covering");
      List<S2Polygon> pieces = Lists.newArrayList();
      for (int i = 0; i < cells.size(); ++i) {
        S2Cell cell = new S2Cell(cells.get(i));
        S2Polygon window = new S2Polygon(cell);
        S2Polygon piece = new S2Polygon();
        piece.initToIntersection(polygon, window);
        pieces.add(piece);
        logger.fine("\nPiece " + i + ":\n Window: " + window + "\n Piece: " + piece);
      }

      // Now we repeatedly remove two random pieces, compute their union, and insert the result as a
      // new piece until only one piece is left.
      //
      // We don't use S2Polygon.union() because it joins the pieces in a mostly deterministic order.
      // We don't just randomly shuffle the pieces and repeatedly join the last two pieces because
      // this always joins a single original piece to the current union rather than doing the unions
      // according to some random tree structure.
      while (pieces.size() > 1) {
        S2Polygon a = choosePiece(pieces);
        S2Polygon b = choosePiece(pieces);
        S2Polygon c = new S2Polygon();
        c.initToUnion(a, b);
        pieces.add(c);
        logger.fine(
            "\nSplitAndAssemble Joined piece a: "
                + S2TextFormat.toString(a)
                + "\n With piece b: "
                + S2TextFormat.toString(b)
                + "\n To get piece c: "
                + S2TextFormat.toString(c));
      }
      S2Polygon result = pieces.get(0);

      // The moment of truth!
      assertTrue(
          "\nActual:\n" + result + "\nExpected:\n" + expected,
          expected.boundaryNear(result, 1e-15));
    }
  }

  /** Removes a random polygon from {@code pieces} and returns it. */
  private S2Polygon choosePiece(List<S2Polygon> pieces) {
    int i = data.random(pieces.size());
    S2Polygon result = pieces.get(i);
    pieces.remove(i);
    return result;
  }

  /**
   * Checks that contains(S2Cell) and mayIntersect(S2Cell) are implemented conservatively, by
   * comparing against the contains/intersect result with the 'cell polygon' defined by the four
   * cell vertices. Please note that the cell polygon is *not* an exact representation of the
   * S2Cell: cell vertices are rounded from their true mathematical positions, which leads to tiny
   * cracks and overlaps between the cell polygons at different cell levels. That is why
   * contains(S2Cell) and mayIntersect(S2Cell) cannot be implemented by simply converting the cell
   * to an S2Polygon. But it is still useful to do this sanity check. In particular:
   *
   * <ul>
   *   <li>If contains(cell) is true, the polygon must contain the cell polygon.
   *   <li>If the polygon intersects the cell polygon, then mayIntersect(cell) must return true.
   * </ul>
   */
  private static void checkCoveringIsConservative(S2Polygon polygon, List<S2CellId> cells) {
    for (int i = 0; i < cells.size(); ++i) {
      S2Cell cell = new S2Cell(cells.get(i));
      S2Polygon cellPoly = new S2Polygon(cell);
      if (polygon.contains(cell)) {
        assertTrue(polygon.contains(cellPoly));
      }
      if (polygon.intersects(cellPoly)) {
        assertTrue(polygon.mayIntersect(cell));
      }
    }
  }

  /** Checks that the index is properly deserialized. */
  @GwtIncompatible("Object serialization")
  @Test
  public void testIndexDeserialization() throws IOException, ClassNotFoundException {
    S2Point center = S2Point.X_POS;
    S1Angle angle = S1Angle.radians(10);
    int numVertices = 10;
    S2Polygon polygon = new S2Polygon(S2Loop.makeRegularLoop(center, angle, numVertices));
    // Initialize the index.
    S2Iterator<S2ShapeIndex.Cell> unused = polygon.index.iterator();
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    ObjectOutputStream out = new ObjectOutputStream(output);

    // Serialize the object.
    out.writeObject(polygon);
    out.close();

    // Deserialize the object.
    ByteArrayInputStream input = new ByteArrayInputStream(output.toByteArray());
    ObjectInputStream in = new ObjectInputStream(input);
    S2Polygon copy = (S2Polygon) in.readObject();
    in.close();

    assertTrue(copy.index != null);
    assertTrue(copy.index.iterator().locate(S2Point.X_POS));
    assertTrue(copy.index.shapes.size() == 1);
    assertTrue(copy.getNumVertices() == numVertices);
    assertTrue(copy.contains(S2Point.X_POS));
  }

  @Test
  public void testEmptyPolygonShape() {
    S2Polygon.Shape shape = empty.shape();
    assertTrue(shape.isEmpty());
    assertFalse(shape.isFull());
    assertEquals(empty, shape.polygon());
    assertTrue(shape.hasInterior());
    assertFalse(shape.containsOrigin());
    assertEquals(0, shape.numEdges());
    assertEquals(0, shape.numChains());
    assertEquals(2, shape.dimension());
  }

  @Test
  public void testFullPolygonShape() {
    S2Polygon.Shape shape = full.shape();
    assertFalse(shape.isEmpty());
    assertTrue(shape.isFull());
    assertEquals(full, shape.polygon());
    assertTrue(shape.hasInterior());
    assertTrue(shape.containsOrigin());
    assertEquals(0, shape.numEdges());
    assertEquals(1, shape.numChains());
    checkFirstNChainStarts(shape, 0);
    checkFirstNChainLengths(shape, 0);
    assertEquals(2, shape.dimension());
  }

  /** Verify that polygons remove their empty loops during construction. */
  @Test
  public void testEmptyLoopsRemoved() throws IOException {
    for (int iter = 0; iter < VALIDITY_ITERS; ++iter) {
      int numLoops = data.uniformInt(1, 6);
      List<List<S2Point>> loops = getConcentricLoops(numLoops, /* minVertices= */ 3);

      // Insert one or more empty loops into the list at random positions, which should be ignored.
      int nEmptyLoops = data.uniformInt(1, 6);
      for (int j = 0; j < nEmptyLoops; ++j) {
        int insertPos = data.uniformInt(0, loops.size());
        loops.add(insertPos, S2Loop.empty().vertices());
      }

      S2Polygon polygon = new S2Polygon(uncheckedLoops(loops));
      assertTrue(polygon.isValid());
      assertEquals(numLoops, polygon.numLoops());
    }
  }

  @Test
  public void testEmptyLoopsRemovedDuringDecoding() throws IOException {
    // Verify that polygons remove serialized canonical empty loops during decoding. Note that
    // S2Polygon.decode guarantees it will not return an invalid polygon, so other kinds of invalid
    // loops cannot be included in this test. That includes having a full loop along with other
    // loops. Also, the remaining loop must have a valid loop depth, but fromExplicitLoops() doesn't
    // set the depth, so we need to fix that ourselves.
    // Reset the loop depth which was modified above.
    farH.loop(0).setDepth(0);
    S2Polygon polygon = S2Polygon.fromExplicitLoops(ImmutableList.of(farH.loop(0), S2Loop.empty()));
    assertTrue(polygon.getLoops().stream().anyMatch(S2Loop::isEmpty));

    ByteArrayOutputStream encoded = new ByteArrayOutputStream();
    polygon.encode(encoded);
    S2Polygon decoded = S2Polygon.decode(new ByteArrayInputStream(encoded.toByteArray()));

    assertEquals(0, farH.loop(0).depth());
    assertEquals(ImmutableSet.of(farH.loop(0)), ImmutableSet.copyOf(decoded.getLoops()));
  }

  @Test
  public void testOneLoopPolygonShape() {
    S2Polygon.Shape shape = near0.shape();
    assertFalse(shape.isEmpty());
    assertFalse(shape.isFull());
    assertEquals(near0, shape.polygon());
    assertTrue(shape.hasInterior());
    assertFalse(shape.containsOrigin());
    assertEquals(4, shape.numEdges());
    checkGetEdge(near0, shape);
    assertEquals(1, shape.numChains());
    checkFirstNChainStarts(shape, 0);
    checkFirstNChainLengths(shape, 4);
    assertEquals(2, shape.dimension());
  }

  @Test
  public void testSeveralLoopPolygonShape() {
    S2Polygon.Shape shape = near3210.shape();
    assertFalse(shape.isEmpty());
    assertFalse(shape.isFull());
    assertEquals(near3210, shape.polygon());
    assertTrue(shape.hasInterior());
    assertFalse(shape.containsOrigin());
    assertEquals(18, shape.numEdges());
    checkGetEdge(near3210, shape);
    assertEquals(4, shape.numChains());
    checkFirstNChainStarts(shape, 0, 3, 6, 14);
    checkFirstNChainLengths(shape, 3, 3, 8, 4);
    assertEquals(2, shape.dimension());
  }

  @Test
  public void testManyLoopPolygonShape() {
    int numLoops = 100;
    int numVerticesPerLoop = 6;
    S2Polygon polygon = concentricLoopsPolygon(S2Point.X_POS, numLoops, numVerticesPerLoop);
    S2Polygon.Shape shape = polygon.shape();
    assertEquals(polygon, shape.polygon());
    assertTrue(shape.hasInterior());
    assertFalse(shape.containsOrigin());
    assertEquals(600, shape.numEdges());
    checkGetEdge(polygon, shape);
    assertEquals(100, shape.numChains());
    checkFirstNChainStarts(shape, 0, 6, 12, 18, 24);
    checkFirstNChainLengths(shape, 6, 6, 6, 6, 6);
    assertEquals(2, shape.dimension());
  }

  private static void checkGetEdge(S2Polygon polygon, S2Shape shape) {
    MutableEdge edge = new MutableEdge();
    ChainPosition position = new ChainPosition();
    for (int e = 0, i = 0; i < polygon.numLoops(); ++i) {
      S2Loop loop = polygon.loop(i);
      for (int j = 0; j < loop.numVertices(); ++j, ++e) {
        if (!loop.isEmptyOrFull()) {
          shape.getEdge(e, edge);
          assertEquals(loop.orientedVertex(j), edge.a);
          assertEquals(loop.orientedVertex(j + 1), edge.b);
          shape.getChainEdge(i, j, edge);
          assertEquals(loop.orientedVertex(j), edge.a);
          assertEquals(loop.orientedVertex(j + 1), edge.b);
          // Assert the consistency of getChainPosition with getChainEdge and getEdge.
          shape.getChainPosition(e, position);
          assertEquals(position.getChainId(), i);
          assertEquals(position.getOffset(), j);
        }
      }
    }
  }

  @Test
  public void testInitRecursion() {
    String loop = "-18.84:-40.96, -18.93:-40.96, -18.93:-40.86, -18.84:-40.86";
    S2Polygon duplicateLoops = uncheckedCreate(() -> makePolygonOrDie(loop + "; " + loop));
    // This could be POLYGON_LOOPS_SHARE_EDGE, but S2ValidQuery checks for duplicate edges by cell
    // without tracking which loops own the edges, and categorizes it as OVERLAPPING_GEOMETRY.
    checkInvalid(duplicateLoops, S2Error.Code.OVERLAPPING_GEOMETRY);
  }

  /** Verifies a bug in S2ShapeIndex has been fixed. */
  @Test
  public void testPointInBigLoop() {
    S2LatLng center = S2LatLng.fromRadians(0.3, 2);
    S2Loop loop = S2Loop.makeRegularLoop(center.toPoint(), S1Angle.degrees(80), 10);
    assertTrue(new S2Polygon(loop).mayIntersect(new S2Cell(S2CellId.fromLatLng(center))));
  }

  /** Verifies S2Polygon.getCentroid() can handle polygons with empty loops. */
  @Test
  public void testEmptyLoopCentroid() {
    S2Polygon polygon = new S2Polygon(S2Loop.empty());
    assertEquals(polygon.getCentroid(), S2Point.ZERO);
  }

  @Test
  public void testEncodeDecodeLossless() throws Exception {
    S2Polygon polygon = makePolygon("1:2, 3:4, 5:6; 7:8, 9:10, 11:12");

    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    LittleEndianOutput encoder = new LittleEndianOutput(baos);
    encoder.writeByte((byte) 1); // S2Polygon.LOSSLESS_ENCODING_VERSION

    // These two bytes are ignored. Respectively encoding that owns loop is false and there are
    // no holes.
    encoder.writeByte((byte) 0);
    encoder.writeByte((byte) 0);

    encoder.writeInt(2);
    polygon.loop(0).encode(encoder);
    polygon.loop(1).encode(encoder);
    polygon.getRectBound().encode(encoder);

    ByteArrayInputStream in = new ByteArrayInputStream(baos.toByteArray());
    S2Polygon decodedPolygon = S2Polygon.decode(in);
    assertEquals(polygon, decodedPolygon);
    assertEquals(0, in.available());
  }

  @Test
  public void testDecodeEncodeSingleLoopLossless() throws IOException {
    String encodedBytesHexString =
        "010100010000000108000000D44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA1BFB4"
            + "825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF91BFB4825F3C81FDEF3F27DCF7C958DE"
            + "913F1EDD892B0BDF913FD44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA13FD44A84"
            + "42C3F9EF3F7EDA2AB341DC91BF27DCF7C958DEA13FB4825F3C81FDEF3F27DCF7C958DE91BF"
            + "1EDD892B0BDF913FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91BFD44A8442C3"
            + "F9EF3F7EDA2AB341DC91BF27DCF7C958DEA1BF0000000000013EFC10E8F8DFA1BF3EFC10E8"
            + "F8DFA13F389D52A246DF91BF389D52A246DF913F013EFC10E8F8DFA1BF3EFC10E8F8DFA13F"
            + "389D52A246DF91BF389D52A246DF913F";
    ByteArrayInputStream bais =
        new ByteArrayInputStream(BaseEncoding.base16().decode(encodedBytesHexString));
    S2Polygon decodedPolygon = S2Polygon.decode(bais);
    assertEquals(1, decodedPolygon.numLoops());
    assertEquals(
        makeLoop("-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1"), decodedPolygon.loop(0));
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    decodedPolygon.encode(bos);
    assertEquals(encodedBytesHexString, BaseEncoding.base16().encode(bos.toByteArray()));
  }

  @Test
  public void testEncodeDecodeEmptyLossless() throws IOException {
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    new S2Polygon().encode(bos);
    assertEquals("041E00", BaseEncoding.base16().encode(bos.toByteArray()));
  }

  @Test
  public void testEncodeDecodeCorruptLossless() {
    IOException e =
        assertThrows(
            IOException.class,
            () ->
                S2Polygon.decode(
                    new ByteArrayInputStream(BaseEncoding.base16().decode("0101FFFFFF"))));
    assertEquals("EOF", e.getMessage());
  }

  @Test
  public void testEncodeDecodeTwoPolygons() throws IOException {
    S2Polygon p1 = makePolygon("10:10,10:0,0:0");
    S2Polygon p2 = makePolygon("0:0,1:0,1:1,0:1");
    ByteArrayOutputStream out = new ByteArrayOutputStream();
    p1.encode(out);
    p2.encode(out);

    byte[] data = out.toByteArray();
    ByteArrayInputStream in = new ByteArrayInputStream(data);
    S2Polygon p1Decoded = S2Polygon.decode(in);
    S2Polygon p2Decoded = S2Polygon.decode(in);
    assertEquals(p1, p1Decoded);
    assertEquals(p2, p2Decoded);
    assertEquals(0, in.available());
  }

  @Test
  public void testDecodeEncodeTwoLoopsLossless() throws IOException {
    // Two counter-clockwise loops, one nested inside the other.
    String encodedBytesHexString =
        "010101020000000108000000D44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA1BFB4"
            + "825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF91BFB4825F3C81FDEF3F27DCF7C958DE"
            + "913F1EDD892B0BDF913FD44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA13FD44A84"
            + "42C3F9EF3F7EDA2AB341DC91BF27DCF7C958DEA13FB4825F3C81FDEF3F27DCF7C958DE91BF"
            + "1EDD892B0BDF913FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91BFD44A8442C3"
            + "F9EF3F7EDA2AB341DC91BF27DCF7C958DEA1BF0000000000013EFC10E8F8DFA1BF3EFC10E8"
            + "F8DFA13F389D52A246DF91BF389D52A246DF913F0104000000C5D7FA4B60FFEF3F1EDD892B"
            + "0BDF813F214C95C437DF81BFC5D7FA4B60FFEF3F1EDD892B0BDF813F214C95C437DF813FC5"
            + "D7FA4B60FFEF3F1EDD892B0BDF81BF214C95C437DF813FC5D7FA4B60FFEF3F1EDD892B0BDF"
            + "81BF214C95C437DF81BF000100000001900C5E3B73DF81BF900C5E3B73DF813F399D52A246"
            + "DF81BF399D52A246DF813F013EFC10E8F8DFA1BF3EFC10E8F8DFA13F389D52A246DF91BF38"
            + "9D52A246DF913F";
    ByteArrayInputStream bais =
        new ByteArrayInputStream(BaseEncoding.base16().decode(encodedBytesHexString));
    S2Polygon decodedPolygon = S2Polygon.decode(bais);
    assertEquals(2, decodedPolygon.numLoops());
    assertEquals(
        makeLoop("-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1"), decodedPolygon.loop(0));
    assertEquals(makeLoop("-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5"), decodedPolygon.loop(1));
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    decodedPolygon.encode(bos);
    assertEquals(encodedBytesHexString, BaseEncoding.base16().encode(bos.toByteArray()));
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2PolygonSerialization() {
    // This used to be a loop from the cell at point S2Point(0.1, 0.2, 0.3) but that causes
    // a POLYGON_LOOPS_CROSS error.
    S2Loop loop1 = new S2Loop(new S2Cell(new S2Point(-0.2, 0.3, 0.4)));
    S2Loop loop2 =
        new S2Loop(
            ImmutableList.of(
                new S2Point(0.1, 0.2, 0.3).normalize(),
                new S2Point(0.5, 0, 0).normalize(),
                new S2Point(4, 4, 4).normalize(),
                new S2Point(6, 7, 8).normalize()));
    doSerializationTest(
        new S2Polygon(Lists.newArrayList(loop1, loop2)), Ordering.<S2Polygon>natural());
  }

  private static S2Point p(double x, double y, double z) {
    return new S2Point(x, y, z);
  }
}
