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

import static com.google.common.geometry.R2EdgeClipper.BOTTOM;
import static com.google.common.geometry.R2EdgeClipper.INSIDE;
import static com.google.common.geometry.R2EdgeClipper.LEFT;
import static com.google.common.geometry.R2EdgeClipper.OUTSIDE;
import static com.google.common.geometry.R2EdgeClipper.RIGHT;
import static com.google.common.geometry.R2EdgeClipper.TOP;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Verifies {@link UVEdgeClipper}. */
@RunWith(JUnit4.class)
public class UVEdgeClipperTest extends GeometryTestCase {
  // Cell's 0 and 1 share a U boundary (left and right).
  // Cell's 0 and 2 share a V boundary (top and bottom).
  private static final S2Cell kCellToken0 = new S2Cell(S2CellId.fromToken("89c25c1"));
  private static final S2Cell kCellToken1 = new S2Cell(S2CellId.fromToken("89c25c3"));
  private static final S2Cell kCellToken2 = new S2Cell(S2CellId.fromToken("89c25c7"));
  private final UVEdgeClipper clipper = new UVEdgeClipper();

  @Before
  public void init() {
    clipper.init(kCellToken0);
  }

  @Test
  public void testAllPermutationsWork() {
    // Convert to UV space for comparison with clipped result.
    R2Vector raw0 = new R2Vector();
    R2Vector raw1 = new R2Vector();
    for (var testCase : buildTestCases()) {
      S2Projections.validFaceXyzToUv(kCellToken0.face(), testCase.v0, raw0);
      S2Projections.validFaceXyzToUv(kCellToken0.face(), testCase.v1, raw1);

      R2Rect boundary = kCellToken0.getBoundUV();

      if (clipper.clipEdge(testCase.v0, testCase.v1)) {
        assertEquals(INSIDE, testCase.outcode0 & testCase.outcode1);

        // Clipped vertices should not land on the original vertices unless they were inside the
        // clipping region to begin with.
        if (testCase.outcode0 != INSIDE) {
          assertNotEquals(raw0, clipper.clippedUvEdge().v0);
          assertTrue(onBoundary(clipper.clippedUvEdge().v0, boundary));
          assertNotEquals(INSIDE, clipper.outcode(0));
        } else {
          assertEquals(raw0, clipper.clippedUvEdge().v0);
          assertEquals(INSIDE, clipper.outcode(0));
        }

        if (testCase.outcode1 != INSIDE) {
          assertNotEquals(raw1, clipper.clippedUvEdge().v1);
          assertTrue(onBoundary(clipper.clippedUvEdge().v1, boundary));
          assertNotEquals(INSIDE, clipper.outcode(1));
        } else {
          assertEquals(raw1, clipper.clippedUvEdge().v1);
          assertEquals(INSIDE, clipper.outcode(1));
        }

        // Clipped edges should be co-linear with the original edge.
        checkColinear(raw0, clipper.clippedUvEdge().v0, clipper.clippedUvEdge().v1);
        checkColinear(raw1, clipper.clippedUvEdge().v0, clipper.clippedUvEdge().v1);
      } else {
        // Not clipped, the two vertices must have been in the same exterior region.
        assertNotEquals(INSIDE, testCase.outcode0 & testCase.outcode1);
        assertEquals(OUTSIDE, clipper.outcode(0));
        assertEquals(OUTSIDE, clipper.outcode(1));
      }
    }
  }

  /** Tests that an edge that crosses a cell boundary clips to the exact same point in each cell. */
  @Test
  public void testSameResultAcrossCellBoundary() {
    // Generate an edge crossing the boundary between the cells. Note that face 4 is rotated so the
    // U coordinate is in the y axis and the V coordinate is in the x axis.
    R2Rect boundary0 = kCellToken0.getBoundUV();
    R2Rect boundary1 = kCellToken1.getBoundUV();
    R2Rect boundary2 = kCellToken2.getBoundUV();
    assertEquals(boundary0.hi().y(), boundary1.lo().y(), 0);
    assertEquals(boundary0.hi().x(), boundary2.lo().x(), 0);

    double kBoundaryU = boundary0.hi().y();
    double kBoundaryV = boundary0.hi().x();

    for (int i = 0; i < 100; ++i) {
      R2Vector uv0 = randPointInBoundary(boundary0);
      R2Vector uv1 = randPointInBoundary(boundary1);
      R2Vector uv2 = randPointInBoundary(boundary2);

      S2Point pnt0 = S2Projections.faceUvToXyz(kCellToken0.face(), uv0.x(), uv0.y());
      S2Point pnt1 = S2Projections.faceUvToXyz(kCellToken0.face(), uv1.x(), uv1.y());
      S2Point pnt2 = S2Projections.faceUvToXyz(kCellToken0.face(), uv2.x(), uv2.y());

      // Test crossing the U boundary for consistency. Note that we check for exact equality here,
      // which should always be true.
      clipper.init(kCellToken0);
      assertTrue(clipper.clipEdge(pnt0, pnt1));
      assertEquals(kBoundaryU, clipper.clippedUvEdge().v1.y(), 0);

      clipper.init(kCellToken1);
      assertTrue(clipper.clipEdge(pnt0, pnt1));
      assertEquals(kBoundaryU, clipper.clippedUvEdge().v0.y(), 0);

      // Test crossing the V boundary for consistency. Note that we check for exact equality here,
      // which should always be true.
      clipper.init(kCellToken0);
      assertTrue(clipper.clipEdge(pnt0, pnt2));
      assertEquals(kBoundaryV, clipper.clippedUvEdge().v1.x(), 0);

      clipper.init(kCellToken2);
      assertTrue(clipper.clipEdge(pnt0, pnt2));
      assertEquals(kBoundaryV, clipper.clippedUvEdge().v0.x(), 0);
    }
  }

  /** Returns a random point within the given cell boundary. */
  private R2Vector randPointInBoundary(R2Rect boundary) {
    double u = data.uniform(boundary.lo().x(), boundary.hi().x());
    double v = data.uniform(boundary.lo().y(), boundary.hi().y());
    return new R2Vector(u, v);
  }

  private static S2Point toXyz(int face, R2Vector uv) {
    return S2Projections.faceUvToXyz(face, uv.x(), uv.y());
  }

  /**
   * Tests multiple iterations of clipping connected edges to ensure that using 'connected' returns
   * the same results as clipping the edges individually.
   */
  @Test
  public void testClipConnectedEdges() {
    // A rectangle to sample edge endpoints from.
    R2Rect parentRect = kCellToken0.id().parent().getBoundUV();
    int face = kCellToken0.face();

    // Two clippers, one of which will use the "connected" optimization.
    clipper.init(kCellToken0);
    UVEdgeClipper connectedClipper = new UVEdgeClipper(kCellToken0);

    // Create and clip the initial edge. */
    S2Point first = toXyz(face, randPointInBoundary(parentRect));
    S2Point prev = toXyz(face, randPointInBoundary(parentRect));

    assertEquals(
        clipper.clipEdge(first, prev, false), connectedClipper.clipEdge(first, prev, false));
    assertClipResultsEqual(clipper, connectedClipper);

    // Repeatedly generate another point and clip the edge from the previous point to the new one.
    for (int i = 0; i < 100; ++i) {
      S2Point next = toXyz(face, randPointInBoundary(parentRect));
      assertEquals(
          clipper.clipEdge(prev, next, false), connectedClipper.clipEdge(prev, next, true));
      assertClipResultsEqual(clipper, connectedClipper);
      prev = next;
    }
  }

  /**
   * Asserts that the two UVEdgeClippers produce the same results when clipping the same edge.
   * The first clipper clips edges individually, while the second clips connected edges, so must be
   * called with a sequence of connected edges.
   */
  private void assertClipResultsEqual(UVEdgeClipper clipper, UVEdgeClipper connectedClipper) {
    R2Edge clipped = clipper.clippedUvEdge();
    R2Edge connectedClipped = connectedClipper.clippedUvEdge();
    assertTrue("Clipped = " + clipped + ", connectedClipped = " + connectedClipped,
        clipped.isEqualTo(connectedClipped));
    assertEquals(clipper.outcode(0), connectedClipper.outcode(0));
    assertEquals(clipper.outcode(1), connectedClipper.outcode(1));
  }

  /** Tests that vertices exactly on the boundary don't get modified. */
  @Test
  public void testPointsOnBoundaryUnchanged() {
    R2Rect rect = kCellToken0.getBoundUV();
    R2Vector uv0 = new R2Vector(rect.lo().x(), rect.lo().y());
    R2Vector uv1 = new R2Vector(rect.hi().x(), rect.lo().y());
    R2Vector uv2 = new R2Vector(rect.hi().x(), rect.hi().y());
    R2Vector uv3 = new R2Vector(rect.lo().x(), rect.hi().y());
    S2Point pnt0 = S2Projections.faceUvToXyz(kCellToken0.face(), uv0.x(), uv0.y());
    S2Point pnt1 = S2Projections.faceUvToXyz(kCellToken0.face(), uv1.x(), uv1.y());
    S2Point pnt2 = S2Projections.faceUvToXyz(kCellToken0.face(), uv2.x(), uv2.y());
    S2Point pnt3 = S2Projections.faceUvToXyz(kCellToken0.face(), uv3.x(), uv3.y());

    clipper.init(kCellToken0);

    assertTrue(clipper.clipEdge(pnt0, pnt1));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv0, uv1));
    assertEdgesEqual(clipper.clippedUvEdge(), clipper.faceUvEdge());
    assertEquals(INSIDE, clipper.outcode(0));
    assertEquals(INSIDE, clipper.outcode(1));

    assertTrue(clipper.clipEdge(pnt1, pnt2));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv1, uv2));
    assertEdgesEqual(clipper.clippedUvEdge(), clipper.faceUvEdge());
    assertEquals(INSIDE, clipper.outcode(0));
    assertEquals(INSIDE, clipper.outcode(1));

    assertTrue(clipper.clipEdge(pnt2, pnt3));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv2, uv3));
    assertEdgesEqual(clipper.clippedUvEdge(), clipper.faceUvEdge());
    assertEquals(INSIDE, clipper.outcode(0));
    assertEquals(INSIDE, clipper.outcode(1));

    assertTrue(clipper.clipEdge(pnt3, pnt0));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv3, uv0));
    assertEdgesEqual(clipper.clippedUvEdge(), clipper.faceUvEdge());
    assertEquals(INSIDE, clipper.outcode(0));
    assertEquals(INSIDE, clipper.outcode(1));
  }

  /** Tests that line segments coincident with the boundary are clipped to it. */
  @Test
  public void testLineOnBoundaryClipped() {
    int kFace = kCellToken0.face();

    R2Rect rect = kCellToken0.getBoundUV();
    R2Vector uv0 = new R2Vector(rect.lo().x(), rect.lo().y());
    R2Vector uv1 = new R2Vector(rect.hi().x(), rect.lo().y());
    R2Vector miduv0 = uv0.add(uv1.sub(uv0).mul(1.0 / 3));
    R2Vector miduv1 = uv0.add(uv1.sub(uv0).mul(2.0 / 3));
    S2Point pnt0 = S2Projections.faceUvToXyz(kFace, uv0.x(), uv0.y());
    S2Point pnt1 = S2Projections.faceUvToXyz(kFace, uv1.x(), uv1.y());
    S2Point midpnt0 = S2Projections.faceUvToXyz(kFace, miduv0.x(), miduv0.y());
    S2Point midpnt1 = S2Projections.faceUvToXyz(kFace, miduv1.x(), miduv1.y());

    S2Point hi = pnt0.add(pnt1.sub(pnt0).mul(2));
    S2Point lo = pnt0.sub(pnt1.sub(pnt0).mul(2));

    // Segment extends past vertex 1.
    assertTrue(clipper.clipEdge(pnt0, hi));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv0, uv1));

    // Segment extends past vertex 0.
    assertTrue(clipper.clipEdge(lo, pnt1));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv0, uv1));

    // Segment extends past both.
    assertTrue(clipper.clipEdge(lo, hi));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv0, uv1));

    // Point in the middle of the boundary to vertex 1.
    assertTrue(clipper.clipEdge(midpnt0, pnt1));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(miduv0, uv1));

    // Vertex 0 to point in the middle of the boundary.
    assertTrue(clipper.clipEdge(pnt0, midpnt0));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(uv0, miduv0));

    // Two points inside the boundary segment.
    assertTrue(clipper.clipEdge(midpnt0, midpnt1));
    assertEdgesEqual(clipper.clippedUvEdge(), edge(miduv0, miduv1));
  }

  /** Asserts that the two edges have the same endpoints. */
  private static void assertEdgesEqual(R2Edge expected, R2Edge actual) {
    assertTrue(expected.isEqualTo(actual));
  }

  /** Tests that a line in two exterior regions but not the cell is still clipped away. */
  @Test
  public void testOutsideBetweenDifferentRegionsClipped() {
    // Each of these lines goes from an exterior region of the cell to a different one, without
    // touching the cell itself. They should all get clipped out.
    S2Point[][] edges = {
        {point(40.708180, -73.937992), point(40.710042, -73.934527)},
        {point(40.707250, -73.914758), point(40.703740, -73.909094)},
        {point(40.686950, -73.909737), point(40.684151, -73.915059)},
        {point(40.691128, -73.939047), point(40.687125, -73.933897)},
        {point(40.671643, -73.924799), point(40.741652, -73.976984)}};
    for (S2Point[] edge : edges) {
      assertFalse(clipper.clipEdge(edge[0], edge[1]));
      assertEquals(OUTSIDE, clipper.outcode(0));
      assertEquals(OUTSIDE, clipper.outcode(1));
    }
  }

  /** Verifies an edge spanning faces but intersecting the cell is detected. */
  @Test
  public void testEdgeHitAcrossFacesWorks() {
    S2Point v0 = S2LatLng.fromDegrees(-1.801, -68.044).toPoint();
    S2Point v1 = S2LatLng.fromDegrees(8.0160, +82.425).toPoint();
    clipper.init(new S2Cell(S2CellId.fromToken("11")));
    assertTrue(clipper.clipEdge(v0, v1));
  }

  /** Verifies an edge spanning faces but missing the face entirely is detected. */
  @Test
  public void testEdgeMissedFaceDetected0() {
    // First edge touches face 0, so it'll get re-used. Second vertex moves to face two.
    S2Point[] edge0 = {point(29, -15), point(-2, -97)};
    S2Point[] edge1 = {point(-2, -97), point(-10, -173)};

    clipper.init(new S2Cell(S2CellId.fromToken("11")));
    assertFalse(clipper.clipEdge(edge0[0], edge0[1]));
    assertFalse(clipper.clipEdge(edge1[0], edge1[1]));
    assertTrue(clipper.missedFace());
  }

  /** As {@link #testEdgeMissedFaceDetected0} but the edge goes from face 0 -> 4 and stays there. */
  @Test
  public void testEdgeMissedFaceDetected1() {
    // First edge touches face 0, so it'll get re-used. Second vertex moves to face two.
    S2Point[] edge0 = {point(29, -15), point(-2, -97)};
    S2Point[] edge1 = {point(-2, -97), point(-10, -110)};
    clipper.init(new S2Cell(S2CellId.fromToken("11")));
    assertFalse(clipper.clipEdge(edge0[0], edge0[1]));
    assertFalse(clipper.clipEdge(edge1[0], edge1[1]));
    assertTrue(clipper.missedFace());
  }

  /** Verifies a single edge, on face 4 missing face 0. */
  @Test
  public void testEdgeMissedFaceDetected2() {
    clipper.init(new S2Cell(S2CellId.fromToken("11")));
    S2Point[] edge = {point(-2, -97), point(-10, -110)};
    assertFalse(clipper.clipEdge(edge[0], edge[1]));
    assertTrue(clipper.missedFace());
  }

  /** Verifies an edge that barely grazes a face. */
  @Test
  public void testAlmostMissesFace() {
    // The cell is on face 4 and the edge is on face 2 but touches face 4 at a point. That's close
    // enough that we have to consider it a hit.
    UVEdgeClipper clipper = new UVEdgeClipper();
    clipper.init(new S2Cell(S2CellId.fromToken("800c")));
    S2Point[] kEdge = {
        new S2Point(-0.535533175298229, -0.597161710994181, 0.597161710994181),
        new S2Point(-0.528388331336849, -0.589982219252971, 0.610513515225011)};

    assertTrue(clipper.clipEdge(kEdge[0], kEdge[1]));
    assertEquals(-1.0, clipper.clippedUvEdge().v0.x(), 0);
    assertEquals(-1.0, clipper.clippedUvEdge().v1.x(), 0);
  }

  /** Verifies that if we init(S2Cell) the face and clip boundary are correct. */
  @Test
  public void testFaceAndBoundarySetFromCell() {
    S2Cell kCell = new S2Cell(S2CellId.fromToken("800c"));

    UVEdgeClipper clipper = new UVEdgeClipper();
    clipper.init(kCell);

    assertEquals(4, clipper.clipFace());
    assertEquals(kCell.getBoundUV(), clipper.clipRect());
  }

  /** Test two edges that graze the boundary between faces 3 and 4. */
  @Test
  public void testBouncesOffFace() {
    // Clip two edges ABC such that B is a corner of the clip cell. This kind of "bounce" should be
    // on the same side of the clip boundary, even though which endpoint intersects is different,
    // and the coordinates in that axis should *exactly* match the clip boundary's.
    clipper.init(new S2Cell(S2CellId.fromToken("81c")));
    S2Cell cell = new S2Cell(S2CellId.fromToken("7f4"));

    S2Point a = cell.getVertex(1);
    S2Point b = cell.getVertex(2);
    S2Point c = cell.getVertex(0);
    assertTrue(clipper.clipEdge(a, b));
    assertEquals(BOTTOM, clipper.outcode(0));
    assertExactly(-1.0, clipper.clippedUvEdge().v0.y());
    assertExactly(-1.0, clipper.clippedUvEdge().v1.y());

    assertTrue(clipper.clipEdge(b, c));
    assertEquals(BOTTOM, clipper.outcode(1));
    assertExactly(-1.0, clipper.clippedUvEdge().v0.y());
    assertExactly(-1.0, clipper.clippedUvEdge().v1.y());
  }

  private static S2Point point(double lat, double lng) {
    return S2LatLng.fromDegrees(lat, lng).toPoint();
  }

  private static R2Edge edge(R2Vector v0, R2Vector v1) {
    R2Edge e = new R2Edge();
    e.init(v0, v1);
    return e;
  }

  /**
   * Asserts that the given points are collinear. This isn't a robust geometric predicate so we
   * allow some error in the determinant.
   */
  private static void checkColinear(R2Vector a, R2Vector b, R2Vector c) {
    R2Vector ab = b.sub(a);
    R2Vector ac = c.sub(a);
    double det = ab.x() * ac.y() - ab.y() * ac.x();

    // The interpolation to an edge has a max absolute error of 2.25 epsilon, and the determinant
    // above 8 epsilon, which is 2.27e-15, round that up a bit to get a rough error bound to allow.
    assertEquals(0, det, 2.5e-15);
  }

  /**
   * Returns true iff 'a' is on the boundary of 'rect'. That is, one of its coordinates equals one
   * of the boundary coordinates.
   */
  private static boolean onBoundary(R2Vector a, R2Rect rect) {
    return (a.x() == rect.lo().x()) || (a.x() == rect.hi().x())
        || (a.y() == rect.lo().y()) || (a.y() == rect.hi().y());
  }

  /**
   * Builds 9 combinations of points in each sector of the 9x9 grid that Cohen-Sutherland operates
   * on:
   * <pre>
   *    L       R
   * T    |   |
   *   ---|---|---
   *      |   |
   *   ---|---|---
   * B    |   |
   * </pre>
   *
   * And then builds all 81 combinations of edges between those points.
   */
  private static List<EdgeClipperTestCase> buildTestCases() {
    TestPoint1D[] latTestPoints = {
        new TestPoint1D(TOP, 40.715, 'T'),
        new TestPoint1D(INSIDE, 40.700, '_'),
        new TestPoint1D(BOTTOM, 40.685, 'B')};
    TestPoint1D[] lonTestPoints = {
        new TestPoint1D(LEFT, -73.940, 'L'),
        new TestPoint1D(INSIDE, -73.925, '_'),
        new TestPoint1D(RIGHT, -73.910, 'R')};

    // Build all 9 combinations of coordinates to form 2D test points
    List<TestPoint2D> points = new ArrayList<>();
    for (int i = 0; i < 3; i++) {
      TestPoint1D latPnt = latTestPoints[i];
      for (int j = 0; j < 3; j++) {
        TestPoint1D lonPnt = lonTestPoints[j];
        points.add(new TestPoint2D(
            (byte) (latPnt.code | lonPnt.code),
            S2LatLng.fromDegrees(latPnt.val, lonPnt.val).toPoint()));
      }
    }

    // Build all 81 combinations of endpoints.
    List<EdgeClipperTestCase> testCases = new ArrayList<>();
    for (int i = 0; i < 9; i++) {
      TestPoint2D v0 = points.get(i);
      for (int j = 0; j < 9; ++j) {
        TestPoint2D v1 = points.get(j);
        testCases.add(new EdgeClipperTestCase(v0.point, v0.code, v1.point, v1.code));
      }
    }

    return testCases;
  }

  private static final class TestPoint1D {
    public final byte code;
    public final double val;
    public final char abbrev;
    private TestPoint1D(byte code, double val, char abbrev) {
      this.code = code;
      this.val = val;
      this.abbrev = abbrev;
    }
    @Override public String toString() {
      return Platform.formatString("TestPoint1D(%d, %f, %s)", code, val, abbrev);
    }
  }

  private static final class TestPoint2D {
    public final byte code; // The combined codes of the two endpoints.
    public final S2Point point;
    TestPoint2D(byte code, S2Point point) {
      this.code = code;
      this.point = point;
    }
  }

  private static final class EdgeClipperTestCase {
    public final S2Point v0;
    public final S2Point v1;

    // Outcodes for each vertex.
    public final byte outcode0;
    public final byte outcode1;

    public EdgeClipperTestCase(S2Point v0, byte code0, S2Point v1, byte code1) {
      this.v0 = v0;
      this.outcode0 = code0;
      this.v1 = v1;
      this.outcode1 = code1;
    }
  }
}
