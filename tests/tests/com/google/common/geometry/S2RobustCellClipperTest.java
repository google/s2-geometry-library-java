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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2Cell.Boundary.BOTTOM_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.LEFT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.RIGHT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.TOP_EDGE;
import static com.google.common.geometry.S2EdgeUtil.interpolate;
import static com.google.common.geometry.S2RobustCellClipper.CrossingType.INCOMING;
import static com.google.common.geometry.S2RobustCellClipper.CrossingType.OUTGOING;
import static com.google.common.geometry.S2RobustCellClipper.RobustClipResult.HIT_BOTH;
import static com.google.common.geometry.S2RobustCellClipper.RobustClipResult.HIT_NONE;
import static com.google.common.geometry.S2RobustCellClipper.RobustClipResult.HIT_V0;
import static com.google.common.geometry.S2RobustCellClipper.RobustClipResult.HIT_V1;
import static com.google.common.geometry.S2RobustCellClipper.RobustClipResult.MISS;
import static java.lang.Math.abs;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Cell.Boundary;
import com.google.common.geometry.S2RobustCellClipper.Crossing;
import com.google.common.geometry.S2RobustCellClipper.CrossingType;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.primitives.PooledList;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for the S2RobustCellClipper. */
@RunWith(JUnit4.class)
public final class S2RobustCellClipperTest extends GeometryTestCase {
  private static MutableEdge e(S2Point p0, S2Point p1) {
    return MutableEdge.of(p0, p1);
  }

  private static S2Cell c(S2CellId id) {
    return new S2Cell(id);
  }

  private static S2Point p(double x, double y, double z) {
    return new S2Point(x, y, z);
  }

  private static Crossing crossing(
      Boundary boundary, CrossingType type, double coord, double intercept, int edgeIndex) {
    Crossing c = new Crossing();
    c.set(boundary, type, coord, intercept, edgeIndex);
    return c;
  }

  private static Crossing crossing(
      Boundary boundary, CrossingType type, int coord, int intercept, int edgeIndex) {
    return crossing(boundary, type, (double) coord, (double) intercept, edgeIndex);
  }

  @Test
  public void testInteriorEdges() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));
    S2Point p0 = cell.getCenter();

    S2Point c0 = cell.getVertex(0);
    S2Point c1 = cell.getVertex(1);
    S2Point c2 = cell.getVertex(2);
    S2Point c3 = cell.getVertex(3);

    // Points that are a tiny bit inside the cell boundaries from the cell corners, enough that
    // they are properly contained in the cell but not far enough to bypass exact edge checks.
    S2Point[] v = {
      p0.sub(c0).normalize().mul(1e-30).add(c0).normalize(),
      p0.sub(c1).normalize().mul(1e-30).add(c1).normalize(),
      p0.sub(c2).normalize().mul(1e-30).add(c2).normalize(),
      p0.sub(c3).normalize().mul(1e-30).add(c3).normalize()
    };

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);

    // Check four edges between the points just inside the cell corners.

    assertEquals(HIT_BOTH, clipper.clipEdge(v[0], v[1]));
    assertTrue(clipper.getCrossings().isEmpty());
    assertEquals(HIT_BOTH, clipper.clipEdge(v[1], v[2]));
    assertTrue(clipper.getCrossings().isEmpty());
    assertEquals(HIT_BOTH, clipper.clipEdge(v[2], v[3]));
    assertTrue(clipper.getCrossings().isEmpty());
    assertEquals(HIT_BOTH, clipper.clipEdge(v[3], v[0]));
    assertTrue(clipper.getCrossings().isEmpty());

    // Check many tiny edges that are just inside the cell corners.
    for (int i = 0; i < 100; ++i) {
      S2Point e0 = v[data.uniform(3)];
      S2Point perturbation = data.getRandomPoint().mul(1e-31);
      S2Point e1 = e0.add(perturbation).normalize();
      assertEquals(HIT_BOTH, clipper.clipEdge(e0, e1));
      assertTrue(clipper.getCrossings().isEmpty());
    }
  }

  @Test
  public void testFaceMissDetected() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));

    // Two points on face 4, this edge misses face 0 entirely, and we should detect that.
    S2Point pnt0 = S2LatLng.fromDegrees(40.6714, -73.9181).toPoint();
    S2Point pnt1 = S2LatLng.fromDegrees(40.6344, -73.9737).toPoint();

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertEquals(MISS, clipper.clipEdge(pnt0, pnt1));
  }

  @Test
  public void testCornerToCorner() {
    // This will trigger the exact crossing path. We form an edge that crosses between the four
    // cells exactly at their shared corner.
    S2Cell cell0 = new S2Cell(S2CellId.fromToken("05"));
    S2Cell cell1 = new S2Cell(S2CellId.fromToken("1b"));
    S2Cell cell2 = new S2Cell(S2CellId.fromToken("11"));
    S2Cell cell3 = new S2Cell(S2CellId.fromToken("0f"));
    MutableEdge edge = e(cell0.getVertex(0), cell2.getVertex(2));

    {
      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(cell0);
      assertEquals(HIT_V0, clipper.clipEdge(edge));

      PooledList<Crossing> crossings = clipper.getCrossings();
      assertEquals(2, crossings.size());
      assertTrue(crossings.get(0).isEqualTo(crossing(RIGHT_EDGE, OUTGOING, 0, 0, 0)));
      assertTrue(crossings.get(1).isEqualTo(crossing(TOP_EDGE, OUTGOING, 0, 0, 0)));
    }

    {
      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(cell1);
      assertEquals(HIT_NONE, clipper.clipEdge(edge));

      PooledList<Crossing> crossings = clipper.getCrossings();
      assertEquals(2, crossings.size());
      assertTrue(crossings.get(0).isEqualTo(crossing(TOP_EDGE, OUTGOING, 0, 0, 0)));
      assertTrue(crossings.get(1).isEqualTo(crossing(LEFT_EDGE, INCOMING, 0, 0, 0)));
    }

    {
      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(cell2);
      assertEquals(HIT_V1, clipper.clipEdge(edge));

      PooledList<Crossing> crossings = clipper.getCrossings();
      assertEquals(2, crossings.size());
      assertTrue(crossings.get(0).isEqualTo(crossing(BOTTOM_EDGE, INCOMING, 0, 0, 0)));
      assertTrue(crossings.get(1).isEqualTo(crossing(LEFT_EDGE, INCOMING, 0, 0, 0)));
    }

    {
      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(cell3);
      assertEquals(HIT_NONE, clipper.clipEdge(edge));

      PooledList<Crossing> crossings = clipper.getCrossings();
      assertEquals(2, crossings.size());
      assertTrue(crossings.get(0).isEqualTo(crossing(BOTTOM_EDGE, INCOMING, 0, 0, 0)));
      assertTrue(crossings.get(1).isEqualTo(crossing(RIGHT_EDGE, OUTGOING, 0, 0, 0)));
    }
  }

  @Test
  public void testCornerGrazingDetected0() {
    // Cell with vertex 0 at (1,0,0.)
    S2Cell cell = new S2Cell(S2CellId.fromToken("14"));
    S2Point corner = cell.getVertex(0);
    assertTrue(corner.equalsPoint(new S2Point(1, 0, 0)));

    // Max error converting to UV and clipping to a face is ~1.5e-15, so we'll move the vertices out
    // a bit more than that so they don't fall within the error margin of the cell. First we'll move
    // them so that the edge does not cross the cell (passes corner at a distance of ~8e-14)
    double kTiny = 2e-15;
    S2Point pnt0 = corner.add(new S2Point(-kTiny, -2 * kTiny, +kTiny)).normalize();
    S2Point pnt1 = corner.add(new S2Point(-kTiny, +kTiny, -2 * kTiny)).normalize();

    // Now we'll displace them so that the edge does cross the cell.
    S2Point pnt2 = corner.add(new S2Point(-kTiny, -kTiny, +kTiny)).normalize();
    S2Point pnt3 = corner.add(new S2Point(-kTiny, +kTiny, -kTiny)).normalize();

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertEquals(MISS, clipper.clipEdge(pnt0, pnt1));
    assertEquals(HIT_NONE, clipper.clipEdge(pnt2, pnt3));

    PooledList<Crossing> crossings = clipper.getCrossings();
    assertEquals(2, crossings.size());
    assertTrue(crossings.get(0).isEqualTo(crossing(BOTTOM_EDGE, OUTGOING, 0, 0, 0)));
    assertTrue(crossings.get(1).isEqualTo(crossing(LEFT_EDGE, INCOMING, 0, 0, 0)));
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // to use exactly the same literals as C++.
  @Test
  public void testFalseMissDetected() {
    // Cell with vertex zero positive x, negative y, negative z.
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));

    // These two points form an edge that intersects the cell but that the UVEdgeClipper misses (a
    // false miss). We should detect this case and fall back to exact edge crossing tests which will
    // return the correct result.
    S2Point pnt0 = new S2Point(0.861698489631126274, -0.359041173122905954, -0.358559826207515031);
    S2Point pnt1 = new S2Point(0.858304494417438169, -0.357623947939909426, -0.367992536232324974);

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertEquals(HIT_NONE, clipper.clipEdge(pnt0, pnt1));

    PooledList<Crossing> crossings = clipper.getCrossings();
    double uv0 = cell.getUVCoordOfBoundary(BOTTOM_EDGE);
    double uv3 = cell.getUVCoordOfBoundary(LEFT_EDGE);

    assertEquals(2, crossings.size());
    assertEquals(BOTTOM_EDGE, crossings.get(0).boundary);
    assertEquals(OUTGOING, crossings.get(0).crossingType);
    assertAlmostEquals(crossings.get(0).intercept, uv0);
    assertAlmostEquals(crossings.get(0).coord, uv3);

    assertEquals(LEFT_EDGE, crossings.get(1).boundary);
    assertEquals(INCOMING, crossings.get(1).crossingType);
    assertDoubleNear(crossings.get(1).intercept, uv3, 1e-12);
    assertAlmostEquals(crossings.get(1).coord, uv0);
  }

  @Test
  public void testTrueHitDetected() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));
    MutableEdge edge = e(cell.getCenter(), cell.getCenter());

    // (degenerate) edge at the center of the cell, definitely contained with no surprises.
    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertEquals(HIT_BOTH, clipper.clipEdge(edge));
  }

  @Test
  public void testFalseHitDetected() {
    // Cell with vertex zero positive x, negative y, negative z.
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));

    // These two points form an edge that does not intersect the cell but that UVEdgeClipper thinks
    // does (a false hit). We should detect this case and fall back to exact edge crossing tests
    // which will return the correct result.
    S2Point pnt0 = new S2Point(0.955698120920362, 0.190765132372670, 0.224164595643751);
    S2Point pnt1 = new S2Point(0.957295555679071, 0.160089511834893, 0.240741702406469);

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertEquals(MISS, clipper.clipEdge(pnt0, pnt1));
    assertTrue(clipper.getCrossings().isEmpty());
  }

  /** Reflects a point across the plane defined by an edge. */
  S2Point reflectAcross(S2Point pnt, S2Point v0, S2Point v1) {
    S2Point normal = S2RobustCrossProd.robustCrossProd(v0, v1).normalize();
    return Matrix.householder(normal).mult(pnt);
  }

  // Check that multiple crossings on a boundary are ordered by intercept.
  @Test
  public void testCrossingsOrderedByIntercept() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));

    S2Point center = cell.getCenter();
    S2Point v0 = cell.getVertex(0);
    S2Point v1 = cell.getVertex(1);
    S2Point v2 = cell.getVertex(2);
    S2Point v3 = cell.getVertex(3);

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);

    // Reflect bottom two vertices across the top and clip edges.
    assertEquals(HIT_V0, clipper.clipEdge(center, reflectAcross(v0, v2, v3)));
    assertEquals(HIT_V0, clipper.clipEdge(center, reflectAcross(v1, v2, v3)));

    PooledList<Crossing> crossings = clipper.getCrossings();
    assertEquals(2, crossings.size());
    assertGreaterOrEqual(crossings.get(0).intercept, crossings.get(1).intercept);
  }

  // Check that disabling crossings produces no results.
  @Test
  public void testNoCrossingsWorks() {
    S2Cell cell = new S2Cell(S2CellId.fromToken("05"));
    S2Point center = cell.getCenter();
    S2Point v0 = cell.getVertex(0);
    S2Point v1 = cell.getVertex(1);
    S2Point v2 = cell.getVertex(2);
    S2Point v3 = cell.getVertex(3);

    S2RobustCellClipper.Options options = new S2RobustCellClipper.Options();
    options.setEnableCrossings(false);

    S2RobustCellClipper clipper = new S2RobustCellClipper(options);
    clipper.startCell(cell);

    // Reflect bottom two vertices across the top and clip edges.
    assertEquals(HIT_V0, clipper.clipEdge(center, reflectAcross(v0, v2, v3)));
    assertEquals(HIT_V0, clipper.clipEdge(center, reflectAcross(v1, v2, v3)));

    PooledList<Crossing> crossings = clipper.getCrossings();
    assertTrue(crossings.isEmpty());
  }

  @Test
  public void testCoplanarEdges() {
    // Parent/child cells that have consecutive coplanar boundaries. The first pair are coplanar on
    // edge 0, the second on edge 1, etc.
    ImmutableList<Pair<S2Cell, S2Cell>> kCells =
        ImmutableList.of(
            Pair.of(c(S2CellId.fromToken("11")), c(S2CellId.fromToken("107"))),
            Pair.of(c(S2CellId.fromToken("0f")), c(S2CellId.fromToken("0e3"))),
            Pair.of(c(S2CellId.fromToken("05")), c(S2CellId.fromToken("057"))),
            Pair.of(c(S2CellId.fromToken("1b")), c(S2CellId.fromToken("1ad"))));

    S2RobustCellClipper clipper = new S2RobustCellClipper();

    // For each pair (a,b) of cells, we'll test each combination of swapping the cells and swapping
    // the edge's vertices.
    for (int rep = 0; rep < 4; ++rep) {
      for (int i = 0; i < 4; ++i) {
        boolean swap = (rep & 1) > 0;
        boolean flip = (rep & 2) > 0;

        S2Cell cell0 = kCells.get(i).first;
        S2Cell cell1 = kCells.get(i).second;

        if (swap) {
          // swap(cell0, cell1): Clip the larger cell to the smaller instead.
          S2Cell tmp = cell0;
          cell0 = cell1;
          cell1 = tmp;
        }

        // Clip the edge.
        clipper.startCell(cell0);
        S2Point v0 = cell1.getVertex(i + 0);
        S2Point v1 = cell1.getVertex(i + 1);
        if (flip) {
          S2Point tmp = v0;
          v0 = v1;
          v1 = tmp;
        }

        // The perturbed sign is based on the cell edge normal's first non-zero component, figure
        // out what that is so we know if the edge _should_
        // test as crossing the boundary or not.
        S2Point normal = cell1.getEdgeRaw(i);
        int sign = 0;
        for (int j = 0; j < 3; ++j) {
          if (normal.get(j) != 0) {
            sign = (normal.get(j) > 0) ? +1 : -1;
            break;
          }
        }

        // When clipping the smaller cell to the larger, the edge will be entirely on one side or
        // the other, with no crossings.
        //
        // When we swap and clip the larger cell to the smaller, the edge extends past the other two
        // boundaries so we should get a crossing on each of
        // them when the sign is positive.
        int expectedCrossings;
        if (swap) {
          if (sign > 0) {
            expectedCrossings = 2;
            assertEquals(HIT_NONE, clipper.clipEdge(v0, v1));
          } else {
            expectedCrossings = 0;
            assertEquals(MISS, clipper.clipEdge(v0, v1));
          }
        } else {
          expectedCrossings = 0;

          // The perturbation should make us miss on boundaries 1 and 2.
          if (i == 1 || i == 2) {
            assertEquals(MISS, clipper.clipEdge(v0, v1));
          } else {
            assertEquals(HIT_BOTH, clipper.clipEdge(v0, v1));
          }
        }
        assertEquals(expectedCrossings, clipper.getCrossings().size());

        S2Cell smaller = cell0.id().lowestOnBit() >= cell1.id().lowestOnBit() ? cell1 : cell0;
        double uv0 = smaller.getBoundUV().getVertex(i + 0).get(i % 2);
        double uv1 = smaller.getBoundUV().getVertex(i + 1).get(i % 2);

        CrossingType crossingType0 = INCOMING;
        CrossingType crossingType1 = OUTGOING;
        if (flip) {
          // swap(type0, type1);
          CrossingType crossingTypeTmp = crossingType0;
          crossingType0 = crossingType1;
          crossingType1 = crossingTypeTmp;

          // swap(uv0, uv1);
          double uvTmp = uv0;
          uv0 = uv1;
          uv1 = uvTmp;
        }

        double coordCurr = smaller.getUVCoordOfBoundary(Boundary.fromValue(i));
        double coordNext = smaller.getUVCoordOfBoundary(Boundary.fromValue(i + 1));
        double coordPrev = smaller.getUVCoordOfBoundary(Boundary.fromValue(i + 3));

        List<Crossing> crossings = new ArrayList<>();
        if (swap && sign > 0) {
          Boundary b1 = Boundary.fromValue(i + 1);
          crossings.add(crossing(b1, crossingType1, coordNext, coordCurr, 0));
          Boundary b2 = Boundary.fromValue(i + 3);
          crossings.add(crossing(b2, crossingType0, coordPrev, coordCurr, 0));
        }

        // EXPECT_THAT(clipper.getCrossings(), UnorderedPointwise(Eq(), crossings));
        assertCrossingsEqual(crossings, clipper.getCrossings());
      }
    }
  }

  /**
   * Asserts that the "expected" List of Crossings has the same Crossings as the "actual" PooledList
   * of Crossings. TODO(torrey): Improve this so that order does not matter. For now, the tests pass
   * anyway.
   */
  private static void assertCrossingsEqual(List<Crossing> expected, PooledList<Crossing> actual) {
    assertEquals(expected.size(), actual.size());
    for (int i = 0; i < expected.size(); ++i) {
      assertTrue(expected.get(i).isEqualTo(actual.get(i)));
    }
  }

  @Test
  public void testCoplanarExterior() {
    // Two cells who's edge 0 are parallel but do not intersect.
    S2Cell kCell0 = new S2Cell(S2CellId.fromToken("107"));
    S2Cell kCell1 = new S2Cell(S2CellId.fromToken("10b"));

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(kCell0);

    S2Point v0 = kCell1.getVertex(0);
    S2Point v1 = kCell1.getVertex(1);
    assertEquals(MISS, clipper.clipEdge(v0, v1));
  }

  @Test
  public void testCoplanarStraddling() {
    // A cell on the equator, and two cells that straddle each end of edge 0.
    S2Cell kCell0 = new S2Cell(S2CellId.fromToken("104"));
    S2Cell[][] kCells = {
      {c(S2CellId.fromToken("0ff")), c(S2CellId.fromToken("101"))},
      {c(S2CellId.fromToken("107")), c(S2CellId.fromToken("109"))}
    };

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    for (int i = 0; i < 2; ++i) {
      clipper.startCell(kCell0);
      S2Point v0 = kCells[i][0].getVertex(0);
      S2Point v1 = kCells[i][1].getVertex(1);
      assertTrue(clipper.clipEdge(v0, v1).hit());

      double coord = kCell0.getUVCoordOfBoundary(BOTTOM_EDGE);
      PooledList<Crossing> crossings = clipper.getCrossings();
      assertEquals(1, crossings.size());

      // We should get a crossing on the boundary we extend past.
      Boundary b = LEFT_EDGE;
      CrossingType crossingType = INCOMING;
      if (i == 1) {
        b = RIGHT_EDGE;
        crossingType = OUTGOING;
      }
      assertTrue(
          crossings
              .get(0)
              .isEqualTo(crossing(b, crossingType, kCell0.getUVCoordOfBoundary(b), coord, 0)));
    }
  }

  @Test
  public void testClipCellToSelf() {
    S2Cell kCell0 = new S2Cell(S2CellId.fromToken("1"));

    // We should get a hit for each edge but they shouldn't add a crossing.
    S2RobustCellClipper clipper = new S2RobustCellClipper();
    for (boolean flip : ImmutableList.of(false, true)) {
      for (int b = 0; b < 4; ++b) {
        clipper.startCell(kCell0);

        S2Point v0 = kCell0.getVertex(b + 0);
        S2Point v1 = kCell0.getVertex(b + 1);
        if (flip) {
          // swap(v0, v1);
          S2Point tmp = v0;
          v0 = v1;
          v1 = tmp;
        }
        assertTrue(clipper.clipEdge(v0, v1).hit());
        assertTrue(clipper.getCrossings().isEmpty());
      }
    }
  }

  @Test
  public void testExactEquatorPointDoesNotCross() {
    // An edge with a vertex that lands exactly on the equator, which is the boundary between two
    // cells. Shouldn't count as a crossing.
    MutableEdge kEdge =
        MutableEdge.of(
            S2LatLng.fromDegrees(0.00048, 120.032482).toPoint(),
            S2LatLng.fromDegrees(0.00000, 120.032743).toPoint());

    S2Cell kCell0 = new S2Cell(S2CellId.fromToken("3275f89d"));
    S2Cell kCell1 = new S2Cell(S2CellId.fromToken("2d8a077"));

    // Test the edge and its reverse.
    for (MutableEdge edge : ImmutableList.of(kEdge, MutableEdge.of(kEdge.b, kEdge.a))) {
      S2RobustCellClipper clipper = new S2RobustCellClipper();

      clipper.startCell(kCell0);
      assertEquals(HIT_BOTH, clipper.clipEdge(edge));
      assertTrue(clipper.getCrossings().isEmpty());

      clipper.startCell(kCell1);
      assertEquals(MISS, clipper.clipEdge(edge));
      assertTrue(clipper.getCrossings().isEmpty());
    }
  }

  @Test
  public void testCloseCrossingsOrderedCorrectly0() {
    // A central cell and four cells around it to use as reference points. We order them by which
    // boundary they're adjacent to, which also lets us enumerate the two corners further from the
    // center point easily:
    //
    //      3   2
    //      ┌───┐
    //      │ 2 │
    // 3┌───┼───┼───┐2
    //  │ 3 │ • │ 1 │
    // 0└───┼───┼───┘1
    //      │ 0 │
    //      └───┘
    //      0   1
    S2Cell kCell = new S2Cell(S2CellId.fromToken("1b"));
    S2Cell[] kCellNeighbor = {
      c(S2CellId.fromToken("1d")), c(S2CellId.fromToken("19")),
      c(S2CellId.fromToken("11")), c(S2CellId.fromToken("05"))
    };

    // Add three evenly spaced crossings on each boundary. One through the middle drawing to the
    // adjacent cell center, and two drawing to the corners of the adjacent cell. This gives us 12
    // crossings total, three per boundary:
    //
    //   [A0, A1, A2, B0, B1, B2, C0, C1, C2, D0, D1, D2]
    //
    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(kCell);

    S2Point v0 = kCell.getCenter();
    for (int k = 0; k < 4; ++k) {
      clipper.clipEdge(v0, kCellNeighbor[k].getVertex(k + 0));
      clipper.clipEdge(v0, kCellNeighbor[k].getCenter());
      clipper.clipEdge(v0, kCellNeighbor[k].getVertex(k + 1));
    }
    assertEquals(12, clipper.getCrossings().size());

    // Now add another crossing near the center, but perturb the vertices a tiny
    // bit so that it crosses very close to the other center crossing.
    S2Point yinc = new S2Point(0, DBL_EPSILON, 0);
    S2Point zinc = new S2Point(0, 0, DBL_EPSILON);

    clipper.clipEdge(v0.add(yinc), kCellNeighbor[0].getCenter().add(yinc));
    clipper.clipEdge(v0.add(zinc), kCellNeighbor[1].getCenter().add(zinc));
    clipper.clipEdge(v0.sub(yinc), kCellNeighbor[2].getCenter().sub(yinc));
    clipper.clipEdge(v0.sub(zinc), kCellNeighbor[3].getCenter().sub(zinc));

    PooledList<Crossing> crossings = clipper.getCrossings();
    assertEquals(16, crossings.size());

    // Check that each edge did indeed cross very close to its neighbor and that
    // we ordered them correctly. We should have an array of 16 crossings now,
    // four per boundary, and the perturbed crossing should come right after the
    // middle crossing of each:
    //
    // k           0               1               2               3
    // i      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    //      [A0, A1, AP, A2, B0, B1, BP, B2, C0, C1, CP, C2, D0, D1, DP, D2]
    // edge   0   1  12   2   3   4  13   5   6   7  14   8   9  10  15  11
    //
    // The tolerance here is just determined empirically.
    for (int k = 0; k < 4; ++k) {
      assertLessOrEqual(
          abs(crossings.get(4 * k + 1).intercept - crossings.get(4 * k + 2).intercept), 2.5e-16);

      assertEquals(3 * k + 1, crossings.get(4 * k + 1).edgeIndex);
      assertEquals(k + 12, crossings.get(4 * k + 2).edgeIndex);
    }
  }

  // Same as above but we only cross once per boundary before perturbing.
  @Test
  public void testCloseCrossingsOrderedCorrectly1() {
    S2Cell kCell = new S2Cell(S2CellId.fromToken("1b"));

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(kCell);

    S2Point v0 = kCell.getCenter();
    for (int k = 0; k < 4; ++k) {
      S2Point v1 = reflectAcross(v0, kCell.getVertex(k), kCell.getVertex(k + 1));
      clipper.clipEdge(v0, v1);
    }
    assertEquals(4, clipper.getCrossings().size());

    // Now add another crossing near the center, but perturb the vertices a tiny
    // bit so that it crosses very close to the other center crossing.
    S2Point yinc = new S2Point(0, DBL_EPSILON, 0);
    S2Point zinc = new S2Point(0, 0, DBL_EPSILON);

    S2Point v1;
    v1 = reflectAcross(v0, kCell.getVertex(0), kCell.getVertex(1));
    clipper.clipEdge(v0.add(yinc), v1.add(yinc));

    v1 = reflectAcross(v0, kCell.getVertex(1), kCell.getVertex(2));
    clipper.clipEdge(v0.add(zinc), v1.add(zinc));

    v1 = reflectAcross(v0, kCell.getVertex(2), kCell.getVertex(3));
    clipper.clipEdge(v0.sub(yinc), v1.sub(yinc));

    v1 = reflectAcross(v0, kCell.getVertex(3), kCell.getVertex(0));
    clipper.clipEdge(v0.sub(zinc), v1.sub(zinc));

    PooledList<Crossing> crossings = clipper.getCrossings();
    assertEquals(8, crossings.size());
    //
    // The tolerance here is just determined empirically.
    for (int k = 0; k < 4; ++k) {
      assertLessOrEqual(
          abs(crossings.get(2 * k).intercept - crossings.get(2 * k + 1).intercept), 2.5e-16);

      assertEquals(k, crossings.get(2 * k + 0).edgeIndex);
      assertEquals(k + 4, crossings.get(2 * k + 1).edgeIndex);
    }
  }

  // Same as above but we perturb the opposite direction to switch the order.
  @Test
  public void testCloseCrossingsOrderedCorrectly2() {
    // Same as above but we only cross once per boundary before perturbing.
    S2Cell kCell = new S2Cell(S2CellId.fromToken("1b"));

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(kCell);

    S2Point v0 = kCell.getCenter();
    for (int k = 0; k < 4; ++k) {
      S2Point v1 = reflectAcross(v0, kCell.getVertex(k), kCell.getVertex(k + 1));
      clipper.clipEdge(v0, v1);
    }
    assertEquals(4, clipper.getCrossings().size());

    // Now add another crossing near the center, but perturb the vertices a tiny
    // bit so that it crosses very close to the other center crossing.
    S2Point yinc = new S2Point(0, DBL_EPSILON, 0);
    S2Point zinc = new S2Point(0, 0, DBL_EPSILON);

    S2Point v1;
    v1 = reflectAcross(v0, kCell.getVertex(0), kCell.getVertex(1));
    clipper.clipEdge(v0.sub(yinc), v1.sub(yinc));

    v1 = reflectAcross(v0, kCell.getVertex(1), kCell.getVertex(2));
    clipper.clipEdge(v0.sub(zinc), v1.sub(zinc));

    v1 = reflectAcross(v0, kCell.getVertex(2), kCell.getVertex(3));
    clipper.clipEdge(v0.add(yinc), v1.add(yinc));

    v1 = reflectAcross(v0, kCell.getVertex(3), kCell.getVertex(0));
    clipper.clipEdge(v0.add(zinc), v1.add(zinc));

    PooledList<Crossing> crossings = clipper.getCrossings();
    assertEquals(8, crossings.size());
    //
    // The tolerance here is just determined empirically.
    for (int k = 0; k < 4; ++k) {
      assertLessOrEqual(
          abs(crossings.get(2 * k).intercept - crossings.get(2 * k + 1).intercept), 2.5e-16);

      assertEquals(k, crossings.get(2 * k + 1).edgeIndex);
      assertEquals(k + 4, crossings.get(2 * k + 0).edgeIndex);
    }
  }

  @Test
  public void testDuplicateCrossingsCancel() {
    S2Cell kCell = new S2Cell(S2CellId.fromToken("1b"));

    for (int i = 0; i < 2; ++i) {
      boolean swap = (i == 1);

      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(kCell);

      // Cross each boundary of a cell in the middle by reflecting its center over the boundary, we
      // should get a crossing for each one.
      S2Point v0 = kCell.getCenter();
      for (int k = 0; k < 4; ++k) {
        S2Point v1 = reflectAcross(v0, kCell.getVertex(k), kCell.getVertex(k + 1));
        clipper.clipEdge(v0, v1);
        assertEquals(k + 1, clipper.getCrossings().size());
      }
      assertEquals(4, clipper.getCrossings().size());

      // Now cross the boundaries again at the same exact point with either an exact duplicate or
      // exact reverse reversed edge. Either way, the duplicate crossings will cancel out and the
      // number of crossings will be reduced by one each iteration.
      for (int k = 0; k < 4; ++k) {
        S2Point v1 = reflectAcross(v0, kCell.getVertex(k), kCell.getVertex(k + 1));
        if (swap) {
          clipper.clipEdge(v1, v0);
        } else {
          clipper.clipEdge(v0, v1);
        }
        assertEquals(3 - k, clipper.getCrossings().size());
      }
    }
  }

  @Test
  public void testFlippedCrossingsCorrectOrder() {

    // Cell who's right and top boundaries are on longitude 0 and its neighbor to
    // the right.
    //
    // 3┌───┬───┐2
    //  │ L │ R │
    // 0└───┴───┘1
    //
    S2Cell kCellL = c(S2CellId.fromToken("05"));
    S2Cell kCellR = c(S2CellId.fromToken("1b"));

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(kCellR);

    // Move corners of the cells toward their respective centers to avoid crossing the far
    // boundaries of the cells as well. We only want to cross the shared boundary segment.
    S2Point v0 = interpolate(kCellL.getVertex(0), kCellL.getCenter(), 1e-15);
    S2Point v1 = interpolate(kCellR.getVertex(1), kCellR.getCenter(), 1e-16);
    S2Point v2 = interpolate(kCellR.getVertex(2), kCellR.getCenter(), 1e-16);
    S2Point v3 = interpolate(kCellL.getVertex(3), kCellL.getCenter(), 1e-15);

    // Displacements for vertices 0 and 3 statically fuzzed to cause the boundary crossings to
    // switch order due to numerical error.
    S2Point[] kDelta =
        new S2Point[] {
          p(8.054597e-17, -4.491415e-17, 3.866609e-17), p(7.969033e-17, -4.919605e-17, 3.505995e-17)
        };

    for (int i = 0; i < kDelta.length; ++i) {
      clipper.reset();

      // Cross (nearly) in the middle of the line segment with some displacement that causes the
      // true crossing order to swap places.
      clipper.clipEdge(v1, v3.add(kDelta[i]).normalize());
      clipper.clipEdge(v2, v0.sub(kDelta[i]).normalize());

      // We crossed boundary 3, so our intercepts should be ordered in descending order except that
      // numerical error caused two crossings to swap , which the exact predicate logic should have
      // put back into the correct place.
      //
      // So check that the first crossing has a smaller intercept value than the second crossing,
      // despite coming before it in the crossing list.
      PooledList<Crossing> crossings = clipper.getCrossings();
      assertEquals(2, crossings.size());
      assertEquals(LEFT_EDGE, crossings.get(0).boundary);
      assertEquals(LEFT_EDGE, crossings.get(1).boundary);
      assertLessThan(crossings.get(0).intercept, crossings.get(1).intercept);
    }
  }

  @Test
  public void testHorizontalAfterUVConversion() {
    // This is a regression encountered in production consisting of an edge that crosses the top
    // boundary of a cell when using exact tests, but, when converted to UV, becomes exactly
    // horizontal, resulting in an interpolation failure.
    //
    // We compute the intercept point using cross products now instead so this should work.
    S2Cell kCell = c(S2CellId.fromToken("1284"));

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(kCell);

    S2Point[] kVertices = {
      new S2Point(0x1.96bf38faa05bp-1, 0x1.a9d02fa65fdf4p-5, 0x1.35d3a866e8255p-1),
      new S2Point(0x1.970f50cbe3162p-1, 0x1.17da878c2c1f3p-5, 0x1.3610aa8b4df9ep-1)
    };

    // The edge should miss but not throw an assertion.
    assertFalse(clipper.clipEdge(kVertices[0], kVertices[1]).hit());
  }

  public void checkCrossingCellBoundaryWorks(int boundaryEdge) {
    // Draw an edge from the center of the cell to the center's image across each edge boundary.
    // This should be very well conditioned so we can check that basic intersection and crossing
    // reporting works.

    S2Cell cell = new S2Cell(S2CellId.fromToken("114"));
    S2Point center = cell.getCenter();
    S2Point v0 = cell.getVertex(boundaryEdge);
    S2Point v1 = cell.getVertex(boundaryEdge + 1);
    S2Point v2 = cell.getVertex(boundaryEdge + 2);
    S2Point v3 = cell.getVertex(boundaryEdge + 3);

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertTrue(clipper.clipEdge(center, reflectAcross(center, v0, v1)).hit());

    // Check that crossing point was close-ish to the midpoint of the edge.
    PooledList<Crossing> crossings = clipper.getCrossings();
    assertEquals(1, crossings.size());

    R2Rect uvbound = cell.getBoundUV();
    // Left and right boundaries (1, and 3) need to test against y component, top and bottom
    // boundaries (0 and 2) need to test against x component.
    double midpnt = ((boundaryEdge % 2 > 0) ? uvbound.y() : uvbound.x()).getCenter();
    assertDoubleNear(crossings.get(0).intercept, midpnt, 0.1);

    // Now reflect across the boundary and the opposite boundary. This should result in two crossing
    // entries.
    clipper.reset();
    assertEquals(
        HIT_NONE, clipper.clipEdge(reflectAcross(center, v0, v1), reflectAcross(center, v2, v3)));
    assertEquals(2, clipper.getCrossings().size());
  }

  @Test
  public void testCrossingCellBoundaryWorks() {
    checkCrossingCellBoundaryWorks(0 /* Bottom */);
    checkCrossingCellBoundaryWorks(1 /* Right */);
    checkCrossingCellBoundaryWorks(2 /* Top */);
    checkCrossingCellBoundaryWorks(3 /* Left */);
  }

  @Test
  public void testCornerGrazingBoundaryContainment() {
    // Cell with vertex 0 at (1,0,0.)
    S2Cell cell = new S2Cell(S2CellId.fromToken("14"));
    S2Point corner = cell.getVertex(0);
    assertTrue(corner.equalsPoint(new S2Point(1, 0, 0)));

    // Make an edge outside the cell that comes within ~8e-14 of the corner.
    double kTiny = 2e-15;
    S2Point pnt0 = corner.add(new S2Point(-kTiny, -2 * kTiny, +kTiny)).normalize();
    S2Point pnt1 = corner.add(new S2Point(-kTiny, +kTiny, -2 * kTiny)).normalize();

    S2RobustCellClipper clipper = new S2RobustCellClipper();
    clipper.startCell(cell);
    assertEquals(MISS, clipper.clipEdge(pnt0, pnt1));

    // The edge should just be ignored and we return the center containment.
    assertTrue(clipper.isBoundaryContained(true));
    assertFalse(clipper.isBoundaryContained(false));
  }

  @Test
  public void testRingAroundCenterFlipsBoundary() {
    S2Point v0 = S2LatLng.fromDegrees(-10, 0).toPoint();
    S2Point v1 = S2LatLng.fromDegrees(0, +10).toPoint();
    S2Point v2 = S2LatLng.fromDegrees(+10, 0).toPoint();
    S2Point v3 = S2LatLng.fromDegrees(0, -10).toPoint();

    {
      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(S2Cell.fromFace(0));
      assertEquals(HIT_BOTH, clipper.clipEdge(v0, v1));
      assertEquals(HIT_BOTH, clipper.clipEdge(v1, v2));
      assertEquals(HIT_BOTH, clipper.clipEdge(v2, v3));
      assertEquals(HIT_BOTH, clipper.clipEdge(v3, v0));

      // The ring should be entirely contained in the face so crossing edges should just flip the
      // center containment bit to get the boundary containment.
      assertFalse(clipper.isBoundaryContained(true));
      assertTrue(clipper.isBoundaryContained(false));
    }

    {
      // Flip the order of vertices to invert the polygon, the center containment bit should still
      // flip.
      S2RobustCellClipper clipper = new S2RobustCellClipper();
      clipper.startCell(S2Cell.fromFace(0));
      assertEquals(HIT_BOTH, clipper.clipEdge(v0, v3));
      assertEquals(HIT_BOTH, clipper.clipEdge(v3, v2));
      assertEquals(HIT_BOTH, clipper.clipEdge(v2, v1));
      assertEquals(HIT_BOTH, clipper.clipEdge(v1, v0));

      assertFalse(clipper.isBoundaryContained(true));
      assertTrue(clipper.isBoundaryContained(false));
    }
  }
}
