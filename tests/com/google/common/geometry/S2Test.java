/*
 * Copyright 2005 Google Inc.
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

import java.util.logging.Logger;

public strictfp class S2Test extends GeometryTestCase {

  private static Logger logger = Logger.getLogger(S2Test.class.getName());

  // Test helper methods for testing the traversal order.
  private static int swapAxes(int ij) {
    return ((ij >> 1) & 1) + ((ij & 1) << 1);
  }

  private static int invertBits(int ij) {
    return ij ^ 3;
  }

  public void testTraversalOrder() {
    for (int r = 0; r < 4; ++r) {
      for (int i = 0; i < 4; ++i) {
        // Check consistency with respect to swapping axes.
        assertEquals(S2.ijToPos(r, i),
            S2.ijToPos(r ^ S2.SWAP_MASK, swapAxes(i)));
        assertEquals(S2.posToIJ(r, i),
            swapAxes(S2.posToIJ(r ^ S2.SWAP_MASK, i)));

        // Check consistency with respect to reversing axis directions.
        assertEquals(S2.ijToPos(r, i),
            S2.ijToPos(r ^ S2.INVERT_MASK, invertBits(i)));
        assertEquals(S2.posToIJ(r, i),
            invertBits(S2.posToIJ(r ^ S2.INVERT_MASK, i)));

        // Check that the two tables are inverses of each other.
        assertEquals(S2.ijToPos(r, S2.posToIJ(r, i)), i);
        assertEquals(S2.posToIJ(r, S2.ijToPos(r, i)), i);
      }
    }
  }

  public void testSTUV() {
    // Check boundary conditions.
    for (double x = -1; x <= 1; ++x) {
      assertEquals(S2Projections.stToUV(x), x);
      assertEquals(S2Projections.uvToST(x), x);
    }
    // Check that UVtoST and STtoUV are inverses.
    for (double x = -1; x <= 1; x += 0.0001) {
      assertDoubleNear(S2Projections.uvToST(S2Projections.stToUV(x)), x);
      assertDoubleNear(S2Projections.stToUV(S2Projections.uvToST(x)), x);
    }
  }

  public void testFaceUVtoXYZ() {
    // Check that each face appears exactly once.
    S2Point sum = new S2Point();
    for (int face = 0; face < 6; ++face) {
      S2Point center = S2Projections.faceUvToXyz(face, 0, 0);
      assertEquals(S2Projections.getNorm(face), center);
      assertEquals(Math.abs(center.get(center.largestAbsComponent())), 1.0);
      sum = S2Point.add(sum, S2Point.fabs(center));
    }
    assertEquals(sum, new S2Point(2, 2, 2));

    // Check that each face has a right-handed coordinate system.
    for (int face = 0; face < 6; ++face) {
      assertEquals(
          S2Point.crossProd(S2Projections.getUAxis(face), S2Projections.getVAxis(face)).dotProd(
              S2Projections.faceUvToXyz(face, 0, 0)), 1.0);
    }

    // Check that the Hilbert curves on each face combine to form a
    // continuous curve over the entire cube.
    for (int face = 0; face < 6; ++face) {
      // The Hilbert curve on each face starts at (-1,-1) and terminates
      // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
      int sign = ((face & S2.SWAP_MASK) != 0) ? -1 : 1;
      assertEquals(S2Projections.faceUvToXyz(face, sign, -sign),
          S2Projections.faceUvToXyz((face + 1) % 6, -1, -1));
    }
  }

  public void testUVNorms() {
    // Check that GetUNorm and GetVNorm compute right-handed normals for
    // an edge in the increasing U or V direction.
    for (int face = 0; face < 6; ++face) {
      for (double x = -1; x <= 1; x += 1 / 1024.) {
        assertDoubleNear(
            S2Point.crossProd(
                S2Projections.faceUvToXyz(face, x, -1), S2Projections.faceUvToXyz(face, x, 1))
                .angle(S2Projections.getUNorm(face, x)), 0);
        assertDoubleNear(
            S2Point.crossProd(
                S2Projections.faceUvToXyz(face, -1, x), S2Projections.faceUvToXyz(face, 1, x))
                .angle(S2Projections.getVNorm(face, x)), 0);
      }
    }
  }

  public void testUVAxes() {
    // Check that axes are consistent with FaceUVtoXYZ.
    for (int face = 0; face < 6; ++face) {
      assertEquals(S2Projections.getUAxis(face), S2Point.sub(
          S2Projections.faceUvToXyz(face, 1, 0), S2Projections.faceUvToXyz(face, 0, 0)));
      assertEquals(S2Projections.getVAxis(face), S2Point.sub(
          S2Projections.faceUvToXyz(face, 0, 1), S2Projections.faceUvToXyz(face, 0, 0)));
    }
  }

  public void testAngleArea() {
    S2Point pz = new S2Point(0, 0, 1);
    S2Point p000 = new S2Point(1, 0, 0);
    S2Point p045 = new S2Point(1, 1, 0);
    S2Point p090 = new S2Point(0, 1, 0);
    S2Point p180 = new S2Point(-1, 0, 0);
    assertDoubleNear(S2.angle(p000, pz, p045), S2.M_PI_4);
    assertDoubleNear(S2.angle(p045, pz, p180), 3 * S2.M_PI_4);
    assertDoubleNear(S2.angle(p000, pz, p180), S2.M_PI);
    assertDoubleNear(S2.angle(pz, p000, pz), 0);
    assertDoubleNear(S2.angle(pz, p000, p045), S2.M_PI_2);

    assertDoubleNear(S2.area(p000, p090, pz), S2.M_PI_2);
    assertDoubleNear(S2.area(p045, pz, p180), 3 * S2.M_PI_4);

    // Make sure that area() has good *relative* accuracy even for
    // very small areas.
    final double eps = 1e-10;
    S2Point pepsx = new S2Point(eps, 0, 1);
    S2Point pepsy = new S2Point(0, eps, 1);
    double expected1 = 0.5 * eps * eps;
    assertDoubleNear(S2.area(pepsx, pepsy, pz), expected1, 1e-14 * expected1);

    // Make sure that it can handle degenerate triangles.
    S2Point pr = new S2Point(0.257, -0.5723, 0.112);
    S2Point pq = new S2Point(-0.747, 0.401, 0.2235);
    assertEquals(S2.area(pr, pr, pr), 0.0);
    // TODO: The following test is not exact in optimized mode because the
    // compiler chooses to mix 64-bit and 80-bit intermediate results.
    assertDoubleNear(S2.area(pr, pq, pr), 0);
    assertEquals(S2.area(p000, p045, p090), 0.0);

    double maxGirard = 0;
    for (int i = 0; i < 10000; ++i) {
      S2Point p0 = randomPoint();
      S2Point d1 = randomPoint();
      S2Point d2 = randomPoint();
      S2Point p1 = S2Point.add(p0, S2Point.mul(d1, 1e-15));
      S2Point p2 = S2Point.add(p0, S2Point.mul(d2, 1e-15));
      // The actual displacement can be as much as 1.2e-15 due to roundoff.
      // This yields a maximum triangle area of about 0.7e-30.
      assertTrue(S2.area(p0, p1, p2) < 0.7e-30);
      maxGirard = Math.max(maxGirard, S2.girardArea(p0, p1, p2));
    }
    logger.info("Worst case Girard for triangle area 1e-30: " + maxGirard);

    // Try a very long and skinny triangle.
    S2Point p045eps = new S2Point(1, 1, eps);
    double expected2 = 5.8578643762690495119753e-11; // Mathematica.
    assertDoubleNear(S2.area(p000, p045eps, p090), expected2, 1e-9 * expected2);

    // Triangles with near-180 degree edges that sum to a quarter-sphere.
    final double eps2 = 1e-10;
    S2Point p000eps2 = new S2Point(1, 0.1 * eps2, eps2);
    double quarterArea1 =
        S2.area(p000eps2, p000, p090) + S2.area(p000eps2, p090, p180) + S2.area(p000eps2, p180, pz)
            + S2.area(p000eps2, pz, p000);
    assertDoubleNear(quarterArea1, S2.M_PI);

    // Four other triangles that sum to a quarter-sphere.
    S2Point p045eps2 = new S2Point(1, 1, eps2);
    double quarterArea2 =
        S2.area(p045eps2, p000, p090) + S2.area(p045eps2, p090, p180) + S2.area(p045eps2, p180, pz)
            + S2.area(p045eps2, pz, p000);
    assertDoubleNear(quarterArea2, S2.M_PI);
  }

  public void testCCW() {
    S2Point a = new S2Point(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
    S2Point b = new S2Point(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
    S2Point c = new S2Point(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
    assertTrue(S2.robustCCW(a, b, c) != 0);
  }

  // Note: obviously, I could have defined a bundle of metrics like this in the
  // S2 class itself rather than just for testing. However, it's not clear that
  // this is useful other than for testing purposes, and I find
  // S2.kMinWidth.GetMaxLevel(width) to be slightly more readable than
  // than S2.kWidth.min().GetMaxLevel(width). Also, there is no fundamental
  // reason that we need to analyze the minimum, maximum, and average values of
  // every metric; it would be perfectly reasonable to just define one of these.

  class MetricBundle {
    public MetricBundle(S2.Metric min, S2.Metric max, S2.Metric avg) {
      this.min_ = min;
      this.max_ = max;
      this.avg_ = avg;
    }

    S2.Metric min_;
    S2.Metric max_;
    S2.Metric avg_;
  }

  public void testMinMaxAvg(MetricBundle bundle) {
    assertTrue(bundle.min_.deriv() < bundle.avg_.deriv());
    assertTrue(bundle.avg_.deriv() < bundle.max_.deriv());
  }

  public void testLessOrEqual(MetricBundle a, MetricBundle b) {
    assertTrue(a.min_.deriv() <= b.min_.deriv());
    assertTrue(a.max_.deriv() <= b.max_.deriv());
    assertTrue(a.avg_.deriv() <= b.avg_.deriv());
  }

  public void testMetrics() {

    MetricBundle angleSpan = new MetricBundle(
        S2Projections.MIN_ANGLE_SPAN, S2Projections.MAX_ANGLE_SPAN, S2Projections.AVG_ANGLE_SPAN);
    MetricBundle width =
        new MetricBundle(S2Projections.MIN_WIDTH, S2Projections.MAX_WIDTH, S2Projections.AVG_WIDTH);
    MetricBundle edge =
        new MetricBundle(S2Projections.MIN_EDGE, S2Projections.MAX_EDGE, S2Projections.AVG_EDGE);
    MetricBundle diag =
        new MetricBundle(S2Projections.MIN_DIAG, S2Projections.MAX_DIAG, S2Projections.AVG_DIAG);
    MetricBundle area =
        new MetricBundle(S2Projections.MIN_AREA, S2Projections.MAX_AREA, S2Projections.AVG_AREA);

    // First, check that min <= avg <= max for each metric.
    testMinMaxAvg(angleSpan);
    testMinMaxAvg(width);
    testMinMaxAvg(edge);
    testMinMaxAvg(diag);
    testMinMaxAvg(area);

    // Check that the maximum aspect ratio of an individual cell is consistent
    // with the global minimums and maximums.
    assertTrue(S2Projections.MAX_EDGE_ASPECT >= 1.0);
    assertTrue(S2Projections.MAX_EDGE_ASPECT
        < S2Projections.MAX_EDGE.deriv() / S2Projections.MIN_EDGE.deriv());
    assertTrue(S2Projections.MAX_DIAG_ASPECT >= 1);
    assertTrue(S2Projections.MAX_DIAG_ASPECT
        < S2Projections.MAX_DIAG.deriv() / S2Projections.MIN_DIAG.deriv());

    // Check various conditions that are provable mathematically.
    testLessOrEqual(width, angleSpan);
    testLessOrEqual(width, edge);
    testLessOrEqual(edge, diag);

    assertTrue(S2Projections.MIN_AREA.deriv()
        >= S2Projections.MIN_WIDTH.deriv() * S2Projections.MIN_EDGE.deriv() - 1e-15);
    assertTrue(S2Projections.MAX_AREA.deriv()
        < S2Projections.MAX_WIDTH.deriv() * S2Projections.MAX_EDGE.deriv() + 1e-15);

    // GetMinLevelForLength() and friends have built-in assertions, we just need
    // to call these functions to test them.
    //
    // We don't actually check that the metrics are correct here, e.g. that
    // GetMinWidth(10) is a lower bound on the width of cells at level 10.
    // It is easier to check these properties in s2cell_unittest, since
    // S2Cell has methods to compute the cell vertices, etc.

    for (int level = -2; level <= S2CellId.MAX_LEVEL + 3; ++level) {
      double dWidth = (2 * S2Projections.MIN_WIDTH.deriv()) * Math.pow(2, -level);
      if (level >= S2CellId.MAX_LEVEL + 3) {
        dWidth = 0;
      }

      // Check boundary cases (exactly equal to a threshold value).
      int expectedLevel = Math.max(0, Math.min(S2CellId.MAX_LEVEL, level));
      assertEquals(S2Projections.MIN_WIDTH.getMinLevel(dWidth), expectedLevel);
      assertEquals(S2Projections.MIN_WIDTH.getMaxLevel(dWidth), expectedLevel);
      assertEquals(S2Projections.MIN_WIDTH.getClosestLevel(dWidth), expectedLevel);

      // Also check non-boundary cases.
      assertEquals(S2Projections.MIN_WIDTH.getMinLevel(1.2 * dWidth), expectedLevel);
      assertEquals(S2Projections.MIN_WIDTH.getMaxLevel(0.8 * dWidth), expectedLevel);
      assertEquals(S2Projections.MIN_WIDTH.getClosestLevel(1.2 * dWidth), expectedLevel);
      assertEquals(S2Projections.MIN_WIDTH.getClosestLevel(0.8 * dWidth), expectedLevel);

      // Same thing for area1.
      double area1 = (4 * S2Projections.MIN_AREA.deriv()) * Math.pow(4, -level);
      if (level <= -3) {
        area1 = 0;
      }
      assertEquals(S2Projections.MIN_AREA.getMinLevel(area1), expectedLevel);
      assertEquals(S2Projections.MIN_AREA.getMaxLevel(area1), expectedLevel);
      assertEquals(S2Projections.MIN_AREA.getClosestLevel(area1), expectedLevel);
      assertEquals(S2Projections.MIN_AREA.getMinLevel(1.2 * area1), expectedLevel);
      assertEquals(S2Projections.MIN_AREA.getMaxLevel(0.8 * area1), expectedLevel);
      assertEquals(S2Projections.MIN_AREA.getClosestLevel(1.2 * area1), expectedLevel);
      assertEquals(S2Projections.MIN_AREA.getClosestLevel(0.8 * area1), expectedLevel);
    }
  }

  public void testExp() {
    for (int i = 0; i < 10; ++i) {
      assertEquals(i + 1, S2.exp(Math.pow(2, i)));
    }

    for (int i = 0; i < 10; ++i) {
      assertEquals(i + 1, S2.exp(-Math.pow(2, i)));
    }

    assertEquals(0, S2.exp(0));
    assertEquals(2, S2.exp(3));
    assertEquals(3, S2.exp(5));
  }
}
