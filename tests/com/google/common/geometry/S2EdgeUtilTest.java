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

import com.google.common.annotations.GwtCompatible;
import com.google.common.annotations.GwtIncompatible;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2EdgeUtil.FaceSegment;
import com.google.common.geometry.S2EdgeUtil.LongitudePruner;
import com.google.common.geometry.S2EdgeUtil.WedgeRelation;

/**
 * Tests for {@link S2EdgeUtil}.
 *
 */
@GwtCompatible(emulated = true)
public strictfp class S2EdgeUtilTest extends GeometryTestCase {
  /**
   * S2.robustCrossing() is allowed to return 0 or -1 if either edge is degenerate. We use this
   * value to represent this possibility.
   */
  private static final int DEGEN = -2;

  /** Maximum allowed error in latlng calculations. */
  private static final S2LatLng RECT_ERROR = S2EdgeUtil.RectBounder.maxErrorForTests();

  /**
   * Helper that asserts actual equals expected, unless we expect a degenerate result, in which case
   * we check that actual is either of the possible degenerate values.
   */
  private static void compareResult(int actual, int expected) {
    if (expected == DEGEN) {
      assertTrue(actual <= 0);
    } else {
      assertEquals(expected, actual);
    }
  }

  private void checkCrossing(S2Point a, S2Point b, S2Point c, S2Point d, int robust,
      boolean edgeOrVertex, boolean simple) {
    a = S2Point.normalize(a);
    b = S2Point.normalize(b);
    c = S2Point.normalize(c);
    d = S2Point.normalize(d);
    compareResult(S2EdgeUtil.robustCrossing(a, b, c, d), robust);
    if (simple) {
      assertEquals(robust > 0, S2EdgeUtil.simpleCrossing(a, b, c, d));
    }
    EdgeCrosser crosser = new EdgeCrosser(a, b, c);
    compareResult(crosser.robustCrossing(d), robust);
    compareResult(crosser.robustCrossing(c), robust);

    assertEquals(edgeOrVertex, S2EdgeUtil.edgeOrVertexCrossing(a, b, c, d));
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(d));
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(c));

    // Check that the crosser can be re-used.
    crosser = new EdgeCrosser(c, d);
    crosser.restartAt(a);
    compareResult(crosser.robustCrossing(b), robust);
    compareResult(crosser.robustCrossing(a), robust);
  }

  private void checkCrossings(S2Point a, S2Point b, S2Point c, S2Point d,
      int robust, boolean edgeOrVertex, boolean simple) {
    checkCrossing(a, b, c, d, robust, edgeOrVertex, simple);
    checkCrossing(b, a, c, d, robust, edgeOrVertex, simple);
    checkCrossing(a, b, d, c, robust, edgeOrVertex, simple);
    checkCrossing(b, a, d, c, robust, edgeOrVertex, simple);
    checkCrossing(a, a, c, d, DEGEN, false, false);
    checkCrossing(a, b, c, c, DEGEN, false, false);
    checkCrossing(a, b, a, b, 0, true, false);
    checkCrossing(c, d, a, b, robust, edgeOrVertex ^ (robust == 0), simple);
  }

  /**
   * The real tests of edge crossings are in loop and polygon unit tests, but we do a few simple
   * tests here.
   */
  public void testCrossings () {
    // Two regular edges that cross.
    checkCrossings(new S2Point(1, 2, 1), new S2Point(1, -3, 0.5),
        new S2Point(1, -0.5, -3), new S2Point(0.1, 0.5, 3), 1, true, true);

    // Two regular edges that cross antipodal points.
    checkCrossings(new S2Point(1, 2, 1), new S2Point(1, -3, 0.5),
        new S2Point(-1, 0.5, 3), new S2Point(-0.1, -0.5, -3), -1, false, true);

    // Two edges on the same great circle.
    checkCrossings(new S2Point(0, 0, -1), new S2Point(0, 1, 0),
        new S2Point(0, 1, 1), new S2Point(0, 0, 1), -1, false, true);

    // Two edges that cross where one vertex is S2.origin().
    checkCrossings(new S2Point(1, 0, 0), S2.origin(),
        new S2Point(1, -0.1, 1), new S2Point(1, 1, -0.1), 1, true, true);

    // Two edges that cross antipodal points where one vertex is S2.origin().
    checkCrossings(new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        new S2Point(0, 0, -1), new S2Point(-1, -1, 1), -1, false, true);

    // Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
    // and edge CD is further CCW around (2,3,4) than AB.
    checkCrossings(new S2Point(2, 3, 4), new S2Point(-1, 2, 5),
        new S2Point(7, -2, 3), new S2Point(2, 3, 4), 0, false, true);

    // Two edges that barely cross each other near the middle of one edge.  The
    // edge AB is approximately in the x=y plane, while CD is approximately
    // perpendicular to it and ends exactly at the x=y plane.
    checkCrossings(new S2Point(1, 1, 1), new S2Point(1, nextAfter(1d, 0), -1),
        new S2Point(11, -12, -1), new S2Point(10, 10, 1), 1, true, false);

    // In this version, the edges are separated by a distance of about 1e-15.
    checkCrossings(new S2Point(1, 1, 1), new S2Point(1, nextAfter(1d, 2), -1),
        new S2Point(1, -1, 0), new S2Point(1, 1, 0), -1, false, false);

    // In this version, the edges are separated by a distance of about 1e-640.
    checkCrossings(new S2Point(0, 0, 1), new S2Point(2, 1e-323, 1),
        new S2Point(1, -1, 1), new S2Point(1e-323, 0, 1), -1, false, false);

    // Two edges that barely cross each other near the middle of one edge.
    // Computing the exact determinant of some of the triangles in this test
    // requires more than 2000 bits of precision.
    checkCrossings(new S2Point(1, -1e-323, -1e-323), new S2Point(1e-323, 1, 1e-323),
        new S2Point(1, -1, 1e-323), new S2Point(1, 1, 0),
        1, true, false);

    /*
     * Tests in the C code that can't succeed in Java, since we don't have better-than-double
     * precision arithmetic:
     *
     * // Two edges that barely cross each other near the end of both edges.  This
     * // example cannot be handled using regular double-precision arithmetic due
     * // to floating-point underflow.
     * checkCrossings(new S2Point(0, 0, 1), new S2Point(2, -1e-323, 1),
     *     new S2Point(1, -1, 1), new S2Point(1e-323, 0, 1), 1, true, false);
     *
     * // In this version, the edges are separated by a distance of about 1e-640.
     * checkCrossings(new S2Point(1, 1e-323, -1e-323), new S2Point(-1e-323, 1, 1e-323),
     *     new S2Point(1, -1, 1e-323), new S2Point(1, 1, 0), -1, false, false);
     */
  }

  private static S2LatLngRect getEdgeBound(S2Point a, S2Point b) {
    S2EdgeUtil.RectBounder bounder = new S2EdgeUtil.RectBounder();
    bounder.addPoint(a);
    bounder.addPoint(b);
    return bounder.getBound();
  }

  private static S2LatLngRect getEdgeBound(
      double x1, double y1, double z1, double x2, double y2, double z2) {
    return getEdgeBound(
        S2Point.normalize(new S2Point(x1, y1, z1)),
        S2Point.normalize(new S2Point(x2, y2, z2)));
  }

  public void testMaxLatitudeSimple() {
    // Check cases where the min/max latitude is attained at a vertex.
    final double kCubeLat = Math.asin(1 / Math.sqrt(3));  // 35.26 degrees
    assertTrue(getEdgeBound(1, 1, 1, 1, -1, -1).approxEquals(
        new S2LatLngRect(new R1Interval(-kCubeLat, kCubeLat),
            new S1Interval(-S2.M_PI_4, S2.M_PI_4)), RECT_ERROR));
    assertTrue(getEdgeBound(1, -1, 1, 1, 1, -1).approxEquals(
        new S2LatLngRect(new R1Interval(-kCubeLat, kCubeLat),
            new S1Interval(-S2.M_PI_4, S2.M_PI_4)), RECT_ERROR));

    // Check cases where the min/max latitude occurs in the edge interior.
    // These tests expect the result to be pretty close to the middle of the
    // allowable error range (i.e., by adding 0.5 * kRectError).

    // Max latitude, CW edge
    assertDoubleNear(S2.M_PI_4 + 0.5 * RECT_ERROR.lat().radians(),
        getEdgeBound(1, 1, 1, 1, -1, 1).lat().hi());
    // Max latitude, CCW edge
    assertDoubleNear(S2.M_PI_4 + 0.5 * RECT_ERROR.lat().radians(),
        getEdgeBound(1, -1, 1, 1, 1, 1).lat().hi());
    // Min latitude, CW edge
    assertDoubleNear(-S2.M_PI_4 - 0.5 * RECT_ERROR.lat().radians(),
        getEdgeBound(1, -1, -1, -1, -1, -1).lat().lo());
    // Min latitude, CCW edge
    assertDoubleNear(-S2.M_PI_4 - 0.5 * RECT_ERROR.lat().radians(),
        getEdgeBound(-1, 1, -1, -1, -1, -1).lat().lo());

    // Check cases where the edge passes through one of the poles.
    assertEquals(S2.M_PI_2, getEdgeBound(.3, .4, 1, -.3, -.4, 1).lat().hi());
    assertEquals(-S2.M_PI_2, getEdgeBound(.3, .4, -1, -.3, -.4, -1).lat().lo());
  }

  public void testMaxLatitudeRandom() {
    // Check that the maximum latitude of edges is computed accurately to within
    // 3 * DBL_EPSILON (the expected maximum error).  We concentrate on maximum
    // latitudes near the equator and north pole since these are the extremes.
    double latRadTolerance = RECT_ERROR.lat().radians();

    final int kIters = 100;
    for (int iter = 0; iter < kIters; ++iter) {
      // Construct a right-handed coordinate frame (U,V,W) such that U points
      // slightly above the equator, V points at the equator, and W is slightly
      // offset from the north pole.
      S2Point u = randomPoint();
      // log is uniform
      u = S2Point.normalize(new S2Point(u.x, u.y,
          S2.DBL_EPSILON * 1e-6 * Math.pow(1e12, rand.nextDouble())));
      S2Point v = S2Point.normalize(S2.robustCrossProd(new S2Point(0, 0, 1), u));
      S2Point w = S2Point.normalize(S2.robustCrossProd(u, v));

      // Construct a line segment AB that passes through U, and check that the
      // maximum latitude of this segment matches the latitude of U.
      S2Point a = S2Point.normalize(S2Point.sub(u, S2Point.mul(v, rand.nextDouble())));
      S2Point b = S2Point.normalize(S2Point.add(u, S2Point.mul(v, rand.nextDouble())));
      S2LatLngRect abBound = getEdgeBound(a, b);
      assertDoubleNear(S2LatLng.latitude(u).radians(), abBound.lat().hi(), latRadTolerance);

      // Construct a line segment CD that passes through W, and check that the
      // maximum latitude of this segment matches the latitude of W.
      S2Point c = S2Point.normalize(S2Point.sub(w, S2Point.mul(v, rand.nextDouble())));
      S2Point d = S2Point.normalize(S2Point.add(w, S2Point.mul(v, rand.nextDouble())));
      S2LatLngRect cdBound = getEdgeBound(c, d);
      assertDoubleNear(S2LatLng.latitude(w).radians(), cdBound.lat().hi(), latRadTolerance);
    }
  }

  private S2Point perturbATowardsB(S2Point a, S2Point b) {
    double choice = rand.nextDouble();
    if (choice < 0.1) {
      return a;
    }
    if (choice < 0.3) {
      // Return a point that is exactly proportional to A and that still
      // satisfies S2.isUnitLength().
      for (;;) {
        S2Point temp = S2Point.mul(a,
            2 - a.norm() + 5 * (rand.nextDouble() - 0.5) * S2.DBL_EPSILON);
        if (!temp.equals(a) && S2.isUnitLength(temp)) {
          return temp;
        }
      }
    }
    if (choice < 0.5) {
      // Return a point such that the distance squared to A will underflow.
      return S2EdgeUtil.interpolateAtDistance(S1Angle.radians(1e-300), a, b);
    }
    // Otherwise return a point whose distance from A is near DBL_EPSILON such
    // that the log of the pdf is uniformly distributed.
    double distance = S2.DBL_EPSILON * 1e-5 * Math.pow(1e6, rand.nextDouble());
    return S2EdgeUtil.interpolateAtDistance(S1Angle.radians(distance), a, b);
  }

  private S2Point randomPole() {
    return new S2Point(0, 0, oneIn(2) ? 1 : -1);
  }

  private S2Point pointNearPole() {
    return perturbATowardsB(randomPole(), randomPoint());
  }

  private S2Point pointNearEquator() {
    return perturbATowardsB(S2Point.normalize(
        new S2Point(rand.nextDouble(), rand.nextDouble(), 0)), randomPole());
  }

  public void testNearlyIdenticalOrAntipodalPoints() {
    // Test pairs of points that are either:
    //  - identical
    //  - nearly or exactly proportional, e.g. (1,0,0) vs. (1+2e-16, 0, 0)
    //  - very close to each other
    // Furthermore we want to test cases where the two points are:
    //  - on a nearly-polar great circle
    //  - on a nearly-equatorial great circle
    //  - near the poles, but on any great circle
    //  - near the equator, but on any great circle
    //  - positioned arbitrarily
    // Also test the corresponding situations for antipodal points, i.e. by
    // negating one of the points so that they are almost 180 degrees apart.
    final int kIters = 10000;
    for (int iter = 0; iter < kIters; ++iter) {
      S2Point a, b;
      switch (rand.nextInt(5)) {
        case 0:
          // Two nearby points on a nearly-polar great circle.
          a = randomPoint();
          b = perturbATowardsB(a, pointNearPole());
          break;
        case 1:
          // Two nearby points on a nearly-equatorial great circle.
          a = pointNearEquator();
          b = perturbATowardsB(a, pointNearEquator());
          break;
        case 2:
          // Two nearby points near a pole, but on any great circle.
          a = pointNearPole();
          b = perturbATowardsB(a, randomPoint());
          break;
        case 3:
          // Two nearby points near the equator, but on any great circle.
          a = pointNearEquator();
          b = perturbATowardsB(a, randomPoint());
          break;
        case 4:
          // Two nearby points anywhere on the sphere.
          a = randomPoint();
          b = perturbATowardsB(a, randomPoint());
          break;
        default:
          throw new IllegalStateException("Can't happen.");
      }

      // The two points are chosen to be so close to each other that the min/max
      // latitudes are nearly always achieved at the edge endpoints.  The only
      // thing we need to watch out for is that the latitude error bound is
      // slightly larger if the min/max latitude occurs in the edge interior.
      S2LatLngRect expectedBound = S2LatLngRect.fromPointPair(new S2LatLng(a), new S2LatLng(b));
      S2LatLngRect bound = getEdgeBound(a, b);
      assertTrue(bound.contains(expectedBound));
      assertTrue(expectedBound.expanded(RECT_ERROR).polarClosure().contains(bound));

      // If the two points are close enough and one point is negated (antipodal
      // points), the bound should be the entire sphere.
      if (S2Point.crossProd(S2Point.sub(a, b), S2Point.add(a, b)).norm()
          <= 6.110 * S2.DBL_EPSILON) {
        assertEquals(S2LatLngRect.full(), getEdgeBound(a, S2Point.neg(b)));
      }
    }
  }

  private S2LatLngRect getSubregionBound(double xLat, double xLng, double yYat, double yLng) {
    S2LatLngRect in = S2LatLngRect.fromPointPair(
        S2LatLng.fromRadians(xLat, xLng),
        S2LatLng.fromRadians(yYat, yLng));
    S2LatLngRect out = S2EdgeUtil.RectBounder.expandForSubregions(in);

    // Test that the bound is actually expanded.
    assertTrue(out.contains(in));
    if (in.lat().equals(S2LatLngRect.fullLat())) {
      assertFalse(in.lat().contains(out.lat()));
    }

    return out;
  }

  public void testExpandForSubregions() {
    // First we check the various situations where the bound contains
    // nearly-antipodal points.  The tests are organized into pairs where the
    // two bounds are similar except that the first bound meets the
    // nearly-antipodal criteria while the second does not.

    // Cases where the bound does not straddle the equator (but almost does),
    // and spans nearly 180 degrees in longitude.
    assertTrue(getSubregionBound(3e-16, 0, 1e-14, S2.M_PI).isFull());
    assertFalse(getSubregionBound(9e-16, 0, 1e-14, S2.M_PI).isFull());
    assertTrue(getSubregionBound(1e-16, 7e-16, 1e-14, S2.M_PI).isFull());
    assertFalse(getSubregionBound(3e-16, 14e-16, 1e-14, S2.M_PI).isFull());
    assertTrue(getSubregionBound(1e-100, 14e-16, 1e-14, S2.M_PI).isFull());
    assertFalse(getSubregionBound(1e-100, 22e-16, 1e-14, S2.M_PI).isFull());

    // Cases where the bound spans at most 90 degrees in longitude, and almost
    // 180 degrees in latitude.  Note that DBL_EPSILON is about 2.22e-16, which
    // implies that the double-precision value just below Pi/2 can be written as
    // (S2.M_PI_2 - 2e-16).
    assertTrue(getSubregionBound(-S2.M_PI_2, -1e-15, S2.M_PI_2 - 7e-16, 0).
        isFull());
    assertFalse(getSubregionBound(-S2.M_PI_2, -1e-15, S2.M_PI_2 - 30e-16, 0).
        isFull());
    assertTrue(getSubregionBound(-S2.M_PI_2 + 4e-16, 0, S2.M_PI_2 - 2e-16, 1e-7).
        isFull());
    assertFalse(getSubregionBound(-S2.M_PI_2 + 30e-16, 0, S2.M_PI_2, 1e-7).
        isFull());
    assertTrue(getSubregionBound(-S2.M_PI_2 + 4e-16, 0, S2.M_PI_2 - 4e-16, S2.M_PI_2).
        isFull());
    assertFalse(getSubregionBound(-S2.M_PI_2, 0, S2.M_PI_2 - 30e-16, S2.M_PI_2).
        isFull());

    // Cases where the bound straddles the equator and spans more than 90
    // degrees in longitude.  These are the cases where the critical distance is
    // between a corner of the bound and the opposite longitudinal edge.  Unlike
    // the cases above, here the bound may contain nearly-antipodal points (to
    // within 3.055 * DBL_EPSILON) even though the latitude and longitude ranges
    // are both significantly less than (Pi - 3.055 * DBL_EPSILON).
    assertTrue(getSubregionBound(-S2.M_PI_2, 0, S2.M_PI_2 - 1e-8, S2.M_PI - 1e-7).isFull());
    assertFalse(getSubregionBound(-S2.M_PI_2, 0, S2.M_PI_2 - 1e-7, S2.M_PI - 1e-7).isFull());
    assertTrue(getSubregionBound(-S2.M_PI_2 + 1e-12, -S2.M_PI + 1e-4, S2.M_PI_2, 0).isFull());
    assertTrue(getSubregionBound(-S2.M_PI_2 + 1e-11, -S2.M_PI + 1e-4, S2.M_PI_2, 0).isFull());

    // Now we test cases where the bound does not contain nearly-antipodal
    // points, but it does contain points that are approximately 180 degrees
    // apart in latitude.
    assertTrue(getSubregionBound(1.5, -S2.M_PI_2, 1.5, S2.M_PI_2 - 2e-16).approxEquals(
        new S2LatLngRect(new R1Interval(1.5, 1.5), S1Interval.full()),
        RECT_ERROR));
    assertTrue(getSubregionBound(1.5, -S2.M_PI_2, 1.5, S2.M_PI_2 - 7e-16).approxEquals(
        new S2LatLngRect(new R1Interval(1.5, 1.5), new S1Interval(-S2.M_PI_2, S2.M_PI_2 - 7e-16)),
        RECT_ERROR));

    // Test the full and empty bounds.
    assertTrue(S2EdgeUtil.RectBounder.expandForSubregions(
        S2LatLngRect.full()).isFull());
    assertTrue(S2EdgeUtil.RectBounder.expandForSubregions(
        S2LatLngRect.empty()).isEmpty());

    // Check for cases where the bound is expanded to include one of the poles.
    assertTrue(getSubregionBound(-S2.M_PI_2 + 1e-15, 0, -S2.M_PI_2 + 1e-15, 0).approxEquals(
        new S2LatLngRect(new R1Interval(-S2.M_PI_2, -S2.M_PI_2 + 1e-15), S1Interval.full()),
        RECT_ERROR));
    assertTrue(getSubregionBound(S2.M_PI_2 - 1e-15, 0, S2.M_PI_2 - 1e-15, 0).approxEquals(
        new S2LatLngRect(new R1Interval(S2.M_PI_2 - 1e-15, S2.M_PI_2), S1Interval.full()),
        RECT_ERROR));
  }

  public void testLongitudePruner () {
    LongitudePruner pruner1 = new LongitudePruner(
        new S1Interval(0.75 * S2.M_PI, -0.75 * S2.M_PI), new S2Point(0, 1, 2));
    assertFalse(pruner1.intersects(new S2Point(1, 1, 3)));
    assertTrue(pruner1.intersects(new S2Point(-1 - 1e-15, -1, 0)));
    assertTrue(pruner1.intersects(new S2Point(-1, 0, 0)));
    assertTrue(pruner1.intersects(new S2Point(-1, 0, 0)));
    assertTrue(pruner1.intersects(new S2Point(1, -1, 8)));
    assertFalse(pruner1.intersects(new S2Point(1, 0, -2)));
    assertTrue(pruner1.intersects(new S2Point(-1, -1e-15, 0)));

    LongitudePruner pruner2 = new LongitudePruner(
        new S1Interval(0.25 * S2.M_PI, 0.25 * S2.M_PI), new S2Point(1, 0, 0));
    assertFalse(pruner2.intersects(new S2Point(2, 1, 2)));
    assertTrue(pruner2.intersects(new S2Point(1, 2, 3)));
    assertFalse(pruner2.intersects(new S2Point(0, 1, 4)));
    assertFalse(pruner2.intersects(new S2Point(-1e-15, -1, -1)));
  }

  private void checkWedge(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2,
      boolean contains, boolean intersects,
      S2EdgeUtil.WedgeRelation wedgeRelation) {
    a0 = S2Point.normalize(a0);
    ab1 = S2Point.normalize(ab1);
    a2 = S2Point.normalize(a2);
    b0 = S2Point.normalize(b0);
    b2 = S2Point.normalize(b2);
    assertEquals(contains, new S2EdgeUtil.WedgeContains().test(a0, ab1, a2, b0, b2) == 1);
    assertEquals(intersects, new S2EdgeUtil.WedgeIntersects().test(a0, ab1, a2, b0, b2) == -1);
    assertEquals(wedgeRelation, S2EdgeUtil.getWedgeRelation(a0, ab1, a2, b0, b2));
  }

  public void testWedges () {
    // For simplicity, all of these tests use an origin of (0, 0, 1).
    // This shouldn't matter as long as the lower-level primitives are
    // implemented correctly.

    // Intersection in one wedge.
    checkWedge(new S2Point(-1, 0, 10), new S2Point(0, 0, 1), new S2Point(1, 2, 10),
        new S2Point(0, 1, 10), new S2Point(1, -2, 10),
        false, true, WedgeRelation.WEDGE_PROPERLY_OVERLAPS);
    // Intersection in two wedges.
    checkWedge(new S2Point(-1, -1, 10), new S2Point(0, 0, 1), new S2Point(1, -1, 10),
        new S2Point(1, 0, 10), new S2Point(-1, 1, 10),
        false, true, WedgeRelation.WEDGE_PROPERLY_OVERLAPS);

    // Normal containment.
    checkWedge(new S2Point(-1, -1, 10), new S2Point(0, 0, 1), new S2Point(1, -1, 10),
        new S2Point(-1, 0, 10), new S2Point(1, 0, 10),
        true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
    // Containment with equality on one side.
    checkWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(-1, -1, 10),
        new S2Point(2, 1, 10), new S2Point(1, -5, 10),
        true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);
    // Containment with equality on the other side.
    checkWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(-1, -1, 10),
        new S2Point(1, -2, 10), new S2Point(-1, -1, 10),
        true, true, WedgeRelation.WEDGE_PROPERLY_CONTAINS);

    // Containment with equality on both sides.
    checkWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10), new S2Point(4, -5, 10),
        true, true, WedgeRelation.WEDGE_EQUALS);

    // Disjoint with equality on one side.
    checkWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(4, -5, 10),
        new S2Point(4, -5, 10), new S2Point(-2, -3, 10),
        false, false, WedgeRelation.WEDGE_IS_DISJOINT);
    // Disjoint with equality on the other side.
    checkWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(0, 5, 10),
        new S2Point(4, -5, 10), new S2Point(-2, 3, 10),
        false, false, WedgeRelation.WEDGE_IS_DISJOINT);
    // Disjoint with equality on both sides.
    checkWedge(new S2Point(-2, 3, 10), new S2Point(0, 0, 1), new S2Point(4, -5, 10),
        new S2Point(4, -5, 10), new S2Point(-2, 3, 10),
        false, false, WedgeRelation.WEDGE_IS_DISJOINT);

    // B contains A with equality on one side.
    checkWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(1, -5, 10),
        new S2Point(2, 1, 10), new S2Point(-1, -1, 10),
        false, true, WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
    // B contains A with equality on the other side.
    checkWedge(new S2Point(2, 1, 10), new S2Point(0, 0, 1), new S2Point(1, -5, 10),
        new S2Point(-2, 1, 10), new S2Point(1, -5, 10),
        false, true, WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
  }

  /**
   * Given a point X and an edge AB, check that the distance from X to AB is "distanceRadians" and
   * the closest point on AB is "expectedClosest".
   */
  private void checkDistance(S2Point x, S2Point a, S2Point b,
      double distanceRadians, S2Point expectedClosest) {
    x = S2Point.normalize(x);
    a = S2Point.normalize(a);
    b = S2Point.normalize(b);
    expectedClosest = S2Point.normalize(expectedClosest);
    // Note the 1e-15 epsilon is not needed by the C++ library, but we need it here due to use of
    // strictfp.
    assertEquals(distanceRadians, S2EdgeUtil.getDistance(x, a, b).radians(), 1E-15);
    S2Point closest = S2EdgeUtil.getClosestPoint(x, a, b);
    if (expectedClosest.equals(S2Point.ORIGIN)) {
      // This special value says that the result should be A or B.
      assertTrue(closest.equals(a) || closest.equals(b));
    } else {
      assertTrue(S2.approxEquals(closest, expectedClosest));
    }
    assertEquals(S1ChordAngle.ZERO, S2EdgeUtil.updateMinDistance(x, a, b, S1ChordAngle.ZERO));
    S1ChordAngle minDistance = S2EdgeUtil.updateMinDistance(x, a, b, S1ChordAngle.INFINITY);
    assertEquals(distanceRadians, minDistance.toAngle().radians(), 1e-15);
  }

  public void testDistance () {
    checkDistance(new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        0, new S2Point(1, 0, 0));
    checkDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        0, new S2Point(0, 1, 0));
    checkDistance(new S2Point(1, 3, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        0, new S2Point(1, 3, 0));
    checkDistance(new S2Point(0, 0, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        S2.M_PI_2, new S2Point(1, 0, 0));
    checkDistance(new S2Point(0, 0, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        S2.M_PI_2, new S2Point(1, 0, 0));
    checkDistance(new S2Point(-1, -1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        0.75 * S2.M_PI, new S2Point(0, 0, 0));

    checkDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(1, 1, 0),
        S2.M_PI_4, new S2Point(1, 1, 0));
    checkDistance(new S2Point(0, -1, 0), new S2Point(1, 0, 0), new S2Point(1, 1, 0),
        S2.M_PI_2, new S2Point(1, 0, 0));

    checkDistance(new S2Point(0, -1, 0), new S2Point(1, 0, 0), new S2Point(-1, 1, 0),
        S2.M_PI_2, new S2Point(1, 0, 0));
    checkDistance(new S2Point(-1, -1, 0), new S2Point(1, 0, 0), new S2Point(-1, 1, 0),
        S2.M_PI_2, new S2Point(-1, 1, 0));

    checkDistance(new S2Point(1, 1, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        Math.asin(Math.sqrt(1. / 3)), new S2Point(1, 1, 0));
    checkDistance(new S2Point(1, 1, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        Math.asin(Math.sqrt(1. / 3)), new S2Point(1, 1, 0));

    checkDistance(new S2Point(-1, 0, 0), new S2Point(1, 1, 0), new S2Point(1, 1, 0),
        0.75 * S2.M_PI, new S2Point(1, 1, 0));
    checkDistance(new S2Point(0, 0, -1), new S2Point(1, 1, 0), new S2Point(1, 1, 0),
        S2.M_PI_2, new S2Point(1, 1, 0));
    checkDistance(new S2Point(-1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0),
        S2.M_PI, new S2Point(1, 0, 0));
  }

  private void checkInterpolate(double t, S2Point a, S2Point b, S2Point expected) {
    a = S2Point.normalize(a);
    b = S2Point.normalize(b);
    expected = S2Point.normalize(expected);
    S2Point actual = S2EdgeUtil.interpolate(t, a, b);

    // We allow a bit more than the usual 1e-15 error tolerance because interpolate() uses trig
    // functions.
    assertTrue(S2.approxEquals(expected, actual, 3e-15));
  }

  public void testInterpolate () {
    // A zero-length edge.
    checkInterpolate(0, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));
    checkInterpolate(1, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));

    // Start, end, and middle of a medium-length edge.
    checkInterpolate(0, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(1, 0, 0));
    checkInterpolate(1, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(0, 1, 0));
    checkInterpolate(0.5, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(1, 1, 0));

    // Test that interpolation is done using distances on the sphere rather than
    // linear distances.
    checkInterpolate(
        1. / 3, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(Math.sqrt(3), 1, 0));
    checkInterpolate(
        2. / 3, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(1, Math.sqrt(3), 0));

    // Test that interpolation is accurate on a long edge (but not so long that
    // the definition of the edge itself becomes too unstable).
    final double kLng = S2.M_PI - 1e-2;
    S2Point a = S2LatLng.fromRadians(0, 0).toPoint();
    S2Point b = S2LatLng.fromRadians(0, kLng).toPoint();
    for (double f = 0.4; f > 1e-15; f *= 0.1) {
      checkInterpolate(f, a, b, S2LatLng.fromRadians(0, f * kLng).toPoint());
      checkInterpolate(1 - f, a, b, S2LatLng.fromRadians(0, (1 - f) * kLng).toPoint());
    }
  }

  public void testInterpolateCanExtrapolate () {
    final S2Point i = new S2Point(1, 0, 0);
    final S2Point j = new S2Point(0, 1, 0);

    // Initial vectors at 90 degrees.
    checkInterpolate(0, i, j, new S2Point(1, 0, 0));
    checkInterpolate(1, i, j, new S2Point(0, 1, 0));
    checkInterpolate(1.5, i, j, new S2Point(-1, 1, 0));
    checkInterpolate(2, i, j, new S2Point(-1, 0, 0));
    checkInterpolate(3, i, j, new S2Point(0, -1, 0));
    checkInterpolate(4, i, j, new S2Point(1, 0, 0));

    // Negative values of t.
    checkInterpolate(-1, i, j, new S2Point(0, -1, 0));
    checkInterpolate(-2, i, j, new S2Point(-1, 0, 0));
    checkInterpolate(-3, i, j, new S2Point(0, 1, 0));
    checkInterpolate(-4, i, j, new S2Point(1, 0, 0));

    // Initial vectors at 45 degrees.
    checkInterpolate(2, i, new S2Point(1, 1, 0), new S2Point(0, 1, 0));
    checkInterpolate(3, i, new S2Point(1, 1, 0), new S2Point(-1, 1, 0));
    checkInterpolate(4, i, new S2Point(1, 1, 0), new S2Point(-1, 0, 0));

    // Initial vectors at 135 degrees.
    checkInterpolate(2, i, new S2Point(-1, 1, 0), new S2Point(0, -1, 0));

    // Take a small fraction along the curve.
    S2Point p = S2EdgeUtil.interpolate(0.001, i, j);
    // We should get back where we started.
    checkInterpolate(1000, i, p, j);
  }

  public void testRepeatedInterpolation () {
    // Check that points do not drift away from unit length when repeated
    // interpolations are done.
    for (int i = 0; i < 100; ++i) {
      S2Point a = randomPoint();
      S2Point b = randomPoint();
      for (int j = 0; j < 1000; ++j) {
        a = S2EdgeUtil.interpolate(0.01, a, b);
      }
      assertTrue(S2.isUnitLength(a));
    }
  }

  public void testIntersectionTolerance () {
    // We repeatedly construct two edges that cross near a random point "p",
    // and measure the distance from the actual intersection point "x" to the
    // the expected intersection point "p" and also to the edges that cross
    // near "p".
    //
    // Note that GetIntersection() does not guarantee that "x" and "p" will be
    // close together (since the intersection point is numerically unstable
    // when the edges cross at a very small angle), but it does guarantee that
    // "x" will be close to both of the edges that cross.
    S1Angle maxPointDist = new S1Angle();
    S1Angle maxEdgeDist = new S1Angle();

    for (int i = 0; i < 1000; ++i) {
      // We construct two edges AB and CD that intersect near "p".  The angle
      // between AB and CD (expressed as a slope) is chosen randomly between
      // 1e-15 and 1.0 such that its logarithm is uniformly distributed.  This
      // implies that small values are much more likely to be chosen.
      //
      // Once the slope is chosen, the four points ABCD must be offset from P
      // by at least (1e-15 / slope) so that the points are guaranteed to have
      // the correct circular ordering around P.  This is the distance from P
      // at which the two edges are separated by about 1e-15, which is
      // approximately the minimum distance at which we can expect computed
      // points on the two lines to be distinct and have the correct ordering.
      //
      // The actual offset distance from P is chosen randomly in the range
      // [1e-15 / slope, 1.0], again uniformly distributing the logarithm.
      // This ensures that we test both long and very short segments that
      // intersect at both large and very small angles.

      Matrix3x3 frame = S2.getFrame(randomPoint());
      S2Point p = frame.getCol(0);
      S2Point d1 = frame.getCol(1);
      S2Point d2 = frame.getCol(2);

      double slope = Math.pow(1e-15, rand.nextDouble());
      d2 = S2Point.add(d1, S2Point.mul(d2, slope));
      S2Point a = S2Point.normalize(
          S2Point.add(p, S2Point.mul(d1, Math.pow(1e-15 / slope, rand.nextDouble()))));
      S2Point b = S2Point.normalize(
          S2Point.sub(p, S2Point.mul(d1, Math.pow(1e-15 / slope, rand.nextDouble()))));
      S2Point c = S2Point.normalize(
          S2Point.add(p, S2Point.mul(d2, Math.pow(1e-15 / slope, rand.nextDouble()))));
      S2Point d = S2Point.normalize(
          S2Point.sub(p, S2Point.mul(d2, Math.pow(1e-15 / slope, rand.nextDouble()))));
      S2Point x = S2EdgeUtil.getIntersection(a, b, c, d);
      S1Angle abDist = S2EdgeUtil.getDistance(x, a, b);
      S1Angle cdDist = S2EdgeUtil.getDistance(x, c, d);
      assertTrue(abDist.lessOrEquals(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
      assertTrue(cdDist.lessOrEquals(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
      maxEdgeDist = max(maxEdgeDist, max(abDist, cdDist));
      maxPointDist = max(maxPointDist, new S1Angle(p, x));
    }
  }

  /** Returns the larger of a or b, or if one of them is null the other one is returned. */
  private static S1Angle max(S1Angle a, S1Angle b) {
    return a.compareTo(b) >= 0 ? a : b;
  }

  private boolean isEdgeBNearEdgeA(String aStr, String bStr, double maxErrorDegrees) {
    S2Polyline a = makePolyline(aStr);
    assertEquals(2, a.numVertices());
    S2Polyline b = makePolyline(bStr);
    assertEquals(2, b.numVertices());
    return isEdgeBNearEdgeA(a.vertex(0), a.vertex(1), b.vertex(0), b.vertex(1),
        S1Angle.degrees(maxErrorDegrees));
  }

  // TODO(eengle): This exists in s2edgeutil.cc, but seems to only be needed for tests?!?
  // If relocated to S2EdgeUtil.java, be sure to comment out assertions.
  private static boolean isEdgeBNearEdgeA(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, S1Angle tolerance) {
    assert (tolerance.radians() < S2.M_PI / 2);
    assert (tolerance.radians() > 0);

    // The point on edge B=b0b1 furthest from edge A=a0a1 is either b0, b1, or
    // some interior point on B.  If it is an interior point on B, then it must be
    // one of the two points where the great circle containing B (circ(B)) is
    // furthest from the great circle containing A (circ(A)).  At these points,
    // the distance between circ(B) and circ(A) is the angle between the planes
    // containing them.
    S2Point aOrtho = S2Point.normalize(S2.robustCrossProd(a0, a1));
    S2Point aNearestB0 = S2EdgeUtil.getClosestPoint(b0, a0, a1, aOrtho);
    S2Point aNearestB1 = S2EdgeUtil.getClosestPoint(b1, a0, a1, aOrtho);

    // If aNearestB0 and aNearestB1 have opposite orientation from a0 and a1,
    // we invert aOrtho so that it points in the same direction as aNearestB0 x
    // aNearestB1.  This helps us handle the case where A and B are oppositely
    // oriented but otherwise might be near each other.  We check orientation and
    // invert rather than computing aNearestB0 x aNearestB1 because those two
    // points might be equal, and have an unhelpful cross product.
    if (S2.robustCCW(aOrtho, aNearestB0, aNearestB1) < 0) {
      aOrtho = S2Point.neg(aOrtho);
    }

    // To check if all points on B are within tolerance of A, we first check to
    // see if the endpoints of B are near A.  If they are not, B is not near A.
    S1Angle b0Distance = new S1Angle(b0, aNearestB0);
    S1Angle b1Distance = new S1Angle(b1, aNearestB1);
    if (b0Distance.greaterThan(tolerance) || b1Distance.greaterThan(tolerance)) {
      return false;
    }

    // If b0 and b1 are both within tolerance of A, we check to see if the angle
    // between the planes containing B and A is greater than tolerance.  If it is
    // not, no point on B can be further than tolerance from A (recall that we
    // already know that b0 and b1 are close to A, and S2Edges are all shorter
    // than 180 degrees).  The angle between the planes containing circ(A) and
    // circ(B) is the angle between their normal vectors.
    S2Point bOrtho = S2Point.normalize(S2.robustCrossProd(b0, b1));
    S1Angle planarAngle = new S1Angle(aOrtho, bOrtho);
    if (planarAngle.lessOrEquals(tolerance)) {
      return true;
    }

    // As planarAngle approaches M_PI, the projection of aOrtho onto the plane
    // of B approaches the null vector, and normalizing it is numerically
    // unstable.  This makes it unreliable or impossible to identify pairs of
    // points where circ(A) is furthest from circ(B).  At this point in the
    // algorithm, this can only occur for two reasons:
    //
    //  1.) b0 and b1 are closest to A at distinct endpoints of A, in which case
    //      the opposite orientation of aOrtho and bOrtho means that A and B are
    //      in opposite hemispheres and hence not close to each other.
    //
    //  2.) b0 and b1 are closest to A at the same endpoint of A, in which case
    //      the orientation of aOrtho was chosen arbitrarily to be that of a0
    //      cross a1.  B must be shorter than 2*tolerance and all points in B are
    //      close to one endpoint of A, and hence to A.
    //
    // The logic applies when planarAngle is robustly greater than M_PI/2, but
    // may be more computationally expensive than the logic beyond, so we choose a
    // value close to M_PI.
    if (planarAngle.radians() >= S2.M_PI - 0.01) {
      return (new S1Angle(b0, a0).lessThan(new S1Angle(b0, a1)))
          == (new S1Angle(b1, a0).lessThan(new S1Angle(b1, a1)));
    }

    // Finally, if either of the two points on circ(B) where circ(B) is furthest
    // from circ(A) lie on edge B, edge B is not near edge A.
    //
    // The normalized projection of aOrtho onto the plane of circ(B) is one of
    // the two points along circ(B) where it is furthest from circ(A).  The other
    // is -1 times the normalized projection.
    S2Point furthest = S2Point.normalize(
        S2Point.sub(aOrtho, S2Point.mul(bOrtho, aOrtho.dotProd(bOrtho))));
    assert (S2.isUnitLength(furthest));
    S2Point furthestInv = S2Point.neg(furthest);

    // A point p lies on B if you can proceed from bOrtho to b0 to p to b1 and
    // back to bOrtho without ever turning right.  We test this for furthest and
    // furthestInv, and return true if neither point lies on B.
    return !((S2.robustCCW(bOrtho, b0, furthest) > 0
            && S2.robustCCW(furthest, b1, bOrtho) > 0)
        || (S2.robustCCW(bOrtho, b0, furthestInv) > 0
            && S2.robustCCW(furthestInv, b1, bOrtho) > 0));
  }

  public void testEdgeBNearEdgeA () {
    // Edge is near itself.
    assertTrue(isEdgeBNearEdgeA("5:5, 10:-5", "5:5, 10:-5", 1e-6));

    // Edge is near its reverse
    assertTrue(isEdgeBNearEdgeA("5:5, 10:-5", "10:-5, 5:5", 1e-6));

    // Short edge is near long edge.
    assertTrue(isEdgeBNearEdgeA("10:0, -10:0", "2:1, -2:1", 1.0));

    // Long edges cannot be near shorter edges.
    assertFalse(isEdgeBNearEdgeA("2:1, -2:1", "10:0, -10:0", 1.0));

    // Orthogonal crossing edges are not near each other...
    assertFalse(isEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 1.0));

    // ... unless all points on B are within tolerance of A.
    assertTrue(isEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 2.0));

    // Very long edges whose endpoints are close may have interior points that are
    // far apart.  An implementation that only considers the vertices of polylines
    // will incorrectly consider such edges as "close" when they are not.
    // Consider, for example, two consecutive lines of longitude.  As they
    // approach the poles, they become arbitrarily close together, but along the
    // equator they bow apart.
    assertFalse(isEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 0.5));
    assertTrue(isEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 1.5));

    // The two arcs here are nearly as long as S2 edges can be (just shy of 180
    // degrees), and their endpoints are less than 1 degree apart.  Their
    // midpoints, however, are at opposite ends of the sphere along its equator.
    assertFalse(isEdgeBNearEdgeA(
        "0:-179.75, 0:-0.25", "0:179.75, 0:0.25", 1.0));

    // At the equator, the second arc here is 9.75 degrees from the first, and
    // closer at all other points.  However, the southern point of the second arc
    // (-1, 9.75) is too far from the first arc for the short-circuiting logic in
    // isEdgeBNearEdgeA to apply.
    assertTrue(isEdgeBNearEdgeA("40:0, -5:0", "39:0.975, -1:0.975", 1.0));

    // Same as above, but B's orientation is reversed, causing the angle between
    // the normal vectors of circ(B) and circ(A) to be (180-9.75) = 170.5 degrees,
    // not 9.75 degrees.  The greatest separation between the planes is still 9.75
    // degrees.
    assertTrue(isEdgeBNearEdgeA("10:0, -10:0", "-.4:0.975, 0.4:0.975", 1.0));

    // A and B are on the same great circle, A and B partially overlap, but the
    // only part of B that does not overlap A is shorter than tolerance.
    assertTrue(isEdgeBNearEdgeA("0:0, 1:0", "0.9:0, 1.1:0", 0.25));

    // A and B are on the same great circle, all points on B are close to A at its
    // second endpoint, (1,0).
    assertTrue(isEdgeBNearEdgeA("0:0, 1:0", "1.1:0, 1.2:0", 0.25));

    // Same as above, but B's orientation is reversed.  This case is special
    // because the projection of the normal defining A onto the plane containing B
    // is the null vector, and must be handled by a special case.
    assertTrue(isEdgeBNearEdgeA("0:0, 1:0", "1.2:0, 1.1:0", 0.25));
  }

  public void testCollinearEdgesThatDontTouch () {
    final int kIters = 500;
    for (int iter = 0; iter < kIters; ++iter) {
      S2Point a = randomPoint();
      S2Point d = randomPoint();
      S2Point b = S2EdgeUtil.interpolate(0.05, a, d);
      S2Point c = S2EdgeUtil.interpolate(0.95, a, d);
      assertTrue(0 > S2EdgeUtil.robustCrossing(a, b, c, d));
      EdgeCrosser crosser = new EdgeCrosser(a, b, c);
      assertTrue(0 > crosser.robustCrossing(d));
      assertTrue(0 > crosser.robustCrossing(c));
    }
  }

  // TODO(eengle): testCoincidentZeroLengthEdgesThatDontTouch requires S2.expensiveCCW use arbitrary
  // precision arithmetic. This has been ported as part of the S2EdgeUtilTest refresh, but commented
  // out until we can update S2.java.
  /*
  // Returns 0 or a power of 2, based on a randomly selected value drawn from a skewed distribution.
  private double randomSkewedCoord() {
    int binaryExp = skewed(11);
    return (binaryExp > 1022) ? 0 : Math.pow(2, -binaryExp);
  }

  public void testCoincidentZeroLengthEdgesThatDontTouch () {
    // It is important that the edge primitives can handle vertices that are
    // exactly proportional to each other, i.e. that are not identical but are
    // nevertheless exactly coincident when projected onto the unit sphere.
    // There are various ways that such points can arise.  For example,
    // normalize() itself is not idempotent: there exist distinct points A,B
    // such that normalize(A) == B  and normalize(B) == A.  Another issue is
    // that sometimes calls to normalize() are skipped when the result of a
    // calculation "should" be unit length mathematically (e.g., when computing
    // the cross product of two orthonormal vectors).
    //
    // This test checks pairs of edges AB and CD where A,B,C,D are exactly
    // coincident on the sphere and the norms of A,B,C,D are monotonically
    // increasing.  Such edge pairs should never intersect.  (This is not
    // obvious, since it depends on the particular symbolic perturbations used
    // by S2.robustCCW().  It would be better to replace this with a test that
    // says that the CCW results must be consistent with each other.)
    final int kIters = 1000;
    for (int iter = 0; iter < kIters; ++iter) {
      // Construct a point P where every component is zero or a power of 2.
      S2Point p = new S2Point(randomSkewedCoord(), randomSkewedCoord(), randomSkewedCoord());

      // If all components were zero, try again.  Note that normalization may
      // convert a non-zero point into a zero one due to underflow (!)
      p = S2Point.normalize(p);
      if (p.equals(S2Point.ORIGIN)) {
        --iter;
        continue;
      }

      // Now every non-zero component should have exactly the same mantissa.
      // This implies that if we scale the point by an arbitrary factor, every
      // non-zero component will still have the same mantissa.  Scale the points
      // so that they are all distinct and are still very likely to satisfy
      // S2.isUnitLength (which allows for a small amount of error in the norm).
      S2Point a = S2Point.mul(p, 1-3e-16);
      S2Point b = S2Point.mul(p, 1-1e-16);
      S2Point c = p;
      S2Point d = S2Point.mul(p, 1+2e-16);
      if (!S2.isUnitLength(a) || !S2.isUnitLength(d)) {
        --iter;
        continue;
      }

      // Verify that the expected edges do not cross.
      assertTrue(0 > S2EdgeUtil.robustCrossing(a, b, c, d));
      EdgeCrosser crosser = new EdgeCrosser(a, b, c);
      assertTrue(0 > crosser.robustCrossing(d));
      assertTrue(0 > crosser.robustCrossing(c));
    }
  }
  */

  private void checkFaceClipping(S2Point aRaw, S2Point bRaw) {
    S2Point a = S2Point.normalize(aRaw);
    S2Point b = S2Point.normalize(bRaw);
    // TODO(user): Remove the following line once S2.robustCrossProd is
    // extended to use simulation of simplicity.
    if (a.equals(S2Point.neg(b))) {
      return;
    }

    // First we test GetFaceSegments.
    FaceSegment[] segments = new FaceSegment[6];
    int n = S2EdgeUtil.getFaceSegments(a, b, segments);
    assertTrue(n >= 1);

    R2Rect biunit = new R2Rect(new R1Interval(-1, 1), new R1Interval(-1, 1));

    // The first and last vertices should approximately equal A and B.
    assertTrue(a.angle(S2Projections.faceUvToXyz(segments[0].face, segments[0].a))
        <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
    assertTrue(b.angle(S2Projections.faceUvToXyz(segments[n - 1].face, segments[n - 1].b))
        <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);

    S2Point norm = S2Point.normalize(S2.robustCrossProd(a, b));
    S2Point aTangent = S2Point.crossProd(norm, a);
    S2Point bTangent = S2Point.crossProd(b, norm);
    for (int i = 0; i < n; ++i) {
      // Vertices may not protrude outside the biunit square.
      assertTrue(biunit.contains(segments[i].a));
      assertTrue(biunit.contains(segments[i].b));
      if (i == 0) {
        continue;
      }

      // The two representations of each interior vertex (on adjacent faces)
      // must correspond to exactly the same S2Point.
      assertNotSame(segments[i - 1].face, segments[i].face);
      assertEquals(S2Projections.faceUvToXyz(segments[i - 1].face, segments[i - 1].b),
          S2Projections.faceUvToXyz(segments[i].face, segments[i].a));

      // Interior vertices should be in the plane containing A and B, and should
      // be contained in the wedge of angles between A and B (i.e., the dot
      // products with aTangent and bTangent should be non-negative).
      S2Point p = S2Point.normalize(S2Projections.faceUvToXyz(segments[i].face, segments[i].a));
      assertTrue(Math.abs(p.dotProd(norm)) <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
      assertTrue(p.dotProd(aTangent) >= -S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
      assertTrue(p.dotProd(bTangent) >= -S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
    }

    // Now we test ClipToPaddedFace (sometimes with a padding of zero).  We do
    // this by defining an (x,y) coordinate system for the plane containing AB,
    // and converting points along the great circle AB to angles in the range
    // [-Pi, Pi].  We then accumulate the angle intervals spanned by each
    // clipped edge; the union over all 6 faces should approximately equal the
    // interval covered by the original edge.
    double padding = oneIn(10) ? 0.0 : 1e-10 * Math.pow(1e-5, rand.nextDouble());
    S2Point xAxis = a, yAxis = aTangent;
    S1Interval expectedAngles = new S1Interval(0, a.angle(b));
    S1Interval maxAngles = expectedAngles.expanded(S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
    S1Interval actualAngles = null;
    for (int face = 0; face < 6; ++face) {
      R2Vector aUv = new R2Vector(), bUv = new R2Vector();
      if (S2EdgeUtil.clipToPaddedFace(a, b, face, padding, aUv, bUv)) {
        S2Point aClip = S2Point.normalize(S2Projections.faceUvToXyz(face, aUv));
        S2Point bClip = S2Point.normalize(S2Projections.faceUvToXyz(face, bUv));
        assertTrue(Math.abs(aClip.dotProd(norm)) <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
        assertTrue(Math.abs(bClip.dotProd(norm)) <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
        if (aClip.angle(a) > S2EdgeUtil.FACE_CLIP_ERROR_RADIANS) {
          assertDoubleNear(1 + padding, Math.max(Math.abs(aUv.x), Math.abs(aUv.y)));
        }
        if (bClip.angle(b) > S2EdgeUtil.FACE_CLIP_ERROR_RADIANS) {
          assertDoubleNear(1 + padding, Math.max(Math.abs(bUv.x), Math.abs(bUv.y)));
        }
        double aAngle = Math.atan2(aClip.dotProd(yAxis), aClip.dotProd(xAxis));
        double bAngle = Math.atan2(bClip.dotProd(yAxis), bClip.dotProd(xAxis));
        // Rounding errors may cause bAngle to be slightly less than aAngle.
        // We handle this by constructing the interval with FromPointPair(),
        // which is okay since the interval length is much less than S2.M_PI.
        S1Interval faceAngles = S1Interval.fromPointPair(aAngle, bAngle);
        assertTrue(maxAngles.contains(faceAngles));
        actualAngles = actualAngles == null ? faceAngles : actualAngles.union(faceAngles);
      }
    }
    assertTrue(actualAngles.expanded(S2EdgeUtil.FACE_CLIP_ERROR_RADIANS).contains(expectedAngles));
  }

  /**
   * Calls {@link #checkFaceClipping(S2Point, S2Point)} with both [a,b] and [b,a] pairs, to ensure
   * symmetry.
   */
  private void checkFaceClippingEdgePair(S2Point a, S2Point b) {
    checkFaceClipping(a, b);
    checkFaceClipping(b, a);
  }

  /**
   * This function is designed to choose line segment endpoints that are difficult to handle
   * correctly. Given two adjacent cube vertices P and Q, it returns either an edge midpoint, face
   * midpoint, or corner vertex that is in the plane of PQ and that has been perturbed slightly. It
   * also sometimes returns a random point from anywhere on the sphere.
   */
  private S2Point perturbedCornerOrMidpoint(S2Point p, S2Point q) {
    S2Point a = S2Point.add(
        S2Point.mul(p, rand.nextInt(3) - 1),
        S2Point.mul(q, rand.nextInt(3) - 1));
    if (oneIn(10)) {
      // This perturbation often has no effect except on coordinates that are
      // zero, in which case the perturbed value is so small that operations on
      // it often result in underflow.
      a = S2Point.add(a, S2Point.mul(randomPoint(), Math.pow(1e-300, rand.nextDouble())));
    } else if (oneIn(2)) {
      // For coordinates near 1 (say > 0.5), this perturbation yields values
      // that are only a few representable values away from the initial value.
      a = S2Point.add(a, S2Point.mul(randomPoint(), 4 * S2.DBL_EPSILON));
    } else {
      // A perturbation whose magnitude is in the range [1e-25, 1e-10].
      a = S2Point.add(a, S2Point.mul(randomPoint(), 1e-10 * Math.pow(1e-15, rand.nextDouble())));
    }
    if (a.norm2() < Double.MIN_VALUE) {
      // If a.norm2() is denormalizedn, normalize() loses too much precision.
      return perturbedCornerOrMidpoint(p, q);
    }
    return a;
  }

  @GwtIncompatible("flaky for reasons to be determined")
  public void testFaceClipping() {
    // Start with a few simple cases.
    // An edge that is entirely contained within one cube face:
    checkFaceClippingEdgePair(new S2Point(1, -0.5, -0.5), new S2Point(1, 0.5, 0.5));
    // An edge that crosses one cube edge:
    checkFaceClippingEdgePair(new S2Point(1, 0, 0), new S2Point(0, 1, 0));
    // An edge that crosses two opposite edges of face 0:
    checkFaceClippingEdgePair(new S2Point(0.75, 0, -1), new S2Point(0.75, 0, 1));
    // An edge that crosses two adjacent edges of face 2:
    checkFaceClippingEdgePair(new S2Point(1, 0, 0.75), new S2Point(0, 1, 0.75));
    // An edges that crosses three cube edges (four faces):
    checkFaceClippingEdgePair(new S2Point(1, 0.9, 0.95), new S2Point(-1, 0.95, 0.9));

    // Comprehensively test edges that are difficult to handle, especially those
    // that nearly follow one of the 12 cube edges.
    R2Rect biunit = new R2Rect(new R1Interval(-1, 1), new R1Interval(-1, 1));
    final int kIters = 1000;  // Test passes with 1e6 iterations
    for (int iter = 0; iter < kIters; ++iter) {
      // Choose two adjacent cube corners P and Q.
      int face = rand.nextInt(6);
      int i = rand.nextInt(4);
      int j = (i + 1) & 3;
      S2Point p = S2Projections.faceUvToXyz(face, biunit.getVertex(i));
      S2Point q = S2Projections.faceUvToXyz(face, biunit.getVertex(j));

      // Now choose two points that are nearly in the plane of PQ, preferring
      // points that are near cube corners, face midpoints, or edge midpoints.
      S2Point a = perturbedCornerOrMidpoint(p, q);
      S2Point b = perturbedCornerOrMidpoint(p, q);
      checkFaceClipping(a, b);
    }
  }

  /**
   * Choose a random point in the rectangle defined by points A and B, sometimes returning a point
   * on the edge AB or the points A and B themselves.
   */
  private R2Vector chooseRectPoint(R2Vector a, R2Vector b) {
    if (oneIn(5)) {
      return oneIn(2) ? a : b;
    } else if (oneIn(3)) {
      return R2Vector.add(a, R2Vector.mul(R2Vector.sub(b, a), rand.nextDouble()));
    } else {
      return new R2Vector(uniform(a.x, b.x), uniform(a.y, b.y));
    }
  }

  /**
   * Given a point X on the line AB (which is checked), return the fraction "t" such that x =
   * (1-t)*a + t*b. Return 0 if A = B.
   */
  private double getFraction(R2Vector x, R2Vector a, R2Vector b) {
    // A bound for the error in edge clipping plus the error in the calculation
    // below (which is similar to IntersectsRect).
    final double kError = S2EdgeUtil.EDGE_CLIP_ERROR_UV_DIST
        + S2EdgeUtil.INTERSECTS_RECT_ERROR_UV_DIST;
    if (a.equals(b)) {
      return 0.0;
    }
    R2Vector dir = R2Vector.normalize(R2Vector.sub(b, a));
    assertTrue(Math.abs(R2Vector.sub(x, a).dotProd(dir.ortho())) <= kError);
    return R2Vector.sub(x, a).dotProd(dir);
  }

  /**
   * Given a point P representing a possibly clipped endpoint A of an edge AB, verify that "clip"
   * contains P, and that if clipping occurred (i.e., P != A) then P is on the boundary of "clip".
   */
  private void checkPointOnBoundary(R2Vector p, R2Vector a, R2Rect clip) {
    assertTrue(clip.contains(p));
    if (!p.equals(a)) {
      assertFalse(clip.contains(new R2Vector(nextAfter(p.x, a.x), nextAfter(p.y, a.y))));
    }
  }

  /**
   * Given an edge AB and a rectangle "clip", verify that IntersectsRect(), ClipEdge(), and
   * clipEdgeBound() produce consistent results.
   */
  private void checkClipEdge(R2Vector a, R2Vector b, R2Rect clip) {
    // A bound for the error in edge clipping plus the error in the
    // IntersectsRect calculation below.
    final double kError = S2EdgeUtil.EDGE_CLIP_ERROR_UV_DIST
        + S2EdgeUtil.INTERSECTS_RECT_ERROR_UV_DIST;
    R2Vector aClipped = new R2Vector();
    R2Vector bClipped = new R2Vector();
    if (!S2EdgeUtil.clipEdge(a, b, clip, aClipped, bClipped)) {
      assertFalse(S2EdgeUtil.intersectsRect(a, b, clip.expanded(-kError)));
    } else {
      assertTrue(S2EdgeUtil.intersectsRect(a, b, clip.expanded(kError)));
      // Check that the clipped points lie on the edge AB, and that the points
      // have the expected order along the segment AB.
      assertTrue(getFraction(aClipped, a, b) <= getFraction(bClipped, a, b));
      // Check that the clipped portion of AB is as large as possible.
      checkPointOnBoundary(aClipped, a, clip);
      checkPointOnBoundary(bClipped, b, clip);
    }

    // Choose a random initial bound to pass to clipEdgeBound.
    R2Rect initialClip = R2Rect.fromPointPair(chooseRectPoint(a, b),
        chooseRectPoint(a, b));
    R2Rect bound = S2EdgeUtil.getClippedEdgeBound(a, b, initialClip);
    if (bound.isEmpty()) {
      // Precondition of clipEdgeBound not met.
      return;
    }
    R2Rect maxBound = bound.intersection(clip);
    if (!S2EdgeUtil.clipEdgeBound(a, b, clip, bound)) {
      assertFalse(S2EdgeUtil.intersectsRect(a, b, maxBound.expanded(-kError)));
    } else {
      assertTrue(S2EdgeUtil.intersectsRect(a, b, maxBound.expanded(kError)));
      // Check that the bound is as large as possible.
      int ai = (a.x > b.x) ? 1 : 0;
      int aj = (a.y > b.y) ? 1 : 0;
      checkPointOnBoundary(bound.getVertex(ai, aj), a, maxBound);
      checkPointOnBoundary(bound.getVertex(1 - ai, 1 - aj), b, maxBound);
    }
  }

  /**
   * Given an interval "clip", randomly choose either a value in the interval, a value outside the
   * interval, or one of the two interval endpoints, ensuring that all cases have reasonable
   * probability for any interval "clip".
   */
  private double chooseEndpoint(R1Interval clip) {
    if (oneIn(5)) {
      return oneIn(2) ? clip.lo() : clip.hi();
    } else {
      switch (rand.nextInt(3)) {
        case 0:  return clip.lo() - rand.nextDouble();
        case 1:  return clip.hi() + rand.nextDouble();
        default: return clip.lo() + rand.nextDouble() * clip.getLength();
      }
    }
  }

  /**
   * Given a rectangle "clip", choose a point that may lie in the rectangle interior, along an
   * extended edge, exactly at a vertex, or in one of the eight regions exterior to "clip" that are
   * separated by its extended edges. Also sometimes return points that are exactly on one of the
   * extended diagonals of "clip". All cases are reasonably likely to occur for any given rectangle
   * "clip".
   */
  private R2Vector chooseEndpoint(R2Rect clip) {
    if (oneIn(10)) {
      // Return a point on one of the two extended diagonals.
      int diag = rand.nextInt(2);
      double t = uniform(-1, 2);
      return R2Vector.add(
          R2Vector.mul(clip.getVertex(diag), 1 - t),
          R2Vector.mul(clip.getVertex(diag + 2), t));
    } else {
      return new R2Vector(chooseEndpoint(clip.x()), chooseEndpoint(clip.y()));
    }
  }

  /**
   * Given a rectangle "clip", test the S2EdgeUtil edge clipping methods using many edges that are
   * randomly constructed to trigger special cases.
   */
  private void checkEdgeClipping(R2Rect clip) {
    final int kIters = 1000;  // Test passes with 1e6 iterations
    for (int iter = 0; iter < kIters; ++iter) {
      checkClipEdge(chooseEndpoint(clip), chooseEndpoint(clip), clip);
    }
  }

  public void testEdgeClipping () {
    // Test clipping against random rectangles.
    for (int i = 0; i < 5; ++i) {
      checkEdgeClipping(R2Rect.fromPointPair(
          new R2Vector(uniform(-1, 1), uniform(-1, 1)),
          new R2Vector(uniform(-1, 1), uniform(-1, 1))));
    }
    // Also clip against one-dimensional, singleton, and empty rectangles.
    checkEdgeClipping(new R2Rect(new R1Interval(-0.7, -0.7), new R1Interval(0.3, 0.35)));
    checkEdgeClipping(new R2Rect(new R1Interval(0.2, 0.5), new R1Interval(0.3, 0.3)));
    checkEdgeClipping(new R2Rect(new R1Interval(-0.7, 0.3), new R1Interval(0, 0)));
    checkEdgeClipping(R2Rect.fromPoint(new R2Vector(0.3, 0.8)));
    checkEdgeClipping(R2Rect.empty());
  }

  /**
   * We implement this in terms of DBL_EPSILON, because we don't have universal support for
   * Math.nextAfter, and it's not worth implementing in Platform.
   */
  private static double nextAfter(double x, double direction) {
    return direction < x ? x - S2.DBL_EPSILON : x + S2.DBL_EPSILON;
  }
}
