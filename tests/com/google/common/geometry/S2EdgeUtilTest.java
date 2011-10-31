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

import com.google.common.collect.ImmutableList;

/**
 * Tests for {@link S2EdgeUtil}.
 *
 */
public strictfp class S2EdgeUtilTest extends GeometryTestCase {

  public static final int DEGENERATE = -2;

  private void compareResult(int actual, int expected) {
    // HACK ALERT: RobustCrossing() is allowed to return 0 or -1 if either edge
    // is degenerate. We use the value kDegen to represent this possibility.
    if (expected == DEGENERATE) {
      assertTrue(actual <= 0);
    } else {
      assertEquals(expected, actual);
    }
  }

  private void assertCrossing(S2Point a,
      S2Point b,
      S2Point c,
      S2Point d,
      int robust,
      boolean edgeOrVertex,
      boolean simple) {
    a = S2Point.normalize(a);
    b = S2Point.normalize(b);
    c = S2Point.normalize(c);
    d = S2Point.normalize(d);

    compareResult(S2EdgeUtil.robustCrossing(a, b, c, d), robust);
    if (simple) {
      assertEquals(robust > 0, S2EdgeUtil.simpleCrossing(a, b, c, d));
    }
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a, b, c);
    compareResult(crosser.robustCrossing(d), robust);
    compareResult(crosser.robustCrossing(c), robust);

    assertEquals(S2EdgeUtil.edgeOrVertexCrossing(a, b, c, d), edgeOrVertex);
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(d));
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(c));
  }

  private void assertCrossings(S2Point a,
      S2Point b,
      S2Point c,
      S2Point d,
      int robust,
      boolean edgeOrVertex,
      boolean simple) {
    assertCrossing(a, b, c, d, robust, edgeOrVertex, simple);
    assertCrossing(b, a, c, d, robust, edgeOrVertex, simple);
    assertCrossing(a, b, d, c, robust, edgeOrVertex, simple);
    assertCrossing(b, a, d, c, robust, edgeOrVertex, simple);
    assertCrossing(a, a, c, d, DEGENERATE, false, false);
    assertCrossing(a, b, c, c, DEGENERATE, false, false);
    assertCrossing(a, b, a, b, 0, true, false);
    assertCrossing(c, d, a, b, robust, (edgeOrVertex ^ (robust == 0)), simple);
  }

  public void testCrossings() {
    // The real tests of edge crossings are in s2{loop,polygon}_unittest,
    // but we do a few simple tests here.

    // Two regular edges that cross.
    assertCrossings(new S2Point(1, 2, 1),
        new S2Point(1, -3, 0.5),
        new S2Point(1, -0.5, -3),
        new S2Point(0.1, 0.5, 3),
        1,
        true,
        true);

    // Two regular edges that cross antipodal points.
    assertCrossings(new S2Point(1, 2, 1),
        new S2Point(1, -3, 0.5),
        new S2Point(-1, 0.5, 3),
        new S2Point(-0.1, -0.5, -3),
        -1,
        false,
        true);

    // Two edges on the same great circle.
    assertCrossings(new S2Point(0, 0, -1),
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 1),
        new S2Point(0, 0, 1),
        -1,
        false,
        true);

    // Two edges that cross where one vertex is S2.Origin().
    assertCrossings(new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 0, 1),
        new S2Point(1, 1, -1),
        1,
        true,
        true);

    // Two edges that cross antipodal points where one vertex is S2.Origin().
    assertCrossings(new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 0, -1),
        new S2Point(-1, -1, 1),
        -1,
        false,
        true);

    // Two edges that share an endpoint. The Ortho() direction is (-4,0,2),
    // and edge CD is further CCW around (2,3,4) than AB.
    assertCrossings(new S2Point(2, 3, 4),
        new S2Point(-1, 2, 5),
        new S2Point(7, -2, 3),
        new S2Point(2, 3, 4),
        0,
        false,
        true);

    // Two edges that barely cross edge other.
    assertCrossings(new S2Point(1, 1, 1),
        new S2Point(1, 1 - 1e-15, -1),
        new S2Point(-1, -1, 0),
        new S2Point(1, 1, 0),
        1,
        true,
        false);
  }

  private S2LatLngRect getEdgeBound(double x1,
      double y1,
      double z1,
      double x2,
      double y2,
      double z2) {
    S2EdgeUtil.RectBounder bounder = new S2EdgeUtil.RectBounder();
    S2Point p1 = S2Point.normalize(new S2Point(x1, y1, z1));
    S2Point p2 = S2Point.normalize(new S2Point(x2, y2, z2));
    bounder.addPoint(p1);
    bounder.addPoint(p2);
    return bounder.getBound();
  }

  public void testRectBounder() {
    // Check cases where min/max latitude is not at a vertex.
    // Max, CW
    assertDoubleNear(getEdgeBound(1, 1, 1, 1, -1, 1).lat().hi(), S2.M_PI_4);
    // Max, CCW
    assertDoubleNear(getEdgeBound(1, -1, 1, 1, 1, 1).lat().hi(), S2.M_PI_4);
    // Min, CW
    assertDoubleNear(getEdgeBound(1, -1, -1, -1, -1, -1).lat().lo(), -S2.M_PI_4);
    // Min, CCW
    assertDoubleNear(getEdgeBound(-1, 1, -1, -1, -1, -1).lat().lo(), -S2.M_PI_4);

    // Check cases where the edge passes through one of the poles.
    assertDoubleNear(getEdgeBound(.3, .4, 1, -.3, -.4, 1).lat().hi(), S2.M_PI_2);
    assertDoubleNear(getEdgeBound(.3, .4, -1, -.3, -.4, -1).lat().lo(), -S2.M_PI_2);

    // Check cases where the min/max latitude is attained at a vertex.
    double kCubeLat = Math.asin(Math.sqrt(1. / 3)); // 35.26 degrees
    assertTrue(
        getEdgeBound(1, 1, 1, 1, -1, -1).lat().approxEquals(new R1Interval(-kCubeLat, kCubeLat)));
    assertTrue(
        getEdgeBound(1, -1, 1, 1, 1, -1).lat().approxEquals(new R1Interval(-kCubeLat, kCubeLat)));
  }

  // Produce a normalized S2Point for testing.
  private S2Point S2NP(double x, double y, double z) {
    return S2Point.normalize(new S2Point(x, y, z));
  }

  public void testXYZPruner() {
    S2EdgeUtil.XYZPruner pruner = new S2EdgeUtil.XYZPruner();

    // We aren't actually normalizing these points but it doesn't
    // matter too much as long as we are reasonably close to unit vectors.
    // This is a simple triangle on the equator.
    pruner.addEdgeToBounds(S2NP(0, 1, 0), S2NP(0.1, 1, 0));
    pruner.addEdgeToBounds(S2NP(0.1, 1, 0), S2NP(0.1, 1, 0.1));
    pruner.addEdgeToBounds(S2NP(0.1, 1, 0.1), S2NP(0, 1, 0));

    // try a loop around the triangle but far enough out to not overlap.
    pruner.setFirstIntersectPoint(S2NP(-0.1, 1.0, 0.0));
    assertFalse(pruner.intersects(S2NP(-0.1, 1.0, 0.2)));
    assertFalse(pruner.intersects(S2NP(0.0, 1.0, 0.2)));
    assertFalse(pruner.intersects(S2NP(0.2, 1.0, 0.2)));
    assertFalse(pruner.intersects(S2NP(0.2, 1.0, 0.05)));
    assertFalse(pruner.intersects(S2NP(0.2, 1.0, -0.1)));
    assertFalse(pruner.intersects(S2NP(-0.1, 1.0, -0.1)));
    assertFalse(pruner.intersects(S2NP(-0.1, 1.0, 0.0)));

    // now we go to a point in the bounding box of the triangle but well
    // out of the loop. This will be a hit even though it really does not
    // need to be.
    assertTrue(pruner.intersects(S2NP(0.02, 1.0, 0.04)));

    // now we zoom out to do an edge *just* below the triangle. This should
    // be a hit because we are within the deformation zone.
    assertTrue(pruner.intersects(S2NP(-0.1, 1.0, -0.03)));
    assertFalse(pruner.intersects(S2NP(0.05, 1.0, -0.03))); // not close
    assertTrue(pruner.intersects(S2NP(0.05, 1.0, -0.01))); // close
    assertTrue(pruner.intersects(S2NP(0.05, 1.0, 0.13)));
    assertFalse(pruner.intersects(S2NP(0.13, 1.0, 0.14)));

    // Create a new pruner with very small area and correspondingly narrow
    // deformation tolerances.
    S2EdgeUtil.XYZPruner spruner = new S2EdgeUtil.XYZPruner();
    spruner.addEdgeToBounds(S2NP(0, 1, 0.000), S2NP(0.001, 1, 0));
    spruner.addEdgeToBounds(S2NP(0.001, 1, 0.000), S2NP(0.001, 1, 0.001));
    spruner.addEdgeToBounds(S2NP(0.001, 1, 0.001), S2NP(0.000, 1, 0));

    spruner.setFirstIntersectPoint(S2NP(0, 1.0, -0.1));
    assertFalse(spruner.intersects(S2NP(0.0005, 1.0, -0.0005)));
    assertFalse(spruner.intersects(S2NP(0.0005, 1.0, -0.0005)));
    assertFalse(spruner.intersects(S2NP(0.0005, 1.0, -0.00001)));
    assertTrue(spruner.intersects(S2NP(0.0005, 1.0, -0.0000001)));
  }

  public void testLongitudePruner() {
    S2EdgeUtil.LongitudePruner pruner1 = new S2EdgeUtil.LongitudePruner(
        new S1Interval(0.75 * S2.M_PI, -0.75 * S2.M_PI), new S2Point(0, 1, 2));

    assertFalse(pruner1.intersects(new S2Point(1, 1, 3)));
    assertTrue(pruner1.intersects(new S2Point(-1 - 1e-15, -1, 0)));
    assertTrue(pruner1.intersects(new S2Point(-1, 0, 0)));
    assertTrue(pruner1.intersects(new S2Point(-1, 0, 0)));
    assertTrue(pruner1.intersects(new S2Point(1, -1, 8)));
    assertFalse(pruner1.intersects(new S2Point(1, 0, -2)));
    assertTrue(pruner1.intersects(new S2Point(-1, -1e-15, 0)));

    S2EdgeUtil.LongitudePruner pruner2 = new S2EdgeUtil.LongitudePruner(
        new S1Interval(0.25 * S2.M_PI, 0.25 * S2.M_PI), new S2Point(1, 0, 0));

    assertFalse(pruner2.intersects(new S2Point(2, 1, 2)));
    assertTrue(pruner2.intersects(new S2Point(1, 2, 3)));
    assertFalse(pruner2.intersects(new S2Point(0, 1, 4)));
    assertFalse(pruner2.intersects(new S2Point(-1e-15, -1, -1)));
  }

  private void assertWedge(S2Point a0,
      S2Point ab1,
      S2Point a2,
      S2Point b0,
      S2Point b2,
      boolean contains,
      boolean intersects,
      boolean crosses) {
    a0 = S2Point.normalize(a0);
    ab1 = S2Point.normalize(ab1);
    a2 = S2Point.normalize(a2);
    b0 = S2Point.normalize(b0);
    b2 = S2Point.normalize(b2);

    assertEquals(new S2EdgeUtil.WedgeContains().test(a0, ab1, a2, b0, b2), contains ? 1 : 0);
    assertEquals(new S2EdgeUtil.WedgeIntersects().test(a0, ab1, a2, b0, b2), intersects ? -1 : 0);
    assertEquals(new S2EdgeUtil.WedgeContainsOrIntersects().test(a0, ab1, a2, b0, b2),
        contains ? 1 : intersects ? -1 : 0);
    assertEquals(new S2EdgeUtil.WedgeContainsOrCrosses().test(a0, ab1, a2, b0, b2),
        contains ? 1 : crosses ? -1 : 0);
  }

  public void testWedges() {
    // For simplicity, all of these tests use an origin of (0, 0, 1).
    // This shouldn't matter as long as the lower-level primitives are
    // implemented correctly.

    // Intersection in one wedge.
    assertWedge(new S2Point(-1, 0, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, 2, 10),
        new S2Point(0, 1, 10),
        new S2Point(1, -2, 10),
        false,
        true,
        true);
    // Intersection in two wedges.
    assertWedge(new S2Point(-1, -1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -1, 10),
        new S2Point(1, 0, 10),
        new S2Point(-1, 1, 10),
        false,
        true,
        true);

    // Normal containment.
    assertWedge(new S2Point(-1, -1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -1, 10),
        new S2Point(-1, 0, 10),
        new S2Point(1, 0, 10),
        true,
        true,
        false);
    // Containment with equality on one side.
    assertWedge(new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(-1, -1, 10),
        new S2Point(2, 1, 10),
        new S2Point(1, -5, 10),
        true,
        true,
        false);
    // Containment with equality on the other side.
    assertWedge(new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(-1, -1, 10),
        new S2Point(1, -2, 10),
        new S2Point(-1, -1, 10),
        true,
        true,
        false);
    // Containment with equality on both sides.
    assertWedge(new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10),
        new S2Point(4, -5, 10),
        true,
        true,
        false);

    // Disjoint with equality on one side.
    assertWedge(new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(4, -5, 10),
        new S2Point(4, -5, 10),
        new S2Point(-2, -3, 10),
        false,
        false,
        false);
    // Disjoint with equality on the other side.
    assertWedge(new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(0, 5, 10),
        new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10),
        false,
        false,
        false);
    // Disjoint with equality on both sides.
    assertWedge(new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(4, -5, 10),
        new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10),
        false,
        false,
        false);

    // B contains A with equality on one side.
    assertWedge(new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -5, 10),
        new S2Point(2, 1, 10),
        new S2Point(-1, -1, 10),
        false,
        true,
        false);
    // B contains A with equality on the other side.
    assertWedge(new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -5, 10),
        new S2Point(-2, 1, 10),
        new S2Point(1, -5, 10),
        false,
        true,
        false);
  }

  public void testGetClosestPoint() {
    final double kMargin = 1e-6;

    S2Point a = S2LatLng.fromDegrees(-0.5, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(+0.5, 0).toPoint();

    // On edge at end points.
    assertEquals(a, S2EdgeUtil.getClosestPoint(a, a, b));
    assertEquals(b, S2EdgeUtil.getClosestPoint(b, a, b));

    // On edge in between.
    S2Point mid = S2LatLng.fromDegrees(0, 0).toPoint();
    assertEquals(mid, S2EdgeUtil.getClosestPoint(mid, a, b));

    // End points are closest
    assertEquals(a, S2EdgeUtil.getClosestPoint(S2LatLng.fromDegrees(-1, 0).toPoint(), a, b));
    assertEquals(b, S2EdgeUtil.getClosestPoint(S2LatLng.fromDegrees(+1, 0).toPoint(), a, b));

    // Intermediate point is closest.
    S2Point x = S2LatLng.fromDegrees(+0.1, 1).toPoint();
    S2Point expectedClosestPoint = S2LatLng.fromDegrees(+0.1, 0).toPoint();

    assertTrue(expectedClosestPoint.aequal(S2EdgeUtil.getClosestPoint(x, a, b), kMargin));
  }

  // Given a point X and an edge AB, check that the distance from X to AB is
  // "distanceRadians" and the closest point on AB is "expectedClosest".
  private static void checkDistance(
      S2Point x, S2Point a, S2Point b, double distanceRadians, S2Point expectedClosest) {
    final double kEpsilon = 1e-10;
    x = S2Point.normalize(x);
    a = S2Point.normalize(a);
    b = S2Point.normalize(b);
    expectedClosest = S2Point.normalize(expectedClosest);

    assertEquals(distanceRadians, S2EdgeUtil.getDistance(x, a, b).radians(), kEpsilon);

    S2Point closest = S2EdgeUtil.getClosestPoint(x, a, b);
    if (expectedClosest.equals(new S2Point(0, 0, 0))) {
      // This special value says that the result should be A or B.
      assertTrue(closest == a || closest == b);
    } else {
      assertTrue(S2.approxEquals(closest, expectedClosest));
    }
  }

  public void testGetDistance() {
    checkDistance(
        new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0), 0, new S2Point(1, 0, 0));
    checkDistance(
        new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0), 0, new S2Point(0, 1, 0));
    checkDistance(
        new S2Point(1, 3, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0), 0, new S2Point(1, 3, 0));
    checkDistance(new S2Point(0, 0, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0), Math.PI / 2,
        new S2Point(1, 0, 0));
    checkDistance(new S2Point(0, 0, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0), Math.PI / 2,
        new S2Point(1, 0, 0));
    checkDistance(new S2Point(-1, -1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        0.75 * Math.PI, new S2Point(0, 0, 0));
    checkDistance(new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(1, 1, 0), Math.PI / 4,
        new S2Point(1, 1, 0));
    checkDistance(new S2Point(0, -1, 0), new S2Point(1, 0, 0), new S2Point(1, 1, 0), Math.PI / 2,
        new S2Point(1, 0, 0));
    checkDistance(new S2Point(0, -1, 0), new S2Point(1, 0, 0), new S2Point(-1, 1, 0), Math.PI / 2,
        new S2Point(1, 0, 0));
    checkDistance(new S2Point(-1, -1, 0), new S2Point(1, 0, 0), new S2Point(-1, 1, 0), Math.PI / 2,
        new S2Point(-1, 1, 0));
    checkDistance(new S2Point(1, 1, 1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        Math.asin(Math.sqrt(1. / 3)), new S2Point(1, 1, 0));
    checkDistance(new S2Point(1, 1, -1), new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        Math.asin(Math.sqrt(1. / 3)), new S2Point(1, 1, 0));
    checkDistance(new S2Point(-1, 0, 0), new S2Point(1, 1, 0), new S2Point(1, 1, 0), 0.75 * Math.PI,
        new S2Point(1, 1, 0));
    checkDistance(new S2Point(0, 0, -1), new S2Point(1, 1, 0), new S2Point(1, 1, 0), Math.PI / 2,
        new S2Point(1, 1, 0));
    checkDistance(new S2Point(-1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0), Math.PI,
        new S2Point(1, 0, 0));
  }

  public void testIntersectionTolerance() {
    // We repeatedly construct two edges that cross near a random point "p",
    // and measure the distance from the actual intersection point "x" to the
    // the expected intersection point "p" and also to the edges that cross
    // near "p".
    //
    // Note that getIntersection() does not guarantee that "x" and "p" will be
    // close together (since the intersection point is numerically unstable
    // when the edges cross at a very small angle), but it does guarantee that
    // "x" will be close to both of the edges that cross.
    S1Angle maxPointDist = new S1Angle();
    S1Angle maxEdgeDist = new S1Angle();

    for (int i = 0; i < 1000; ++i) {
      // We construct two edges AB and CD that intersect near "p". The angle
      // between AB and CD (expressed as a slope) is chosen randomly between
      // 1e-15 and 1.0 such that its logarithm is uniformly distributed. This
      // implies that small values are much more likely to be chosen.
      //
      // Once the slope is chosen, the four points ABCD must be offset from P
      // by at least (1e-15 / slope) so that the points are guaranteed to have
      // the correct circular ordering around P. This is the distance from P
      // at which the two edges are separated by about 1e-15, which is
      // approximately the minimum distance at which we can expect computed
      // points on the two lines to be distinct and have the correct ordering.
      //
      // The actual offset distance from P is chosen randomly in the range
      // [1e-15 / slope, 1.0], again uniformly distributing the logarithm.
      // This ensures that we test both long and very short segments that
      // intersect at both large and very small angles.

      ImmutableList<S2Point> points = getRandomFrame();
      S2Point p = points.get(0);
      S2Point d1 = points.get(1);
      S2Point d2 = points.get(2);
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
      S1Angle distAb = S2EdgeUtil.getDistance(x, a, b);
      S1Angle distCd = S2EdgeUtil.getDistance(x, c, d);

      assertTrue(distAb.lessThan(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
      assertTrue(distCd.lessThan(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));

      // test getIntersection() post conditions
      assertTrue(S2.orderedCCW(a, x, b, S2Point.normalize(S2.robustCrossProd(a, b))));
      assertTrue(S2.orderedCCW(c, x, d, S2Point.normalize(S2.robustCrossProd(c, d))));

      maxEdgeDist = S1Angle.max(maxEdgeDist, S1Angle.max(distAb, distCd));
      maxPointDist = S1Angle.max(maxPointDist, new S1Angle(p, x));
    }
  }
}
