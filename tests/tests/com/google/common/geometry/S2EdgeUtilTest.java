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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2.M_PI_4;
import static com.google.common.geometry.S2EdgeUtil.INTERSECTION_ERROR;
import static com.google.common.geometry.S2Point.ORIGIN;
import static com.google.common.geometry.S2Point.X_POS;
import static com.google.common.geometry.S2Point.Y_POS;
import static com.google.common.geometry.S2Point.Z_NEG;
import static com.google.common.geometry.S2Point.Z_POS;
import static com.google.common.geometry.S2RobustCrossProd.robustCrossProd;
import static com.google.common.geometry.S2TextFormat.makePolyline;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashMultiset;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2EdgeUtil.FaceSegment;
import com.google.common.geometry.S2EdgeUtil.LongitudePruner;
import com.google.common.geometry.S2EdgeUtil.ResultError;
import com.google.common.geometry.S2EdgeUtil.WedgeRelation;
import java.math.BigDecimal;
import java.util.Collection;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Tests for {@link S2EdgeUtil}.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@RunWith(JUnit4.class)
public class S2EdgeUtilTest extends GeometryTestCase {

  /** Maximum allowed error in latlng calculations. */
  private static final S2LatLng RECT_ERROR = S2EdgeUtil.RectBounder.maxErrorForTests();

  /** The approximate maximum error in S2EdgeUtil.getDistance() for small distances. */
  private static final double GET_DISTANCE_ABS_ERROR = 3 * DBL_EPSILON;

  // The tests below verify that default or NaN S2Point arguments don't cause crashes when
  // assertions are disabled, as would be the case in production. Normally, assertions are expected
  // to be thrown on those invalid inputs, and these tests also verify that.

  private void checkCrossingSignInvalid(S2Point point, int expected) {
    // Calling the constructor with assertions enabled should throw an AssertionError.
    Assert.assertThrows(
        AssertionError.class,
        () -> {
          EdgeCrosser unused = new EdgeCrosser(point, point);
        });

    // Calling the constructor and robustCrossing() with assertions disabled should not crash or
    // throw an exception, but when assertions are enabled, robustCrossing() should throw an
    // AssertionError.
    EdgeCrosser crosser = uncheckedCreate(() -> new EdgeCrosser(point, point));
    Assert.assertThrows(AssertionError.class, () -> crosser.robustCrossing(point, point));
    assertEquals(expected, (long) uncheckedCreate(() -> crosser.robustCrossing(point, point)));
  }

  private void checkEdgeOrVertexCrossingInvalid(S2Point point, boolean expected) {
    Assert.assertThrows(
        AssertionError.class,
        () -> {
          EdgeCrosser unused = new EdgeCrosser(point, point);
        });

    EdgeCrosser crosser = uncheckedCreate(() -> new EdgeCrosser(point, point));
    Assert.assertThrows(AssertionError.class, () -> crosser.edgeOrVertexCrossing(point, point));
    assertEquals(expected, uncheckedCreate(() -> crosser.edgeOrVertexCrossing(point, point)));
  }

  private void checkSignedEdgeOrVertexCrossingInvalid(S2Point point, int expected) {
    Assert.assertThrows(
        AssertionError.class,
        () -> {
          EdgeCrosser unused = new EdgeCrosser(point, point);
        });

    EdgeCrosser crosser = uncheckedCreate(() -> new EdgeCrosser(point, point));
    Assert.assertThrows(
        AssertionError.class, () -> crosser.signedEdgeOrVertexCrossing(point, point));
    assertEquals(
        expected, (long) uncheckedCreate(() -> crosser.signedEdgeOrVertexCrossing(point, point)));
  }

  @Test // Regression test for b/393637044
  @SuppressWarnings("FloatingPointLiteralPrecision") // to be identical to the C++ test.
  public void distanceOptimizationIsConservative() {
    // Verifies that updateMinInteriorDistance() computes the lower bound on the true distance
    // conservatively. (This test used to fail.)
    S2Point x = new S2Point(-0.017952729194524016, -0.30232422079175203, 0.95303607751077712);
    S2Point a = new S2Point(-0.017894725505830295, -0.30229974986194175, 0.95304493075220664);
    S2Point b = new S2Point(-0.017986591360900289, -0.30233851195954353, 0.95303090543659963);
    S1ChordAngle minDistance = S2EdgeUtil.updateMinDistance(x, a, b, S1ChordAngle.INFINITY);

    assertTrue(minDistance.equals(
                   S2EdgeUtil.updateMinDistance(x, a, b, minDistance.successor())));
  }

  @Test
  public void testInvalidDefaultPoints() throws Exception {
    // Check that default-constructed S2Point arguments do throw assertions, when assertions are
    // enabled, but don't cause crashes when assertions are disabled.
    S2Point point = new S2Point();
    checkCrossingSignInvalid(point, 0);
    checkEdgeOrVertexCrossingInvalid(point, false);
    checkSignedEdgeOrVertexCrossingInvalid(point, 0);
  }

  @Test
  public void testInvalidNanPoints() throws Exception {
    // Check that NaN S2Point arguments do throw assertions, when assertions are  enabled, but
    // don't cause crashes when assertions are disabled.
    double nan = Double.NaN;  // In c++ this test uses std::numeric_limits<double>::quiet_NaN();
    S2Point point = new S2Point(nan, nan, nan);
    checkCrossingSignInvalid(point, -1);
    checkEdgeOrVertexCrossingInvalid(point, false);
    checkSignedEdgeOrVertexCrossingInvalid(point, 0);
  }

  private void checkCrossing(
      S2Point a,
      S2Point b,
      S2Point c,
      S2Point d,
      int crossingSign,
      int signedCrossingSign,
      boolean checkSimple) {
    // For degenerate edges, robustCrossing() is documented to return 0 if two vertices from
    // different edges are the same and -1 otherwise. The checkCrossings() function below uses
    // various argument permutations that can sometimes create this case, so we fix it now if
    // necessary.
    if (a.equalsPoint(c) || a.equalsPoint(d) || b.equalsPoint(c) || b.equalsPoint(d)) {
      crossingSign = 0;
    }

    // As a sanity check, make sure that the expected value of "signedCrossingSign" is consistent
    // with its documented properties.
    if (crossingSign == 1) {
      assertEquals(signedCrossingSign, S2Predicates.sign(a, b, c));
    } else if (crossingSign == 0 && signedCrossingSign != 0) {
      assertEquals(signedCrossingSign, (a.equalsPoint(c) || b.equalsPoint(d)) ? 1 : -1);
    }

    assertEquals(crossingSign, S2EdgeUtil.robustCrossing(a, b, c, d));
    boolean edgeOrVertex = signedCrossingSign != 0;
    assertEquals(edgeOrVertex, S2EdgeUtil.edgeOrVertexCrossing(a, b, c, d));

    EdgeCrosser crosser = new EdgeCrosser(a, b, c);
    assertEquals(crossingSign, crosser.robustCrossing(d));
    assertEquals(crossingSign, crosser.robustCrossing(c));
    assertEquals(crossingSign, crosser.robustCrossing(d, c));
    assertEquals(crossingSign, crosser.robustCrossing(c, d));

    crosser.restartAt(c);
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(d));
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(c));
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(d, c));
    assertEquals(edgeOrVertex, crosser.edgeOrVertexCrossing(c, d));

    crosser.restartAt(c);
    assertEquals(signedCrossingSign, crosser.signedEdgeOrVertexCrossing(d));
    assertEquals(-signedCrossingSign, crosser.signedEdgeOrVertexCrossing(c));
    assertEquals(-signedCrossingSign, crosser.signedEdgeOrVertexCrossing(d, c));
    assertEquals(signedCrossingSign, crosser.signedEdgeOrVertexCrossing(c, d));

    // Check that the crosser can be re-used.
    crosser = new EdgeCrosser(c, d);
    crosser.restartAt(a);
    assertEquals(crossingSign, crosser.robustCrossing(b));
    assertEquals(crossingSign, crosser.robustCrossing(a));

    if (checkSimple) {
      assertEquals(crossingSign > 0, S2EdgeUtil.simpleCrossing(a, b, c, d));
    }
  }

  private void checkCrossings(
      S2Point a,
      S2Point b,
      S2Point c,
      S2Point d,
      int crossingSign,
      int signedCrossingSign,
      boolean checkSimple) {
    a = a.normalize();
    b = b.normalize();
    c = c.normalize();
    d = d.normalize();

    checkCrossing(a, b, c, d, crossingSign, signedCrossingSign, checkSimple);
    checkCrossing(b, a, c, d, crossingSign, -signedCrossingSign, checkSimple);
    checkCrossing(a, b, d, c, crossingSign, -signedCrossingSign, checkSimple);
    checkCrossing(b, a, d, c, crossingSign, signedCrossingSign, checkSimple);
    checkCrossing(a, a, c, d, -1, 0, false);
    checkCrossing(a, b, c, c, -1, 0, false);
    checkCrossing(a, a, c, c, -1, 0, false);
    checkCrossing(a, b, a, b, 0, 1, false);
    if (crossingSign == 0) {
      // For vertex crossings, if AB crosses CD then CD does not cross AB. In order to get the
      // crossing sign right in both cases, all tests are specified such that AB crosses CD. The
      // other case is tested here.
      Preconditions.checkArgument(signedCrossingSign != 0);
      checkCrossing(c, d, a, b, crossingSign, 0, checkSimple);
    } else {
      checkCrossing(c, d, a, b, crossingSign, -signedCrossingSign, checkSimple);
    }
  }

  /**
   * The real tests of edge crossings are in loop and polygon unit tests, but we do a few simple
   * tests here.
   */
  @SuppressWarnings("FloatingPointLiteralPrecision") // Can't be helped for values like 1e-323.
  @Test
  public void testCrossings() {
    // 1. Two regular edges that cross.
    checkCrossings(
        new S2Point(1, 2, 1), new S2Point(1, -3, 0.5),
        new S2Point(1, -0.5, -3), new S2Point(0.1, 0.5, 3),
        1, 1, true);

    // 2. Two regular edges that intersect antipodal points.
    checkCrossings(
        new S2Point(1, 2, 1), new S2Point(1, -3, 0.5),
        new S2Point(-1, 0.5, 3), new S2Point(-0.1, -0.5, -3),
        -1, 0, true);

    // 3. Two edges on the same great circle that start at antipodal points.
    checkCrossings(
        new S2Point(0, 0, -1), new S2Point(0, 1, 0),
        new S2Point(0, 1, 1), new S2Point(0, 0, 1),
        -1, 0, true);

    // 4. Two edges that cross where one vertex is S2.origin().
    checkCrossings(
        new S2Point(1, 0, 0), S2.origin(),
        new S2Point(1, -0.1, 1), new S2Point(1, 1, -0.1),
        1, 1, true);

    // 5. Two edges that intersect antipodal points where one vertex is S2.origin().
    checkCrossings(
        new S2Point(1, 0, 0), S2.origin(),
        new S2Point(-1, 0.1, -1), new S2Point(-1, -1, 0.1),
        -1, 0, true);

    // 6. Two edges that share an endpoint. The Ortho() direction is (-4,0,2), and edge AB is
    // further CCW around (2,3,4) than CD.
    checkCrossings(
        new S2Point(7, -2, 3), new S2Point(2, 3, 4),
        new S2Point(2, 3, 4), new S2Point(-1, 2, 5),
        0, -1, true);

    // 7. Two edges that barely cross each other near the middle of one edge. The edge AB is
    // approximately in the x=y plane, while CD is approximately perpendicular to it and ends
    // exactly at the x=y plane.
    checkCrossings(
        new S2Point(1, 1, 1), new S2Point(1, nextAfter(1d, 0), -1),
        new S2Point(11, -12, -1), new S2Point(10, 10, 1),
        1, -1, false);

    // 8. In this version, the edges are separated by a distance of about 1e-15.
    checkCrossings(
        new S2Point(1, 1, 1), new S2Point(1, nextAfter(1d, 2), -1),
        new S2Point(1, -1, 0), new S2Point(1, 1, 0),
        -1, 0, false);

    // 9. Two edges that barely cross each other near the end of both edges. This example cannot be
    // handled using regular double-precision arithmetic due to floating-point underflow.
    checkCrossings(
        new S2Point(0, 0, 1), new S2Point(2, -1e-323, 1),
        new S2Point(1, -1, 1), new S2Point(1e-323, 0, 1),
        1, -1, false);

    // 10. In this version, the edges are separated by a distance of about 1e-640.
    checkCrossings(
        new S2Point(0, 0, 1), new S2Point(2, 1e-323, 1),
        new S2Point(1, -1, 1), new S2Point(1e-323, 0, 1),
        -1, 0, false);

    // 11. Two edges that barely cross each other near the middle of one edge. Computing the exact
    // determinant of some of the triangles in this test requires more than 2000 bits of precision.
    checkCrossings(
        new S2Point(1, -1e-323, -1e-323), new S2Point(1e-323, 1, 1e-323),
        new S2Point(1, -1, 1e-323), new S2Point(1, 1, 0),
        1, 1, false);

    // 12. In this version, the edges are separated by a distance of about 1e-640.
    checkCrossings(
        new S2Point(1, 1e-323, -1e-323), new S2Point(-1e-323, 1, 1e-323),
        new S2Point(1, -1, 1e-323), new S2Point(1, 1, 0),
        -1, 0, false);
  }

  private static S2LatLngRect getEdgeBound(S2Point a, S2Point b) {
    S2EdgeUtil.RectBounder bounder = new S2EdgeUtil.RectBounder();
    bounder.addPoint(a);
    bounder.addPoint(b);
    return bounder.getBound();
  }

  private static S2LatLngRect getEdgeBound(
      double x1, double y1, double z1, double x2, double y2, double z2) {
    return getEdgeBound(new S2Point(x1, y1, z1).normalize(), new S2Point(x2, y2, z2).normalize());
  }

  @Test
  public void testMaxLatitudeSimple() {
    // Check cases where the min/max latitude is attained at a vertex.
    final double kCubeLat = asin(1 / sqrt(3)); // 35.26 degrees
    assertTrue(
        getEdgeBound(1, 1, 1, 1, -1, -1)
            .approxEquals(
                new S2LatLngRect(
                    new R1Interval(-kCubeLat, kCubeLat), new S1Interval(-M_PI_4, M_PI_4)),
                RECT_ERROR));
    assertTrue(
        getEdgeBound(1, -1, 1, 1, 1, -1)
            .approxEquals(
                new S2LatLngRect(
                    new R1Interval(-kCubeLat, kCubeLat), new S1Interval(-M_PI_4, M_PI_4)),
                RECT_ERROR));

    // Check cases where the min/max latitude occurs in the edge interior. These tests expect the
    // result to be pretty close to the middle of the allowable error range (i.e., by adding 0.5 *
    // kRectError).

    // Max latitude, CW edge
    assertAlmostEquals(
        M_PI_4 + 0.5 * RECT_ERROR.lat().radians(), getEdgeBound(1, 1, 1, 1, -1, 1).lat().hi());
    // Max latitude, CCW edge
    assertAlmostEquals(
        M_PI_4 + 0.5 * RECT_ERROR.lat().radians(), getEdgeBound(1, -1, 1, 1, 1, 1).lat().hi());
    // Min latitude, CW edge
    assertAlmostEquals(
        -M_PI_4 - 0.5 * RECT_ERROR.lat().radians(), getEdgeBound(1, -1, -1, -1, -1, -1).lat().lo());
    // Min latitude, CCW edge
    assertAlmostEquals(
        -M_PI_4 - 0.5 * RECT_ERROR.lat().radians(), getEdgeBound(-1, 1, -1, -1, -1, -1).lat().lo());

    // Check cases where the edge passes through one of the poles.
    assertExactly(M_PI_2, getEdgeBound(.3, .4, 1, -.3, -.4, 1).lat().hi());
    assertExactly(-M_PI_2, getEdgeBound(.3, .4, -1, -.3, -.4, -1).lat().lo());
  }

  @Test
  public void testMaxLatitudeRandom() {
    // Check that the maximum latitude of edges is computed accurately to within 3 * DBL_EPSILON
    // (the expected maximum error). We concentrate on maximum latitudes near the equator and
    // north pole since these are the extremes.
    double latRadTolerance = RECT_ERROR.lat().radians();

    final int kIters = 100;
    for (int iter = 0; iter < kIters; ++iter) {
      // Construct a right-handed coordinate frame (U,V,W) such that U points slightly above the
      // equator, V points at the equator, and W is slightly offset from the north pole.
      S2Point u = data.getRandomPoint();
      // log is uniform
      u = new S2Point(u.x, u.y, S2.DBL_EPSILON * 1e-6 * pow(1e12, data.nextDouble())).normalize();
      S2Point v = robustCrossProd(new S2Point(0, 0, 1), u).normalize();
      S2Point w = robustCrossProd(u, v).normalize();

      // Construct a line segment AB that passes through U, and check that the maximum latitude of
      // this segment matches the latitude of U.
      S2Point a = u.sub(v.mul(data.nextDouble())).normalize();
      S2Point b = u.add(v.mul(data.nextDouble())).normalize();
      S2LatLngRect abBound = getEdgeBound(a, b);
      assertDoubleNear(S2LatLng.latitude(u).radians(), abBound.lat().hi(), latRadTolerance);

      // Construct a line segment CD that passes through W, and check that the maximum latitude of
      // this segment matches the latitude of W.
      S2Point c = w.sub(v.mul(data.nextDouble())).normalize();
      S2Point d = w.add(v.mul(data.nextDouble())).normalize();
      S2LatLngRect cdBound = getEdgeBound(c, d);
      assertDoubleNear(S2LatLng.latitude(w).radians(), cdBound.lat().hi(), latRadTolerance);
    }
  }

  private S2Point perturbATowardsB(S2Point a, S2Point b) {
    double choice = data.nextDouble();
    if (choice < 0.1) {
      return a;
    }
    if (choice < 0.3) {
      // Return a point that is exactly proportional to A and that still satisfies
      // S2.isUnitLength().
      for (; ; ) {
        S2Point temp = a.mul(2 - a.norm() + 5 * (data.nextDouble() - 0.5) * S2.DBL_EPSILON);
        if (!temp.equals(a) && S2.isUnitLength(temp)) {
          return temp;
        }
      }
    }
    if (choice < 0.5) {
      // Return a point such that the distance squared to A will underflow.
      return S2EdgeUtil.getPointOnLine(a, b, S1Angle.radians(1e-300));
    }
    // Otherwise return a point whose distance from A is near DBL_EPSILON such that the log of the
    // pdf is uniformly distributed.
    double distance = S2.DBL_EPSILON * 1e-5 * pow(1e6, data.nextDouble());
    return S2EdgeUtil.getPointOnLine(a, b, S1Angle.radians(distance));
  }

  private S2Point randomPole() {
    return new S2Point(0, 0, data.oneIn(2) ? 1 : -1);
  }

  private S2Point pointNearPole() {
    return perturbATowardsB(randomPole(), data.getRandomPoint());
  }

  private S2Point pointNearEquator() {
    return perturbATowardsB(
        new S2Point(data.nextDouble(), data.nextDouble(), 0).normalize(), randomPole());
  }

  @Test
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
    // Also test the corresponding situations for antipodal points, i.e. by negating one of the
    // points so that they are almost 180 degrees apart.
    final int kIters = 10000;
    for (int iter = 0; iter < kIters; ++iter) {
      S2Point a;
      S2Point b;
      switch (data.nextInt(5)) {
        case 0:
          // Two nearby points on a nearly-polar great circle.
          a = data.getRandomPoint();
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
          b = perturbATowardsB(a, data.getRandomPoint());
          break;
        case 3:
          // Two nearby points near the equator, but on any great circle.
          a = pointNearEquator();
          b = perturbATowardsB(a, data.getRandomPoint());
          break;
        case 4:
          // Two nearby points anywhere on the sphere.
          a = data.getRandomPoint();
          b = perturbATowardsB(a, data.getRandomPoint());
          break;
        default:
          throw new IllegalStateException("Can't happen.");
      }

      // The two points are chosen to be so close to each other that the min/max latitudes are
      // nearly always achieved at the edge endpoints. The only thing we need to watch out for is
      // that the latitude error bound is slightly larger if the min/max latitude occurs in the edge
      // interior.
      S2LatLngRect expectedBound = S2LatLngRect.fromPointPair(new S2LatLng(a), new S2LatLng(b));
      S2LatLngRect bound = getEdgeBound(a, b);
      assertTrue(bound.contains(expectedBound));
      assertTrue(expectedBound.expanded(RECT_ERROR).polarClosure().contains(bound));

      // If the two points are close enough and one point is negated (antipodal points), the bound
      // should be the entire sphere.
      if (a.sub(b).crossProd(a.add(b)).norm() <= 6.110 * S2.DBL_EPSILON) {
        assertEquals(S2LatLngRect.full(), getEdgeBound(a, b.neg()));
      }
    }
  }

  private S2LatLngRect getSubregionBound(double xLat, double xLng, double yYat, double yLng) {
    S2LatLngRect in =
        S2LatLngRect.fromPointPair(
            S2LatLng.fromRadians(xLat, xLng), S2LatLng.fromRadians(yYat, yLng));
    S2LatLngRect out = S2EdgeUtil.RectBounder.expandForSubregions(in);

    // Test that the bound is actually expanded.
    assertTrue(out.contains(in));
    if (in.lat().equals(S2LatLngRect.fullLat())) {
      assertFalse(in.lat().contains(out.lat()));
    }

    return out;
  }

  @Test
  public void testExpandForSubregions() {
    // First we check the various situations where the bound contains nearly-antipodal points. The
    // tests are organized into pairs where the two bounds are similar except that the first bound
    // meets the nearly-antipodal criteria while the second does not.

    // Cases where the bound does not straddle the equator (but almost does), and spans nearly 180
    // degrees in longitude.
    assertTrue(getSubregionBound(3e-16, 0, 1e-14, PI).isFull());
    assertFalse(getSubregionBound(9e-16, 0, 1e-14, PI).isFull());
    assertTrue(getSubregionBound(1e-16, 7e-16, 1e-14, PI).isFull());
    assertFalse(getSubregionBound(3e-16, 14e-16, 1e-14, PI).isFull());
    assertTrue(getSubregionBound(1e-100, 14e-16, 1e-14, PI).isFull());
    assertFalse(getSubregionBound(1e-100, 22e-16, 1e-14, PI).isFull());

    // Cases where the bound spans at most 90 degrees in longitude, and almost 180 degrees in
    // latitude. Note that DBL_EPSILON is about 2.22e-16, which implies that the double-precision
    // value just below Pi/2 can be written as (M_PI_2 - 2e-16).
    assertTrue(getSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 7e-16, 0).isFull());
    assertFalse(getSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 30e-16, 0).isFull());
    assertTrue(getSubregionBound(-M_PI_2 + 4e-16, 0, M_PI_2 - 2e-16, 1e-7).isFull());
    assertFalse(getSubregionBound(-M_PI_2 + 30e-16, 0, M_PI_2, 1e-7).isFull());
    assertTrue(getSubregionBound(-M_PI_2 + 4e-16, 0, M_PI_2 - 4e-16, M_PI_2).isFull());
    assertFalse(getSubregionBound(-M_PI_2, 0, M_PI_2 - 30e-16, M_PI_2).isFull());

    // Cases where the bound straddles the equator and spans more than 90 degrees in longitude.
    // These are the cases where the critical distance is between a corner of the bound and the
    // opposite longitudinal edge. Unlike the cases above, here the bound may contain
    // nearly-antipodal points (to within 3.055 * DBL_EPSILON) even though the latitude and
    // longitude ranges are both significantly less than (Pi - 3.055 * DBL_EPSILON).
    assertTrue(getSubregionBound(-M_PI_2, 0, M_PI_2 - 1e-8, PI - 1e-7).isFull());
    assertFalse(getSubregionBound(-M_PI_2, 0, M_PI_2 - 1e-7, PI - 1e-7).isFull());
    assertTrue(getSubregionBound(-M_PI_2 + 1e-12, -PI + 1e-4, M_PI_2, 0).isFull());
    assertTrue(getSubregionBound(-M_PI_2 + 1e-11, -PI + 1e-4, M_PI_2, 0).isFull());

    // Now we test cases where the bound does not contain nearly-antipodal points, but it does
    // contain points that are approximately 180 degrees apart in latitude.
    assertTrue(
        getSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 2e-16)
            .approxEquals(
                new S2LatLngRect(new R1Interval(1.5, 1.5), S1Interval.full()), RECT_ERROR));
    assertTrue(
        getSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 7e-16)
            .approxEquals(
                new S2LatLngRect(new R1Interval(1.5, 1.5), new S1Interval(-M_PI_2, M_PI_2 - 7e-16)),
                RECT_ERROR));

    // Test the full and empty bounds.
    assertTrue(S2EdgeUtil.RectBounder.expandForSubregions(S2LatLngRect.full()).isFull());
    assertTrue(S2EdgeUtil.RectBounder.expandForSubregions(S2LatLngRect.empty()).isEmpty());

    // Check for cases where the bound is expanded to include one of the poles.
    assertTrue(
        getSubregionBound(-M_PI_2 + 1e-15, 0, -M_PI_2 + 1e-15, 0)
            .approxEquals(
                new S2LatLngRect(new R1Interval(-M_PI_2, -M_PI_2 + 1e-15), S1Interval.full()),
                RECT_ERROR));
    assertTrue(
        getSubregionBound(M_PI_2 - 1e-15, 0, M_PI_2 - 1e-15, 0)
            .approxEquals(
                new S2LatLngRect(new R1Interval(M_PI_2 - 1e-15, M_PI_2), S1Interval.full()),
                RECT_ERROR));
  }

  @Test
  public void testLongitudePruner() {
    LongitudePruner pruner1 =
        new LongitudePruner(new S1Interval(0.75 * PI, -0.75 * PI), new S2Point(0, 1, 2));
    assertFalse(pruner1.intersects(new S2Point(1, 1, 3)));
    assertTrue(pruner1.intersects(new S2Point(-1 - 1e-15, -1, 0)));
    assertTrue(pruner1.intersects(new S2Point(-1, 0, 0)));
    assertTrue(pruner1.intersects(new S2Point(-1, 0, 0)));
    assertTrue(pruner1.intersects(new S2Point(1, -1, 8)));
    assertFalse(pruner1.intersects(new S2Point(1, 0, -2)));
    assertTrue(pruner1.intersects(new S2Point(-1, -1e-15, 0)));

    LongitudePruner pruner2 =
        new LongitudePruner(new S1Interval(0.25 * PI, 0.25 * PI), new S2Point(1, 0, 0));
    assertFalse(pruner2.intersects(new S2Point(2, 1, 2)));
    assertTrue(pruner2.intersects(new S2Point(1, 2, 3)));
    assertFalse(pruner2.intersects(new S2Point(0, 1, 4)));
    assertFalse(pruner2.intersects(new S2Point(-1e-15, -1, -1)));
  }

  private void checkWedge(
      S2Point a0,
      S2Point ab1,
      S2Point a2,
      S2Point b0,
      S2Point b2,
      boolean contains,
      boolean intersects,
      S2EdgeUtil.WedgeRelation wedgeRelation) {
    a0 = a0.normalize();
    ab1 = ab1.normalize();
    a2 = a2.normalize();
    b0 = b0.normalize();
    b2 = b2.normalize();
    assertEquals(contains, new S2EdgeUtil.WedgeContains().test(a0, ab1, a2, b0, b2) == 1);
    assertEquals(intersects, new S2EdgeUtil.WedgeIntersects().test(a0, ab1, a2, b0, b2) == -1);
    assertEquals(wedgeRelation, S2EdgeUtil.getWedgeRelation(a0, ab1, a2, b0, b2));
  }

  @Test
  public void testWedges() {
    // For simplicity, all of these tests use an origin of (0, 0, 1). This shouldn't matter as long
    // as the lower-level primitives are implemented correctly.

    // Intersection in one wedge.
    checkWedge(
        new S2Point(-1, 0, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, 2, 10),
        new S2Point(0, 1, 10),
        new S2Point(1, -2, 10),
        false,
        true,
        WedgeRelation.WEDGE_PROPERLY_OVERLAPS);
    // Intersection in two wedges.
    checkWedge(
        new S2Point(-1, -1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -1, 10),
        new S2Point(1, 0, 10),
        new S2Point(-1, 1, 10),
        false,
        true,
        WedgeRelation.WEDGE_PROPERLY_OVERLAPS);

    // Normal containment.
    checkWedge(
        new S2Point(-1, -1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -1, 10),
        new S2Point(-1, 0, 10),
        new S2Point(1, 0, 10),
        true,
        true,
        WedgeRelation.WEDGE_PROPERLY_CONTAINS);
    // Containment with equality on one side.
    checkWedge(
        new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(-1, -1, 10),
        new S2Point(2, 1, 10),
        new S2Point(1, -5, 10),
        true,
        true,
        WedgeRelation.WEDGE_PROPERLY_CONTAINS);
    // Containment with equality on the other side.
    checkWedge(
        new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(-1, -1, 10),
        new S2Point(1, -2, 10),
        new S2Point(-1, -1, 10),
        true,
        true,
        WedgeRelation.WEDGE_PROPERLY_CONTAINS);

    // Containment with equality on both sides.
    checkWedge(
        new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10),
        new S2Point(4, -5, 10),
        true,
        true,
        WedgeRelation.WEDGE_EQUALS);

    // Disjoint with equality on one side.
    checkWedge(
        new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(4, -5, 10),
        new S2Point(4, -5, 10),
        new S2Point(-2, -3, 10),
        false,
        false,
        WedgeRelation.WEDGE_IS_DISJOINT);
    // Disjoint with equality on the other side.
    checkWedge(
        new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(0, 5, 10),
        new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10),
        false,
        false,
        WedgeRelation.WEDGE_IS_DISJOINT);
    // Disjoint with equality on both sides.
    checkWedge(
        new S2Point(-2, 3, 10),
        new S2Point(0, 0, 1),
        new S2Point(4, -5, 10),
        new S2Point(4, -5, 10),
        new S2Point(-2, 3, 10),
        false,
        false,
        WedgeRelation.WEDGE_IS_DISJOINT);

    // B contains A with equality on one side.
    checkWedge(
        new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -5, 10),
        new S2Point(2, 1, 10),
        new S2Point(-1, -1, 10),
        false,
        true,
        WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
    // B contains A with equality on the other side.
    checkWedge(
        new S2Point(2, 1, 10),
        new S2Point(0, 0, 1),
        new S2Point(1, -5, 10),
        new S2Point(-2, 1, 10),
        new S2Point(1, -5, 10),
        false,
        true,
        WedgeRelation.WEDGE_IS_PROPERLY_CONTAINED);
  }

  /**
   * Given a point X and an edge AB, check that the distance from X to AB is "distanceRadians" and
   * the closest point on AB is "expectedClosest".
   */
  private void checkDistance(
      S2Point x, S2Point a, S2Point b, double distanceRadians, S2Point expectedClosest) {
    x = x.normalize();
    a = a.normalize();
    b = b.normalize();
    expectedClosest = expectedClosest.normalize();
    // Note the 1e-15 epsilon is not needed by the C++ library, but we need it here.
    // TODO(torrey): Investigate why this is needed, fix if possible.
    assertEquals(distanceRadians, S2EdgeUtil.getDistance(x, a, b).radians(), 1E-15);
    S2Point closest = S2EdgeUtil.getClosestPoint(x, a, b);
    if (expectedClosest.equals(ORIGIN)) {
      // This special value says that the result should be A or B.
      assertTrue(closest.equals(a) || closest.equals(b));
    } else {
      assertTrue(S2.approxEquals(closest, expectedClosest));
    }
    assertEquals(S1ChordAngle.ZERO, S2EdgeUtil.updateMinDistance(x, a, b, S1ChordAngle.ZERO));
    S1ChordAngle minDistance = S2EdgeUtil.updateMinDistance(x, a, b, S1ChordAngle.INFINITY);
    assertEquals(distanceRadians, minDistance.toAngle().radians(), 1e-15);
  }

  @Test
  public void testDistance() {
    checkDistance(
        new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0), 0, new S2Point(1, 0, 0));
    checkDistance(
        new S2Point(0, 1, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0), 0, new S2Point(0, 1, 0));
    checkDistance(
        new S2Point(1, 3, 0), new S2Point(1, 0, 0), new S2Point(0, 1, 0), 0, new S2Point(1, 3, 0));
    checkDistance(
        new S2Point(0, 0, 1),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        M_PI_2,
        new S2Point(1, 0, 0));
    checkDistance(
        new S2Point(0, 0, -1),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        M_PI_2,
        new S2Point(1, 0, 0));
    checkDistance(
        new S2Point(-1, -1, 0),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        0.75 * PI,
        new S2Point(0, 0, 0));

    checkDistance(
        new S2Point(0, 1, 0),
        new S2Point(1, 0, 0),
        new S2Point(1, 1, 0),
        M_PI_4,
        new S2Point(1, 1, 0));
    checkDistance(
        new S2Point(0, -1, 0),
        new S2Point(1, 0, 0),
        new S2Point(1, 1, 0),
        M_PI_2,
        new S2Point(1, 0, 0));

    checkDistance(
        new S2Point(0, -1, 0),
        new S2Point(1, 0, 0),
        new S2Point(-1, 1, 0),
        M_PI_2,
        new S2Point(1, 0, 0));
    checkDistance(
        new S2Point(-1, -1, 0),
        new S2Point(1, 0, 0),
        new S2Point(-1, 1, 0),
        M_PI_2,
        new S2Point(-1, 1, 0));

    checkDistance(
        new S2Point(1, 1, 1),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        asin(sqrt(1. / 3)),
        new S2Point(1, 1, 0));
    checkDistance(
        new S2Point(1, 1, -1),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        asin(sqrt(1. / 3)),
        new S2Point(1, 1, 0));

    checkDistance(
        new S2Point(-1, 0, 0),
        new S2Point(1, 1, 0),
        new S2Point(1, 1, 0),
        0.75 * PI,
        new S2Point(1, 1, 0));
    checkDistance(
        new S2Point(0, 0, -1),
        new S2Point(1, 1, 0),
        new S2Point(1, 1, 0),
        M_PI_2,
        new S2Point(1, 1, 0));
    checkDistance(
        new S2Point(-1, 0, 0),
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0),
        PI,
        new S2Point(1, 0, 0));
  }

  private static void checkMaxDistance(S2Point x, S2Point a, S2Point b, double distanceRadians) {
    x = x.normalize();
    a = a.normalize();
    b = b.normalize();

    assertEquals(
        S2EdgeUtil.updateMaxDistance(x, a, b, S1ChordAngle.STRAIGHT), S1ChordAngle.STRAIGHT);

    S1ChordAngle maxDistance = S2EdgeUtil.updateMaxDistance(x, a, b, S1ChordAngle.NEGATIVE);
    assertEquals(distanceRadians, maxDistance.toAngle().radians(), 1e-15);
  }

  @Test
  public void testUpdateMaxDistance() {
    checkMaxDistance(new S2Point(1, 0, 1), X_POS, Y_POS, M_PI_2);
    checkMaxDistance(new S2Point(1, 0, -1), X_POS, Y_POS, M_PI_2);
    checkMaxDistance(new S2Point(0, 1, 1), X_POS, Y_POS, M_PI_2);
    checkMaxDistance(new S2Point(0, 1, -1), X_POS, Y_POS, M_PI_2);

    double expectedRadians = asin(sqrt(2. / 3));
    checkMaxDistance(new S2Point(1, 1, 1), X_POS, Y_POS, expectedRadians);
    checkMaxDistance(new S2Point(1, 1, -1), X_POS, Y_POS, expectedRadians);

    checkMaxDistance(X_POS, new S2Point(1, 1, 0), new S2Point(1, -1, 0), M_PI_4);
    checkMaxDistance(Y_POS, new S2Point(1, 1, 0), new S2Point(-1, 1, 0), M_PI_4);
    checkMaxDistance(Z_POS, new S2Point(0, 1, 1), new S2Point(0, -1, 1), M_PI_4);

    checkMaxDistance(Z_POS, X_POS, new S2Point(1, 0, -1), 3 * M_PI_4);
    checkMaxDistance(Z_POS, X_POS, new S2Point(1, 1, -sqrt(2)), 3 * M_PI_4);

    checkMaxDistance(Z_POS, Z_NEG, Z_NEG, PI);
  }

  private void checkInterpolate(double t, S2Point a, S2Point b, S2Point expected) {
    a = a.normalize();
    b = b.normalize();
    expected = expected.normalize();

    // We allow a bit more than the usual 1e-15 error tolerance because interpolate() uses trig
    // functions.
    final double kErrorRadians = 3e-15;
    assertApproxEquals(expected, S2EdgeUtil.interpolate(t, a, b), kErrorRadians);

    // Now test the other interpolation functions.
    S1Angle r = new S1Angle(a, b).mul(t);
    assertLessOrEqual(
        new S1Angle(S2EdgeUtil.getPointOnLine(a, b, r), expected).radians(), kErrorRadians);
    if (a.dotProd(b) == 0) {  // Common in the test cases below.
      assertLessOrEqual(
          new S1Angle(S2EdgeUtil.getPointOnRay(a, b, r), expected).radians(), kErrorRadians);
    }
    if (r.radians() >= 0 && r.radians() < 0.99 * PI) {
      S1ChordAngle rCA = S1ChordAngle.fromS1Angle(r);
      assertLessOrEqual(
          new S1Angle(S2EdgeUtil.getPointOnLine(a, b, rCA), expected).radians(), kErrorRadians);
      if (a.dotProd(b) == 0) {
        assertLessOrEqual(
            new S1Angle(S2EdgeUtil.getPointOnRay(a, b, rCA), expected).radians(), kErrorRadians);
      }
    }
  }

  @Test
  public void testInterpolate() {
    // A zero-length edge, getting the end points with 't' == 0 or 1.
    checkInterpolate(0, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));
    checkInterpolate(1, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));
    checkInterpolate(0.5, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));
    checkInterpolate(
        Double.MIN_VALUE, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));

    // Cases where A and B are the same but actually trying to interpolate between them.
    checkInterpolate(0.5, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));
    checkInterpolate(
        Double.MIN_VALUE, new S2Point(1, 0, 0), new S2Point(1, 0, 0), new S2Point(1, 0, 0));

    // Start, end, and middle of a medium-length edge.
    checkInterpolate(0, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(1, 0, 0));
    checkInterpolate(1, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(0, 1, 0));
    checkInterpolate(0.5, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(1, 1, 0));

    // Choose test points designed to expose floating-point errors.
    S2Point p1 = new S2Point(0.1, 1e-30, 0.3).normalize();
    S2Point p2 = new S2Point(-0.7, -0.55, -1e30).normalize();

    // Another zero-length edge.
    checkInterpolate(0, p1, p1, p1);
    checkInterpolate(1, p1, p1, p1);
    checkInterpolate(0.5, p1, p1, p1);
    checkInterpolate(Double.MIN_VALUE, p1, p1, p1);

    // Start, end, and middle of a medium-length edge.
    checkInterpolate(0, p1, p2, p1);
    checkInterpolate(1, p1, p2, p2);
    checkInterpolate(0.5, p1, p2, p1.add(p2).mul(0.5));

    // Test that interpolation is done using distances on the sphere rather than linear distances.
    checkInterpolate(
        1. / 3, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(sqrt(3), 1, 0));
    checkInterpolate(
        2. / 3, new S2Point(1, 0, 0), new S2Point(0, 1, 0), new S2Point(1, sqrt(3), 0));

    // Test that interpolation is accurate on a long edge (but not so long that the definition of
    // the edge itself becomes too unstable).
    final double kLng = PI - 1e-2;
    S2Point a = S2LatLng.fromRadians(0, 0).toPoint();
    S2Point b = S2LatLng.fromRadians(0, kLng).toPoint();
    for (double f = 0.4; f > 1e-15; f *= 0.1) {
      checkInterpolate(f, a, b, S2LatLng.fromRadians(0, f * kLng).toPoint());
      checkInterpolate(1 - f, a, b, S2LatLng.fromRadians(0, (1 - f) * kLng).toPoint());
    }

    // Test that interpolation on a 180 degree edge (antipodal endpoints) yields a result with the
    // correct distance from each endpoint.
    for (double t = 0; t <= 1; t += 0.125) {
      S2Point actual = S2EdgeUtil.interpolate(t, p1, p1.neg());
      assertDoubleNear(new S1Angle(actual, p1).radians(), t * PI, 3e-15);
    }
  }

  @Test
  public void testInterpolateCanExtrapolate() {
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

  @Test
  public void testRepeatedInterpolation() {
    // Check that points do not drift away from unit length when repeated interpolations are done.
    for (int i = 0; i < 100; ++i) {
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      for (int j = 0; j < 1000; ++j) {
        a = S2EdgeUtil.interpolate(0.01, a, b);
      }
      assertTrue(S2.isUnitLength(a));
    }
  }

  @Test
  public void testGetPointToLeftS1Angle() {
    S2Point a = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(0, 5).toPoint();  // east
    final S1Angle kDistance = metersToAngle(10);

    S2Point c = S2EdgeUtil.getPointToLeft(a, b, kDistance);
    assertDoubleNear(new S1Angle(a, c).radians(), kDistance.radians(), 1e-15);
    // CAB must be a right angle with C to the left of AB.
    assertDoubleNear(S2.turnAngle(c, a, b), M_PI_2 /*radians*/, 1e-15);
  }

  @Test
  public void testGetPointToLeftS1ChordAngle() {
    S2Point a = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(0, 5).toPoint();  // east
    final S1Angle kDistance = metersToAngle(10);

    S2Point c = S2EdgeUtil.getPointToLeft(a, b, S1ChordAngle.fromS1Angle(kDistance));
    assertDoubleNear(new S1Angle(a, c).radians(), kDistance.radians(), 1e-15);
    // CAB must be a right angle with C to the left of AB.
    assertDoubleNear(S2.turnAngle(c, a, b),  M_PI_2 /*radians*/, 1e-15);
  }

  @Test
  public void testGetPointToRightS1Angle() {
    S2Point a = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(0, 5).toPoint();  // east
    final S1Angle kDistance = metersToAngle(10);

    S2Point c = S2EdgeUtil.getPointToRight(a, b, kDistance);
    assertDoubleNear(new S1Angle(a, c).radians(), kDistance.radians(), 1e-15);
    // CAB must be a right angle with C to the right of AB.
    assertDoubleNear(S2.turnAngle(c, a, b),  -M_PI_2 /*radians*/, 1e-15);
  }

  @Test
  public void testGetPointToRightS1ChordAngle() {
    S2Point a = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(0, 5).toPoint();  // east
    final S1Angle kDistance = metersToAngle(10);

    S2Point c = S2EdgeUtil.getPointToRight(a, b, S1ChordAngle.fromS1Angle(kDistance));
    assertDoubleNear(new S1Angle(a, c).radians(), kDistance.radians(), 1e-15);
    // CAB must be a right angle with C to the right of AB.
    assertDoubleNear(S2.turnAngle(c, a, b), -M_PI_2 /*radians*/, 1e-15);
  }

  /** Used to denote which method was used to calculate an intersection. */
  private enum IntersectionMethod {
    /**
     * Denotes that {@link S2EdgeUtil#getIntersectionApprox(S2Point, S2Point, S2Point, S2Point,
     * ResultError)} was used to calculate an intersection.
     */
    APPROX,
    /**
     * Denotes that {@link S2EdgeUtil#getIntersectionExact(S2Point, S2Point, S2Point, S2Point)} was
     * used to calculate an intersection.
     */
    EXACT
  }

  @Test
  public void testIntersectionError() {
    // We repeatedly construct two edges that cross near a random point "p", and measure the
    // distance from the actual intersection point "x" to the expected intersection point "p" and
    // also to the edges that cross near "p".
    //
    // Note that GetIntersection() does not guarantee that "x" and "p" will be close together (since
    // the intersection point is numerically unstable when the edges cross at a very small angle),
    // but it does guarantee that "x" will be close to both of the edges that cross.
    Collection<IntersectionMethod> methodsUsed = HashMultiset.create();
    checkCrossingEdges(
        (S2Point a, S2Point b, S2Point c, S2Point d, S2Point p, double slope) -> {
          // Each constructed edge should be at most 1.5 * DBL_EPSILON away from the original point.
          S1Angle abDist = S2EdgeUtil.getDistance(p, a, b);
          S1Angle cdDist = S2EdgeUtil.getDistance(p, c, d);
          assertTrue(
              abDist.lessOrEquals(S1Angle.radians(1.5 * DBL_EPSILON + GET_DISTANCE_ABS_ERROR)));
          assertTrue(
              cdDist.lessOrEquals(S1Angle.radians(1.5 * DBL_EPSILON + GET_DISTANCE_ABS_ERROR)));

          // Verify that the expected intersection point is close to both edges and also close to
          // the original point P. (It might not be very close to P if the angle between the edges
          // is very small.)
          S2Point expected = getIntersectionExact(a, b, c, d);
          abDist = S2EdgeUtil.getDistance(expected, a, b);
          cdDist = S2EdgeUtil.getDistance(expected, c, d);
          S1Angle pDist = new S1Angle(expected, p);
          assertTrue(abDist.lessOrEquals(S1Angle.radians(2 * GET_DISTANCE_ABS_ERROR)));
          assertTrue(cdDist.lessOrEquals(S1Angle.radians(2 * GET_DISTANCE_ABS_ERROR)));
          assertTrue(
              pDist.lessOrEquals(
                  S1Angle.radians(GET_DISTANCE_ABS_ERROR / slope + INTERSECTION_ERROR)));

          // Now we actually test the getIntersection method.
          ResultError resultError = new ResultError();
          S2Point actual = S2EdgeUtil.getIntersection(a, b, c, d, resultError);
          methodsUsed.add(
              (resultError.error > INTERSECTION_ERROR)
                  ? IntersectionMethod.EXACT
                  : IntersectionMethod.APPROX);
          abDist = S2EdgeUtil.getDistance(actual, a, b);
          cdDist = S2EdgeUtil.getDistance(actual, c, d);
          pDist = new S1Angle(expected, actual);
          assertTrue(abDist.lessOrEquals(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
          assertTrue(cdDist.lessOrEquals(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE));
          assertTrue(
              pDist.lessOrEquals(S1Angle.radians(INTERSECTION_ERROR + GET_DISTANCE_ABS_ERROR)));
        });
    System.out.println(methodsUsed);
  }

  @Test
  public void testGetIntersectionApproxError() {
    checkCrossingEdges(
        (S2Point a, S2Point b, S2Point c, S2Point d, S2Point p, double slope) -> {
          ResultError resultError = new ResultError();
          S2Point expected = getIntersectionExact(a, b, c, d);
          S2Point actual = getIntersectionApprox(a, b, c, d, resultError);
          S1Angle dist = new S1Angle(expected, actual);
          // GET_DISTANCE_ABS_ERROR and DBL_EPSILON are added to account for S1Angle distance
          // calculation error and BigDecimal rounding error, respectively.
          assertTrue(
              dist.lessOrEquals(
                  S1Angle.radians(resultError.error + GET_DISTANCE_ABS_ERROR + DBL_EPSILON)));
        });
  }

  /** Returns {@link S2EdgeUtil#getIntersectionExact} result with sign corrected (if necessary). */
  private static S2Point getIntersectionExact(S2Point a0, S2Point a1, S2Point b0, S2Point b1) {
    return S2EdgeUtil.correctIntersectionSign(
        a0, a1, b0, b1, S2EdgeUtil.getIntersectionExact(a0, a1, b0, b1));
  }

  /** Returns {@link S2EdgeUtil#getIntersectionApprox} result with sign corrected (if necessary). */
  private static S2Point getIntersectionApprox(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, ResultError resultError) {
    return S2EdgeUtil.correctIntersectionSign(
        a0, a1, b0, b1, S2EdgeUtil.getIntersectionApprox(a0, a1, b0, b1, resultError));
  }

  @Test
  public void testGetProjectionError() {
    checkCrossingEdges(
        (S2Point a, S2Point b, S2Point c, S2Point d, S2Point p, double slope) -> {
          // Precalculate abNorm and its length.
          S2Point abNorm = robustCrossProd(a, b);
          double abNormLen = abNorm.norm();

          // Test projection of c on (a, b). Will be twice dot(c, cross(a, b)) given how
          // S2EdgeUtil.getProjection() is implemented.
          ResultError resultError = new ResultError();
          double expected = new BigDecimal(2).multiply(getDotCrossExact(c, a, b)).doubleValue();
          double actual = S2EdgeUtil.getProjection(c, abNorm, abNormLen, a, b, resultError);
          // DBL_EPSILON is added to result error due to rounding error of actual from BigDecimal.
          assertDoubleNear(expected, actual, resultError.error + DBL_EPSILON);

          // Test projection of d on (a, b). Will be twice dot(d, cross(a, b)) given how
          // S2EdgeUtil.getProjection() is implemented.
          expected = new BigDecimal(2).multiply(getDotCrossExact(d, a, b)).doubleValue();
          actual = S2EdgeUtil.getProjection(d, abNorm, abNormLen, a, b, resultError);
          // DBL_EPSILON is added to result error due to rounding error of actual from BigDecimal.
          assertDoubleNear(expected, actual, resultError.error + DBL_EPSILON);
        });
  }

  /** Return the dot product between x and the cross product of (a, b) using exact arithmetic. */
  private static BigDecimal getDotCrossExact(S2Point x, S2Point a, S2Point b) {
    return new BigPoint(x).dotProd(new BigPoint(a).crossProd(new BigPoint(b)));
  }

  /**
   * Given two edges a0a1 and b0b1, check that the minimum distance between them is {@code
   * distanceRadians}, and that {@link S2EdgeUtil#getEdgePairClosestPoints} returns {@code
   * expectedA} and {@code expectedB} as the points that achieve this distance. {@link
   * S2Point#ORIGIN} may be passed for {@code expectedA} or {@code expectedB} to indicate that both
   * endpoints of the corresponding edge are equally distant, and therefore either one might be
   * returned. Parameters are passed by value so that this function can normalize them.
   */
  private void checkEdgePairDistance(
      S2Point a0,
      S2Point a1,
      S2Point b0,
      S2Point b1,
      double distanceRadians,
      S2Point expectedA,
      S2Point expectedB) {
    a0 = a0.normalize();
    a1 = a1.normalize();
    b0 = b0.normalize();
    b1 = b1.normalize();
    expectedA = expectedA.normalize();
    expectedB = expectedB.normalize();
    S2Shape.MutableEdge result = new S2Shape.MutableEdge();
    S2EdgeUtil.getEdgePairClosestPoints(a0, a1, b0, b1, result);
    if (result.a.equalsPoint(ORIGIN)) {
      // This special value says that the result should be a0 or a1.
      assertTrue(result.a.equalsPoint(a0) || result.a.equalsPoint(a1));
    } else {
      assertTrue(S2.approxEquals(expectedA, result.a));
    }
    if (expectedB.equalsPoint(ORIGIN)) {
      // This special value says that the result should be b0 or b1.
      assertTrue(result.b.equalsPoint(b0) || result.b.equalsPoint(b1));
    } else {
      assertTrue(S2.approxEquals(expectedB, result.b));
    }
    assertEquals(
        S1ChordAngle.ZERO, S2EdgeUtil.getEdgePairMinDistance(a0, a1, b0, b1, S1ChordAngle.ZERO));
    assertNotSame(
        S1ChordAngle.INFINITY,
        S2EdgeUtil.getEdgePairMinDistance(a0, a1, b0, b1, S1ChordAngle.INFINITY));
    S1ChordAngle actual = S2EdgeUtil.getEdgePairMinDistance(a0, a1, b0, b1, S1ChordAngle.INFINITY);
    assertDoubleNear(distanceRadians, actual.toAngle().radians(), 1e-15);
  }

  @Test
  public void testEdgePairMinDistance() {
    // One edge is degenerate.
    checkEdgePairDistance(
        new S2Point(1, 0, 1),
        new S2Point(1, 0, 1),
        new S2Point(1, -1, 0),
        new S2Point(1, 1, 0),
        M_PI_4,
        new S2Point(1, 0, 1),
        new S2Point(1, 0, 0));

    checkEdgePairDistance(
        new S2Point(1, -1, 0),
        new S2Point(1, 1, 0),
        new S2Point(1, 0, 1),
        new S2Point(1, 0, 1),
        M_PI_4,
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 1));

    // Both edges are degenerate.
    checkEdgePairDistance(
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 0),
        M_PI_2,
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0));

    // Both edges are degenerate and antipodal.
    checkEdgePairDistance(
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0),
        new S2Point(-1, 0, 0),
        new S2Point(-1, 0, 0),
        PI,
        new S2Point(1, 0, 0),
        new S2Point(-1, 0, 0));

    // Two identical edges.
    checkEdgePairDistance(
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        0,
        new S2Point(0, 0, 0),
        new S2Point(0, 0, 0));

    // Both edges are degenerate and identical.
    checkEdgePairDistance(
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0),
        0,
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0));

    // Edges that share exactly one vertex (all 4 possibilities).
    checkEdgePairDistance(
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 1),
        0,
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 0));

    checkEdgePairDistance(
        new S2Point(0, 1, 0),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 1),
        0,
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 0));

    checkEdgePairDistance(
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 1),
        new S2Point(0, 1, 0),
        0,
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 0));

    checkEdgePairDistance(
        new S2Point(0, 1, 0),
        new S2Point(1, 0, 0),
        new S2Point(0, 1, 1),
        new S2Point(0, 1, 0),
        0,
        new S2Point(0, 1, 0),
        new S2Point(0, 1, 0));

    // Two edges whose interiors cross.
    checkEdgePairDistance(
        new S2Point(1, -1, 0),
        new S2Point(1, 1, 0),
        new S2Point(1, 0, -1),
        new S2Point(1, 0, 1),
        0,
        new S2Point(1, 0, 0),
        new S2Point(1, 0, 0));

    // The closest distance occurs between two edge endpoints, but more than one endpoint pair is
    // equally distant.
    checkEdgePairDistance(
        new S2Point(1, -1, 0),
        new S2Point(1, 1, 0),
        new S2Point(-1, 0, 0),
        new S2Point(-1, 0, 1),
        acos(-0.5),
        new S2Point(0, 0, 0),
        new S2Point(-1, 0, 1));

    checkEdgePairDistance(
        new S2Point(-1, 0, 0),
        new S2Point(-1, 0, 1),
        new S2Point(1, -1, 0),
        new S2Point(1, 1, 0),
        acos(-0.5),
        new S2Point(-1, 0, 1),
        new S2Point(0, 0, 0));

    checkEdgePairDistance(
        new S2Point(1, -1, 0),
        new S2Point(1, 1, 0),
        new S2Point(-1, 0, -1),
        new S2Point(-1, 0, 1),
        acos(-0.5),
        new S2Point(0, 0, 0),
        new S2Point(0, 0, 0));
  }

  /**
   * Given two edges a0,a1 and b0,b1, check that the maximum distance between them is within 1e-15
   * of "distanceRadians". Parameters are normalized before checking them.
   */
  private static void checkEdgePairMaxDistance(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, double distanceRadians) {
    a0 = a0.normalize();
    a1 = a1.normalize();
    b0 = b0.normalize();
    b1 = b1.normalize();

    S1ChordAngle maxed = S2EdgeUtil.getEdgePairMaxDistance(a0, a1, b0, b1, S1ChordAngle.STRAIGHT);
    assertSame(maxed, S1ChordAngle.STRAIGHT);
    S1ChordAngle value = S2EdgeUtil.getEdgePairMaxDistance(a0, a1, b0, b1, S1ChordAngle.NEGATIVE);
    assertEquals(distanceRadians, value.toAngle().radians(), 1e-15);
  }

  /** Simple tests for edge pair max distances. */
  @Test
  public void testEdgePairMaxDistanceSimple() {
    // Some test points.
    S2Point center = ll(12, 12);
    S2Point north = ll(13, 12);
    S2Point south = ll(11, 12);
    S2Point east = ll(12, 13);
    S2Point west = ll(12, 11);
    S2Point nw = ll(13, 11);
    S2Point sw = ll(11, 11);
    S2Point ss = ll(10, 12); // Further south

    // Edges with a common endpoint, where the maximum distance is between the other endpoints:
    // Straight line north-south
    checkEdgePairMaxDistance(north, center, center, south, north.angle(south));
    // Straight line east-west
    checkEdgePairMaxDistance(east, center, center, west, east.angle(west));
    // 90 degree turn north to west
    checkEdgePairMaxDistance(north, center, center, west, north.angle(west));

    // Lines that cross:
    // 90 degree crossing at 'center'
    checkEdgePairMaxDistance(north, south, east, west, south.angle(west));
    // Acute angle crossing at the midpoint between "west" and "center".
    checkEdgePairMaxDistance(north, sw, south, nw, south.angle(north));

    // Lines that don't cross. Max distance is between the most distant endpoints.
    checkEdgePairMaxDistance(north, nw, south, ss, ss.angle(nw));
  }

  @Test
  public void testEdgePairMaxDistance() {
    // Standard situation. Same hemisphere, not degenerate.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        new S2Point(1, 1, 0), new S2Point(1, 1, 1),
        acos(1 / sqrt(3)));

    // One edge is degenerate.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 1), new S2Point(1, 0, 1),
        new S2Point(1, -1, 0), new S2Point(1, 1, 0),
        acos(0.5));
    checkEdgePairMaxDistance(
        new S2Point(1, -1, 0), new S2Point(1, 1, 0),
        new S2Point(1, 0, 1), new S2Point(1, 0, 1),
        acos(0.5));

    // Both edges are degenerate.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 0), new S2Point(1, 0, 0),
        new S2Point(0, 1, 0), new S2Point(0, 1, 0),
        M_PI_2);

    // Both edges are degenerate and antipodal.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 0), new S2Point(1, 0, 0),
        new S2Point(-1, 0, 0), new S2Point(-1, 0, 0),
        PI);

    // Two identical edges.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        new S2Point(1, 0, 0), new S2Point(0, 1, 0),
        M_PI_2);

    // Both edges are degenerate and identical.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 0), new S2Point(1, 0, 0),
        new S2Point(1, 0, 0), new S2Point(1, 0, 0),
        0);

    // Antipodal reflection of one edge crosses the other edge.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 1),  new S2Point(1, 0, -1),
        new S2Point(-1, -1, 0), new S2Point(-1, 1, 0),
        PI);

    // One vertex of one edge touches the interior of the antipodal reflection of the other edge.
    checkEdgePairMaxDistance(
        new S2Point(1, 0, 1), new S2Point(1, 0, 0),
        new S2Point(-1, -1, 0), new S2Point(-1, 1, 0),
        PI);
  }


  private boolean isEdgeBNearEdgeA(String aStr, String bStr, double maxErrorDegrees) {
    S2Polyline a = makePolyline(aStr);
    assertEquals(2, a.numVertices());
    S2Polyline b = makePolyline(bStr);
    assertEquals(2, b.numVertices());
    return S2EdgeUtil.isEdgeBNearEdgeA(
        a.vertex(0), a.vertex(1), b.vertex(0), b.vertex(1), S1Angle.degrees(maxErrorDegrees));
  }

  @Test
  public void testEdgeBNearEdgeA() {
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

    // Very long edges whose endpoints are close may have interior points that are far apart. An
    // implementation that only considers the vertices of polylines will incorrectly consider such
    // edges as "close" when they are not. Consider, for example, two consecutive lines of
    // longitude. As they approach the poles, they become arbitrarily close together, but along the
    // equator they bow apart.
    assertFalse(isEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 0.5));
    assertTrue(isEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 1.5));

    // The two arcs here are nearly as long as S2 edges can be (just shy of 180 degrees), and their
    // endpoints are less than 1 degree apart. Their midpoints, however, are at opposite ends of the
    // sphere along its equator.
    assertFalse(isEdgeBNearEdgeA("0:-179.75, 0:-0.25", "0:179.75, 0:0.25", 1.0));

    // At the equator, the second arc here is 9.75 degrees from the first, and closer at all other
    // points. However, the southern point of the second arc (-1, 9.75) is too far from the first
    // arc for the short-circuiting logic in isEdgeBNearEdgeA to apply.
    assertTrue(isEdgeBNearEdgeA("40:0, -5:0", "39:0.975, -1:0.975", 1.0));

    // Same as above, but B's orientation is reversed, causing the angle between the normal vectors
    // of circ(B) and circ(A) to be (180-9.75) = 170.5 degrees, not 9.75 degrees. The greatest
    // separation between the planes is still 9.75 degrees.
    assertTrue(isEdgeBNearEdgeA("10:0, -10:0", "-.4:0.975, 0.4:0.975", 1.0));

    // A and B are on the same great circle, A and B partially overlap, but the only part of B that
    // does not overlap A is shorter than tolerance.
    assertTrue(isEdgeBNearEdgeA("0:0, 1:0", "0.9:0, 1.1:0", 0.25));

    // A and B are on the same great circle, all points on B are close to A at its second endpoint,
    // (1,0).
    assertTrue(isEdgeBNearEdgeA("0:0, 1:0", "1.1:0, 1.2:0", 0.25));

    // Same as above, but B's orientation is reversed. This case is special because the projection
    // of the normal defining A onto the plane containing B is the null vector, and must be handled
    // by a special case.
    assertTrue(isEdgeBNearEdgeA("0:0, 1:0", "1.2:0, 1.1:0", 0.25));
  }

  @Test
  public void testCollinearEdgesThatDontTouch() {
    final int kIters = 500;
    for (int iter = 0; iter < kIters; ++iter) {
      S2Point a = data.getRandomPoint();
      S2Point d = data.getRandomPoint();
      S2Point b = S2EdgeUtil.interpolate(0.05, a, d);
      S2Point c = S2EdgeUtil.interpolate(0.95, a, d);
      assertTrue(S2EdgeUtil.robustCrossing(a, b, c, d) < 0);
      EdgeCrosser crosser = new EdgeCrosser(a, b, c);
      assertTrue(crosser.robustCrossing(d) < 0);
      assertTrue(crosser.robustCrossing(c) < 0);
    }
  }

  // Returns 0 or a power of 2, based on a randomly selected value drawn from a skewed distribution.
  private double randomSkewedCoord() {
    int binaryExp = data.skewed(11);
    return (binaryExp > 1022) ? 0 : pow(2, -binaryExp);
  }

  @GwtIncompatible("fails due to inadequate robustness in J2CL's S2Predicates.sign method")
  @Test
  public void testCoincidentZeroLengthEdgesThatDontTouch() {
    // It is important that the edge primitives can handle vertices that are exactly proportional to
    // each other, i.e. that are not identical but are nevertheless exactly coincident when
    // projected onto the unit sphere. There are various ways that such points can arise. For
    // example, normalize() itself is not idempotent: there exist distinct points A,B such that
    // normalize(A) == B  and normalize(B) == A. Another issue is that sometimes calls to
    // normalize() are skipped when the result of a calculation "should" be unit length
    // mathematically (e.g., when computing the cross product of two orthonormal vectors).
    //
    // This test checks pairs of edges AB and CD where A,B,C,D are exactly coincident on the sphere
    // and the norms of A,B,C,D are monotonically increasing. Such edge pairs should never
    // intersect. (This is not obvious, since it depends on the particular symbolic perturbations
    // used by sign(). It would be better to replace this with a test that says that the CCW results
    // must be consistent with each other.)
    final int kIters = 1000;
    for (int iter = 0; iter < kIters; ++iter) {
      // Construct a point P where every component is zero or a power of 2.
      S2Point p = new S2Point(randomSkewedCoord(), randomSkewedCoord(), randomSkewedCoord());

      // If all components were zero, try again. Note that normalization may convert a non-zero
      // point into a zero one due to underflow (!)
      p = p.normalize();
      if (p.equals(ORIGIN)) {
        --iter;
        continue;
      }

      // Now every non-zero component should have exactly the same mantissa. This implies that if we
      // scale the point by an arbitrary factor, every non-zero component will still have the same
      // mantissa. Scale the points so that they are all distinct and are still very likely to
      // satisfy S2.isUnitLength (which allows for a small amount of error in the norm).
      S2Point a = p.mul(1 - 3e-16);
      S2Point b = p.mul(1 - 1e-16);
      S2Point c = p;
      S2Point d = p.mul(1 + 2e-16);
      if (!S2.isUnitLength(a) || !S2.isUnitLength(d)) {
        --iter;
        continue;
      }

      // Verify that the expected edges do not cross.
      assertTrue(S2EdgeUtil.robustCrossing(a, b, c, d) < 0);
      EdgeCrosser crosser = new EdgeCrosser(a, b, c);
      assertTrue(crosser.robustCrossing(d) < 0);
      assertTrue(crosser.robustCrossing(c) < 0);
    }
  }

  private void checkFaceClipping(S2Point aRaw, S2Point bRaw) {
    S2Point a = aRaw.normalize();
    S2Point b = bRaw.normalize();

    // First we test GetFaceSegments.
    FaceSegment[] segments = FaceSegment.allFaces();
    int n = S2EdgeUtil.getFaceSegments(a, b, segments);
    assertTrue(n >= 1);

    R2Rect biunit = new R2Rect(new R1Interval(-1, 1), new R1Interval(-1, 1));

    // The first and last vertices should approximately equal A and B.
    assertTrue(
        a.angle(S2Projections.faceUvToXyz(segments[0].face, segments[0].a))
            <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
    assertTrue(
        b.angle(S2Projections.faceUvToXyz(segments[n - 1].face, segments[n - 1].b))
            <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);

    S2Point norm = robustCrossProd(a, b).normalize();
    S2Point aTangent = norm.crossProd(a);
    S2Point bTangent = b.crossProd(norm);
    for (int i = 0; i < n; ++i) {
      // Vertices may not protrude outside the biunit square.
      assertTrue(biunit.contains(segments[i].a));
      assertTrue(biunit.contains(segments[i].b));
      if (i == 0) {
        continue;
      }

      // The two representations of each interior vertex (on adjacent faces) must correspond to
      // exactly the same S2Point.
      assertNotSame(segments[i - 1].face, segments[i].face);
      assertEquals(
          S2Projections.faceUvToXyz(segments[i - 1].face, segments[i - 1].b),
          S2Projections.faceUvToXyz(segments[i].face, segments[i].a));

      // Interior vertices should be in the plane containing A and B, and should be contained in the
      // wedge of angles between A and B (i.e., the dot products with aTangent and bTangent should
      // be non-negative).
      S2Point p = S2Projections.faceUvToXyz(segments[i].face, segments[i].a).normalize();
      assertTrue(abs(p.dotProd(norm)) <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
      assertTrue(p.dotProd(aTangent) >= -S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
      assertTrue(p.dotProd(bTangent) >= -S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
    }

    // Now we test ClipToPaddedFace (sometimes with a padding of zero). We do this by defining an
    // (x,y) coordinate system for the plane containing AB, and converting points along the great
    // circle AB to angles in the range [-Pi, Pi]. We then accumulate the angle intervals spanned by
    // each clipped edge; the union over all 6 faces should approximately equal the interval covered
    // by the original edge.
    double padding = data.oneIn(10) ? 0.0 : 1e-10 * pow(1e-5, data.nextDouble());
    S2Point xAxis = a;
    S2Point yAxis = aTangent;
    S1Interval expectedAngles = new S1Interval(0, a.angle(b));
    S1Interval maxAngles = expectedAngles.expanded(S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
    S1Interval actualAngles = null;
    for (int face = 0; face < 6; ++face) {
      R2Vector aUv = new R2Vector();
      R2Vector bUv = new R2Vector();
      if (S2EdgeUtil.clipToPaddedFace(a, b, face, padding, aUv, bUv)) {
        S2Point aClip = S2Projections.faceUvToXyz(face, aUv).normalize();
        S2Point bClip = S2Projections.faceUvToXyz(face, bUv).normalize();
        assertTrue(abs(aClip.dotProd(norm)) <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
        assertTrue(abs(bClip.dotProd(norm)) <= S2EdgeUtil.FACE_CLIP_ERROR_RADIANS);
        if (aClip.angle(a) > S2EdgeUtil.FACE_CLIP_ERROR_RADIANS) {
          assertDoubleNear(1 + padding, max(abs(aUv.x), abs(aUv.y)));
        }
        if (bClip.angle(b) > S2EdgeUtil.FACE_CLIP_ERROR_RADIANS) {
          assertDoubleNear(1 + padding, max(abs(bUv.x), abs(bUv.y)));
        }
        double aAngle = atan2(aClip.dotProd(yAxis), aClip.dotProd(xAxis));
        double bAngle = atan2(bClip.dotProd(yAxis), bClip.dotProd(xAxis));
        // Rounding errors may cause bAngle to be slightly less than aAngle. We handle this by
        // constructing the interval with fromPointPair(), which is okay since the interval length
        // is much less than PI. .
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
    S2Point a = p.mul(data.nextInt(3) - 1).add(q.mul(data.nextInt(3) - 1));
    if (data.oneIn(10)) {
      // This perturbation often has no effect except on coordinates that are zero, in which case
      // the perturbed value is so small that operations on it often result in underflow.
      a = a.add(data.getRandomPoint().mul(pow(1e-300, data.nextDouble())));
    } else if (data.oneIn(2)) {
      // For coordinates near 1 (say > 0.5), this perturbation yields values that are only a few
      // representable values away from the initial value.
      a = a.add(data.getRandomPoint().mul(4 * S2.DBL_EPSILON));
    } else {
      // A perturbation whose magnitude is in the range [1e-25, 1e-10].
      a = a.add(data.getRandomPoint().mul(1e-10 * pow(1e-15, data.nextDouble())));
    }
    if (a.norm2() < Double.MIN_VALUE) {
      // If a.norm2() is denormalized, normalize() loses too much precision.
      return perturbedCornerOrMidpoint(p, q);
    }
    return a;
  }

  @Test
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

    // Comprehensively test edges that are difficult to handle, especially those that nearly follow
    // one of the 12 cube edges.
    R2Rect biunit = new R2Rect(new R1Interval(-1, 1), new R1Interval(-1, 1));
    final int kIters = 1000; // Test passes with 1e6 iterations
    for (int iter = 0; iter < kIters; ++iter) {
      // Choose two adjacent cube corners P and Q.
      int face = data.nextInt(6);
      int i = data.nextInt(4);
      int j = (i + 1) & 3;
      S2Point p = S2Projections.faceUvToXyz(face, biunit.getVertex(i));
      S2Point q = S2Projections.faceUvToXyz(face, biunit.getVertex(j));

      // Now choose two points that are nearly in the plane of PQ, preferring points that are near
      // cube corners, face midpoints, or edge midpoints.
      S2Point a = perturbedCornerOrMidpoint(p, q);
      S2Point b = perturbedCornerOrMidpoint(p, q);
      checkFaceClipping(a, b);
    }
  }

  /**
   * Return a uniformly distributed "double" in the range [min, max), where min is the lower of 'a'
   * and 'b', and max is the higher.
   */
  private double unorderedUniform(double a, double b) {
    return (a < b) ? data.uniform(a, b) : data.uniform(b, a);
  }

  /**
   * Choose a random point in the rectangle defined by points A and B, sometimes returning a point
   * on the edge AB or the points A and B themselves.
   */
  private R2Vector chooseRectPoint(R2Vector a, R2Vector b) {
    if (data.oneIn(5)) {
      return data.oneIn(2) ? a : b;
    } else if (data.oneIn(3)) {
      return R2Vector.add(a, R2Vector.mul(R2Vector.sub(b, a), data.nextDouble()));
    } else {
      return new R2Vector(unorderedUniform(a.x, b.x), unorderedUniform(a.y, b.y));
    }
  }

  /**
   * Given a point X on the line AB (which is checked), return the fraction "t" such that x =
   * (1-t)*a + t*b. Return 0 if A = B.
   */
  private double getFraction(R2Vector x, R2Vector a, R2Vector b) {
    // A bound for the error in edge clipping plus the error in the calculation below (which is
    // similar to IntersectsRect).
    final double kError =
        S2EdgeUtil.EDGE_CLIP_ERROR_UV_DIST + S2EdgeUtil.INTERSECTS_RECT_ERROR_UV_DIST;
    if (a.equals(b)) {
      return 0.0;
    }
    R2Vector dir = R2Vector.normalize(R2Vector.sub(b, a));
    assertTrue(abs(R2Vector.sub(x, a).dotProd(dir.ortho())) <= kError);
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
   * Given an edge AB and a rectangle "clip", verify that intersectsRect(), clipEdge(), and
   * clipEdgeBound() produce consistent results.
   */
  private void checkClipEdge(R2Vector a, R2Vector b, R2Rect clip) {
    // A bound for the error in edge clipping plus the error in the intersectsRect calculation
    // below.
    final double kError =
        S2EdgeUtil.EDGE_CLIP_ERROR_UV_DIST + S2EdgeUtil.INTERSECTS_RECT_ERROR_UV_DIST;
    R2Vector aClipped = new R2Vector();
    R2Vector bClipped = new R2Vector();
    if (!S2EdgeUtil.clipEdge(a, b, clip, aClipped, bClipped)) {
      assertFalse(S2EdgeUtil.intersectsRect(a, b, clip.expanded(-kError)));
    } else {
      assertTrue(S2EdgeUtil.intersectsRect(a, b, clip.expanded(kError)));
      // Check that the clipped points lie on the edge AB, and that the points have the expected
      // order along the segment AB.
      assertTrue(getFraction(aClipped, a, b) <= getFraction(bClipped, a, b));
      // Check that the clipped portion of AB is as large as possible.
      checkPointOnBoundary(aClipped, a, clip);
      checkPointOnBoundary(bClipped, b, clip);
    }

    // Choose a random initial bound to pass to clipEdgeBound.
    R2Rect initialClip = R2Rect.fromPointPair(chooseRectPoint(a, b), chooseRectPoint(a, b));
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
    if (data.oneIn(5)) {
      return data.oneIn(2) ? clip.lo() : clip.hi();
    } else {
      switch (data.nextInt(3)) {
        case 0:
          return clip.lo() - data.nextDouble();
        case 1:
          return clip.hi() + data.nextDouble();
        default:
          return clip.lo() + data.nextDouble() * clip.getLength();
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
    if (data.oneIn(10)) {
      // Return a point on one of the two extended diagonals.
      int diag = data.nextInt(2);
      double t = data.uniform(-1, 2);
      return R2Vector.add(
          R2Vector.mul(clip.getVertex(diag), 1 - t), R2Vector.mul(clip.getVertex(diag + 2), t));
    } else {
      return new R2Vector(chooseEndpoint(clip.x()), chooseEndpoint(clip.y()));
    }
  }

  /**
   * Given a rectangle "clip", test the S2EdgeUtil edge clipping methods using many edges that are
   * randomly constructed to trigger special cases.
   */
  private void checkEdgeClipping(R2Rect clip) {
    final int kIters = 1000; // Test passes with 1e6 iterations
    for (int iter = 0; iter < kIters; ++iter) {
      checkClipEdge(chooseEndpoint(clip), chooseEndpoint(clip), clip);
    }
  }

  @Test
  public void testEdgeClipping() {
    // Test clipping against random rectangles.
    for (int i = 0; i < 5; ++i) {
      checkEdgeClipping(
          R2Rect.fromPointPair(
              new R2Vector(data.uniform(-1, 1), data.uniform(-1, 1)),
              new R2Vector(data.uniform(-1, 1), data.uniform(-1, 1))));
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

  /** A visitor of crossing edge pairs AB and CD. */
  private interface CrossingEdgesVisitor {
    /** Edges (a, b) and (c, d) intersect near p with given slope. */
    void visit(S2Point a, S2Point b, S2Point c, S2Point d, S2Point p, double slope);
  }

  /** Visits many (5000) randomly generated crossing edge pairs. */
  private void checkCrossingEdges(CrossingEdgesVisitor visitor) {
    for (int i = 0; i < 5000; ++i) {
      // Construct two random edges that intersect. The angle between these edges (expressed as a
      // slope) is chosen randomly between 1e-15 and 1e15 such that its logarithm is uniformly
      // distributed. Similarly, two edge lengths approximately between 1e-15 and 1 are chosen. The
      // edge endpoints are chosen such that they are often very close to the other edge (i.e.,
      // barely crossing). Taken together this ensures that we test both long and very short edges
      // that intersect at both large and very small angles.
      Matrix frame = S2.getFrame(data.getRandomPoint());
      S2Point p = frame.getCol(0);
      S2Point d1 = frame.getCol(1);
      S2Point d2 = frame.getCol(2);

      double slope = 1e-15 * pow(1e30, data.nextDouble());
      d2 = d1.add(d2.mul(slope)).normalize();

      S2Point a;
      S2Point b;
      S2Point c;
      S2Point d;
      do {
        double abLen = pow(1e-15, data.nextDouble());
        double cdLen = pow(1e-15, data.nextDouble());
        double aFraction = pow(1e-5, data.nextDouble());
        aFraction = data.nextBoolean() ? aFraction : 1 - aFraction;
        double cFraction = pow(1e-5, data.nextDouble());
        cFraction = data.nextBoolean() ? cFraction : 1 - cFraction;
        a = p.sub(d1.mul(aFraction * abLen)).normalize();
        b = p.add(d1.mul((1 - aFraction) * abLen)).normalize();
        c = p.sub(d2.mul(cFraction * cdLen)).normalize();
        d = p.add(d2.mul((1 - cFraction) * cdLen)).normalize();
      } while (S2EdgeUtil.robustCrossing(a, b, c, d) <= 0);
      visitor.visit(a, b, c, d, p, slope);
    }
  }
}
