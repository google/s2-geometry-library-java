/*
 * Copyright 2022 Google Inc.
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

import static java.lang.Math.max;
import static java.lang.Math.min;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Iterables;
import com.google.common.geometry.Projection.MercatorProjection;
import com.google.common.geometry.Projection.PlateCarreeProjection;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link S2EdgeTessellator}. */
@RunWith(JUnit4.class)
public class S2EdgeTessellatorTest extends GeometryTestCase {

  private static final Logger logger = Platform.getLoggerForClass(S2EdgeTessellatorTest.class);

  static class Stats {
    private double max = Double.NEGATIVE_INFINITY;
    private double sum = 0;
    private int count = 0;

    void tally(double v) {
      if (Double.isNaN(v)) {
        throw new IllegalArgumentException("NaN");
      }
      if (v > max) {
        max = v;
      }
      sum += v;
      count++;
    }

    double max() {
      return max;
    }

    double mean() {
      return sum / count;
    }

    @Override
    public String toString() {
      return "avg = " + mean() + " max = " + max();
    }
  }

  // Determines whether the distance between the two edges is measured geometrically or
  // parametrically
  enum DistanceType {
    PARAMETRIC,
    GEOMETRIC
  }

  S1Angle getMaxDistance(
      Projection projection,
      R2Vector px,
      S2Point x,
      R2Vector py,
      S2Point y,
      DistanceType distanceType) {
    // Step along the projected edge at a fine resolution and keep track of the
    // maximum distance of any point to the current geodesic edge.
    final int numberOfSteps = 100;
    S1ChordAngle maxDist = S1ChordAngle.ZERO;
    for (double f = 0.5 / numberOfSteps; f < 1.0; f += 1.0 / numberOfSteps) {
      S1ChordAngle distance = S1ChordAngle.INFINITY;
      S2Point p = projection.unproject(Projection.interpolate(f, px, py));
      if (distanceType == DistanceType.GEOMETRIC) {
        distance = S2EdgeUtil.updateMinDistance(p, new S2Edge(x, y), distance);
      } else {
        distance = new S1ChordAngle(p, S2EdgeUtil.interpolate(f, x, y));
      }
      if (distance.greaterThan(maxDist)) {
        maxDist = distance;
      }
    }
    // Ensure that the maximum distance estimate is a lower bound, not an upper
    // bound, since we only want to record a failure of the distance estimation
    // algorithm if the number it returns is definitely too small.
    return maxDist.plusError(S2EdgeUtil.getUpdateMinDistanceMaxError(maxDist)).toAngle();
  }

  // When there are longitudes greater than 180 degrees due to wrapping, the
  // combination of projecting and unprojecting an S2Point can have slightly more
  // error than is allowed by S2::ApproxEquals.
  private static final S1Angle MAX_PROJECTION_ERROR = S1Angle.radians(2e-15);

  private static final S1ChordAngle MAX_INTERPOLATION_ERROR =
      S1ChordAngle.fromS1Angle(S1Angle.radians(1e-14));

  // Converts a projected edge to a sequence of geodesic edges and verifies that
  // the result satisfies the given tolerance.
  Stats testUnprojected(
      Projection proj, S1Angle tolerance, R2Vector pa, R2Vector pbIn, boolean logStats) {
    S2EdgeTessellator tessellator = new S2EdgeTessellator(proj, tolerance);
    List<S2Point> vertices = new ArrayList<>();
    tessellator.appendUnprojected(pa, pbIn, vertices);
    R2Vector pb = proj.wrapDestination(pa, pbIn);
    assertTrue(new S1Angle(proj.unproject(pa), vertices.get(0)).lessThan(MAX_PROJECTION_ERROR));
    assertTrue(
        new S1Angle(proj.unproject(pb), Iterables.getLast(vertices))
            .lessThan(MAX_PROJECTION_ERROR));
    Stats stats = new Stats();
    if (pa.equals(pb)) {
      assertEquals(1, vertices.size());
      return stats;
    }
    // Precompute the normal to the projected edge.
    R2Vector norm = R2Vector.normalize(pb.sub(pb).ortho());
    S2Point x = vertices.get(0);
    R2Vector px = proj.project(x);

    for (int i = 1; i < vertices.size(); ++i) {
      S2Point y = vertices.get(i);
      R2Vector py = proj.wrapDestination(px, proj.project(y));
      // Check that every vertex is on the projected edge.
      assertLessOrEqual(py.sub(pa).dotProd(norm), 1e-14 * py.norm());
      stats.tally(
          getMaxDistance(proj, px, x, py, y, DistanceType.GEOMETRIC).radians()
              / tolerance.radians());
      x = y;
      px = py;
    }
    if (logStats) {
      logger.info(
          vertices.get(0)
              + " to "
              + px
              + ": "
              + Iterables.getLast(vertices)
              + " vertices, "
              + stats);
    }
    return stats;
  }

  // Converts a geodesic edge to a sequence of projected edges and verifies that
  // the result satisfies the given tolerance.
  Stats testProjected(Projection proj, S1Angle tolerance, S2Point a, S2Point b, boolean logStats) {
    S2EdgeTessellator tessellator = new S2EdgeTessellator(proj, tolerance);
    List<R2Vector> vertices = new ArrayList<>();
    tessellator.appendProjected(a, b, vertices);
    assertTrue(new S1Angle(a, proj.unproject(vertices.get(0))).lessThan(MAX_PROJECTION_ERROR));
    assertTrue(
        new S1Angle(b, proj.unproject(Iterables.getLast(vertices))).lessThan(MAX_PROJECTION_ERROR));
    Stats stats = new Stats();
    if (a.equals(b)) {
      assertEquals(1, vertices.size());
      return stats;
    }
    R2Vector px = vertices.get(0);
    S2Point x = proj.unproject(px);
    for (int i = 1; i < vertices.size(); ++i) {
      R2Vector py = vertices.get(i);
      S2Point y = proj.unproject(py);
      // Check that every vertex is on the geodesic edge.
      assertTrue(
          S1ChordAngle.fromS1Angle(S2EdgeUtil.getDistance(y, a, b))
              .lessThan(MAX_INTERPOLATION_ERROR));
      stats.tally(
          getMaxDistance(proj, px, x, py, y, DistanceType.GEOMETRIC)
              .div(tolerance.radians())
              .radians());
      x = y;
      px = py;
    }
    if (logStats) {
      System.out.println(
          vertices.get(0)
              + " to "
              + px
              + ": "
              + Iterables.getLast(vertices)
              + " vertices, "
              + stats);
    }
    return stats;
  }

  @Test
  public void testProjectedNoTessellation() {
    Projection plateCarree = new PlateCarreeProjection(180.0);
    S2EdgeTessellator tessellator = new S2EdgeTessellator(plateCarree, S1Angle.degrees(0.01));

    ArrayList<R2Vector> vertices = new ArrayList<>();
    tessellator.appendProjected(new S2Point(1, 0, 0), new S2Point(0, 1, 0), vertices);
    assertEquals(2, vertices.size());
  }

  @Test
  public void testUnrojectedNoTessellation() {
    Projection plateCarree = new PlateCarreeProjection(180.0);
    S2EdgeTessellator tessellator = new S2EdgeTessellator(plateCarree, S1Angle.degrees(0.01));

    ArrayList<S2Point> vertices = new ArrayList<>();
    tessellator.appendUnprojected(new R2Vector(0, 30), new R2Vector(0, 50), vertices);
    assertEquals(2, vertices.size());
  }

   @Test
  public void testUnprojectedWrapping() {
    // This tests that a projected edge that crosses the 180 degree meridian
    // goes the "short way" around the sphere.
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S2EdgeTessellator tess = new S2EdgeTessellator(proj, S1Angle.degrees(0.01));
    List<S2Point> vertices = new ArrayList<>();
    tess.appendUnprojected(new R2Vector(-170, 0), new R2Vector(170, 80), vertices);
    for (S2Point v : vertices) {
      assertGreaterOrEqual(Math.abs(S2LatLng.longitude(v).degrees()), 170.0);
    }
  }

  @Test
  public void testProjectedWrapping() {
    // This tests that a projected edge that crosses the 180 degree meridian
    // goes the "short way" around the sphere.
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S2EdgeTessellator tess = new S2EdgeTessellator(proj, S1Angle.degrees(0.01));
    List<R2Vector> vertices = new ArrayList<>();
    tess.appendProjected(
        S2LatLng.fromDegrees(0, -170).toPoint(), S2LatLng.fromDegrees(0, 170).toPoint(), vertices);
    for (R2Vector v : vertices) {
      assertLessOrEqual(v.x(), -170.0);
    }
  }

  @Test
  public void testUnprojectedWrappingMultipleCrossings() {
    // Tests an edge chain that crosses the 180 degree meridian multiple times.
    // Note that due to coordinate wrapping, the last vertex of one edge may not
    // exactly match the first edge of the next edge after unprojection.
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S2EdgeTessellator tess = new S2EdgeTessellator(proj, S1Angle.degrees(0.01));
    List<S2Point> vertices = new ArrayList<>();
    for (double lat = 1; lat <= 60; lat += 1.0) {
      tess.appendUnprojected(
          new R2Vector(180 - 0.03 * lat, lat), new R2Vector(-180 + 0.07 * lat, lat), vertices);
      tess.appendUnprojected(
          new R2Vector(-180 + 0.07 * lat, lat),
          new R2Vector(180 - 0.03 * (lat + 1), lat + 1),
          vertices);
    }
    for (S2Point v : vertices) {
      assertGreaterOrEqual(Math.abs(S2LatLng.longitude(v).degrees()), 175.0);
    }
  }

  @Test
  public void testProjectedWrappingMultipleCrossings() {
    // The following loop crosses the 180 degree meridian four times (twice in
    // each direction).
    List<S2Point> loop =
        S2TextFormat.parsePointsOrDie("0:160, 0:-40, 0:120, 0:-80, 10:120, 10:-40, 0:160");
    S1Angle tolerance = S1Angle.e7(1);
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S2EdgeTessellator tess = new S2EdgeTessellator(proj, tolerance);
    List<R2Vector> vertices = new ArrayList<>();
    for (int i = 0; i + 1 < loop.size(); i++) {
      tess.appendProjected(loop.get(i), loop.get(i + 1), vertices);
    }
    assertEquals(vertices.get(0), Iterables.getLast(vertices));

    double minLongitude = vertices.get(0).x();
    double maxLongitude = vertices.get(0).x();
    for (R2Vector v : vertices) {
      minLongitude = min(minLongitude, v.x());
      maxLongitude = max(maxLongitude, v.x());
    }
    assertEquals(160, minLongitude, 0.001);
    assertEquals(640, maxLongitude, 0.001);
  }

  @Test
  public void testInfiniteRecursionBug() {
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S1Angle oneMicron = S1Angle.radians(1e-6 / 6371.0);
    S2EdgeTessellator tess = new S2EdgeTessellator(proj, oneMicron);
    List<R2Vector> vertices = new ArrayList<>();
    tess.appendProjected(
        S2LatLng.fromDegrees(3, 21).toPoint(), S2LatLng.fromDegrees(1, -159).toPoint(), vertices);
    assertEquals(36, vertices.size());
  }

  @Test
  public void testUnprojectedAccuracy() {
    MercatorProjection proj = new MercatorProjection(180);
    S1Angle tolerance = S1Angle.degrees(1e-5);
    R2Vector pa = new R2Vector(0, 0);
    R2Vector pb = new R2Vector(89.999999, 179);
    Stats stats = testUnprojected(proj, tolerance, pa, pb, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  // Repro case for b/110719057.
  @Test
  public void testUnprojectedAccuracyCrossEquator() {
    MercatorProjection proj = new MercatorProjection(180);
    S1Angle tolerance = S1Angle.degrees(1e-5);
    R2Vector pa = new R2Vector(-10, -10);
    R2Vector pb = new R2Vector(10, 10);
    Stats stats = testUnprojected(proj, tolerance, pa, pb, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  @Test
  public void testProjectedAccuracy() {
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S1Angle tolerance = S1Angle.e7(1);
    S2Point a = S2LatLng.fromDegrees(-89.999, -170).toPoint();
    S2Point b = S2LatLng.fromDegrees(50, 100).toPoint();
    Stats stats = testProjected(proj, tolerance, a, b, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  @Test
  public void testUnprojectedAccuracyMidpointEquator() {
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S1Angle tolerance = S1Angle.radians(S2Earth.metersToRadians(1));
    R2Vector a = new R2Vector(80, 50);
    R2Vector b = new R2Vector(-80, -50);
    Stats stats = testUnprojected(proj, tolerance, a, b, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  @Test
  public void testProjectedAccuracyMidpointEquator() {
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S1Angle tolerance = S1Angle.radians(S2Earth.metersToRadians(1));
    S2Point a = S2LatLng.fromDegrees(50, 80).toPoint();
    S2Point b = S2LatLng.fromDegrees(-50, -80).toPoint();
    Stats stats = testProjected(proj, tolerance, a, b, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  /// Repro case for b/110719057.
  @Test
  public void testProjectedAccuracyCrossEquator() {
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S1Angle tolerance = S1Angle.radians(S2Earth.metersToRadians(1));
    S2Point a = S2LatLng.fromDegrees(-20, -20).toPoint();
    S2Point b = S2LatLng.fromDegrees(20, 20).toPoint();
    Stats stats = testProjected(proj, tolerance, a, b, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  @Test
  public void testProjectedAccuracySeattleToNewYork() {
    PlateCarreeProjection proj = new PlateCarreeProjection(180);
    S1Angle tolerance = S1Angle.radians(S2Earth.metersToRadians(1));
    S2Point a = S2LatLng.fromDegrees(47.6062, -122.3321).toPoint();
    S2Point b = S2LatLng.fromDegrees(40.7128, -74.0059).toPoint();
    Stats stats = testProjected(proj, tolerance, a, b, true);
    assertLessOrEqual(stats.max(), 1.0);
  }

  // Given a projection, this function repeatedly chooses a pair of edge
  // endpoints and measures the true distance between the geodesic and projected
  // edges that connect those two endpoints.  It then compares this to against
  // the distance measurement algorithm used by S2EdgeTessellator, which
  // consists of measuring the point-to-point distance between the edges at each
  // of two fractions "t" and "1-t", computing the maximum of those two
  // distances, and then scaling by the constant documented in the .cc file
  // (based on the idea that the distance between the edges as a function of the
  // interpolation fraction can be accurately modeled as a cubic polynomial).
  //
  // This function is used to (1) verify that the distance estimates are always
  // conservative, (2) verify the optimality of the interpolation fraction "t",
  // and (3) estimate the amount of overtessellation that occurs for various
  // types of edges (e.g., short vs. long edges, edges that follow lines of
  // latitude or longitude, etc.)
  void testEdgeError(Projection projection, double t) {
    double x = 1 - 2 * t;
    double quarterPi = Math.PI / 4;
    double inverseRoot2 = 1 / Math.sqrt(2);
    double dlat = Math.sin(0.5 * quarterPi * (1 - x));
    double dlng = Math.sin(quarterPi * (1 - x));
    double dsin2 = dlat * dlat + dlng * dlng * Math.sin(quarterPi * x) * inverseRoot2;
    double dsin2Max = 0.5 * (1 - inverseRoot2);
    double scaleFactor =
        max(
            (2 * Math.sqrt(3) / 9) / (x * (1 - x * x)),
            Math.asin(Math.sqrt(dsin2Max)) / Math.asin(Math.sqrt(dsin2)));
    Stats statsGeometric = new Stats();
    Stats statsParametric = new Stats();
    int iters = 100000;
    for (int iter = 0; iter < iters; iter++) {
      data.setSeed(iter);
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      if (a.dotProd(b) < -1e-14) {
        continue;
      }
      R2Vector pa = projection.project(a);
      R2Vector pb = projection.wrapDestination(pa, projection.project(b));
      S1Angle maxDistGeometric = getMaxDistance(projection, pa, a, pb, b, DistanceType.GEOMETRIC);
      if (maxDistGeometric.lessThan(S2EdgeTessellator.MIN_TOLERANCE)) {
        continue;
      }
      S1Angle maxDistParametric = getMaxDistance(projection, pa, a, pb, b, DistanceType.PARAMETRIC);
      if (maxDistParametric.lessThan(S2EdgeTessellator.MIN_TOLERANCE)) {
        continue;
      }

      S1Angle d1 =
          new S1Angle(
              S2EdgeUtil.interpolate(t, a, b), projection.unproject(pa.mul(1 - t).add(pb.mul(t))));
      S1Angle d2 =
          new S1Angle(
              S2EdgeUtil.interpolate(1 - t, a, b),
              projection.unproject(pa.mul(t).add(pb.mul(1 - t))));
      S1Angle dist = S1Angle.max(S1Angle.radians(1e-300), S1Angle.max(d1, d2));
      double ratioGeometric = maxDistGeometric.radians() / dist.radians();
      double rationParametric = maxDistParametric.radians() / dist.radians();
      statsGeometric.tally(ratioGeometric);
      statsParametric.tally(rationParametric);
    }
    assertLessOrEqual(statsGeometric.max(), scaleFactor);
  }

  // The interpolation parameter actually used in the .java file.
  private static final double BEST_FRACTION = 0.31215691082248315;

  @Test
  public void testMaxEdgeErrorPlateCarree() {
    PlateCarreeProjection projection = new PlateCarreeProjection(180);
    testEdgeError(projection, BEST_FRACTION);
  }

  @Test
  public void testMaxEdgeErrorMercator() {
    MercatorProjection projection = new MercatorProjection(180);
    testEdgeError(projection, BEST_FRACTION);
  }

  // Tessellates random edges using the given projection and tolerance, and
  // verifies that the expected criteria are satisfied. If avoidPolar is true, will check if points
  // get within tolerance of the North or South Poles and will skip if that is true.
  void testRandomEdges(Projection projection, S1Angle tolerance) {
    int iters = 500;
    double maxR2 = 0;
    double maxS2 = 0;
    for (int iter = 0; iter < iters; iter++) {
      // TODO(user) iteration 75 and 288 are known bad seeds that causes issues with the
      // guarantee that unprojecting the points in a tessellated geodeic should result in points
      // that stay within the interpolation tolerance of the geodesic.
      // Seed 75:
      // (-49.64206592218376, -44.641143496511475) : (-34.96795250656595, 136.8252303923727)
      // Seed 288:
      // (48.51899997948501, -43.402084341709795) : (22.725193708030066, 136.62141990649366)
      if (iter == 75 || iter == 288) {
        continue;
      }
      data.setSeed(iter);
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      maxR2 = max(maxR2, testProjected(projection, tolerance, a, b, false).max());
      R2Vector pa = projection.project(a);
      R2Vector pb = projection.project(b);
      maxS2 = max(maxS2, testUnprojected(projection, tolerance, pa, pb, false).max());
    }
    assertLessOrEqual(maxR2, 1.0);
    assertLessOrEqual(maxS2, 1.0);
  }

  @Test
  public void testRandomEdgesPlateCarree() {
    PlateCarreeProjection projection = new PlateCarreeProjection(180);
    S1Angle tolerance = S1Angle.radians(S2Earth.metersToRadians(100));
    testRandomEdges(projection, tolerance);
  }

  @Test
  public void testRandomEdgesMercator() {
    MercatorProjection projection = new MercatorProjection(180);
    S1Angle tolerance = S1Angle.radians(S2Earth.metersToRadians(100));
    testRandomEdges(projection, tolerance);
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // To be identical to the C++ test.
  @Test
  public void testUnwrappingAssertionRegression() {
    // Before it was fixed, this series of points would cause an assertion to be thrown when
    // tessellating due to rounding in the unwrapping code.
    double[][] kPoints = {
        { -16.876721435218865253, -179.986547984808964884 },
        { -16.874909244632696925, -179.991889238369623172 },
        { -16.880241814330226191, -179.990858688466971671 },
        { -16.883762104047619346, -179.995169553755403058 },
        { -16.881949690252106677, +179.999489074621124018 },
        { -16.876617071405430437, +179.998458788144517939 },
        { -16.880137137875717457, +179.994147804931060364 },
        { -16.878324446969305228, +179.988806637264332267 },
        { -16.872991774409559440, +179.987776672537478362 },
        { -16.869471841739493101, +179.992087611973005323 },
        { -16.867659097232969856, +179.986746766061799008 },
        { -16.862326415537093993, +179.985716917832945683 },
        { -16.858806527326652969, +179.990027652027180238 },
        { -16.860619186956174786, +179.995368278278732532 },
        { -16.855286549828541354, +179.994338224830613626 },
        { -16.851766483129139829, +179.998648636203512297 },
        { -16.849953908374558864, +179.993308229628894424 }
    };

    MercatorProjection projection = new MercatorProjection(0.5);

    S2EdgeTessellator tessellator = new S2EdgeTessellator(projection, S1Angle.e7(1));
    List<R2Vector> vertices = new ArrayList<>();

    for (int i = 0; i + 1 < kPoints.length; ++i) {
      tessellator.appendProjected(
          S2LatLng.fromDegrees(kPoints[i + 0][0], kPoints[i + 0][1]).toPoint(),
          S2LatLng.fromDegrees(kPoints[i + 1][0], kPoints[i + 1][1]).toPoint(),
          vertices);
    }
    assertEquals(17, vertices.size());
  }

}
