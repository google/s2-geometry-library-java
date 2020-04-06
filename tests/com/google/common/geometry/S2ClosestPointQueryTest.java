/*
 * Copyright 2015 Google Inc.
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
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.geometry.S2ClosestPointQuery.Result;
import java.util.List;

/** Verifies {@link S2ClosestPointQuery}. */
@GwtCompatible
public class S2ClosestPointQueryTest extends GeometryTestCase {
  /** The approximate radius of S2Cap from which query points are chosen. */
  protected static final S1Angle QUERY_RADIUS = kmToAngle(10);

  /**
   * An approximate bound on the distance measurement error for "reasonable" distances (say, less
   * than Pi/2) due to using S1ChordAngle.
   */
  private static final double MAX_CHORD_ANGLE_ERROR = 1e-15;

  public void testNoPoints() {
    S2PointIndex<Integer> index = new S2PointIndex<>();
    S2ClosestPointQuery<Integer> query = new S2ClosestPointQuery<>(index);
    assertNull(query.findClosestPoint(S2Point.X_POS));
    assertEquals(0, query.findClosestPoints(S2Point.X_POS).size());
  }

  public void testManyDuplicatePoints() {
    int numPoints = 10000;
    S2PointIndex<Integer> index = new S2PointIndex<>();
    for (int i = 0; i < numPoints; ++i) {
      index.add(S2Point.X_POS, i);
    }
    S2ClosestPointQuery<Integer> query = new S2ClosestPointQuery<>(index);
    assertEquals(numPoints, new PointTarget(S2Point.X_POS).findClosestPoints(query).size());
  }

  @GwtIncompatible("Edge tests require Real.java, which is incompatible")
  public void testPoints() {
    int numIndexes = 10;
    int numPoints = 1000;
    int numQueries = 50;
    for (PointFactory factory : PointFactory.values()) {
      checkFactory(factory, numIndexes, numPoints, numQueries);
    }
  }

  /**
   * Check that result set "x" contains all the expected results from "y", and does not include any
   * duplicate results.
   */
  private static void checkResultSet(
      List<Result<Integer>> x,
      List<Result<Integer>> y,
      int maxSize,
      S1Angle maxDistance,
      S1Angle maxError,
      S1Angle maxPruningError,
      String label) {
    // Results should be sorted by distance.
    assertTrue(Ordering.natural().reverse().isOrdered(x));

    // Make sure there are no duplicate values.
    assertEquals("Result set contains duplicates", Sets.newHashSet(x).size(), x.size());

    // Result set X should contain all the items from U whose distance is less
    // than "limit" computed below.
    double limit = 0;
    if (x.size() < maxSize) {
      // Result set X was not limited by "max_size", so it should contain all
      // the items up to "max_distance", except that a few items right near the
      // distance limit may be missed because the distance measurements used for
      // pruning S2Cells are not conservative.
      limit = maxDistance.radians() - maxPruningError.radians();
    } else if (!x.isEmpty()) {
      // Result set X contains only the closest "max_size" items, to within a
      // tolerance of "max_error + max_pruning_error".
      limit =
          x.get(x.size() - 1).distance().toAngle().radians()
              - maxError.radians()
              - maxPruningError.radians();
    }
    for (int i = 0; i < y.size(); ++i) {
      Result<Integer> item = y.get(i);
      if (item.distance().toAngle().radians() < limit) {
        assertEquals(label + " " + item, 1, Iterables.frequency(x, item));
      }
    }
  }

  /**
   * Compares two sets of "closest" items, where "expected" is computed via brute force (i.e.,
   * considering every possible candidate) and "actual" is computed using a spatial data structure.
   * Here "maxSize" is a bound on the maximum number of items, "maxDistance" is a limit on the
   * distance to any item, and "maxError" is the maximum error allowed when selecting which items
   * are closest.
   */
  private static void compareResults(
      List<Result<Integer>> expected,
      List<Result<Integer>> actual,
      int maxSize,
      S1Angle maxDistance,
      S1Angle maxError) {
    S1Angle maxPruningError = S1Angle.radians(1e-15);
    checkResultSet(actual, expected, maxSize, maxDistance, maxError, maxPruningError, "Missing");
    checkResultSet(expected, actual, maxSize, maxDistance, maxError, S1Angle.ZERO, "Extra");
  }

  /** An abstract class that adds points to an S2PointIndex for benchmarking. */
  protected enum PointFactory {
    /**
     * Generator for points regularly spaced along a circle. The circle is centered within the query
     * cap and occupies 25% of its area, so that random query points have a 25% chance of being
     * inside the circle.
     *
     * <p>Points along a circle are nearly the worst case for distance calculations, since many
     * points are nearly equidistant from any query point that is not immediately adjacent to the
     * circle.
     */
    CIRCLE {
      @Override
      public List<S2Point> createPoints(
          S2ClosestPointQueryTest helper, S2Cap queryCap, int numPoints) {
        return S2Loop.makeRegularVertices(
            queryCap.axis(), S1Angle.radians(0.5 * queryCap.angle().radians()), numPoints);
      }
    },
    /** Generator for points of a fractal whose convex hull approximately matches the query cap. */
    FRACTAL {
      @Override
      public List<S2Point> createPoints(
          S2ClosestPointQueryTest helper, S2Cap queryCap, int numPoints) {
        S2FractalBuilder builder = new S2FractalBuilder(helper.rand);
        builder.setLevelForApproxMaxEdges(numPoints);
        builder.setFractalDimension(1.5);
        return builder.makeVertices(helper.getRandomFrameAt(queryCap.axis()), queryCap.angle());
      }
    },
    /** Generator for points on a square grid that includes the entire query cap. */
    GRID {
      @Override
      public List<S2Point> createPoints(
          S2ClosestPointQueryTest helper, S2Cap queryCap, int numPoints) {
        int sqrtNumPoints = (int) Math.ceil(Math.sqrt(numPoints));
        Matrix3x3 frame = helper.getRandomFrameAt(queryCap.axis());
        double radius = queryCap.angle().radians();
        double spacing = 2 * radius / sqrtNumPoints;
        List<S2Point> points = Lists.newArrayList();
        for (int i = 0; i < sqrtNumPoints; ++i) {
          for (int j = 0; j < sqrtNumPoints; ++j) {
            points.add(
                S2.fromFrame(
                    frame,
                    S2Point.normalize(
                        new S2Point(
                            Math.tan((i + 0.5) * spacing - radius),
                            Math.tan((j + 0.5) * spacing - radius),
                            1.0))));
          }
        }
        return points;
      }
    };

    /**
     * Returns a list of approximately {@code numPoints} random points sampled from {@code queryCap}
     * by some geometric strategy. Typically the indexed points will occupy some fraction of this
     * cap.)
     */
    protected abstract List<S2Point> createPoints(
        S2ClosestPointQueryTest helper, S2Cap queryCap, int numPoints);
  }

  private interface Target {
    /** Returns the distance from the target to the given point. */
    S1Angle getDistance(S2Point x);
    /** Returns the queue of results, after verifying both methods produce the same queue. */
    List<Result<Integer>> findClosestPoints(S2ClosestPointQuery<Integer> query);
  }

  private static final class PointTarget implements Target {
    final S2Point point;

    PointTarget(S2Point point) {
      this.point = point;
    }

    @Override
    public S1Angle getDistance(S2Point x) {
      return new S1Angle(x, point);
    }

    @Override
    public List<Result<Integer>> findClosestPoints(S2ClosestPointQuery<Integer> query) {
      // Fill 'x' with the query results.
      List<Result<Integer>> x = Lists.newArrayList();
      query.findClosestPoints(x, point);
      // Let the query allocate a queue for us.
      List<Result<Integer>> y = query.findClosestPoints(point);
      // Verify the results are the same and return one of them.
      assertEquals(x, y);
      return y;
    }
  }

  private static final class EdgeTarget implements Target {
    final S2Point a;
    final S2Point b;

    EdgeTarget(S2Point a, S2Point b) {
      this.a = a;
      this.b = b;
    }

    @Override
    public S1Angle getDistance(S2Point x) {
      return S2EdgeUtil.getDistance(x, a, b);
    }

    @Override
    public List<Result<Integer>> findClosestPoints(S2ClosestPointQuery<Integer> query) {
      // Fill 'x' with the query results.
      List<Result<Integer>> x = Lists.newArrayList();
      query.findClosestPointsToEdge(x, a, b);
      // Let the query allocate a queue for us.
      List<Result<Integer>> y = query.findClosestPointsToEdge(a, b);
      // Verify the results are the same and return one of them.
      assertEquals(x, y);
      return y;
    }
  }

  /**
   * Use "query" to find the closest point(s) to the given target, and extract the query results.
   * Also verify that the results satisfy the search criteria.
   */
  private static List<Result<Integer>> getClosestPoints(
      Target target, S2ClosestPointQuery<Integer> query) {
    List<Result<Integer>> actual = target.findClosestPoints(query);
    assertTrue(actual.size() <= query.getMaxPoints());
    if (query.getRegion() != null && query.getMaxDistance() == S1Angle.INFINITY) {
      // We can predict exactly how many points should be returned.
      assertEquals(Math.min(query.getMaxPoints(), query.index().numPoints()), actual.size());
    }
    for (Result<Integer> result : actual) {
      // Check that query.distance() is approximately equal to the angle between point and target.
      // They may be slightly different because query.distance() is computed using S1ChordAngle.
      // Note that the error gets considerably larger (1e-7) as the angle approaches Pi.
      S2Point p = result.entry().point();
      S1Angle angle = result.distance().toAngle();
      assertEquals(target.getDistance(p).radians(), angle.radians(), MAX_CHORD_ANGLE_ERROR);
      // Check that the point satisfies the region() criteria.
      if (query.getRegion() != null) {
        assertTrue(query.getRegion().contains(p));
      }
      // Check that it satisfies the maxDistance() criteria.
      assertTrue(angle.compareTo(query.getMaxDistance()) <= 0);
    }
    return actual;
  }

  private static void checkFindClosestPoints(Target target, S2ClosestPointQuery<Integer> query) {
    query.useBruteForce(true);
    List<Result<Integer>> expected = getClosestPoints(target, query);
    query.useBruteForce(false);
    List<Result<Integer>> actual = getClosestPoints(target, query);
    compareResults(expected, actual, query.getMaxPoints(), query.getMaxDistance(), S1Angle.ZERO);
  }

  /** Validates a variety of random queries. */
  private void checkFactory(PointFactory factory, int numIndexes, int numPoints, int numQueries) {
    S2PointIndex<Integer> index = new S2PointIndex<>();
    S2ClosestPointQuery<Integer> query = new S2ClosestPointQuery<>(index);
    for (int i = 0; i < numIndexes; i++) {
      // Generate a point set and index it.
      S2Cap queryCap = S2Cap.fromAxisAngle(randomPoint(), QUERY_RADIUS);
      index.reset();
      addPoints(index, factory.createPoints(this, queryCap, numPoints));
      query.reset();
      for (int j = 0; j < numQueries; j++) {
        query.setMaxPoints((int) uniform(1, 100));
        if (oneIn(2)) {
          query.setMaxDistance(randomAngle());
        }
        S2LatLngRect rect =
            S2LatLngRect.fromCenterSize(
                new S2LatLng(samplePoint(queryCap)), new S2LatLng(randomAngle(), randomAngle()));
        if (oneIn(5)) {
          query.setRegion(rect);
        }
        if (oneIn(2)) {
          S2Point p = samplePoint(queryCap);
          checkFindClosestPoints(new PointTarget(p), query);
        } else {
          S2Point a = samplePoint(queryCap);
          S2Point b = samplePoint(S2Cap.fromAxisAngle(a, randomDecimilliAngle()));
          checkFindClosestPoints(new EdgeTarget(a, b), query);
        }
      }
    }
  }

  protected static void addPoints(S2PointIndex<Integer> index, List<S2Point> points) {
    for (int i = 0; i < points.size(); i++) {
      index.add(points.get(i), i);
    }
  }

  protected S1Angle randomDecimilliAngle() {
    return S1Angle.radians(Math.pow(1e-4, randomAngle().radians()));
  }

  private S1Angle randomAngle() {
    return S1Angle.radians(uniform(0, QUERY_RADIUS.radians()));
  }
}
