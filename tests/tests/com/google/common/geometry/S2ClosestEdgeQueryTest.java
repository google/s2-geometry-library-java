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

import static com.google.common.geometry.BestEdgesTestUtils.resultsToString;
import static com.google.common.geometry.S2Projections.PROJ;
import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static java.lang.Math.PI;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import com.google.common.collect.ImmutableSet;
import com.google.common.geometry.BestEdgesTestUtils.MinOrMax;
import com.google.common.geometry.S2BestEdgesQueryBase.Result;
import com.google.common.geometry.S2ClosestEdgeQuery.Query;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

/**
 * Tests for S2ClosestEdgeQuery and the underlying S2BestEdgesQueryBase which supplies the
 * interesting parts of the implementation.
 */
public final class S2ClosestEdgeQueryTest extends GeometryTestCase {
  /**
   * For the given S2ClosestEdgeQuery.Target and S2ShapeIndex, visits and returns up to maxShapes
   * S2Shapes from the index containing a point on a connected component of the target.
   */
  public static Set<S2Shape> getContainingShapes(
      S2ClosestEdgeQuery.Target<S1ChordAngle> target, S2ShapeIndex index, int maxShapes) {
    S2ClosestEdgeQuery.Query closestEdgeQuery = S2ClosestEdgeQuery.builder().build(index);

    final Set<S2Shape> visitedShapes = new HashSet<>();
    closestEdgeQuery.visitContainingShapes(
        target,
        (shape) -> {
          visitedShapes.add(shape);
          return visitedShapes.size() < maxShapes;
        });

    assertLessOrEqual(visitedShapes.size(), maxShapes);
    return visitedShapes;
  }

  public void testNoEdges() {
    S2ShapeIndex index = makeIndexOrDie("# #");
    Query query = S2ClosestEdgeQuery.builder().build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.PointTarget<>(new S2Point(1, 0, 0));
    Optional<Result<S1ChordAngle>> result = query.findClosestEdge(target);
    assertFalse(result.isPresent());
  }

  public void testOptionsNotModified() {
    // Tests that findClosestEdge(), getDistance(), and isDistanceLess() do not modify
    // query.options(), even though all of these methods have their own specific options
    // requirements.
    S2ClosestEdgeQuery.Builder builder = S2ClosestEdgeQuery.builder();
    builder
        .setMaxDistance(S1ChordAngle.fromDegrees(3))
        .setMaxError(S1ChordAngle.fromDegrees(0.001))
        .setMaxResults(3);

    // First, verify that the builder is storing the expected values.
    assertEquals(3, builder.maxResults());
    assertExactly(
        S1ChordAngle.fromDegrees(3).getLength2(), builder.maxDistance().getLength2());
    assertExactly(
        S1ChordAngle.fromDegrees(0.001).getLength2(), builder.maxError().getLength2());

    // An index with three points.
    S2ShapeIndex index = makeIndexOrDie("1:1 | 1:2 | 1:3 # #");
    Query query = builder.build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.PointTarget<>(makePointOrDie("2:2"));

    assertEquals(1, query.findClosestEdge(target).get().edgeId());
    assertDoubleNear(1d, query.getDistance(target).degrees(), 1e-15);
    assertTrue(query.isDistanceLess(target, S1ChordAngle.fromDegrees(1.5)));

    // Verify that none of the query options were modified.
    assertEquals(builder.maxResults(), query.options().maxResults());
    assertEquals(builder.maxDistance(), query.options().distanceLimit());
    assertEquals(builder.maxError(), query.options().maxError());
  }

  /** Verifies that visitContainingShapes visits the correct shapes for a point target. */
  public void testVisitContainingShapeForPointTarget() {
    // Builds an index containing 1 point, 1 edge, and 3 polygons, each with a single triangular
    // loop. Only shapes 2 and 4 (the first and third polygons) should contain the target point.
    // (Note that the S2ContainsPointQuery vertex model is SEMI_OPEN, therefore no point or edge
    // contains anything, not even edges containing their endpoints.)
    S2ShapeIndex index =
        makeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.PointTarget<>(makePointOrDie("1:1"));

    // Visit up to one containing shape. Either shape 2 or 4 should be the only one.
    Set<S2Shape> shapeSet = getContainingShapes(target, index, 1);
    assertEquals(1, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // If up to 5 containing shapes may be visited, exactly shapes 2 and 4 actually should be.
    assertEquals(getContainingShapes(target, index, 5),
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)));
  }

  /** Verifies that visitContainingShapes visits the correct shapes for an edge target. */
  public void testVisitContainingShapeForEdgeTarget() {
    // Builds an index containing 1 point, 1 edge, and 3 polygons, each with a single triangular
    // loop. Only shapes 2 and 4 (the first and third polygons) should contain the target edge.
    S2ShapeIndex index =
        makeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
    S2ClosestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.EdgeTarget<>(makePointOrDie("1:2"), makePointOrDie("2:1"));

    // Visit up to one containing shape. Either shape 2 or 4 should be the only one.
    Set<S2Shape> shapeSet = getContainingShapes(target, index, 1);
    assertEquals(1, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // If up to 5 containing shapes may be visited, exactly shapes 2 and 4 actually should be.
    assertEquals(getContainingShapes(target, index, 5),
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)));
  }

  /** Verifies that visitContainingShapes visits the correct shapes for a cell target. */
  public void testVisitContainingShapeForCellTarget() {
    // Builds an index containing 1 point, 1 edge, and 3 polygons, each with a single triangular
    // loop. Only shapes 2 and 4 (the first and third polygons) should contain the very small (leaf)
    // target cell near 1:1.
    S2CellId cellId1 = S2CellId.fromPoint(makePointOrDie("1:1"));
    S2ShapeIndex index =
        makeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
    S2ClosestEdgeQuery.CellTarget<S1ChordAngle> target1 =
        new S2ClosestEdgeQuery.CellTarget<>(new S2Cell(cellId1));

    // Visit up to one containing shape. Either shape 2 or 4 should be the only one.
    Set<S2Shape> shapeSet = getContainingShapes(target1, index, 1);
    assertEquals(1, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // If up to 5 containing shapes may be visited, only shapes 2 and 4 actually should be.
    assertEquals(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)),
        getContainingShapes(target1, index, 5));

    // For a larger cell that properly contains one or more index cells, all shapes that intersect
    // the first such cell in S2CellId order are returned. In the test below, this happens to again
    // be the first and third polygons, whose shape ids are 2 and 4.
    S2CellId cellId2 = cellId1.parent(5);
    S2ClosestEdgeQuery.CellTarget<S1ChordAngle> target2 =
        new S2ClosestEdgeQuery.CellTarget<>(new S2Cell(cellId2));
    assertEquals(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)),
        getContainingShapes(target2, index, 5));
  }

  public void testDistanceEqualToLimit() {
    // Tests the behavior of isDistanceLess, isDistanceLessOrEqual, and
    // isConservativeDistanceLessOrEqual (and the corresponding Options) when the distance to the
    // target exactly equals the chosen limit.
    S2Point p0 = makePointOrDie("23:12");
    S2Point p1 = makePointOrDie("47:11");
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(S2Point.Shape.singleton(p0));
    Query query = S2ClosestEdgeQuery.builder().build(index);

    // Start with two identical points and a zero distance.
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target0 = new S2ClosestEdgeQuery.PointTarget<>(p0);
    S1ChordAngle dist0 = S1ChordAngle.ZERO;
    assertFalse(query.isDistanceLess(target0, dist0));
    assertTrue(query.isDistanceLessOrEqual(target0, dist0));
    assertTrue(query.isConservativeDistanceLessOrEqual(target0, dist0));

    // Now try two points separated by a non-zero distance.
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target1 = new S2ClosestEdgeQuery.PointTarget<>(p1);
    S1ChordAngle dist1 = new S1ChordAngle(p0, p1);
    assertFalse(query.isDistanceLess(target1, dist1));
    assertTrue(query.isDistanceLessOrEqual(target1, dist1));
  }

  @SuppressWarnings("FloatingPointLiteralPrecision")
  public void testTrueDistanceLessThanS1ChordAngleDistance() {
    // Tests that isConservativeDistanceLessOrEqual returns points where the true distance is
    // slightly less than the one computed by S1ChordAngle.
    //
    // The points below had the worst error from among 100,000 random pairs.
    S2Point p0 = new S2Point(0.78516762584829192, -0.50200400690845970, -0.36263449417782678);
    S2Point p1 = new S2Point(0.78563011732429433, -0.50187655940493503, -0.36180828883938054);

    // For the points above, the S1ChordAngle distance is ~4 ulps greater than the true distance.
    S1ChordAngle dist1 = new S1ChordAngle(p0, p1);
    S1ChordAngle limit = dist1.predecessor().predecessor().predecessor().predecessor();
    assertLessThan(S2Predicates.compareDistance(p0, p1, limit.getLength2()), 0);

    // Verify that isConservativeDistanceLessOrEqual() still returns "p1".
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(S2Point.Shape.singleton(p0));
    Query query = S2ClosestEdgeQuery.builder().build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target1 = new S2ClosestEdgeQuery.PointTarget<>(p1);
    assertFalse(query.isDistanceLess(target1, limit));
    assertFalse(query.isDistanceLessOrEqual(target1, limit));
    assertTrue(query.isConservativeDistanceLessOrEqual(target1, limit));
  }

  public void testReuseOfQuery() {
    // Tests that between queries, the internal mechanism for de-duplicating results is re-set.
    S2ShapeIndex index = makeIndexOrDie("2:2 # #");
    Query query = S2ClosestEdgeQuery.builder().setMaxError(S1ChordAngle.fromDegrees(1)).build(index);
    S2ShapeIndex targetIndex = makeIndexOrDie("## 0:0, 0:5, 5:5, 5:0");
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2ClosestEdgeQuery.createShapeIndexTarget(targetIndex);
    List<Result<S1ChordAngle>> results1 = query.findClosestEdges(target);
    List<Result<S1ChordAngle>> results2 = query.findClosestEdges(target);
    assertEquals(results1.size(), results2.size());
  }

  public void testTargetPointInsideIndexedPolygon() {
    // Tests a target point in the interior of an indexed polygon.
    // (The index also includes a polyline with no interior.)
    S2ShapeIndex index = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
    S2ClosestEdgeQuery.Builder builder = S2ClosestEdgeQuery.builder();
    builder
        .setMaxDistance(S1Angle.degrees(1))
        .setIncludeInteriors(true);

    Query query = builder.build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.PointTarget<>(makePointOrDie("2:12"));
    List<Result<S1ChordAngle>> results = query.findClosestEdges(target);
    assertEquals(1, results.size());
    assertEquals(S1ChordAngle.ZERO, results.get(0).distance());
    assertEquals(index.getShapes().get(1), results.get(0).shape());
    assertEquals(-1, results.get(0).edgeId());
    assertTrue(results.get(0).isInterior());
  }

  public void testTargetPointOutsideIndexedPolygon() {
    // Tests a target point partly surrounded by a polyline of three edges around four corners of a
    // square. The index also includes a nearby polygon.
    S2ShapeIndex index = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
    S2ClosestEdgeQuery.Builder builder = S2ClosestEdgeQuery.builder();
    builder
        .setMaxDistance(S1Angle.degrees(1))
        .setIncludeInteriors(true);

    Query query = builder.build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2ClosestEdgeQuery.PointTarget<>(makePointOrDie("2:2"));
    // The target point is two degrees away from the closest edges, so no results are expected.
    List<Result<S1ChordAngle>> results = query.findClosestEdges(target);
    assertEquals(0, results.size());
  }

  public void testTargetPolygonContainingIndexedPoints() {
    // Four query points, where the first two are partly surrounded by the three edges of a target
    // polyline through the corners of a square. The second two are contained by a square polygon
    // target.
    // The point 2:2 is 2 degrees from the 0:0,0:5 polyline edge.
    // The point 3:3 is 2 degrees from the 0:5,5:5 and 5:5,5:0 polyline edges.
    // The point 1:11 is 1 degree from the 0:10, 1:15 edge of the containing polygon.
    // The point 3:13 is 2 degrees from the 5:15, 5:10 edge of the containing polygon.
    S2ShapeIndex queryIndex = makeIndexOrDie("2:2 | 3:3 | 1:11 | 3:13 # #");
    // A three-edge line, and a square polygon.
    S2ShapeIndex targetIndex = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");

    // With a max distance of 1, only one point is within the max distance of an edge, but two of
    // them are contained by a target, so have distance zero.
    Query query = S2ClosestEdgeQuery.builder().setMaxDistance(S1Angle.degrees(1)).build(queryIndex);
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2ClosestEdgeQuery.createShapeIndexTarget(targetIndex);

    // Including interiors for the target will cause the query to find the two query points
    // contained by the target.
    target.setIncludeInteriors(true);
    List<Result<S1ChordAngle>> results = query.findClosestEdges(target);

    // There should be two results at distance zero. They should be edges (points, really) #2 and #3
    // on the same single shape, which is an S2Point.Shape containing multiple points.
    assertEquals(2, results.size());
    assertEquals(S1ChordAngle.ZERO, results.get(0).distance());
    assertEquals(queryIndex.getShapes().get(0), results.get(0).shape());
    assertEquals(2, results.get(0).edgeId()); // 1:11
    assertEquals(S1ChordAngle.ZERO, results.get(1).distance());
    assertEquals(queryIndex.getShapes().get(0), results.get(1).shape());
    assertEquals(3, results.get(1).edgeId()); // 3:13

    // Run the same query again, but with setIncludeInteriors of the target false. This should a
    // find a single index point within the distance limit of 1, at 1:11
    target.setIncludeInteriors(false);
    results = query.findClosestEdges(target);
    assertEquals(resultsToString(results), 1, results.size());
    assertEquals(resultsToString(results), 2, results.get(0).edgeId()); // 1:11

    // Changing maxDistance to zero and keeping includeInteriors false, there should be no results.
    Query query2 = query.toBuilder().setMaxDistance(S1ChordAngle.ZERO).build(queryIndex);
    results = query2.findClosestEdges(target);
    assertEquals(0, results.size());
  }

  public void testEmptyTargetOptimized() {
    // Ensure that the optimized algorithm handles empty targets when a distance limit is specified.
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(
        S2LaxPolygonShape.create(
            new S2Polygon(
                S2Loop.makeRegularLoop(new S2Point(1, 0, 0), S1Angle.radians(0.1), 1000))));
    Query query = S2ClosestEdgeQuery.builder().setMaxDistance(S1Angle.radians(1e-5)).build(index);
    S2ShapeIndex targetIndex = new S2ShapeIndex();
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2ClosestEdgeQuery.createShapeIndexTarget(targetIndex);
    assertEquals(0, query.findClosestEdges(target).size());
  }

  public void testEmptyPolygonTarget() {
    // Verifies that distances are measured correctly to empty polygon targets.
    S2ShapeIndex emptyPolygonIndex = makeIndexOrDie("# # empty");
    S2ShapeIndex fullPolygonIndex = makeIndexOrDie("# # full");

    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> emptyPolygonTarget =
        S2ClosestEdgeQuery.createShapeIndexTarget(emptyPolygonIndex);
    emptyPolygonTarget.setIncludeInteriors(true);

    Query emptyQuery = S2ClosestEdgeQuery.builder().build(emptyPolygonIndex);
    assertEquals(S1ChordAngle.INFINITY, emptyQuery.getDistance(emptyPolygonTarget));

    Query pointQuery = S2ClosestEdgeQuery.builder().build(makeIndexOrDie("1:1 # #"));
    assertEquals(S1ChordAngle.INFINITY, pointQuery.getDistance(emptyPolygonTarget));

    Query fullQuery =
        S2ClosestEdgeQuery.builder().setIncludeInteriors(true).build(fullPolygonIndex);
    assertEquals(S1ChordAngle.INFINITY, fullQuery.getDistance(emptyPolygonTarget));
  }

  public void testFullLaxPolygonTarget() {
    // Verifies that distances are measured correctly to full LaxPolygon targets.
    S2ShapeIndex emptyPolygonIndex = makeIndexOrDie("# # empty");
    S2ShapeIndex pointIndex = makeIndexOrDie("1:1 # #");
    S2ShapeIndex fullPolygonIndex = makeIndexOrDie("# # full");
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2ClosestEdgeQuery.createShapeIndexTarget(fullPolygonIndex);
    target.setIncludeInteriors(true);

    Query emptyQuery = S2ClosestEdgeQuery.builder().build(emptyPolygonIndex);
    assertEquals(S1ChordAngle.INFINITY, emptyQuery.getDistance(target));

    Query pointQuery = S2ClosestEdgeQuery.builder().build(pointIndex);
    assertEquals(S1ChordAngle.ZERO, pointQuery.getDistance(target));

    Query fullQuery =
        S2ClosestEdgeQuery.builder().setIncludeInteriors(true).build(fullPolygonIndex);
    assertEquals(S1ChordAngle.ZERO, fullQuery.getDistance(target));
  }

  public void testIsConservativeDistanceLessOrEqual() {
    int numTested = 0;
    int numConservativeNeeded = 0;

    for (int iter = 0; iter < 1000; ++iter) {
      data.setSeed(iter + 1); // Easier to reproduce a specific case.
      S2Point x = data.getRandomPoint();
      S2Point dir = data.getRandomPoint();
      S1Angle r = S1Angle.radians(PI * pow(1e-30, data.nextDouble()));
      S2Point y = S2EdgeUtil.interpolateAtDistance(r, x, dir);
      S1ChordAngle limit = S1ChordAngle.fromS1Angle(r);
      if (S2Predicates.compareDistance(x, y, limit.getLength2()) <= 0) {
        S2ShapeIndex index = new S2ShapeIndex();
        index.add(S2Point.Shape.singleton(x));
        Query query = S2ClosestEdgeQuery.builder().build(index);
        S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
            new S2ClosestEdgeQuery.PointTarget<>(y);
        assertTrue(query.isConservativeDistanceLessOrEqual(target, limit));
        ++numTested;
        if (!query.isDistanceLess(target, limit)) {
          ++numConservativeNeeded;
        }
      }
    }

    // Verify that in most test cases, the distance between the target points was close to the
    // desired value. Also verify that at least in some test cases, the conservative distance test
    // was actually necessary.
    assertGreaterOrEqual(numTested, 300);
    assertLessOrEqual(numTested, 700);
    // Note that in the C++ unit test, this value is 25, but Java uses a different random sequence.
    assertGreaterOrEqual(numConservativeNeeded, 20);
  }

  //////////////////////////////////////////////////////////////////////////////
  //  General query testing by comparing with brute force method.
  //////////////////////////////////////////////////////////////////////////////

  // The approximate radius of S2Cap from which query edges are chosen.
  static S1Angle kTestCapRadius = kmToAngle(10);

  // An approximate bound on the distance measurement error for "reasonable" distances (say, less
  // than Pi/2) due to using S1ChordAngle.
  static double kTestChordAngleError = 1e-15;

  // Use "query" to find the closest edge(s) to the given target. Verify that the results satisfy
  // the search criteria.
  void getClosestEdges(
      S2ClosestEdgeQuery.Target<S1ChordAngle> target,
      Query query,
      List<Result<S1ChordAngle>> edges) {
    // For half of the tests, have the query return a list of results. For the other half of
    // tests, use the visitor, although this visitor just puts new Results into the list for
    // verification anyway.
    if (data.oneIn(2)) {
      edges.addAll(query.findClosestEdges(target));
    } else {
      query.findClosestEdges(
          target,
          (S1ChordAngle distance, S2Shape shape, int edgeId) ->
              edges.add(new Result<>(distance, shape, edgeId)));
    }

    assertLessOrEqual(edges.size(), query.options().maxResults());
    if (query.options().distanceLimit().isInfinity()) {
      int minExpected = min(query.options().maxResults(), S2ShapeUtil.countEdges(query.index()));
      if (!query.options().includeInteriors()) {
        // We can predict exactly how many edges should be returned.
        assertEquals(minExpected, edges.size());
      } else {
        // All edges should be returned, and possibly some shape interiors.
        assertLessOrEqual(minExpected, edges.size());
      }
    }
    for (Result<S1ChordAngle> edge : edges) {
      // Check that the edge satisfies the maxDistance() condition.
      assertLessThan(edge.distance(), query.options().distanceLimit());
    }
  }

  /**
   * Runs the closest edge query in both brute force and optimized modes, and compares results.
   * Verifies upper and lower bounds on the closest optimized result w.r.t. maxError. Returns the
   * true closest result.
   */
  @CanIgnoreReturnValue
  Result<S1ChordAngle> testFindClosestEdges(
      S2ShapeIndex index,
      S2ClosestEdgeQuery.Target<S1ChordAngle> target,
      S2ClosestEdgeQuery.Builder queryBuilder) {
    List<Result<S1ChordAngle>> bruteForceResults = new ArrayList<>();
    List<Result<S1ChordAngle>> optimizedResults = new ArrayList<>();

    // Get the closest edges using brute force, and verify w.r.t the search criteria.
    queryBuilder.setUseBruteForce(true);
    S2ClosestEdgeQuery.Query bruteForceQuery = queryBuilder.build(index);
    getClosestEdges(target, bruteForceQuery, bruteForceResults);

    // Get the closest edges using the optimized algorithm, and verify w.r.t the search criteria.
    queryBuilder.setUseBruteForce(false);
    S2ClosestEdgeQuery.Query optimizedQuery = queryBuilder.build(index);
    getClosestEdges(target, optimizedQuery, optimizedResults);

    assertTrue(BestEdgesTestUtils.resultsAreEquivalent(
        MinOrMax.MIN,
        bruteForceResults,
        optimizedResults,
        queryBuilder.maxResults(),
        queryBuilder.distanceLimit(),
        queryBuilder.maxError()));

    if (bruteForceResults.isEmpty()) {
      return new Result<>(S1ChordAngle.INFINITY, null, -1);
    }

    // Note that when options.maxError() > 0, expected.get(0).distance() may not be the minimum
    // distance. It is never larger by more than maxError(), but the actual value also depends on
    // maxResults().
    //
    // Here we verify that getDistance() and isDistanceLess() return results that are consistent
    // with the maxError() setting.
    S1ChordAngle maxError = queryBuilder.maxError();
    // The distance of the truly closest result.
    S1ChordAngle bruteForceMinDistance = bruteForceResults.get(0).distance();
    // The worst acceptable result distance, if using maxError.
    S1ChordAngle upperBound = S1ChordAngle.add(bruteForceMinDistance, maxError);
    S1ChordAngle lowerBound = S1ChordAngle.sub(bruteForceMinDistance, maxError);

    assertLessOrEqual(optimizedQuery.getDistance(target), upperBound);

    // Test isDistanceLess().
    assertFalse(optimizedQuery.isDistanceLess(target, lowerBound));
    assertTrue(optimizedQuery.isConservativeDistanceLessOrEqual(target, bruteForceMinDistance));

    // Return the closest edge result so that we can also test project().
    return bruteForceResults.get(0);
  }

  /**
   * Does 'numQueries' ClosestEdgeQuery tests. Each test runs a query against a randomly selected
   * index from one of 'numIndexes' randomly generated ShapeIndexes, each of which contains one
   * pseudo-random shape generated by the given ShapeFactory. The ShapeFactory is initialized from
   * the TestDataGenerator's internal Random. For each query, random parameters and a random target
   * are used, also using the TestDataGenerator's internal Random.
   *
   * <p>To obtain repeatable test sequences, callers can use data.setSeed() before calling
   * testWithShapeFactory().
   *
   * <p>The running time of this test is proportional to (numIndexes + numQueries) * numEdges.
   * (Note that every query is checked using the brute force algorithm.)
   */
  void testWithShapeFactory(
      TestDataGenerator.ShapeFactory factory,
      int numIndexes,
      int numEdges,
      int numQueries,
      boolean checkConservativeCellDistance) {
    // Build a set of S2ShapeIndexes containing a pseudo-random shape from the provided factory.
    List<S2Cap> indexCaps = new ArrayList<>();
    List<S2ShapeIndex> indexes = new ArrayList<>();

    // Test indexes are placed at random locations across the whole sphere.
    for (int i = 0; i < numIndexes; ++i) {
      S2Cap cap = S2Cap.fromAxisAngle(data.getRandomPoint(), kTestCapRadius);
      indexCaps.add(cap);
      indexes.add(factory.getShape(cap, numEdges));
    }

    for (int i = 0; i < numQueries; ++i) {
      // Select an index at random.
      int iIndex = data.uniform(numIndexes);
      S2Cap indexCap = indexCaps.get(iIndex);
      S2ShapeIndex index = indexes.get(iIndex);

      // Choose query points from an area approximately 4x larger than the geometry being tested.
      S1Angle queryRadius = indexCap.angle().mul(2);
      S2Cap queryCap = S2Cap.fromAxisAngle(indexCap.axis(), queryRadius);
      S2ClosestEdgeQuery.Builder builder = S2ClosestEdgeQuery.builder();

      // Occasionally we don't set any limit on the number of result edges. (This may return all
      // edges if we also don't set a distance limit.)
      if (!data.oneIn(5)) {
        builder.setMaxResults(1 + data.uniform(10));
      }
      // We set a distance limit 2/3 of the time.
      if (!data.oneIn(3)) {
        builder.setMaxDistance(queryRadius.mul(data.nextDouble()));
      }
      if (checkConservativeCellDistance || data.oneIn(2)) {
        // Choose a maximum error whose logarithm is uniformly distributed over a reasonable range,
        // except that it is sometimes zero.
        builder.setMaxError(
            S1ChordAngle.fromRadians(pow(1e-4, data.nextDouble()) * queryRadius.radians()));
      }
      builder.setIncludeInteriors(data.oneIn(2));

      int targetType = data.uniform(4);
      // For checking conservative cell distance, we need a target that uses maxError.
      if (checkConservativeCellDistance) {
        targetType = 3;
      }
      if (targetType == 0) {
        // Find the edges closest to a given point.
        S2Point point = data.samplePoint(queryCap);
        S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target =
            new S2ClosestEdgeQuery.PointTarget<>(point);
        Result<S1ChordAngle> closest = testFindClosestEdges(index, target, builder);

        if (!closest.distance().isInfinity()) {
          // Also test the Project method.
          Query query = builder.build(index);
          Optional<Result<S1ChordAngle>> toProject = query.findClosestEdge(target);
          assertTrue(toProject.isPresent());
          assertDoubleNear(
              toProject.get().distance().toAngle().radians(),
              new S1Angle(point, query.project(point, toProject.get())).radians(),
              kTestChordAngleError);
        }
      } else if (targetType == 1) {
        // Find the edges closest to a given edge.
        S2Point a = data.samplePoint(queryCap);
        S2Point b =
            data.samplePoint(S2Cap.fromAxisAngle(a, queryRadius.mul(pow(1e-4, data.nextDouble()))));
        S2ClosestEdgeQuery.EdgeTarget<S1ChordAngle> target =
            new S2ClosestEdgeQuery.EdgeTarget<>(a, b);
        testFindClosestEdges(index, target, builder);
      } else if (targetType == 2) {
        // Find the edges closest to a given cell.
        int minLevel = PROJ.maxDiag.getMinLevel(queryRadius.radians());
        int level = minLevel + data.uniform(S2CellId.MAX_LEVEL - minLevel + 1);
        S2Point a = data.samplePoint(queryCap);
        S2CellId targetId = S2CellId.fromPoint(a).parent(level);
        S2Cell cell = new S2Cell(targetId);
        S2ClosestEdgeQuery.CellTarget<S1ChordAngle> target =
            new S2ClosestEdgeQuery.CellTarget<>(cell);
        testFindClosestEdges(index, target, builder);
      } else {
        assertEquals(3, targetType);
        // Use another one of the pre-built indexes as the target.
        int jIndex = data.uniform(numIndexes);
        S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
            S2ClosestEdgeQuery.createShapeIndexTarget(indexes.get(jIndex));
        target.setIncludeInteriors(data.oneIn(2));
        testFindClosestEdges(index, target, builder);
      }
    }
  }

  static int kNumIndexes = 50;
  static int kNumEdges = 100;
  static int kNumQueries = 200;

  public void testCircleEdges() {
    testWithShapeFactory(
        data.new RegularLoopShapeFactory(), kNumIndexes, kNumEdges, kNumQueries, false);
  }

  public void testFractalEdges() {
    testWithShapeFactory(
        data.new FractalLoopShapeFactory(), kNumIndexes, kNumEdges, kNumQueries, false);
  }

  public void testPointCloudEdges() {
    testWithShapeFactory(
        data.new PointCloudShapeFactory(), kNumIndexes, kNumEdges, kNumQueries, false);
  }

  public void testConservativeCellDistanceIsUsed() {
    // These specific test cases happen to fail if maxError() is not properly taken into account
    // when measuring distances to S2ShapeIndex cells.
    int[] testSeeds = {126, 397, 564 };
    for (int seed : testSeeds) {
      data.setSeed(seed);
      testWithShapeFactory(data.new FractalLoopShapeFactory(), 5, 100, 10, true);
    }
  }
}
