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

import static com.google.common.geometry.S2Projections.MAX_DIAG;
import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.parsePointsOrDie;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.geometry.BestEdgesTestUtils.MinOrMax;
import com.google.common.geometry.S2BestEdgesQueryBase.Result;
import com.google.common.geometry.S2FurthestEdgeQuery.Query;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Tests for S2FurthestEdgeQuery and the underlying S2BestEdgesQueryBase which supplies the
 * interesting parts of the implementation.
 */
@RunWith(JUnit4.class)
public final class S2FurthestEdgeQueryTest extends GeometryTestCase {
  /**
   * For the given S2FurthestEdgeQuery.Target and S2ShapeIndex, visits and returns up to maxShapes
   * S2Shapes from the index containing a point antipodal to a connected component of the target.
   */
  public static Set<S2Shape> getAntipodalShapes(
      S2FurthestEdgeQuery.Target<S1ChordAngle> target, S2ShapeIndex index, int maxShapes) {
    S2FurthestEdgeQuery.Query furthestEdgeQuery = S2FurthestEdgeQuery.builder().build(index);

    // Visits shapes that contain a point antipodal to a connected component of the target.
    final Set<S2Shape> visitedShapes = new HashSet<>();
    furthestEdgeQuery.visitAntipodalShapes(
        target,
        shapeId -> {
          visitedShapes.add(index.getShapes().get(shapeId));
          return visitedShapes.size() < maxShapes;
        });

    assertLessOrEqual(visitedShapes.size(), maxShapes);
    return visitedShapes;
  }

  @Test
  public void testNoEdges() {
    S2ShapeIndex index = makeIndexOrDie("# #");
    Query query = S2FurthestEdgeQuery.builder().build(index);
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(new S2Point(1, 0, 0));
    Optional<Result<S1ChordAngle>> result = query.findFurthestEdge(target);
    assertFalse(result.isPresent());
  }

  @Test
  public void testOptionsNotModified() {
    // Tests that findFurthestEdge(), getDistance(), and isDistanceGreater() do not modify
    // query.options(), even though all of these methods have their own specific options
    // requirements.
    S2FurthestEdgeQuery.Builder builder = S2FurthestEdgeQuery.builder()
        .setMinDistance(S1ChordAngle.fromDegrees(1))
        .setMaxError(S1ChordAngle.fromDegrees(0.001))
        .setMaxResults(3);

    // First, verify that the builder is storing the expected values.
    assertEquals(3, builder.maxResults());
    assertExactly(
        S1ChordAngle.fromDegrees(1).getLength2(), builder.minDistance().getLength2());
    assertExactly(
        S1ChordAngle.fromDegrees(0.001).getLength2(), builder.maxError().getLength2());

    S2ShapeIndex index = makeIndexOrDie("0:1 | 0:2 | 0:3 # #");
    Query query = builder.build(index);
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(makePointOrDie("0:4"));

    assertEquals(0, query.findFurthestEdge(target).get().edgeId());
    assertDoubleNear(3.0, query.getDistance(target).degrees(), 1e-15);
    assertTrue(query.isDistanceGreater(target, S1ChordAngle.fromDegrees(1.5)));

    // Verify that none of the query options were modified.
    assertEquals(builder.maxResults(), query.options().maxResults());
    assertEquals(builder.minDistance(), query.options().distanceLimit());
    assertEquals(builder.maxError(), query.options().maxError());
  }

  /** Verifies that visitAntipodalShapes visits the correct shapes for a point target. */
  @Test
  public void testVisitAntipodalShapesForPointTarget() {
    // Builds an index containing 1 point, 1 edge, and 3 polygons, each with a single triangular
    // loop. Only shapes 2 and 4 (the first and third polygons) are at maximum distance to the
    // target.
    // (Note that the S2ContainsPointQuery vertex model is SEMI_OPEN, therefore no point or edge
    // contains anything, not even edges containing their endpoints.)
    S2ShapeIndex index =
        makeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
    S2Point p = makePointOrDie("1:1");
    // Make the target antipodal to 1:1 so the polygons at max distance are those that contain 1:1.
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(p.neg());

    // Visit up to one shape at best (maximum) distance from the target. Either shape 2 or 4 should
    // be the only one.
    Set<S2Shape> shapeSet = getAntipodalShapes(target, index, 1);
    assertEquals(1, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // If up to 5 antipodal shapes may be visited, exactly shapes 2 and 4 actually should be.
    shapeSet = getAntipodalShapes(target, index, 5);
    assertEquals(2, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));
  }

  /** Verifies that visitAntipodalShapes visits the correct shapes for an edge target. */
  @Test
  public void testVisitAntipodalShapesForEdgeTarget() {
    // Builds an index containing 1 point, 1 edge, and 3 polygons, each with a single triangular
    // loop. Only shapes 2 and 4 (the first and third polygons) should be antipodal to the target.
    S2ShapeIndex index =
        makeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");

    // Make the target antipodal to the 1:1,2:1 edge, so the polygons at max distance are the
    // two triangles above.
    List<S2Point> edge = parsePointsOrDie("1:2, 2:1");
    S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.EdgeTarget<>(edge.get(0).neg(), edge.get(1).neg());

    // Visit up to one antipodal shape. Either shape 2 or 4 should be the only one.
    Set<S2Shape> shapeSet = getAntipodalShapes(target, index, 1);
    assertEquals(1, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // If up to 5 antipodal shapes may be visited, exactly shapes 2 and 4 actually should be.
    shapeSet = getAntipodalShapes(target, index, 5);
    assertEquals(2, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));
  }

  /** Verifies that visitAntipodalShapes visits the correct shapes for cell target. */
  @Test
  public void testVisitAntipodalShapesForCellTarget() {
    // Builds an index containing 1 point, 1 edge, and 3 polygons, each with a single triangular
    // loop. Only shapes 2 and 4 (the first and third polygons) should contain the very small (leaf)
    // target cell near the antipode of 1:1.
    S2ShapeIndex index =
        makeIndexOrDie("1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");

    S2CellId cellId1 = S2CellId.fromPoint(makePointOrDie("1:1").neg());
    S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target1 =
        new S2FurthestEdgeQuery.CellTarget<>(new S2Cell(cellId1));

    // Visit up to one antipodal shape. Either shape 2 or 4 should be the only one.
    Set<S2Shape> shapeSet = getAntipodalShapes(target1, index, 1);
    assertEquals(1, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // If up to 5 antipodal shapes may be visited, only shapes 2 and 4 actually should be.
    shapeSet = getAntipodalShapes(target1, index, 5);
    assertEquals(2, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));

    // For a larger cell that properly contains the antipodal reflection of one or more index cells,
    // all shapes that intersect the first such cell in S2CellId order are returned. In the test
    // below, this happens to again be the first and third polygons, whose shape ids are 2 and 4.
    S2CellId cellId2 = cellId1.parent(5);
    S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target2 =
        new S2FurthestEdgeQuery.CellTarget<>(new S2Cell(cellId2));
    shapeSet = getAntipodalShapes(target2, index, 5);
    assertEquals(2, shapeSet.size());
    assertTrue(
        ImmutableSet.of(index.getShapes().get(2), index.getShapes().get(4)).containsAll(shapeSet));
  }

  // In furthest edge queries, the following distance computation is used when updating max
  // distances.
  S1ChordAngle getMaxDistanceToEdge(S2Point x, S2Point y0, S2Point y1) {
    return S2EdgeUtil.updateMaxDistance(x, y0, y1, S1ChordAngle.ZERO);
  }

  @Test
  public void testDistanceEqualToLimit() {
    // Tests the behavior of isDistanceGreater, isDistanceGreaterOrEqual, and
    // isConservativeDistanceGreaterOrEqual (and the corresponding Options) when the distance to the
    // target exactly equals the chosen limit.
    S2Point p0 = makePointOrDie("23:12");
    S2Point p1 = makePointOrDie("47:11");
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(S2Point.Shape.singleton(p0));
    Query query = S2FurthestEdgeQuery.builder().build(index);

    // Start with antipodal points and a maximum (180 degrees) distance.
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target0 =
        new S2FurthestEdgeQuery.PointTarget<>(p0.neg());
    S1ChordAngle distMax = S1ChordAngle.STRAIGHT;
    assertFalse(query.isDistanceGreater(target0, distMax));
    assertTrue(query.isDistanceGreaterOrEqual(target0, distMax));
    assertTrue(query.isConservativeDistanceGreaterOrEqual(target0, distMax));

    // Now try two points separated by a non-maximal distance.
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target1 =
        new S2FurthestEdgeQuery.PointTarget<>(p1.neg());
    S1ChordAngle dist1 = getMaxDistanceToEdge(p0, p1.neg(), p1.neg());
    assertFalse("dist1 = " + dist1, query.isDistanceGreater(target1, dist1));
    assertTrue(query.isDistanceGreaterOrEqual(target1, dist1));
    assertTrue(query.isConservativeDistanceGreaterOrEqual(target1, dist1));
  }

  @SuppressWarnings("FloatingPointLiteralPrecision")
  @Test
  public void testTrueDistanceGreaterThanS1ChordAngleDistance() {
    // Tests that isConservativeDistanceGreaterOrEqual returns points where the true distance is
    // slightly greater than the one computed by S1ChordAngle.
    //
    // The points below had the worst error from among 1x10^6 random pairs.
    S2Point p0 = new S2Point(0.72362949088190598, -0.39019820403414807, -0.56930283812266336);
    S2Point p1 = new S2Point(0.54383822931548842, 0.758981734255934404, 0.35803171284238039);

    // Verify that the S1ChordAngle distance is ~3 ulps greater than the true distance.
    S1ChordAngle dist1 = getMaxDistanceToEdge(p0, p1, p1);
    S1ChordAngle limit = dist1.successor().successor().successor();
    assertGreaterThan(S2Predicates.compareDistance(p0, p1, limit.getLength2()), 0);

    // Verify that isConservativeDistanceGreaterOrEqual() still returns "p1".
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(S2Point.Shape.singleton(p0));
    Query query = S2FurthestEdgeQuery.builder().build(index);
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target1 =
        new S2FurthestEdgeQuery.PointTarget<>(p1);
    assertFalse(query.isDistanceGreater(target1, limit));
    assertFalse(query.isDistanceGreaterOrEqual(target1, limit));
    assertTrue(query.isConservativeDistanceGreaterOrEqual(target1, limit));
  }

  @Test
  public void testAntipodalPointInsideIndexedPolygon() {
    // Tests a target point antipodal to the interior of an indexed polygon.
    // (The index also includes a polyline loop with no interior.)
    S2ShapeIndex index = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
    S2FurthestEdgeQuery.Builder builder = S2FurthestEdgeQuery.builder()
        .setMinDistance(S1Angle.degrees(178))
        .setIncludeInteriors(true);

    // First check that with includeInteriors set to true, the distance is maximum (180 degrees).
    Query query = builder.build(index);
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(makePointOrDie("2:12").neg());
    List<Result<S1ChordAngle>> results = query.findFurthestEdges(target);
    assertGreaterThan(results.size(), 0);
    assertEquals(S1ChordAngle.STRAIGHT, results.get(0).distance());
    assertEquals(1, results.get(0).shapeId());
    assertEquals(-1, results.get(0).edgeId());

    // Next check that with includeInteriors set to false, the distance is less than 180 for the
    // same target and index.
    Query query2 = query.toBuilder().setIncludeInteriors(false).build(index);
    results = query2.findFurthestEdges(target);
    assertGreaterThan(results.size(), 0);
    assertLessOrEqual(results.get(0).distance(), S1ChordAngle.STRAIGHT);
    // The same polygon shape, with id 1.
    assertEquals(1, results.get(0).shapeId());
    // Found a specific edge, so id should be positive.
    assertGreaterThan(results.get(0).edgeId(), 0);
  }

  @Test
  public void testAntipodalPointOutsideIndexedPolygon() {
    // Tests a target point antipodal to the interior of a polyline loop with no interior. The index
    // also includes a polygon almost antipodal to the target, but with all edges closer than the
    // minDistance threshold.
    S2ShapeIndex index = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
    S2FurthestEdgeQuery.Builder builder = S2FurthestEdgeQuery.builder()
        .setMinDistance(S1Angle.degrees(179))
        .setIncludeInteriors(true);

    Query query = builder.build(index);
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(makePointOrDie("2:2").neg());
    List<Result<S1ChordAngle>> results = query.findFurthestEdges(target);
    assertEquals(0, results.size());
  }

  @Test
  public void testTargetPolygonContainingIndexedPoints() {
    // Two points are contained within a target polyline loop (no interior) and two points are
    // contained within a target polygon.
    S2ShapeIndex index = makeIndexOrDie("2:2 | 4:4 | 1:11 | 3:12 # #");
    Query query = S2FurthestEdgeQuery.builder().setUseBruteForce(false).build(index);
    S2ShapeIndex targetIndex = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
    S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2FurthestEdgeQuery.createShapeIndexTarget(targetIndex);
    target.setIncludeInteriors(true);

    List<Result<S1ChordAngle>> results1 = query.findFurthestEdges(target);
    // All points should be returned since we did not specify maxResults.
    assertEquals(4, results1.size());
    assertGreaterThan(results1.get(0).distance(), S1ChordAngle.ZERO);
    assertEquals(0, results1.get(0).shapeId());
    assertEquals(0, results1.get(0).edgeId()); // 2:2 (to 5:15)
    assertGreaterThan(results1.get(1).distance(), S1ChordAngle.ZERO);
    assertEquals(0, results1.get(1).shapeId());
    assertEquals(3, results1.get(1).edgeId()); // 3:12 (to 0:0)
  }

  @Test
  public void testAntipodalPolygonContainingIndexedPoints() {
    // Two antipodal points are contained within a polyline loop (no interior)
    // and two antipodal points are contained within a polygon.
    List<S2Point> points = parsePointsOrDie("2:2, 3:3, 1:11, 3:13");
    S2ShapeIndex index = new S2ShapeIndex();
    List<S2Point> antipodalPoints = new ArrayList<>();
    for (S2Point p : points) {
      antipodalPoints.add(p.neg());
    }
    index.add(S2Point.Shape.fromList(antipodalPoints));

    // With a minimum distance of 179, only the antipodal points contained by the polygon will be
    // found, at distance STRAIGHT.
    Query query = S2FurthestEdgeQuery.builder().setMinDistance(S1Angle.degrees(179)).build(index);
    S2ShapeIndex targetIndex = makeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
    S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2FurthestEdgeQuery.createShapeIndexTarget(targetIndex);
    target.setIncludeInteriors(true);
    List<Result<S1ChordAngle>> results = query.findFurthestEdges(target);
    assertEquals(2, results.size());

    assertEquals(S1ChordAngle.STRAIGHT, results.get(0).distance());
    // The two points could be returned in either order.
    boolean foundShape0Edge2 = false;
    boolean foundShape0Edge3 = false;
    for (Result<S1ChordAngle> result : results) {
      assertEquals(0, result.shapeId());
      assertEquals(S1ChordAngle.STRAIGHT, result.distance());
      if (result.edgeId() == 2) { // 1:11
        assertFalse(foundShape0Edge2);
        foundShape0Edge2 = true;
      }
      if (result.edgeId() == 3) { // 3:13
        assertFalse(foundShape0Edge3);
        foundShape0Edge3 = true;
      }
    }
    assertTrue(foundShape0Edge2);
    assertTrue(foundShape0Edge3);
  }

  @Test
  public void testEmptyPolygonTarget() {
    // Verifies that distances are measured correctly to empty polygon targets.
    S2ShapeIndex emptyPolygonIndex = makeIndexOrDie("# # empty");
    S2ShapeIndex fullPolygonIndex = makeIndexOrDie("# # full");

    S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> emptyPolygonTarget =
        S2FurthestEdgeQuery.createShapeIndexTarget(emptyPolygonIndex);
    emptyPolygonTarget.setIncludeInteriors(true);

    Query pointIndexQuery = S2FurthestEdgeQuery.builder().build(makeIndexOrDie("1:1 # #"));
    assertEquals(S1ChordAngle.NEGATIVE, pointIndexQuery.getDistance(emptyPolygonTarget));

    Query emptyIndexQuery = S2FurthestEdgeQuery.builder().build(emptyPolygonIndex);
    assertEquals(S1ChordAngle.NEGATIVE, emptyIndexQuery.getDistance(emptyPolygonTarget));

    Query fullIndexQuery =
        S2FurthestEdgeQuery.builder().setIncludeInteriors(true).build(fullPolygonIndex);
    assertEquals(S1ChordAngle.NEGATIVE, fullIndexQuery.getDistance(emptyPolygonTarget));
  }

  @Test
  public void testFullLaxPolygonTarget() {
    // Verifies that distances are measured correctly to full LaxPolygon targets.
    S2ShapeIndex emptyPolygonIndex = makeIndexOrDie("# # empty");
    S2ShapeIndex pointIndex = makeIndexOrDie("1:1 # #");
    S2ShapeIndex fullPolygonIndex = makeIndexOrDie("# # full");
    S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
        S2FurthestEdgeQuery.createShapeIndexTarget(fullPolygonIndex);
    target.setIncludeInteriors(true);

    Query emptyQuery = S2FurthestEdgeQuery.builder().build(emptyPolygonIndex);
    assertEquals(S1ChordAngle.NEGATIVE, emptyQuery.getDistance(target));

    Query pointQuery = S2FurthestEdgeQuery.builder().build(pointIndex);
    assertEquals(S1ChordAngle.STRAIGHT, pointQuery.getDistance(target));

    Query fullQuery =
        S2FurthestEdgeQuery.builder().setIncludeInteriors(true).build(fullPolygonIndex);
    assertEquals(S1ChordAngle.STRAIGHT, fullQuery.getDistance(target));
  }

  //////////////////////////////////////////////////////////////////////////////
  //  General query testing by comparing with brute force method.
  //////////////////////////////////////////////////////////////////////////////

  // The approximate radius of S2Cap from which query edges are chosen.
  static S1Angle kTestCapRadius = kmToAngle(10);

  // Use "query" to find the furthest edge(s) to the given target. Verify that the results satisfy
  // the search criteria.
  void getFurthestEdges(
      S2FurthestEdgeQuery.Target<S1ChordAngle> target,
      Query query,
      List<Result<S1ChordAngle>> edges) {
    // For half of the tests, have the query return a list of results. For the other half of
    // tests, use the visitor, although this visitor just puts new Results into the list for
    // verification anyway.
    if (data.oneIn(2)) {
      edges.addAll(query.findFurthestEdges(target));
    } else {
      query.findFurthestEdges(target,
          (distance, shapeId, edgeId) -> edges.add(new Result<>(distance, shapeId, edgeId)));
    }

    assertLessOrEqual(edges.size(), query.options().maxResults());
    if (query.options().distanceLimit().isNegative()) {
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
      // Check that the edge satisfies the minDistance() condition.
      assertGreaterOrEqual(edge.distance(), query.options().distanceLimit());
    }
  }

  void testFindFurthestEdges(
      S2ShapeIndex index,
      S2FurthestEdgeQuery.Target<S1ChordAngle> target,
      S2FurthestEdgeQuery.Builder queryBuilder) {
    List<Result<S1ChordAngle>> bruteForceResults = new ArrayList<>();
    List<Result<S1ChordAngle>> optimizedResults = new ArrayList<>();

    queryBuilder.setUseBruteForce(true);
    S2FurthestEdgeQuery.Query bruteForceQuery = queryBuilder.build(index);
    getFurthestEdges(target, bruteForceQuery, bruteForceResults);

    queryBuilder.setUseBruteForce(true);
    S2FurthestEdgeQuery.Query optimizedQuery = queryBuilder.build(index);
    getFurthestEdges(target, optimizedQuery, optimizedResults);

     assertTrue(BestEdgesTestUtils.resultsAreEquivalent(
        MinOrMax.MAX,
        bruteForceResults,
        optimizedResults,
        queryBuilder.maxResults(),
        queryBuilder.distanceLimit(),
        queryBuilder.maxError()));

    if (bruteForceResults.isEmpty()) {
      return;
    }

    // Note that when options.maxError() > 0, expected.get(0).distance may not be the maximum
    // distance. It is never smaller by more than maxError(), but the actual value also depends on
    // maxResults().
    //
    // Here we verify that getDistance() and isDistanceGreater() return results that are consistent
    // with the maxError() setting.
    S1ChordAngle maxError = queryBuilder.maxError();

    // The distance of the truly furthest result.
    S1ChordAngle bruteForceMaxDistance = bruteForceResults.get(0).distance();
    // The closest acceptable result distance, if using maxError.
    S1ChordAngle lowerBound = S1ChordAngle.sub(bruteForceMaxDistance, maxError);
    // isDistanceGreater may not return a distance
    S1ChordAngle upperBound = S1ChordAngle.add(bruteForceMaxDistance, maxError);

    // getDistance() is permitted to underestimate by as much as maxError, but no more.
    assertGreaterOrEqual(optimizedQuery.getDistance(target), lowerBound);

    // Test isDistanceGreater().
    assertFalse(optimizedQuery.isDistanceGreater(target, upperBound));
    assertTrue(optimizedQuery.isDistanceGreater(target, bruteForceMaxDistance.predecessor()));
    assertTrue(optimizedQuery.isConservativeDistanceGreaterOrEqual(target, bruteForceMaxDistance));
  }

  /**
   * Puts shape A in a shape index, edge B in an edge target, and uses a FurthestEdge query to get
   * the maximum distance between them, after setting 'includeInteriors' for shape A on the query.
   * Compares the distance to "expectedDistanceRadians".
   */
  private void checkShapeToEdgeMaxDistance(
    // TODO(user): Port checkShapeToEdgeMaxDistance and the unit tests that use it to C++.
      S2Shape a, boolean includeInteriors, S2Point b0, S2Point b1, double expectedDistanceRadians) {
    b0 = b0.normalize();
    b1 = b1.normalize();

    S2ShapeIndex index = new S2ShapeIndex();
    index.add(a);

    S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.EdgeTarget<>(b0, b1);

    // Get the single furthest edge.
    S2FurthestEdgeQuery.Builder builder = S2FurthestEdgeQuery.builder();
    builder.setIncludeInteriors(includeInteriors);
    Query query = builder.build(index);
    Optional<Result<S1ChordAngle>> opt = query.findFurthestEdge(target);
    assertTrue(opt.isPresent());
    Result<S1ChordAngle> result = opt.get();

    // Get the indexed edge that's farthest from the target.
    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
    index.getShapes().get(result.shapeId()).getEdge(result.edgeId(), edge);

    assertEquals("Furthest Indexed Edge from target is " + edge.toDegreesString(),
        expectedDistanceRadians, result.distance().radians(), 1e-15);
  }

  /**
   * Puts edge A in a shape index, edge B in an edge target, and uses a FurthestEdge query to get
   * the maximum distance between them. Compares to "expectedDistanceRadians".
   */
  private void checkEdgePairMaxDistance(
      S2Point a0, S2Point a1, S2Point b0, S2Point b1, double expectedDistanceRadians) {
    checkShapeToEdgeMaxDistance(
        new S2Edge(a0.normalize(), a1.normalize()), false, b0, b1, expectedDistanceRadians);
  }

  /**
   * Exactly the same test cases as {@link S2EdgeUtilTest#testEdgePairMaxDistanceSimple()} but using
   * FurthestEdgeQuery. This demonstrates that for simple edge-to-edge queries, the maximum distance
   * between the index and the target is the maximum distance between the edges.
   */
  @Test
  public void testFurthestEdgeDistanceEdgePairs() {
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

  /** Tests simple cases of maximum distance from an edge to a polygon. */
  @Test
  public void testFurthestEdgeDistanceEdgeToPolygon() {
    S2Point center = ll(12, 12);
    // Four points around the center
    S2Point ne = ll(13, 13);
    S2Point nw = ll(13, 11);
    S2Point sw = ll(11, 11);
    S2Point se = ll(11, 13);
    // A test square, two degrees on a side, centered at 'center'. (The interior of the polygon is
    // on the left when walking the vertices in the given order).
    S2LaxPolygonShape square = S2LaxPolygonShape.fromLoop(ImmutableList.of(ne, nw, sw, se));
    // The same square, but reversed, so representing the sphere with a square hole cut out of it.
    S2LaxPolygonShape hole = S2LaxPolygonShape.fromLoop(ImmutableList.of(se, sw, nw, ne));

    // A point just inside the square, 0.1 degrees from the west edge (i.e. (13,11) to (11,11))
    S2Point inside = ll(12, 11.1);

    // Two points outside the square, half a degree to the west and east.
    S2Point west = ll(12, 10.5);
    S2Point east = ll(12, 13.5);

    // Check the maximum distance between the square and an edge inside the square, with one edge
    // endpoint at the center of the square and one endpoint close to a side of the square. The
    // maximum distance is from the endpoint near the western edge of the square, to the south east
    // corner of the square.
    checkShapeToEdgeMaxDistance(square, false, inside, center, inside.angle(se));
    // The edge same inside the square hole has the same maximum distance.
    checkShapeToEdgeMaxDistance(hole, false, inside, center, inside.angle(se));

    // Check the maximum distance between the square and an edge crossing one square edge, with one
    // edge endpoint half a degree outside the square on the west, and one endpoint close to the
    // west edge inside the square. The maximum distance is
    // attained between the outside endpoint at 'west' and the opposite (east) side of the square:
    // the south east corner specifically.
    checkShapeToEdgeMaxDistance(square, false, west, inside, west.angle(se));
    checkShapeToEdgeMaxDistance(hole, false, west, inside, west.angle(se));

    // Check the maximum distance between the square and an edge crossing two edges of the square,
    // so both endpoints are outside the square. In some intuitive sense, the maximum distance
    // between the edge and the square is between the center of the edge and a corner of the square,
    // but that is the point on the edge furthest from ALL points on the square. The query returns
    // the maximum distance between the edge and ANY point on the square, which is the same as the
    // previous case: the distance from an edge endpoint and a corner of the square on the opposite
    // side.
    checkShapeToEdgeMaxDistance(square, false, east, west, west.angle(se));
    // Here, both edge endpoints are inside the polygon, but it makes no difference.
    checkShapeToEdgeMaxDistance(hole, false, east, west, west.angle(se));
  }

  /**
   * Does 'numQueries' FurthestEdgeQuery tests. Each test runs a query against a randomly selected
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
      TestDataGenerator.ShapeFactory factory, int numIndexes, int numEdges, int numQueries) {
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
      // Exercise the opposite-hemisphere code 1/5 of the time.
      S2Point center = data.oneIn(5) ? indexCap.axis().neg() : indexCap.axis();
      S2Cap queryCap = S2Cap.fromAxisAngle(center, queryRadius);
      S2FurthestEdgeQuery.Builder builder = S2FurthestEdgeQuery.builder();

      // Occasionally we don't set any limit on the number of result edges. (This may return all
      // edges if we also don't set a distance limit.)
      if (!data.oneIn(5)) {
        builder.setMaxResults(1 + data.uniform(10));
      }
      // We set a distance limit 2/3 of the time.
      if (!data.oneIn(3)) {
        builder.setMinDistance(queryRadius.mul(data.nextDouble()));
      }
      if (data.oneIn(2)) {
        // Choose a maximum error whose logarithm is uniformly distributed over a reasonable range,
        // except that it is sometimes zero.
        builder.setMaxError(
            S1ChordAngle.fromRadians(pow(1e-4, data.nextDouble()) * queryRadius.radians()));
      }
      builder.setIncludeInteriors(data.oneIn(2));

      int targetType = data.uniform(4);
      if (targetType == 0) {
        // Find the edges furthest from a given point.
        S2Point point = data.samplePoint(queryCap);
        S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
            new S2FurthestEdgeQuery.PointTarget<>(point);
        testFindFurthestEdges(index, target, builder);
      } else if (targetType == 1) {
        // Find the edges furthest from a given edge.
        S2Point a = data.samplePoint(queryCap);
        S2Point b =
            data.samplePoint(S2Cap.fromAxisAngle(a, queryRadius.mul(pow(1e-4, data.nextDouble()))));
        S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
            new S2FurthestEdgeQuery.EdgeTarget<>(a, b);
        testFindFurthestEdges(index, target, builder);
      } else if (targetType == 2) {
        // Find the edges furthest from a given cell.
        int minLevel = MAX_DIAG.getMinLevel(queryRadius.radians());
        int level = minLevel + data.uniform(S2CellId.MAX_LEVEL - minLevel + 1);
        S2Point a = data.samplePoint(queryCap);
        S2CellId targetId = S2CellId.fromPoint(a).parent(level);
        S2Cell cell = new S2Cell(targetId);
        S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target =
            new S2FurthestEdgeQuery.CellTarget<>(cell);
        testFindFurthestEdges(index, target, builder);
      } else {
        assertEquals(3, targetType);
        // Use another one of the pre-built indexes as the target.
        int jIndex = data.uniform(numIndexes);
        S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
            S2FurthestEdgeQuery.createShapeIndexTarget(indexes.get(jIndex));
        boolean targetIncludeInteriors = data.oneIn(2);
        target.setIncludeInteriors(targetIncludeInteriors);
        testFindFurthestEdges(index, target, builder);
      }
    }
  }

  static int kNumIndexes = 50;
  static int kNumEdges = 100;
  static int kNumQueries = 200;

  @Test
  public void testCircleEdges() {
    testWithShapeFactory(
        data.new RegularLoopShapeFactory(), kNumIndexes, kNumEdges, kNumQueries);
  }

  @Test
  public void testFractalEdges() {
    testWithShapeFactory(
        data.new FractalLoopShapeFactory(), kNumIndexes, kNumEdges, kNumQueries);
  }

  @Test
  public void testPointCloudEdges() {
    testWithShapeFactory(
        data.new PointCloudShapeFactory(), kNumIndexes, kNumEdges, kNumQueries);
  }
}
