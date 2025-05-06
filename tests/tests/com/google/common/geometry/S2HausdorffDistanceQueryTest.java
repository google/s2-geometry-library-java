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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2HausdorffDistanceQuery.DirectedResult;
import com.google.common.geometry.S2HausdorffDistanceQuery.Options;
import com.google.common.geometry.S2HausdorffDistanceQuery.Result;
import java.util.List;
import java.util.Optional;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for S2HausdorffDistanceQuery. */
@RunWith(JUnit4.class)
public final class S2HausdorffDistanceQueryTest extends GeometryTestCase {

  /** Test the constructors and accessors of Result and DirectedResult. */
  @Test
  public void testResultConstructorsAndAccessorsWork() {
    S2Point point1 = S2LatLng.fromDegrees(3, 4).toPoint();
    S2Point point2 = S2LatLng.fromDegrees(5, 6).toPoint();
    S1ChordAngle distance1 = S1ChordAngle.fromDegrees(5);
    S1ChordAngle distance2 = S1ChordAngle.fromDegrees(5);

    DirectedResult directedResult1 = new DirectedResult(distance1, point1);
    DirectedResult directedResult2 = new DirectedResult(distance2, point2);
    Result result12 = new Result(directedResult1, directedResult2);

    assertEquals(directedResult1.targetPoint(), point1);
    assertEquals(directedResult1.distance(), distance1);
    assertEquals(directedResult2.targetPoint(), point2);
    assertEquals(directedResult2.distance(), distance2);

    assertEquals(result12.targetToSource().targetPoint(), point1);
    assertEquals(result12.sourceToTarget().targetPoint(), point2);
    assertEquals(result12.distance(), directedResult2.distance());
  }

  /** Test the constructors and accessors of the Options & Builder. */
  @Test
  public void testOptionsConstructorsAndAccessorsWork() {
    S2HausdorffDistanceQuery.Builder builder = new S2HausdorffDistanceQuery.Builder();
    assertTrue(builder.includeInteriors());

    Options defaultOptions = new Options(builder);
    assertTrue(defaultOptions.includeInteriors());

    builder.setIncludeInteriors(false);
    assertFalse(builder.includeInteriors());

    Options modifiedOptions = new Options(builder);
    assertFalse(modifiedOptions.includeInteriors());
  }

  /** Test the options accessors of the query. */
  @Test
  public void testQueryOptionsAccessorsWorks() {
    S2HausdorffDistanceQuery defaultQuery = S2HausdorffDistanceQuery.builder().build();
    boolean defaultIncludeInteriors = defaultQuery.options().includeInteriors();

    S2HausdorffDistanceQuery modifiedQuery =
        defaultQuery.options().toBuilder().setIncludeInteriors(!defaultIncludeInteriors).build();
    boolean modifiedIncludeInteriors = modifiedQuery.options().includeInteriors();

    assertTrue(defaultIncludeInteriors);
    assertFalse(modifiedIncludeInteriors);
  }

  /** Tests involving single-vertex S2Polylines, which are valid, but edge cases. */
  @Test
  public void testSingleVertexPolylineQueriesSucceed() {
    S2Polyline s0 = new S2Polyline(S2TextFormat.parsePointsOrDie("0.5:0.5"));
    S2Polyline s1 = new S2Polyline(S2TextFormat.parsePointsOrDie("2.5:2.5"));

    // Shape indexes containing the polylines
    S2ShapeIndex s0I = new S2ShapeIndex();
    s0I.add(s0.shape());
    S2ShapeIndex s1I = new S2ShapeIndex();
    s1I.add(s1.shape());

    // Directed distances between single-vertex polylines should just be the distance between those
    // points.
    Optional<DirectedResult> directedS0ToS1 =
        S2HausdorffDistanceQuery.builder().build().getDirectedResult(s0I, s1I);
    S1ChordAngle expected = new S1ChordAngle(s0.vertex(0), s1.vertex(0));
    assertTrue(directedS0ToS1.isPresent());
    assertAlmostEquals(directedS0ToS1.get().distance().degrees(), expected.degrees());

    // Undirected Hausdorff distance between two degenerate edges is the same.
    Optional<Result> undirectedS0ToS1 =
        S2HausdorffDistanceQuery.builder().build().getResult(s0I, s1I);
    assertTrue(undirectedS0ToS1.isPresent());
    assertAlmostEquals(undirectedS0ToS1.get().distance().degrees(), expected.degrees());

    // Check directed Hausdorff distance from an edge to a single-vertex polyline.
    S2Polyline edge = new S2Polyline(S2TextFormat.parsePointsOrDie("1:1, 2:2"));
    // Should be the distance from the further vertex of the edge, which is 2:2.
    expected = new S1ChordAngle(edge.vertex(1), s0.vertex(0));
    S2ShapeIndex edgeI = new S2ShapeIndex();
    edgeI.add(edge.shape());
    Optional<DirectedResult> directedEdgeToSingle =
        S2HausdorffDistanceQuery.builder().build().getDirectedResult(edgeI, s0I);
    assertTrue(directedEdgeToSingle.isPresent());
    assertAlmostEquals(directedEdgeToSingle.get().distance().degrees(), expected.degrees());

    // Check directed Hausdorff distance from a single-vertex polyline to an edge.
    // Should be the distance from the single vertex to the closest vertex of the edge, 1:1
    expected = new S1ChordAngle(s0.vertex(0), edge.vertex(0));
    Optional<DirectedResult> directedSingleToEdge =
        S2HausdorffDistanceQuery.builder().build().getDirectedResult(s0I, edgeI);
    assertTrue(directedSingleToEdge.isPresent());
    assertAlmostEquals(directedSingleToEdge.get().distance().degrees(), expected.degrees());
  }

  /** Test involving 2 simple polyline shape indexes. */
  @Test
  public void testSimplePolylineQueriesSucceed() {
    List<S2Point> a0 = S2TextFormat.parsePointsOrDie("0:0, 0:1, 0:1.5");
    List<S2Point> a1 = S2TextFormat.parsePointsOrDie("0:2, 0:1.5, -10:1");
    List<S2Point> b0 = S2TextFormat.parsePointsOrDie("1:0, 1:1, 3:2");

    // Setup the shape indexes.
    S2ShapeIndex emptyIndex = new S2ShapeIndex();

    // Shape index a consists of 2 polylines, a0 and a1.
    S2ShapeIndex a = new S2ShapeIndex();
    a.add(S2LaxPolylineShape.create(a0));
    a.add(S2LaxPolylineShape.create(a1));

    // Shape index b consists of 1 polyline: b0.
    S2ShapeIndex b = new S2ShapeIndex();
    b.add(S2LaxPolylineShape.create(b0));

    // Calculate expected distances.
    // Directed a to b HD is achieved at the vertex 2 of a1 and vertex 1 of b0.
    S1ChordAngle expectedAToB = new S1ChordAngle(a1.get(2), b0.get(1));
    // Directed b to a HD is achieved at the vertex 2 of b0 and vertex 0 of a1.
    S1ChordAngle expectedBToA = new S1ChordAngle(b0.get(2), a1.get(0));

    S2HausdorffDistanceQuery query = S2HausdorffDistanceQuery.builder().build();

    Optional<DirectedResult> directedEmptyToA = query.getDirectedResult(emptyIndex, a);
    Optional<DirectedResult> directedAToEmpty = query.getDirectedResult(a, emptyIndex);
    S1ChordAngle directedAToEmptyDistance = query.getDirectedDistance(a, emptyIndex);

    // These two should be false since an empty set is an infinite distance away.
    boolean emptyToADirectedDistanceLess =
        query.isDirectedDistanceLess(emptyIndex, a, S1ChordAngle.fromDegrees(360));
    boolean aToEmptyDirectedDistanceLess =
        query.isDirectedDistanceLess(a, emptyIndex, S1ChordAngle.fromDegrees(360));

    assertFalse(directedEmptyToA.isPresent());
    assertFalse(directedAToEmpty.isPresent());
    assertTrue(directedAToEmptyDistance.isInfinity());
    assertFalse(emptyToADirectedDistanceLess);
    assertFalse(aToEmptyDirectedDistanceLess);

    Optional<DirectedResult> directedAToB = query.getDirectedResult(a, b);
    Optional<DirectedResult> directedBToA = query.getDirectedResult(b, a);
    S1ChordAngle directedAToBDistance = query.getDirectedDistance(a, b);

    // Tests for isDirectedDistanceLess with limits near the Hausdorff distance.
    boolean aToBDirectedDistanceLessThanDistancePlus =
        query.isDirectedDistanceLess(
            a, b, S1ChordAngle.add(directedAToBDistance, S1ChordAngle.fromDegrees(1.0)));
    boolean aToBDirectedDistanceLessThanDistanceMinus =
        query.isDirectedDistanceLess(
            a, b, S1ChordAngle.sub(directedAToBDistance, S1ChordAngle.fromDegrees(1.0)));

    assertTrue(directedAToB.isPresent());
    assertTrue(directedBToA.isPresent());
    assertAlmostEquals(directedAToB.get().distance().degrees(), expectedAToB.degrees());
    assertAlmostEquals(directedAToBDistance.degrees(), expectedAToB.degrees());
    assertAlmostEquals(directedBToA.get().distance().degrees(), expectedBToA.degrees());
    assertTrue(aToBDirectedDistanceLessThanDistancePlus);
    assertFalse(aToBDirectedDistanceLessThanDistanceMinus);

    // Tests for undirected cases.
    Optional<Result> aToB = query.getResult(a, b);
    Optional<Result> bToA = query.getResult(b, a);
    S1ChordAngle bToADistance = query.getDistance(b, a);
    Optional<Result> bb = query.getResult(b, b);

    assertTrue(aToB.isPresent());
    assertTrue(bToA.isPresent());
    assertTrue(bb.isPresent());

    // Tests for isDistanceLess with limits near the Hausdorff distance and the average of the two
    // directed Hausdorff distances.
    double largerAAndBDistance =
        max(directedAToB.get().distance().radians(), directedBToA.get().distance().radians());
    double smallerAAndBDistance =
        min(directedAToB.get().distance().radians(), directedBToA.get().distance().radians());
    double averageAAndBDistance = (largerAAndBDistance + smallerAAndBDistance) / 2.0;

    // This should be true if we add a small epsilon upwards to account for any floating point
    // error.
    boolean distanceLessLargerDistance =
        query.isDistanceLess(a, b, S1ChordAngle.fromRadians(largerAAndBDistance + 0.001));
    // The average should cause one direction to succeed and the other to fail so overall this
    // should return false.
    boolean distanceLessAverageDistance =
        query.isDistanceLess(a, b, S1ChordAngle.fromRadians(averageAAndBDistance));
    // The isWithin(Directed)DistanceLimit methods are inclusive so subtract some small epsilon to
    // cause the method to return false.
    boolean distanceLessSmallerDistance =
        query.isDistanceLess(a, b, S1ChordAngle.fromRadians(smallerAAndBDistance - 0.001));

    boolean bbAlwaysWithin = query.isDistanceLess(b, b, S1ChordAngle.fromDegrees(0));

    assertAlmostEquals(aToB.get().distance().degrees(), bToA.get().distance().degrees());
    assertAlmostEquals(bb.get().distance().degrees(), 0);
    assertExactly(
        aToB.get().distance().degrees(),
        max(aToB.get().distance().degrees(), bToA.get().distance().degrees()));
    assertExactly(bToADistance.degrees(), bToA.get().distance().degrees());

    assertTrue(distanceLessLargerDistance);
    assertFalse(distanceLessAverageDistance);
    assertFalse(distanceLessSmallerDistance);
    assertTrue(bbAlwaysWithin);
  }

  /** Test involving a polyline shape (dimension == 1) and a point shape (dimension == 0). */
  @Test
  public void testPointVectorShapeQueriesSucceed() {
    // Points for the polyline shape.
    List<S2Point> aPoints = S2TextFormat.parsePointsOrDie("2:0, 0:1, 1:2, 0:3, 0:4");
    // Points for the point vector shape.
    List<S2Point> bPoints = S2TextFormat.parsePointsOrDie("-1:2, -0.5:0.5, -0.5:3.5");

    S2ShapeIndex a = new S2ShapeIndex();
    a.add(S2LaxPolylineShape.create(aPoints));

    S2ShapeIndex b = new S2ShapeIndex();
    b.add(S2Point.Shape.fromList(bPoints));

    S2HausdorffDistanceQuery query = S2HausdorffDistanceQuery.builder().build();

    // Directed Hausdorff distance from a to b is achieved at the vertex 0 of a and vertex 1 of b.
    S1ChordAngle expectedAToB = new S1ChordAngle(aPoints.get(0), bPoints.get(1));

    // Directed Hausdorff distance from b to a is achieved at the vertex 0 of b and vertex 3 of a.
    S1ChordAngle expectedBToA = new S1ChordAngle(bPoints.get(0), aPoints.get(3));

    // Undirected Hausdorff distance between a and b is the maximum of the two directed Hausdorff
    // distances.
    S1ChordAngle expectedAB = S1ChordAngle.max(expectedAToB, expectedBToA);

    Optional<DirectedResult> directedAToB = query.getDirectedResult(a, b);
    Optional<DirectedResult> directedBToA = query.getDirectedResult(b, a);
    S1ChordAngle undirectedAB = query.getDistance(a, b);

    assertTrue(directedAToB.isPresent());
    assertTrue(directedBToA.isPresent());
    assertFalse(undirectedAB.isInfinity());
    assertExactly(undirectedAB.degrees(), expectedAB.degrees());

    assertAlmostEquals(directedAToB.get().distance().degrees(), expectedAToB.degrees());
    assertEquals(directedAToB.get().targetPoint(), aPoints.get(0));

    assertAlmostEquals(directedBToA.get().distance().degrees(), expectedBToA.degrees());
    assertEquals(directedBToA.get().targetPoint(), bPoints.get(0));

    boolean aToBDirectedDistanceLessPlus =
        query.isDirectedDistanceLess(a, b, S1ChordAngle.fromDegrees(expectedAToB.degrees() + 0.01));
    boolean bToADirectedDistanceLessPlus =
        query.isDirectedDistanceLess(b, a, S1ChordAngle.fromDegrees(expectedBToA.degrees() + 0.01));
    boolean aToBDirectedDistanceLessMinus =
        query.isDirectedDistanceLess(a, b, S1ChordAngle.fromDegrees(expectedAToB.degrees() - 0.01));
    boolean bToADirectedDistanceLessMinus =
        query.isDirectedDistanceLess(b, a, S1ChordAngle.fromDegrees(expectedBToA.degrees() - 0.01));

    boolean aBDistanceLessPlus =
        query.isDistanceLess(a, b, S1ChordAngle.fromDegrees(expectedAB.degrees() + 0.01));
    boolean bADistanceLessMinus =
        query.isDistanceLess(b, a, S1ChordAngle.fromDegrees(expectedBToA.degrees() - 0.01));

    assertTrue(aToBDirectedDistanceLessPlus);
    assertTrue(bToADirectedDistanceLessPlus);
    assertFalse(aToBDirectedDistanceLessMinus);
    assertFalse(bToADirectedDistanceLessMinus);
    assertTrue(aBDistanceLessPlus);
    assertFalse(bADistanceLessMinus);
  }

  /** Test involving partially overlapping polygons. */
  @Test
  public void testOverlappingPolygons() {
    // The first polygon is a triangle. Its first two vertices are inside the quadrangle b (defined
    // below), and the last vertex is outside of b.
    S2ShapeIndex a = new S2ShapeIndex();
    a.add(S2TextFormat.makeLaxPolygonOrDie("1:1, 1:2, 3.5:1.5"));

    // The other polygon is a quadrangle.
    S2ShapeIndex b = new S2ShapeIndex();
    b.add(S2TextFormat.makeLaxPolygonOrDie("0:0, 0:3, 3:3, 3:0"));

    // A triangle.
    S2ShapeIndex c = new S2ShapeIndex();
    c.add(S2TextFormat.makeLaxPolygonOrDie("0:0, 0:2, 3:0"));

    // Error tolerance to account for the difference between the northern edge of the quadrangle,
    // which is a geodesic line, and the parallel lat=3 connecting the vertices of that edge.
    final double kEpsilon = 3.0e-3;

    // The first query does not include the interiors.
    S2HausdorffDistanceQuery query1 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(false).build();

    // The directed Hausdorff distance from the first query is achieved on the vertex of the
    // triangle that is inside the quadrangle, and is approximately 1 degree away from the nearest
    // edge of the quadrangle.
    S2Point expectedTargetPoint1 = S2LatLng.fromDegrees(1, 2).toPoint();

    Optional<DirectedResult> aToB1 = query1.getDirectedResult(a, b);

    boolean cToBLessThan =
        query1.isDirectedDistanceLess(c, b, S1ChordAngle.fromDegrees(1.0 + kEpsilon));

    assertTrue(aToB1.isPresent());
    assertDoubleNear(aToB1.get().distance().degrees(), 1, kEpsilon);
    assertEquals(expectedTargetPoint1, aToB1.get().targetPoint());
    assertTrue(cToBLessThan);

    // The second query has includeInteriors set to true.
    S2HausdorffDistanceQuery query2 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(true).build();
    Optional<DirectedResult> aToB2 = query2.getDirectedResult(a, b);

    // The directed Hausdorff distance from the second query is achieved on the last vertex of the
    // triangle that is outside the quadrangle, and is about 0.5 degrees away from the nearest edge
    // of the quadrangle.
    S2Point expectedTargetPoint2 = S2LatLng.fromDegrees(3.5, 1.5).toPoint();
    // C is fully contained in B so all points are 0 distance to B.
    boolean cToBLessThan2 = query2.isDirectedDistanceLess(c, b, S1ChordAngle.fromDegrees(kEpsilon));

    assertTrue(aToB2.isPresent());
    assertDoubleNear(aToB2.get().distance().degrees(), 0.5, kEpsilon);
    assertEquals(expectedTargetPoint2, aToB2.get().targetPoint());
    assertTrue(cToBLessThan2);
  }

  // Test involving full geometries.
  @Test
  public void testWholeWorld() {
    S2ShapeIndex a = new S2ShapeIndex();
    a.add(S2Point.Shape.fromList(S2TextFormat.parsePointsOrDie("1:1")));

    S2ShapeIndex b = S2TextFormat.makeIndexOrDie("# # full");

    S2HausdorffDistanceQuery query1 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(true).build();
    Optional<DirectedResult> aToB1 = query1.getDirectedResult(a, b);

    assertTrue(aToB1.isPresent());
    assertExactly(aToB1.get().distance().degrees(), 0.0);

    // Going from full geometry to non-full geometry should return empty option.
    Optional<DirectedResult> bToA1 = query1.getDirectedResult(b, a);

    assertFalse(bToA1.isPresent());

    Optional<Result> undirectedAToB = query1.getResult(b, a);
    Optional<Result> undirectedBToA = query1.getResult(a, b);

    assertFalse(undirectedAToB.isPresent());
    assertFalse(undirectedBToA.isPresent());

    // A point to the whole world should always work.
    boolean aToBDirectedDistanceLessZero = query1.isDirectedDistanceLess(a, b, S1ChordAngle.ZERO);
    // The whole world to a point should always fail.
    boolean bToADirectedDistanceLessInf =
        query1.isDirectedDistanceLess(b, a, S1ChordAngle.INFINITY);
    // In undirected case, must consider distance from full geometry to single point which is
    // infinite.
    boolean aBDistanceLessZero = query1.isDistanceLess(a, b, S1ChordAngle.INFINITY);

    assertTrue(aToBDirectedDistanceLessZero);
    assertFalse(bToADirectedDistanceLessInf);
    assertFalse(aBDistanceLessZero);
  }

  @Test
  public void testWholeWorldSameReference() {
    S2ShapeIndex a = S2TextFormat.makeIndexOrDie("# # full");
    S2ShapeIndex b = S2TextFormat.makeIndexOrDie("# # full");

    S2HausdorffDistanceQuery query1 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(true).build();
    Optional<Result> aToB = query1.getResult(a, b);
    assertFalse(aToB.isPresent());

    S2HausdorffDistanceQuery query2 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(true).build();
    Optional<Result> aToA = query1.getResult(a, a);
    assertFalse(aToA.isPresent());

    boolean aToBDistanceLessInf = query2.isDistanceLess(a, b, S1ChordAngle.INFINITY);
    boolean aToADistanceLessInf = query2.isDistanceLess(a, a, S1ChordAngle.INFINITY);

    assertFalse(aToBDistanceLessInf);
    assertFalse(aToADistanceLessInf);
  }
}
