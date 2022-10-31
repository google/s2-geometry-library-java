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

import com.google.common.geometry.S2HausdorffDistanceQuery.DirectedResult;
import com.google.common.geometry.S2HausdorffDistanceQuery.Options;
import com.google.common.geometry.S2HausdorffDistanceQuery.Result;
import java.util.List;
import java.util.Optional;

/** Tests for S2HausdorffDistanceQuery. */
public final class S2HausdorffDistanceQueryTest extends GeometryTestCase {

  /** Test the constructors and accessors of Result and DirectedResult. */
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
  public void testQueryOptionsAccessorsWorks() {
    S2HausdorffDistanceQuery defaultQuery = S2HausdorffDistanceQuery.builder().build();
    boolean defaultIncludeInteriors = defaultQuery.options().includeInteriors();

    S2HausdorffDistanceQuery modifiedQuery =
        defaultQuery.options().toBuilder().setIncludeInteriors(!defaultIncludeInteriors).build();
    boolean modifiedIncludeInteriors = modifiedQuery.options().includeInteriors();

    assertTrue(defaultIncludeInteriors);
    assertFalse(modifiedIncludeInteriors);
  }

  /** Test involving 2 simple polyline shape indexes. */
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

    // Shape index b consists of 1 polylines: b0.
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

    Optional<DirectedResult> directedAToB = query.getDirectedResult(a, b);
    Optional<DirectedResult> directedBToA = query.getDirectedResult(b, a);
    S1ChordAngle directedAToBDistance = query.getDirectedDistance(a, b);

    Optional<Result> aToB = query.getResult(a, b);
    Optional<Result> bToA = query.getResult(b, a);
    S1ChordAngle bToADistance = query.getDistance(b, a);
    Optional<Result> bb = query.getResult(b, b);

    assertFalse(directedEmptyToA.isPresent());
    assertFalse(directedAToEmpty.isPresent());
    assertTrue(directedAToEmptyDistance.isInfinity());

    assertTrue(directedAToB.isPresent());
    assertTrue(directedBToA.isPresent());

    assertDoubleEquals(directedAToB.get().distance().degrees(), expectedAToB.degrees());
    assertDoubleEquals(directedAToBDistance.degrees(), expectedAToB.degrees());
    assertDoubleEquals(directedBToA.get().distance().degrees(), expectedBToA.degrees());

    assertTrue(aToB.isPresent());
    assertTrue(bToA.isPresent());
    assertTrue(bb.isPresent());

    assertDoubleEquals(aToB.get().distance().degrees(), bToA.get().distance().degrees());
    assertDoubleEquals(bb.get().distance().degrees(), 0);
    assertExactly(
        aToB.get().distance().degrees(),
        max(aToB.get().distance().degrees(), bToA.get().distance().degrees()));
    assertExactly(bToADistance.degrees(), bToA.get().distance().degrees());
  }

  /** Test involving a polyline shape (dimension == 1) and a point shape (dimension == 0). */
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

    assertDoubleEquals(directedAToB.get().distance().degrees(), expectedAToB.degrees());
    assertEquals(directedAToB.get().targetPoint(), aPoints.get(0));

    assertDoubleEquals(directedBToA.get().distance().degrees(), expectedBToA.degrees());
    assertEquals(directedBToA.get().targetPoint(), bPoints.get(0));

    assertExactly(undirectedAB.degrees(), expectedAB.degrees());
  }

  /** Test involving partially overlapping polygons. */
  public void testOverlappingPolygons() {
    // The first polygon is a triangle. Its first two vertices are inside the quadrangle b (defined
    // below), and the last vertex is outside of b.
    S2ShapeIndex a = new S2ShapeIndex();
    a.add(S2TextFormat.makeLaxPolygonOrDie("1:1, 1:2, 3.5:1.5"));

    // The other polygon is a quadrangle.
    S2ShapeIndex b = new S2ShapeIndex();
    b.add(S2TextFormat.makeLaxPolygonOrDie("0:0, 0:3, 3:3, 3:0"));

    // The first query does not include the interiors.
    S2HausdorffDistanceQuery query1 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(false).build();
    Optional<DirectedResult> aToB1 = query1.getDirectedResult(a, b);

    // The second query has includeInteriors set to true.
    S2HausdorffDistanceQuery query2 =
        S2HausdorffDistanceQuery.builder().setIncludeInteriors(true).build();
    Optional<DirectedResult> aToB2 = query2.getDirectedResult(a, b);

    // Error tolerance to account for the difference between the northern edge of the quadrangle,
    // which is a geodesic line, and the parallel lat=3 connecting the vertices of that edge.
    final double kEpsilon = 3.0e-3;

    // The directed Hausdorff distance from the first query is achieved on the vertex of the
    // triangle that is inside the quadrangle, and is approximately 1 degree away from the nearest
    // edge of the quadrangle.
    S2Point expectedTargetPoint1 = S2LatLng.fromDegrees(1, 2).toPoint();

    // The directed Hausdorff distance from the second query is achieved on the last vertex of the
    // triangle that is outside the quadrangle, and is about 0.5 degrees away from the nearest edge
    // of the quadrangle.
    S2Point expectedTargetPoint2 = S2LatLng.fromDegrees(3.5, 1.5).toPoint();

    assertTrue(aToB1.isPresent());
    assertDoubleNear(aToB1.get().distance().degrees(), 1, kEpsilon);
    assertEquals(expectedTargetPoint1, aToB1.get().targetPoint());

    assertTrue(aToB2.isPresent());
    assertDoubleNear(aToB2.get().distance().degrees(), 0.5, kEpsilon);
    assertEquals(expectedTargetPoint2, aToB2.get().targetPoint());
  }
}
