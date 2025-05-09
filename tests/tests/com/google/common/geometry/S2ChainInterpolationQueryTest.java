/*
 * Copyright 2021 Google Inc.
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

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for S2ChainInterpolationQuery. */
@RunWith(JUnit4.class)
public class S2ChainInterpolationQueryTest extends GeometryTestCase {
  static final double EPSILON = 1.e-8;
  static final S1Angle TOLERANCE = S1Angle.degrees(0.01);

  private static class TestResult {
    S2Point point;
    int edgeId;
    S1Angle distance;

    TestResult(S2Point point, int edgeId, S1Angle distance) {
      this.point = point;
      this.edgeId = edgeId;
      this.distance = distance;
    }
  }

  @Test
  public void testS1AngleQuerySimplePolylines() {
    // Set up the test inputs.
    final double latitudeB = 1.0;
    final double latitudeC = 2.5;
    final double totalLengthABC = latitudeC;

    S2Point a = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(latitudeB, 0).toPoint();
    S2Point c = S2LatLng.fromDegrees(latitudeC, 0).toPoint();

    S2LaxPolylineShape emptyShape = S2LaxPolylineShape.EMPTY;
    S2LaxPolylineShape shapeAC = S2LaxPolylineShape.create(ImmutableList.of(a, c));
    S2LaxPolylineShape shapeABC = S2LaxPolylineShape.create(ImmutableList.of(a, b, c));
    S2LaxPolylineShape shapeBB = S2LaxPolylineShape.create(ImmutableList.of(b, b));
    S2Polyline polylineCC = new S2Polyline(ImmutableList.of(c));

    S2ChainInterpolationQuery queryEmpty = new S2ChainInterpolationQuery(emptyShape);
    S2ChainInterpolationQuery queryAC = new S2ChainInterpolationQuery(shapeAC);
    S2ChainInterpolationQuery queryABC = new S2ChainInterpolationQuery(shapeABC);
    S2ChainInterpolationQuery queryBB = new S2ChainInterpolationQuery(shapeBB);
    S2ChainInterpolationQuery queryCC = new S2ChainInterpolationQuery(polylineCC);

    ArrayList<TestResult> groundTruth = new ArrayList<>();
    ImmutableList<Double> distances =
        ImmutableList.of(
            -1.0,
            0.0,
            1.0e-8,
            latitudeB / 2,
            latitudeB - 1.0e-7,
            latitudeB,
            latitudeB + 1.0e-5,
            latitudeB + 0.5,
            latitudeC - 10.e-7, //
            latitudeC,
            latitudeC + 10.e-16,
            1.e6);

    for (double distance : distances) {
      double lat = max(0.0, min(totalLengthABC, distance));
      final S2Point point = S2LatLng.fromDegrees(lat, 0).toPoint();
      int edgeId = distance < latitudeB ? 0 : 1;
      groundTruth.add(new TestResult(point, edgeId, S1Angle.degrees(lat)));
    }

    // Run the tests.
    final double lengthEmpty = queryEmpty.getLength().degrees();
    final double lengthABC = queryABC.getLength().degrees();
    final double lengthAC = queryAC.getLength().degrees();
    final double lengthBB = queryBB.getLength().degrees();
    final double lengthCC = queryCC.getLength().degrees();

    boolean acResultAtInfinity = queryAC.findPoint(S1Angle.INFINITY);
    final S2Point acPointAtInfinity = queryAC.resultPoint();

    ArrayList<TestResult> resultsAC = new ArrayList<>();
    ArrayList<TestResult> resultsABC = new ArrayList<>();
    ArrayList<TestResult> resultsBB = new ArrayList<>();
    ArrayList<TestResult> resultsCC = new ArrayList<>();

    boolean emptyQueryResult = false;
    for (final double distance : distances) {
      final double totalFraction = distance / totalLengthABC;
      emptyQueryResult |= queryEmpty.findPointAtFraction(totalFraction);
      if (queryAC.findPointAtFraction(totalFraction)) {
        resultsAC.add(
            new TestResult(
                queryAC.resultPoint(), queryAC.resultEdgeId(), queryAC.resultDistance()));
      }
      if (queryABC.findPointAtFraction(totalFraction)) {
        resultsABC.add(
            new TestResult(
                queryABC.resultPoint(), queryABC.resultEdgeId(), queryABC.resultDistance()));
      }
      if (queryBB.findPointAtFraction(totalFraction)) {
        resultsBB.add(
            new TestResult(
                queryBB.resultPoint(), queryBB.resultEdgeId(), queryBB.resultDistance()));
      }
      if (queryCC.findPointAtFraction(totalFraction)) {
        resultsCC.add(
            new TestResult(
                queryCC.resultPoint(), queryCC.resultEdgeId(), queryCC.resultDistance()));
      }
    }

    // Check the test results.
    assertFalse(emptyQueryResult);
    assertTrue(acResultAtInfinity);

    assertTrue(lengthEmpty <= EPSILON);
    assertEquals(totalLengthABC, lengthAC, EPSILON);
    assertEquals(totalLengthABC, lengthABC, EPSILON);
    assertTrue(lengthBB <= EPSILON);
    assertTrue(lengthCC <= EPSILON);

    assertEquals(groundTruth.size(), resultsAC.size());
    assertEquals(groundTruth.size(), resultsABC.size());
    assertEquals(groundTruth.size(), resultsBB.size());
    assertEquals(groundTruth.size(), resultsCC.size());

    assertTrue(new S1Angle(acPointAtInfinity, c).degrees() <= EPSILON);

    for (int i = 0; i < groundTruth.size(); ++i) {
      assertTrue(
          new S1Angle(resultsAC.get(i).point, groundTruth.get(i).point).degrees() <= EPSILON);
      assertTrue(
          new S1Angle(resultsABC.get(i).point, groundTruth.get(i).point).degrees() <= EPSILON);
      assertTrue(new S1Angle(resultsBB.get(i).point, shapeBB.vertex(0)).degrees() <= EPSILON);

      assertEquals(0, resultsAC.get(i).edgeId);
      assertEquals(0, resultsBB.get(i).edgeId);
      assertEquals(0, resultsCC.get(i).edgeId);
      assertEquals(groundTruth.get(i).edgeId, resultsABC.get(i).edgeId);
    }
  }

  @Test
  public void testS2ParametricQueryDistance() {
    // Set up the test inputs.
    final ImmutableList<Double> distances =
        ImmutableList.of(
            -1.0,
            -1.0e-8,
            0.0,
            1.0e-8,
            0.2,
            0.5,
            1.0 - 1.0e-8,
            1.0,
            1.0 + 1.e-8,
            1.2,
            1.2,
            1.2 + 1.0e-10,
            1.5,
            1.999999,
            2.0,
            2.00000001,
            1.e6);

    final List<S2Point> vertices =
        S2TextFormat.parsePointsOrDie(
            "0:0, 0:0, 1.0e-7:0, 0.1:0, 0.2:0, 0.2:0, 0.6:0, 0.999999:0, 0.999999:0, "
                + "1:0, 1:0, 1.000001:0, 1.000001:0, 1.1:0, 1.2:0, 1.2000001:0, 1.7:0, "
                + "1.99999999:0, 2:0");

    final double totalLength = new S1Angle(vertices.get(0), Iterables.getLast(vertices)).degrees();

    S2LaxPolylineShape shape = S2LaxPolylineShape.create(vertices);
    S2ChainInterpolationQuery query = new S2ChainInterpolationQuery(shape);

    // Run the tests.
    final double length = query.getLength().degrees();
    ArrayList<TestResult> results = new ArrayList<>();
    for (final double d : distances) {
      boolean findResult = query.findPoint(S1Angle.degrees(d));
      assertTrue(findResult);
      results.add(
          new TestResult(query.resultPoint(), query.resultEdgeId(), query.resultDistance()));
    }

    // Check the test results.
    assertDoubleNear(totalLength, length, EPSILON);
    assertEquals(distances.size(), results.size());

    for (int i = 0; i < distances.size(); ++i) {
      final double queryDistance = distances.get(i);
      final double lat = S2LatLng.latitude(results.get(i).point).degrees();
      final int edgeId = results.get(i).edgeId;
      final S1Angle distance = results.get(i).distance;

      if (queryDistance < 0) {
        assertAlmostEquals(lat, 0);
        assertEquals(0, edgeId);
        assertAlmostEquals(distance.degrees(), 0.0);
      } else if (queryDistance > 2) {
        assertEquals(2, lat, EPSILON);
        assertEquals(shape.numEdges() - 1, edgeId);
        assertAlmostEquals(distance.degrees(), totalLength);
      } else {
        assertEquals(queryDistance, lat, EPSILON);
        assertTrue(edgeId >= 0);
        assertTrue(edgeId <= shape.numEdges());
        S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
        shape.getEdge(edgeId, edge);
        assertTrue(lat >= S2LatLng.latitude(edge.a).degrees());
        assertTrue(lat <= S2LatLng.latitude(edge.b).degrees());
        assertEquals(queryDistance, distance.degrees(), EPSILON);
      }
    }
  }

  @Test
  public void testS2ParametricQueryChains() {
    // Set up the test input. Two closed degenerate loops, each with two edges.
    ArrayList<List<S2Point>> loops = new ArrayList<>();
    loops.add(makeInvalidLoop("0:0, 1:0").vertices());
    loops.add(makeInvalidLoop("2:0, 3:0").vertices());
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(loops);

    S2ChainInterpolationQuery query = new S2ChainInterpolationQuery(shape);
    S2ChainInterpolationQuery query0 = new S2ChainInterpolationQuery(shape, 0);
    S2ChainInterpolationQuery query1 = new S2ChainInterpolationQuery(shape, 1);

    // Run the tests.
    boolean queryResult = query.findPointAtFraction(0.25);
    boolean query0Result = query0.findPointAtFraction(0.25);
    boolean query1Result = query1.findPointAtFraction(0.25);

    // Check the test results.
    assertTrue(queryResult);
    assertTrue(query0Result);
    assertTrue(query1Result);

    assertEquals(1, S2LatLng.latitude(query.resultPoint()).degrees(), EPSILON);
    assertEquals(0.5, S2LatLng.latitude(query0.resultPoint()).degrees(), EPSILON);
    assertEquals(2.5, S2LatLng.latitude(query1.resultPoint()).degrees(), EPSILON);
  }

  @Test
  public void testGetLengthAtEdgeEmpty() {
    S2LaxPolylineShape shape = S2LaxPolylineShape.EMPTY;
    S2ChainInterpolationQuery query = new S2ChainInterpolationQuery(shape);
    assertEquals(query.getLengthAtEdgeEnd(0), S1Angle.ZERO);
  }

  @Test
  public void testGetLengthAtEdgePolyline() {
    S2Point vertex0 = S2LatLng.fromDegrees(0.0, 0.0).toPoint();
    S2Point vertex1 = S2LatLng.fromDegrees(0.0, 1.0).toPoint();
    S2Point vertex2 = S2LatLng.fromDegrees(0.0, 3.0).toPoint();
    S2Point vertex3 = S2LatLng.fromDegrees(0.0, 6.0).toPoint();

    S2LaxPolylineShape shape =
        S2LaxPolylineShape.create(ImmutableList.of(vertex0, vertex1, vertex2, vertex3));

    S2ChainInterpolationQuery query = new S2ChainInterpolationQuery(shape);

    assertEquals(query.getLength(), S1Angle.degrees(6.0));
    assertEquals(query.getLengthAtEdgeEnd(-100), S1Angle.INFINITY);
    assertEquals(query.getLengthAtEdgeEnd(0), S1Angle.degrees(1.0));
    assertEquals(query.getLengthAtEdgeEnd(1), S1Angle.degrees(3.0));
    assertEquals(query.getLengthAtEdgeEnd(2), S1Angle.degrees(6.0));
    assertEquals(query.getLengthAtEdgeEnd(100), S1Angle.INFINITY);
  }

  @Test
  public void testGetLengthAtEdgePolygon() {
    S2Point vertex0 = S2LatLng.fromDegrees(1.0, 1.0).toPoint();
    S2Point vertex1 = S2LatLng.fromDegrees(2.0, 1.0).toPoint();
    S2Point vertex2 = S2LatLng.fromDegrees(2.0, 3.0).toPoint();
    S2Point vertex3 = S2LatLng.fromDegrees(1.0, 3.0).toPoint();
    S2Point vertex4 = S2LatLng.fromDegrees(0.0, 0.0).toPoint();
    S2Point vertex5 = S2LatLng.fromDegrees(0.0, 4.0).toPoint();
    S2Point vertex6 = S2LatLng.fromDegrees(3.0, 4.0).toPoint();
    S2Point vertex7 = S2LatLng.fromDegrees(3.0, 0.0).toPoint();

    List<List<S2Point>> loops = new ArrayList<>();
    loops.add(ImmutableList.of(vertex0, vertex1, vertex2, vertex3));
    loops.add(ImmutableList.of(vertex4, vertex5, vertex6, vertex7));

    S2LaxPolygonShape shape = S2LaxPolygonShape.create(loops);

    S2ChainInterpolationQuery query0 = new S2ChainInterpolationQuery(shape, 0);

    assertDoubleNear(query0.getLength().degrees(), 6.0, TOLERANCE.degrees());
    assertEquals(query0.getLengthAtEdgeEnd(-100), S1Angle.INFINITY);
    assertDoubleNear(query0.getLengthAtEdgeEnd(0).degrees(), 1.0, TOLERANCE.degrees());
    assertDoubleNear(query0.getLengthAtEdgeEnd(1).degrees(), 3.0, TOLERANCE.degrees());
    assertDoubleNear(query0.getLengthAtEdgeEnd(2).degrees(), 4.0, TOLERANCE.degrees());
    assertDoubleNear(query0.getLengthAtEdgeEnd(3).degrees(), 6.0, TOLERANCE.degrees());
    assertEquals(query0.getLengthAtEdgeEnd(4), S1Angle.INFINITY);
    assertEquals(query0.getLengthAtEdgeEnd(5), S1Angle.INFINITY);
    assertEquals(query0.getLengthAtEdgeEnd(6), S1Angle.INFINITY);
    assertEquals(query0.getLengthAtEdgeEnd(7), S1Angle.INFINITY);
    assertEquals(query0.getLengthAtEdgeEnd(100), S1Angle.INFINITY);

    S2ChainInterpolationQuery query1 = new S2ChainInterpolationQuery(shape, 1);

    assertDoubleNear(query1.getLength().degrees(), 14.0, TOLERANCE.degrees());
    assertEquals(query1.getLengthAtEdgeEnd(-100), S1Angle.INFINITY);
    assertEquals(query1.getLengthAtEdgeEnd(0), S1Angle.INFINITY);
    assertEquals(query1.getLengthAtEdgeEnd(1), S1Angle.INFINITY);
    assertEquals(query1.getLengthAtEdgeEnd(2), S1Angle.INFINITY);
    assertEquals(query1.getLengthAtEdgeEnd(3), S1Angle.INFINITY);
    assertDoubleNear(query1.getLengthAtEdgeEnd(4).degrees(), 4.0, TOLERANCE.degrees());
    assertDoubleNear(query1.getLengthAtEdgeEnd(5).degrees(), 7.0, TOLERANCE.degrees());
    assertDoubleNear(query1.getLengthAtEdgeEnd(6).degrees(), 11.0, TOLERANCE.degrees());
    assertDoubleNear(query1.getLengthAtEdgeEnd(7).degrees(), 14.0, TOLERANCE.degrees());
    assertEquals(query1.getLengthAtEdgeEnd(100), S1Angle.INFINITY);
  }

  @Test
  public void testSlice() {
    S2ChainInterpolationQuery emptyQuery = new S2ChainInterpolationQuery(S2LaxPolylineShape.EMPTY);
    assertEquals("", S2TextFormat.toString(emptyQuery.slice(0, 1)));

    S2LaxPolylineShape polyline = S2TextFormat.makeLaxPolylineOrDie("0:0, 0:1, 0:2");
    S2ChainInterpolationQuery query = new S2ChainInterpolationQuery(polyline);
    assertEquals("0:0, 0:1, 0:2", S2TextFormat.toString(query.slice(0, 1)));
    assertEquals("0:0, 0:1", S2TextFormat.toString(query.slice(0, 0.5)));
    assertEquals("0:2, 0:1", S2TextFormat.toString(query.slice(1, 0.5)));
    assertEquals("0:0.5, 0:1, 0:1.5", S2TextFormat.toString(query.slice(0.25, 0.75)));
  }
}
