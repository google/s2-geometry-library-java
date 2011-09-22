/*
 * Copyright 2005 Google Inc.
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

public strictfp class S2LatLngRectTest extends GeometryTestCase {

  public void testIntervalOps(S2LatLngRect x, S2LatLngRect y, String expectedRelation,
      S2LatLngRect expectedUnion, S2LatLngRect expectedIntersection) {
    // Test all of the interval operations on the given pair of intervals.
    // "expected_relation" is a sequence of "T" and "F" characters corresponding
    // to the expected results of Contains(), InteriorContains(), Intersects(),
    // and InteriorIntersects() respectively.

    assertEquals(x.contains(y), expectedRelation.charAt(0) == 'T');
    assertEquals(x.interiorContains(y), expectedRelation.charAt(1) == 'T');
    assertEquals(x.intersects(y), expectedRelation.charAt(2) == 'T');
    assertEquals(x.interiorIntersects(y), expectedRelation.charAt(3) == 'T');

    assertEquals(x.contains(y), x.union(y).equals(x));
    assertEquals(x.intersects(y), !x.intersection(y).isEmpty());

    assertTrue(x.union(y).equals(expectedUnion));
    assertTrue(x.intersection(y).equals(expectedIntersection));

    if (y.getSize() == S2LatLng.fromRadians(0, 0)) {
      S2LatLngRect r = x.addPoint(y.lo());
      assertTrue(r == expectedUnion);
    }
  }

  public void testCellOps(S2LatLngRect r, S2Cell cell, int level) {
    // Test the relationship between the given rectangle and cell:
    // 0 == no intersection, 1 == MayIntersect, 2 == Intersects,
    // 3 == Vertex Containment, 4 == Contains

    boolean vertexContained = false;
    for (int i = 0; i < 4; ++i) {
      if (r.contains(cell.getVertexRaw(i))
          || (!r.isEmpty() && cell.contains(r.getVertex(i).toPoint()))) {
        vertexContained = true;
      }
    }
    assertEquals(r.mayIntersect(cell), level >= 1);
    assertEquals(r.intersects(cell), level >= 2);
    assertEquals(vertexContained, level >= 3);
    assertEquals(r.contains(cell), level >= 4);
  }

  public void testBasic() {
    // Most of the S2LatLngRect methods have trivial implementations that
    // use the R1Interval and S1Interval classes, so most of the testing
    // is done in those unit tests.

    // Test basic properties of empty and full caps.
    S2LatLngRect empty = S2LatLngRect.empty();
    S2LatLngRect full = S2LatLngRect.full();
    assertTrue(empty.isValid());
    assertTrue(empty.isEmpty());
    assertTrue(full.isValid());
    assertTrue(full.isFull());

    // assertTrue various constructors and accessor methods.
    S2LatLngRect d1 = rectFromDegrees(-90, 0, -45, 180);
    assertDoubleNear(d1.latLo().degrees(), -90);
    assertDoubleNear(d1.latHi().degrees(), -45);
    assertDoubleNear(d1.lngLo().degrees(), 0);
    assertDoubleNear(d1.lngHi().degrees(), 180);
    assertTrue(d1.lat().equals(new R1Interval(-S2.M_PI_2, -S2.M_PI_4)));
    assertTrue(d1.lng().equals(new S1Interval(0, S2.M_PI)));

    // FromCenterSize()
    assertTrue(
        S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(80, 170), S2LatLng.fromDegrees(40, 60))
            .approxEquals(rectFromDegrees(60, 140, 90, -160)));
    assertTrue(S2LatLngRect
        .fromCenterSize(S2LatLng.fromDegrees(10, 40), S2LatLng.fromDegrees(210, 400)).isFull());
    assertTrue(
        S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(-90, 180), S2LatLng.fromDegrees(20, 50))
            .approxEquals(rectFromDegrees(-90, 155, -80, -155)));

    // FromPoint(), FromPointPair()
    assertEquals(S2LatLngRect.fromPoint(d1.lo()), new S2LatLngRect(d1.lo(), d1.lo()));
    assertEquals(
        S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(-35, -140), S2LatLng.fromDegrees(15, 155)),
        rectFromDegrees(-35, 155, 15, -140));
    assertEquals(
        S2LatLngRect.fromPointPair(S2LatLng.fromDegrees(25, -70), S2LatLng.fromDegrees(-90, 80)),
        rectFromDegrees(-90, -70, 25, 80));

    // GetCenter(), GetVertex(), Contains(S2LatLng), InteriorContains(S2LatLng).
    S2LatLng eqM180 = S2LatLng.fromRadians(0, -S2.M_PI);
    S2LatLng northPole = S2LatLng.fromRadians(S2.M_PI_2, 0);
    S2LatLngRect r1 = new S2LatLngRect(eqM180, northPole);

    assertEquals(r1.getCenter(), S2LatLng.fromRadians(S2.M_PI_4, -S2.M_PI_2));
    assertEquals(r1.getVertex(0), S2LatLng.fromRadians(0, S2.M_PI));
    assertEquals(r1.getVertex(1), S2LatLng.fromRadians(0, 0));
    assertEquals(r1.getVertex(2), S2LatLng.fromRadians(S2.M_PI_2, 0));
    assertEquals(r1.getVertex(3), S2LatLng.fromRadians(S2.M_PI_2, S2.M_PI));
    assertTrue(r1.contains(S2LatLng.fromDegrees(30, -45)));
    assertTrue(!r1.contains(S2LatLng.fromDegrees(30, 45)));
    assertTrue(!r1.interiorContains(eqM180) && !r1.interiorContains(northPole));
    assertTrue(r1.contains(new S2Point(0.5, -0.3, 0.1)));
    assertTrue(!r1.contains(new S2Point(0.5, 0.2, 0.1)));

    // Make sure that GetVertex() returns vertices in CCW order.
    for (int i = 0; i < 4; ++i) {
      double lat = S2.M_PI_4 * (i - 2);
      double lng = S2.M_PI_2 * (i - 2) + 0.2;
      S2LatLngRect r = new S2LatLngRect(new R1Interval(lat, lat + S2.M_PI_4), new S1Interval(
          Math.IEEEremainder(lng, 2 * S2.M_PI), Math.IEEEremainder(lng + S2.M_PI_2, 2 * S2.M_PI)));
      for (int k = 0; k < 4; ++k) {
        assertTrue(
            S2.simpleCCW(r.getVertex((k - 1) & 3).toPoint(), r.getVertex(k).toPoint(),
                r.getVertex((k + 1) & 3).toPoint()));
      }
    }

    // Contains(S2LatLngRect), InteriorContains(S2LatLngRect),
    // Intersects(), InteriorIntersects(), Union(), Intersection().
    //
    // Much more testing of these methods is done in s1interval_unittest
    // and r1interval_unittest.

    S2LatLngRect r1Mid = rectFromDegrees(45, -90, 45, -90);
    S2LatLngRect reqM180 = new S2LatLngRect(eqM180, eqM180);
    S2LatLngRect rNorthPole = new S2LatLngRect(northPole, northPole);

    testIntervalOps(r1, r1Mid, "TTTT", r1, r1Mid);
    testIntervalOps(r1, reqM180, "TFTF", r1, reqM180);
    testIntervalOps(r1, rNorthPole, "TFTF", r1, rNorthPole);

    assertTrue(r1.equals(rectFromDegrees(0, -180, 90, 0)));
    testIntervalOps(r1, rectFromDegrees(-10, -1, 1, 20), "FFTT", rectFromDegrees(-10, -180, 90, 20),
        rectFromDegrees(0, -1, 1, 0));
    testIntervalOps(r1, rectFromDegrees(-10, -1, 0, 20), "FFTF", rectFromDegrees(-10, -180, 90, 20),
        rectFromDegrees(0, -1, 0, 0));
    testIntervalOps(r1, rectFromDegrees(-10, 0, 1, 20), "FFTF", rectFromDegrees(-10, -180, 90, 20),
        rectFromDegrees(0, 0, 1, 0));

    testIntervalOps(rectFromDegrees(-15, -160, -15, -150), rectFromDegrees(20, 145, 25, 155),
        "FFFF", rectFromDegrees(-15, 145, 25, -150), empty);
    testIntervalOps(rectFromDegrees(70, -10, 90, -140), rectFromDegrees(60, 175, 80, 5), "FFTT",
        rectFromDegrees(60, -180, 90, 180), rectFromDegrees(70, 175, 80, 5));

    // assertTrue that the intersection of two rectangles that overlap in
    // latitude
    // but not longitude is valid, and vice versa.
    testIntervalOps(rectFromDegrees(12, 30, 60, 60), rectFromDegrees(0, 0, 30, 18), "FFFF",
        rectFromDegrees(0, 0, 60, 60), empty);
    testIntervalOps(rectFromDegrees(0, 0, 18, 42), rectFromDegrees(30, 12, 42, 60), "FFFF",
        rectFromDegrees(0, 0, 42, 60), empty);

    // AddPoint()
    S2LatLngRect p = S2LatLngRect.empty();
    p = p.addPoint(S2LatLng.fromDegrees(0, 0));
    p = p.addPoint(S2LatLng.fromRadians(0, -S2.M_PI_2));
    p = p.addPoint(S2LatLng.fromRadians(S2.M_PI_4, -S2.M_PI));
    p = p.addPoint(new S2Point(0, 0, 1));
    assertTrue(p.equals(r1));

    // Expanded()
    assertTrue(
        rectFromDegrees(70, 150, 80, 170).expanded(S2LatLng.fromDegrees(20, 30)).approxEquals(
            rectFromDegrees(50, 120, 90, -160)));
    assertTrue(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(20, 30)).isEmpty());
    assertTrue(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(20, 30)).isFull());
    assertTrue(
        rectFromDegrees(-90, 170, 10, 20).expanded(S2LatLng.fromDegrees(30, 80)).approxEquals(
            rectFromDegrees(-90, -180, 40, 180)));

    // ConvolveWithCap()
    S2LatLngRect llr1 =
        new S2LatLngRect(S2LatLng.fromDegrees(0, 170), S2LatLng.fromDegrees(0, -170))
            .convolveWithCap(S1Angle.degrees(15));
    S2LatLngRect llr2 =
        new S2LatLngRect(S2LatLng.fromDegrees(-15, 155), S2LatLng.fromDegrees(15, -155));
    assertTrue(llr1.approxEquals(llr2));

    llr1 = new S2LatLngRect(S2LatLng.fromDegrees(60, 150), S2LatLng.fromDegrees(80, 10))
        .convolveWithCap(S1Angle.degrees(15));
    llr2 = new S2LatLngRect(S2LatLng.fromDegrees(45, -180), S2LatLng.fromDegrees(90, 180));
    assertTrue(llr1.approxEquals(llr2));

    // GetCapBound(), bounding cap at center is smaller:
    assertTrue(new S2LatLngRect(S2LatLng.fromDegrees(-45, -45), S2LatLng.fromDegrees(45, 45))
        .getCapBound().approxEquals(S2Cap.fromAxisHeight(new S2Point(1, 0, 0), 0.5)));
    // GetCapBound(), bounding cap at north pole is smaller:
    assertTrue(new S2LatLngRect(S2LatLng.fromDegrees(88, -80), S2LatLng.fromDegrees(89, 80))
        .getCapBound().approxEquals(S2Cap.fromAxisAngle(new S2Point(0, 0, 1), S1Angle.degrees(2))));
    // GetCapBound(), longitude span > 180 degrees:
    assertTrue(
        new S2LatLngRect(S2LatLng.fromDegrees(-30, -150), S2LatLng.fromDegrees(-10, 50))
            .getCapBound()
            .approxEquals(S2Cap.fromAxisAngle(new S2Point(0, 0, -1), S1Angle.degrees(80))));

    // Contains(S2Cell), MayIntersect(S2Cell), Intersects(S2Cell)

    // Special cases.
    testCellOps(empty, S2Cell.fromFacePosLevel(3, (byte) 0, 0), 0);
    testCellOps(full, S2Cell.fromFacePosLevel(2, (byte) 0, 0), 4);
    testCellOps(full, S2Cell.fromFacePosLevel(5, (byte) 0, 25), 4);

    // This rectangle includes the first quadrant of face 0. It's expanded
    // slightly because cell bounding rectangles are slightly conservative.
    S2LatLngRect r4 = rectFromDegrees(-45.1, -45.1, 0.1, 0.1);
    testCellOps(r4, S2Cell.fromFacePosLevel(0, (byte) 0, 0), 3);
    testCellOps(r4, S2Cell.fromFacePosLevel(0, (byte) 0, 1), 4);
    testCellOps(r4, S2Cell.fromFacePosLevel(1, (byte) 0, 1), 0);

    // This rectangle intersects the first quadrant of face 0.
    S2LatLngRect r5 = rectFromDegrees(-10, -45, 10, 0);
    testCellOps(r5, S2Cell.fromFacePosLevel(0, (byte) 0, 0), 3);
    testCellOps(r5, S2Cell.fromFacePosLevel(0, (byte) 0, 1), 3);
    testCellOps(r5, S2Cell.fromFacePosLevel(1, (byte) 0, 1), 0);

    // Rectangle consisting of a single point.
    testCellOps(rectFromDegrees(4, 4, 4, 4), S2Cell.fromFacePosLevel(0, (byte) 0, 0), 3);

    // Rectangles that intersect the bounding rectangle of a face
    // but not the face itself.
    testCellOps(rectFromDegrees(41, -87, 42, -79), S2Cell.fromFacePosLevel(2, (byte) 0, 0), 1);
    testCellOps(rectFromDegrees(-41, 160, -40, -160), S2Cell.fromFacePosLevel(5, (byte) 0, 0), 1);
    {
      // This is the leaf cell at the top right hand corner of face 0.
      // It has two angles of 60 degrees and two of 120 degrees.
      S2Cell cell0tr = new S2Cell(new S2Point(1 + 1e-12, 1, 1));
      S2LatLngRect bound0tr = cell0tr.getRectBound();
      S2LatLng v0 = new S2LatLng(cell0tr.getVertexRaw(0));
      testCellOps(
          rectFromDegrees(v0.lat().degrees() - 1e-8, v0.lng().degrees() - 1e-8,
              v0.lat().degrees() - 2e-10, v0.lng().degrees() + 1e-10), cell0tr, 1);
    }

    // Rectangles that intersect a face but where no vertex of one region
    // is contained by the other region. The first one passes through
    // a corner of one of the face cells.
    testCellOps(rectFromDegrees(-37, -70, -36, -20), S2Cell.fromFacePosLevel(5, (byte) 0, 0), 2);
    {
      // These two intersect like a diamond and a square.
      S2Cell cell202 = S2Cell.fromFacePosLevel(2, (byte) 0, 2);
      S2LatLngRect bound202 = cell202.getRectBound();
      testCellOps(
          rectFromDegrees(bound202.lo().lat().degrees() + 3, bound202.lo().lng().degrees() + 3,
              bound202.hi().lat().degrees() - 3, bound202.hi().lng().degrees() - 3), cell202, 2);
    }
  }

  public void testArea() {
    assertEquals(0.0, S2LatLngRect.empty().area());
    assertDoubleNear(4 * Math.PI, S2LatLngRect.full().area());
    assertDoubleNear(Math.PI / 2, rectFromDegrees(0, 0, 90, 90).area());
  }

  public void testEdgeBound() {
    // assertTrue cases where min/max latitude is not at a vertex.
    assertDoubleNear(getEdgeBound(1, 1, 1, 1, -1, 1).lat().hi(), S2.M_PI_4); // Max,
    // CW
    assertDoubleNear(getEdgeBound(1, -1, 1, 1, 1, 1).lat().hi(), S2.M_PI_4); // Max,
    // CCW
    assertDoubleNear(getEdgeBound(1, -1, -1, -1, -1, -1).lat().lo(), -S2.M_PI_4); // Min,
                                                                                  // CW
    assertDoubleNear(getEdgeBound(-1, 1, -1, -1, -1, -1).lat().lo(), -S2.M_PI_4); // Min,
                                                                                  // CCW

    // assertTrue cases where the edge passes through one of the poles.
    assertDoubleNear(getEdgeBound(.3, .4, 1, -.3, -.4, 1).lat().hi(), S2.M_PI_2);
    assertDoubleNear(getEdgeBound(.3, .4, -1, -.3, -.4, -1).lat().lo(), -S2.M_PI_2);

    // assertTrue cases where the min/max latitude is attained at a vertex.
    final double kCubeLat = Math.asin(Math.sqrt(1. / 3)); // 35.26 degrees
    assertTrue(
        getEdgeBound(1, 1, 1, 1, -1, -1).lat().approxEquals(new R1Interval(-kCubeLat, kCubeLat)));
    assertTrue(
        getEdgeBound(1, -1, 1, 1, 1, -1).lat().approxEquals(new R1Interval(-kCubeLat, kCubeLat)));
  }

  public void testGetDistanceOverlapping() {
    // Check pairs of rectangles that overlap: (should all return 0):
    S2LatLngRect a = rectFromDegrees(0, 0, 2, 2);
    S2LatLngRect b = pointRectFromDegrees(0, 0);
    S1Angle zero = S1Angle.radians(0);
    assertEquals(zero, a.getDistance(a));
    assertEquals(zero, a.getDistance(b));
    assertEquals(zero, b.getDistance(b));
    assertEquals(zero, a.getDistance(S2LatLng.fromDegrees(0, 0)));
    assertEquals(zero, a.getDistance(rectFromDegrees(0, 1, 2, 3)));
    assertEquals(zero, a.getDistance(rectFromDegrees(0, 2, 2, 4)));
    assertEquals(zero, a.getDistance(rectFromDegrees(1, 0, 3, 2)));
    assertEquals(zero, a.getDistance(rectFromDegrees(2, 0, 4, 2)));
    assertEquals(zero, a.getDistance(rectFromDegrees(1, 1, 3, 3)));
    assertEquals(zero, a.getDistance(rectFromDegrees(2, 2, 4, 4)));
  }

  public void testGetDistanceRectVsPoint() {
    // Rect that spans 180.
    S2LatLngRect a = rectFromDegrees(-1, -1, 2, 1);
    verifyGetDistance(a, pointRectFromDegrees(-2, -1));
    verifyGetDistance(a, pointRectFromDegrees(1, 2));

    verifyGetDistance(pointRectFromDegrees(-2, -1), a);
    verifyGetDistance(pointRectFromDegrees(1, 2), a);

    verifyGetRectPointDistance(a, S2LatLng.fromDegrees(-2, -1));
    verifyGetRectPointDistance(a, S2LatLng.fromDegrees(1, 2));

    // Tests near the north pole.
    S2LatLngRect b = rectFromDegrees(86, 0, 88, 2);
    verifyGetDistance(b, pointRectFromDegrees(87, 3));
    verifyGetDistance(b, pointRectFromDegrees(87, -1));
    verifyGetDistance(b, pointRectFromDegrees(89, 1));
    verifyGetDistance(b, pointRectFromDegrees(89, 181));
    verifyGetDistance(b, pointRectFromDegrees(85, 1));
    verifyGetDistance(b, pointRectFromDegrees(85, 181));
    verifyGetDistance(b, pointRectFromDegrees(90, 0));

    verifyGetDistance(pointRectFromDegrees(87, 3), b);
    verifyGetDistance(pointRectFromDegrees(87, -1), b);
    verifyGetDistance(pointRectFromDegrees(89, 1), b);
    verifyGetDistance(pointRectFromDegrees(89, 181), b);
    verifyGetDistance(pointRectFromDegrees(85, 1), b);
    verifyGetDistance(pointRectFromDegrees(85, 181), b);
    verifyGetDistance(pointRectFromDegrees(90, 0), b);

    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(87, 3));
    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(87, -1));
    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(89, 1));
    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(89, 181));
    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(85, 1));
    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(85, 181));
    verifyGetRectPointDistance(b, S2LatLng.fromDegrees(90, 0));

    // Rect that touches the north pole.
    S2LatLngRect c = rectFromDegrees(88, 0, 90, 2);
    verifyGetDistance(c, pointRectFromDegrees(89, 3));
    verifyGetDistance(c, pointRectFromDegrees(89, 90));
    verifyGetDistance(c, pointRectFromDegrees(89, 181));
    verifyGetDistance(pointRectFromDegrees(89, 3), c);
    verifyGetDistance(pointRectFromDegrees(89, 90), c);
    verifyGetDistance(pointRectFromDegrees(89, 181), c);
  }

  public void testGetDistanceRectVsRect() {
    // Rect that spans 180.
    S2LatLngRect a = rectFromDegrees(-1, -1, 2, 1);
    verifyGetDistance(a, rectFromDegrees(0, 2, 1, 3));
    verifyGetDistance(a, rectFromDegrees(-2, -3, -1, -2));

    // Tests near the south pole.
    S2LatLngRect b = rectFromDegrees(-87, 0, -85, 3);
    verifyGetDistance(b, rectFromDegrees(-89, 1, -88, 2));
    verifyGetDistance(b, rectFromDegrees(-84, 1, -83, 2));
    verifyGetDistance(b, rectFromDegrees(-88, 90, -86, 91));
    verifyGetDistance(b, rectFromDegrees(-84, -91, -83, -90));
    verifyGetDistance(b, rectFromDegrees(-90, 181, -89, 182));
    verifyGetDistance(b, rectFromDegrees(-84, 181, -83, 182));
  }

  public void testGetDistanceRandomPairs() {
    // Test random pairs.
    for (int i = 0; i < 10000; ++i) {
      S2LatLngRect a =
          S2LatLngRect.fromPointPair(new S2LatLng(randomPoint()), new S2LatLng(randomPoint()));
      S2LatLngRect b =
          S2LatLngRect.fromPointPair(new S2LatLng(randomPoint()), new S2LatLng(randomPoint()));
      verifyGetDistance(a, b);


      S2LatLng c = new S2LatLng(randomPoint());
      verifyGetRectPointDistance(a, c);
      verifyGetRectPointDistance(b, c);
    }
  }

  private static S1Angle bruteForceDistance(S2LatLngRect a, S2LatLngRect b) {
    if (a.intersects(b)) {
      return S1Angle.radians(0);
    }

    // Compare every point in 'a' against every latitude edge and longitude edge
    // in 'b', and vice-versa, for a total of 16 point-vs-latitude-edge tests
    // and 16 point-vs-longitude-edge tests.
    S2LatLng pntA[] =
        {new S2LatLng(a.latLo(), a.lngLo()), new S2LatLng(a.latLo(), a.lngHi()),
            new S2LatLng(a.latHi(), a.lngHi()), new S2LatLng(a.latHi(), a.lngLo())};
    S2LatLng pntB[] =
        {new S2LatLng(b.latLo(), b.lngLo()), new S2LatLng(b.latLo(), b.lngHi()),
            new S2LatLng(b.latHi(), b.lngHi()), new S2LatLng(b.latHi(), b.lngLo())};

    // Make arrays containing the lo/hi latitudes and the lo/hi longitude edges.
    S1Angle latA[] = {a.latLo(), a.latHi()};
    S1Angle latB[] = {b.latLo(), b.latHi()};
    S2Point lng_edge_a[][] =
        { {pntA[0].toPoint(), pntA[3].toPoint()}, {pntA[1].toPoint(), pntA[2].toPoint()}};
    S2Point lng_edge_b[][] =
        { {pntB[0].toPoint(), pntB[3].toPoint()}, {pntB[1].toPoint(), pntB[2].toPoint()}};

    S1Angle minDistance = S1Angle.degrees(180.0);
    for (int i = 0; i < 4; ++i) {
      // For each point in a and b.
      S2LatLng currentA = pntA[i];
      S2LatLng currentB = pntB[i];

      for (int j = 0; j < 2; ++j) {
        // Get distances to latitude and longitude edges.
        S1Angle aToLat = getDistance(currentA, latB[j], b.lng());
        S1Angle bToLat = getDistance(currentB, latA[j], a.lng());
        S1Angle aToLng =
            S2EdgeUtil.getDistance(currentA.toPoint(), lng_edge_b[j][0], lng_edge_b[j][1]);
        S1Angle bToLng =
            S2EdgeUtil.getDistance(currentB.toPoint(), lng_edge_a[j][0], lng_edge_a[j][1]);

        minDistance = S1Angle.min(
            minDistance, S1Angle.min(aToLat, S1Angle.min(bToLat, S1Angle.min(aToLng, bToLng))));
      }
    }
    return minDistance;
  }

  private static S1Angle bruteForceRectPointDistance(S2LatLngRect a, S2LatLng b) {
    if (a.contains(b)) {
      return S1Angle.radians(0);
    }

    S1Angle bToLoLat = getDistance(b, a.latLo(), a.lng());
    S1Angle bToHiLat = getDistance(b, a.latHi(), a.lng());
    S1Angle bToLoLng =
        S2EdgeUtil.getDistance(b.toPoint(), new S2LatLng(a.latLo(), a.lngLo()).toPoint(),
            new S2LatLng(a.latHi(), a.lngLo()).toPoint());
    S1Angle bToHiLng =
        S2EdgeUtil.getDistance(b.toPoint(), new S2LatLng(a.latLo(), a.lngHi()).toPoint(),
            new S2LatLng(a.latHi(), a.lngHi()).toPoint());
    return S1Angle.min(bToLoLat, S1Angle.min(bToHiLat, S1Angle.min(bToLoLng, bToHiLng)));
  }

  /**
   * Returns the minimum distance from X to the latitude line segment defined by
   * the given latitude and longitude interval.
   */
  private static S1Angle getDistance(S2LatLng x, S1Angle lat, S1Interval interval) {
    assertTrue(x.isValid());
    assertTrue(interval.isValid());

    // Is X inside the longitude interval?
    if (interval.contains(x.lng().radians()))
      return S1Angle.radians(Math.abs(x.lat().radians() - lat.radians()));

    // Return the distance to the closer endpoint.
    return S1Angle.min(x.getDistance(new S2LatLng(lat, S1Angle.radians(interval.lo()))),
        x.getDistance(new S2LatLng(lat, S1Angle.radians(interval.hi()))));
  }

  private static S2LatLngRect getEdgeBound(double x1,
      double y1,
      double z1,
      double x2,
      double y2,
      double z2) {
    return S2LatLngRect.fromEdge(
        S2Point.normalize(new S2Point(x1, y1, z1)), S2Point.normalize(new S2Point(x2, y2, z2)));
  }

  private static S2LatLngRect pointRectFromDegrees(double lat, double lng) {
    return S2LatLngRect.fromPoint(S2LatLng.fromDegrees(lat, lng).normalized());
  }

  private static S2LatLngRect rectFromDegrees(
      double latLo, double lngLo, double latHi, double lngHi) {
    // Convenience method to construct a rectangle. This method is
    // intentionally *not* in the S2LatLngRect interface because the
    // argument order is ambiguous, but hopefully it's not too confusing
    // within the context of this unit test.

    return new S2LatLngRect(S2LatLng.fromDegrees(latLo, lngLo).normalized(),
        S2LatLng.fromDegrees(latHi, lngHi).normalized());
  }

  /**
   * This method verifies a.getDistance(b), where b is a S2LatLng, by comparing
   * its result against a.getDistance(c), c being the point rectangle created
   * from b.
   */
  private static void verifyGetRectPointDistance(S2LatLngRect a, S2LatLng p) {
    S1Angle distance1 = bruteForceRectPointDistance(a, p.normalized());
    S1Angle distance2 = a.getDistance(p.normalized());
    assertEquals(distance1.radians(), distance2.radians(), 1e-10);
  }

  /**
   * This method verifies a.getDistance(b) by comparing its result against a
   * brute-force implementation. The correctness of the brute-force version is
   * much easier to verify by inspection.
   */
  private static void verifyGetDistance(S2LatLngRect a, S2LatLngRect b) {
    S1Angle distance1 = bruteForceDistance(a, b);
    S1Angle distance2 = a.getDistance(b);
    assertEquals(distance1.radians(), distance2.radians(), 1e-10);
  }
}
