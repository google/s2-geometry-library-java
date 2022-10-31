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

import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2.M_PI_4;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/** Tests for {@link S2LatLngRect}. */
public strictfp class S2LatLngRectTest extends GeometryTestCase {
  public void testIntervalOps(
      S2LatLngRect x,
      S2LatLngRect y,
      String expectedRelation,
      S2LatLngRect expectedUnion,
      S2LatLngRect expectedIntersection) {
    // Test all of the interval operations on the given pair of intervals. "expected_relation" is a
    // sequence of "T" and "F" characters corresponding to the expected results of Contains(),
    // InteriorContains(), Intersects(), and InteriorIntersects() respectively.

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
    // Test the relationship between the given rectangle and cell: 0 == no intersection, 1 ==
    // MayIntersect, 2 == Intersects, 3 == Vertex Containment, 4 == Contains

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
    // Most of the S2LatLngRect methods have trivial implementations that use the R1Interval and
    // S1Interval classes, so most of the testing is done in those unit tests.

    // Test basic properties of empty and full caps.
    S2LatLngRect empty = S2LatLngRect.empty();
    S2LatLngRect full = S2LatLngRect.full();
    assertTrue(empty.isValid());
    assertTrue(empty.isEmpty());
    assertFalse(empty.isPoint());
    assertTrue(full.isValid());
    assertTrue(full.isFull());
    assertFalse(full.isPoint());

    // assertTrue various constructors and accessor methods.
    S2LatLngRect d1 = rectFromDegrees(-90, 0, -45, 180);
    assertDoubleEquals(d1.latLo().degrees(), -90);
    assertDoubleEquals(d1.latHi().degrees(), -45);
    assertDoubleEquals(d1.lngLo().degrees(), 0);
    assertDoubleEquals(d1.lngHi().degrees(), 180);
    assertEquals(new R1Interval(-M_PI_2, -M_PI_4), d1.lat());
    assertEquals(new S1Interval(0, PI), d1.lng());

    // FromCenterSize()
    assertTrue(
        S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(80, 170), S2LatLng.fromDegrees(40, 60))
            .approxEquals(rectFromDegrees(60, 140, 90, -160)));
    assertTrue(
        S2LatLngRect.fromCenterSize(S2LatLng.fromDegrees(10, 40), S2LatLng.fromDegrees(210, 400))
            .isFull());
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
    S2LatLng eqM180 = S2LatLng.fromRadians(0, -PI);
    S2LatLng northPole = S2LatLng.fromRadians(M_PI_2, 0);
    S2LatLngRect r1 = new S2LatLngRect(eqM180, northPole);

    assertEquals(r1.getCenter(), S2LatLng.fromRadians(M_PI_4, -M_PI_2));
    assertEquals(r1.getVertex(0), S2LatLng.fromRadians(0, PI));
    assertEquals(r1.getVertex(1), S2LatLng.fromRadians(0, 0));
    assertEquals(r1.getVertex(2), S2LatLng.fromRadians(M_PI_2, 0));
    assertEquals(r1.getVertex(3), S2LatLng.fromRadians(M_PI_2, PI));
    // Arguments are reduced modulo 4.
    assertEquals(r1.getVertex(4), S2LatLng.fromRadians(0, PI));
    assertEquals(r1.getVertex(5), S2LatLng.fromRadians(0, 0));
    assertEquals(r1.getVertex(6), S2LatLng.fromRadians(M_PI_2, 0));
    assertEquals(r1.getVertex(7), S2LatLng.fromRadians(M_PI_2, PI));

    assertTrue(r1.contains(S2LatLng.fromDegrees(30, -45)));
    assertTrue(!r1.contains(S2LatLng.fromDegrees(30, 45)));
    assertTrue(!r1.interiorContains(eqM180) && !r1.interiorContains(northPole));
    assertTrue(r1.contains(new S2Point(0.5, -0.3, 0.1)));
    assertTrue(!r1.contains(new S2Point(0.5, 0.2, 0.1)));

    // Make sure that GetVertex() returns vertices in CCW order.
    for (int i = 0; i < 4; ++i) {
      double lat = M_PI_4 * (i - 2);
      double lng = M_PI_2 * (i - 2) + 0.2;
      S2LatLngRect r =
          new S2LatLngRect(
              new R1Interval(lat, lat + M_PI_4),
              new S1Interval(
                  Platform.IEEEremainder(lng, 2 * PI),
                  Platform.IEEEremainder(lng + M_PI_2, 2 * PI)));
      for (int k = 0; k < 4; ++k) {
        assertTrue(
            S2Predicates.sign(
                    r.getVertex((k - 1) & 3).toPoint(),
                    r.getVertex(k).toPoint(),
                    r.getVertex((k + 1) & 3).toPoint())
                > 0);
      }
    }

    // Contains(S2LatLngRect), InteriorContains(S2LatLngRect), Intersects(), InteriorIntersects(),
    // Union(), Intersection().
    //
    // Much more testing of these methods is done in s1interval_unittest and r1interval_unittest.

    S2LatLngRect r1Mid = rectFromDegrees(45, -90, 45, -90);
    S2LatLngRect reqM180 = new S2LatLngRect(eqM180, eqM180);
    S2LatLngRect rNorthPole = new S2LatLngRect(northPole, northPole);

    testIntervalOps(r1, r1Mid, "TTTT", r1, r1Mid);
    testIntervalOps(r1, reqM180, "TFTF", r1, reqM180);
    testIntervalOps(r1, rNorthPole, "TFTF", r1, rNorthPole);

    assertTrue(r1.equals(rectFromDegrees(0, -180, 90, 0)));
    testIntervalOps(
        r1,
        rectFromDegrees(-10, -1, 1, 20),
        "FFTT",
        rectFromDegrees(-10, -180, 90, 20),
        rectFromDegrees(0, -1, 1, 0));
    testIntervalOps(
        r1,
        rectFromDegrees(-10, -1, 0, 20),
        "FFTF",
        rectFromDegrees(-10, -180, 90, 20),
        rectFromDegrees(0, -1, 0, 0));
    testIntervalOps(
        r1,
        rectFromDegrees(-10, 0, 1, 20),
        "FFTF",
        rectFromDegrees(-10, -180, 90, 20),
        rectFromDegrees(0, 0, 1, 0));

    testIntervalOps(
        rectFromDegrees(-15, -160, -15, -150),
        rectFromDegrees(20, 145, 25, 155),
        "FFFF",
        rectFromDegrees(-15, 145, 25, -150),
        empty);
    testIntervalOps(
        rectFromDegrees(70, -10, 90, -140),
        rectFromDegrees(60, 175, 80, 5),
        "FFTT",
        rectFromDegrees(60, -180, 90, 180),
        rectFromDegrees(70, 175, 80, 5));

    // assertTrue that the intersection of two rectangles that overlap in latitude but not longitude
    // is valid, and vice versa.
    testIntervalOps(
        rectFromDegrees(12, 30, 60, 60),
        rectFromDegrees(0, 0, 30, 18),
        "FFFF",
        rectFromDegrees(0, 0, 60, 60),
        empty);
    testIntervalOps(
        rectFromDegrees(0, 0, 18, 42),
        rectFromDegrees(30, 12, 42, 60),
        "FFFF",
        rectFromDegrees(0, 0, 42, 60),
        empty);

    // AddPoint()
    S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();
    builder.addPoint(S2LatLng.fromDegrees(0, 0));
    assertTrue(builder.isPoint());
    builder.addPoint(S2LatLng.fromRadians(0, -M_PI_2));
    assertFalse(builder.isPoint());
    builder.addPoint(S2LatLng.fromRadians(M_PI_4, -PI));
    builder.addPoint(new S2Point(0, 0, 1));
    assertTrue(builder.build().equals(r1));

    // Expanded()
    assertTrue(
        rectFromDegrees(70, 150, 80, 170)
            .expanded(S2LatLng.fromDegrees(20, 30))
            .approxEquals(rectFromDegrees(50, 120, 90, -160)));
    assertTrue(S2LatLngRect.empty().expanded(S2LatLng.fromDegrees(20, 30)).isEmpty());
    assertTrue(S2LatLngRect.full().expanded(S2LatLng.fromDegrees(20, 30)).isFull());
    assertTrue(
        rectFromDegrees(-90, 170, 10, 20)
            .expanded(S2LatLng.fromDegrees(30, 80))
            .approxEquals(rectFromDegrees(-90, -180, 40, 180)));

    // ConvolveWithCap()
    S2LatLngRect.Builder llr1 =
        new S2LatLngRect.Builder(S2LatLng.fromDegrees(0, 170), S2LatLng.fromDegrees(0, -170));
    llr1.convolveWithCap(S1Angle.degrees(15));
    S2LatLngRect llr2 =
        new S2LatLngRect(S2LatLng.fromDegrees(-15, 155), S2LatLng.fromDegrees(15, -155));
    assertTrue(llr1.approxEquals(llr2));

    llr1 = new S2LatLngRect.Builder(S2LatLng.fromDegrees(60, 150), S2LatLng.fromDegrees(80, 10));
    llr1.convolveWithCap(S1Angle.degrees(15));
    llr2 = new S2LatLngRect(S2LatLng.fromDegrees(45, -180), S2LatLng.fromDegrees(90, 180));
    assertTrue(llr1.approxEquals(llr2));

    // GetCapBound(), bounding cap at center is smaller:
    assertTrue(
        new S2LatLngRect(S2LatLng.fromDegrees(-45, -45), S2LatLng.fromDegrees(45, 45))
            .getCapBound()
            .approxEquals(S2Cap.fromAxisHeight(new S2Point(1, 0, 0), 0.5)));
    // GetCapBound(), bounding cap at north pole is smaller:
    assertTrue(
        new S2LatLngRect(S2LatLng.fromDegrees(88, -80), S2LatLng.fromDegrees(89, 80))
            .getCapBound()
            .approxEquals(S2Cap.fromAxisAngle(new S2Point(0, 0, 1), S1Angle.degrees(2))));
    // GetCapBound(), longitude span > 180 degrees:
    assertTrue(
        new S2LatLngRect(S2LatLng.fromDegrees(-30, -150), S2LatLng.fromDegrees(-10, 50))
            .getCapBound()
            .approxEquals(S2Cap.fromAxisAngle(new S2Point(0, 0, -1), S1Angle.degrees(80))));

    // GetCapBound(), ensure hemispheres are bounded conservatively:
    assertTrue(
        new S2LatLngRect(S2LatLng.fromDegrees(-10, -100), S2LatLng.fromDegrees(0, 100))
            .getCapBound().radius().getLength2() >= S1ChordAngle.RIGHT.getLength2());

    // Contains(S2Cell), MayIntersect(S2Cell), Intersects(S2Cell)

    // Special cases.
    testCellOps(empty, S2Cell.fromFace(3), 0);
    testCellOps(full, S2Cell.fromFacePosLevel(2, 0, 0), 4);
    testCellOps(full, S2Cell.fromFacePosLevel(5, 0, 25), 4);

    // This rectangle includes the first quadrant of face 0. It's expanded slightly because cell
    // bounding rectangles are slightly conservative.
    S2LatLngRect r4 = rectFromDegrees(-45.1, -45.1, 0.1, 0.1);
    testCellOps(r4, S2Cell.fromFace(0), 3);
    testCellOps(r4, S2Cell.fromFacePosLevel(0, 0, 1), 4);
    testCellOps(r4, S2Cell.fromFacePosLevel(1, 0, 1), 0);

    // This rectangle intersects the first quadrant of face 0.
    S2LatLngRect r5 = rectFromDegrees(-10, -45, 10, 0);
    testCellOps(r5, S2Cell.fromFace(0), 3);
    testCellOps(r5, S2Cell.fromFacePosLevel(0, 0, 1), 3);
    testCellOps(r5, S2Cell.fromFacePosLevel(1, 0, 1), 0);

    // Rectangle consisting of a single point.
    testCellOps(rectFromDegrees(4, 4, 4, 4), S2Cell.fromFace(0), 3);

    // Rectangles that intersect the bounding rectangle of a face but not the face itself.
    testCellOps(rectFromDegrees(41, -87, 42, -79), S2Cell.fromFace(2), 1);
    testCellOps(rectFromDegrees(-41, 160, -40, -160), S2Cell.fromFace(5), 1);
    {
      // This is the leaf cell at the top right hand corner of face 0. It has two angles of 60
      // degrees and two of 120 degrees.
      S2Cell cell0tr = new S2Cell(new S2Point(1 + 1e-12, 1, 1));
      S2LatLng v0 = new S2LatLng(cell0tr.getVertexRaw(0));
      testCellOps(
          rectFromDegrees(
              v0.lat().degrees() - 1e-8,
              v0.lng().degrees() - 1e-8,
              v0.lat().degrees() - 2e-10,
              v0.lng().degrees() + 1e-10),
          cell0tr,
          1);
    }

    // Rectangles that intersect a face but where no vertex of one region is contained by the other
    // region. The first one passes through a corner of one of the face cells.
    testCellOps(rectFromDegrees(-37, -70, -36, -20), S2Cell.fromFace(5), 2);
    {
      // These two intersect like a diamond and a square.
      S2Cell cell202 = S2Cell.fromFacePosLevel(2, 0, 2);
      S2LatLngRect bound202 = cell202.getRectBound();
      testCellOps(
          rectFromDegrees(
              bound202.lo().lat().degrees() + 3,
              bound202.lo().lng().degrees() + 3,
              bound202.hi().lat().degrees() - 3,
              bound202.hi().lng().degrees() - 3),
          cell202,
          2);
    }
  }

  public void testBoundaryIntersects_EmptyRectangle() {
    S2LatLngRect rect = S2LatLngRect.empty();
    S2Point lo = rect.lo().toPoint();
    S2Point hi = rect.hi().toPoint();
    assertFalse(rect.boundaryIntersects(lo, lo));
    assertFalse(rect.boundaryIntersects(lo, hi));
  }

  public void testBoundaryIntersects_FullRectangle() {
    S2LatLngRect rect = S2LatLngRect.full();
    S2Point lo = rect.lo().toPoint();
    S2Point hi = rect.hi().toPoint();
    assertFalse(rect.boundaryIntersects(lo, lo));
    assertFalse(rect.boundaryIntersects(lo, hi));
  }

  public void testBoundaryIntersects_SphericalLune() {
    // This rectangle only has two non-degenerate sides.
    S2LatLngRect rect = rectFromDegrees(-90, 100, 90, 120);
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("60:60"), S2TextFormat.makePointOrDie("90:60")));
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("-60:110"), S2TextFormat.makePointOrDie("60:110")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("-60:95"), S2TextFormat.makePointOrDie("60:110")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("60:115"), S2TextFormat.makePointOrDie("80:125")));
  }

  public void testBoundaryIntersects_NorthHemisphere() {
    // This rectangle only has only one non-degenerate side.
    S2LatLngRect rect = rectFromDegrees(0, -180, 90, 180);
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("60:-180"), S2TextFormat.makePointOrDie("90:-180")));
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("60:-170"), S2TextFormat.makePointOrDie("60:170")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("-10:-180"), S2TextFormat.makePointOrDie("10:-180")));
  }

  public void testBoundaryIntersects_SouthHemisphere() {
    // This rectangle only has only one non-degenerate side.
    S2LatLngRect rect = rectFromDegrees(-90, -180, 0, 180);
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("-90:-180"), S2TextFormat.makePointOrDie("-60:-180")));
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("-60:-170"), S2TextFormat.makePointOrDie("-60:170")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("-10:-180"), S2TextFormat.makePointOrDie("10:-180")));
  }

  public void testBoundaryIntersects_RectCrossingAntiMeridian() {
    S2LatLngRect rect = rectFromDegrees(20, 170, 40, -170);
    assertTrue(rect.contains(S2TextFormat.makePointOrDie("30:180")));

    // Check that crossings of all four sides are detected.
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("25:160"), S2TextFormat.makePointOrDie("25:180")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("25:-160"), S2TextFormat.makePointOrDie("25:-180")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("15:175"), S2TextFormat.makePointOrDie("30:175")));
    assertTrue(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("45:175"), S2TextFormat.makePointOrDie("30:175")));

    // Check that the edges on the opposite side of the sphere but at the same
    // latitude do not intersect the rectangle boundary.
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("25:-20"), S2TextFormat.makePointOrDie("25:0")));
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("25:20"), S2TextFormat.makePointOrDie("25:0")));
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("15:-5"), S2TextFormat.makePointOrDie("30:-5")));
    assertFalse(
        rect.boundaryIntersects(
            S2TextFormat.makePointOrDie("45:-5"), S2TextFormat.makePointOrDie("30:-5")));
  }

  public void testExpandedByDistancePositive() {
    assertTrue(rectFromDegrees(0, 170, 0, -170).expandedByDistance(S1Angle.degrees(15))
        .approxEquals(rectFromDegrees(-15, 155, 15, -155)));
    assertTrue(rectFromDegrees(60, 150, 80, 10).expandedByDistance(S1Angle.degrees(15))
        .approxEquals(rectFromDegrees(45, -180, 90, 180)));
  }

  public void testExpandedByDistanceContractNorthEast() {
    S2LatLngRect inRect = rectFromDegrees(0, 0, 30, 90);
    S1Angle distance = S1Angle.degrees(5);
    S2LatLngRect outRect = inRect.expandedByDistance(distance).expandedByDistance(distance.mul(-1));
    assertTrue(outRect.toString(), outRect.approxEquals(inRect));
  }

  public void testExpandedByDistanceContractSouthWest() {
    S2LatLngRect inRect = rectFromDegrees(-30, -90, 0, 0);
    S1Angle distance = S1Angle.degrees(5);
    S2LatLngRect outRect = inRect.expandedByDistance(distance).expandedByDistance(distance.mul(-1));
    assertTrue(outRect.toString(), outRect.approxEquals(inRect));
  }

  public void testExpandedByDistanceContractLatWithNorthPole() {
    S2LatLngRect rect = rectFromDegrees(0, -90, 90, 180).expandedByDistance(S1Angle.degrees(-5));
    assertTrue(rect.toString(), rect.approxEquals(rectFromDegrees(5, 0, 85, 90)));
  }

  public void testExpandedByDistanceContractLatWithNorthPoleAndLngFull() {
    S2LatLngRect rect = rectFromDegrees(0, -180, 90, 180).expandedByDistance(S1Angle.degrees(-5));
    assertTrue(rect.toString(), rect.approxEquals(rectFromDegrees(5, -180, 90, 180)));
  }

  public void testExpandedByDistanceContractLatWithSouthPole() {
    S2LatLngRect rect = rectFromDegrees(-90, -90, 0, 180).expandedByDistance(S1Angle.degrees(-5));
    assertTrue(rect.toString(), rect.approxEquals(rectFromDegrees(-85, 0, -5, 90)));
  }

  public void testExpandedByDistanceContractLatWithSouthPoleAndLngFull() {
    S2LatLngRect rect = rectFromDegrees(-90, -180, 0, 180).expandedByDistance(S1Angle.degrees(-5));
    assertTrue(rect.toString(), rect.approxEquals(rectFromDegrees(-90, -180, -5, 180)));
  }

  public void testExpandedByDistanceContractLngFull() {
    S2LatLngRect rect = rectFromDegrees(0, -180, 30, 180).expandedByDistance(S1Angle.degrees(-5));
    assertTrue(rect.toString(), rect.approxEquals(rectFromDegrees(5, -180, 25, 180)));
  }

  public void testExpandedByDistanceContractLatResultEmpty() {
    S2LatLngRect rect = rectFromDegrees(0, 0, 9.9, 90).expandedByDistance(S1Angle.degrees(-5));
    assertTrue(rect.toString(), rect.isEmpty());
  }

  public void testExpandedByDistanceContractLngResultEmpty() {
    S2LatLngRect rect = rectFromDegrees(0, 0, 30, 11).expandedByDistance(S1Angle.degrees(-5));

    // The cap center is at latitude 30 - 5 = 25 degrees. The length of the latitude 25 degree line
    // is 0.906 times the length of the equator. Thus the cap whose radius is 5 degrees covers the
    // rectangle whose latitude interval is 11 degrees.
    assertTrue(rect.toString(), rect.isEmpty());
  }

  public void testPolarClosure() {
    assertEquals(rectFromDegrees(-89, 0, 89, 1), rectFromDegrees(-89, 0, 89, 1).polarClosure());
    assertEquals(
        rectFromDegrees(-90, -180, -45, 180), rectFromDegrees(-90, -30, -45, 100).polarClosure());
    assertEquals(
        rectFromDegrees(89, -180, 90, 180), rectFromDegrees(89, 145, 90, 146).polarClosure());
    assertEquals(S2LatLngRect.full(), rectFromDegrees(-90, -145, 90, -144).polarClosure());
  }

  public void testApproxEquals() {
    // S1Interval and R1Interval have additional testing.

    assertTrue(S2LatLngRect.empty().approxEquals(rectFromDegrees(1, 5, 1, 5)));
    assertTrue(rectFromDegrees(1, 5, 1, 5).approxEquals(S2LatLngRect.empty()));
    assertFalse(rectFromDegrees(1, 5, 1, 5).approxEquals(rectFromDegrees(2, 7, 2, 7)));

    // Test the maxError (double) parameter.
    assertTrue(
        rectFromDegrees(10, 10, 20, 20)
            .approxEquals(rectFromDegrees(11, 11, 19, 19), S1Angle.degrees(1.001).radians()));
    assertFalse(
        rectFromDegrees(10, 10, 20, 20)
            .approxEquals(rectFromDegrees(11, 11, 19, 19), S1Angle.degrees(0.999).radians()));

    // Test the maxError (S2LatLng) parameter.
    assertTrue(
        rectFromDegrees(0, 10, 20, 30)
            .approxEquals(rectFromDegrees(-1, 8, 21, 32), S2LatLng.fromDegrees(1.001, 2.001)));
    assertFalse(
        rectFromDegrees(0, 10, 20, 30)
            .approxEquals(rectFromDegrees(-1, 8, 21, 32), S2LatLng.fromDegrees(0.999, 1.999)));
  }

  public void testArea() {
    assertExactly(0.0, S2LatLngRect.empty().area());
    assertDoubleEquals(4 * PI, S2LatLngRect.full().area());
    assertDoubleEquals(M_PI_2, rectFromDegrees(0, 0, 90, 90).area());
  }

  public void testEdgeBound() {
    // assertTrue cases where min/max latitude is not at a vertex.
    assertDoubleEquals(getEdgeBound(1, 1, 1, 1, -1, 1).lat().hi(), M_PI_4); // Max, CW
    assertDoubleEquals(getEdgeBound(1, -1, 1, 1, 1, 1).lat().hi(), M_PI_4); // Max, CCW
    assertDoubleEquals(getEdgeBound(1, -1, -1, -1, -1, -1).lat().lo(), -M_PI_4); // Min, CW
    assertDoubleEquals(getEdgeBound(-1, 1, -1, -1, -1, -1).lat().lo(), -M_PI_4); // Min, CCW

    // assertTrue cases where the edge passes through one of the poles.
    assertDoubleEquals(getEdgeBound(.3, .4, 1, -.3, -.4, 1).lat().hi(), M_PI_2);
    assertDoubleEquals(getEdgeBound(.3, .4, -1, -.3, -.4, -1).lat().lo(), -M_PI_2);

    // assertTrue cases where the min/max latitude is attained at a vertex.
    final double kCubeLat = asin(sqrt(1. / 3)); // 35.26 degrees
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
          S2LatLngRect.fromPointPair(
              new S2LatLng(data.getRandomPoint()), new S2LatLng(data.getRandomPoint()));
      S2LatLngRect b =
          S2LatLngRect.fromPointPair(
              new S2LatLng(data.getRandomPoint()), new S2LatLng(data.getRandomPoint()));
      verifyGetDistance(a, b);

      S2LatLng c = new S2LatLng(data.getRandomPoint());
      verifyGetRectPointDistance(a, c);
      verifyGetRectPointDistance(b, c);
    }
  }

  /**
   * This function assumes that GetDirectedHausdorffDistance() always returns a distance from some
   * point in a to b. So the function mainly tests whether the returned distance is large enough,
   * and only does a weak test on whether it is small enough.
   */
  private static void verifyGetDirectedHausdorffDistance(S2LatLngRect a, S2LatLngRect b) {
    S1Angle hausdorffDistance = a.getDirectedHausdorffDistance(b);

    final double resolution = 0.1;
    S1Angle maxDistance = new S1Angle();

    int sampleSizeOnLat = (int) (a.lat().getLength() / resolution) + 1;
    int sampleSizeOnLng = (int) (a.lng().getLength() / resolution) + 1;
    double deltaOnLat = a.lat().getLength() / sampleSizeOnLat;
    double deltaOnLng = a.lng().getLength() / sampleSizeOnLng;

    double lng = a.lng().lo();
    for (int i = 0; i <= sampleSizeOnLng; ++i, lng += deltaOnLng) {
      double lat = a.lat().lo();
      for (int j = 0; j <= sampleSizeOnLat; ++j, lat += deltaOnLat) {
        S2LatLng latlng = S2LatLng.fromRadians(lat, lng).normalized();
        S1Angle distanceToB = b.getDistance(latlng);

        if (distanceToB.greaterOrEquals(maxDistance)) {
          maxDistance = distanceToB;
        }
      }
    }

    assertTrue(a + ":" + b, maxDistance.radians() <= hausdorffDistance.radians() + 1e-10);
    assertTrue(a + ":" + b, maxDistance.radians() >= hausdorffDistance.radians() - resolution);
  }

  public void testGetDirectedHausdorffDistanceRandomPairs() {
    // Test random pairs.
    final int iterations = 1000;
    for (int i = 0; i < iterations; ++i) {
      S2LatLngRect a =
          S2LatLngRect.fromPointPair(
              new S2LatLng(data.getRandomPoint()), new S2LatLng(data.getRandomPoint()));
      S2LatLngRect b =
          S2LatLngRect.fromPointPair(
              new S2LatLng(data.getRandomPoint()), new S2LatLng(data.getRandomPoint()));

      // a and b are *minimum* bounding rectangles of two random points, in particular, their
      // Voronoi diagrams are always of the same topology. We take the "complements" of a and b for
      // more thorough testing.
      S2LatLngRect a2 = new S2LatLngRect(a.lat(), a.lng().complement());
      S2LatLngRect b2 = new S2LatLngRect(b.lat(), b.lng().complement());

      // Note that "a" and "b" come from the same distribution, so there is no need to test pairs
      // such as (b, a), (b, a2), etc.
      verifyGetDirectedHausdorffDistance(a, b);
      verifyGetDirectedHausdorffDistance(a, b2);
      verifyGetDirectedHausdorffDistance(a2, b);
      verifyGetDirectedHausdorffDistance(a2, b2);
    }
  }

  public void testGetDirectedHausdorffDistanceContained() {
    // Caller rect is contained in callee rect. Should return 0.
    S2LatLngRect a = rectFromDegrees(-10, 20, -5, 90);
    assertEquals(
        S1Angle.radians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-10, 20, -5, 90)));
    assertEquals(
        S1Angle.radians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-10, 19, -5, 91)));
    assertEquals(
        S1Angle.radians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-11, 20, -4, 90)));
    assertEquals(
        S1Angle.radians(0), a.getDirectedHausdorffDistance(rectFromDegrees(-11, 19, -4, 91)));
  }

  public void testGetDirectHausdorffDistancePointToRect() {
    // The Hausdorff distance from a point to a rect should be the same as its distance to the rect.
    S2LatLngRect a1 = pointRectFromDegrees(5, 8);
    S2LatLngRect a2 = pointRectFromDegrees(90, 10); // north pole

    S2LatLngRect b = rectFromDegrees(-85, -50, -80, 10);
    assertDoubleEquals(a1.getDirectedHausdorffDistance(b).radians(), a1.getDistance(b).radians());
    assertDoubleEquals(a2.getDirectedHausdorffDistance(b).radians(), a2.getDistance(b).radians());

    b = rectFromDegrees(4, -10, 80, 10);
    assertDoubleEquals(a1.getDirectedHausdorffDistance(b).radians(), a1.getDistance(b).radians());
    assertDoubleEquals(a2.getDirectedHausdorffDistance(b).radians(), a2.getDistance(b).radians());

    b = rectFromDegrees(70, 170, 80, -170);
    assertDoubleEquals(a1.getDirectedHausdorffDistance(b).radians(), a1.getDistance(b).radians());
    assertDoubleEquals(a2.getDirectedHausdorffDistance(b).radians(), a2.getDistance(b).radians());
  }

  public void testGetDirectedHausdorffDistanceRectToPoint() {
    S2LatLngRect a = rectFromDegrees(1, -8, 10, 20);
    verifyGetDirectedHausdorffDistance(a, pointRectFromDegrees(5, 8));
    verifyGetDirectedHausdorffDistance(a, pointRectFromDegrees(-6, -100));
    // south pole
    verifyGetDirectedHausdorffDistance(a, pointRectFromDegrees(-90, -20));
    // north pole
    verifyGetDirectedHausdorffDistance(a, pointRectFromDegrees(90, 0));
  }

  public void testGetDirectedHausdorffDistanceRectToRectNearPole() {
    // Tests near south pole.
    S2LatLngRect a = rectFromDegrees(-87, 0, -85, 3);
    verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-89, 1, -88, 2));
    verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, 1, -83, 2));
    verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-88, 90, -86, 91));
    verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, -91, -83, -90));
    verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-90, 181, -89, 182));
    verifyGetDirectedHausdorffDistance(a, rectFromDegrees(-84, 181, -83, 182));
  }

  public void testGetDirectedHausdorffDistanceRectToRectDegenerateCases() {
    // Rectangles that contain poles.
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(0, 10, 90, 20), rectFromDegrees(-4, -10, 4, 0));
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(-4, -10, 4, 0), rectFromDegrees(0, 10, 90, 20));

    // Two rectangles share same or complement longitudinal intervals.
    S2LatLngRect a = rectFromDegrees(-50, -10, 50, 10);
    S2LatLngRect b = rectFromDegrees(30, -10, 60, 10);
    verifyGetDirectedHausdorffDistance(a, b);
    S2LatLngRect c = new S2LatLngRect(a.lat(), a.lng().complement());
    verifyGetDirectedHausdorffDistance(c, b);

    // rectangle a touches b_opposite_lng.
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(10, 170, 30, 180), rectFromDegrees(-50, -10, 50, 10));
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(10, -180, 30, -170), rectFromDegrees(-50, -10, 50, 10));

    // rectangle b's Voronoi diagram is degenerate (lng interval spans 180 degrees), and a touches
    // the degenerate Voronoi vertex.
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(-30, 170, 30, 180), rectFromDegrees(-10, -90, 10, 90));
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(-30, -180, 30, -170), rectFromDegrees(-10, -90, 10, 90));

    // rectangle a touches a voronoi vertex of rectangle b.
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(-20, 105, 20, 110), rectFromDegrees(-30, 5, 30, 15));
    verifyGetDirectedHausdorffDistance(
        rectFromDegrees(-20, 95, 20, 105), rectFromDegrees(-30, 5, 30, 15));
  }

  public void testConvolveWithCap() {
    // ConvolveWithCap()
    S2LatLngRect.Builder llr1 =
        new S2LatLngRect.Builder(S2LatLng.fromDegrees(0, 170), S2LatLng.fromDegrees(0, -170));
    llr1.convolveWithCap(S1Angle.degrees(15));
    S2LatLngRect llr2 =
        new S2LatLngRect(S2LatLng.fromDegrees(-15, 155), S2LatLng.fromDegrees(15, -155));
    assertTrue(llr1.approxEquals(llr2));

    llr1 = new S2LatLngRect.Builder(S2LatLng.fromDegrees(60, 150), S2LatLng.fromDegrees(80, 10));
    llr1.convolveWithCap(S1Angle.degrees(15));
    llr2 = new S2LatLngRect(S2LatLng.fromDegrees(45, -180), S2LatLng.fromDegrees(90, 180));
    assertTrue(llr1.approxEquals(llr2));
  }

  public void testEncodeDecode() throws IOException {
    S2LatLngRect r = rectFromDegrees(-20, -80, 10, 20);
    ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
    r.encode(outputStream);
    byte[] encodedBytes = outputStream.toByteArray();
    assertTrue(encodedBytes.length > 0);
    S2LatLngRect decoded = S2LatLngRect.decode(new ByteArrayInputStream(encodedBytes));
    assertEquals(r, decoded);
  }

  public void testCoder() {
    S2LatLngRect value = S2LatLngRect.fromEdge(S2Point.X_POS, S2Point.Y_POS);
    assertEquals(value, roundtrip(S2LatLngRect.CODER, value));
  }

  // Recursively verify that when a rectangle is split into two pieces, the centroids of the
  // children sum to give the centroid of the parent.
  void checkCentroidSplitting(S2LatLngRect r, int splitsLeft) {
    S2LatLngRect child0;
    S2LatLngRect child1;
    if (data.oneIn(2)) {
      double lat = data.uniform(r.lat().lo(), r.lat().hi());
      child0 = new S2LatLngRect(new R1Interval(r.lat().lo(), lat), r.lng());
      child1 = new S2LatLngRect(new R1Interval(lat, r.lat().hi()), r.lng());
    } else {
      assertTrue(r.lng().lo() < r.lng().hi());
      double lng = data.uniform(r.lng().lo(), r.lng().hi());
      child0 = new S2LatLngRect(r.lat(), new S1Interval(r.lng().lo(), lng));
      child1 = new S2LatLngRect(r.lat(), new S1Interval(lng, r.lng().hi()));
    }
    assertTrue(
        (r.getCentroid().sub(child0.getCentroid()).sub(child1.getCentroid()).norm() <= 1e-15));
    if (splitsLeft > 0) {
      checkCentroidSplitting(child0, splitsLeft - 1);
      checkCentroidSplitting(child1, splitsLeft - 1);
    }
  }

  public void testGetCentroid() {
    // Empty and full rectangles.
    assertEquals(S2Point.ORIGIN, S2LatLngRect.empty().getCentroid());
    assertTrue(S2LatLngRect.full().getCentroid().norm() <= 1e-15);

    // Rectangles that cover the full longitude range.
    for (int i = 0; i < 100; ++i) {
      double lat1 = data.uniform(-M_PI_2, M_PI_2);
      double lat2 = data.uniform(-M_PI_2, M_PI_2);
      S2LatLngRect r = new S2LatLngRect(R1Interval.fromPointPair(lat1, lat2), S1Interval.full());
      S2Point centroid = r.getCentroid();
      assertEquals(0.5 * (sin(lat1) + sin(lat2)) * r.area(), centroid.getZ(), 1e-15);
      assertTrue(new R2Vector(centroid.getX(), centroid.getY()).norm() <= 1e-15);
    }

    // Rectangles that cover the full latitude range.
    for (int i = 0; i < 100; ++i) {
      double lng1 = data.uniform(-PI, PI);
      double lng2 = data.uniform(-PI, PI);
      S2LatLngRect r =
          new S2LatLngRect(S2LatLngRect.fullLat(), S1Interval.fromPointPair(lng1, lng2));
      S2Point centroid = r.getCentroid();
      assertTrue(abs(centroid.getZ()) <= 1e-15);
      assertEquals(r.lng().getCenter(), new S2LatLng(centroid).lng().radians(), 1e-15);
      double alpha = 0.5 * r.lng().getLength();
      assertEquals(
          0.25 * PI * sin(alpha) / alpha * r.area(),
          new R2Vector(centroid.getX(), centroid.getY()).norm(),
          1e-15);
    }

    // Finally, verify that when a rectangle is recursively split into pieces, the centroids of the
    // pieces add to give the centroid of their parent. To make the code simpler we avoid rectangles
    // that cross the 180 degree line of longitude.
    checkCentroidSplitting(
        new S2LatLngRect(S2LatLngRect.fullLat(), new S1Interval(-3.14, 3.14)), 10 /*splits_left*/);
  }

  public void testBuilder_empty() {
    assertEquals(S2LatLngRect.empty(), S2LatLngRect.Builder.empty().build());
  }

  public void testBuilder_full() {
    assertEquals(
        S2LatLngRect.full(),
        S2LatLngRect.Builder.empty().addPoint(new S2Point(0.3, 0.2, 0.1)).setFull().build());
  }

  public void testBuilder_fromPoints() {
    S2LatLngRect.Builder builder =
        new S2LatLngRect.Builder(S2LatLng.fromDegrees(-1, 2), S2LatLng.fromDegrees(1, 3));
    builder.addPoint(S2LatLng.fromDegrees(0.25, 6.5)); // São Tomé
    builder.addPoint(S2LatLng.fromDegrees(5.6, -0.17)); // Kotoka airport in Accra
    assertEquals(rectFromDegrees(-1, -0.17, 5.6, 6.5), builder.build());
  }

  public void testBuilder_union() {
    S2LatLngRect.Builder builder =
        S2LatLngRect.Builder.empty()
            .union(rectFromDegrees(-2, -176, -1, -175))
            .union(rectFromDegrees(3, 170, 4, 171));
    S2LatLngRect rect1 = builder.build();
    builder.union(rectFromDegrees(10, 179, 11, -179));
    S2LatLngRect rect2 = builder.build();
    assertEquals(rectFromDegrees(-2, 170, 4, -175), rect1);
    assertEquals(rectFromDegrees(-2, 170, 11, -175), rect2);
  }

  public void testBuilder_intersection() {
    S2LatLngRect.Builder builder =
        rectFromDegrees(14.7, -17.5, 31.8, 34.2) // Sahara
            .toBuilder()
            .intersection(rectFromDegrees(-2.7, 29.3, 31.5, 38.6)); // Nile
    S2LatLngRect rect1 = builder.build();
    builder.intersection(rectFromDegrees(15.4, 32.4, 15.7, 32.7)); // Khartoum
    S2LatLngRect rect2 = builder.build();
    assertEquals(rectFromDegrees(14.7, 29.3, 31.5, 34.2), rect1);
    assertEquals(rectFromDegrees(15.4, 32.4, 15.7, 32.7), rect2);
  }

  private static S1Angle bruteForceDistance(S2LatLngRect a, S2LatLngRect b) {
    if (a.intersects(b)) {
      return S1Angle.radians(0);
    }

    // Compare every point in 'a' against every latitude edge and longitude edge in 'b', and
    // vice-versa, for a total of 16 point-vs-latitude-edge tests and 16 point-vs-longitude-edge
    // tests.
    S2LatLng[] pntA = {
      new S2LatLng(a.latLo(), a.lngLo()),
      new S2LatLng(a.latLo(), a.lngHi()),
      new S2LatLng(a.latHi(), a.lngHi()),
      new S2LatLng(a.latHi(), a.lngLo())
    };
    S2LatLng[] pntB = {
      new S2LatLng(b.latLo(), b.lngLo()),
      new S2LatLng(b.latLo(), b.lngHi()),
      new S2LatLng(b.latHi(), b.lngHi()),
      new S2LatLng(b.latHi(), b.lngLo())
    };

    // Make arrays containing the lo/hi latitudes and the lo/hi longitude edges.
    S1Angle[] latA = {a.latLo(), a.latHi()};
    S1Angle[] latB = {b.latLo(), b.latHi()};
    S2Point[][] lngEdgeA = {
      {pntA[0].toPoint(), pntA[3].toPoint()}, {pntA[1].toPoint(), pntA[2].toPoint()}
    };
    S2Point[][] lngEdgeB = {
      {pntB[0].toPoint(), pntB[3].toPoint()}, {pntB[1].toPoint(), pntB[2].toPoint()}
    };

    S1Angle minDistance = S1Angle.degrees(180.0);
    for (int i = 0; i < 4; ++i) {
      // For each point in a and b.
      S2LatLng currentA = pntA[i];
      S2LatLng currentB = pntB[i];

      for (int j = 0; j < 2; ++j) {
        // Get distances to latitude and longitude edges.
        S1Angle aToLat = getDistance(currentA, latB[j], b.lng());
        S1Angle bToLat = getDistance(currentB, latA[j], a.lng());
        S1Angle aToLng = S2EdgeUtil.getDistance(currentA.toPoint(), lngEdgeB[j][0], lngEdgeB[j][1]);
        S1Angle bToLng = S2EdgeUtil.getDistance(currentB.toPoint(), lngEdgeA[j][0], lngEdgeA[j][1]);

        minDistance =
            S1Angle.min(
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
        S2EdgeUtil.getDistance(
            b.toPoint(),
            new S2LatLng(a.latLo(), a.lngLo()).toPoint(),
            new S2LatLng(a.latHi(), a.lngLo()).toPoint());
    S1Angle bToHiLng =
        S2EdgeUtil.getDistance(
            b.toPoint(),
            new S2LatLng(a.latLo(), a.lngHi()).toPoint(),
            new S2LatLng(a.latHi(), a.lngHi()).toPoint());
    return S1Angle.min(bToLoLat, S1Angle.min(bToHiLat, S1Angle.min(bToLoLng, bToHiLng)));
  }

  /**
   * Returns the minimum distance from X to the latitude line segment defined by the given latitude
   * and longitude interval.
   */
  private static S1Angle getDistance(S2LatLng x, S1Angle lat, S1Interval interval) {
    assertTrue(x.isValid());
    assertTrue(interval.isValid());

    // Is X inside the longitude interval?
    if (interval.contains(x.lng().radians())) {
      return S1Angle.radians(abs(x.lat().radians() - lat.radians()));
    }

    // Return the distance to the closer endpoint.
    return S1Angle.min(
        x.getDistance(new S2LatLng(lat, S1Angle.radians(interval.lo()))),
        x.getDistance(new S2LatLng(lat, S1Angle.radians(interval.hi()))));
  }

  /** Returns the S2LatLngRect from the edge (x1,y1,z1).normalize() to (x2,y2,z2).normalize. */
  private static S2LatLngRect getEdgeBound(
      double x1, double y1, double z1, double x2, double y2, double z2) {
    return S2LatLngRect.fromEdge(
        new S2Point(x1, y1, z1).normalize(), new S2Point(x2, y2, z2).normalize());
  }

  private static S2LatLngRect pointRectFromDegrees(double lat, double lng) {
    return S2LatLngRect.fromPoint(S2LatLng.fromDegrees(lat, lng).normalized());
  }

  /**
   * Convenience method to construct a rectangle. This method is intentionally *not* in the
   * S2LatLngRect interface because the argument order is ambiguous, but hopefully it's not too
   * confusing within the context of this unit test.
   */
  private static S2LatLngRect rectFromDegrees(
        double latLo, double lngLo, double latHi, double lngHi) {

    return new S2LatLngRect(
        S2LatLng.fromDegrees(latLo, lngLo).normalized(),
        S2LatLng.fromDegrees(latHi, lngHi).normalized());
  }

  /**
   * This method verifies a.getDistance(b), where b is a S2LatLng, by comparing its result against
   * a.getDistance(c), c being the point rectangle created from b.
   */
  private static void verifyGetRectPointDistance(S2LatLngRect a, S2LatLng p) {
    S1Angle distance1 = bruteForceRectPointDistance(a, p.normalized());
    S1Angle distance2 = a.getDistance(p.normalized());
    assertEquals(distance1.radians(), distance2.radians(), 1e-10);
  }

  /**
   * This method verifies a.getDistance(b) by comparing its result against a brute-force
   * implementation. The correctness of the brute-force version is much easier to verify by
   * inspection.
   */
  private static void verifyGetDistance(S2LatLngRect a, S2LatLngRect b) {
    S1Angle distance1 = bruteForceDistance(a, b);
    S1Angle distance2 = a.getDistance(b);
    assertEquals(distance1.radians(), distance2.radians(), 1e-10);
  }
}
