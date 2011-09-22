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

public strictfp class S2CapTest extends GeometryTestCase {

  public S2Point getLatLngPoint(double latDegrees, double lngDegrees) {
    return S2LatLng.fromDegrees(latDegrees, lngDegrees).toPoint();
  }

  // About 9 times the double-precision roundoff relative error.
  public static final double EPS = 1e-15;

  public void testBasic() {
    // Test basic properties of empty and full caps.
    S2Cap empty = S2Cap.empty();
    S2Cap full = S2Cap.full();
    assertTrue(empty.isValid());
    assertTrue(empty.isEmpty());
    assertTrue(empty.complement().isFull());
    assertTrue(full.isValid());
    assertTrue(full.isFull());
    assertTrue(full.complement().isEmpty());
    assertEquals(full.height(), 2.0);
    assertDoubleNear(full.angle().degrees(), 180);

    // Containment and intersection of empty and full caps.
    assertTrue(empty.contains(empty));
    assertTrue(full.contains(empty));
    assertTrue(full.contains(full));
    assertTrue(!empty.interiorIntersects(empty));
    assertTrue(full.interiorIntersects(full));
    assertTrue(!full.interiorIntersects(empty));

    // Singleton cap containing the x-axis.
    S2Cap xaxis = S2Cap.fromAxisHeight(new S2Point(1, 0, 0), 0);
    assertTrue(xaxis.contains(new S2Point(1, 0, 0)));
    assertTrue(!xaxis.contains(new S2Point(1, 1e-20, 0)));
    assertEquals(xaxis.angle().radians(), 0.0);

    // Singleton cap containing the y-axis.
    S2Cap yaxis = S2Cap.fromAxisAngle(new S2Point(0, 1, 0), S1Angle.radians(0));
    assertTrue(!yaxis.contains(xaxis.axis()));
    assertEquals(xaxis.height(), 0.0);

    // Check that the complement of a singleton cap is the full cap.
    S2Cap xcomp = xaxis.complement();
    assertTrue(xcomp.isValid());
    assertTrue(xcomp.isFull());
    assertTrue(xcomp.contains(xaxis.axis()));

    // Check that the complement of the complement is *not* the original.
    assertTrue(xcomp.complement().isValid());
    assertTrue(xcomp.complement().isEmpty());
    assertTrue(!xcomp.complement().contains(xaxis.axis()));

    // Check that very small caps can be represented accurately.
    // Here "kTinyRad" is small enough that unit vectors perturbed by this
    // amount along a tangent do not need to be renormalized.
    final double kTinyRad = 1e-10;
    S2Cap tiny =
        S2Cap.fromAxisAngle(S2Point.normalize(new S2Point(1, 2, 3)), S1Angle.radians(kTinyRad));
    S2Point tangent = S2Point.normalize(S2Point.crossProd(tiny.axis(), new S2Point(3, 2, 1)));
    assertTrue(tiny.contains(S2Point.add(tiny.axis(), S2Point.mul(tangent, 0.99 * kTinyRad))));
    assertTrue(!tiny.contains(S2Point.add(tiny.axis(), S2Point.mul(tangent, 1.01 * kTinyRad))));

    // Basic tests on a hemispherical cap.
    S2Cap hemi = S2Cap.fromAxisHeight(S2Point.normalize(new S2Point(1, 0, 1)), 1);
    assertEquals(hemi.complement().axis(), S2Point.neg(hemi.axis()));
    assertEquals(hemi.complement().height(), 1.0);
    assertTrue(hemi.contains(new S2Point(1, 0, 0)));
    assertTrue(!hemi.complement().contains(new S2Point(1, 0, 0)));
    assertTrue(hemi.contains(S2Point.normalize(new S2Point(1, 0, -(1 - EPS)))));
    assertTrue(!hemi.interiorContains(S2Point.normalize(new S2Point(1, 0, -(1 + EPS)))));

    // A concave cap.
    S2Cap concave = S2Cap.fromAxisAngle(getLatLngPoint(80, 10), S1Angle.degrees(150));
    assertTrue(concave.contains(getLatLngPoint(-70 * (1 - EPS), 10)));
    assertTrue(!concave.contains(getLatLngPoint(-70 * (1 + EPS), 10)));
    assertTrue(concave.contains(getLatLngPoint(-50 * (1 - EPS), -170)));
    assertTrue(!concave.contains(getLatLngPoint(-50 * (1 + EPS), -170)));

    // Cap containment tests.
    assertTrue(!empty.contains(xaxis));
    assertTrue(!empty.interiorIntersects(xaxis));
    assertTrue(full.contains(xaxis));
    assertTrue(full.interiorIntersects(xaxis));
    assertTrue(!xaxis.contains(full));
    assertTrue(!xaxis.interiorIntersects(full));
    assertTrue(xaxis.contains(xaxis));
    assertTrue(!xaxis.interiorIntersects(xaxis));
    assertTrue(xaxis.contains(empty));
    assertTrue(!xaxis.interiorIntersects(empty));
    assertTrue(hemi.contains(tiny));
    assertTrue(hemi.contains(
        S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.radians(S2.M_PI_4 - EPS))));
    assertTrue(!hemi.contains(
        S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.radians(S2.M_PI_4 + EPS))));
    assertTrue(concave.contains(hemi));
    assertTrue(concave.interiorIntersects(hemi.complement()));
    assertTrue(!concave.contains(S2Cap.fromAxisHeight(S2Point.neg(concave.axis()), 0.1)));
  }

  public void testRectBound() {
    // Empty and full caps.
    assertTrue(S2Cap.empty().getRectBound().isEmpty());
    assertTrue(S2Cap.full().getRectBound().isFull());

    final double kDegreeEps = 1e-13;
    // Maximum allowable error for latitudes and longitudes measured in
    // degrees. (assertDoubleNear uses a fixed tolerance that is too small.)

    // Cap that includes the south pole.
    S2LatLngRect rect =
        S2Cap.fromAxisAngle(getLatLngPoint(-45, 57), S1Angle.degrees(50)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), -90, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 5, kDegreeEps);
    assertTrue(rect.lng().isFull());

    // Cap that is tangent to the north pole.
    rect = S2Cap.fromAxisAngle(S2Point.normalize(new S2Point(1, 0, 1)), S1Angle.radians(S2.M_PI_4))
        .getRectBound();
    assertDoubleNear(rect.lat().lo(), 0);
    assertDoubleNear(rect.lat().hi(), S2.M_PI_2);
    assertTrue(rect.lng().isFull());

    rect = S2Cap
        .fromAxisAngle(S2Point.normalize(new S2Point(1, 0, 1)), S1Angle.degrees(45)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), 0, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 90, kDegreeEps);
    assertTrue(rect.lng().isFull());

    // The eastern hemisphere.
    rect = S2Cap
        .fromAxisAngle(new S2Point(0, 1, 0), S1Angle.radians(S2.M_PI_2 + 5e-16)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), -90, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 90, kDegreeEps);
    assertTrue(rect.lng().isFull());

    // A cap centered on the equator.
    rect = S2Cap.fromAxisAngle(getLatLngPoint(0, 50), S1Angle.degrees(20)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), -20, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 20, kDegreeEps);
    assertDoubleNear(rect.lngLo().degrees(), 30, kDegreeEps);
    assertDoubleNear(rect.lngHi().degrees(), 70, kDegreeEps);

    // A cap centered on the north pole.
    rect = S2Cap.fromAxisAngle(getLatLngPoint(90, 123), S1Angle.degrees(10)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), 80, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 90, kDegreeEps);
    assertTrue(rect.lng().isFull());
  }

  public void testCells() {
    // For each cube face, we construct some cells on
    // that face and some caps whose positions are relative to that face,
    // and then check for the expected intersection/containment results.

    // The distance from the center of a face to one of its vertices.
    final double kFaceRadius = Math.atan(S2.M_SQRT2);

    for (int face = 0; face < 6; ++face) {
      // The cell consisting of the entire face.
      S2Cell rootCell = S2Cell.fromFacePosLevel(face, (byte) 0, 0);

      // A leaf cell at the midpoint of the v=1 edge.
      S2Cell edgeCell = new S2Cell(S2Projections.faceUvToXyz(face, 0, 1 - EPS));

      // A leaf cell at the u=1, v=1 corner.
      S2Cell cornerCell = new S2Cell(S2Projections.faceUvToXyz(face, 1 - EPS, 1 - EPS));

      // Quick check for full and empty caps.
      assertTrue(S2Cap.full().contains(rootCell));
      assertTrue(!S2Cap.empty().mayIntersect(rootCell));

      // Check intersections with the bounding caps of the leaf cells that are
      // adjacent to 'corner_cell' along the Hilbert curve. Because this corner
      // is at (u=1,v=1), the curve stays locally within the same cube face.
      S2CellId first = cornerCell.id().prev().prev().prev();
      S2CellId last = cornerCell.id().next().next().next().next();
      for (S2CellId id = first; id.lessThan(last); id = id.next()) {
        S2Cell cell = new S2Cell(id);
        assertEquals(cell.getCapBound().contains(cornerCell), id.equals(cornerCell.id()));
        assertEquals(
            cell.getCapBound().mayIntersect(cornerCell), id.parent().contains(cornerCell.id()));
      }

      int antiFace = (face + 3) % 6; // Opposite face.
      for (int capFace = 0; capFace < 6; ++capFace) {
        // A cap that barely contains all of 'cap_face'.
        S2Point center = S2Projections.getNorm(capFace);
        S2Cap covering = S2Cap.fromAxisAngle(center, S1Angle.radians(kFaceRadius + EPS));
        assertEquals(covering.contains(rootCell), capFace == face);
        assertEquals(covering.mayIntersect(rootCell), capFace != antiFace);
        assertEquals(covering.contains(edgeCell), center.dotProd(edgeCell.getCenter()) > 0.1);
        assertEquals(covering.contains(edgeCell), covering.mayIntersect(edgeCell));
        assertEquals(covering.contains(cornerCell), capFace == face);
        assertEquals(
            covering.mayIntersect(cornerCell), center.dotProd(cornerCell.getCenter()) > 0);

        // A cap that barely intersects the edges of 'cap_face'.
        S2Cap bulging = S2Cap.fromAxisAngle(center, S1Angle.radians(S2.M_PI_4 + EPS));
        assertTrue(!bulging.contains(rootCell));
        assertEquals(bulging.mayIntersect(rootCell), capFace != antiFace);
        assertEquals(bulging.contains(edgeCell), capFace == face);
        assertEquals(bulging.mayIntersect(edgeCell), center.dotProd(edgeCell.getCenter()) > 0.1);
        assertTrue(!bulging.contains(cornerCell));
        assertTrue(!bulging.mayIntersect(cornerCell));

        // A singleton cap.
        S2Cap singleton = S2Cap.fromAxisAngle(center, S1Angle.radians(0));
        assertEquals(singleton.mayIntersect(rootCell), capFace == face);
        assertTrue(!singleton.mayIntersect(edgeCell));
        assertTrue(!singleton.mayIntersect(cornerCell));
      }
    }
  }
}
