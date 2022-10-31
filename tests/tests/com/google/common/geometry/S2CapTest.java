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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2.M_PI_4;
import static com.google.common.geometry.S2.M_SQRT2;
import static java.lang.Math.atan;

import com.google.common.io.BaseEncoding;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/** Tests for S2Cap. */
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
    assertExactly(2.0, full.height());
    assertDoubleEquals(full.angle().degrees(), 180);

    // Test the S1Angle constructor using out-of-range arguments.
    assertTrue(S2Cap.fromAxisAngle(S2Point.X_POS, S1Angle.radians(-20)).isEmpty());
    assertTrue(S2Cap.fromAxisAngle(S2Point.X_POS, S1Angle.radians(5)).isFull());
    assertTrue(S2Cap.fromAxisAngle(S2Point.X_POS, S1Angle.INFINITY).isFull());

    // Containment and intersection of empty and full caps.
    assertTrue(empty.contains(empty));
    assertTrue(full.contains(empty));
    assertTrue(full.contains(full));
    assertFalse(empty.interiorIntersects(empty));
    assertTrue(full.interiorIntersects(full));
    assertFalse(full.interiorIntersects(empty));

    // Singleton cap containing the x-axis.
    S2Cap xaxis = S2Cap.fromAxisHeight(new S2Point(1, 0, 0), 0);
    assertTrue(xaxis.contains(new S2Point(1, 0, 0)));
    assertFalse(xaxis.contains(new S2Point(1, 1e-20, 0)));
    assertExactly(0.0, xaxis.angle().radians());

    // Singleton cap containing the y-axis.
    S2Cap yaxis = S2Cap.fromAxisAngle(new S2Point(0, 1, 0), S1Angle.radians(0));
    assertFalse(yaxis.contains(xaxis.axis()));
    assertExactly(0.0, xaxis.height());

    // Check that the complement of a singleton cap is the full cap.
    S2Cap xcomp = xaxis.complement();
    assertTrue(xcomp.isValid());
    assertTrue(xcomp.isFull());
    assertTrue(xcomp.contains(xaxis.axis()));

    // Check that the complement of the complement is *not* the original.
    assertTrue(xcomp.complement().isValid());
    assertTrue(xcomp.complement().isEmpty());
    assertFalse(xcomp.complement().contains(xaxis.axis()));

    // Check that very small caps can be represented accurately. Here "kTinyRad" is small enough
    // that unit vectors perturbed by this amount along a tangent do not need to be renormalized.
    final double kTinyRad = 1e-10;
    S2Cap tiny = S2Cap.fromAxisAngle(new S2Point(1, 2, 3).normalize(), S1Angle.radians(kTinyRad));
    S2Point tangent = tiny.axis().crossProd(new S2Point(3, 2, 1)).normalize();
    assertTrue(tiny.contains(tiny.axis().add(tangent.mul(0.99 * kTinyRad))));
    assertFalse(tiny.contains(tiny.axis().add(tangent.mul(1.01 * kTinyRad))));

    // Basic tests on a hemispherical cap.
    S2Cap hemi = S2Cap.fromAxisHeight(new S2Point(1, 0, 1).normalize(), 1);
    assertEquals(hemi.complement().axis(), hemi.axis().neg());
    assertEquals(1.0, hemi.complement().height(), 0.0);
    assertTrue(hemi.contains(new S2Point(1, 0, 0)));
    assertFalse(hemi.complement().contains(new S2Point(1, 0, 0)));
    assertTrue(hemi.contains(new S2Point(1, 0, -(1 - EPS)).normalize()));
    assertFalse(hemi.interiorContains(new S2Point(1, 0, -(1 + EPS)).normalize()));

    // A concave cap. Note that the error bounds for point containment tests increase with the cap
    // angle, so we need to use a larger error bound here. (It would be painful to do this
    // everywhere, but this at least gives an example of how to compute the maximum error.)
    S2Point center = getLatLngPoint(80, 10);
    S1ChordAngle radius = S1ChordAngle.fromS1Angle(S1Angle.degrees(150));
    double maxError =
        radius.getS2PointConstructorMaxError()
            + radius.getS1AngleConstructorMaxError()
            + 3 * DBL_EPSILON; // getLatLngPoint() error
    S2Cap concave = S2Cap.fromAxisChord(center, radius);
    S2Cap concaveMin = S2Cap.fromAxisChord(center, radius.plusError(-maxError));
    S2Cap concaveMax = S2Cap.fromAxisChord(center, radius.plusError(maxError));
    assertTrue(concaveMax.contains(getLatLngPoint(-70, 10)));
    assertFalse(concaveMin.contains(getLatLngPoint(-70, 10)));
    assertTrue(concaveMax.contains(getLatLngPoint(-50, -170)));
    assertFalse(concaveMin.contains(getLatLngPoint(-50, -170)));

    // Cap containment tests.
    assertFalse(empty.contains(xaxis));
    assertFalse(empty.interiorIntersects(xaxis));
    assertTrue(full.contains(xaxis));
    assertTrue(full.interiorIntersects(xaxis));
    assertFalse(xaxis.contains(full));
    assertFalse(xaxis.interiorIntersects(full));
    assertTrue(xaxis.contains(xaxis));
    assertFalse(xaxis.interiorIntersects(xaxis));
    assertTrue(xaxis.contains(empty));
    assertFalse(xaxis.interiorIntersects(empty));
    assertTrue(hemi.contains(tiny));
    assertTrue(
        hemi.contains(S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.radians(M_PI_4 - EPS))));
    assertFalse(
        hemi.contains(S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.radians(M_PI_4 + EPS))));
    assertTrue(concave.contains(hemi));
    assertTrue(concave.interiorIntersects(hemi.complement()));
    assertFalse(concave.contains(S2Cap.fromAxisHeight(concave.axis().neg(), 0.1)));
  }

  public void testAddEmptyCapToNonEmptyCap() {
    S2Cap nonEmptyCap = S2Cap.fromAxisAngle(S2Point.X_POS, S1Angle.degrees(10));
    assertEquals(nonEmptyCap, nonEmptyCap.addCap(S2Cap.empty()));
  }

  public void testAddNonEmptyCapToEmptyCap() {
    S2Cap nonEmptyCap = S2Cap.fromAxisAngle(S2Point.X_POS, S1Angle.degrees(10));
    assertEquals(nonEmptyCap, S2Cap.empty().addCap(nonEmptyCap));
  }

  public void testCustomEmpty() {
    // Verifies that clients can still create custom negative-height empty caps.
    S2Cap empty = S2Cap.fromAxisHeight(S2Point.X_POS, -1);
    assertExactly(-1.0, empty.height());
    assertTrue(empty.isEmpty());
  }

  public void testRectBound() {
    // Empty and full caps.
    assertTrue(S2Cap.empty().getRectBound().isEmpty());
    assertTrue(S2Cap.full().getRectBound().isFull());

    final double kDegreeEps = 1e-13;
    // Maximum allowable error for latitudes and longitudes measured in degrees. (assertDoubleNear
    // uses a fixed tolerance that is too relaxed, and assertDoubleEquals is too strict.)

    // Cap that includes the south pole.
    S2LatLngRect rect =
        S2Cap.fromAxisAngle(getLatLngPoint(-45, 57), S1Angle.degrees(50)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), -90, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 5, kDegreeEps);
    assertTrue(rect.lng().isFull());

    // Cap that is tangent to the north pole.
    rect =
        S2Cap.fromAxisAngle(new S2Point(1, 0, 1).normalize(), S1Angle.radians(M_PI_4))
            .getRectBound();
    assertExactly(rect.lat().lo(), 0);
    assertDoubleEquals(rect.lat().hi(), M_PI_2);
    assertTrue(rect.lng().isFull());

    rect =
        S2Cap.fromAxisAngle(new S2Point(1, 0, 1).normalize(), S1Angle.degrees(45)).getRectBound();
    assertDoubleNear(rect.latLo().degrees(), 0, kDegreeEps);
    assertDoubleNear(rect.latHi().degrees(), 90, kDegreeEps);
    assertTrue(rect.lng().isFull());

    // The eastern hemisphere.
    rect =
        S2Cap.fromAxisAngle(new S2Point(0, 1, 0), S1Angle.radians(M_PI_2 + 5e-16)).getRectBound();
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
    // For each cube face, we construct some cells on that face and some caps whose positions are
    // relative to that face, and then check for the expected intersection/containment results.

    // The distance from the center of a face to one of its vertices.
    final double kFaceRadius = atan(M_SQRT2);

    for (int face = 0; face < 6; ++face) {
      // The cell consisting of the entire face.
      S2Cell rootCell = S2Cell.fromFace(face);

      // A leaf cell at the midpoint of the v=1 edge.
      S2Cell edgeCell = new S2Cell(S2Projections.faceUvToXyz(face, 0, 1 - EPS));

      // A leaf cell at the u=1, v=1 corner.
      S2Cell cornerCell = new S2Cell(S2Projections.faceUvToXyz(face, 1 - EPS, 1 - EPS));

      // Quick check for full and empty caps.
      assertTrue(S2Cap.full().contains(rootCell));
      assertFalse(S2Cap.empty().mayIntersect(rootCell));

      // Check intersections with the bounding caps of the leaf cells that are adjacent to
      // 'cornerCell' along the Hilbert curve. Because this corner is at (u=1,v=1), the curve stays
      // locally within the same cube face.
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
        // A cap that barely contains all of 'capFace'.
        S2Point center = S2Projections.getNorm(capFace);
        S2Cap covering = S2Cap.fromAxisAngle(center, S1Angle.radians(kFaceRadius + EPS));
        assertEquals(covering.contains(rootCell), capFace == face);
        assertEquals(covering.mayIntersect(rootCell), capFace != antiFace);
        assertEquals(covering.contains(edgeCell), center.dotProd(edgeCell.getCenter()) > 0.1);
        assertEquals(covering.contains(edgeCell), covering.mayIntersect(edgeCell));
        assertEquals(covering.contains(cornerCell), capFace == face);
        assertEquals(covering.mayIntersect(cornerCell), center.dotProd(cornerCell.getCenter()) > 0);

        // A cap that barely intersects the edges of 'capFace'.
        S2Cap bulging = S2Cap.fromAxisAngle(center, S1Angle.radians(M_PI_4 + EPS));
        assertFalse(bulging.contains(rootCell));
        assertEquals(bulging.mayIntersect(rootCell), capFace != antiFace);
        assertEquals(bulging.contains(edgeCell), capFace == face);
        assertEquals(bulging.mayIntersect(edgeCell), center.dotProd(edgeCell.getCenter()) > 0.1);
        assertFalse(bulging.contains(cornerCell));
        assertFalse(bulging.mayIntersect(cornerCell));

        // A singleton cap.
        S2Cap singleton = S2Cap.fromAxisAngle(center, S1Angle.radians(0));
        assertEquals(singleton.mayIntersect(rootCell), capFace == face);
        assertFalse(singleton.mayIntersect(edgeCell));
        assertFalse(singleton.mayIntersect(cornerCell));
      }
    }
  }

  public void testExpanded() {
    assertTrue(S2Cap.empty().expanded(S1Angle.radians(2)).isEmpty());
    assertTrue(S2Cap.full().expanded(S1Angle.radians(2)).isFull());
    S2Cap cap50 = S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.degrees(50));
    S2Cap cap51 = S2Cap.fromAxisAngle(new S2Point(1, 0, 0), S1Angle.degrees(51));
    assertTrue(cap50.expanded(S1Angle.radians(0)).approxEquals(cap50));
    assertTrue(cap50.expanded(S1Angle.degrees(1)).approxEquals(cap51));
    assertFalse(cap50.expanded(S1Angle.degrees(129.99)).isFull());
    assertTrue(cap50.expanded(S1Angle.degrees(130.01)).isFull());
  }

  public void testGetCentroid() {
    // Empty and full caps.
    assertEquals(new S2Point(), S2Cap.empty().getCentroid());
    assertLessOrEqual(S2Cap.full().getCentroid().norm(), 1e-15);

    // Random caps.
    for (int i = 0; i < 100; ++i) {
      S2Point center = data.getRandomPoint();
      double height = data.uniform(0.0, 2.0);
      S2Cap cap = S2Cap.fromAxisHeight(center, height);
      S2Point centroid = cap.getCentroid();
      S2Point expected = center.mul((1.0 - height / 2.0) * cap.area());
      assertLessOrEqual(expected.sub(centroid).norm(), 1e-15);
    }
  }

  public void testUnion() {
    // Two caps which have the same center but one has a larger radius.
    S2Cap a = S2Cap.fromAxisAngle(getLatLngPoint(50.0, 10.0), S1Angle.degrees(0.2));
    S2Cap b = S2Cap.fromAxisAngle(getLatLngPoint(50.0, 10.0), S1Angle.degrees(0.3));
    assertTrue(b.contains(a));
    assertEquals(b, a.union(b));

    // Two caps where one is the full cap.
    assertTrue(a.union(S2Cap.full()).isFull());

    // Two caps where one is the empty cap.
    assertEquals(a, a.union(S2Cap.empty()));

    // Two caps which have different centers, one entirely encompasses the other.
    S2Cap c = S2Cap.fromAxisAngle(getLatLngPoint(51.0, 11.0), S1Angle.degrees(1.5));
    assertTrue(c.contains(a));
    assertEquals(a.union(c).axis(), c.axis());
    assertEquals(a.union(c).angle(), c.angle());

    // Two entirely disjoint caps.
    S2Cap d = S2Cap.fromAxisAngle(getLatLngPoint(51.0, 11.0), S1Angle.degrees(0.1));
    assertFalse(d.contains(a));
    assertFalse(d.intersects(a));
    assertTrue(a.union(d).approxEquals(d.union(a)));
    assertDoubleNear(50.4588, new S2LatLng(a.union(d).axis()).lat().degrees(), 0.001);
    assertDoubleNear(10.4525, new S2LatLng(a.union(d).axis()).lng().degrees(), 0.001);
    assertDoubleNear(0.7425, a.union(d).angle().degrees(), 0.001);

    // Two partially overlapping caps.
    S2Cap e = S2Cap.fromAxisAngle(getLatLngPoint(50.3, 10.3), S1Angle.degrees(0.2));
    assertFalse(e.contains(a));
    assertTrue(e.intersects(a));
    assertTrue(a.union(e).approxEquals(e.union(a)));
    assertDoubleNear(50.1500, new S2LatLng(a.union(e).axis()).lat().degrees(), 0.001);
    assertDoubleNear(10.1495, new S2LatLng(a.union(e).axis()).lng().degrees(), 0.001);
    assertDoubleNear(0.3781, a.union(e).angle().degrees(), 0.001);

    // Two very large caps, whose radius sums to in excess of 180 degrees, and whose centers are not
    // antipodal.
    S2Cap f = S2Cap.fromAxisAngle(new S2Point(0, 0, 1).normalize(), S1Angle.degrees(150));
    S2Cap g = S2Cap.fromAxisAngle(new S2Point(0, 1, 0).normalize(), S1Angle.degrees(150));
    assertTrue(f.union(g).isFull());

    // Two non-overlapping hemisphere caps with antipodal centers.
    S2Cap hemi = S2Cap.fromAxisHeight(new S2Point(0, 0, 1).normalize(), 1);
    assertTrue(hemi.union(hemi.complement()).isFull());
  }

  public void testSerialization() throws IOException {
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    S2Cap testCap = S2Cap.fromAxisHeight(S2Point.X_NEG, 0.123);
    testCap.encode(bos);
    assertEquals(testCap, S2Cap.decode(new ByteArrayInputStream(bos.toByteArray())));
  }

  public void testDecodeEmptyCap() throws IOException {
    checkCoder("000000000000F03F00000000000000000000000000000000000000000000F0BF", S2Cap.empty());
  }

  public void testDecodeFullCap() throws IOException {
    checkCoder("000000000000F03F000000000000000000000000000000000000000000001040", S2Cap.full());
  }

  public void testDecodeCapWithHeight() throws IOException {
    checkCoder(
        "00000000000000000000000000000000000000000000F03F0000000000001040",
        S2Cap.fromAxisHeight(new S2Point(0, 0, 1), 5.0));
  }

  private static void checkCoder(String hex, S2Cap expected) throws IOException {
    BaseEncoding testFormat = BaseEncoding.base16();
    ByteArrayInputStream bais = new ByteArrayInputStream(testFormat.decode(hex));
    assertEquals(expected, S2Cap.decode(bais));
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    expected.encode(baos);
    assertEquals(hex, testFormat.encode(baos.toByteArray()));
  }

  public void testCoder() throws IOException {
    S2Cap cap = S2Cap.fromAxisChord(S2Point.Y_POS, S1ChordAngle.RIGHT);
    assertEquals(cap, roundtrip(S2Cap.CODER, cap));
  }
}
