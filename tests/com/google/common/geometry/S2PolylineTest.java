/*
 * Copyright 2006 Google Inc.
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

import com.google.common.collect.Lists;

import junit.framework.Assert;

import java.util.List;

/**
 * Tests for {@link S2Polyline}.
 *
 */
public strictfp class S2PolylineTest extends GeometryTestCase {

  @Override
  public void setUp() {
    super.setUp();
  }

  public void testBasic() {
    List<S2Point> vertices = Lists.newArrayList();
    S2Polyline empty = new S2Polyline(vertices);
    assertEquals(empty.getRectBound(), S2LatLngRect.empty());
  }

  public void testGetLengthCentroid() {
    // Construct random great circles and divide them randomly into segments.
    // Then make sure that the length and centroid are correct. Note that
    // because of the way the centroid is computed, it does not matter how
    // we split the great circle into segments.

    for (int i = 0; i < 100; ++i) {
      // Choose a coordinate frame for the great circle.
      S2Point x = randomPoint();
      S2Point y = S2Point.normalize(S2Point.crossProd(x, randomPoint()));
      S2Point z = S2Point.normalize(S2Point.crossProd(x, y));

      List<S2Point> vertices = Lists.newArrayList();
      for (double theta = 0; theta < 2 * S2.M_PI; theta += Math.pow(rand.nextDouble(), 10)) {
        S2Point p = S2Point.add(S2Point.mul(x, Math.cos(theta)), S2Point.mul(y, Math.sin(theta)));
        if (vertices.isEmpty() || !p.equals(vertices.get(vertices.size() - 1))) {
          vertices.add(p);
        }
      }
      // Close the circle.
      vertices.add(vertices.get(0));
      S2Polyline line = new S2Polyline(vertices);
      S1Angle length = line.getArclengthAngle();
      assertTrue(Math.abs(length.radians() - 2 * S2.M_PI) < 2e-14);
    }
  }

  public void testMayIntersect() {
    List<S2Point> vertices = Lists.newArrayList();
    vertices.add(S2Point.normalize(new S2Point(1, -1.1, 0.8)));
    vertices.add(S2Point.normalize(new S2Point(1, -0.8, 1.1)));
    S2Polyline line = new S2Polyline(vertices);
    for (int face = 0; face < 6; ++face) {
      S2Cell cell = S2Cell.fromFacePosLevel(face, (byte) 0, 0);
      assertEquals(line.mayIntersect(cell), (face & 1) == 0);
    }
  }

  public void testInterpolate() {
    List<S2Point> vertices = Lists.newArrayList();
    vertices.add(new S2Point(1, 0, 0));
    vertices.add(new S2Point(0, 1, 0));
    vertices.add(S2Point.normalize(new S2Point(0, 1, 1)));
    vertices.add(new S2Point(0, 0, 1));
    S2Polyline line = new S2Polyline(vertices);

    assertEquals(line.interpolate(-0.1), vertices.get(0));
    assertTrue(S2.approxEquals(
        line.interpolate(0.1), S2Point.normalize(new S2Point(1, Math.tan(0.2 * S2.M_PI / 2), 0))));
    assertTrue(S2.approxEquals(line.interpolate(0.25), S2Point.normalize(new S2Point(1, 1, 0))));

    assertEquals(line.interpolate(0.5), vertices.get(1));
    assertEquals(line.interpolate(0.75), vertices.get(2));
    assertEquals(line.interpolate(1.1), vertices.get(3));
  }

  public void testEqualsAndHashCode() {
    List<S2Point> vertices = Lists.newArrayList();
    vertices.add(new S2Point(1, 0, 0));
    vertices.add(new S2Point(0, 1, 0));
    vertices.add(S2Point.normalize(new S2Point(0, 1, 1)));
    vertices.add(new S2Point(0, 0, 1));


    S2Polyline line1 = new S2Polyline(vertices);
    S2Polyline line2 = new S2Polyline(vertices);

    checkEqualsAndHashCodeMethods(line1, line2, true);

    List<S2Point> moreVertices = Lists.newLinkedList(vertices);
    moreVertices.remove(0);

    S2Polyline line3 = new S2Polyline(moreVertices);

    checkEqualsAndHashCodeMethods(line1, line3, false);
    checkEqualsAndHashCodeMethods(line1, null, false);
    checkEqualsAndHashCodeMethods(line1, "", false);
  }

  public void testProject() {
    List<S2Point> latLngs = Lists.newArrayList();
    latLngs.add(S2LatLng.fromDegrees(0, 0).toPoint());
    latLngs.add(S2LatLng.fromDegrees(0, 1).toPoint());
    latLngs.add(S2LatLng.fromDegrees(0, 2).toPoint());
    latLngs.add(S2LatLng.fromDegrees(1, 2).toPoint());
    S2Polyline line = new S2Polyline(latLngs);

    int edgeIndex = -1;
    S2Point testPoint = null;

    testPoint = S2LatLng.fromDegrees(0.5, -0.5).toPoint();
    edgeIndex = line.getNearestEdgeIndex(testPoint);
    assertTrue(S2.approxEquals(
        line.projectToEdge(testPoint, edgeIndex), S2LatLng.fromDegrees(0, 0).toPoint()));
    assertEquals(0, edgeIndex);

    testPoint = S2LatLng.fromDegrees(0.5, 0.5).toPoint();
    edgeIndex = line.getNearestEdgeIndex(testPoint);
    assertTrue(S2.approxEquals(
        line.projectToEdge(testPoint, edgeIndex), S2LatLng.fromDegrees(0, 0.5).toPoint()));
    assertEquals(0, edgeIndex);

    testPoint = S2LatLng.fromDegrees(0.5, 1).toPoint();
    edgeIndex = line.getNearestEdgeIndex(testPoint);
    assertTrue(S2.approxEquals(
        line.projectToEdge(testPoint, edgeIndex), S2LatLng.fromDegrees(0, 1).toPoint()));
    assertEquals(0, edgeIndex);

    testPoint = S2LatLng.fromDegrees(-0.5, 2.5).toPoint();
    edgeIndex = line.getNearestEdgeIndex(testPoint);
    assertTrue(S2.approxEquals(
        line.projectToEdge(testPoint, edgeIndex), S2LatLng.fromDegrees(0, 2).toPoint()));
    assertEquals(1, edgeIndex);

    testPoint = S2LatLng.fromDegrees(2, 2).toPoint();
    edgeIndex = line.getNearestEdgeIndex(testPoint);
    assertTrue(S2.approxEquals(
        line.projectToEdge(testPoint, edgeIndex), S2LatLng.fromDegrees(1, 2).toPoint()));
    assertEquals(2, edgeIndex);
  }

  /**
   * Utility for testing equals() and hashCode() results at once.
   * Tests that lhs.equals(rhs) matches expectedResult, as well as
   * rhs.equals(lhs).  Also tests that hashCode() return values are
   * equal if expectedResult is true.  (hashCode() is not tested if
   * expectedResult is false, as unequal objects can have equal hashCodes.)
   *
   * @param lhs An Object for which equals() and hashCode() are to be tested.
   * @param rhs As lhs.
   * @param expectedResult True if the objects should compare equal,
   *   false if not.
   */
  private static void checkEqualsAndHashCodeMethods(Object lhs, Object rhs,
                                             boolean expectedResult) {
    if ((lhs == null) && (rhs == null)) {
      Assert.assertTrue(
          "Your check is dubious...why would you expect null != null?",
          expectedResult);
      return;
    }

    if ((lhs == null) || (rhs == null)) {
      Assert.assertFalse(
          "Your check is dubious...why would you expect an object "
          + "to be equal to null?", expectedResult);
    }

    if (lhs != null) {
      assertEquals(expectedResult, lhs.equals(rhs));
    }
    if (rhs != null) {
      assertEquals(expectedResult, rhs.equals(lhs));
    }

    if (expectedResult) {
      String hashMessage =
          "hashCode() values for equal objects should be the same";
      Assert.assertTrue(hashMessage, lhs.hashCode() == rhs.hashCode());
    }
  }
}
