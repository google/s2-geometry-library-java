/*
 * Copyright 2013 Google Inc.
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
import static com.google.common.geometry.S2TextFormat.makePoint;
import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Verifies S2Point. */
@RunWith(JUnit4.class)
public class S2PointTest extends GeometryTestCase {
  @Test
  public void testAddition() {
    S2Point a = new S2Point(1, 2, 3);
    S2Point b = new S2Point(1, 1, 1);
    S2Point aPlusB = a.add(b);
    assertEquals(new S2Point(2, 3, 4), aPlusB);
  }

  @Test
  public void testSubtraction() {
    S2Point a = new S2Point(1, 2, 3);
    S2Point b = new S2Point(1, 1, 1);
    S2Point aMinusB = a.sub(b);
    assertEquals(new S2Point(0, 1, 2), aMinusB);
  }

  @Test
  public void testScalarMultiplication() {
    S2Point a = new S2Point(1, 2, 3);
    assertEquals(new S2Point(5, 10, 15), a.mul(5.0));
  }

  @Test
  public void testScalarDivision() {
    S2Point a = new S2Point(3, 6, 9);
    assertEquals(new S2Point(1, 2, 3), a.div(3));
  }

  @Test
  public void testNegation() {
    S2Point a = new S2Point(3, 6, 9);
    assertEquals(new S2Point(-3, -6, -9), a.neg());
  }

  @Test
  public void testComponentWiseAbs() {
    S2Point a = new S2Point(-3, 6, -9);
    assertEquals(new S2Point(3, 6, 9), a.fabs());
  }

  @Test
  public void testVectorMethods() {
    S2Point a = new S2Point(1, 2, 3);
    S2Point b = new S2Point(0, 1, 0);
    assertEquals(new S2Point(1, 3, 3), a.add(b));
    assertEquals(new S2Point(1, 1, 3), a.sub(b));
    assertEquals(new S2Point(2, 4, 6), a.mul(2));
    assertEquals(new S2Point(0.5, 1, 1.5), a.div(2));
    assertEquals(2.0, a.dotProd(b), 0);
    assertEquals(new S2Point(-3, 0, 1), a.crossProd(b));
    assertEquals(b, b.normalize());
    assertEquals(a, a.neg().fabs());
    assertEquals(S2Point.Y_NEG, S2Point.X_POS.ortho());
  }

  @Test
  public void testRotateSimple() {
    // Check simple axial cases.
    S2Point p = makePoint("0:0");
    assertPointsWithinDistance(makePoint("0:90"), p.rotate(S2Point.Z_POS, M_PI_2), 1e-15);
    assertPointsWithinDistance(makePoint("-90:0"), p.rotate(S2Point.Y_POS, M_PI_2), 1e-15);
    assertPointsWithinDistance(makePoint("0:0"), p.rotate(S2Point.X_POS, M_PI_2), 1e-15);
    // Verify rotations in the plane containing the origin and two random points.
    for (int i = 0; i < 1000; i++) {
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      S2Point axis = a.crossProd(b);
      double angle = a.angle(b);
      // Rotate 'a' onto 'b'.
      assertPointsWithinDistance(b, a.rotate(axis, angle), 1e-15);
      // Rotate 'b' onto 'a'.
      assertPointsWithinDistance(a, b.rotate(axis, -angle), 1e-15);
      // Rotate 'a' to the midpoint of 'a' and 'b'.
      // TODO(user): This test is fragile.
      assertPointsWithinDistance(a.add(b).normalize(), a.rotate(axis, angle / 2), 1e-14);
      // Rotate 'b' to the antipodal point of 'a'.
      assertPointsWithinDistance(a.neg(), b.rotate(axis, PI - angle), 1e-14);
    }
  }

  @Test
  public void testRotateDirection() {
    // A point at the corner of Oregon, California, and Nevada that we will rotate around.
    S2Point axis = makePoint("42:-120");
    // A point one degree east of 'axis' that will be rotated clockwise and counter-clockwise.
    S2Point p = makePoint("42:-119");

    // Verify that rotation by a positive radians is counter-clockwise. Starting a degree east of
    // 'axis', rotating +90 degrees or Pi/2 radians should end up approximately a degree north.
    assertPointsWithinDegrees(makePoint("43:-120"), p.rotate(axis, PI / 2), 0.26);

    // Verify that rotation by a negative radians is clockwise. Starting a degree east of 'axis',
    // rotating -90 degrees or -Pi/2 radians should end up approximately a degree south of 'axis'.
    assertPointsWithinDegrees(makePoint("41:-120.0"), p.rotate(axis, -PI / 2), 0.26);
  }

  private static void checkRotate(S2Point p, S2Point axis, double angle) {
    S2Point result = S2Point.rotate(p, axis, angle);

    // "result" should be unit length.
    assertTrue(S2.isUnitLength(result));

    // "result" and "p" should be the same distance from "axis".
    double kMaxPositionError = 1e-15;
    assertLessOrEqual(
        new S1Angle(result, axis).sub(new S1Angle(p, axis)).abs().radians(), kMaxPositionError);

    // Check that the rotation angle is correct. We allow a fixed error in the *position* of the
    // result, so we need to convert this into a rotation angle. The allowable error can be very
    // large as "p" approaches "axis".
    double axisDistance = p.crossProd(axis).norm();
    double maxRotationError = (axisDistance < kMaxPositionError)
        ? 2 * PI
        : asin(kMaxPositionError / axisDistance);

    double actualRotation = S2.turnAngle(p, axis, result) + PI;
    double rotationError = Platform.IEEEremainder(angle - actualRotation, 2 * PI);
    assertLessOrEqual(rotationError, maxRotationError);
  }

  @Test
  public void testRotate() {
    for (int iter = 0; iter < 1000; ++iter) {
      S2Point axis = data.getRandomPoint();
      S2Point target = data.getRandomPoint();
      // Choose a distance whose logarithm is uniformly distributed.
      double distance = PI * Math.pow(1e-15, data.nextDouble());
      // Sometimes choose points near the far side of the axis.
      if (data.oneIn(5)) {
        distance = PI - distance;
      }
      S2Point p = S2EdgeUtil.getPointOnLine(axis, target, S1Angle.radians(distance));
      // Choose the rotation angle.
      double angle = 2 * PI * Math.pow(1e-15, data.nextDouble());
      if (data.oneIn(3)) {
        angle = -angle;
      }
      if (data.oneIn(10)) {
        angle = 0;
      }
      checkRotate(p, axis, angle);
    }
  }

  @Test
  public void testScalarTripleProduct() {
    for (int i = 0; i < 1000; i++) {
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      S2Point c = data.getRandomPoint();
      assertEquals(a.dotProd(b.crossProd(c)), S2Point.scalarTripleProduct(a, b, c), 1e-14);
    }
  }

  @Test
  public void testCrossProdNorm() {
    for (int i = 0; i < 1000; i++) {
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      assertEquals(a.crossProd(b).norm(), a.crossProdNorm(b), 1e-14);
    }
  }

  @Test
  public void testS2Region() {
    S2Point point = new S2Point(1, 0, 0);

    S2Cap expectedCapBound = S2Cap.fromAxisHeight(point, 0);
    assertEquals(expectedCapBound, point.getCapBound());

    S2LatLng ll = new S2LatLng(point);
    S2LatLngRect expectedRect = new S2LatLngRect(ll, ll);
    assertEquals(expectedRect, point.getRectBound());

    // The leaf cell containing a point is still much larger than the point.
    S2Cell cell = new S2Cell(point);
    assertFalse(point.contains(cell));
    assertTrue(point.mayIntersect(cell));
  }

  @Test
  public void testShape() {
    List<S2Point> points = ImmutableList.of(makePoint("0:0"), makePoint("1:1"), makePoint("2:2"));
    for (int j = 0; j < points.size(); j++) {
      List<S2Point> subset = points.subList(0, j);
      S2Shape shape = S2Point.Shape.fromList(subset);
      assertFalse(shape.containsOrigin());
      assertFalse(shape.hasInterior());
      assertEquals(subset.size(), shape.numEdges());
      assertEquals(subset.size(), shape.numChains());

      if (j == 0) {
        assertTrue(shape.isEmpty());
        assertFalse(shape.isFull());
      } else {
        assertFalse(shape.isEmpty());
        assertFalse(shape.isFull());
      }

      MutableEdge edge = new MutableEdge();
      for (int i = 0; i < subset.size(); i++) {
        shape.getEdge(i, edge);
        assertEquals(edge.a, edge.b);
        assertEquals(subset.get(i), edge.a);
        assertEquals(1, shape.getChainLength(i));
        assertEquals(edge.a, shape.getChainVertex(i, 0));
        assertEquals(edge.b, shape.getChainVertex(i, 1));
      }
    }
  }

  @Test
  public void testSerialization() throws IOException {
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    S2Point.X_NEG.encode(bos);
    assertEquals(S2Point.X_NEG, S2Point.decode(new ByteArrayInputStream(bos.toByteArray())));
  }

  @Test
  public void testCoder() throws IOException {
    assertEquals(S2Point.Y_NEG, roundtrip(S2Point.CODER, S2Point.Y_NEG));
  }

  @Test
  public void testShapeCoderFast() {
    S2Point.Shape shape = S2Point.Shape.singleton(S2Point.Z_POS);
    assertShapesEqual(shape, roundtrip(S2Point.Shape.FAST_CODER, shape));
  }

  @Test
  public void testShapeCoderCompact() {
    S2Point.Shape shape = S2Point.Shape.singleton(S2Point.Z_POS);
    assertShapesEqual(shape, roundtrip(S2Point.Shape.COMPACT_CODER, shape));
  }

  @Test
  public void testBuilder_origin() {
    assertEquals(S2Point.ORIGIN, new S2Point.Builder().build());
  }

  @Test
  public void testBuilder_single() {
    S2Point base = new S2Point(0.1, 0.2, 0.3);
    S2Point built = new S2Point.Builder().add(base).build();
    assertEquals(base, built);
  }

  @Test
  public void testBuilder_multiple() {
    assertEquals(
        new S2Point(0.51, -2.23, 2.73),
        new S2Point.Builder()
            .add(new S2Point(1, -2, 3))
            .add(new S2Point(-0.5, -0.25, -0.3))
            .add(new S2Point(0.01, 0.02, 0.03))
            .build());
  }

  @Test
  public void testToBuilder() {
    S2Point base = new S2Point(0.25, 0.3, 0.5);
    S2Point.Builder builder = base.toBuilder();
    builder.add(new S2Point(0.1, 0.25, -0.5));
    assertEquals(new S2Point(0.35, 0.55, 0.0), builder.build());
    assertEquals(new S2Point(0.25, 0.3, 0.5), base); // Builder has a copy, not original
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2PointSerialization() {
    doSerializationTest(new S2Point(0.1, 0.2, 0.3));
  }
}
