/*
 * Copyright 2014 Google Inc.
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

import static com.google.common.geometry.S1ChordAngle.add;
import static com.google.common.geometry.S1ChordAngle.sub;
import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.M_PI_2;
import static java.lang.Math.PI;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.tan;


/** Tests for {@link S1ChordAngle}. */
public strictfp class S1ChordAngleTest extends GeometryTestCase {

  public void testTwoPointConstructor() {
    for (int iter = 0; iter < 100; ++iter) {
      Matrix frame = data.getRandomFrame();
      S2Point x = frame.getCol(0);
      S2Point y = frame.getCol(1);
      S2Point z = frame.getCol(2);
      assertEquals(S1Angle.ZERO, new S1ChordAngle(z, z).toAngle());
      assertDoubleNear(PI, new S1ChordAngle(z.neg(), z).toAngle().radians(), 1e-7);
      assertDoubleEquals(PI / 2, new S1ChordAngle(x, z).toAngle().radians());
      S2Point w = y.add(z).normalize();
      assertDoubleEquals(PI / 4, new S1ChordAngle(w, z).toAngle().radians());
    }
  }

  public void testFromLength2() {
    assertExactly(0.0, S1ChordAngle.fromLength2(0).degrees());
    assertDoubleEquals(60.0, S1ChordAngle.fromLength2(1).degrees());
    assertDoubleEquals(90.0, S1ChordAngle.fromLength2(2).degrees());
    assertExactly(180.0, S1ChordAngle.fromLength2(4).degrees());
    assertExactly(180.0, S1ChordAngle.fromLength2(5).degrees());
  }

  public void testZero() {
    assertEquals(S1Angle.ZERO, S1ChordAngle.ZERO.toAngle());
  }

  public void testStraight() {
    assertEquals(S1Angle.degrees(180), S1ChordAngle.STRAIGHT.toAngle());
  }

  public void testRight() {
    assertEquals(90.0, S1ChordAngle.RIGHT.degrees(), Platform.ulp(90D));
  }

  public void testInfinity() {
    assertTrue(S1ChordAngle.STRAIGHT.compareTo(S1ChordAngle.INFINITY) < 0);
    assertEquals(S1Angle.INFINITY, S1ChordAngle.INFINITY.toAngle());
  }

  public void testNegative() {
    assertTrue(S1ChordAngle.NEGATIVE.compareTo(S1ChordAngle.ZERO) < 0);
    assertTrue(S1ChordAngle.NEGATIVE.toAngle().compareTo(S1Angle.ZERO) < 0);
  }

  public void testEquals() {
    S1ChordAngle[] angles = {
      S1ChordAngle.NEGATIVE, S1ChordAngle.ZERO, S1ChordAngle.STRAIGHT, S1ChordAngle.INFINITY
    };
    for (int i = 0; i < angles.length; ++i) {
      for (int j = 0; j < angles.length; ++j) {
        assertEquals(angles[i] == angles[j], angles[i].equals(angles[j]));
      }
    }
  }

  public void testPredicates() {
    assertTrue(S1ChordAngle.ZERO.isZero());
    assertFalse(S1ChordAngle.ZERO.isNegative());
    assertFalse(S1ChordAngle.ZERO.isSpecial());
    assertFalse(S1ChordAngle.STRAIGHT.isSpecial());
    assertTrue(S1ChordAngle.NEGATIVE.isNegative());
    assertTrue(S1ChordAngle.NEGATIVE.isSpecial());
    assertTrue(S1ChordAngle.INFINITY.isInfinity());
    assertTrue(S1ChordAngle.INFINITY.isSpecial());
  }

  public void testToFromS1Angle() {
    assertExactly(0.0, S1ChordAngle.fromS1Angle(S1Angle.ZERO).toAngle().radians());
    assertExactly(4.0, S1ChordAngle.fromS1Angle(S1Angle.radians(Math.PI)).getLength2());
    assertExactly(PI, S1ChordAngle.fromS1Angle(S1Angle.radians(Math.PI)).toAngle().radians());
    assertEquals(S1Angle.INFINITY, S1ChordAngle.fromS1Angle(S1Angle.INFINITY).toAngle());
    assertEquals(
        S1Angle.INFINITY,
        S1ChordAngle.fromS1Angle(S1Angle.radians(Double.POSITIVE_INFINITY)).toAngle());
    assertTrue(S1ChordAngle.fromS1Angle(S1Angle.radians(-1)).toAngle().radians() < 0.0);
    assertDoubleEquals(1.0, S1ChordAngle.fromS1Angle(S1Angle.radians(1.0)).toAngle().radians());
  }

  public void testSuccessor() {
    assertEquals(S1ChordAngle.ZERO, S1ChordAngle.NEGATIVE.successor());
    assertEquals(S1ChordAngle.INFINITY, S1ChordAngle.STRAIGHT.successor());
    assertEquals(S1ChordAngle.INFINITY, S1ChordAngle.INFINITY.successor());
    S1ChordAngle x = S1ChordAngle.NEGATIVE;
    for (int i = 0; i < 10; ++i) {
      assertTrue(x.compareTo(x.successor()) < 0);
      x = x.successor();
    }
  }

  public void testPredecessor() {
    assertEquals(S1ChordAngle.STRAIGHT, S1ChordAngle.INFINITY.predecessor());
    assertEquals(S1ChordAngle.NEGATIVE, S1ChordAngle.ZERO.predecessor());
    assertEquals(S1ChordAngle.NEGATIVE, S1ChordAngle.NEGATIVE.predecessor());
    S1ChordAngle x = S1ChordAngle.INFINITY;
    for (int i = 0; i < 10; ++i) {
      assertTrue(x.compareTo(x.predecessor()) > 0);
      x = x.predecessor();
    }
  }

  public void testArithmetic() {
    S1ChordAngle zero = S1ChordAngle.ZERO;
    S1ChordAngle degree30 = S1ChordAngle.fromS1Angle(S1Angle.degrees(30));
    S1ChordAngle degree60 = S1ChordAngle.fromS1Angle(S1Angle.degrees(60));
    S1ChordAngle degree90 = S1ChordAngle.fromS1Angle(S1Angle.degrees(90));
    S1ChordAngle degree120 = S1ChordAngle.fromS1Angle(S1Angle.degrees(120));
    S1ChordAngle degree180 = S1ChordAngle.STRAIGHT;
    assertExactly(0.0, add(zero, zero).degrees());
    assertExactly(0.0, sub(zero, zero).degrees());
    assertExactly(0.0, sub(degree60, degree60).degrees());
    assertExactly(0.0, sub(degree180, degree180).degrees());
    assertExactly(0.0, sub(zero, degree60).degrees());
    assertExactly(0.0, sub(degree30, degree90).degrees());
    assertDoubleEquals(60.0, add(degree60, zero).degrees());
    assertDoubleEquals(60.0, sub(degree60, zero).degrees());
    assertDoubleEquals(60.0, add(zero, degree60).degrees());
    assertDoubleEquals(90.0, add(degree30, degree60).degrees());
    assertDoubleEquals(90.0, add(degree60, degree30).degrees());
    assertDoubleEquals(60.0, sub(degree90, degree30).degrees());
    assertDoubleEquals(30.0, sub(degree90, degree60).degrees());
    assertExactly(180.0, add(degree180, zero).degrees());
    assertExactly(180.0, sub(degree180, zero).degrees());
    assertExactly(180.0, add(degree90, degree90).degrees());
    assertExactly(180.0, add(degree120, degree90).degrees());
    assertExactly(180.0, add(degree120, degree120).degrees());
    assertExactly(180.0, add(degree30, degree180).degrees());
    assertExactly(180.0, add(degree180, degree180).degrees());
  }

  // Verifies that S1ChordAngle is capable of adding and subtracting angles extremely accurately up
  // to Pi/2 radians. (Accuracy continues to be good well beyond this value but degrades as angles
  // approach Pi.)
  public void testArithmeticPrecision() {
    S1ChordAngle kEps = S1ChordAngle.fromRadians(1e-15);
    S1ChordAngle k90 = S1ChordAngle.RIGHT;
    S1ChordAngle k90MinusEps = sub(k90, kEps);
    S1ChordAngle k90PlusEps = add(k90, kEps);

    double kMaxError = 2 * DBL_EPSILON;
    assertDoubleNear(k90MinusEps.radians(), M_PI_2 - kEps.radians(), kMaxError);
    assertDoubleNear(k90PlusEps.radians(), M_PI_2 + kEps.radians(), kMaxError);
    assertDoubleNear(sub(k90, k90MinusEps).getLength2(), kEps.getLength2(), kMaxError);
    assertDoubleNear(sub(k90, k90MinusEps).radians(), kEps.radians(), kMaxError);
    assertDoubleNear(sub(k90PlusEps, k90).radians(), kEps.radians(), kMaxError);
    assertDoubleNear(add(k90MinusEps, kEps).radians(), M_PI_2, kMaxError);
  }

  public void testTrigonometry() {
    final int iters = 20;
    for (int iter = 0; iter <= iters; ++iter) {
      double radians = PI * iter / iters;
      S1ChordAngle angle = S1ChordAngle.fromS1Angle(S1Angle.radians(radians));
      assertEquals(sin(radians), S1ChordAngle.sin(angle), 1e-15);
      assertEquals(cos(radians), S1ChordAngle.cos(angle), 1e-15);
      // Since the tan(x) is unbounded near Pi/4, we map the result back to an angle before
      // comparing. (The assertion is that the result is equal to the tangent of a nearby angle.)
      assertEquals(atan(tan(radians)), atan(S1ChordAngle.tan(angle)), 1e-15);
    }

    // Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
    S1ChordAngle angle90 = S1ChordAngle.fromLength2(2);
    S1ChordAngle angle180 = S1ChordAngle.fromLength2(4);
    assertExactly(1.0, S1ChordAngle.sin(angle90));
    assertExactly(0.0, S1ChordAngle.cos(angle90));
    assertExactly(Double.POSITIVE_INFINITY, S1ChordAngle.tan(angle90));
    assertExactly(0.0, S1ChordAngle.sin(angle180));
    assertExactly(-1.0, S1ChordAngle.cos(angle180));
    assertExactly(-0.0, S1ChordAngle.tan(angle180));
  }

  public void testPlusError() {
    assertEquals(S1ChordAngle.NEGATIVE, S1ChordAngle.NEGATIVE.plusError(5));
    assertEquals(S1ChordAngle.INFINITY, S1ChordAngle.INFINITY.plusError(-5));
    assertEquals(S1ChordAngle.STRAIGHT, S1ChordAngle.STRAIGHT.plusError(5));
    assertEquals(S1ChordAngle.ZERO, S1ChordAngle.ZERO.plusError(-5));
    assertEquals(S1ChordAngle.fromLength2(1.25), S1ChordAngle.fromLength2(1).plusError(0.25));
    assertEquals(S1ChordAngle.fromLength2(0.75), S1ChordAngle.fromLength2(1).plusError(-0.25));
  }

  /**
   * Verifies that the error bound returned by {@link S1ChordAngle#getS2PointConstructorMaxError} is
   * large enough.
   */
  public void testGetS2PointConstructorMaxError() {
    for (int iter = 0; iter < 10000; ++iter) {
      data.setSeed(iter);
      S2Point x = data.getRandomPoint();
      S2Point y = data.getRandomPoint();
      if (data.oneIn(10)) {
        // Occasionally test a point pair that is nearly identical or antipodal.
        S1Angle r = S1Angle.radians(1e-15 * data.nextDouble());
        y = S2EdgeUtil.interpolateAtDistance(r, x, y);
        if (data.oneIn(2)) {
          y = y.neg();
        }
      }
      S1ChordAngle dist = new S1ChordAngle(x, y);
      double error = dist.getS2PointConstructorMaxError();
      String msg = "angle=" + dist + ", iter=" + iter;
      assertTrue(msg, S2Predicates.compareDistance(x, y, dist.plusError(error).getLength2()) <= 0);
      assertTrue(msg, S2Predicates.compareDistance(x, y, dist.plusError(-error).getLength2()) >= 0);
    }
  }

  public void testHashCodeZero() {
    // Check that hashCode() and equals(...) work consistently for the ±0 edge case.
    S1ChordAngle positive0 = S1ChordAngle.fromLength2(0);
    S1ChordAngle negative0 = S1ChordAngle.fromLength2(-0.0);

    assertTrue(positive0.equals(negative0));
    assertEquals(positive0.hashCode(), negative0.hashCode());
  }

  public void testHashCodeDifferent() {
    S1ChordAngle zero = S1ChordAngle.fromLength2(0);
    S1ChordAngle nonZero = S1ChordAngle.fromLength2(1);

    assertFalse(zero.equals(nonZero));
    // ant_test uses a different (old) version of JUnit without assertNotEquals.
    assertTrue(zero.hashCode() != nonZero.hashCode());
  }

  public void testS1AngleConsistency() {
    // This test checks that the error bounds in the S1ChordAngle constructors are consistent with
    // the maximum error in S1Angle(x, y).
    double maxS1AngleError = 3.25 * DBL_EPSILON;
    for (int iter = 0; iter < 10000; ++iter) {
      S2Point x = data.getRandomPoint();
      S2Point y = data.getRandomPoint();
      S1ChordAngle dist1 = S1ChordAngle.fromS1Angle(new S1Angle(x, y));
      S1ChordAngle dist2 = new S1ChordAngle(x, y);
      double maxError =
          (maxS1AngleError
              + dist1.getS1AngleConstructorMaxError()
              + dist2.getS2PointConstructorMaxError());
      assertTrue(dist1.compareTo(dist2.plusError(maxError)) <= 0);
      assertTrue(dist1.compareTo(dist2.plusError(-maxError)) >= 0);
    }
  }
}
