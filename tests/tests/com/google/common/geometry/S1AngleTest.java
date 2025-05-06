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

import static com.google.common.geometry.Platform.nextAfter;
import static java.lang.Math.PI;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.annotations.GwtIncompatible;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S1Angle}. */
@RunWith(JUnit4.class)
public class S1AngleTest extends GeometryTestCase {

  @Test
  public void testBasic() {
    // Check that the conversion between Pi radians and 180 degrees is exact.
    assertExactly(PI, S1Angle.radians(PI).radians());
    assertExactly(180.0, S1Angle.radians(PI).degrees());
    assertExactly(PI, S1Angle.degrees(180).radians());
    assertExactly(180.0, S1Angle.degrees(180).degrees());

    assertExactly(90.0, S1Angle.radians(PI / 2).degrees());

    // Check negative angles.
    assertExactly(-90.0, S1Angle.radians(-PI / 2).degrees());
    assertExactly(-PI / 4, S1Angle.degrees(-45).radians());

    // Check that E5/E6/E7 representations work as expected.
    assertEquals(S1Angle.e5(2000000), S1Angle.degrees(20));
    assertEquals(S1Angle.e6(-60000000), S1Angle.degrees(-60));
    assertEquals(S1Angle.e7(750000000), S1Angle.degrees(75));
    assertEquals(S1Angle.e5(2000000), S1Angle.degrees(20));
    assertEquals(S1Angle.e6(-60000000), S1Angle.degrees(-60));
    assertEquals(S1Angle.e7(750000000), S1Angle.degrees(75));
    assertEquals(1234567, S1Angle.degrees(12.34567).e5());
    assertEquals(12345678, S1Angle.degrees(12.345678).e6());
    assertEquals(-123456789, S1Angle.degrees(-12.3456789).e7());
  }

  @Test
  public void testMath() {
    // 29.999999999999996
    assertAlmostEquals(30.0, S1Angle.degrees(10).add(S1Angle.degrees(20)).degrees());
    assertExactly(-10, S1Angle.degrees(10).sub(S1Angle.degrees(20)).degrees());
    assertExactly(20, S1Angle.degrees(10).mul(2.0).degrees());
    assertExactly(5, S1Angle.degrees(10).div(2.0).degrees());
    assertExactly(1.0, S1Angle.degrees(0).cos());
    assertExactly(1.0, S1Angle.degrees(90).sin());
    assertAlmostEquals(1.0, S1Angle.degrees(45).tan()); // 0.9999999999999999
  }

  @SuppressWarnings("deprecation") // earthDistance and fromEarthDistance are deprecated but tested.
  @Test
  public void testDistance() {
    // Check distance accessor for arbitrary sphere
    assertExactly(100.0 * PI, S1Angle.radians(PI).distance(100.0));
    assertExactly(50.0 * PI, S1Angle.radians(PI / 2).distance(100.0));
    assertExactly(25.0 * PI, S1Angle.radians(PI / 4).distance(100.0));
  }

  @Test
  public void testE7Overflow() {
    // Normalized angles should never overflow.
    assertEquals(-1800000000, S1Angle.degrees(-180.0).e7());
    assertEquals(1800000000, S1Angle.degrees(180.0).e7());

    // Overflow starts at 214.7483648 degrees.
    // Don't use assertThrows since GWT doesn't understand that.
    try {
      S1Angle.degrees(215).e7();
      fail("IllegalArgumentException expected, not thrown.");
    } catch (IllegalArgumentException unused) {
      // Good, test passed.
    }
    try {
      S1Angle.degrees(-215).e7();
      fail("IllegalArgumentException expected, not thrown.");
    } catch (IllegalArgumentException unused) {
      // Good, test passed.
    }
  }

  @Test
  public void testNormalize() {
    assertExactly(0.0, S1Angle.degrees(360.0).normalize().degrees());
    assertExactly(-90.0, S1Angle.degrees(-90.0).normalize().degrees());
    assertExactly(180.0, S1Angle.degrees(-180.0).normalize().degrees());
    assertExactly(90.0, S1Angle.degrees(90.0).normalize().degrees());
    assertExactly(180.0, S1Angle.degrees(180.0).normalize().degrees());
    assertExactly(-90.0, S1Angle.degrees(270.0).normalize().degrees());
    assertExactly(180.0, S1Angle.degrees(540.0).normalize().degrees());
    assertExactly(90.0, S1Angle.degrees(-270.0).normalize().degrees());

    // PI is unchanged.
    assertExactly(PI, S1Angle.radians(PI).normalize().radians());

    // nextAfter PI should wrap around to _exactly_ nextAfter -Pi.
    assertExactly(
        nextAfter(-PI, 0.0),
        S1Angle.radians(nextAfter(PI, 4.0)).normalize().radians());

    // -PI is wrapped around to _exactly_ PI.
    assertExactly(PI, S1Angle.radians(-PI).normalize().radians());

    // nextAfter (downwards) -PI should wrap around to _exactly_ nextAfter (downwards) PI.
    assertExactly(nextAfter(PI, 0.0), S1Angle.radians(nextAfter(-PI, -4.0)).normalize().radians());
  }

  @Test
  public void testAlreadyNormalizedIsIdentity() {
    final S1Angle angle = S1Angle.degrees(90.0);
    assertSame(angle.normalize(), angle);
  }

  @SuppressWarnings("SelfEquals")
  @Test
  public void testInfinity() {
    assertTrue(S1Angle.radians(1e30).compareTo(S1Angle.INFINITY) < 0);
    assertTrue(S1Angle.INFINITY.neg().compareTo(S1Angle.ZERO) < 0);
    assertTrue(S1Angle.INFINITY.equals(S1Angle.INFINITY));
  }

  @Test
  public void testBuilder_zero() {
    assertEquals(S1Angle.ZERO, new S1Angle.Builder().build());
  }

  @Test
  public void testBuilder_radians() {
    assertExactly(1.5, new S1Angle.Builder().add(0.5).add(0.75).add(0.25).build().radians());
  }

  @Test
  public void testBuilder_angles() {
    S1Angle angle =
        new S1Angle.Builder()
            .add(S1Angle.degrees(90))
            .add(S1Angle.degrees(45))
            .add(S1Angle.degrees(30))
            .add(S1Angle.degrees(15))
            .build();
    assertExactly(180.0, angle.degrees());
  }

  @Test
  public void testBuilder_mixed() {
    S1Angle.Builder builder = new S1Angle.Builder();
    builder.add(PI / 2);
    builder.add(S1Angle.degrees(45));
    S1Angle angle1 = builder.build();
    assertExactly(0.75 * PI, angle1.radians());
    // Add more to the builder, make sure it didn't mutate the already-built value
    builder.add(PI);
    S1Angle angle2 = builder.build();
    assertExactly(1.75 * PI, angle2.radians());
    assertExactly(0.75 * PI, angle1.radians());
  }

  @Test
  public void testToBuilder() {
    S1Angle base = S1Angle.radians(PI);
    S1Angle.Builder builder = base.toBuilder();
    builder.add(S1Angle.degrees(90));
    assertExactly(1.5 * PI, builder.build().radians());
    assertExactly(PI, base.radians()); // Builder uses a copy, not original
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS1AngleSerialization() {
    doSerializationTest(S1Angle.radians(0.234));
  }
}
