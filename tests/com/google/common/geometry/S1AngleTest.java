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

import com.google.common.annotations.GwtCompatible;
import junit.framework.TestCase;

/** Tests for {@link S1Angle}. */
@GwtCompatible
public strictfp class S1AngleTest extends TestCase {

  private static final double EPSILON = 1e-12;

  public void testBasic() {
    // Check that the conversion between Pi radians and 180 degrees is exact.
    assertEquals(Math.PI, S1Angle.radians(Math.PI).radians());
    assertEquals(180.0, S1Angle.radians(Math.PI).degrees());
    assertEquals(Math.PI, S1Angle.degrees(180).radians());
    assertEquals(180.0, S1Angle.degrees(180).degrees());

    assertEquals(90.0, S1Angle.radians(Math.PI / 2).degrees());

    // Check negative angles.
    assertEquals(-90.0, S1Angle.radians(-Math.PI / 2).degrees());
    assertEquals(-Math.PI / 4, S1Angle.degrees(-45).radians());

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

  public void testMath() {
    assertEquals(30, S1Angle.degrees(10).add(S1Angle.degrees(20)).degrees(), EPSILON);
    assertEquals(-10, S1Angle.degrees(10).sub(S1Angle.degrees(20)).degrees(), EPSILON);
    assertEquals(20, S1Angle.degrees(10).mul(2.0).degrees(), EPSILON);
    assertEquals(5, S1Angle.degrees(10).div(2.0).degrees(), EPSILON);
    assertEquals(1.0, S1Angle.degrees(0).cos(), EPSILON);
    assertEquals(1.0, S1Angle.degrees(90).sin(), EPSILON);
    assertEquals(1.0, S1Angle.degrees(45).tan(), EPSILON);
  }

  public void testDistance() {
    // Check distance accessor for arbitrary sphere
    assertEquals(100.0 * Math.PI, S1Angle.radians(Math.PI).distance(100.0), EPSILON);
    assertEquals(50.0 * Math.PI, S1Angle.radians(Math.PI / 2).distance(100.0), EPSILON);
    assertEquals(25.0 * Math.PI, S1Angle.radians(Math.PI / 4).distance(100.0), EPSILON);

  }

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

  public void testNormalize() {
    assertEquals(0.0, S1Angle.degrees(360.0).normalize().degrees(), EPSILON);
    assertEquals(-90.0, S1Angle.degrees(-90.0).normalize().degrees(), EPSILON);
    assertEquals(180.0, S1Angle.degrees(-180.0).normalize().degrees(), EPSILON);
    assertEquals(90.0, S1Angle.degrees(90.0).normalize().degrees(), EPSILON);
    assertEquals(180.0, S1Angle.degrees(180.0).normalize().degrees(), EPSILON);
    assertEquals(-90.0, S1Angle.degrees(270.0).normalize().degrees(), EPSILON);
    assertEquals(180.0, S1Angle.degrees(540.0).normalize().degrees(), EPSILON);
    assertEquals(90.0, S1Angle.degrees(-270.0).normalize().degrees(), EPSILON);

    assertEquals(Math.PI, S1Angle.radians(Math.PI).normalize().radians(), EPSILON);
    // Should wrap around.
    assertEquals(
        -Math.PI, S1Angle.radians(Platform.nextAfter(Math.PI, 4.0)).normalize().radians(), EPSILON);
    assertEquals(
        Math.PI, S1Angle.radians(Platform.nextAfter(Math.PI, 0.0)).normalize().radians(), EPSILON);

    // -PI maps to PI.
    assertEquals(Math.PI, S1Angle.radians(-Math.PI).normalize().radians(), EPSILON);
    // Should wrap around.
    assertEquals(
        Math.PI,
        S1Angle.radians(Platform.nextAfter(-Math.PI, -4.0)).normalize().radians(),
        EPSILON);
    assertEquals(
        -Math.PI,
        S1Angle.radians(Platform.nextAfter(-Math.PI, 0.0)).normalize().radians(),
        EPSILON);
  }

  public void testAlreadyNormalizedIsIdentity() {
    final S1Angle angle = S1Angle.degrees(90.0);
    assertSame(angle.normalize(), angle);
  }

  @SuppressWarnings("SelfEquals")
  public void testInfinity() {
    assertTrue(S1Angle.radians(1e30).compareTo(S1Angle.INFINITY) < 0);
    assertTrue(S1Angle.INFINITY.neg().compareTo(S1Angle.ZERO) < 0);
    assertTrue(S1Angle.INFINITY.equals(S1Angle.INFINITY));
  }

  public void testBuilder_zero() {
    assertEquals(S1Angle.ZERO, new S1Angle.Builder().build());
  }

  public void testBuilder_radians() {
    assertEquals(1.5, new S1Angle.Builder().add(0.5).add(0.75).add(0.25).build().radians());
  }

  public void testBuilder_angles() {
    S1Angle angle =
        new S1Angle.Builder()
            .add(S1Angle.degrees(90))
            .add(S1Angle.degrees(45))
            .add(S1Angle.degrees(30))
            .add(S1Angle.degrees(15))
            .build();
    assertEquals(180.0, angle.degrees());
  }

  public void testBuilder_mixed() {
    S1Angle.Builder builder = new S1Angle.Builder();
    builder.add(Math.PI / 2);
    builder.add(S1Angle.degrees(45));
    S1Angle angle1 = builder.build();
    assertEquals(0.75 * Math.PI, angle1.radians());
    // Add more to the builder, make sure it didn't mutate the already-built value
    builder.add(Math.PI);
    S1Angle angle2 = builder.build();
    assertEquals(1.75 * Math.PI, angle2.radians());
    assertEquals(0.75 * Math.PI, angle1.radians());
  }

  public void testToBuilder() {
    S1Angle base = S1Angle.radians(Math.PI);
    S1Angle.Builder builder = base.toBuilder();
    builder.add(S1Angle.degrees(90));
    assertEquals(1.5 * Math.PI, builder.build().radians());
    assertEquals(Math.PI, base.radians()); // Builder uses a copy, not original
  }
}
