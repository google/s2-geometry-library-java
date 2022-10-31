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

import static java.lang.Math.PI;
import static java.lang.Math.round;

import com.google.common.primitives.Ints;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;
import java.io.Serializable;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * This class represents a one-dimensional angle (as opposed to a two-dimensional solid angle). It
 * has methods for converting angles to or from radians, degrees, and the E6/E7 representations
 * (i.e. degrees multiplied by 1e6/1e7 and rounded to the nearest integer).
 *
 * <p>The main differences between this class and a naive Angle implementation are:
 *
 * <ul>
 *   <li>More accurate conversion between radians and degrees. The {@code Angle<>} classes use a
 *       rational approximation to Pi/180 that results in errors in the two least-significant bits.
 *   <li>Built-in support for the E6 and E7 representations.
 *   <li>A simple interface and implementation.
 * </ul>
 *
 * @author danieldanciu@google.com (Daniel Danciu) ported from util/geometry
 */
@JsType
public final strictfp class S1Angle implements S1Distance<S1Angle>, Serializable {

  /** An angle larger than any finite angle. */
  public static final S1Angle INFINITY = new S1Angle(Double.POSITIVE_INFINITY);

  /** An explicit shorthand for the default constructor. */
  public static final S1Angle ZERO = new S1Angle();

  private final double radians;

  /** Returns a new S1Angle specified in radians. */
  public static S1Angle radians(double radians) {
    return new S1Angle(radians);
  }

  /** Returns the angle in radians. */
  public double radians() {
    return radians;
  }

  /**
   * Returns a new S1Angle converted from degrees. Note that <code>degrees(x).degrees() == x</code>
   * may not hold due to inexact arithmetic.
   */
  public static S1Angle degrees(double degrees) {
    return new S1Angle(degrees * (PI / 180));
  }

  /** Returns the angle in degrees. */
  public double degrees() {
    return radians * (180 / PI);
  }

  /** Returns a new S1Angle converted from tens of microdegrees. */
  public static S1Angle e5(int e5) {
    return degrees(e5 * 1e-5);
  }

  /**
   * Returns angle in tens of microdegrees, rounded to the nearest ten microdegrees.
   *
   * <p>Normalized angles will never overflow an int.
   *
   * @throws IllegalArgumentException if the result overflows an int
   */
  public int e5() {
    return Ints.checkedCast(round(degrees() * 1e5));
  }

  /** Returns a new S1Angle converted from microdegrees. */
  public static S1Angle e6(int e6) {
    // Multiplying by 1e-6 isn't quite as accurate as dividing by 1e6, but it's about 10 times
    // faster and more than accurate enough.
    return degrees(e6 * 1e-6);
  }

  /**
   * Returns angle in microdegrees, rounded to the nearest microdegree.
   *
   * <p>Normalized angles will never overflow an int.
   *
   * @throws IllegalArgumentException if the result overflows an int
   */
  public int e6() {
    return Ints.checkedCast(round(degrees() * 1e6));
  }

  /** Returns a new S1Angle converted from tenths of a microdegree. */
  public static S1Angle e7(int e7) {
    return degrees(e7 * 1e-7);
  }

  /**
   * Returns angle in tenths of a microdegree, rounded to the nearest tenth of a microdegree.
   *
   * <p>Normalized angles will never overflow an int.
   *
   * @throws IllegalArgumentException if the result overflows an int
   */
  public int e7() {
    return Ints.checkedCast(round(degrees() * 1e7));
  }

  /** The default constructor yields a zero angle. */
  @JsIgnore
  public S1Angle() {
    this(0);
  }

  @JsConstructor
  private S1Angle(double radians) {
    this.radians = radians;
  }

  /**
   * Return the angle between two points, which is also equal to the distance between these points
   * on the unit sphere. The points do not need to be normalized.
   */
  @JsIgnore
  public S1Angle(S2Point x, S2Point y) {
    this(x.angle(y));
  }

  @Override
  public boolean equals(Object that) {
    if (that instanceof S1Angle) {
      return this.radians == ((S1Angle) that).radians;
    }
    return false;
  }

  @Override
  public boolean isZero() {
    return radians == 0;
  }

  @Override
  public int hashCode() {
    long value = Double.doubleToLongBits(radians);
    return (int) (value ^ (value >>> 32));
  }

  @Override
  public boolean lessThan(S1Angle that) {
    return this.radians < that.radians;
  }

  @Override
  public boolean greaterThan(S1Angle that) {
    return this.radians > that.radians;
  }

  @Override
  public boolean lessOrEquals(S1Angle that) {
    return this.radians <= that.radians;
  }

  @Override
  public boolean greaterOrEquals(S1Angle that) {
    return this.radians >= that.radians;
  }

  public static S1Angle max(S1Angle left, S1Angle right) {
    return right.greaterThan(left) ? right : left;
  }

  public static S1Angle min(S1Angle left, S1Angle right) {
    return right.greaterThan(left) ? left : right;
  }

  /** Returns the distance along the surface of a sphere of the given radius. */
  public double distance(double radius) {
    return radians * radius;
  }

  public S1Angle neg() {
    return new S1Angle(-radians);
  }

  /** Returns an {@link S1Angle} whose angle is <code>(this + a)</code>. */
  @CheckReturnValue
  public S1Angle add(S1Angle a) {
    return new S1Angle(radians + a.radians);
  }

  /** Returns an {@link S1Angle} whose angle is <code>(this - a)</code>. */
  @CheckReturnValue
  public S1Angle sub(S1Angle a) {
    return new S1Angle(radians - a.radians);
  }

  /** Returns an {@link S1Angle} whose angle is <code>(this * m)</code>. */
  @CheckReturnValue
  public S1Angle mul(double m) {
    return new S1Angle(radians * m);
  }

  /** Returns an {@link S1Angle} whose angle is <code>(this / d)</code>. */
  @CheckReturnValue
  public S1Angle div(double d) {
    return new S1Angle(radians / d);
  }

  /**
   * Returns the trigonometric cosine of the angle.
   */
  public double cos() {
    return Math.cos(radians);
  }

  /**
   * Returns the trigonometric sine of the angle.
   */
  public double sin() {
    return Math.sin(radians);
  }

  /**
   * Returns the trigonometric tangent of the angle.
   */
  public double tan() {
    return Math.tan(radians);
  }

  /**
   * Returns the angle normalized to the range (-180, 180] degrees.
   */
  @CheckReturnValue
  public S1Angle normalize() {
    final boolean isNormalized = radians > -PI && radians <= PI;
    if (isNormalized) {
      return this;
    }
    double normalized = Platform.IEEEremainder(radians, 2.0 * PI);
    if (normalized <= -PI) {
      normalized = PI;
    }
    assert normalized > -PI;
    assert normalized <= PI;
    return new S1Angle(normalized);
  }

  /**
   * Writes the angle in degrees with a "d" suffix, e.g. "17.3745d". By default 6 digits are
   * printed; this can be changed using setprecision(). Up to 17 digits are required to distinguish
   * one angle from another.
   */
  @Override
  public String toString() {
    return degrees() + "d";
  }

  @Override
  public int compareTo(S1Angle that) {
    return this.radians < that.radians ? -1 : this.radians > that.radians ? 1 : 0;
  }

  /** Creates a new Builder initialized to a copy of this angle. */
  public Builder toBuilder() {
    return new Builder().add(this);
  }

  /** A builder of {@link S1Angle} instances. */
  @JsType
  public static final class Builder {
    private double radians;

    /** Constructs a new builder initialized to {@link #ZERO}. */
    public Builder() {}

    /** Adds angle. */
    @CanIgnoreReturnValue
    public Builder add(S1Angle angle) {
      radians += angle.radians;
      return this;
    }

    /** Adds radians. */
    @CanIgnoreReturnValue
    @JsMethod(name = "addRadians")
    public Builder add(double radians) {
      this.radians += radians;
      return this;
    }

    /** Returns a new {@link S1Angle} copied from the current state of this builder. */
    public S1Angle build() {
      return new S1Angle(radians);
    }
  }
}
