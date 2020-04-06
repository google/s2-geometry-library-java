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

import static com.google.common.base.Preconditions.checkArgument;
import static java.lang.Math.asin;
import static java.lang.Math.sqrt;

import com.google.common.annotations.GwtCompatible;
import com.google.common.primitives.Doubles;
import java.io.Serializable;

/**
 * S1ChordAngle represents the angle subtended by a chord (i.e., the straight 3D Cartesian line
 * segment connecting two points on the unit sphere). Its representation makes it very efficient for
 * computing and comparing distances, but unlike S1Angle it is only capable of representing angles
 * between 0 and Pi radians. Generally, S1ChordAngle should only be used in loops where many angles
 * need to be calculated and compared. Otherwise it is simpler to use S1Angle.
 *
 * <p>S1ChordAngle also loses some accuracy as the angle approaches Pi radians. Specifically, the
 * representation of (Pi - x) radians can be expected to have an error of about (1e-15 / x), with a
 * maximum error of about 1e-7.
 */
@GwtCompatible(serializable = true)
public final strictfp class S1ChordAngle implements Comparable<S1ChordAngle>, Serializable {

  /** Max value that can be returned from {@link #getLength2()}. */
  public static final double MAX_LENGTH2 = 4.0;

  /** The zero chord angle. */
  public static final S1ChordAngle ZERO = new S1ChordAngle(0);

  /** The chord angle of 90 degrees (a "right angle"). */
  public static final S1ChordAngle RIGHT = new S1ChordAngle(2);

  /** The chord angle of 180 degrees (a "straight angle"). This is the max finite chord angle. */
  public static final S1ChordAngle STRAIGHT = new S1ChordAngle(MAX_LENGTH2);

  /**
   * A chord angle larger than any finite chord angle. The only valid operations on {@code INFINITY}
   * are comparisons and {@link S1Angle} conversions.
   */
  public static final S1ChordAngle INFINITY = new S1ChordAngle(Double.POSITIVE_INFINITY);

  /**
   * A chord angle smaller than {@link #ZERO}. The only valid operations on {@code NEGATIVE} are
   * comparisons and {@link S1Angle} conversions.
   */
  public static final S1ChordAngle NEGATIVE = new S1ChordAngle(-1);

  private final double length2;

  /**
   * Constructs the S1ChordAngle corresponding to the distance between the two given points. The
   * points must be unit length.
   */
  public S1ChordAngle(S2Point x, S2Point y) {
    checkArgument(S2.isUnitLength(x));
    checkArgument(S2.isUnitLength(y));
    // The distance may slightly exceed 4.0 due to roundoff errors.
    length2 = Math.min(MAX_LENGTH2, x.getDistance2(y));
    checkArgument(isValid());
  }

  /**
   * Returns a new chord angle approximated from {@code angle} (see {@link
   * #getS1AngleConstructorMaxError()} for the max magnitude of the error).
   *
   * <p>Angles outside the range [0, Pi] are handled as follows:
   *
   * <ul>
   *   <li>{@link S1Angle#INFINITY} is mapped to {@link #INFINITY}
   *   <li>negative angles are mapped to {@link #NEGATIVE}
   *   <li>finite angles larger than Pi are mapped to {@link #STRAIGHT}
   * </ul>
   *
   * <p>Note that this operation is relatively expensive and should be avoided. To use {@link
   * S1ChordAngle} effectively, you should structure your code so that input arguments are converted
   * to S1ChordAngles at the beginning of your algorithm, and results are converted back to {@link
   * S1Angle}s only at the end.
   */
  public static S1ChordAngle fromS1Angle(S1Angle angle) {
    if (angle.radians() < 0) {
      return NEGATIVE;
    } else if (angle.equals(S1Angle.INFINITY)) {
      return INFINITY;
    } else {
      // The chord length is 2 * sin(angle / 2).
      double length = 2 * Math.sin(0.5 * Math.min(Math.PI, angle.radians()));
      return new S1ChordAngle(length * length);
    }
  }

  /**
   * S1ChordAngles are represented by the squared chord length, which can range from 0 to {@code
   * MAX_LENGTH2}. {@link #INFINITY} uses an infinite squared length.
   */
  private S1ChordAngle(double length2) {
    this.length2 = length2;
    checkArgument(isValid());
  }

  /**
   * Construct an S1ChordAngle from the squared chord length. Note that the argument is
   * automatically clamped to a maximum of {@code MAX_LENGTH2} to handle possible roundoff errors.
   * The argument must be non-negative.
   */
  public static S1ChordAngle fromLength2(double length2) {
    return new S1ChordAngle(Math.min(MAX_LENGTH2, length2));
  }

  /** Returns whether the chord distance is exactly 0. */
  public boolean isZero() {
    return length2 == 0;
  }

  /** Returns whether the chord distance is negative. */
  public boolean isNegative() {
    return length2 < 0;
  }

  /** Returns whether the chord distance is exactly (positive) infinity. */
  public boolean isInfinity() {
    return length2 == Double.POSITIVE_INFINITY;
  }

  /** Returns true if the angle is negative or infinity. */
  public boolean isSpecial() {
    return isNegative() || isInfinity();
  }

  /**
   * Returns true if getLength2() is within the normal range of 0 to 4 (inclusive) or the angle is
   * special.
   */
  public boolean isValid() {
    return (length2 >= 0 && length2 <= MAX_LENGTH2) || isNegative() || isInfinity();
  }

  /**
   * Convert the chord angle to an {@link S1Angle}. {@link #INFINITY} is converted to {@link
   * S1Angle#INFINITY}, and {@link #NEGATIVE} is converted to a negative {@link S1Angle}. This
   * operation is relatively expensive.
   */
  public S1Angle toAngle() {
    if (isNegative()) {
      return S1Angle.radians(-1);
    } else if (isInfinity()) {
      return S1Angle.INFINITY;
    } else {
      return S1Angle.radians(2 * asin(0.5 * sqrt(length2)));
    }
  }

  /** The squared length of the chord. (Most clients will not need this.) */
  public double getLength2() {
    return length2;
  }

  /**
   * Returns the smallest representable S1ChordAngle larger than this object. This can be used to
   * convert a "<" comparison to a "<=" comparison.
   *
   * <p>Note the following special cases:
   *
   * <ul>
   *   <li>NEGATIVE.successor() == ZERO
   *   <li>STRAIGHT.successor() == INFINITY
   *   <li>INFINITY.Successor() == INFINITY
   * </ul>
   */
  public S1ChordAngle successor() {
    if (length2 >= MAX_LENGTH2) {
      return INFINITY;
    }
    if (length2 < 0.0) {
      return ZERO;
    }
    return new S1ChordAngle(Platform.nextAfter(length2, 10.0));
  }

  /**
   * As {@link #successor}, but returns the largest representable S1ChordAngle less than this
   * object.
   *
   * <p>Note the following special cases:
   *
   * <ul>
   *   <li>INFINITY.predecessor() == STRAIGHT
   *   <li>ZERO.predecessor() == NEGATIVE
   *   <li>NEGATIVE.predecessor() == NEGATIVE
   * </ul>
   */
  public S1ChordAngle predecessor() {
    if (length2 <= 0.0) {
      return NEGATIVE;
    }
    if (length2 > MAX_LENGTH2) {
      return STRAIGHT;
    }
    return new S1ChordAngle(Platform.nextAfter(length2, -10.0));
  }

  /**
   * Returns a new S1ChordAngle whose chord distance represents the sum of the angular distances
   * represented by the 'a' and 'b' chord angles.
   *
   * <p>Note that this method is much more efficient than converting the chord angles to S1Angles
   * and adding those. It requires only one square root plus a few additions and multiplications.
   */
  public static S1ChordAngle add(S1ChordAngle a, S1ChordAngle b) {
    checkArgument(!a.isSpecial());
    checkArgument(!b.isSpecial());

    // Optimization for the common case where "b" is an error tolerance parameter that happens to be
    // set to zero.
    double a2 = a.length2;
    double b2 = b.length2;
    if (b2 == 0) {
      return a;
    }

    // Clamp the angle sum to at most 180 degrees.
    if (a2 + b2 >= MAX_LENGTH2) {
      return S1ChordAngle.STRAIGHT;
    }

    // Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
    // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
    // Then the formula below can be derived from c = 2 * sin(A+B) and the relationships
    //   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
    //   cos(X) = sqrt(1 - sin^2(X)) .
    double x = a2 * (1 - 0.25 * b2); // isValid() => non-negative
    double y = b2 * (1 - 0.25 * a2); // isValid() => non-negative
    return new S1ChordAngle(Math.min(MAX_LENGTH2, x + y + 2 * sqrt(x * y)));
  }

  /**
   * Subtract one S1ChordAngle from another.
   *
   * <p>Note that this method is much more efficient than converting the chord angles to S1Angles
   * and adding those. It requires only one square root plus a few additions and multiplications.
   */
  public static S1ChordAngle sub(S1ChordAngle a, S1ChordAngle b) {
    // See comments in add(S1ChordAngle, S1ChordAngle).
    checkArgument(!a.isSpecial());
    checkArgument(!b.isSpecial());
    double a2 = a.length2;
    double b2 = b.length2;
    if (b2 == 0) {
      return a;
    }
    if (a2 <= b2) {
      return S1ChordAngle.ZERO;
    }
    double x = a2 * (1 - 0.25 * b2);
    double y = b2 * (1 - 0.25 * a2);
    return new S1ChordAngle(Math.max(0.0, x + y - 2 * sqrt(x * y)));
  }

  /** Returns the smaller of the given instances. */
  public static S1ChordAngle min(S1ChordAngle a, S1ChordAngle b) {
    return a.length2 <= b.length2 ? a : b;
  }

  /** Returns the larger of the given instances. */
  public static S1ChordAngle max(S1ChordAngle a, S1ChordAngle b) {
    return a.length2 > b.length2 ? a : b;
  }

  /** Returns the square of Math.sin(toAngle().radians()), but computed more efficiently. */
  public static double sin2(S1ChordAngle a) {
    checkArgument(!a.isSpecial());
    // Let "a" be the (non-squared) chord length, and let A be the corresponding half-angle
    // (a = 2*sin(A)). The formula below can be derived from:
    //   sin(2*A) = 2 * sin(A) * cos(A)
    //   cos^2(A) = 1 - sin^2(A)
    // This is much faster than converting to an angle and computing its sine.
    return a.length2 * (1 - 0.25 * a.length2);
  }

  /** Returns Math.sin(toAngle().radians()), but computed more efficiently. */
  public static double sin(S1ChordAngle a) {
    return sqrt(sin2(a));
  }

  /** Returns Math.cos(toAngle().radians()), but computed more efficiently. */
  public static double cos(S1ChordAngle a) {
    // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
    checkArgument(!a.isSpecial());
    return 1 - 0.5 * a.length2;
  }

  /** Returns Math.tan(toAngle().radians()), but computed more efficiently. */
  public static double tan(S1ChordAngle a) {
    return sin(a) / cos(a);
  }

  /**
   * Returns a new S1ChordAngle that has been adjusted by the given error bound (which can be
   * positive or negative). {@code error} should be the value returned by one of the error bound
   * methods below. For example:
   *
   * <pre>
   *    {@code S1ChordAngle a = new S1ChordAngle(x, y);}
   *    {@code S1ChordAngle a1 = a.plusError(a.getS2PointConstructorMaxError());}
   * </pre>
   *
   * <p>If this {@link #isSpecial}, we return {@code this}.
   */
  public S1ChordAngle plusError(double error) {
    return isSpecial() ? this : fromLength2(Math.max(0.0, Math.min(MAX_LENGTH2, length2 + error)));
  }

  /** Returns the error in {@link #fromS1Angle}. */
  public double getS1AngleConstructorMaxError() {
    return S2.DBL_EPSILON * length2;
  }

  /**
   * There is a relative error of {@code 2.5 * DBL_EPSILON} when computing the squared distance,
   * plus a relative error of {@code 2 * DBL_EPSILON} and an absolute error of {@code 16 *
   * DBL_EPSILON^2} because the lengths of the input points may differ from 1 by up to {@code 2 *
   * DBL_EPSILON} each. (This is the maximum length error in {@link S2Point#normalize}).
   */
  public double getS2PointConstructorMaxError() {
    return (4.5 * S2.DBL_EPSILON * length2) + (16 * S2.DBL_EPSILON * S2.DBL_EPSILON);
  }

  /** Returns the string of the closest {@link S1Angle} to this chord distance. */
  @Override
  public String toString() {
    return toAngle().toString();
  }

  @Override
  public int compareTo(S1ChordAngle that) {
    return Double.compare(this.length2, that.length2);
  }

  @Override
  public boolean equals(Object other) {
    return (other instanceof S1ChordAngle) && length2 == ((S1ChordAngle) other).length2;
  }

  @Override
  public int hashCode() {
    return length2 == 0.0 ? 0 : Doubles.hashCode(length2);
  }
}
