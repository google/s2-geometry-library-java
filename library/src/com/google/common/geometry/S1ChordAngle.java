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
import static com.google.common.geometry.S2.DBL_EPSILON;
import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static java.lang.Math.sqrt;

import com.google.common.primitives.Doubles;
import java.io.Serializable;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * S1ChordAngle represents the angle subtended by a chord (i.e., the straight 3D Cartesian line
 * segment connecting two points on the unit sphere). Its representation makes it very efficient for
 * computing and comparing distances, but unlike S1Angle it is only capable of representing angles
 * between 0 and Pi radians. Generally, S1ChordAngle should only be used in loops where many angles
 * need to be calculated and compared. Otherwise it is simpler to use S1Angle.
 *
 * <p>S1ChordAngle also loses some accuracy as the angle approaches Pi radians. There are several
 * different ways to measure this error, including the representational error (i.e., how accurately
 * S1ChordAngle can represent angles near Pi radians), the conversion error (i.e., how much
 * precision is lost when an S1Angle is converted to an S1ChordAngle), and the measurement error
 * (i.e., how accurate the S1ChordAngle(a, b) constructor is when the points A and B are separated
 * by angles close to Pi radians). All of these errors differ by a small constant factor.
 *
 * <p>For the measurement error (which is the largest of these errors and also the most important in
 * practice), let the angle between A and B be (Pi - x) radians, i.e. A and B are within "x" radians
 * of being antipodal. The corresponding chord length is
 *
 * <pre>r = 2 * sin((Pi - x) / 2) = 2 * cos(x / 2) .</pre>
 *
 * <p>For values of x not close to Pi the relative error in the squared chord length is at most
 * {@code 4.5 * DBL_EPSILON} (see {@link getS2PointConstructorMaxError}). The relative error in "r"
 * is thus at most {@code 2.25 * DBL_EPSILON ~= 5e-16}. To convert this error into an equivalent
 * angle, we have
 *
 * <pre>|dr / dx| = sin(x / 2)</pre>
 *
 * <p>and therefore
 *
 * <pre>
 *    |dx| = dr / sin(x / 2)
 *         = 5e-16 * (2 * cos(x / 2)) / sin(x / 2)
 *         = 1e-15 / tan(x / 2)
 * </pre>
 *
 * <p>The maximum error is attained when
 *
 * <pre>
 *    x  = |dx|
 *       = 1e-15 / tan(x / 2)
 *      ~= 1e-15 / (x / 2)
 *      ~= sqrt(2e-15)
 * </pre>
 *
 * <p>In summary, the measurement error for an angle (Pi - x) is at most
 *
 * <pre>
 *    dx  = min(1e-15 / tan(x / 2), sqrt(2e-15))
 *      (~= min(2e-15 / x, sqrt(2e-15)) when x is small).
 * </pre>
 *
 * <p>On the Earth's surface (assuming a radius of 6371km), this corresponds to the following
 * worst-case measurement errors:
 *
 * <pre>
 *     Accuracy:             Unless antipodal to within:
 *     ---------             ---------------------------
 *     6.4 nanometers        10,000 km (90 degrees)
 *     1 micrometer          81.2 kilometers
 *     1 millimeter          81.2 meters
 *     1 centimeter          8.12 meters
 *     28.5 centimeters      28.5 centimeters
 * </pre>
 *
 * <p>The representational and conversion errors referred to earlier are somewhat smaller than this.
 * For example, maximum distance between adjacent representable S1ChordAngle values is only 13.5 cm
 * rather than 28.5 cm. To see this, observe that the closest representable value to r^2 = 4 is
 * {@code r^2 = 4 * (1 - DBL_EPSILON / 2)}. Thus {@code r = 2 * (1 - DBL_EPSILON / 4)} and the angle
 * between these two representable values is
 *
 * <pre>
 *    x  = 2 * acos(r / 2)
 *       = 2 * acos(1 - DBL_EPSILON / 4)
 *      ~= 2 * asin(sqrt(DBL_EPSILON / 2)
 *      ~= sqrt(2 * DBL_EPSILON)
 *      ~= 2.1e-8
 * </pre>
 *
 * <p>which is 13.5 cm on the Earth's surface.
 *
 * <p>The worst case rounding error occurs when the value halfway between these two representable
 * values is rounded up to 4. This halfway value is {@code r^2 = (4 * (1 - DBL_EPSILON / 4))}, thus
 * {@code r = 2 * (1 - DBL_EPSILON / 8)} and the worst case rounding error is
 *
 * <pre>
 *    x  = 2 * acos(r / 2)
 *       = 2 * acos(1 - DBL_EPSILON / 8)
 *      ~= 2 * asin(sqrt(DBL_EPSILON / 4)
 *      ~= sqrt(DBL_EPSILON)
 *      ~= 1.5e-8
 * </pre>
 *
 * which is 9.5 cm on the Earth's surface.
 */
@JsType
public final class S1ChordAngle implements S1Distance<S1ChordAngle>, Serializable {
  /**
   * Maximum relative error when summing two S1ChordAngles together. Absolute error is length2() of
   * the summed S1ChordAngle times this value. See {@link add(S1ChordAngle, S1ChordAngle)} for the
   * error analysis.
   */
  public static final double RELATIVE_SUM_ERROR = 2.02 * DBL_EPSILON;

  /**
   * Returns a new DistanceCollector for finding a minimum S1ChordAngle, starting from an initial
   * value of infinity.
   */
  public static final DistanceCollector<S1ChordAngle> minCollector() {
    return new DistanceCollector<S1ChordAngle>() {
      private S1ChordAngle min = S1ChordAngle.INFINITY;

      @Override
      public S1ChordAngle distance() {
        return min;
      }

      @Override
      public void reset() {
        min = S1ChordAngle.INFINITY;
      }

      @Override
      public void set(S1ChordAngle value) {
        min = value;
      }

      @Override
      public boolean update(S1ChordAngle other) {
        if (other.lessThan(min)) {
          min = other;
          return true;
        }
        return false;
      }

      @Override
      public boolean update(S2Point p1, S2Point p2) {
        // TODO(user): Benchmark using S2Predicates.compareDistance() to reject distances
        // further than the current best before accurately computing the actual distance.
        return update(new S1ChordAngle(p1, p2));
      }

      @Override
      public boolean update(S2Point p, S2Point v0, S2Point v1) {
        // TODO(user): Move or copy the S2EdgeUtil code into this class.
        return update(S2EdgeUtil.updateMinDistance(p, v0, v1, min));
      }

      @Override
      public boolean update(S2Point v0, S2Point v1, S2Point w0, S2Point w1) {
        // TODO(user): Move or copy the S2EdgeUtil code into this class.
        return update(S2EdgeUtil.getEdgePairMinDistance(v0, v1, w0, w1, min));
      }

      @Override
      public boolean update(S2Point p, S2Cell c) {
        return update(c.getDistance(p));
      }

      @Override
      public boolean update(S2Point v0, S2Point v1, S2Cell c) {
        return update(c.getDistanceToEdge(v0, v1));
      }

      @Override
      public boolean update(S2Cell c1, S2Cell c2) {
        return update(c1.getDistance(c2));
      }
    };
  }

  /**
   * Returns a new DistanceCollector for finding a maximum S1ChordAngle, starting from an initial
   * value of negative.
   */
  public static final DistanceCollector<S1ChordAngle> maxCollector() {
    return new DistanceCollector<S1ChordAngle>() {
      private S1ChordAngle max = S1ChordAngle.NEGATIVE;

      @Override
      public S1ChordAngle distance() {
        return max;
      }

      @Override
      public void reset() {
        max = S1ChordAngle.NEGATIVE;
      }

      @Override
      public void set(S1ChordAngle value) {
        max = value;
      }

      @Override
      public boolean update(S1ChordAngle other) {
        if (other.greaterThan(max)) {
          max = other;
          return true;
        }
        return false;
      }

      @Override
      public boolean update(S2Point p1, S2Point p2) {
        S1ChordAngle dist = new S1ChordAngle(p1, p2);
        // For greater accuracy and consistency with update(S2Point, S2Point, S2Point),
        // recompute large chord distances by subtracting the antipodal distance from STRAIGHT.
        if (dist.compareTo(S1ChordAngle.RIGHT) > 0) {
          dist = new S1ChordAngle(p1.neg(), p2);
          dist = S1ChordAngle.sub(S1ChordAngle.STRAIGHT, dist);
        }
        return update(dist);
      }

      @Override
      public boolean update(S2Point p, S2Point v0, S2Point v1) {
        // TODO(user): Move or copy the S2EdgeUtil code into this class.
        return update(S2EdgeUtil.updateMaxDistance(p, v0, v1, max));
      }

      @Override
      public boolean update(S2Point v0, S2Point v1, S2Point w0, S2Point w1) {
        // TODO(user): Move or copy the S2EdgeUtil code into this class.
        return update(S2EdgeUtil.getEdgePairMaxDistance(v0, v1, w0, w1, max));
      }

      @Override
      public boolean update(S2Point p, S2Cell c) {
        return update(c.getMaxDistance(p));
      }

      @Override
      public boolean update(S2Point v0, S2Point v1, S2Cell c) {
        return update(c.getMaxDistance(v0, v1));
      }

      @Override
      public boolean update(S2Cell c1, S2Cell c2) {
        return update(c1.getMaxDistance(c2));
      }
    };
  }

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
   * S1ChordAngles are represented by the squared chord length, which can range from 0 to {@code
   * MAX_LENGTH2}. {@link #INFINITY} uses an infinite squared length.
   */
  @JsConstructor
  private S1ChordAngle(double length2) {
    this.length2 = length2;
    if (!isValid()) {
      throw new IllegalArgumentException("Invalid length2: " + length2);
    }
  }

  /**
   * Constructs the S1ChordAngle corresponding to the distance between the two given points. The
   * points must be unit length.
   */
  @JsIgnore
  public S1ChordAngle(S2Point x, S2Point y) {
    // The distance may slightly exceed 4.0 due to roundoff errors.
    this(Math.min(MAX_LENGTH2, x.getDistance2(y)));
    checkArgument(S2.isUnitLength(x));
    checkArgument(S2.isUnitLength(y));
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
      double length = 2 * Math.sin(0.5 * Math.min(PI, angle.radians()));
      return new S1ChordAngle(length * length);
    }
  }

  /** Returns a new chord angle approximated from an angle specified in radians. */
  public static S1ChordAngle fromRadians(double radians) {
    return fromS1Angle(S1Angle.radians(radians));
  }

  /** Returns a new chord angle approximated from an angle specified in degrees. */
  public static S1ChordAngle fromDegrees(double degrees) {
    return fromS1Angle(S1Angle.degrees(degrees));
  }

  /** Returns a new chord angle approximated from an angle specified in tens of microdegrees. */
  public static S1ChordAngle fromE5(int e5) {
    return fromS1Angle(S1Angle.e5(e5));
  }

  /** Returns a new chord angle approximated from an angle specified in microdegrees. */
  public static S1ChordAngle fromE6(int e6) {
    return fromS1Angle(S1Angle.e6(e6));
  }

  /** Returns a new chord angle approximated from an angle specified in tenths of microdegrees. */
  public static S1ChordAngle fromE7(int e7) {
    return fromS1Angle(S1Angle.e7(e7));
  }

  /**
   * Returns an approximation of this chord angle in radians, implemented by first converting to an
   * S1Angle.
   */
  public double radians() {
    return toAngle().radians();
  }

  /**
   * Returns an approximation of this chord angle in degrees, implemented by first converting to an
   * S1Angle.
   */
  public double degrees() {
    return toAngle().degrees();
  }

  /**
   * Returns an approximation of this chord angle in tens of microdegrees, implemented by first
   * converting to an S1Angle.
   */
  public int e5() {
    return toAngle().e5();
  }

  /**
   * Returns an approximation of this chord angle in microdegrees, implemented by first converting
   * to an S1Angle.
   */
  public int e6() {
    return toAngle().e6();
  }

  /**
   * Returns an approximation of this chord angle in tenths of microdegrees, implemented by first
   * converting to an S1Angle.
   */
  public int e7() {
    return toAngle().e7();
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
  @Override
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

  /** Returns true iff this S1ChordAngle represents a distance less than 'other'. */
  @Override
  public boolean lessThan(S1ChordAngle other) {
    return length2 < other.length2;
  }

  /** Returns true iff this S1ChordAngle represents a distance greater than 'other'. */
  @Override
  public boolean greaterThan(S1ChordAngle other) {
    return length2 > other.length2;
  }

  /** Returns true iff this S1ChordAngle represents a distance less than or equal to 'other'. */
  @Override
  public boolean lessOrEquals(S1ChordAngle other) {
    return length2 <= other.length2;
  }

  /** Returns true iff this S1ChordAngle represents a distance greater than or equal to 'other'. */
  @Override
  public boolean greaterOrEquals(S1ChordAngle other) {
    return length2 >= other.length2;
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
   * convert a {@code <} comparison to a {@code <=} comparison.
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
   * represented by the 'a' and 'b' chord angles, with a ceiling of 180 degrees, the maximum value
   * of S1ChordAngle.
   *
   * <p>Note that this method is much more efficient than converting the chord angles to S1Angles
   * and adding those. It requires only one square root plus a few additions and multiplications.
   */
  public static S1ChordAngle add(S1ChordAngle a, S1ChordAngle b) {
    checkArgument(!a.isSpecial());
    checkArgument(!b.isSpecial());

    // Error Analysis:
    //
    //   u is the unit round-off, equal to ε/2 = 2^-53
    //   x and y below are both computed with (1+u)² error.  So we have
    //     length2 = x + y + 2*sqrt(x * y)
    //     length2 = ((x + y)(1+u)³ + 2*sqrt((x * y)((1+u)³))(1+u))(1+u)
    //
    //   Taking (1+u)^1.5 out of the square root and rounding (1+u)^2.5 to (1+u)³:
    //     length2 = (x + y + 2*sqrt(x * y))(1+u)⁴
    //
    //   Bounding (1+u)⁴ as 4*1.01u = 4.04*ε/2 = 2.02ε, total relative error is:
    //     length2()*2.02ε.

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
   * Subtract one S1ChordAngle from another, with a floor of zero as S1ChordAngles are non-negative.
   *
   * <p>Note that this method is much more efficient than converting the chord angles to S1Angles
   * and subtracting those. It requires only two square roots plus a few additions and
   * multiplications.
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

    // The calculation below is formulated differently (with two square roots rather than one) to
    // avoid excessive cancellation error when two nearly equal values are subtracted.
    double c = Math.max(0.0, sqrt(x) - sqrt(y));
    return new S1ChordAngle(c * c);
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
    // The sin() call and the multiply each have a relative error of 0.5 * DBL_EPSILON. However, the
    // sin() error is squared.
    return 1.5 * DBL_EPSILON * length2;
  }

  /**
   * There is a relative error of {@code 2.5 * DBL_EPSILON} when computing the squared distance,
   * plus a relative error of {@code 2 * DBL_EPSILON} and an absolute error of {@code 16 *
   * DBL_EPSILON^2} because the lengths of the input points may differ from 1 by up to {@code 2 *
   * DBL_EPSILON} each. (This is the maximum length error in {@link S2Point#normalize}).
   */
  public double getS2PointConstructorMaxError() {
    return (4.5 * DBL_EPSILON * length2) + (16 * DBL_EPSILON * DBL_EPSILON);
  }

  /** Returns the string of the closest {@link S1Angle} to this chord distance. */
  @Override
  public String toString() {
    if (length2 == NEGATIVE.length2) {
      return "NEGATIVE";
    }
    if (length2 == ZERO.length2) {
      return "ZERO";
    }
    if (length2 == STRAIGHT.length2) {
      return "STRAIGHT";
    }
    if (length2 == INFINITY.length2) {
      return "INFINITY";
    }
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
