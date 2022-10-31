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
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.max;

import com.google.errorprone.annotations.CheckReturnValue;
import java.io.Serializable;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * An S1Interval represents a closed interval on a unit circle (also known as a 1-dimensional
 * sphere). It is capable of representing the empty interval (containing no points), the full
 * interval (containing all points), and zero-length intervals (containing a single point).
 *
 * <p>Points are represented by the angle they make with the positive x-axis in the range [-Pi, Pi].
 * An interval is represented by its lower and upper bounds (both inclusive, since the interval is
 * closed). The lower bound may be greater than the upper bound, in which case the interval is
 * "inverted" (i.e. it passes through the point (-1, 0)).
 *
 * <p>Note that the point (-1, 0) has two valid representations, Pi and -Pi. The normalized
 * representation of this point internally is Pi, so that endpoints of normal intervals are in the
 * range (-Pi, Pi]. However, we take advantage of the point -Pi to construct two special intervals:
 * the full() interval is [-Pi, Pi], and the Empty() interval is [Pi, -Pi].
 *
 * @author danieldanciu@google.com (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
public final strictfp class S1Interval implements Cloneable, Serializable {
  private double lo;
  private double hi;

  @JsIgnore
  public S1Interval() {
    this(PI, -PI);
  }

  /**
   * Both endpoints must be in the range -Pi to Pi inclusive. The value -Pi is converted internally
   * to Pi except for the full() and empty() intervals.
   */
  @JsConstructor
  public S1Interval(double lo, double hi) {
    set(lo, hi, false);
  }

  /** Copy constructor. */
  @JsIgnore
  public S1Interval(S1Interval interval) {
    this(interval.lo, interval.hi);
  }

  /**
   * Creates an S1Interval and passes the arguments down to {@link #set(double, double, boolean)}.
   */
  private static S1Interval fromEndpointsChecked(double lo, double hi, boolean checked) {
    S1Interval newInterval = new S1Interval(lo, hi);
    newInterval.set(lo, hi, checked);
    return newInterval;
  }

  /**
   * Assigns the range of this interval, assuming both arguments are in the correct range, i.e.
   * normalization from -Pi to Pi is already done. If {@code checked} is false, endpoints at -Pi
   * will be moved to +Pi unless the other endpoint is already there.
   *
   * <p>Note that because S1Interval has invariants to maintain after each update, values cannot be
   * set singly, both endpoints must be set together.
   */
  void set(double newLo, double newHi, boolean checked) {
    this.lo = newLo;
    this.hi = newHi;
    if (!checked) {
      if (newLo == -PI && newHi != PI) {
        lo = PI;
      }
      if (newHi == -PI && newLo != PI) {
        hi = PI;
      }
    }
  }

  /**
   * Sets the range of this interval to the empty interval.
   *
   * <p>Package private since only S2 code needs to mutate S1Intervals for now.
   */
  void setEmpty() {
    lo = PI;
    hi = -PI;
  }

  /**
   * Sets the range of this interval to the full interval.
   *
   * <p>Package private since only S2 code needs to mutate S1Intervals for now.
   */
  void setFull() {
    lo = -PI;
    hi = PI;
  }

  public static S1Interval empty() {
    S1Interval result = new S1Interval();
    result.setEmpty();
    return result;
  }

  public static S1Interval full() {
    S1Interval result = new S1Interval();
    result.setFull();
    return result;
  }

  /** Convenience method to construct an interval containing a single point. */
  public static S1Interval fromPoint(double radians) {
    if (radians == -PI) {
      radians = PI;
    }
    return fromEndpointsChecked(radians, radians, true);
  }

  /**
   * Convenience method to construct the minimal interval containing the two given points. This is
   * equivalent to starting with an empty interval and calling addPoint() twice, but it is more
   * efficient.
   */
  public static S1Interval fromPointPair(double p1, double p2) {
    // assert abs(p1) <= PI && abs(p2) <= PI;
    S1Interval result = new S1Interval();
    result.initFromPointPair(p1, p2);
    return result;
  }

  void initFromPointPair(double p1, double p2) {
    if (p1 == -PI) {
      p1 = PI;
    }
    if (p2 == -PI) {
      p2 = PI;
    }
    if (positiveDistance(p1, p2) <= PI) {
      this.lo = p1;
      this.hi = p2;
    } else {
      this.lo = p2;
      this.hi = p1;
    }
  }

  public double lo() {
    return lo;
  }

  public double hi() {
    return hi;
  }

  /**
   * An interval is valid if neither bound exceeds Pi in absolute value, and the value -Pi appears
   * only in the Empty() and full() intervals.
   */
  public boolean isValid() {
    return (abs(lo) <= PI && abs(hi) <= PI && !(lo == -PI && hi != PI) && !(hi == -PI && lo != PI));
  }

  /** Returns true if the interval contains all points on the unit circle. */
  public boolean isFull() {
    return hi - lo == 2 * PI;
  }

  /** Returns true if the interval is empty, i.e. it contains no points. */
  public boolean isEmpty() {
    return lo - hi == 2 * PI;
  }

  /** Returns true if {@code lo() > hi()}. (This is true for empty intervals.) */
  public boolean isInverted() {
    return lo > hi;
  }

  /**
   * Returns the midpoint of the interval. For full and empty intervals, the result is arbitrary.
   */
  public double getCenter() {
    double center = 0.5 * (lo + hi);
    if (!isInverted()) {
      return center;
    }
    // Return the center in the range (-Pi, Pi].
    return (center <= 0) ? (center + PI) : (center - PI);
  }

  /** Returns the length of the interval. The length of an empty interval is negative. */
  public double getLength() {
    double length = hi - lo;
    if (length >= 0) {
      return length;
    }
    length += 2 * PI;
    // Empty intervals have a negative length.
    return (length > 0) ? length : -1;
  }

  /**
   * Return the complement of the interior of the interval. An interval and its complement have the
   * same boundary but do not share any interior values. The complement operator is not a bijection,
   * since the complement of a singleton interval (containing a single value) is the same as the
   * complement of an empty interval.
   */
  public S1Interval complement() {
    if (lo == hi) {
      return full(); // Singleton.
    }
    return fromEndpointsChecked(hi, lo, true); // Handles empty and full.
  }

  /**
   * Return the midpoint of the complement of the interval. For full and empty intervals, the result
   * is arbitrary. For a singleton interval (containing a single point), the result is its antipodal
   * point on S1.
   */
  public double getComplementCenter() {
    if (lo() != hi()) {
      return complement().getCenter();
    } else { // Singleton.
      return (hi() <= 0) ? (hi() + PI) : (hi() - PI);
    }
  }

  /** Returns true if the interval (which is closed) contains the point 'p'. */
  @JsMethod(name = "containsPoint")
  public boolean contains(double p) {
    // Works for empty, full, and singleton intervals.
    // assert abs(p) <= PI;
    if (p == -PI) {
      p = PI;
    }
    return fastContains(p);
  }

  /**
   * Returns true if the interval contains the interval {@code y}. Works for empty, full, and
   * singleton intervals.
   */
  public boolean contains(final S1Interval y) {
    // It might be helpful to compare the structure of these tests to the simpler
    // fastContains(double) method below.

    if (isInverted()) {
      if (y.isInverted()) {
        return y.lo >= lo && y.hi <= hi;
      }
      return (y.lo >= lo || y.hi <= hi) && !isEmpty();
    } else {
      if (y.isInverted()) {
        return isFull() || y.isEmpty();
      }
      return y.lo >= lo && y.hi <= hi;
    }
  }

  /**
   * Returns true if the interval (which is closed) contains the point 'p'. Skips the normalization
   * of 'p' from -Pi to Pi.
   */
  public boolean fastContains(double p) {
    if (isInverted()) {
      return (p >= lo || p <= hi) && !isEmpty();
    } else {
      return p >= lo && p <= hi;
    }
  }

  /** Returns true if the interior of the interval contains the point 'p'. */
  @JsMethod(name = "interiorContainsPoint")
  public boolean interiorContains(double p) {
    // Works for empty, full, and singleton intervals.
    // assert abs(p) <= PI;
    if (p == -PI) {
      p = PI;
    }

    if (isInverted()) {
      return p > lo || p < hi;
    } else {
      return (p > lo && p < hi) || isFull();
    }
  }

  /**
   * Returns true if the interior of this interval contains the entire interval 'y'. Note that
   * x.interiorContains(x) is true only when x is the empty or full interval, and
   * x.interiorContains(S1Interval(p,p)) is equivalent to x.InteriorContains(p).
   */
  public boolean interiorContains(final S1Interval y) {
    if (isInverted()) {
      if (!y.isInverted()) {
        return y.lo > lo || y.hi < hi;
      }
      return (y.lo > lo && y.hi < hi) || y.isEmpty();
    } else {
      if (y.isInverted()) {
        return isFull() || y.isEmpty();
      }
      return (y.lo > lo && y.hi < hi) || isFull();
    }
  }

  /**
   * Returns true if the two intervals contain any points in common. Note that the point +/-Pi has
   * two representations, so the intervals [-Pi,-3] and [2,Pi] intersect, for example.
   */
  public boolean intersects(final S1Interval y) {
    if (isEmpty() || y.isEmpty()) {
      return false;
    }
    if (isInverted()) {
      // Every non-empty inverted interval contains Pi.
      return y.isInverted() || y.lo <= hi || y.hi >= lo;
    } else {
      if (y.isInverted()) {
        return y.lo <= hi || y.hi >= lo;
      }
      return y.lo <= hi && y.hi >= lo;
    }
  }

  /**
   * Returns true if the interior of this interval contains any point of the interval {@code y}
   * (including its boundary). Works for empty, full, and singleton intervals.
   */
  public boolean interiorIntersects(final S1Interval y) {
    if (isEmpty() || y.isEmpty() || lo == hi) {
      return false;
    }
    if (isInverted()) {
      return y.isInverted() || y.lo < hi || y.hi > lo;
    } else {
      if (y.isInverted()) {
        return y.lo < hi || y.hi > lo;
      }
      return (y.lo < hi && y.hi > lo) || isFull();
    }
  }

  /**
   * Return the Hausdorff distance to the given interval {@code y}. For two S1Intervals x and y,
   * this distance is defined by {@code h(x, y) = max_{p in x} min_{q in y} d(p, q)}, where {@code
   * d(.,.)} is measured along S1.
   */
  public double getDirectedHausdorffDistance(final S1Interval y) {
    if (y.contains(this)) {
      return 0.0; // this includes the case this is empty
    }
    if (y.isEmpty()) {
      return PI; // maximum possible distance on S1
    }

    double yComplementCenter = y.getComplementCenter();
    if (contains(yComplementCenter)) {
      return positiveDistance(y.hi(), yComplementCenter);
    } else {
      // The Hausdorff distance is realized by either two hi() endpoints or two lo() endpoints,
      // whichever is farther apart.
      double hiHi =
          new S1Interval(y.hi(), yComplementCenter).contains(hi())
              ? positiveDistance(y.hi(), hi())
              : 0;
      double loLo =
          new S1Interval(yComplementCenter, y.lo()).contains(lo())
              ? positiveDistance(lo(), y.lo())
              : 0;
      // assert hiHi > 0 || loLo > 0;
      return max(hiHi, loLo);
    }
  }

  /**
   * Expands the interval by the minimum amount necessary so that it contains the point {@code p}
   * (an angle in the range [-Pi, Pi]).
   */
  @CheckReturnValue
  public S1Interval addPoint(double p) {
    // assert abs(p) <= PI;
    if (p == -PI) {
      p = PI;
    }

    if (fastContains(p)) {
      return new S1Interval(this);
    }

    if (isEmpty()) {
      return S1Interval.fromPoint(p);
    } else {
      // Compute distance from p to each endpoint.
      double dlo = positiveDistance(p, lo);
      double dhi = positiveDistance(hi, p);
      if (dlo < dhi) {
        return new S1Interval(p, hi);
      } else {
        return new S1Interval(lo, p);
      }
      // Adding a point can never turn a non-full interval into a full one.
    }
  }

  /**
   * Returns the closest point in the interval to the point {@code p}. The interval must be
   * non-empty.
   */
  public double clampPoint(double p) {
    // assert !isEmpty();
    // assert abs(p) <= PI;
    if (p == -PI) {
      p = PI;
    }

    if (fastContains(p)) {
      return p;
    }

    // Compute distance from p to each endpoint.
    double dlo = positiveDistance(p, lo);
    double dhi = positiveDistance(hi, p);
    return (dlo < dhi) ? lo : hi;
  }

  /**
   * Returns a new interval that has been expanded on each side by the distance {@code margin}. If
   * "margin" is negative, then shrink the interval on each side by "margin" instead. The resulting
   * interval may be empty or full. Any expansion (positive or negative) of a full interval remains
   * full, and any expansion of an empty interval remains empty.
   */
  @CheckReturnValue
  public S1Interval expanded(double margin) {
    S1Interval copy = new S1Interval(this);
    copy.expandedInternal(margin);
    return copy;
  }

  /**
   * Expands this interval on each side by the distance {@code margin}. If "margin" is negative,
   * then shrink the interval on each side by "margin" instead. The resulting interval may be empty
   * or full. Any expansion (positive or negative) of a full interval remains full, and any
   * expansion of an empty interval remains empty.
   *
   * <p>Package private since only S2 code should be mutating S1Intervals for now.
   */
  void expandedInternal(double margin) {
    if (margin >= 0) {
      if (isEmpty()) {
        return;
      }
      // Check whether this interval will be full after expansion, allowing for a 1-bit rounding
      // error when computing each endpoint.
      if (getLength() + 2 * margin + 2 * DBL_EPSILON >= 2 * PI) {
        setFull();
        return;
      }
    } else {
      if (isFull()) {
        return;
      }
      // Check whether this interval will be empty after expansion, allowing for a 1-bit rounding
      // error when computing each endpoint.
      if (getLength() + 2 * margin - 2 * DBL_EPSILON <= 0) {
        setEmpty();
        return;
      }
    }

    set(
        Platform.IEEEremainder(lo - margin, 2 * PI),
        Platform.IEEEremainder(hi + margin, 2 * PI),
        false);
    if (lo <= -PI) {
      lo = PI;
    }
  }

  /** Returns the smallest interval that contains this interval and the interval {@code y}. */
  @CheckReturnValue
  public S1Interval union(S1Interval y) {
    S1Interval result = new S1Interval(this);
    result.unionInternal(y);
    return result;
  }

  /**
   * Sets this interval to the union of the current interval and {@code y}.
   *
   * <p>Package private since only S2 classes are intended to mutate S1Intervals for now.
   */
  void unionInternal(S1Interval y) {
    // The y.isFull() case is handled correctly in all cases by the code below, but can follow three
    // separate code paths depending on whether this interval is inverted, is non-inverted but
    // contains Pi, or neither.
    if (!y.isEmpty()) {
      if (fastContains(y.lo)) {
        if (fastContains(y.hi)) {
          // Either this interval contains y, or the union of the two intervals is the full
          // interval.
          if (!contains(y)) {
            setFull();
          }
        } else {
          hi = y.hi;
        }
      } else if (fastContains(y.hi)) {
        lo = y.lo;
      } else if (isEmpty() || y.fastContains(lo)) {
        // This interval contains neither endpoint of y. This means that either y contains all of
        // this interval, or the two intervals are disjoint.
        lo = y.lo;
        hi = y.hi;
      } else {
        // Check which pair of endpoints are closer together.
        double dlo = positiveDistance(y.hi, lo);
        double dhi = positiveDistance(hi, y.lo);
        if (dlo < dhi) {
          lo = y.lo;
        } else {
          hi = y.hi;
        }
      }
    }
  }

  /**
   * Returns the value of the given endpoint in this interval, which must be 0 for the low end, or 1
   * for the high end.
   */
  public double get(int endpoint) {
    if (endpoint < 0 || endpoint > 1) {
      throw new ArrayIndexOutOfBoundsException();
    }
    return endpoint == 0 ? lo : hi;
  }

  /**
   * Returns the smallest interval that contains the intersection of this interval with {@code y}.
   * Note that the region of intersection may consist of two disjoint intervals.
   */
  @CheckReturnValue
  public S1Interval intersection(final S1Interval y) {
    S1Interval result = new S1Interval(this);
    result.intersectionInternal(y);
    return result;
  }

  /**
   * Sets this interval to the intersection of the current interval and {@code y}.
   *
   * <p>Package private since only S2 classes are intended to mutate S1Intervals for now.
   */
  void intersectionInternal(final S1Interval y) {
    // The y.isFull() case is handled correctly in all cases by the code below, but can follow three
    // separate code paths depending on whether this interval is inverted, is non-inverted but
    // contains Pi, or neither.

    if (y.isEmpty()) {
      this.setEmpty();
    } else if (fastContains(y.lo)) {
      if (fastContains(y.hi)) {
        // Either this interval contains y, or the region of intersection consists of two disjoint
        // subintervals. In either case, we want to set the interval to the shorter of the two
        // original intervals.
        if (y.getLength() < getLength()) {
          this.set(y.lo, y.hi, true); // isFull() code path
        }
      } else {
        this.set(y.lo, hi, true);
      }
    } else if (fastContains(y.hi)) {
      this.set(lo, y.hi, true);
    } else {
      // This interval contains neither endpoint of y. This means that either y contains all of this
      // interval, or the two intervals are disjoint.
      if (!y.fastContains(lo)) {
        // assert !intersects(y);
        this.setEmpty();
      }
    }
  }

  /**
   * Returns true if this interval can be transformed into the interval {@code y} by moving each
   * endpoint by at most "maxError" (and without the endpoints crossing, which would invert the
   * interval). Empty and full intervals are considered to start at an arbitrary point on the unit
   * circle, thus any interval with {@code length <= 2*maxError} matches the empty interval, and any
   * interval with {@code length >= 2*Pi - 2*maxError} matches the full interval.
   */
  @JsMethod(name = "approxEqualsWithMaxError")
  public boolean approxEquals(S1Interval y, double maxError) {
    // Full and empty intervals require special cases because the "endpoints" are considered to be
    // positioned arbitrarily.
    if (isEmpty()) {
      return y.getLength() <= 2 * maxError;
    }
    if (y.isEmpty()) {
      return getLength() <= 2 * maxError;
    }
    if (isFull()) {
      return y.getLength() >= 2 * (PI - maxError);
    }
    if (y.isFull()) {
      return getLength() >= 2 * (PI - maxError);
    }

    // The purpose of the last test below is to verify that moving the endpoints does not invert the
    // interval, e.g. [-1e20, 1e20] vs. [1e20, -1e20].
    return (abs(Platform.IEEEremainder(y.lo - lo, 2 * PI)) <= maxError
        && abs(Platform.IEEEremainder(y.hi - hi, 2 * PI)) <= maxError
        && abs(getLength() - y.getLength()) <= 2 * maxError);
  }

  /** As {@link #approxEquals(S1Interval, double)}, with a default maxError of 1e-15. */
  public boolean approxEquals(final S1Interval y) {
    return approxEquals(y, 1e-15);
  }

  /** Returns true if two intervals contains the same set of points. */
  @Override
  public boolean equals(Object that) {
    if (that instanceof S1Interval) {
      S1Interval thatInterval = (S1Interval) that;
      return lo == thatInterval.lo && hi == thatInterval.hi;
    }
    return false;
  }

  @Override
  public int hashCode() {
    long value = 17;
    value = 37 * value + Double.doubleToLongBits(lo);
    value = 37 * value + Double.doubleToLongBits(hi);
    return (int) ((value >>> 32) ^ value);
  }

  @Override
  public String toString() {
    return "[" + this.lo + ", " + this.hi + "]";
  }

  /**
   * Computes the distance from {@code a} to {@code b} in the range [0, 2*Pi). This is equivalent to
   * {@code drem(b - a - PI, 2 * PI) + PI}, except that it is more numerically stable (it does not
   * lose precision for very small positive distances).
   */
  public static double positiveDistance(double a, double b) {
    double d = b - a;
    if (d >= 0) {
      return d;
    }
    // We want to ensure that if b == Pi and a == (-Pi + eps), the return result is approximately
    // 2*Pi and not zero.
    return (b + PI) - (a - PI);
  }
}
