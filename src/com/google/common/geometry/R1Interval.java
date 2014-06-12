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

import java.io.Serializable;

import javax.annotation.CheckReturnValue;

/**
 * An R1Interval represents a closed, bounded interval on the real line. It is capable of
 * representing the empty interval (containing no points) and zero-length intervals (containing a
 * single point).
 *
 */
@GwtCompatible(serializable = true)
public final strictfp class R1Interval implements Serializable {
  private double lo;
  private double hi;

  /**
   * Default constructor, contains the empty interval.
   *
   * <p>Package private since only the S2 library needs to mutate R1Intervals. External code that
   * needs an empty interval should call {@link #empty()}.
   */
  R1Interval() {
    lo = 1;
    hi = 0;
  }

  /** Interval constructor. If lo > hi, the interval is empty. */
  public R1Interval(double lo, double hi) {
    this.lo = lo;
    this.hi = hi;
  }

  /** Copy constructor. */
  public R1Interval(R1Interval interval) {
    this.lo = interval.lo;
    this.hi = interval.hi;
  }

  /** Returns an empty interval. (Any interval where lo > hi is considered empty.) */
  public static R1Interval empty() {
    return new R1Interval(1, 0);
  }

  /** Convenience method to construct an interval containing a single point. */
  public static R1Interval fromPoint(double p) {
    return new R1Interval(p, p);
  }

  /**
   * Convenience method to construct the minimal interval containing the two given points. This is
   * equivalent to starting with an empty interval and calling AddPoint() twice, but it is more
   * efficient.
   */
  public static R1Interval fromPointPair(double p1, double p2) {
    R1Interval result = new R1Interval();
    result.initFromPointPair(p1, p2);
    return result;
  }

  void initFromPointPair(double p1, double p2) {
    if (p1 <= p2) {
      lo = p1;
      hi = p2;
    } else {
      lo = p2;
      hi = p1;
    }
  }

  public double lo() {
    return lo;
  }

  public double hi() {
    return hi;
  }

  /** Designates which end of the interval to work with. */
  enum Endpoint {
    /** The low end of the interval. */
    LO {
      @Override
      public double getValue(R1Interval interval) {
        return interval.lo;
      }
      @Override
      public void setValue(R1Interval interval, double value) {
        interval.lo = value;
      }
      @Override
      public Endpoint opposite() {
        return HI;
      }
    },
    /** The high end of the interval. */
    HI {
      @Override
      public double getValue(R1Interval interval) {
        return interval.hi;
      }
      @Override
      public void setValue(R1Interval interval, double value) {
        interval.hi = value;
      }
      @Override
      public Endpoint opposite() {
        return LO;
      }
    };
    public abstract double getValue(R1Interval interval);
    public abstract void setValue(R1Interval interval, double value);
    public abstract Endpoint opposite();
  }

  /** Returns the value at the given Endpoint, which must not be null. */
  double getValue(Endpoint endpoint) {
    return endpoint.getValue(this);
  }

  /** Sets the value of the given Endpoint, which must not be null. */
  void setValue(Endpoint endpoint, double value) {
    endpoint.setValue(this, value);
  }

  /** Returns true if the interval is empty, i.e. it contains no points. */
  public boolean isEmpty() {
    return lo > hi;
  }

  /** Returns the center of the interval. For empty intervals, the result is arbitrary. */
  public double getCenter() {
    return 0.5 * (lo + hi);
  }

  /** Returns the length of the interval. The length of an empty interval is negative. */
  public double getLength() {
    return hi - lo;
  }

  public boolean contains(double p) {
    return p >= lo && p <= hi;
  }

  public boolean interiorContains(double p) {
    return p > lo && p < hi;
  }

  /** Returns true if this interval contains the interval {@code y}. */
  public boolean contains(R1Interval y) {
    if (y.isEmpty()) {
      return true;
    }
    return y.lo >= lo && y.hi <= hi;
  }

  /**
   * Returns true if the interior of this interval contains the entire interval {@code y}
   * (including its boundary).
   */
  public boolean interiorContains(R1Interval y) {
    if (y.isEmpty()) {
      return true;
    }
    return y.lo > lo && y.hi < hi;
  }

  /**
   * Returns true if this interval intersects {@code y}, i.e. if they have any points in
   * common.
   */
  public boolean intersects(R1Interval y) {
    if (lo <= y.lo) {
      return y.lo <= hi && y.lo <= y.hi;
    } else {
      return lo <= y.hi && lo <= hi;
    }
  }

  /**
   * Returns true if the interior of this interval intersects any point of {@code y} (including its
   * boundary).
   */
  public boolean interiorIntersects(R1Interval y) {
    return y.lo < hi && lo < y.hi && lo < hi && y.lo <= y.hi;
  }

  /**
   * Sets the minimum and maximum value of this interval. If {@code lo} is greater than {@code hi}
   * this interval will become empty.
   *
   * <p>Package private since only the S2 libraries have a current need to mutate R1Intervals.
   */
  void set(double lo, double hi) {
    this.lo = lo;
    this.hi = hi;
  }

  /**
   * Sets the minimum value of this interval. If {@code lo} is greater than {@code hi()} this
   * interval will become empty.
   *
   * <p>Package private since only the S2 libraries have a current need to mutate R1Intervals.
   */
  void setLo(double lo) {
    this.lo = lo;
  }

  /**
   * Sets the maximum value of this interval. If {@code hi} is less than {@code lo()} this interval
   * will become empty.
   *
   * <p>Package private since only the S2 libraries have a current need to mutate R1Intervals.
   */
  void setHi(double hi) {
    this.hi = hi;
  }

  /**
   * Sets the current interval to the empty interval.
   *
   * <p>Package private since only the S2 libraries have a current need to mutate R1Intervals.
   */
  void setEmpty() {
    this.lo = 1;
    this.hi = 0;
  }

  /**
   * Expands this interval so that it contains the point {@code p}.
   *
   * <p>Package private since only the S2 library needs to mutate R1Intervals.
   */
  void unionInternal(double p) {
    if (isEmpty()) {
      lo = p;
      hi = p;
    } else if (p < lo) {
      lo = p;
    } else if (p > hi) {
      hi = p;
    }
  }

  /**
   * Returns the closest point in the interval to the point {@code p}. The interval must be
   * non-empty.
   */
  public double clampPoint(double p) {
    // assert (!isEmpty());
    return Math.max(lo, Math.min(hi, p));
  }

  /**
   * Return an interval that contains all points with a distance "radius" of a
   * point in this interval. Note that the expansion of an empty interval is
   * always empty.
   */
  @CheckReturnValue
  public R1Interval expanded(double radius) {
    // assert (radius >= 0);
    if (isEmpty()) {
      return this;
    }
    return new R1Interval(lo - radius, hi + radius);
  }

  /**
   * Expands this interval to contain all points within a distance "radius" of a point in this
   * interval.
   *
   * <p>Package private since only S2 classes are intended to mutate R1Intervals for now.
   */
  void expandedInternal(double radius) {
    lo -= radius;
    hi += radius;
  }

  /** Returns the smallest interval that contains this interval and {@code y}. */
  @CheckReturnValue
  public R1Interval union(R1Interval y) {
    if (isEmpty()) {
      return y;
    }
    if (y.isEmpty()) {
      return this;
    }
    return new R1Interval(Math.min(lo, y.lo), Math.max(hi, y.hi));
  }

  /**
   * Sets this interval to the union of this interval and {@code y}.
   *
   * <p>Package private since only S2 classes are intended to mutate R11Intervals for now.
   */
  void unionInternal(R1Interval y) {
    if (isEmpty()) {
      lo = y.lo;
      hi = y.hi;
    } else if (!y.isEmpty()) {
      lo = Math.min(lo, y.lo);
      hi = Math.max(hi, y.hi);
    }
  }

  /**
   * Returns the intersection of this interval with {@code y}. Empty intervals do not need to be
   * special-cased.
   */
  @CheckReturnValue
  public R1Interval intersection(R1Interval y) {
    return new R1Interval(Math.max(lo, y.lo), Math.min(hi, y.hi));
  }

  /**
   * Sets this interval to the intersection of the current interval and {@code y}.
   *
   * <p>Package private since only S2 classes are intended to mutate R1 intervals for now.
   */
  void intersectionInternal(R1Interval y) {
    lo = Math.max(lo, y.lo);
    hi = Math.min(hi, y.hi);
  }

  /** Returns the smallest interval that contains this interval and the point {@code p}. */
  @CheckReturnValue
  public R1Interval addPoint(double p) {
    if (isEmpty()) {
      return R1Interval.fromPoint(p);
    } else if (p < lo) {
      return new R1Interval(p, hi);
    } else if (p > hi) {
      return new R1Interval(lo, p);
    } else {
      return new R1Interval(lo, hi);
    }
  }

  @Override
  public boolean equals(Object that) {
    if (that instanceof R1Interval) {
      R1Interval y = (R1Interval) that;
      // Return true if two intervals contain the same set of points.
      return (lo == y.lo && hi == y.hi) || (isEmpty() && y.isEmpty());

    }
    return false;
  }

  @Override
  public int hashCode() {
    if (isEmpty()) {
      return 17;
    }

    long value = 17;
    value = 37 * value + Double.doubleToLongBits(lo);
    value = 37 * value + Double.doubleToLongBits(hi);
    return (int) (value ^ (value >>> 32));
  }

  /**
   * As {@link #approxEquals(R1Interval, double)}, with a default value for maxError just larger
   * than typical rounding errors in computing intervals.
   */
  public boolean approxEquals(R1Interval y) {
    return approxEquals(y, 1e-15);
  }

  /**
   * Returns true if this interval can be transformed into {@code y} by moving each endpoint by at
   * most {@code maxError}. The empty interval is considered to be positioned arbitrarily on the
   * real line, thus any interval for which {@code length <= 2*maxError} is true matches the empty
   * interval.
   */
  public boolean approxEquals(R1Interval y, double maxError) {
    if (isEmpty()) {
      return y.getLength() <= maxError;
    }
    if (y.isEmpty()) {
      return getLength() <= maxError;
    }
    return Math.abs(y.lo - lo) <= maxError && Math.abs(y.hi - hi) <= maxError;
  }

  @Override
  public String toString() {
    return "[" + lo + ", " + hi + "]";
  }
}
