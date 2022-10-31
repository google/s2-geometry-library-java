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

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

import java.io.Serializable;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * R2Vector represents a vector in the two-dimensional space. It defines the basic geometrical
 * operations for 2D vectors, e.g. cross product, addition, norm, comparison, etc.
 *
 * @author danieldanciu@google.com (Daniel Danciu)
 */
@SuppressWarnings("AmbiguousMethodReference")
@JsType
public final strictfp class R2Vector implements Serializable {
  double x;
  double y;

  /** Constructs a new R2Vector at the origin [0,0] of the R2 coordinate system. */
  @JsIgnore
  public R2Vector() {
    this(0, 0);
  }

  /** Constructs a new R2 vector from the given x and y coordinates. */
  public R2Vector(double x, double y) {
    this.x = x;
    this.y = y;
  }

  /** Constructs a new R2 vector from the given coordinates array, which must have length 2. */
  @JsIgnore
  public R2Vector(double[] coord) {
    this(coord[0], coord[1]);
    if (coord.length != 2) {
      throw new IllegalStateException("Points must have exactly 2 coordinates");
    }
  }

  /** Returns the x coordinate of this R2 vector. */
  public double x() {
    return x;
  }

  /** Returns the y coordinate of this R2 vector. */
  public double y() {
    return y;
  }

  /**
   * Returns the coordinate of the given axis, which will be the x axis if index is 0, and the y
   * axis if index is 1.
   *
   * @throws ArrayIndexOutOfBoundsException Thrown if the given index is not 0 or 1.
   */
  public double get(int index) {
    if (index < 0 || index > 1) {
      throw new ArrayIndexOutOfBoundsException(index);
    }
    return index == 0 ? this.x : this.y;
  }

  /**
   * Sets the position of this vector from the given other vector. Package private since this is
   * only mutable for S2.
   */
  void set(R2Vector v) {
    this.x = v.x();
    this.y = v.y();
  }

  /**
   * Sets the position of this vector from the given values. Package private since this is only
   * mutable for S2.
   */
  @JsMethod(name = "setCoordinates")
  void set(double x, double y) {
    this.x = x;
    this.y = y;
  }

  /** Returns the vector result of {@code p1 - p2}. */
  public static R2Vector add(final R2Vector p1, final R2Vector p2) {
    return new R2Vector(p1.x + p2.x, p1.y + p2.y);
  }

  /** Returns add(this, p) */
  public R2Vector add(R2Vector p) {
    return add(this, p);
  }

  /** Returns the vector result of {@code p1 - p2}. */
  public static R2Vector sub(final R2Vector p1, final R2Vector p2) {
    return new R2Vector(p1.x - p2.x, p1.y - p2.y);
  }

  /** Returns sub(this, p) */
  public R2Vector sub(R2Vector p) {
    return sub(this, p);
  }

  /**
   * Returns the element-wise multiplication of p1 and p2, e.g. {@code vector [p1.x*p2.x,
   * p1.y*p2.y]}.
   */
  public static R2Vector mul(final R2Vector p, double m) {
    return new R2Vector(m * p.x, m * p.y);
  }

  /** Returns mul(this, m) */
  public R2Vector mul(double m) {
    return mul(this, m);
  }

  /** Returns the vector magnitude. */
  public double norm() {
    return sqrt(norm2());
  }

  /** Returns the square of the vector magnitude. */
  public double norm2() {
    return (x * x) + (y * y);
  }

  /**
   * Returns a new vector scaled to magnitude 1, or a copy of the original vector if magnitude was
   * 0.
   */
  public static R2Vector normalize(R2Vector vector) {
    double n = vector.norm();
    if (n != 0) {
      return mul(vector, 1.0 / n);
    } else {
      return new R2Vector(vector.x, vector.y);
    }
  }

  /**
   * Returns a new R2 vector orthogonal to the current one with the same norm and counterclockwise
   * to it.
   */
  public R2Vector ortho() {
    return new R2Vector(-y, x);
  }

  /** Returns the dot product of the given vectors. */
  @JsIgnore
  public static double dotProd(final R2Vector p1, final R2Vector p2) {
    return (p1.x * p2.x) + (p1.y * p2.y);
  }

  /** Returns the dot product of this vector with that vector. */
  public double dotProd(R2Vector that) {
    return dotProd(this, that);
  }

  /** Returns the cross product of this vector with that vector. */
  public double crossProd(final R2Vector that) {
    return this.x * that.y - this.y * that.x;
  }

  /**
   * Returns true if this vector is less than that vector, with the x-axis as the primary sort key
   * and the y-axis as the secondary sort key.
   */
  public boolean lessThan(R2Vector that) {
    if (x < that.x) {
      return true;
    }
    if (that.x < x) {
      return false;
    }
    if (y < that.y) {
      return true;
    }
    return false;
  }

  /** Returns true if that object is an R2Vector with exactly the same x and y coordinates. */
  @Override
  public boolean equals(Object that) {
    if (!(that instanceof R2Vector)) {
      return false;
    }
    R2Vector thatPoint = (R2Vector) that;
    return this.x == thatPoint.x && this.y == thatPoint.y;
  }

  /**
   * Calcualates hashcode based on stored coordinates. Since we want +0.0 and -0.0 to be treated the
   * same, we ignore the sign of the coordinates.
   */
  @Override
  public int hashCode() {
    long value = 17;
    value += 37 * value + Double.doubleToLongBits(abs(x));
    value += 37 * value + Double.doubleToLongBits(abs(y));
    return (int) (value ^ (value >>> 32));
  }

  @Override
  public String toString() {
    return "(" + x + ", " + y + ")";
  }
}
