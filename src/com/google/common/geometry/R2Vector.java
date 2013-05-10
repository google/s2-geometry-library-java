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

/**
 * R2Vector represents a vector in the two-dimensional space. It defines the
 * basic geometrical operations for 2D vectors, e.g. cross product, addition,
 * norm, comparison, etc.
 *
 */
public final strictfp class R2Vector {
  private final double x;
  private final double y;

  /** Constructs a new R2Vector at the origin [0,0] of the R2 coordinate system. */
  public R2Vector() {
    this(0, 0);
  }

  /** Constructs a new R2 vector from the given x and y coordinates. */
  public R2Vector(double x, double y) {
    this.x = x;
    this.y = y;
  }

  /** Constructs a new R2 vector from the given coordinates array, which must have length 2. */
  public R2Vector(double[] coord) {
    if (coord.length != 2) {
      throw new IllegalStateException("Points must have exactly 2 coordinates");
    }
    x = coord[0];
    y = coord[1];
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
    if (index > 1) {
      throw new ArrayIndexOutOfBoundsException(index);
    }
    return index == 0 ? this.x : this.y;
  }

  /** Returns the vector result of {@code p1 - p2}. */
  public static R2Vector add(final R2Vector p1, final R2Vector p2) {
    return new R2Vector(p1.x + p2.x, p1.y + p2.y);
  }

  /** Returns the vector result of {@code p1 - p2}. */
  public static R2Vector sub(final R2Vector p1, final R2Vector p2) {
    return new R2Vector(p1.x - p2.x, p1.y - p2.y);
  }

  /**
   * Returns the element-wise multiplication of p1 and p2, e.g.
   * {@code vector [p1.x*p2.x, p1.y*p2.y]}.
   */
  public static R2Vector mul(final R2Vector p, double m) {
    return new R2Vector(m * p.x, m * p.y);
  }

  /** Returns the square of the vector magnitude. */
  public double norm2() {
    return (x * x) + (y * y);
  }

  /**
   * Returns a new R2 vector orthogonal to the current one with the same norm and counterclockwise
   * to it.
   */
  public R2Vector ortho() {
    return new R2Vector(-y, x);
  }

  /** Returns the dot product of the given vectors. */
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
   * Calcualates hashcode based on stored coordinates. Since we want +0.0 and
   * -0.0 to be treated the same, we ignore the sign of the coordinates.
   */
  @Override
  public int hashCode() {
    long value = 17;
    value += 37 * value + Double.doubleToLongBits(Math.abs(x));
    value += 37 * value + Double.doubleToLongBits(Math.abs(y));
    return (int) (value ^ (value >>> 32));
  }

  @Override
  public String toString() {
    return "(" + x + ", " + y + ")";
  }
}
