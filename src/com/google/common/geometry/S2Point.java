/*
 * Copyright 2006 Google Inc.
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
 * An S2Point represents a point on the unit sphere as a 3D vector. Usually
 * points are normalized to be unit length, but some methods do not require
 * this.
 *
 */
@GwtCompatible(serializable = true)
public strictfp class S2Point implements S2Region, Comparable<S2Point>, Serializable {
  /** Origin of the coordinate system, [0,0,0]. */
  public static final S2Point ORIGIN = new S2Point(0, 0, 0);

  /** Direction of the x-axis. */
  public static final S2Point X_POS = new S2Point(1, 0, 0);

  /** Opposite direction of the x-axis. */
  public static final S2Point X_NEG = new S2Point(-1, 0, 0);

  /** Direction of the y-axis. */
  public static final S2Point Y_POS = new S2Point(0, 1, 0);

  /** Opposite direction of the y-axis. */
  public static final S2Point Y_NEG = new S2Point(0, -1, 0);

  /** Direction of the z-axis. */
  public static final S2Point Z_POS = new S2Point(0, 0, 1);

  /** Opposite direction of the z-axis. */
  public static final S2Point Z_NEG = new S2Point(0, 0, -1);

  // coordinates of the points
  final double x;
  final double y;
  final double z;

  public S2Point() {
    x = y = z = 0;
  }

  public S2Point(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public double getX() {
    return x;
  }

  public double getY() {
    return y;
  }

  public double getZ() {
    return z;
  }

  public final static S2Point minus(S2Point p1, S2Point p2) {
    return sub(p1, p2);
  }

  public final static S2Point neg(S2Point p) {
    return new S2Point(-p.x, -p.y, -p.z);
  }

  public final double norm2() {
    return x * x + y * y + z * z;
  }

  public double norm() {
    return Math.sqrt(norm2());
  }

  public static final S2Point crossProd(final S2Point p1, final S2Point p2) {
    return new S2Point(
        p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x);
  }

  public static final S2Point add(final S2Point p1, final S2Point p2) {
    return new S2Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
  }

  public static final S2Point sub(final S2Point p1, final S2Point p2) {
    return new S2Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  }

  public final double dotProd(S2Point that) {
    return this.x * that.x + this.y * that.y + this.z * that.z;
  }

  /**
   * Returns the scalar triple product,
   * {@code a.dotProd(S2Point.crossProd(b, c))}. This is a faster
   * implementation than calling the dotProd and crossProd methods directly.
   */
  public static final double scalarTripleProduct(S2Point a, S2Point b, S2Point c) {
    double x = b.y * c.z - b.z * c.y;
    double y = b.z * c.x - b.x * c.z;
    double z = b.x * c.y - b.y * c.x;
    double result = a.x * x + a.y * y + a.z * z;
    // assert result == a.dotProd(S2Point.crossProd(b, c));
    return result;
  }

  /**
   * Returns the distance in 3D coordinates from this to that. Equivalent to
   * {@code S2Point.sub(a, b).norm()}, but significantly faster.
   */
  public double getDistance(S2Point that) {
    return Math.sqrt(getDistance2(that));
  }

  /**
   * Returns the square of the distance in 3D coordinates from this to that. Equivalent to
   * {@code S2Point.sub(a, b).norm2()}, but significantly faster.
   */
  public double getDistance2(S2Point that) {
    double dx = this.x - that.x;
    double dy = this.y - that.y;
    double dz = this.z - that.z;
    return dx * dx + dy * dy + dz * dz;
  }

  public static final S2Point mul(final S2Point p, double m) {
    return new S2Point(m * p.x, m * p.y, m * p.z);
  }

  public static final S2Point div(final S2Point p, double m) {
    return new S2Point(p.x / m, p.y / m, p.z / m);
  }

  /** return a vector orthogonal to this one */
  @CheckReturnValue
  public final S2Point ortho() {
    int k = largestAbsComponent();
    S2Point temp;
    if (k == 1) {
      temp = X_POS;
    } else if (k == 2) {
      temp = Y_POS;
    } else {
      temp = Z_POS;
    }
    return S2Point.normalize(crossProd(this, temp));
  }

  /** Return the index of the largest component fabs */
  public final int largestAbsComponent() {
    final double absX = Math.abs(x);
    final double absY = Math.abs(y);
    final double absZ = Math.abs(z);
    if (absX > absY) {
      if (absX > absZ) {
        return 0;
      } else {
        return 2;
      }
    } else {
      if (absY > absZ) {
        return 1;
      } else {
        return 2;
      }
    }
  }

  public static final S2Point fabs(S2Point p) {
    return new S2Point(Math.abs(p.x), Math.abs(p.y), Math.abs(p.z));
  }

  public static final S2Point normalize(S2Point p) {
    double norm = p.norm();
    if (norm != 0) {
      norm = 1.0 / norm;
    }
    return S2Point.mul(p, norm);
  }

  public final double get(int axis) {
    return (axis == 0) ? x : (axis == 1) ? y : z;
  }

  /**
   * Returns the norm of the cross product,
   * {@code S2Point.crossProd(this, va).norm()}. This is more efficient
   * than calling crossProd() followed by norm().
   */
  public final double crossProdNorm(S2Point va) {
    double x = this.y * va.z - this.z * va.y;
    double y = this.z * va.x - this.x * va.z;
    double z = this.x * va.y - this.y * va.x;
    double result = Math.sqrt(x * x + y * y + z * z);
    // assert result == S2Point.crossProd(this, va).norm();
    return result;
  }

  /**
   * Rotates this point around an arbitrary axis.
   *
   * @param axis point around which rotation should be performed.
   * @param radians radians to rotate the point counterclockwise around the given axis.
   */
  @CheckReturnValue
  public S2Point rotate(S2Point axis, double radians) {
    S2Point point = S2Point.normalize(this);
    S2Point normAxis = S2Point.normalize(axis);
    S2Point pointOnAxis = S2Point.mul(normAxis, point.dotProd(normAxis));
    S2Point axisToPoint = S2Point.sub(point, pointOnAxis);
    S2Point axisToPointNormal = S2Point.crossProd(normAxis, axisToPoint);
    axisToPoint = S2Point.mul(axisToPoint, Math.cos(radians));
    axisToPointNormal = S2Point.mul(axisToPointNormal, Math.sin(radians));
    return S2Point.add(S2Point.add(axisToPoint, axisToPointNormal), pointOnAxis);
  }

  /** Return the angle between two vectors in radians */
  public final double angle(S2Point va) {
    return Math.atan2(this.crossProdNorm(va), this.dotProd(va));
  }

  /**
   * Compare two vectors, return true if all their components are within a
   * difference of margin.
   */
  boolean aequal(S2Point that, double margin) {
    return (Math.abs(x - that.x) < margin) && (Math.abs(y - that.y) < margin)
        && (Math.abs(z - that.z) < margin);
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof S2Point)) {
      return false;
    }
    S2Point thatPoint = (S2Point) that;
    return this.x == thatPoint.x && this.y == thatPoint.y && this.z == thatPoint.z;
  }

  /**
   * Returns true if this point is equal to {@code that}. Slightly faster than
   * {@link #equals(Object)}.
   */
  public boolean equalsPoint(S2Point that) {
    return this.x == that.x && this.y == that.y && this.z == that.z;
  }

  public boolean lessThan(S2Point vb) {
    if (x < vb.x) {
      return true;
    }
    if (vb.x < x) {
      return false;
    }
    if (y < vb.y) {
      return true;
    }
    if (vb.y < y) {
      return false;
    }
    if (z < vb.z) {
      return true;
    }
    return false;
  }

  // Required for Comparable
  @Override
  public int compareTo(S2Point other) {
    return (lessThan(other) ? -1 : (equalsPoint(other) ? 0 : 1));
  }

  @Override
  public String toString() {
    return "(" + x + ", " + y + ", " + z + ")";
  }

  public String toDegreesString() {
    S2LatLng s2LatLng = new S2LatLng(this);
    return "(" + Double.toString(s2LatLng.latDegrees()) + ", "
        + Double.toString(s2LatLng.lngDegrees()) + ")";
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
    value += 37 * value + Double.doubleToLongBits(Math.abs(z));
    return (int) (value ^ (value >>> 32));
  }

  // S2Region implementation.
  @Override
  public boolean contains(S2Cell cell) {
    return false;
  }

  @Override
  public S2Cap getCapBound() {
    return S2Cap.fromAxisHeight(this, 0);
  }

  @Override
  public S2LatLngRect getRectBound() {
    S2LatLng latLng = new S2LatLng(this);
    return S2LatLngRect.fromPoint(latLng);
  }

  @Override
  public boolean mayIntersect(S2Cell cell) {
    return cell.contains(this);
  }
}
