/*
 * Copyright 2018 Google Inc.
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
import static java.lang.Math.max;

import java.math.BigDecimal;
import java.util.Objects;

/** A point consisting of BigDecimal coordinates. */
final class BigPoint implements Comparable<BigPoint> {
  final BigDecimal x;
  final BigDecimal y;
  final BigDecimal z;

  /** Creates a point of BigDecimal coordinates from a point of double coordinates. */
  BigPoint(S2Point p) {
    this(Platform.newBigDecimal(p.x), Platform.newBigDecimal(p.y), Platform.newBigDecimal(p.z));
  }

  /** Creates a point from the given BigDecimal coordinates. */
  BigPoint(BigDecimal x, BigDecimal y, BigDecimal z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  /** Subtracts 'p' from this. */
  public BigPoint sub(BigPoint p) {
    return new BigPoint(x.subtract(p.x), y.subtract(p.y), z.subtract(p.z));
  }

  /** Returns the negative of this point. */
  public BigPoint neg() {
    return new BigPoint(x.negate(), y.negate(), z.negate());
  }

  /**
   * Returns an S2Point by rounding 'this' to double precision.
   *
   * <p>This can cause non-zero points to underflow to zero. See {@link #normalizedPoint()} for an
   * implementation that avoids this underflow.
   */
  S2Point toS2Point() {
    return new S2Point(x.doubleValue(), y.doubleValue(), z.doubleValue());
  }

  /**
   * Returns a normalized S2Point by rounding 'this' to double precision, taking care to avoid
   * underflow in case all components are very small.
   */
  S2Point normalizedPoint() {
    // The fastest path is to convert all BigDecimals to doubles and normalize. This is safe as long
    // as *at least one* component is large enough to avoid underflow when squared. For IEEE 754
    // doubles, where 2^-1022 is the smallest positive normal number, this requires that the
    // absolute value of the component is at least 2^-511.
    double dx = x.doubleValue();
    double dy = y.doubleValue();
    double dz = z.doubleValue();
    final double safeMinToSquare = 0x1.0p-511; // compile time representation of 2^-511
    if (max(abs(dx), max(abs(dy), abs(dz))) >= safeMinToSquare) {
      return new S2Point(dx, dy, dz).normalize();
    }

    BigDecimal maxAbs = x.abs().max(y.abs()).max(z.abs());
    if (maxAbs.signum() == 0) {
      return S2Point.ZERO;
    }

    // If all components of {@code this} are too small, they will underflow when squared, resulting
    // in an incorrectly normalized S2Point. We uniformly scale the components so that the largest
    // component is in the range [0.1, 1).
    int numLeadingZeroesAfterDecimal = maxAbs.scale() - maxAbs.precision();
    return new S2Point(
            x.scaleByPowerOfTen(numLeadingZeroesAfterDecimal).doubleValue(),
            y.scaleByPowerOfTen(numLeadingZeroesAfterDecimal).doubleValue(),
            z.scaleByPowerOfTen(numLeadingZeroesAfterDecimal).doubleValue())
        .normalize();
  }

  /** Returns the vector cross product of 'this' with 'that'. */
  BigPoint crossProd(BigPoint that) {
    return new BigPoint(
        y.multiply(that.z).subtract(z.multiply(that.y)),
        z.multiply(that.x).subtract(x.multiply(that.z)),
        x.multiply(that.y).subtract(y.multiply(that.x)));
  }

  /** Returns the vector dot product of 'this' with 'that'. */
  BigDecimal dotProd(BigPoint that) {
    return x.multiply(that.x).add(y.multiply(that.y)).add(z.multiply(that.z));
  }

  /** Returns the vector dot product of 'this' with 'that'. */
  BigDecimal dotProd(S2Point that) {
    return dotProd(new BigPoint(that));
  }

  /** Returns true iff this and 'p' are exactly parallel or anti-parallel. */
  boolean isLinearlyDependent(BigPoint p) {
    BigPoint n = crossProd(p);
    return n.x.signum() == 0 && n.y.signum() == 0 && n.z.signum() == 0;
  }

  /** Returns true iff this and 'p' are exactly anti-parallel, antipodal points. */
  boolean isAntipodal(BigPoint p) {
    return isLinearlyDependent(p) && dotProd(p).signum() < 0;
  }

  /** Returns the square of the magnitude of this vector. */
  BigDecimal norm2() {
    return this.dotProd(this);
  }

  @Override
  public int compareTo(BigPoint p) {
    int result = x.compareTo(p.x);
    if (result != 0) {
      return result;
    }
    result = y.compareTo(p.y);
    if (result != 0) {
      return result;
    }
    return z.compareTo(p.z);
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof BigPoint)) {
      return false;
    }
    BigPoint thatPoint = (BigPoint) that;
    return x.equals(thatPoint.x) && y.equals(thatPoint.y) && z.equals(thatPoint.z);
  }

  @Override
  public int hashCode() {
    return Objects.hash(x, y, z);
  }
}
