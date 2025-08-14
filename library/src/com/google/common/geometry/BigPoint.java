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

  /** Returns an S2Point by rounding 'this' to double precision. */
  S2Point toS2Point() {
    return new S2Point(x.doubleValue(), y.doubleValue(), z.doubleValue());
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
