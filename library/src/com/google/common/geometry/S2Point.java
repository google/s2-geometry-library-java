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

import static java.lang.Math.abs;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import com.google.common.base.Preconditions;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.List;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * An S2Point represents a point on the unit sphere as a 3D vector. Usually points are normalized to
 * be unit length, but some methods do not require this.
 *
 * @author danieldanciu@google.com (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@SuppressWarnings("AmbiguousMethodReference")
@JsType
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

  /** An {@link S2Coder} of points that uses {@link #encode} and {@link #decode}. */
  public static final S2Coder<S2Point> CODER =
      new S2Coder<S2Point>() {
        @Override
        @JsIgnore // OutputStream
        public void encode(S2Point value, OutputStream output) throws IOException {
          value.encode(output);
        }

        @Override
        public S2Point decode(Bytes data, Cursor cursor) throws IOException {
          return S2Point.decode(data.toInputStream(cursor));
        }

        @Override public boolean isLazy() {
          return false;
        }
      };

  // Coordinates of the point.
  final double x;
  final double y;
  final double z;

  /** Constructs an S2Point with all coordinates zero, i.e. at ORIGIN. */
  @JsIgnore
  public S2Point() {
    this(0, 0, 0);
  }

  /** Constructs an S2Point from the given coordinates. */
  @JsConstructor
  public S2Point(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  /** Constructs an S2Point from the given vector of doubles in order x, y, z. */
  @JsIgnore
  public S2Point(double[] vec) {
    this(vec[0], vec[1], vec[2]);
  }

  /** Fills the x, y, and z values of this S2Point in to the first three doubles in 'vec'. */
  public void fill(double[] vec) {
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
  }

  /** Returns the value of one of the components by axis index. The axes are ordered x, y, z. */
  public final double get(int axis) {
    return (axis == 0) ? x : (axis == 1) ? y : z;
  }

  /** Returns the 'x' component. */
  public double getX() {
    return x;
  }

  /** Returns the 'y' component. */
  public double getY() {
    return y;
  }

  /** Returns the 'z' component. */
  public double getZ() {
    return z;
  }

  /** Returns add(this,p). */
  public S2Point add(S2Point p) {
    return add(this, p);
  }

  /** Returns the component-wise addition of 'p1' and 'p2'. */
  @JsIgnore
  public static final S2Point add(final S2Point p1, final S2Point p2) {
    return new S2Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
  }

  /** Returns the component-wise sum of all the given points. */
  @JsIgnore
  public static final S2Point sum(S2Point... points) {
    double x = 0;
    double y = 0;
    double z = 0;
    for (S2Point p : points) {
      x += p.x;
      y += p.y;
      z += p.z;
    }
    return new S2Point(x, y, z);
  }

  /** Returns sub(this,p). */
  public S2Point sub(S2Point p) {
    return sub(this, p);
  }

  /** Returns the component-wise subtraction of 'p1' and 'p2'. */
  @JsIgnore
  public static final S2Point sub(final S2Point p1, final S2Point p2) {
    return new S2Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  }

  /** Returns sub(this,p). */
  public static final S2Point minus(S2Point p1, S2Point p2) {
    return sub(p1, p2);
  }

  /** Returns mul(this,scale). */
  public S2Point mul(double scale) {
    return S2Point.mul(this, scale);
  }

  /** Returns the component-wise multiplication of 'p' with 'm'. */
  @JsIgnore
  public static final S2Point mul(final S2Point p, double m) {
    return new S2Point(m * p.x, m * p.y, m * p.z);
  }

  /** Returns div(this,scale). */
  public S2Point div(double scale) {
    return S2Point.div(this, scale);
  }

  /** Returns the component-wise division of 'p' by 'm'. */
  @JsIgnore
  public static final S2Point div(final S2Point p, double m) {
    return new S2Point(p.x / m, p.y / m, p.z / m);
  }

  /** Returns the vector dot product of 'this' with 'that'. */
  public final double dotProd(S2Point that) {
    return this.x * that.x + this.y * that.y + this.z * that.z;
  }

  /**
   * Returns crossProd(this,p). Note that S2.robustCrossProd() may be a better choice when the two
   * vectors are nearly parallel, are equal, or one is the negation of the other.
   */
  public S2Point crossProd(S2Point p) {
    return crossProd(this, p);
  }

  /**
   * Returns the R3 vector cross product of 'p1' and 'p2'. Note that S2.robustCrossProd(p1, p2) may
   * be a better choice when the two vectors are nearly parallel, are equal, or one is the negation
   * of the other.
   */
  @JsIgnore
  public static final S2Point crossProd(final S2Point p1, final S2Point p2) {
    return new S2Point(
        p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x);
  }

  /** Returns the component-size negation of 'this'. */
  public S2Point neg() {
    return S2Point.neg(this);
  }

  /** Returns the component-wise negation of 'p', i.e. its antipodal point. */
  @JsIgnore
  public static final S2Point neg(S2Point p) {
    return new S2Point(-p.x, -p.y, -p.z);
  }

  /** Returns the component-wise absolute value of this. */
  public S2Point fabs() {
    return S2Point.fabs(this);
  }

  /** Returns the component-wise absolute point from 'p'. */
  @JsIgnore
  public static final S2Point fabs(S2Point p) {
    return new S2Point(abs(p.x), abs(p.y), abs(p.z));
  }

  /** Returns a copy of 'this' rescaled to be unit-length. */
  public S2Point normalize() {
    return S2Point.normalize(this);
  }

  /** Returns a copy of 'p' rescaled to be unit-length. */
  @JsIgnore
  public static final S2Point normalize(S2Point p) {
    double norm = p.norm();
    if (norm != 0) {
      norm = 1.0 / norm;
    }
    return S2Point.mul(p, norm);
  }

  /** Returns the vector magnitude {@code sqrt(x*x+y*y+z*z)}. */
  public double norm() {
    return sqrt(norm2());
  }

  /** Returns the square of the vector magnitude {@code x*x+y*y+z*z}. */
  public final double norm2() {
    return x * x + y * y + z * z;
  }

  /**
   * Returns the scalar triple product, {@code a.dotProd(b.crossProd(c))}.
   *
   * <p>This is a faster implementation than calling the dotProd and crossProd methods directly.
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
   * Returns the distance in 3D coordinates from this to that.
   *
   * <p>Equivalent to {@code a.sub(b).norm()}, but significantly faster.
   *
   * <p>If ordering points by angle, this is faster than {@link #norm}, and much faster than {@link
   * #angle}, but consider using {@link S1ChordAngle}.
   *
   * <p>Returns the 3D Cartesian distance (also called the slant distance) between this and that,
   * which are normal vectors on the surface of a unit sphere. If the S2Points represent points on
   * Earth, use {@link S2Earth#getDistanceMeters(S2Point, S2Point)} to get distance in meters.
   */
  public double getDistance(S2Point that) {
    return sqrt(getDistance2(that));
  }

  /**
   * Returns the square of the distance in 3D coordinates from this to that.
   *
   * <p>Equivalent to {@code getDistance(that)<sup>2</sup>}, but significantly faster.
   *
   * <p>If ordering points by angle, this is much faster than {@link #angle}, but consider using
   * {@link S1ChordAngle}.
   */
  public double getDistance2(S2Point that) {
    double dx = this.x - that.x;
    double dy = this.y - that.y;
    double dz = this.z - that.z;
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * Return a vector orthogonal to this one, on one of the X-Y, Y-Z, or X-Z planes. Consider using
   * the alternate implementation in S2.ortho(S2Point), which avoids result coordinates that are
   * zero. These may cause degeneracies and consequently poor performance in functions that handle
   * degeneracies.
   */
  public final S2Point ortho() {
    switch (largestAbsComponent()) {
      case 1:
        return crossProd(X_POS).normalize();
      case 2:
        return crossProd(Y_POS).normalize();
      default:
        return crossProd(Z_POS).normalize();
    }
  }

  /** Return the index of the largest component by absolute value. */
  public final int largestAbsComponent() {
    return largestAbsComponent(x, y, z);
  }

  /** Return the index of the largest of x, y, or z by absolute value, as 0, 1 or 2 respectively. */
  static final int largestAbsComponent(double x, double y, double z) {
    final double absX = abs(x);
    final double absY = abs(y);
    final double absZ = abs(z);
    return (absX > absY) ? ((absX > absZ) ? 0 : 2) : ((absY > absZ) ? 1 : 2);
  }

  /**
   * Returns the norm of the cross product, {@code S2Point.crossProd(this, va).norm()}. This is more
   * efficient than calling crossProd() followed by norm().
   */
  public final double crossProdNorm(S2Point va) {
    double x = this.y * va.z - this.z * va.y;
    double y = this.z * va.x - this.x * va.z;
    double z = this.x * va.y - this.y * va.x;
    double result = sqrt(x * x + y * y + z * z);
    // assert result == S2Point.crossProd(this, va).norm();
    return result;
  }

  /**
   * Rotates this point around an arbitrary axis. This point and the provided axis are not
   * required to be normalized. However, if you already have a normalized point and axis, you can
   * use the static method {@link #rotate(S2Point, S2Point, double)} below to avoid the cost of
   * normalizing again. The returned result is normalized.
   *
   * @param axis point around which rotation should be performed
   * @param radians radians to rotate the point counterclockwise around the given axis
   */
  public S2Point rotate(S2Point axis, double radians) {
    return rotate(normalize(), axis.normalize(), radians);
  }

  /**
   * Rotates the given point around an arbitrary axis. The point and the provided axis must be
   * normalized. The result is normalized.
   *
   * @param point which should be rotated
   * @param axis point around which rotation should be performed
   * @param radians radians to rotate the point counterclockwise around the given axis
   */
  public static S2Point rotate(S2Point point, S2Point axis, double radians) {
    assert S2.isUnitLength(point);
    assert S2.isUnitLength(axis);

    S2Point center = axis.mul(point.dotProd(axis));
    S2Point axisToCenter = point.sub(center);

    // Explicitly normalize the result because there are cases where the accumulated error is
    // a bit larger than the tolerance of isUnitLength().
    return S2Point.sum(
        axisToCenter.mul(cos(radians)),
        axis.crossProd(point).mul(sin(radians)),
        center)
      .normalize();
  }

  /**
   * Returns the angle between two vectors in radians. Note that this technique is fast but can lose
   * about half its precision when the two vectors are nearly (anti-)parallel. If that is a concern,
   * use the {@link S2#stableAngle(S2Point, S2Point)} method instead, which is somewhat slower, but
   * more numerically stable. See the comments there for more information.
   */
  public final double angle(S2Point va) {
    return atan2(crossProdNorm(va), dotProd(va));
  }

  /**
   * Compare two vectors, return true if all their components are within a difference of margin. See
   * also {@link S2#approxEquals(S2Point, S2Point, double)} and its variants that check the angle
   * between two points.
   */
  boolean aequal(S2Point that, double margin) {
    return (abs(x - that.x) < margin) && (abs(y - that.y) < margin) && (abs(z - that.z) < margin);
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
   * Returns true if this point is equal to {@code that}. Slightly faster than {@link
   * #equals(Object)}.
   */
  public boolean equalsPoint(S2Point that) {
    return this.x == that.x && this.y == that.y && this.z == that.z;
  }

  /**
   * Returns true if this point is less than the given point 'vb'. Comparison is first by X,
   * breaking ties in X with Y, and then breaking ties in X and Y with Z.
   */
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

  /** Returns a string withe the representation of this point as a lat,lng in degrees. */
  public String toDegreesString() {
    S2LatLng s2LatLng = new S2LatLng(this);
    return "(" + s2LatLng.latDegrees() + ", " + s2LatLng.lngDegrees() + ")";
  }

  /** Returns a new Builder initialized to a copy of this point. */
  public Builder toBuilder() {
    return new Builder().add(this);
  }

  /**
   * Calculates hashcode based on stored coordinates. Since we want +0.0 and -0.0 to be treated the
   * same, we ignore the sign of the coordinates.
   */
  @Override
  public int hashCode() {
    long value = 17;
    value += 37 * value + Double.doubleToLongBits(abs(x));
    value += 37 * value + Double.doubleToLongBits(abs(y));
    value += 37 * value + Double.doubleToLongBits(abs(z));
    return (int) (value ^ (value >>> 32));
  }

  // S2Region implementation.
  @Override
  @JsMethod(name = "containsCell")
  public boolean contains(S2Cell cell) {
    return false;
  }

  @Override
  public boolean contains(S2Point other) {
    return equalsPoint(other);
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

  /** Writes this point to the given output stream. */
  @JsIgnore
  public void encode(OutputStream os) throws IOException {
    encode(new LittleEndianOutput(os));
  }

  /** Writes this point to the given little endian output stream. */
  @JsIgnore
  void encode(LittleEndianOutput os) throws IOException {
    os.writeDouble(x);
    os.writeDouble(y);
    os.writeDouble(z);
  }

  /** Returns a new S2Point decoded from the given input stream. */
  @JsIgnore
  public static S2Point decode(InputStream is) throws IOException {
    return decode(new LittleEndianInput(is));
  }

  /** Returns a new S2Point decoded from the given little endian input stream. */
  @JsIgnore
  static S2Point decode(LittleEndianInput is) throws IOException {
    return new S2Point(is.readDouble(), is.readDouble(), is.readDouble());
  }

  /**
   * An S2Shape representing a list of S2Points. Each point is represented as a degenerate edge with
   * the same starting and ending vertices.
   *
   * <p>This class is useful for adding a collection of points to an S2ShapeIndex.
   */
  public abstract static class Shape extends AbstractList<S2Point>
      implements S2Shape, Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * A fast {@link S2Coder} of a shape of multiple points. Decodes points on demand, i.e.
     * {@link S2Coder#isLazy()} is true.
     */
    public static final S2Coder<Shape> FAST_CODER =
        S2PointVectorCoder.FAST.delegating(s -> s, Shape::fromList);

    /**
     * A compact {@link S2Coder} of a shape of multiple points.  Decodes points on demand, i.e.
     * {@link S2Coder#isLazy()} is true.
     */
    public static final S2Coder<Shape> COMPACT_CODER =
        S2PointVectorCoder.COMPACT.delegating(s -> s, Shape::fromList);

    /** Returns an S2Shape view of the given single point. */
    public static Shape singleton(final S2Point point) {
      return new Shape() {
        private static final long serialVersionUID = 1L;

        @Override
        public int size() {
          return 1;
        }

        @Override
        public S2Point get(int index) {
          if (index != 0) {
            throw new IndexOutOfBoundsException();
          }
          return point;
        }
      };
    }

    /** Returns an S2Shape view of the given list. The points are *not* copied. */
    public static Shape fromList(final List<S2Point> points) {
      return new Shape() {
        private static final long serialVersionUID = 1L;

        @Override
        public int size() {
          return points.size();
        }

        @Override
        public S2Point get(int index) {
          return points.get(index);
        }
      };
    }

    @Override
    public boolean hasInterior() {
      return false;
    }

    @Override
    public boolean containsOrigin() {
      return false;
    }

    @Override
    public int numEdges() {
      return size();
    }

    @Override
    public void getEdge(int index, MutableEdge result) {
      result.a = result.b = get(index);
    }

    @Override
    public int numChains() {
      return size();
    }

    @Override
    public int getChainStart(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      return chainId;
    }

    @Override
    public int getChainLength(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      return 1;
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge result) {
      result.a = result.b = getChainVertex(chainId, offset);
    }

    @Override
    public S2Point getChainVertex(int chainId, int edgeOffset) {
      Preconditions.checkElementIndex(edgeOffset, 2);
      return get(chainId);
    }

    @Override
    public void getChainPosition(int edgeId, ChainPosition result) {
      // Each edge is in its own chain.
      Preconditions.checkElementIndex(edgeId, numEdges());
      result.set(edgeId, 0);
    }

    @Override
    public int dimension() {
      return 0;
    }

    /** An encoder/decoder of {@link S2Point.Shape}s. */
    public static class Coder implements S2Coder<S2Point.Shape> {
      /**
       * An instance of {@link Coder} which encodes/decodes {@link S2Point.Shape}s in the {@code
       * FAST} format. Decodes on demand, i.e. {@link S2Coder#isLazy()} is true.
       */
      static final Coder FAST = new Coder(S2PointVectorCoder.FAST);

      /**
       * An instance of {@link Coder} which encodes/decodes {@link S2Point.Shape}s in the {@code
       * COMPACT} format.  Decodes on demand, i.e. {@link S2Coder#isLazy()} is true.
       */
      static final Coder COMPACT = new Coder(S2PointVectorCoder.COMPACT);

      private final S2PointVectorCoder coder;

      private Coder(S2PointVectorCoder coder) {
        this.coder = coder;
      }

      @Override
      public void encode(Shape shape, OutputStream output) throws IOException {
        coder.encode(shape, output);
      }

      @Override
      public Shape decode(Bytes data, Cursor cursor) throws IOException {
        return S2Point.Shape.fromList(coder.decode(data, cursor));
      }

      @Override
      public boolean isLazy() {
        return coder.isLazy();
      }
    }
  }

  /** A builder of {@link S2Point} instances. */
  @JsType
  public static final class Builder {
    private double x;
    private double y;
    private double z;

    /** Constructs a new builder initialized to {@link #ORIGIN}. */
    public Builder() {}

    /** Adds point. */
    @CanIgnoreReturnValue
    public Builder add(S2Point point) {
      x += point.x;
      y += point.y;
      z += point.z;
      return this;
    }

    /** Returns a new {@link S2Point} copied from the current state of this builder. */
    public S2Point build() {
      return new S2Point(x, y, z);
    }
  }
}
