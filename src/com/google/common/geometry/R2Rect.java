/*
 * Copyright 2013 Google Inc.
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
 * An R2Rect represents a closed axis-aligned rectangle in the (x,y) plane. This class is mutable
 * to allow iteratively constructing bounds via e.g. {@link #addPoint(R2Vector)}.
 *
 * <p>This class is package private, since it is not intended for external usage.
 */
final strictfp class R2Rect {
  private final R1Interval x;
  private final R1Interval y;

  /** Creates an empty R2Rect. */
  public R2Rect() {
    // The default R1Interval constructor creates an empty interval.
    this.x = new R1Interval();
    this.y = new R1Interval();
    // assert (isValid());
  }

  /** Constructs a rectangle from the given lower-left and upper-right points. */
  public R2Rect(R2Vector lo, R2Vector hi) {
    x = new R1Interval(lo.x(), hi.x());
    y = new R1Interval(lo.y(), hi.y());
    // assert (isValid());
  }

  /**
   * Constructs a rectangle from the given intervals in x and y.  The two intervals must either be
   * both empty or both non-empty.
   */
  public R2Rect(R1Interval x, R1Interval y) {
    this.x = x;
    this.y = y;
    // assert (isValid());
  }

  /**
   * Returns a new instance of the canonical empty rectangle.  Use isEmpty() to test for empty
   * rectangles, since they have more than one representation.
   */
  public static R2Rect empty() {
    return new R2Rect(R1Interval.empty(), R1Interval.empty());
  }

  /**
   * Returns a new rectangle from a center point and size in each dimension. Both components of size
   * should be non-negative, i.e. this method cannot be used to create an empty rectangle.
   */
  public static R2Rect fromCenterSize(R2Vector center, R2Vector size) {
    return new R2Rect(
        new R1Interval(center.x() - 0.5 * size.x(), center.x() + 0.5 * size.x()),
        new R1Interval(center.y() - 0.5 * size.y(), center.y() + 0.5 * size.y()));
  }

  /** Returns a rectangle containing a single point. */
  public static R2Rect fromPoint(R2Vector p) {
    return new R2Rect(p, p);
  }

  /**
   * Returns the minimal bounding rectangle containing the two given points. This is equivalent to
   * starting with an empty rectangle and calling addPoint() twice.  Note that it is different than
   * the R2Rect(lo, hi) constructor, where the first point is always used as the lower-left corner
   * of the resulting rectangle.
   */
  public static R2Rect fromPointPair(R2Vector p1, R2Vector p2) {
    return new R2Rect(R1Interval.fromPointPair(p1.x(), p2.x()),
        R1Interval.fromPointPair(p1.y(), p2.y()));
  }

  /** Returns the interval along the x-axis. */
  public R1Interval x() {
    return x;
  }

  /** Returns the interval along the y-axis. */
  public R1Interval y() {
    return y;
  }

  /** Returns the point in this rectangle with the minimum x and y values. */
  public R2Vector lo() {
    return new R2Vector(x().lo(), y().lo());
  }

  /** Returns the point in this rectangle with the maximum x and y values. */
  public R2Vector hi() {
    return new R2Vector(x().hi(), y().hi());
  }

  /**
   * Returns true if this rectangle is valid, which essentially just means that if the bound for
   * either axis is empty then both must be.
   */
  public boolean isValid() {
    // The x/y ranges must either be both empty or both non-empty.
    return x().isEmpty() == y().isEmpty();
  }

  /** Return true if this rectangle is empty, i.e. it contains no points at all. */
  public boolean isEmpty() {
    return x().isEmpty();
  }

  /**
   * Returns the k<super>th</super> vertex of this rectangle (k = 0,1,2,3) in CCW order. Vertex 0 is
   * in the lower-left corner.
   */
  public R2Vector getVertex(int k) {
    // Twiddle bits to return the points in CCW order (lower left, lower right,
    // upper right, upper left).
    // assert (k >= 0 && k <= 3);
    return getVertex((k >> 1) ^ (k & 1), k >> 1);
  }

  /**
   * Returns the vertex in direction "i" along the x-axis (0=left, 1=right) and direction "j" along
   * the y-axis (0=down, 1=up). Equivalently, returns the vertex constructed by selecting endpoint
   * "i" of the x-interval (0=lo, 1=hi) and vertex "j" of the y-interval.
   */
  public R2Vector getVertex(int i, int j) {
    return new R2Vector(i == 0 ? x.lo() : x.hi(), j == 0 ? y.lo() : y.hi());
  }

  /** Returns the center of this rectangle in (x,y)-space. */
  public R2Vector getCenter() {
    return new R2Vector(x().getCenter(), y().getCenter());
  }

  /**
   * Return the width and height of this rectangle in (x,y)-space.  Empty rectangles have a negative
   * width and height.
   */
  public R2Vector getSize() {
    return new R2Vector(x().getLength(), y().getLength());
  }

  /** Valid axes. */
  public enum Axis {
    X {
      @Override
      public R1Interval getInterval(R2Rect rect) {
        return rect.x;
      }
    },
    Y {
      @Override
      public R1Interval getInterval(R2Rect rect) {
        return rect.y;
      }
      
    };
    public abstract R1Interval getInterval(R2Rect rect);
  }

  /** Returns the interval for the given axis, which must not be null. */
  public R1Interval getInterval(Axis axis) {
    return axis.getInterval(this);
  }

  /**
   * Returns true if this rectangle contains the given point.  Note that rectangles are closed
   * regions, i.e. they contain their boundary.
   */
  public boolean contains(R2Vector p) {
    return x().contains(p.x()) && y().contains(p.y());
  }

  /**
   * Returns true if and only if the given point is contained in the interior of the region (i.e.
   * the region excluding its boundary).
   */
  public boolean interiorContains(R2Vector p) {
    return x().interiorContains(p.x()) && y().interiorContains(p.y());
  }

  /** Returns true if and only if this rectangle contains the given other rectangle. */
  public boolean contains(R2Rect other) {
    return x().contains(other.x()) && y().contains(other.y());
  }

  /**
   * Returns true if and only if the interior of this rectangle contains all points of the given
   * other rectangle (including its boundary).
   */
  public boolean interiorContains(R2Rect other) {
    return x().interiorContains(other.x()) && y().interiorContains(other.y());
  }

  /** Returns true if this rectangle and the given other rectangle have any points in common. */
  public boolean intersects(R2Rect other) {
    return x().intersects(other.x()) && y().intersects(other.y());
  }

  /**
   * Return true if and only if the interior of this rectangle intersects any point (including the
   * boundary) of the given other rectangle.
   */
  public boolean interiorIntersects(R2Rect other) {
    return x().interiorIntersects(other.x()) && y().interiorIntersects(other.y());
  }

  /**
   * Increase the size of the bounding rectangle to include the given point. This rectangle is
   * expanded by the minimum amount possible.
   */
  public void addPoint(R2Vector p) {
    x.unionInternal(p.x());
    y.unionInternal(p.y());
  }

  /**
   * Return the closest point in this rectangle to the given point "p". This rectangle must be
   * non-empty.
   */
  public R2Vector clampPoint(R2Vector p) {
    return new R2Vector(x().clampPoint(p.x()), y().clampPoint(p.y()));
  }

  /**
   * Return a rectangle that has been expanded on each side in the x-direction by margin.x(), and on
   * each side in the y-direction by margin.y().  If either margin is empty, then shrink the
   * interval on the corresponding sides instead.  The resulting rectangle may be empty.  Any
   * expansion of an empty rectangle remains empty.
   */
  public R2Rect expanded(R2Vector margin) {
    R1Interval xx = x().expanded(margin.x());
    R1Interval yy = y().expanded(margin.y());
    if (xx.isEmpty() || yy.isEmpty()) {
      return empty();
    } else {
      return new R2Rect(xx, yy);
    }
  }

  /**
   * Returns a rectangle that has been expanded on both sides by the given margin.  Any expansion of
   * an empty rectangle remains empty.
   */
  public R2Rect expanded(double margin) {
    return expanded(new R2Vector(margin, margin));
  }

  /**
   * Returns the smallest rectangle containing the union of this rectangle and the given rectangle.
   */
  public R2Rect union(R2Rect other) {
    return new R2Rect(x().union(other.x()), y().union(other.y()));
  }

  /**
   * Returns the smallest rectangle containing the intersection of this rectangle and the given
   * rectangle.
   */
  public R2Rect intersection(R2Rect other) {
    R1Interval xx = x().intersection(other.x());
    R1Interval yy = y().intersection(other.y());
    if (xx.isEmpty() || yy.isEmpty()) {
      return empty();
    }
    return new R2Rect(xx, yy);
  }

  /** Returns a simple convolution hashcodes from the x and y internals. */
  @Override
  public int hashCode() {
    return x.hashCode() * 701 + y.hashCode();
  }

  /** Returns true if two rectangles contains the same set of points. */
  @Override
  public boolean equals(Object other) {
    if (other instanceof R2Rect) {
      R2Rect r2 = (R2Rect) other;
      return x().equals(r2.x()) && y().equals(r2.y());
    } else {
      return false;
    }
  }

  /**
   * Returns true if the x- and y-intervals of the two rectangles are the same up to the given
   * tolerance.  See {@link R1Interval} for details on approximate interval equality.
   */
  public boolean approxEquals(R2Rect other) {
    return approxEquals(other, 1e-15);
  }

  /**
   * Returns true if the given rectangles are equal to within {@code maxError}. See
   * {@link R1Interval} for details on approximate interval equality.
   */
  public boolean approxEquals(R2Rect other, double maxError) {
    return x().approxEquals(other.x(), maxError) && y().approxEquals(other.y(), maxError);
  }

  /** Returns a simple string representation of this rectangle's lower and upper corners. */
  @Override
  public String toString() {
    return "[Lo" + lo() + ", Hi" + hi() + "]";
  }
}

