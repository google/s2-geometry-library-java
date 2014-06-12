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

import javax.annotation.CheckReturnValue;

/**
 * S2LatLngRect represents a latitude-longitude rectangle. It is capable of representing the
 * empty and full rectangles as well as single points.
 *
 * <p>Note that the latitude-longitude space is considered to have a <strong>cylindrical</strong>
 * topology rather than a spherical one, i.e. the poles have multiple lat/lng representations. An
 * S2LatLngRect may be defined so that it includes some representations of a pole but not
 * others. Use the polarClosure() method if you want to expand a rectangle so that it contains all
 * possible representations of any contained poles.
 * 
 * <p>Because S2LatLngRect uses S1Interval to store the longitude range, longitudes of -180 degrees
 * are treated specially.  Except for empty and full longitude spans, -180 degree longitudes will
 * turn into +180 degrees.  This sign flip causes lng.lo() to be greater than lng.hi(), indicating
 * that the rectangle will wrap around through -180 instead of through +179.  Thus the math is
 * consistent with the library, but the sign flip can be surprising, especially when working with
 * map projections where -180 and +180 are at opposite ends of the flattened map.  See
 * {@link S1Interval} for more details.
 * 
 * <p>S2LatLngRect is immutable, so all methods that change the boundary return a new instance.  To
 * efficiently make numerous alterations to the bounds, use a {@link Builder} instead,
 * mutate it to compute the desired bounds, and then call {@link Builder#build()} to convert it to
 * an immutable S2LatLngRect.
 *
 */
@GwtCompatible(serializable = true)
public final strictfp class S2LatLngRect extends S2LatLngRectBase {
  /**
   * The canonical empty rectangle, as derived from the empty R1 and S1 intervals.
   * Empty: lat.lo = 1, lat.hi = 0, lng.lo = Pi, lng.hi = -Pi (radians).
   */
  public static S2LatLngRect empty() {
    return new S2LatLngRect(R1Interval.empty(), S1Interval.empty());
  }

  /** The canonical full rectangle. */
  public static S2LatLngRect full() {
    return new S2LatLngRect(fullLat(), fullLng());
  }

  /** The full allowable range of latitudes. */
  public static R1Interval fullLat() {
    return new R1Interval(-S2.M_PI_2, S2.M_PI_2);
  }

  /** The full allowable range of longitudes. */
  public static S1Interval fullLng() {
    return S1Interval.full();
  }

  /**
   * Constructs a rectangle of the given size centered around the given point.  {@code center} needs
   * to be normalized, but {@code size} does not. The latitude interval of the result is clamped to 
   * [-90, 90] degrees, and the longitude interval of the result is full() if and only if the
   * longitude size is 360 degrees or more.  Examples of clamping (in degrees):
   * 
   * <p>center = (80, 170),   size = (40, 60)   ->  lat = [60, 100],   lng = [140, -160]
   * <p>center = (10, 40),    size = (210, 400) ->  lat = [-90, 90],   lng = [-180, 180]
   * <p>center = (-90, 180),  size = (20, 50)   ->  lat = [-90, -80],  lng = [155, -155]
   */
  public static S2LatLngRect fromCenterSize(S2LatLng center, S2LatLng size) {
    return fromPoint(center).expanded(size.mul(0.5));
  }

  /** Convenience method to construct a rectangle containing a single point. */
  public static S2LatLngRect fromPoint(S2LatLng p) {
    // assert (p.isValid());
    return new S2LatLngRect(p, p);
  }

  /**
   * Convenience method to construct the minimal bounding rectangle containing the two given 
   * normalized points. This is equivalent to starting with an empty rectangle and calling
   * addPoint() twice. Note that it is different than the {@link #S2LatLngRect(S2LatLng, S2LatLng)}
   * constructor, where the first point is always used as the lower-left corner of the resulting
   * rectangle.
   */
  public static S2LatLngRect fromPointPair(S2LatLng p1, S2LatLng p2) {
    // assert (p1.isValid() && p2.isValid());
    return new S2LatLngRect(R1Interval.fromPointPair(p1.lat().radians(), p2.lat().radians()),
        S1Interval.fromPointPair(p1.lng().radians(), p2.lng().radians()));
  }

  /**
   * Returns a latitude-longitude rectangle that contains the edge from "a" to "b". Both points must
   * be unit-length. Note that the bounding rectangle of an edge can be larger than the bounding
   * rectangle of its endpoints.
   */
  public static S2LatLngRect fromEdge(S2Point a, S2Point b) {
    // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
    S2LatLngRect r = fromPointPair(new S2LatLng(a), new S2LatLng(b));

    // Check whether the min/max latitude occurs in the edge interior. We find the normal to the
    // plane containing AB, and then a vector "dir" in this plane that also passes through the
    // equator. We use RobustCrossProd to ensure that the edge normal is accurate even when the two
    // points are very close together.
    S2Point ab = S2.robustCrossProd(a, b);
    S2Point dir = S2Point.crossProd(ab, S2Point.Z_POS);
    double da = dir.dotProd(a);
    double db = dir.dotProd(b);
    if (da * db >= 0) {
      // Minimum and maximum latitude are attained at the vertices.
      return r;
    }
    // Minimum/maximum latitude occurs in the edge interior. This affects the latitude bounds but
    // not the longitude bounds.
    double absLat = Math.acos(Math.abs(ab.z / ab.norm()));
    if (da < 0) {
      return new S2LatLngRect(new R1Interval(r.lat().lo(), absLat), r.lng());
    } else {
      return new S2LatLngRect(new R1Interval(-absLat, r.lat().hi()), r.lng());
    }
  }

  /**
   * Constructs a rectangle from minimum and maximum latitudes and longitudes. If 
   * {@code lo.lng() > hi.lng()}, the rectangle spans the 180 degree longitude line. Both points
   * must be normalized, with {@code lo.lat() <= hi.lat()}.  The rectangle contains all the points p
   * such that {@code lo <= p && p <= hi}, where '<=' is defined in the obvious way.
   */
  public S2LatLngRect(final S2LatLng lo, final S2LatLng hi) {
    super(lo, hi);
    // assert (isValid());
  }

  /** Constructs a rectangle from latitude and longitude intervals. */
  public S2LatLngRect(R1Interval lat, S1Interval lng) {
    super(lat, lng);
    // assert (isValid());
  }

  @Override
  public final R1Interval lat() {
    // It is OK to return the instance field because S2LatLngRect won't mutate its 'lat' field.
    return lat;
  }

  @Override
  public final S1Interval lng() {
    // It is OK to return the instance field because S2LatLngRect won't mutate its 'lng' field.
    return lng;
  }
  
  /**
   * Returns a new rectangle that includes this rectangle and the given point, expanding this
   * rectangle to include the point by the minimum amount possible.
   */
  @CheckReturnValue
  public S2LatLngRect addPoint(S2Point p) {
    return addPoint(new S2LatLng(p));
  }

  /**
   * Returns a new rectangle that includes this rectangle and the given S2LatLng, expanding this
   * rectangle to include the point by the minimum amount possible. The S2LatLng argument must be
   * normalized.
   */
  @CheckReturnValue
  public S2LatLngRect addPoint(S2LatLng ll) {
    // assert (ll.isValid());
    R1Interval newLat = lat.addPoint(ll.lat().radians());
    S1Interval newLng = lng.addPoint(ll.lng().radians());
    return new S2LatLngRect(newLat, newLng);
  }

  /**
   * Returns a rectangle that contains all points whose latitude distance from this rectangle is at
   * most margin.lat(), and whose longitude distance from this rectangle is at most margin.lng(). In
   * particular, latitudes are clamped while longitudes are wrapped. Note that any expansion of an
   * empty interval remains empty, and both components of the given margin must be non-negative.
   * 
   * <p>Note that if an expanded rectangle contains a pole, it may not contain all possible lat/lng
   * representations of that pole.  Use polarClosure() if you do not want this behavior.
   *
   * <p>NOTE: If you are trying to grow a rectangle by a certain *distance* on the sphere
   * (e.g. 5km), use the convolveWithCap() method instead.
   */
  @CheckReturnValue
  public S2LatLngRect expanded(S2LatLng margin) {
    // assert (margin.lat().radians() >= 0 && margin.lng().radians() >= 0);
    return new S2LatLngRect(
        lat.expanded(margin.lat().radians()).intersection(fullLat()),
        lng.expanded(margin.lng().radians()));
  }
 
  /**
   * If the rectangle does not include either pole, return it unmodified. Otherwise expand the
   * longitude range to full() so that the rectangle contains all possible representations of the
   * contained pole(s).
   */
  @CheckReturnValue
  public S2LatLngRect polarClosure() {
    if (lat.lo() == -S2.M_PI_2 || lat.hi() == S2.M_PI_2) {
      return new S2LatLngRect(lat, S1Interval.full());
    } else {
      return this;
    }
  }

  /**
   * Returns the smallest rectangle containing the union of this rectangle and the given rectangle.
   */
  @CheckReturnValue
  public S2LatLngRect union(S2LatLngRectBase other) {
    return new S2LatLngRect(lat.union(other.lat), lng.union(other.lng));
  }

  /** 
   * Returns the smallest rectangle containing the intersection of this rectangle and the given
   * rectangle.  Note that the region of intersection may consist of two disjoint rectangles, in
   * which case a single rectangle spanning both of them is returned.
   */ 
  @CheckReturnValue
  public S2LatLngRect intersection(S2LatLngRectBase other) {
    R1Interval intersectLat = lat.intersection(other.lat);
    S1Interval intersectLng = lng.intersection(other.lng);
    if (intersectLat.isEmpty() || intersectLng.isEmpty()) {
      // The lat/lng ranges must either be both empty or both non-empty.
      return empty();
    }
    return new S2LatLngRect(intersectLat, intersectLng);
  }

  /**
   * Returns a rectangle that contains the convolution of this rectangle with a cap of the given
   * angle. This expands the rectangle by a fixed distance (as opposed to growing the rectangle in
   * latitude-longitude space). The returned rectangle includes all points whose minimum distance to
   * the original rectangle is at most the given angle.
   */
  @CheckReturnValue
  public S2LatLngRect convolveWithCap(S1Angle angle) {
    // The most straightforward approach is to build a cap centered on each vertex and take the
    // union of all the bounding rectangles (including the original rectangle; this is necessary for
    // very large rectangles).

    // Optimization: convert the angle to a height exactly once.
    S2Cap cap = S2Cap.fromAxisAngle(S2Point.X_POS, angle);

    S2LatLngRect r = this;
    for (int k = 0; k < 4; ++k) {
      S2Cap vertexCap = S2Cap.fromAxisHeight(getVertex(k).toPoint(), cap.height());
      r = r.union(vertexCap.getRectBound());
    }
    return r;
  }

  // NOTE: This should be marked as @Override, but clone() isn't present in GWT's version of Object,
  // so we can't mark it as such.
  public S2Region clone() {
    return new S2LatLngRect(this.lo(), this.hi());
  }

  @Override
  public S2LatLngRect getRectBound() {
    return this;
  }

  /**
   * This class is a builder for S2LatLngRect instances.  This is much more efficient when creating
   * the bounds from numerous operations, as it ensures that the S2LatLngRect is only created once.
   * 
   * <p>Example usage:
   * 
   * {@code
   *   S2LatLngRect union(List<S2LatLng> points) {
   *     S2LatLngRect.Builder builder = new S2LatLngRect.Builder();
   *     for (S2LatLng point : points) {
   *       builder.addPoint(point);
   *     }
   *     return builder.build();
   *   }
   * }
   */
  public static final strictfp class Builder extends S2LatLngRectBase {
    public Builder(final S2LatLng lo, final S2LatLng hi) {
      super(lo, hi);
    }

    public Builder(R1Interval lat, S1Interval lng) {
      super(lat, lng);
    }

    @Override
    public final R1Interval lat() {
      // 'lat' is copied here to avoid further changes in the builder being visible in the returned
      // object.
      return new R1Interval(lat);
    }
    
    @Override
    public final S1Interval lng() {
      // 'lng' is copied here to avoid further changes in the builder being visible in the returned
      // object.
      return new S1Interval(lng);
    }
    
    /** Returns a new immutable S2LatLngRect copied from the current state of this builder. */
    public S2LatLngRect build() {
      return new S2LatLngRect(new R1Interval(lat), new S1Interval(lng));
    }

    /** A builder initialized to be empty (such that it doesn't contain anything). */
    public static Builder empty() {
      return new Builder(R1Interval.empty(), S1Interval.empty());
    }
   
    /** Sets the rectangle to the full rectangle. */
    public void setFull() {
      lat.set(-S2.M_PI_2, S2.M_PI_2);
      lng.setFull();    
    }

    public void addPoint(S2Point p) {
      addPoint(new S2LatLng(p));
    }

    /** 
     * Increases the size of the bounding rectangle to include the given point. The rectangle is
     * expanded by the minimum amount possible.
     */
    public void addPoint(S2LatLng ll) {
      // assert (ll.isValid());
      lat.unionInternal(ll.lat().radians());
      lng.unionInternal(S1Interval.fromPoint(ll.lng().radians()));
    }

    /**
     * Mutates the rectangle to contain all points whose latitude distance from this rectangle is at
     * most margin.lat(), and whose longitude distance from this rectangle is at most margin.lng().
     * In particular, latitudes are clamped while longitudes are wrapped. Note that any expansion of
     * an empty interval remains empty, and both components of the given margin must be
     * non-negative.
     *
     * <p>NOTE: If you are trying to grow a rectangle by a certain *distance* on the sphere
     * (e.g. 5km), use the convolveWithCap() method instead.
     */
    public void expanded(S2LatLng margin) {
      // assert (margin.lat().radians() >= 0 && margin.lng().radians() >= 0);
      lat.expandedInternal(margin.lat().radians());
      lat.intersectionInternal(fullLat());
      lng.expandedInternal(margin.lng().radians());
    }

    /**
     * If the rectangle does not include either pole, leave it unmodified. Otherwise expand the
     * longitude range to full() so that the rectangle contains all possible representations of the
     * contained pole(s).
     */
    public void polarClosure() {
      if (lat.lo() == -S2.M_PI_2 || lat.hi() == S2.M_PI_2) {
        lng.setFull();
      }
    }

    /**
     * Mutates this rectangle to be the smallest rectangle containing the union of the current and
     * given rectangles.
     */
    public void union(S2LatLngRect other) {
      lat.unionInternal(other.lat);
      lng.unionInternal(other.lng);
    }

    /**
     * Mutates this rectangle to be the smallest rectangle containing the intersection of the
     * current and given rectangles. Note that the region of intersection may consist of two
     * disjoint rectangles, in which case we set the rectangle to be a single rectangle spanning
     * both of them.
     */
    public void intersection(S2LatLngRect other) {
      lat.intersectionInternal(other.lat);
      lng.intersectionInternal(other.lng);
      // The lat/lng ranges must either be both empty or both non-empty.
      if (lat.isEmpty() && !lng.isEmpty()) {
        lng.setEmpty();
      } else if (lng.isEmpty() && !lat.isEmpty()) {
        lat.setEmpty();
      }
    }

    /**
     * Mutates the current rectangle to contain the convolution of this rectangle with a cap of the
     * given angle. This expands the rectangle by a fixed distance (as opposed to growing the
     * rectangle in latitude-longitude space). The new rectangle includes all points whose minimum
     * distance to the original rectangle is at most the given angle.
     */
    public void convolveWithCap(S1Angle angle) {
      S2Cap cap = S2Cap.fromAxisAngle(S2Point.X_POS, angle);
      double latLo = lat.lo();
      double latHi = lat.hi();
      double lngLo = lng.lo();
      double lngHi = lng.hi();

      union(S2Cap.fromAxisHeight(
          S2LatLng.fromRadians(latLo, lngLo).toPoint(), cap.height()).getRectBound());
      union(S2Cap.fromAxisHeight(
          S2LatLng.fromRadians(latLo, lngHi).toPoint(), cap.height()).getRectBound());
      union(S2Cap.fromAxisHeight(
          S2LatLng.fromRadians(latHi, lngLo).toPoint(), cap.height()).getRectBound());
      union(S2Cap.fromAxisHeight(
          S2LatLng.fromRadians(latHi, lngHi).toPoint(), cap.height()).getRectBound());
    }

    // NOTE: This should be marked as @Override, but clone() isn't present in GWT's version of
    // Object, so we can't mark it as such.
    public S2Region clone() {
      return new S2LatLngRect(this.lo(), this.hi());
    }

    @Override
    public S2LatLngRect getRectBound() {
      return build();
    }
  }
}