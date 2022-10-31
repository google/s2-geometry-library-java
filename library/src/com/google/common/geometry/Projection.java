/*
 * Copyright 2022 Google Inc.
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

import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * For the purposes of the S2 library, a projection is a function that maps between S2Points and
 * R2Points. It can also define the coordinate wrapping behavior along each axis.
 */
@JsType
public interface Projection {
  /* Converts a point on the sphere to a projected 2D point. */
  default R2Vector project(S2Point p) {
    return fromLatLng(new S2LatLng(p));
  }

  /**
   * Converts a projected 2D point to a point on the sphere.
   *
   * <p>If wrapping is defined for a given axis (see below), then this method should accept any real
   * number for the corresponding coordinate.
   */
  default S2Point unproject(R2Vector p) {
    return toLatLng(p).toPoint();
  }

  /**
   * Convenience function equivalent to Project(ll.ToPoint()), but the implementation may be more
   * efficient.
   */
  R2Vector fromLatLng(S2LatLng ll);

  /**
   * Convenience function equivalent to S2LatLng(Unproject(p)), but the implementation may be more
   * efficient.
   */
  S2LatLng toLatLng(R2Vector p);

  /**
   * Returns the point obtained by interpolating the given fraction of the distance along the line
   * from A to B. Almost all projections should use the default implementation of this method, which
   * simply interpolates linearly in R2 space. Fractions < 0 or > 1 result in extrapolation instead.
   *
   * <p>The only reason to override this method is if you want edges to be defined as something
   * other than straight lines in the 2D projected coordinate system. For example, using a
   * third-party library such as GeographicLib you could define edges as geodesics over an ellipsoid
   * model of the Earth. (Note that very few data sets define edges this way.)
   *
   * <p>Also note that there is no reason to define a projection where edges are geodesics over the
   * sphere, because this is the native S2 interpretation.
   */
  static R2Vector interpolate(double f, R2Vector a, R2Vector b) {
    return a.mul(1 - f).add(b.mul(f));
  }

  /**
   * Defines the coordinate wrapping distance along each axis. If this value is non-zero for a given
   * axis, the coordinates are assumed to "wrap" with the given period. For example, if
   * wrap_distance.y() == 360 then (x, y) and (x, y + 360) should map to the same S2Point.
   *
   * <p>This information is used to ensure that edges takes the shortest path between two given
   * points. For example, if coordinates represent (latitude, longitude) pairs in degrees and
   * wrap_distance().y() == 360, then the edge (5:179, 5:-179) would be interpreted as spanning 2
   * degrees of longitude rather than 358 degrees.
   *
   * <p>If a given axis does not wrap, its wrap_distance should be set to zero.
   */
  R2Vector wrapDistance();

  /**
   * Helper function that wraps the coordinates of B if necessary in order to obtain the shortest
   * edge AB. For example, suppose that A = [170, 20], B = [-170, 20], and the projection wraps so
   * that [x, y] == [x + 360, y]. Then this function would return [190, 20] for point B (reducing
   * the edge length in the "x" direction from 340 to 20).
   */
  default R2Vector wrapDestination(R2Vector a, R2Vector b) {
    R2Vector wrap = wrapDistance();
    double x = b.x();
    double y = b.y();
    // The code below ensures that "b" is unmodified unless wrapping is required.
    if (wrap.x() > 0 && Math.abs(x - a.x()) > 0.5 * wrap.x()) {
      x = a.x() + Platform.IEEEremainder(x - a.x(), wrap.x());
    }
    if (wrap.y() > 0 && Math.abs(y - a.y()) > 0.5 * wrap.y()) {
      y = a.y() + Platform.IEEEremainder(y - a.y(), wrap.y());
    }
    return new R2Vector(x, y);
  }

  /**
   * MercatorProjection defines the spherical Mercator projection. Google Maps uses this projection
   * together with WGS84 coordinates, in which case it is known as the "Web Mercator" projection
   * (see Wikipedia). This class makes no assumptions regarding the coordinate system of its input
   * points, but simply applies the spherical Mercator projection to them.
   *
   * <p>The Mercator projection is finite in width (x) but infinite in height (y). "x" corresponds
   * to longitude, and spans a finite range such as [-180, 180] (with coordinate wrapping), while
   * "y" is a function of latitude and spans an infinite range. (As "y" coordinates get larger,
   * points get closer to the north pole but never quite reach it.) The north and south poles have
   * infinite "y" values. (Note that this will cause problems if you tessellate a Mercator edge
   * where one endpoint is a pole. If you need to do this, clip the edge first so that the "y"
   * coordinate is no more than about 5 * max_x.)
   */
  @JsType
  public final class MercatorProjection implements Projection {

    private final double xWrap;
    private final double toRadians;
    private final double fromRadians;

    /**
     * Default constructor with the projected 'x' coordinate in [-pi, pi]. Note 'y' coordinate is
     * unbounded from -infinity to infinity.
     */
    @JsIgnore
    public static MercatorProjection createInRadians() {
      return new MercatorProjection(Math.PI);
    }

    /**
     * Constructs a Mercator projection where "x" corresponds to longitude in the range [-maxX,
     * maxX] and "y" corresponds to latitude and can be any real number. The horizontal and vertical
     * scales are equal locally.
     *
     * @param maxX defines the maximal values for x in the positive and negative direction. The
     *     projected 'x' coordinate will be in [-maxX, maxX]. Note that the y coordinate for
     *     projections is unbound.
     */
    public MercatorProjection(double maxX) {
      xWrap = 2 * maxX;
      toRadians = Math.PI / maxX;
      fromRadians = maxX / Math.PI;
    }

    @Override
    public R2Vector fromLatLng(S2LatLng ll) {
      double sinPhi = Math.sin(ll.latRadians());
      double y = 0.5 * Math.log((1 + sinPhi) / (1 - sinPhi));
      return new R2Vector(fromRadians * ll.lngRadians(), fromRadians * y);
    }

    @Override
    public S2LatLng toLatLng(R2Vector p) {
      double x = toRadians * Platform.IEEEremainder(p.x(), xWrap);
      double k = Math.exp(2 * toRadians * p.y());
      double y = Double.isInfinite(k) ? Math.PI / 2 : Math.asin((k - 1) / (k + 1));
      return S2LatLng.fromRadians(y, x);
    }

    @Override
    public R2Vector wrapDistance() {
      return new R2Vector(xWrap, 0);
    }
  }

  /**
   * PlateCarreeProjection defines the "plate carree" (square plate) projection, which converts
   * points on the sphere to (longitude, latitude) pairs. Coordinates can be scaled so that they
   * represent radians, degrees, etc, but the projection is always centered around (latitude=0,
   * longitude=0).
   *
   * <p>Note that (x, y) coordinates are backwards compared to the usual (latitude, longitude)
   * ordering, in order to match the usual convention for graphs in which "x" is horizontal and "y"
   * is vertical.
   */
  @JsType
  public final class PlateCarreeProjection implements Projection {

    private final double xWrap;
    private final double toRadians;
    private final double fromRadians;

    /**
     * Constructor which by default sets the scale to PI. This will cause points in the projection
     * to equal the latitudes and longitudes and vice versa with some accommodation for wrapping the
     * x coordinate in Unproject.
     */
    @JsIgnore
    public static PlateCarreeProjection createInRadians() {
      return new PlateCarreeProjection(Math.PI);
    }

    /**
     * Constructs the plate carree projection where the x coordinates (longitude) span [-x_scale,
     * x_scale] and the y coordinates (latitude) span [-x_scale/2, x_scale/2]. For example if
     * x_scale==180 then the x range is [-180, 180] and the y range is [-90, 90].
     *
     * @param scale defines the ranges of the output from the projection.
     */
    public PlateCarreeProjection(double scale) {
      xWrap = 2 * scale;
      toRadians = Math.PI / scale;
      fromRadians = scale / Math.PI;
    }

    @Override
    public R2Vector fromLatLng(S2LatLng ll) {
      return new R2Vector(fromRadians * ll.lngRadians(), fromRadians * ll.latRadians());
    }

    @Override
    public S2LatLng toLatLng(R2Vector p) {
      return S2LatLng.fromRadians(
          toRadians * p.y(), toRadians * Platform.IEEEremainder(p.x(), xWrap));
    }

    @Override
    public R2Vector wrapDistance() {
      return new R2Vector(xWrap, 0);
    }
  }
}
