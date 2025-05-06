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

import static com.google.common.geometry.S2.M_PI_2;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.primitives.Doubles;
import com.google.errorprone.annotations.CheckReturnValue;
import com.google.errorprone.annotations.Immutable;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * An immutable representation of a point on the unit sphere, as a pair of latitude-longitude
 * coordinates.
 *
 * <p>Like the rest of the "geometry" package, the intent is to represent spherical geometry as a
 * mathematical abstraction, so functions that are specifically related to the Earth's geometry
 * (e.g. easting/northing conversions) are elsewhere. Note that the serialized form of this class is
 * not stable and should not be relied upon for long-term persistence. Avoid using Java
 * serialization, use the provided S2Coder instead.
 *
 * @author danieldanciu@google.com (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@Immutable
@JsType
@SuppressWarnings("Assertion")
public final class S2LatLng implements Serializable {
  /** The center point of the lat/lng coordinate system. */
  public static final S2LatLng CENTER = new S2LatLng(0.0, 0.0);

  /** A {@link S2Coder} of latlngs. */
  public static final S2Coder<S2LatLng> CODER =
      new S2Coder<S2LatLng>() {
        @Override
        public void encode(S2LatLng value, OutputStream output) throws IOException {
          LittleEndianOutput.writeDouble(output, value.latRadians);
          LittleEndianOutput.writeDouble(output, value.lngRadians);
        }

        @Override
        public S2LatLng decode(Bytes data, Cursor cursor) throws IOException {
          try {
            double lat = data.readLittleEndianDouble(cursor.position);
            double lng = data.readLittleEndianDouble(cursor.position + Doubles.BYTES);
            cursor.position += 2 * Doubles.BYTES;
            return S2LatLng.fromRadians(lat, lng);
          } catch (ArrayIndexOutOfBoundsException e) {
            throw new IOException("Insufficient or invalid input bytes: ", e);
          }
        }

        @Override
        public boolean isLazy() {
          return false;
        }
      };

  private final double latRadians;
  private final double lngRadians;

  /**
   * Returns a new S2LatLng specified in radians. Some methods on S2LatLng require that the latitude
   * and longitude are valid. If you are unsure, use {@link #isValid(double, double)} to check that
   * the lat and lng are valid before attempting to construct the S2LatLng, or use {@link
   * #isValid()} after construction.
   */
  public static S2LatLng fromRadians(double latRadians, double lngRadians) {
    return new S2LatLng(latRadians, lngRadians);
  }

  /** Returns a new S2LatLng converted from degrees. */
  public static S2LatLng fromDegrees(double latDegrees, double lngDegrees) {
    return new S2LatLng(S1Angle.degrees(latDegrees), S1Angle.degrees(lngDegrees));
  }

  /** Returns a new S2LatLng converted from tens of microdegrees. */
  public static S2LatLng fromE5(int latE5, int lngE5) {
    return new S2LatLng(S1Angle.e5(latE5), S1Angle.e5(lngE5));
  }

  /** Returns a new S2LatLng converted from microdegrees. */
  public static S2LatLng fromE6(int latE6, int lngE6) {
    return new S2LatLng(S1Angle.e6(latE6), S1Angle.e6(lngE6));
  }

  /** Returns a new S2LatLng converted from tenths of a microdegree. */
  public static S2LatLng fromE7(int latE7, int lngE7) {
    return new S2LatLng(S1Angle.e7(latE7), S1Angle.e7(lngE7));
  }

  /** Returns a new S2LatLng converted from an S2Point (not necessarily normalized). */
  public static S2LatLng fromPoint(S2Point p) {
    return new S2LatLng(p);
  }

  /**
   * Returns the latitude in radians of the given point, which is not required to be unit length.
   */
  public static S1Angle latitude(S2Point p) {
    // We use atan2 rather than asin because the input vector is not necessarily unit length, and
    // atan2 is much more accurate than asin near the poles. The "+ 0.0" is to ensure that points
    // with coordinates of -0.0 and +0.0 (which compare equal) are converted to identical S2LatLng
    // values, since even though -0.0 == +0.0 they can be formatted differently.
    return S1Angle.radians(atan2(p.z + 0.0, sqrt(p.x * p.x + p.y * p.y)));
  }

  /**
   * Returns the longitude in radians of the given point, which is not required to be unit length.
   */
  public static S1Angle longitude(S2Point p) {
    // The "+ 0.0" is to ensure that points with coordinates of -0.0 and +0.0 (which compare equal)
    // are converted to identical S2LatLng values, since even though -0.0 == +0.0 and -180 == 180
    // degrees, they can be formatted differently. Also note that atan2(0, 0) is defined to be zero.
    return S1Angle.radians(atan2(p.y + 0.0, p.x + 0.0));
  }

  /** This is internal to avoid ambiguity about which units are expected. */
  @JsConstructor
  private S2LatLng(double latRadians, double lngRadians) {
    this.latRadians = latRadians;
    this.lngRadians = lngRadians;
  }

  /**
   * Basic constructor. The latitude and longitude must be within the ranges allowed by isValid()
   * below.
   */
  @JsIgnore
  public S2LatLng(S1Angle lat, S1Angle lng) {
    this(lat.radians(), lng.radians());
  }

  /** Default constructor for convenience when declaring arrays, etc. */
  @JsIgnore
  public S2LatLng() {
    this(0, 0);
  }

  /** Convert a point (not necessarily normalized) to an S2LatLng. */
  @JsIgnore
  public S2LatLng(S2Point p) {
    // The "+ 0.0" is to ensure that points with coordinates of -0.0 and +0.0 (which compare equal)
    // are converted to identical S2LatLng values, since even though -0.0 == +0.0 they can be
    // formatted differently.
    this(atan2(p.z + 0.0, sqrt(p.x * p.x + p.y * p.y)), atan2(p.y + 0.0, p.x + 0.0));
    // The latitude and longitude are already normalized. We use atan2 to compute the latitude
    // because the input vector is not necessarily unit length, and atan2 is much more accurate than
    // asin near the poles. Note that atan2(0, 0) is defined to be zero.
  }

  /** Returns the latitude of this point as a new S1Angle. */
  public S1Angle lat() {
    return S1Angle.radians(latRadians);
  }

  /** Returns the latitude of this point as radians. */
  public double latRadians() {
    return latRadians;
  }

  /** Returns the latitude of this point as degrees. */
  public double latDegrees() {
    return 180.0 / PI * latRadians;
  }

  /** Returns the longitude of this point as a new S1Angle. */
  public S1Angle lng() {
    return S1Angle.radians(lngRadians);
  }

  /** Returns the longitude of this point as radians. */
  public double lngRadians() {
    return lngRadians;
  }

  /** Returns the longitude of this point as degrees. */
  public double lngDegrees() {
    return 180.0 / PI * lngRadians;
  }

  /**
   * Return true if this LatLng is {@link normalized()}, i.e. the latitude is between -90 and 90
   * degrees inclusive and the longitude is between -180 and 180 degrees inclusive.
   */
  public boolean isValid() {
    return isValid(latRadians, lngRadians);
  }

  /**
   * Return true if the given lat and lng would be a valid LatLng, i.e. the latitude is between -90
   * and 90 degrees inclusive and the longitude is between -180 and 180 degrees inclusive.
   */
  public static boolean isValid(double latRadians, double lngRadians) {
    return abs(latRadians) <= M_PI_2 && abs(lngRadians) <= PI;
  }

  /**
   * Returns a new S2LatLng based on this instance for which {@link #isValid()} will be {@code
   * true}.
   *
   * <ul>
   *   <li>Latitude is clipped to the range {@code [-90, 90]}
   *   <li>Longitude is normalized to be in the range {@code [-180, 180]}
   * </ul>
   *
   * <p>If the current point is valid then the returned point will have the same coordinates.
   */
  @CheckReturnValue
  public S2LatLng normalized() {
    // IEEEremainder(x, 2 * PI) reduces its argument to the range [-PI, PI] inclusive, which is what
    // we want here.
    return new S2LatLng(
        max(-M_PI_2, min(M_PI_2, latRadians)), Platform.IEEEremainder(lngRadians, 2 * PI));
  }

  /**
   * Convert an S2LatLng to the equivalent unit-length vector (S2Point). Unnormalized values (see
   * {@link normalized()}) are wrapped around the sphere as would be expected based on their
   * definition as spherical angles. So for example the following pairs yield equivalent points
   * (modulo numerical error):
   *
   * <pre>
   *   (90.5, 10) =~ (89.5, -170)
   *   (a, b) =~ (a + 360 * n, b)
   * </pre>
   *
   * The maximum error in the result is 1.5 * DBL_EPSILON. (This does not include the error of
   * converting degrees, E5, E6, or E7 to radians.)
   *
   * <p>The latitude and longitude must be finite.
   */
  public S2Point toPoint() {
    assert Double.isFinite(latRadians) : latRadians;
    assert Double.isFinite(lngRadians) : lngRadians;
    // TODO(torrey): Consider enabling this assertion, even though it isn't actually required by
    // the toPoint() code. Clients calling toPoint() with unnormalized LatLngs may not actually be
    // consciously relying on this code to correctly wrap latitudes > 90 around the sphere, but are
    // instead using invalid data without knowing it.
    // assert isValid();
    double phi = latRadians;
    double theta = lngRadians;
    double cosphi = cos(phi);
    return new S2Point(cos(theta) * cosphi, sin(theta) * cosphi, sin(phi));
  }

  /**
   * Return the distance (measured along the surface of the sphere) to the given point.
   *
   * <p>This implements the Haversine formula, which is numerically stable for small distances but
   * only gets about 8 digits of precision for very large distances (e.g. antipodal points). Note
   * that 8 digits is still accurate to within about 10cm for a sphere the size of the Earth.
   *
   * <p>This could be fixed with another sin() and cos() below, but at that point you might as well
   * just convert both arguments to S2Points and compute the distance that way (which gives about 15
   * digits of accuracy for all distances).
   */
  public S1Angle getDistance(final S2LatLng o) {
    // TODO(torrey): As above, consider enabling these assertions even though they are not required
    // by the getDistance() code, because most violations are likely garbage data.
    // assert isValid();
    // assert o.isValid();
    double lat1 = latRadians;
    double lat2 = o.latRadians;
    double lng1 = lngRadians;
    double lng2 = o.lngRadians;
    double dlat = sin(0.5 * (lat2 - lat1));
    double dlng = sin(0.5 * (lng2 - lng1));
    double x = dlat * dlat + dlng * dlng * cos(lat1) * cos(lat2);
    return S1Angle.radians(2 * asin(sqrt(min(1.0, x))));
  }

  /** Returns the surface distance to the given point assuming a constant radius. */
  @JsMethod(name = "getDistanceWithRadius")
  public double getDistance(final S2LatLng o, double radius) {
    return getDistance(o).distance(radius);
  }

  /**
   * Adds the given point to this point. Note that there is no guarantee that the new point will be
   * <em>valid</em>.
   */
  @CheckReturnValue
  public S2LatLng add(final S2LatLng o) {
    return new S2LatLng(latRadians + o.latRadians, lngRadians + o.lngRadians);
  }

  /**
   * Subtracts the given point from this point. Note that there is no guarantee that the new point
   * will be <em>valid</em>.
   */
  @CheckReturnValue
  public S2LatLng sub(final S2LatLng o) {
    return new S2LatLng(latRadians - o.latRadians, lngRadians - o.lngRadians);
  }

  /**
   * Scales this point by the given scaling factor. Note that there is no guarantee that the new
   * point will be <em>valid</em>.
   */
  @CheckReturnValue
  public S2LatLng mul(final double m) {
    return new S2LatLng(latRadians * m, lngRadians * m);
  }

  @Override
  public boolean equals(Object that) {
    if (that instanceof S2LatLng) {
      S2LatLng o = (S2LatLng) that;
      return (latRadians == o.latRadians) && (lngRadians == o.lngRadians);
    }
    return false;
  }

  @Override
  public int hashCode() {
    long value = 17;
    value += 37 * value + Double.doubleToLongBits(latRadians);
    value += 37 * value + Double.doubleToLongBits(lngRadians);
    return (int) (value ^ (value >>> 32));
  }

  /**
   * Returns true if both the latitude and longitude of the given point are within {@code maxError}
   * radians of this point.
   */
  @JsMethod(name = "approxEqualsWithMaxError")
  public boolean approxEquals(S2LatLng o, double maxError) {
    return (abs(latRadians - o.latRadians) < maxError)
        && (abs(lngRadians - o.lngRadians) < maxError);
  }

  /**
   * Returns true if the given point is within {@code 1e-9} radians of this point. This corresponds
   * to a distance of less than {@code 1cm} at the surface of the Earth.
   */
  @CheckReturnValue
  public boolean approxEquals(S2LatLng o) {
    return approxEquals(o, 1e-9);
  }

  /** Returns a String with the lat and lng in radians. See also {@link #toStringDegrees()}. */
  @Override
  public String toString() {
    return "(" + latRadians + ", " + lngRadians + ")";
  }

  /** Returns a String with the lat and lng in degrees. */
  public String toStringDegrees() {
    return "(" + latDegrees() + ", " + lngDegrees() + ")";
  }
}
