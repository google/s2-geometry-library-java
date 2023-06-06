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
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * The Earth modeled as a sphere.
 *
 * <p>Provides many convenience functions so that it doesn't take 2 lines of code just to do a
 * single conversion. Note that the conversions between angles and distances on the Earth's surface
 * provided here rely on modeling Earth as a sphere; otherwise a given angle would correspond to a
 * range of distances depending on where the corresponding line segment was located.
 *
 * <p>More sophisticated Earth models (such as WSG84) should be used if required for accuracy or
 * interoperability.
 *
 * @author ericv@google.com (Eric Veach)
 * @author norris@google.com (Norris Boyd)
 */
@JsType
@GwtCompatible
public class S2Earth {
  private S2Earth() {}

  /**
   * Returns the Earth's mean radius, which is the radius of the equivalent sphere with the same
   * surface area. According to NASA, this value is 6371.01 +/- 0.02 km. The equatorial radius is
   * 6378.136 km, and the polar radius is 6356.752 km. They differ by one part in 298.257.
   *
   * <p>Reference: http://ssd.jpl.nasa.gov/phys_props_earth.html, which quotes Yoder, C.F. 1995,
   * "Astrometric and Geodetic Properties of Earth and the Solar System" in Global Earth Physics, A
   * Handbook of Physical Constants, AGU Reference Shelf 1, American Geophysical Union, Table 2.
   */
  public static double getRadiusMeters() {
    return 6371010.0;
  }

  /** Returns the Earth's mean radius as above, but in kilometers. */
  public static double getRadiusKm() {
    return 0.001 * getRadiusMeters();
  }

  /**
   * Returns the altitude of the lowest known point on Earth.
   *
   * <p>The lowest known point on Earth is the Challenger Deep with an altitude of -10898 meters
   * relative to the mean radius of the Earth.
   *
   * @return The altitude of the lowest known point on Earth.
   */
  @JsIgnore
  public static double getLowestAltitudeMeters() {
    return -10898;
  }

  /**
   * Returns the altitude of the highest known point on Earth.
   *
   * <p>The highest known point on Earth is Mount Everest with an altitude of 8846 meters above the
   * mean radius of the Earth.
   *
   * @return The altitude of the highest known point on Earth.
   */
  @JsIgnore
  public static double getHighestAltitudeMeters() {
    return 8846;
  }

  /** Converts the given distance in meters to an angle. */
  @JsIgnore
  public static S1Angle metersToAngle(double distanceMeters) {
    return S1Angle.radians(metersToRadians(distanceMeters));
  }

  /** Converts the given distance in meters to an S1ChordAngle. */
  @JsIgnore
  public static S1ChordAngle metersToChordAngle(double distanceMeters) {
    return S1ChordAngle.fromRadians(metersToRadians(distanceMeters));
  }

  /** Converts the given angle to meters. */
  public static double toMeters(S1Angle angle) {
    return angle.radians() * getRadiusMeters();
  }

  /** Converts the given S1ChordAngle to meters. */
  @JsIgnore
  public static double toMeters(S1ChordAngle chordAngle) {
    return chordAngle.radians() * getRadiusMeters();
  }

  /** Converts the given angle to kilometers. */
  public static double toKm(S1Angle angle) {
    return angle.radians() * getRadiusKm();
  }

  /** Converts the given kilometers to radians. */
  public static double kmToRadians(double km) {
    return km / getRadiusKm();
  }

  /** Converts the given radians to kilometers. */
  public static double radiansToKm(double radians) {
    return radians * getRadiusKm();
  }

  /** Converts the given meters to radians. */
  public static double metersToRadians(double meters) {
    return meters / getRadiusMeters();
  }

  /** Converts the given radians to meters. */
  public static double radiansToMeters(double radians) {
    return radians * getRadiusMeters();
  }

  // Functions that convert between areas on the unit sphere (as returned by the S2 library) and
  // areas on the Earth's surface. Note that the area of a region on the unit sphere is equal to the
  // solid angle it subtends from the sphere's center (measured in steradians).

  /** Converts the given square kilometers to steradians. */
  public static double squareKmToSteradians(double km2) {
    return km2 / (getRadiusKm() * getRadiusKm());
  }

  /** Converts the given square meters to steradians. */
  public static double squareMetersToSteradians(double m2) {
    return m2 / (getRadiusMeters() * getRadiusMeters());
  }

  /** Converts the given steradians to square kilometers. */
  public static double steradiansToSquareKm(double steradians) {
    return steradians * getRadiusKm() * getRadiusKm();
  }

  /** Converts the given steradians to square kilometers. */
  public static double steradiansToSquareMeters(double steradians) {
    return steradians * getRadiusMeters() * getRadiusMeters();
  }

  /** Returns the distance between two {@link S2Point}s on the globe in kilometers. */
  @JsMethod(name = "getDistanceBetweenPointsKm")
  public static double getDistanceKm(S2Point a, S2Point b) {
    return radiansToKm(a.angle(b));
  }

  /** Returns the distance between two {@link S2LatLng}s on the globe in kilometers. */
  @JsMethod(name = "getDistanceBetweenLatLngsKm")
  public static double getDistanceKm(S2LatLng a, S2LatLng b) {
    return toKm(a.getDistance(b));
  }

  /** Returns the distance between two {@link S2Point}s on the globe in meters. */
  @JsMethod(name = "getDistanceBetweenPointsMeters")
  public static double getDistanceMeters(S2Point a, S2Point b) {
    return radiansToMeters(a.angle(b));
  }

  /** Returns the distance between two {@link S2LatLng}s on the globe in meters. */
  @JsMethod(name = "getDistanceBetweenLatLngsMeters")
  public static double getDistanceMeters(S2LatLng a, S2LatLng b) {
    return toMeters(a.getDistance(b));
  }

  /**
   * Returns the Haversine of the angle. The versine is 1-cos(theta) and the haversine is "half" of
   * the versine or (1-cos(theta))/2.
   *
   * <p>http://en.wikipedia.org/wiki/Haversine_formula Haversine(x) has very good numerical
   * stability around zero. Haversine(x) == (1-cos(x))/2 == sin(x/2)^2; must be implemented with the
   * second form to reap the numerical benefits.
   *
   * @return The haversine of the angle.
   */
  public static double haversine(double angle) {
    double halfSin = Math.sin(angle / 2);
    return halfSin * halfSin;
  }

  /**
   * Calculates the bearing angle between two S2LatLngs.
   *
   * <p>Sourced from http://www.movable-type.co.uk/scripts/latlong.html.
   *
   * @return The bearing angle between two S2LatLngs.
   */
  public static S1Angle getInitialBearing(S2LatLng a, S2LatLng b) {
    double lat1 = a.latRadians();
    double cosLat2 = Math.cos(b.latRadians());
    double latDiff = b.latRadians() - a.latRadians();
    double lngDiff = b.lngRadians() - a.lngRadians();

    double x = Math.sin(latDiff) + Math.sin(lat1) * cosLat2 * 2 * haversine(lngDiff);
    double y = Math.sin(lngDiff) * cosLat2;
    return S1Angle.radians(Math.atan2(y, x));
  }
}
