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

import static org.junit.Assert.assertEquals;

import com.google.common.annotations.GwtIncompatible;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Testing S2Earth.
 * @author ericv@google.com (Eric Veach)
 * @author norris@google.com (Norris Boyd)
 */
@RunWith(JUnit4.class)
@GwtIncompatible
public class S2EarthTest {

  @Test
  public void testAngleConversion() {
    assertEquals(S2Earth.toMeters(S1Angle.degrees(180)), S2Earth.getRadiusMeters() * Math.PI, 0);
    assertEquals(S2Earth.toKm(S1Angle.radians(0.5)), 0.5 * S2Earth.getRadiusKm(), 0);
    assertEquals(1, S2Earth.kmToRadians(S2Earth.getRadiusMeters() / 1000), 0);
    assertEquals(S2Earth.radiansToKm(0.5), 0.5 * S2Earth.getRadiusKm(), 0);
    assertEquals(0.3, S2Earth.metersToRadians(S2Earth.radiansToKm(0.3) * 1000), 0);
    assertEquals(2500, S2Earth.radiansToMeters(S2Earth.kmToRadians(2.5)), 0);
  }

  @Test
  public void testSolidAngleConversion() {
    assertEquals(1, S2Earth.squareKmToSteradians(Math.pow(S2Earth.getRadiusMeters() / 1000, 2)), 0);
    assertEquals(
        Math.pow(0.5 * S2Earth.getRadiusKm(), 2),
        S2Earth.steradiansToSquareKm(Math.pow(0.5, 2)),
        0);
    assertEquals(
        Math.pow(0.3, 2),
        S2Earth.squareMetersToSteradians(Math.pow(S2Earth.radiansToKm(0.3) * 1000, 2)),
        0);

    // This one test doesn't equal exactly. C++ tests allow an epsilon of 4 ULPs, so this is not an
    // issue of correctness; use 1 ULP here.
    final double expected = Math.pow(2500, 2);
    assertEquals(
        expected,
        S2Earth.steradiansToSquareMeters(Math.pow(S2Earth.kmToRadians(2.5), 2)),
        Math.ulp(expected));
  }

  @Test
  public void testGetDistance() {
    S2Point north = new S2Point(0, 0, 1);
    S2Point south = new S2Point(0, 0, -1);
    S2Point west = new S2Point(0, -1, 0);

    assertEquals(0, S2Earth.getDistanceKm(west, west), 0);
    assertEquals(
        S2Earth.getDistanceMeters(north, west), (Math.PI / 2) * S2Earth.getRadiusMeters(), 0);
    assertEquals(
        S2Earth.getDistanceKm(S2LatLng.fromRadians(0, 0.6), S2LatLng.fromRadians(0, -0.4)),
        S2Earth.getRadiusKm(),
        0);

    // This one test doesn't equal exactly. C++ tests allow an epsilon of 4 ULPs, so this is not an
    // issue of correctness; use 1 ULP here.
    final double expected = 1000 * S2Earth.getRadiusKm() * Math.PI / 4;
    assertEquals(
        S2Earth.getDistanceMeters(S2LatLng.fromDegrees(80, 27), S2LatLng.fromDegrees(55, -153)),
        expected,
        Math.ulp(expected));
  }

  @Test
  public void testGetInitialBearings() {
    S2LatLng equator0 = S2LatLng.fromDegrees(0, 0);
    S2LatLng equator50 = S2LatLng.fromDegrees(0, 50);
    S2LatLng equator100 = S2LatLng.fromDegrees(0, 100);

    // Westward on Equator.
    assertEquals(S2Earth.getInitialBearing(equator50, equator100), S1Angle.degrees(90));

    // Eastwards on Equator.
    assertEquals(S2Earth.getInitialBearing(equator100, equator0), S1Angle.degrees(-90));

    // Northward from Meridian.
    assertEquals(
        S2Earth.getInitialBearing(S2LatLng.fromDegrees(16, 28), S2LatLng.fromDegrees(81, 28)),
        S1Angle.degrees(0));

    // Southward from Meridian.
    assertEquals(
        S2Earth.getInitialBearing(S2LatLng.fromDegrees(24, 64), S2LatLng.fromDegrees(-27, 64)),
        S1Angle.degrees(180));

    // Towards the north pole.
    assertEquals(
        S2Earth.getInitialBearing(S2LatLng.fromDegrees(12, 75), S2LatLng.fromDegrees(90, 50))
            .degrees(),
        0.0,
        1e-7);

    // Towards the south pole.
    assertEquals(
        S2Earth.getInitialBearing(S2LatLng.fromDegrees(-35, 105), S2LatLng.fromDegrees(-90, -120))
            .degrees(),
        180.0,
        1e-7);

    S2LatLng spain = S2LatLng.fromDegrees(40.4379332, -3.749576);
    S2LatLng japan = S2LatLng.fromDegrees(35.6733227, 139.6403486);
    assertEquals(S2Earth.getInitialBearing(spain, japan).degrees(), 29.2, 1e-2);
    assertEquals(S2Earth.getInitialBearing(japan, spain).degrees(), -27.2, 1e-2);
  }
}
