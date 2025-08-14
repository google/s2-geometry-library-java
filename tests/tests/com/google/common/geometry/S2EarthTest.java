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

import static com.google.common.truth.Truth.assertThat;
import static org.junit.Assert.assertEquals;

import com.google.common.annotations.GwtIncompatible;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Testing S2Earth.
 *
 * @author ericv@google.com (Eric Veach)
 * @author norris@google.com (Norris Boyd)
 */
@RunWith(JUnit4.class)
@GwtIncompatible
public class S2EarthTest {

  @Test
  public void testAngleConversion() {
    assertThat(S2Earth.getRadiusMeters() * Math.PI)
        .isWithin(0)
        .of(S2Earth.toMeters(S1Angle.degrees(180)));
    assertThat(0.5 * S2Earth.getRadiusKm()).isWithin(0).of(S2Earth.toKm(S1Angle.radians(0.5)));
    assertThat(S2Earth.kmToRadians(S2Earth.getRadiusMeters() / 1000)).isWithin(0).of(1);
    assertThat(0.5 * S2Earth.getRadiusKm()).isWithin(0).of(S2Earth.radiansToKm(0.5));
    assertThat(S2Earth.metersToRadians(S2Earth.radiansToKm(0.3) * 1000)).isWithin(0).of(0.3);
    assertThat(S2Earth.radiansToMeters(S2Earth.kmToRadians(2.5))).isWithin(0).of(2500);
  }

  @Test
  public void testChordAngleConversion() {
    double chordAngleEpsilon = 1e-15;
    double quarterCircumferenceMeters = S2Earth.getRadiusMeters() * Math.PI / 2.0;
    assertThat(S2Earth.metersToChordAngle(quarterCircumferenceMeters).getLength2())
        .isWithin(chordAngleEpsilon)
        .of(2);
    assertThat(2 * S2Earth.getRadiusMeters())
        .isWithin(0)
        .of(S2Earth.toMeters(S1ChordAngle.fromRadians(2)));
    assertThat(S2Earth.getRadiusMeters() * Math.PI)
        .isWithin(chordAngleEpsilon)
        .of(S2Earth.toMeters(S1ChordAngle.fromDegrees(180)));
  }

  @Test
  public void testSolidAngleConversion() {
    assertThat(S2Earth.squareKmToSteradians(Math.pow(S2Earth.getRadiusMeters() / 1000, 2)))
        .isWithin(0)
        .of(1);
    assertThat(S2Earth.steradiansToSquareKm(Math.pow(0.5, 2)))
        .isWithin(0)
        .of(Math.pow(0.5 * S2Earth.getRadiusKm(), 2));
    assertThat(S2Earth.squareMetersToSteradians(Math.pow(S2Earth.radiansToKm(0.3) * 1000, 2)))
        .isWithin(0)
        .of(Math.pow(0.3, 2));

    // This one test doesn't equal exactly. C++ tests allow an epsilon of 4 ULPs, so this is not an
    // issue of correctness; use 1 ULP here.
    final double expected = Math.pow(2500, 2);
    assertThat(S2Earth.steradiansToSquareMeters(Math.pow(S2Earth.kmToRadians(2.5), 2)))
        .isWithin(Math.ulp(expected))
        .of(expected);
  }

  @Test
  public void testGetDistance() {
    S2Point north = new S2Point(0, 0, 1);
    S2Point south = new S2Point(0, 0, -1);
    S2Point west = new S2Point(0, -1, 0);

    assertThat(S2Earth.getDistanceKm(west, west)).isWithin(0).of(0);
    assertThat((Math.PI / 2) * S2Earth.getRadiusMeters())
        .isWithin(0)
        .of(S2Earth.getDistanceMeters(north, west));
    assertThat(S2Earth.getRadiusKm())
        .isWithin(0)
        .of(S2Earth.getDistanceKm(S2LatLng.fromRadians(0, 0.6), S2LatLng.fromRadians(0, -0.4)));

    // This one test doesn't equal exactly. C++ tests allow an epsilon of 4 ULPs, so this is not an
    // issue of correctness; use 1 ULP here.
    final double expected = 1000 * S2Earth.getRadiusKm() * Math.PI / 4;
    assertThat(
            S2Earth.getDistanceMeters(S2LatLng.fromDegrees(80, 27), S2LatLng.fromDegrees(55, -153)))
        .isWithin(Math.ulp(expected))
        .of(expected);
  }

  @Test
  public void testGetInitialBearing() {
    S2LatLng equator0 = S2LatLng.fromDegrees(0, 0);
    S2LatLng equator50 = S2LatLng.fromDegrees(0, 50);
    S2LatLng equator100 = S2LatLng.fromDegrees(0, 100);

    // Eastward on Equator.
    assertEquals(S2Earth.getInitialBearing(equator50, equator100), S1Angle.degrees(90));

    // Westward on Equator.
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
    assertThat(0.0)
        .isWithin(1e-7)
        .of(
            S2Earth.getInitialBearing(S2LatLng.fromDegrees(12, 75), S2LatLng.fromDegrees(90, 50))
                .degrees());

    // Towards the south pole.
    assertThat(180.0)
        .isWithin(1e-7)
        .of(
            S2Earth.getInitialBearing(
                    S2LatLng.fromDegrees(-35, 105), S2LatLng.fromDegrees(-90, -120))
                .degrees());

    S2LatLng spain = S2LatLng.fromDegrees(40.4379332, -3.749576);
    S2LatLng japan = S2LatLng.fromDegrees(35.6733227, 139.6403486);
    assertThat(29.2).isWithin(1e-2).of(S2Earth.getInitialBearing(spain, japan).degrees());
    assertThat(-27.2).isWithin(1e-2).of(S2Earth.getInitialBearing(japan, spain).degrees());
  }
}
