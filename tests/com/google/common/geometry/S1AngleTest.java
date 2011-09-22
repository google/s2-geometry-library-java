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

import junit.framework.TestCase;

public strictfp class S1AngleTest extends TestCase {


  public void testBasic() {
    // Check that the conversion between Pi radians and 180 degrees is exact.
    assertEquals(S1Angle.radians(Math.PI).radians(), Math.PI);
    assertEquals(S1Angle.radians(Math.PI).degrees(), 180.0);
    assertEquals(S1Angle.degrees(180).radians(), Math.PI);
    assertEquals(S1Angle.degrees(180).degrees(), 180.0);

    assertEquals(S1Angle.radians(Math.PI / 2).degrees(), 90.0);

    // Check negative angles.
    assertEquals(S1Angle.radians(-Math.PI / 2).degrees(), -90.0);
    assertEquals(S1Angle.degrees(-45).radians(), -Math.PI / 4);

    // Check that E5/E6/E7 representations work as expected.
    assertEquals(S1Angle.e5(2000000), S1Angle.degrees(20));
    assertEquals(S1Angle.e6(-60000000), S1Angle.degrees(-60));
    assertEquals(S1Angle.e7(750000000), S1Angle.degrees(75));
    assertEquals(S1Angle.degrees(12.34567).e5(), 1234567);
    assertEquals(S1Angle.degrees(12.345678).e6(), 12345678);
    assertEquals(S1Angle.degrees(-12.3456789).e7(), -123456789);
  }
}
