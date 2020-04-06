/*
 * Copyright 2017 Google Inc.
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
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/** Verifies S2PointRegion. */
@GwtCompatible
public class S2PointRegionTest extends GeometryTestCase {
  public void testS2Region() {
    S2PointRegion pointRegion = new S2PointRegion(1, 0, 0);

    S2Cap expectedCapBound = S2Cap.fromAxisHeight(pointRegion.getPoint(), 0);
    assertEquals(expectedCapBound, pointRegion.getCapBound());

    S2LatLng ll = new S2LatLng(pointRegion.getPoint());
    S2LatLngRect expectedRect = new S2LatLngRect(ll, ll);
    assertEquals(expectedRect, pointRegion.getRectBound());

    // The leaf cell containing a point is still much larger than the point region.
    S2Cell cell = new S2Cell(pointRegion.getPoint());
    assertFalse(pointRegion.contains(cell));
    assertTrue(pointRegion.mayIntersect(cell));
  }

  public void testSerialization() throws IOException {
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    new S2PointRegion(S2Point.X_NEG).encode(bos);
    assertEquals(
        new S2PointRegion(S2Point.X_NEG),
        S2PointRegion.decode(new ByteArrayInputStream(bos.toByteArray())));
  }
}
