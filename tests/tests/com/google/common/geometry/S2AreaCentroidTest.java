/*
 * Copyright 2023 Google Inc.
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
import java.util.Comparator;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2AreaCentroid}. */
@RunWith(JUnit4.class)
public class S2AreaCentroidTest extends GeometryTestCase {
  @Test
  public void testSimple() {
    S2Point point = new S2Point(0.1, 10.3, 7.5);
    S2AreaCentroid ac = new S2AreaCentroid(5.0, point);
    assertExactly(5.0, ac.getArea());
    assertEquals(point, ac.getCentroid());
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2AreaCentroidSerialization() {
    Comparator<S2AreaCentroid> comparator = (a, b) -> {
      if (a.equals(b)) {
        return 0;
      }
      // Only equality matters for the serialization test.
      return Double.compare(a.getArea(), b.getArea());
    };

    doSerializationTest(new S2AreaCentroid(5.0, new S2Point(0.1, 10.3, 7.5)), comparator);
  }
}
