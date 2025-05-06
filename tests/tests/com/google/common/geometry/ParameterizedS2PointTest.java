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
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link ParametrizedS2Point}. */
@RunWith(JUnit4.class)
public class ParameterizedS2PointTest extends GeometryTestCase {
  @Test
  public void testSimple() {
    S2Point point = new S2Point(1, 1e7, 1e-9);
    ParametrizedS2Point p = new ParametrizedS2Point(0.123, point);
    assertExactly(0.123, p.getTime());
    assertEquals(point, p.getPoint());
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testParametrizedS2PointSerialization() {
    doSerializationTest(new ParametrizedS2Point(0.123, new S2Point(1, 1e7, 1e-9)));
  }
}
