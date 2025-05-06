/*
 * Copyright 2013 Google Inc.
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

/** Verifies R2Vector. */
@RunWith(JUnit4.class)
public class R2VectorTest extends GeometryTestCase {
  @Test
  public void testConstructors() {
    double[] coordinates = new double[2];
    coordinates[0] = 1.5;
    coordinates[1] = 2.5;
    R2Vector v = new R2Vector(coordinates);
    assertEquals(v, new R2Vector(1.5, 2.5));
    assertExactly(1.5, v.x());
    assertExactly(2.5, v.y());
  }

  @Test
  public void testOrtho() {
    assertEquals(new R2Vector(1, 1), new R2Vector(1, -1).ortho());
    assertEquals(new R2Vector(1, -1), new R2Vector(-1, -1).ortho());
    assertEquals(new R2Vector(-1, -1), new R2Vector(-1, 1).ortho());
    assertEquals(new R2Vector(1, 1), new R2Vector(1, -1).ortho());
  }

  @Test
  public void testAdd() {
    assertEquals(new R2Vector(5, 5), R2Vector.add(new R2Vector(4, 3), new R2Vector(1, 2)));
    assertEquals(new R2Vector(5, 5), new R2Vector(4, 3).add(new R2Vector(1, 2)));
  }

  @Test
  public void testSub() {
    assertEquals(new R2Vector(3, 1), R2Vector.sub(new R2Vector(4, 3), new R2Vector(1, 2)));
    assertEquals(new R2Vector(3, 1), new R2Vector(4, 3).sub(new R2Vector(1, 2)));
  }

  @Test
  public void testMul() {
    assertAlmostEquals(12.0, R2Vector.mul(new R2Vector(4, 3), 3.0).x());
    assertAlmostEquals(9.0, R2Vector.mul(new R2Vector(4, 3), 3.0).y());
    assertAlmostEquals(12.0, new R2Vector(4, 3).mul(3.0).x());
    assertAlmostEquals(9.0, new R2Vector(4, 3).mul(3.0).y());
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testR2VectorSerialization() {
    doSerializationTest(new R2Vector(0.025, 1e6));
  }
}
