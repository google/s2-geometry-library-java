/*
 * Copyright 2024 Google Inc.
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

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit Tests for {@link S2Shape}. */
@RunWith(JUnit4.class)
public class S2ShapeTest extends GeometryTestCase {
  @Test
  public void testNextPrevEdgePointDoesNotWrap() {
    var index = S2TextFormat.makeIndexOrDie("1:1 | 2:2 ##");
    S2Shape shape = index.getShapes().get(0);

    // Points have one chain per point so we should always get -1.
    assertEquals(-1, shape.prevEdgeWrap(0));
    assertEquals(-1, shape.nextEdgeWrap(0));

    assertEquals(-1, shape.prevEdgeWrap(1));
    assertEquals(-1, shape.nextEdgeWrap(1));
  }

  @Test
  public void testNextPrevEdgeOpenPolylineDoesNotWrap() {
    var index = S2TextFormat.makeIndexOrDie("# 1:1, 2:2, 3:3 #");
    S2Shape shape = index.getShapes().get(0);

    // Open polylines should not wrap around.
    assertEquals(-1, shape.prevEdgeWrap(0));
    assertEquals(1, shape.nextEdgeWrap(0));

    assertEquals(0, shape.prevEdgeWrap(1));
    assertEquals(-1, shape.nextEdgeWrap(1));
  }

  @Test
  public void testNextPrevEdgeClosedPolylineWraps() {
    var index = S2TextFormat.makeIndexOrDie("# 0:0, 1:1, 0:2, -1:1, 0:0 #");
    S2Shape shape = index.getShapes().get(0);

    // Closed polylines should wrap around.
    assertEquals(3, shape.prevEdgeWrap(0));
    assertEquals(1, shape.nextEdgeWrap(0));

    assertEquals(2, shape.prevEdgeWrap(3));
    assertEquals(0, shape.nextEdgeWrap(3));
  }

  @Test
  public void testNextPrevEdgePolygonWraps() {
    var index = S2TextFormat.makeIndexOrDie("## 0:0, 1:1, 0:2, -1:1");
    S2Shape shape = index.getShapes().get(0);

    // Polygons should always wrap.
    assertEquals(3, shape.prevEdgeWrap(0));
    assertEquals(1, shape.nextEdgeWrap(0));

    assertEquals(2, shape.prevEdgeWrap(3));
    assertEquals(0, shape.nextEdgeWrap(3));
  }
}