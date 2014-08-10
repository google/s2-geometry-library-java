/*
 * Copyright 2014 Google Inc.
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
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;

import java.util.List;

/**
 * Verifies S2ShapeUtil.
 */
@GwtCompatible
public class S2ShapeUtilTest extends GeometryTestCase {
  public void testEdgeAccess() {
    S2EdgeVectorShape shape = new S2EdgeVectorShape();
    final int kNumEdges = 100;
    List<S2Edge> edges = Lists.newArrayList();
    for (int i = 0; i < kNumEdges; ++i) {
      S2Point a = randomPoint();
      S2Point b = randomPoint();
      shape.add(a, b);
      edges.add(new S2Edge(a, b));
    }
    assertEquals(kNumEdges, shape.numEdges());
    MutableEdge edge = new MutableEdge();
    for (int i = 0; i < kNumEdges; ++i) {
      shape.getEdge(i, edge);
      assertEquals(edges.get(i), new S2Edge(edge.getStart(), edge.getEnd()));
    }
  }

  public void testSingletonConstructor() {
    S2Point a = S2Point.X_POS;
    S2Point b = S2Point.Y_POS;
    S2EdgeVectorShape shape = new S2EdgeVectorShape(a, b);
    assertEquals(1, shape.numEdges());
    MutableEdge edge = new MutableEdge();
    shape.getEdge(0, edge);
    assertEquals(a, edge.getStart());
    assertEquals(b, edge.getEnd());
  }
}
