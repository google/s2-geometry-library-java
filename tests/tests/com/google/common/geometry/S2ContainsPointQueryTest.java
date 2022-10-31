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

import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static com.google.common.geometry.TestDataGenerator.makeIndex;
import static com.google.common.geometry.TestDataGenerator.makePoint;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2ContainsPointQuery.EdgeVisitor;
import com.google.common.geometry.S2ContainsPointQuery.Options;
import com.google.common.geometry.S2ContainsPointQuery.ShapeVisitor;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.primitives.Ints;
import java.util.ArrayList;
import java.util.List;

public class S2ContainsPointQueryTest extends GeometryTestCase {
  public void testVertexModelOpen() {
    S2ShapeIndex index = makeIndex("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
    S2ContainsPointQuery q = new S2ContainsPointQuery(index, Options.OPEN);
    assertFalse(q.contains(makePoint("0:0")));
    assertFalse(q.contains(makePoint("0:1")));
    assertFalse(q.contains(makePoint("0:2")));
    assertFalse(q.contains(makePoint("0:5")));
    assertFalse(q.contains(makePoint("0:7")));
    assertFalse(q.contains(makePoint("2:6")));
    assertTrue(q.contains(makePoint("1:6")));
    assertFalse(q.contains(makePoint("10:10")));
    S2ContainsPointQuery q2 = new S2ContainsPointQuery(index, Options.OPEN);
    assertFalse(q2.shapeContains(index.getShapes().get(1), makePoint("1:6")));
    assertTrue(q2.shapeContains(index.getShapes().get(2), makePoint("1:6")));
    assertFalse(q2.shapeContains(index.getShapes().get(2), makePoint("0:5")));
    assertFalse(q2.shapeContains(index.getShapes().get(2), makePoint("0:7")));
  }

  public void testVertexModelSemiOpen() {
    S2ShapeIndex index = makeIndex("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
    S2ContainsPointQuery q = new S2ContainsPointQuery(index, Options.SEMI_OPEN);
    assertFalse(q.contains(makePoint("0:0")));
    assertFalse(q.contains(makePoint("0:1")));
    assertFalse(q.contains(makePoint("0:2")));
    assertFalse(q.contains(makePoint("0:5")));
    assertTrue(q.contains(makePoint("0:7"))); // Contained vertex.
    assertFalse(q.contains(makePoint("2:6")));
    assertTrue(q.contains(makePoint("1:6")));
    assertFalse(q.contains(makePoint("10:10")));
    S2ContainsPointQuery q2 = new S2ContainsPointQuery(index, Options.SEMI_OPEN);
    assertFalse(q2.shapeContains(index.getShapes().get(1), makePoint("1:6")));
    assertTrue(q2.shapeContains(index.getShapes().get(2), makePoint("1:6")));
    assertFalse(q2.shapeContains(index.getShapes().get(2), makePoint("0:5")));
    assertTrue(q2.shapeContains(index.getShapes().get(2), makePoint("0:7")));
  }

  public void testVertexModelClosed() {
    S2ShapeIndex index = makeIndex("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
    S2ContainsPointQuery q = new S2ContainsPointQuery(index, Options.CLOSED);
    assertTrue(q.contains(makePoint("0:0")));
    assertTrue(q.contains(makePoint("0:1")));
    assertTrue(q.contains(makePoint("0:2")));
    assertTrue(q.contains(makePoint("0:5")));
    assertTrue(q.contains(makePoint("0:7")));
    assertTrue(q.contains(makePoint("2:6")));
    assertTrue(q.contains(makePoint("1:6")));
    assertFalse(q.contains(makePoint("10:10")));
    S2ContainsPointQuery q2 = new S2ContainsPointQuery(index, Options.CLOSED);
    assertFalse(q2.shapeContains(index.getShapes().get(1), makePoint("1:6")));
    assertTrue(q2.shapeContains(index.getShapes().get(2), makePoint("1:6")));
    assertTrue(q2.shapeContains(index.getShapes().get(2), makePoint("0:5")));
    assertTrue(q2.shapeContains(index.getShapes().get(2), makePoint("0:7")));
  }

  public void testGetContainingShapes() {
    // Also tests shapeContains().
    int numVerticesPerLoop = 10;
    S1Angle maxLoopRadius = kmToAngle(10);
    S2Cap centerCap = S2Cap.fromAxisAngle(data.getRandomPoint(), maxLoopRadius);
    S2ShapeIndex index = new S2ShapeIndex();
    for (int i = 0; i < 100; ++i) {
      S2Point center = data.samplePoint(centerCap);
      S1Angle radius = S1Angle.radians(maxLoopRadius.radians() * data.nextDouble());
      index.add(S2Loop.makeRegularLoop(center, radius, numVerticesPerLoop));
    }
    S2ContainsPointQuery query = new S2ContainsPointQuery(index);
    for (int i = 0; i < 100; i++) {
      S2Point p = data.samplePoint(centerCap);
      List<S2Shape> expected = new ArrayList<>();
      for (S2Shape shape : index.getShapes()) {
        S2Loop loop = (S2Loop) shape;
        if (loop.contains(p)) {
          assertTrue(query.shapeContains(shape, p));
          expected.add(shape);
        } else {
          assertFalse(query.shapeContains(shape, p));
        }
      }
      assertEquals(expected, ImmutableList.copyOf(query.getContainingShapes(p)));
      final List<S2Shape> visited = new ArrayList<>();
      final S2Shape end = expected.get(expected.size() / 2);
      query.visitContainingShapes(
          p,
          new ShapeVisitor() {
            @Override
            public boolean apply(S2Shape shape) {
              if (shape == end) {
                return false;
              } else {
                visited.add(shape);
                return true;
              }
            }
          });
      assertEquals(expected.subList(0, expected.size() / 2), visited);
    }
  }

  public void testVisitIncidentEdges() {
    S2ShapeIndex index = makeIndex("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2");
    expectIncidentEdgeIds(index, makePoint("0:0"), 0, 0);
    expectIncidentEdgeIds(index, makePoint("1:1"), 0, 1, 1, 0);
    expectIncidentEdgeIds(index, makePoint("1:2"), 1, 0, 2, 0, 2, 2);
    expectIncidentEdgeIds(index, makePoint("1:3"), 2, 0, 2, 1);
    expectIncidentEdgeIds(index, makePoint("2:2"), 2, 1, 2, 2);
  }

  private void expectIncidentEdgeIds(S2ShapeIndex index, S2Point p, int... expectedPairs) {
    List<Integer> actual = new ArrayList<>();
    S2ContainsPointQuery q = new S2ContainsPointQuery(index);
    assertTrue(
        q.visitIncidentEdges(
            p,
            new EdgeVisitor() {
              @Override
              public boolean test(S2Shape shape, int edgeId, S2Point a, S2Point b) {
                actual.add(index.getShapes().indexOf(shape));
                actual.add(edgeId);
                return true;
              }
            },
            new MutableEdge()));
    assertEquals(Ints.asList(expectedPairs), actual);
  }
}
