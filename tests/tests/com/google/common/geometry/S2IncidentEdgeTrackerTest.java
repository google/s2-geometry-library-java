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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableSet;
import com.google.common.geometry.S2IncidentEdgeTracker.IncidentEdgeKey;
import com.google.common.geometry.S2IncidentEdgeTracker.IncidentEdgeMap;
import com.google.common.geometry.S2Shape.MutableEdge;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Unit tests for {@link S2IncidentEdgeTracker}. This class is mainly tested through the tests of
 * its users, such as S2ValidationQuery, so only basic functionality is covered here.
 */
@RunWith(JUnit4.class)
public class S2IncidentEdgeTrackerTest extends GeometryTestCase {
  private final S2Point a = S2LatLng.fromDegrees(1, 1).toPoint();
  private final S2Point b = S2LatLng.fromDegrees(1, 2).toPoint();
  private final S2Point c = S2LatLng.fromDegrees(2, 1).toPoint();
  private final S2Point d = S2LatLng.fromDegrees(2, 2).toPoint();

  @Test
  public void testVerticesWithTwoIncidentEdges() {
    S2IncidentEdgeTracker tracker = new S2IncidentEdgeTracker();

    // A square shape of four consecutive edges: ab, bc, cd, da. With only two edges incident on
    // each vertex, the tracker won't bother keeping any of them.
    tracker.startShape(0);
    tracker.addEdge(0, e(a, b));
    tracker.addEdge(1, e(b, c));
    tracker.addEdge(2, e(c, d));
    tracker.addEdge(3, e(d, a));
    tracker.finishShape();

    IncidentEdgeMap map = tracker.incidentEdges();
    assertTrue(map.isEmpty());
  }

  @Test
  public void testVerticesWithFourIncidentEdges() {
    S2IncidentEdgeTracker tracker = new S2IncidentEdgeTracker();

    // A square shape with edges in each direction.
    tracker.startShape(0);
    tracker.addEdge(0, e(a, b));
    tracker.addEdge(5, e(b, a));
    tracker.addEdge(1, e(b, c));
    tracker.addEdge(6, e(c, b));
    tracker.addEdge(2, e(c, d));
    tracker.addEdge(7, e(d, c));
    tracker.addEdge(3, e(d, a));
    tracker.addEdge(4, e(a, d));
    tracker.finishShape();

    IncidentEdgeMap map = tracker.incidentEdges();
    assertFalse(map.isEmpty());

    // Edges 0, 5, 3, and 4 are incident on vertex a.
    assertIncidentAt(map, new IncidentEdgeKey(0, a), 0, 5, 4, 3);
    assertIncidentAt(map, new IncidentEdgeKey(0, b), 0, 5, 1, 6);
    assertIncidentAt(map, new IncidentEdgeKey(0, c), 7, 2, 1, 6);
    assertIncidentAt(map, new IncidentEdgeKey(0, d), 2, 4, 3, 7);
  }

  @Test
  public void testTwoShapes() {
    S2IncidentEdgeTracker tracker = new S2IncidentEdgeTracker();

    // A triangle shape around a,b,c with edges in each direction.
    tracker.startShape(0);
    tracker.addEdge(0, e(a, b));
    tracker.addEdge(1, e(b, c));
    tracker.addEdge(2, e(c, a));
    tracker.addEdge(3, e(a, c));
    tracker.addEdge(4, e(c, b));
    tracker.addEdge(5, e(b, a));
    tracker.finishShape();

    // Another triangle shape around b,c,d with edges in each direction.
    tracker.startShape(1);
    tracker.addEdge(0, e(b, c));
    tracker.addEdge(1, e(c, d));
    tracker.addEdge(2, e(d, b));
    tracker.addEdge(3, e(b, d));
    tracker.addEdge(4, e(d, c));
    tracker.addEdge(5, e(c, b));
    tracker.finishShape();

    IncidentEdgeMap map = tracker.incidentEdges();
    assertFalse(map.isEmpty());

    // Check edge map for shape 0.
    assertIncidentAt(map, new IncidentEdgeKey(0, a), 0, 2, 3, 5);
    assertIncidentAt(map, new IncidentEdgeKey(0, b), 0, 1, 4, 5);
    assertIncidentAt(map, new IncidentEdgeKey(0, c), 1, 2, 3, 4);
    assertIncidentAt(map, new IncidentEdgeKey(0, d));

    // Check edge map for shape 1.
    assertIncidentAt(map, new IncidentEdgeKey(1, a));
    assertIncidentAt(map, new IncidentEdgeKey(1, b), 0, 2, 3, 5);
    assertIncidentAt(map, new IncidentEdgeKey(1, c), 0, 1, 4, 5);
    assertIncidentAt(map, new IncidentEdgeKey(1, d), 1, 2, 3, 4);
  }

  // Just like the previous test, but the edges of each shape are added in multiple sets, working
  // through the edges vertex by vertex. Each edge is added twice, once for each endpoint.
  @Test
  public void testTwoShapesAsTwoAddEdgesSequences() {
    S2IncidentEdgeTracker tracker = new S2IncidentEdgeTracker();

    // Shape 0, all edges adjacent to a.
    tracker.startShape(0);
    tracker.addEdge(0, e(a, b));
    tracker.addEdge(5, e(b, a));
    tracker.addEdge(2, e(c, a));
    tracker.addEdge(3, e(a, c));
    tracker.finishShape();

    // Shape 1, all edges adjacent to b.
    tracker.startShape(1);
    tracker.addEdge(0, e(b, c));
    tracker.addEdge(5, e(c, b));
    tracker.addEdge(2, e(d, b));
    tracker.addEdge(3, e(b, d));
    tracker.finishShape();

    // Shape 0, all edges adjacent to b.
    tracker.startShape(0);
    tracker.addEdge(0, e(a, b));
    tracker.addEdge(5, e(b, a));
    tracker.addEdge(1, e(b, c));
    tracker.addEdge(4, e(c, b));
    tracker.finishShape();

    // Shape 1, all edges adjacent to c.
    tracker.startShape(1);
    tracker.addEdge(0, e(b, c));
    tracker.addEdge(5, e(c, b));
    tracker.addEdge(1, e(c, d));
    tracker.addEdge(4, e(d, c));
    tracker.finishShape();

    // Shape 0, all edges adjacent to c.
    tracker.startShape(0);
    tracker.addEdge(1, e(b, c));
    tracker.addEdge(4, e(c, b));
    tracker.addEdge(2, e(c, a));
    tracker.addEdge(3, e(a, c));
    tracker.finishShape();

    // Shape 1, all edges adjacent to d.
    tracker.startShape(1);
    tracker.addEdge(1, e(c, d));
    tracker.addEdge(4, e(d, c));
    tracker.addEdge(2, e(d, b));
    tracker.addEdge(3, e(b, d));
    tracker.finishShape();

    IncidentEdgeMap map = tracker.incidentEdges();
    assertFalse(map.isEmpty());

    // Check edge map for shape 0.
    assertIncidentAt(map, new IncidentEdgeKey(0, a), 0, 2, 3, 5);
    assertIncidentAt(map, new IncidentEdgeKey(0, b), 0, 1, 4, 5);
    assertIncidentAt(map, new IncidentEdgeKey(0, c), 1, 2, 3, 4);
    assertIncidentAt(map, new IncidentEdgeKey(0, d));

    // Check edge map for shape 1.
    assertIncidentAt(map, new IncidentEdgeKey(1, a));
    assertIncidentAt(map, new IncidentEdgeKey(1, b), 0, 2, 3, 5);
    assertIncidentAt(map, new IncidentEdgeKey(1, c), 0, 1, 4, 5);
    assertIncidentAt(map, new IncidentEdgeKey(1, d), 1, 2, 3, 4);
  }

  public void assertIncidentAt(
      IncidentEdgeMap map, IncidentEdgeKey key, Integer... expectedEdgeIds) {
    ImmutableSet<Integer> expected = ImmutableSet.copyOf(expectedEdgeIds);
    if (expected.isEmpty()) {
      assertNull(map.get(key));
    } else {
      assertEquals(expected, map.get(key));
    }
  }

  /** MutableEdge factory */
  public static MutableEdge e(S2Point a, S2Point b) {
    MutableEdge edge = new MutableEdge();
    edge.a = a;
    edge.b = b;
    return edge;
  }
}

