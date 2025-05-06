/*
 * Copyright 2019 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makePoint;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.OptionalInt;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link S2ContainsVertexQuery}. */
@RunWith(JUnit4.class)
public class S2ContainsVertexQueryTest extends GeometryTestCase {
  @Test
  public void testUndetermined() {
    S2ContainsVertexQuery q = new S2ContainsVertexQuery(makePoint("1:2"));
    q.addOutgoing(makePoint("3:4"));
    q.addIncoming(makePoint("3:4"));
    assertEquals(0, q.safeContainsSign().getAsInt());
    assertFalse(q.duplicateEdges());
  }

  @Test
  public void testContainedWithDuplicates() {
    // The S2.ortho reference direction points approximately due west.
    // Containment is determined by the unmatched edge immediately clockwise.
    S2ContainsVertexQuery q = new S2ContainsVertexQuery(makePoint("0:0"));
    q.addIncoming(makePoint("3:-3"));
    q.addOutgoing(makePoint("1:-5"));
    q.addOutgoing(makePoint("2:-4"));
    q.addIncoming(makePoint("1:-5"));
    assertEquals(1, q.containsSign());
    assertFalse(q.duplicateEdges());

    // Incoming and outgoing edges to 1:-5 cancel, so one more incoming isn't a duplicate.
    q.addIncoming(makePoint("1:-5"));
    assertFalse(q.duplicateEdges());

    // 3:3 has only been seen once incoming, another incoming is a duplicate.
    q.addIncoming(makePoint("3:-3"));
    assertTrue(q.duplicateEdges());
  }

  @Test
  public void testNotContainedWithDuplicates() {
    // The S2.ortho reference direction points approximately due west.
    // Containment is determined by the unmatched edge immediately clockwise.
    S2ContainsVertexQuery q = new S2ContainsVertexQuery(makePoint("1:1"));
    q.addOutgoing(makePoint("1:-5"));
    q.addIncoming(makePoint("2:-4"));
    q.addOutgoing(makePoint("3:-3"));
    q.addIncoming(makePoint("1:-5"));
    OptionalInt result = q.safeContainsSign();
    assertTrue(result.isPresent());
    assertEquals(-1, result.getAsInt());

    // Incoming and outgoing edges to 1:-5 cancel, so one more outgoing isn't a duplicate.
    q.addOutgoing(makePoint("1:-5"));
    result = q.safeContainsSign();
    assertTrue(result.isPresent());

    // 3:3 has only been seen once outgoing, another outgoing is a duplicate.
    q.addOutgoing(makePoint("3:-3"));
    result  = q.safeContainsSign();
    assertFalse(result.isPresent());
  }

  @Test
  // Tests that S2ContainsVertexQuery is compatible with S2Predicates.angleContainsVertex().
  public void testCompatibleWithAngleContainsVertex() {
    List<S2Point> points = S2Loop.makeRegularVertices(
        makePoint("89:1"), S1Angle.degrees(5), 10);
    for (int i = 0; i < points.size(); ++i) {
      S2Point a = points.get(i);
      S2Point b = points.get((i + 1) % points.size());
      S2Point c = points.get((i + 2) % points.size());
      S2ContainsVertexQuery query = new S2ContainsVertexQuery(b);
      query.addIncoming(a);
      query.addOutgoing(c);
      assertEquals(query.containsSign() > 0, S2Predicates.angleContainsVertex(a, b, c));
      assertFalse(query.duplicateEdges());
    }
  }

  // Tests compatibility with S2Predicates.angleContainsVertex() for a degenerate edge.
  @Test
  public void testCompatibleWithAngleContainsVertexDegenerate() {
    S2Point a = new S2Point(1, 0, 0);
    S2Point b = new S2Point(0, 1, 0);
    S2ContainsVertexQuery query = new S2ContainsVertexQuery(b);
    query.addIncoming(a);
    query.addOutgoing(a);
    assertEquals(query.containsSign() > 0, S2Predicates.angleContainsVertex(a, b, a));
    assertFalse(query.duplicateEdges());
  }
}
