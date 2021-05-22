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

import com.google.common.annotations.GwtCompatible;

@GwtCompatible
public class S2ContainsVertexQueryTest extends GeometryTestCase {
  public void testUndetermined() {
    S2ContainsVertexQuery q = new S2ContainsVertexQuery(makePoint("1:2"));
    q.addOutgoing(makePoint("3:4"));
    q.addIncoming(makePoint("3:4"));
    assertEquals(0, q.containsSign());
  }

  public void testContainedWithDuplicates() {
    // The S2.ortho reference direction points approximately due west.
    // Containment is determined by the unmatched edge immediately clockwise.
    S2ContainsVertexQuery q = new S2ContainsVertexQuery(makePoint("0:0"));
    q.addIncoming(makePoint("3:-3"));
    q.addOutgoing(makePoint("1:-5"));
    q.addOutgoing(makePoint("2:-4"));
    q.addIncoming(makePoint("1:-5"));
    assertEquals(1, q.containsSign());
  }

  public void testNotContainedWithDuplicates() {
    // The S2.ortho reference direction points approximately due west.
    // Containment is determined by the unmatched edge immediately clockwise.
    S2ContainsVertexQuery q = new S2ContainsVertexQuery(makePoint("1:1"));
    q.addOutgoing(makePoint("1:-5"));
    q.addIncoming(makePoint("2:-4"));
    q.addOutgoing(makePoint("3:-3"));
    q.addIncoming(makePoint("1:-5"));
    assertEquals(-1, q.containsSign());
  }

  public void testMatchesLoopContainment() {
    // Check that the containment function defined is compatible with S2Loop
    // (which at least currently does not use this class).
    S2Loop loop = S2Loop.makeRegularLoop(makePoint("89:-179"), S1Angle.degrees(10), 1000);
    for (int i = 1; i <= loop.numVertices(); ++i) {
      S2ContainsVertexQuery q = new S2ContainsVertexQuery(loop.vertex(i));
      q.addIncoming(loop.vertex(i - 1));
      q.addOutgoing(loop.vertex(i + 1));
      assertEquals(q.containsSign() > 0, loop.contains(loop.vertex(i)));
    }
  }
}
