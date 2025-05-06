/*
 * Copyright 2022 Google Inc.
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
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertThrows;

import com.google.common.annotations.GwtIncompatible;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for S2Edge. */
@RunWith(JUnit4.class)
public final class S2EdgeTest extends GeometryTestCase {

  private final List<S2Point> vertices = S2TextFormat.parsePointsOrDie("0:0, 1:1");

  @Test
  public void testBasic() {
    final S2Edge edge = new S2Edge(vertices.get(0), vertices.get(1));
    assertSame(vertices.get(0), edge.getStart());
    assertSame(vertices.get(1), edge.getEnd());
    assertFalse(edge.hasInterior());
    assertFalse(edge.containsOrigin());
    assertEquals(1, edge.numChains());
  }

  @Test
  public void testGetChainMethods() {
    final S2Edge edge = new S2Edge(vertices.get(0), vertices.get(1));
    assertEquals(0, edge.getChainStart(0));
    assertEquals(1, edge.getChainLength(0));

    S2Shape.MutableEdge mutableEdge = new S2Shape.MutableEdge();
    edge.getChainEdge(0, 0, mutableEdge);
    assertSame(vertices.get(0), mutableEdge.getStart());
    assertSame(vertices.get(1), mutableEdge.getEnd());

    // getChainVertex can get the start and end of the single edge in the single chain.
    assertSame(vertices.get(0), edge.getChainVertex(0, 0));
    assertSame(vertices.get(1), edge.getChainVertex(0, 1));
  }

  @Test
  public void testGetChainArgumentValidation() {
    final S2Edge edge = new S2Edge(vertices.get(0), vertices.get(1));

    // S2Edge has a single chain 0 so chain 1 must fail a precondition check.
    assertThrows("Index out of range", IndexOutOfBoundsException.class,
        () -> edge.getChainStart(1));

    assertThrows("Index out of range", IndexOutOfBoundsException.class,
        () -> edge.getChainLength(1));

    S2Shape.MutableEdge mutableEdge = new S2Shape.MutableEdge();
    assertThrows("Index out of range", IndexOutOfBoundsException.class,
        () -> edge.getChainEdge(1, 0, mutableEdge));
    // The chain has length 0 so offset 1 must also fail.
    assertThrows("Index out of range", IndexOutOfBoundsException.class,
        () -> edge.getChainEdge(0, 1, mutableEdge));

    // Offsets 0 and 1 are valid offsets into the single chain, as tested above, but anything higher
    // is invalid.
    assertThrows("Index out of range", IndexOutOfBoundsException.class,
        () -> assertEquals(vertices.get(0), edge.getChainVertex(1, 0)));
    assertThrows("Index out of range", IndexOutOfBoundsException.class,
        () -> assertEquals(vertices.get(0), edge.getChainVertex(0, 2)));
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2EdgeSerialization() {
    doSerializationTest(new S2Edge(new S2Point(0.1, 0.2, 0.3), new S2Point(0.5, 0.6, -100)));
  }
}
