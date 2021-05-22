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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeAspect.ChainAspect.Multi;
import com.google.common.geometry.S2ShapeAspect.ChainAspect.Simple;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Exercises all mixes of the S2ShapeAspects components. */
@GwtIncompatible("Insufficient support for generics")
public class S2ShapeAspectTest extends GeometryTestCase {
  private static S2Point A = makePoint("0:0");
  private static S2Point B = makePoint("0:1");
  private static S2Point C = makePoint("1:1");
  private static S2Point D = makePoint("1:0");
  
  // Open, Simple, Array and Packed
  
  public void testOpenSimpleEmpty() {
    S2Point[][] chains = {{}};
    S2Point[][] edges = {};
    check(chains, edges, new LineSimpleArray(chains[0]));
    check(chains, edges, new LineSimplePacked(chains[0]));
  }
  
  public void testOpenSimpleSingleton() {
    S2Point[][] chains = {{A}};
    S2Point[][] edges = {};
    check(chains, edges, new LineSimpleArray(chains[0]));
    check(chains, edges, new LineSimplePacked(chains[0]));
  }
  
  public void testOpenSimpleTriangle() {
    S2Point[][] chains = {{A, B, C}};
    S2Point[][] edges = {{A, B}, {B, C}};
    check(chains, edges, new LineSimpleArray(chains[0]));
    check(chains, edges, new LineSimplePacked(chains[0]));
  }
  
  private static class LineSimpleArray extends Simple.Array implements Line {
    LineSimpleArray(S2Point ... vertices) {
      super(Arrays.asList(vertices));
    }
  }
  
  private static class LineSimplePacked extends Simple.Packed implements Line {
    LineSimplePacked(S2Point ... vertices) {
      super(Arrays.asList(vertices));
    }
  }
  
  // Closed, Simple, Array and Packed
  
  public void testClosedSimpleEmpty() {
    S2Point[][] chains = {{}};
    S2Point[][] edges = {};
    check(chains, edges, new RingSimpleArray(chains[0]));
    check(chains, edges, new RingSimplePacked(chains[0]));
  }
  
  public void testClosedSimpleSingleton() {
    S2Point[][] chains = {{A}};
    S2Point[][] edges = {{A, A}};
    check(chains, edges, new RingSimpleArray(chains[0]));
    check(chains, edges, new RingSimplePacked(chains[0]));
  }
  
  public void testClosedSimpleTriangle() {
    S2Point[][] chains = {{A, B, C}};
    S2Point[][] edges = {{A, B}, {B, C}, {C, A}};
    check(chains, edges, new RingSimpleArray(chains[0]));
    check(chains, edges, new RingSimplePacked(chains[0]));
  }
  
  private static class RingSimpleArray extends Simple.Array implements Ring {
    RingSimpleArray(S2Point ... vertices) {
      super(Arrays.asList(vertices));
    }
  }
  
  private static class RingSimplePacked extends Simple.Packed implements Ring {
    RingSimplePacked(S2Point ... vertices) {
      super(Arrays.asList(vertices));
    }
  }
  
  // Open, Multi, Array and Packed
  
  public void testOpenMultiEmpty() {
    S2Point[][] chains = {};
    S2Point[][] edges = {};
    check(chains, edges, new LineMultiArray(chains));
    check(chains, edges, new LineMultiPacked(chains));
  }
  
  public void testOpenMultiEmptyChain() {
    // JUnit's assertThrows would be nice, but it isn't supported on GWT.
    try {
      new LineMultiArray(new S2Point[][]{{}});
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("Must have at least 1 edge.", e.getMessage());
    }
    try {
      new LineMultiPacked(new S2Point[][]{{A, B}, {}});
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("Must have at least 1 edge.", e.getMessage());
    }
  }
  
  public void testOpenMultiOneChain() {
    S2Point[][] chains = {{A, B, C}};
    S2Point[][] edges = {{A, B}, {B, C}};
    check(chains, edges, new LineMultiArray(chains));
    check(chains, edges, new LineMultiPacked(chains));
  }
  
  public void testOpenMultiTwoChains() {
    S2Point[][] chains = {{A, B, C}, {D, C}};
    S2Point[][] edges = {{A, B}, {B, C}, {D, C}};
    check(chains, edges, new LineMultiArray(chains));
    check(chains, edges, new LineMultiPacked(chains));
  }
  
  public void testOpenMultiPartialChain() {
    // Note that while we don't allow empty chains, degenerate chains have no edges.
    S2Point[][] chains = {{D}, {A, B, C}, {D}, {D, C}, {D}};
    S2Point[][] edges = {{A, B}, {B, C}, {D, C}};
    check(chains, edges, new LineMultiArray(chains));
    check(chains, edges, new LineMultiPacked(chains));
  }
  
  private static class LineMultiArray extends Multi.Array implements Line {
    LineMultiArray(S2Point[][] chains) {
      super(toLists(chains));
    }
  }
  
  private static class LineMultiPacked extends Multi.Packed implements Line {
    LineMultiPacked(S2Point[][] chains) {
      super(toLists(chains));
    }
  }
  
  // Closed, Multi, Array and Packed
  
  public void testClosedMultiEmpty() {
    S2Point[][] chains = {};
    S2Point[][] edges = {};
    check(chains, edges, new RingMultiArray(chains));
    check(chains, edges, new RingMultiPacked(chains));
  }
  
  public void testClosedMultiEmptyChain() {
    S2Point[][] chains = {{}};
    S2Point[][] edges = {};
    check(chains, edges, new RingMultiArray(chains));
    check(chains, edges, new RingMultiPacked(chains));
  }
  
  public void testClosedMultiOneChain() {
    S2Point[][] chains = {{A, B, C}};
    S2Point[][] edges = {{A, B}, {B, C}, {C, A}};
    check(chains, edges, new RingMultiArray(chains));
    check(chains, edges, new RingMultiPacked(chains));
  }
  
  public void testClosedMultiTwoChainsWithEmpties() {
    S2Point[][] chains = {{}, {A, B, C}, {}, {D, C}, {}};
    S2Point[][] edges = {{A, B}, {B, C}, {C, A}, {D, C}, {C, D}};
    check(chains, edges, new RingMultiArray(chains));
    check(chains, edges, new RingMultiPacked(chains));
  }
  
  public void testClosedMultiPartialChain() {
    S2Point[][] chains = {{D}, {A, B, C}, {D}, {D, C}, {D}};
    S2Point[][] edges = {{D, D}, {A, B}, {B, C}, {C, A}, {D, D}, {D, C}, {C, D}, {D, D}};
    check(chains, edges, new RingMultiArray(chains));
    check(chains, edges, new RingMultiPacked(chains));
  }
  
  private static class RingMultiArray extends Multi.Array implements Ring {
    RingMultiArray(S2Point[][] chains) {
      super(toLists(chains));
    }
  }
  
  public void testClosedMultiPacked() {
    
  }
  
  private static class RingMultiPacked extends Multi.Packed implements Ring {
    RingMultiPacked(S2Point[][] chains) {
      super(toLists(chains));
    }
  }
  
  // Snapped
  
  public void testSnapped() {
    // Do a sanity test of the [open,closed]x[simple,multi] mixes, having snapped ABC to abc.
    S2CellId ia = S2CellId.fromPoint(A);
    S2CellId ib = S2CellId.fromPoint(B);
    S2CellId ic = S2CellId.fromPoint(C);
    S2Point a = ia.toPoint();
    S2Point b = ib.toPoint();
    S2Point c = ic.toPoint();
    S2CellId[] simpleCells = {ia, ib, ic};
    S2Point[][] simpleChains = {{a, b, c}};
    S2Point[][] simpleOpenEdges = {{a, b}, {b, c}};
    S2Point[][] simpleClosedEdges = {{a, b}, {b, c}, {c, a}};
    check(simpleChains, simpleOpenEdges, new LineSimpleSnapped(simpleCells));
    check(simpleChains, simpleClosedEdges, new RingSimpleSnapped(simpleCells));
    S2CellId[][] multiCells = {{ia, ib, ic}, {ic, ib, ia}};
    S2Point[][] multiChains = {{a, b, c}, {c, b, a}};
    S2Point[][] multiOpenEdges = {{a, b}, {b, c}, {c, b}, {b, a}};
    S2Point[][] multiClosedEdges = {{a, b}, {b, c}, {c, a}, {c, b}, {b, a}, {a, c}};
    check(multiChains, multiOpenEdges, new LineMultiSnapped(multiCells));
    check(multiChains, multiClosedEdges, new RingMultiSnapped(multiCells));
  }
  
  private static class LineSimpleSnapped extends Simple.Snapped implements Line {
    LineSimpleSnapped(S2CellId ... vertices) {
      super(Arrays.asList(vertices));
    }
  }
  
  private static class RingSimpleSnapped extends Simple.Snapped implements Ring {
    RingSimpleSnapped(S2CellId ... vertices) {
      super(Arrays.asList(vertices));
    }
  }
  
  private static class LineMultiSnapped extends Multi.Snapped implements Line {
    LineMultiSnapped(S2CellId[][] chains) {
      super(toLists(chains));
    }
  }

  private static class RingMultiSnapped extends Multi.Snapped implements Ring {
    RingMultiSnapped(S2CellId[][] chains) {
      super(toLists(chains));
    }
  }
  
  private interface Line extends Undef, S2ShapeAspect.EdgeAspect.Open {
    @Override default int dimension() {
      return 1;
    }
  }
  
  private interface Ring extends Undef, S2ShapeAspect.EdgeAspect.Closed {
    @Override default int dimension() {
      return 2;
    }
  }

  /** Leaves dimension methods effectively unimplemented. */
  private interface Undef extends S2Shape, S2ShapeAspect.TopoAspect {
    @Override default boolean hasInterior() {
      throw new UnsupportedOperationException();
    }
    @Override default boolean containsOrigin() {
      throw new UnsupportedOperationException();
    }
  }
  
  /** Verifies the expected vertices, chains, and edges. */
  private static void check(
      S2Point[][] chains,
      S2Point[][] edges,
      S2ShapeAspect.Mixed shape) {
    List<S2Point> vertices = ImmutableList.copyOf(Iterables.concat(toLists(chains)));
    assertEquals("Unexpected vertices", vertices, shape.vertices());
    assertEquals("Unexpected chains", toLists(chains), shape.chains());
    assertEquals("Unexpected edges", toLists(edges), edges(shape));
  }
  
  private static <T> List<List<T>> toLists(T[][] edges) {
    return Lists.transform(Arrays.asList(edges), Arrays::asList);
  }
  
  private static List<List<S2Point>> edges(S2ShapeAspect.Mixed shape) {
    List<List<S2Point>> edges = new ArrayList<>();
    MutableEdge edge = new MutableEdge();
    for (int i = 0; i < shape.numEdges(); i++) {
      shape.getEdge(i, edge);
      S2Point a = edge.getStart();
      S2Point b = edge.getEnd();
      int chainId = shape.chainId(i);
      int offset = i - shape.getChainStart(chainId);
      shape.getChainEdge(chainId, offset, edge);
      assertEquals(a, edge.a);
      assertEquals(b, edge.b);
      edges.add(Arrays.asList(a, b));
    }
    return edges;
  }
}
