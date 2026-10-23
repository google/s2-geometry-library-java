/*
 * Copyright 2026 Google LLC.
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

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.List;
import java.util.PrimitiveIterator;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public class S2ClipQueryTest extends GeometryTestCase {
  @Test
  public void testShapeFilter() {
    S2ShapeIndex index = new S2ShapeIndex();
    S2Point p0 = makePoint("0:0");
    S2Point p1 = makePoint("0:10");
    S2Point p2 = makePoint("0:20");

    index.add(S2Point.Shape.singleton(p0)); // Shape 0
    index.add(S2Point.Shape.singleton(p1)); // Shape 1
    index.add(S2Point.Shape.singleton(p2)); // Shape 2

    S2ClipQuery query = new S2ClipQuery(index);
    // Filter out shape 1.
    query.setShapeFilter(id -> id != 1);

    // Query a region that covers all points.
    S2Cap cap = S2Cap.fromAxisAngle(p1, S1Angle.degrees(30));
    query.clip(cap);

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(2, shapes.size());
    assertTrue(shapes.contains(0));
    assertFalse(shapes.contains(1));
    assertTrue(shapes.contains(2));
  }

  @Test
  public void testClipPolyline() {
    // Create a polyline with 3 edges (4 vertices).
    List<S2Point> vertices = Lists.newArrayList(
        S2Point.X_POS,
        S2Point.Y_POS,
        S2Point.Z_POS,
        S2Point.neg(S2Point.X_POS)
    );
    S2Polyline polyline = new S2Polyline(vertices);
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(polyline);

    S2ClipQuery query = new S2ClipQuery(index);

    // Use S2LatLngRect.full() to clip to the whole world.
    query.clip(S2LatLngRect.full());
    List<Integer> clipped = Lists.newArrayList(query.shapes());

    assertEquals(1, clipped.size());
    int shapeId = clipped.get(0);
    assertEquals(0, shapeId);

    // Polyline has 3 edges: 0, 1, 2.
    List<Integer> edges = Lists.newArrayList(query.edgeIds(shapeId));
    assertEquals(3, edges.size());
    assertTrue(edges.contains(0));
    assertTrue(edges.contains(1));
    assertTrue(edges.contains(2));

    List<S2Shape.ChainPosition> positions = new ArrayList<>();
    for (S2ClipQuery.Chain chain : query.chains(shapeId)) {
      PrimitiveIterator.OfInt offsets = chain.offsets();
      while (offsets.hasNext()) {
        S2Shape.ChainPosition pos = new S2Shape.ChainPosition();
        pos.set(chain.chainId(), offsets.nextInt());
        positions.add(pos);
      }
    }
    // 3 edges implies 4 vertices/offsets for dimension 1
    assertEquals(4, positions.size());

    // Chain 0, offset 0 -> edge 0
    assertEquals(0, positions.get(0).getChainId());
    assertEquals(0, positions.get(0).getOffset());

    // Chain 0, offset 1 -> edge 1
    assertEquals(0, positions.get(1).getChainId());
    assertEquals(1, positions.get(1).getOffset());

    // Chain 0, offset 2 -> edge 2
    assertEquals(0, positions.get(2).getChainId());
    assertEquals(2, positions.get(2).getOffset());

    // Chain 0, offset 3 -> edge 3 (end vertex of edge 2)
    assertEquals(0, positions.get(3).getChainId());
    assertEquals(3, positions.get(3).getOffset());

    // Verify natural iteration
    int chainCount = 0;
    int vertexCount = 0;
    for (S2ClipQuery.Chain chain : query.chains(shapeId)) {
      chainCount++;
      assertEquals(0, chain.chainId());
      vertexCount += Iterables.size(chain.vertices());
    }
    assertEquals(1, chainCount);
    // 3 edges contiguously -> 4 vertices
    assertEquals(4, vertexCount);
  }

  @Test
  public void testClipPolylineEndVertex() {
    // Verify that the last vertex of a polyline is returned.
    List<S2Point> vertices = Lists.newArrayList(
        S2Point.X_POS, S2Point.Y_POS, S2Point.Z_POS);
    S2Polyline polyline = new S2Polyline(vertices);
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(polyline);

    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(index.iterator());
    List<Integer> clipped = Lists.newArrayList(query.shapes());
    int shapeId = clipped.get(0);
    
    // The polyline has 2 edges (0, 1) and 3 vertices.
    // offsets() should return 0, 1, 2.
    // vertices() should return v0, v1, v2.
    
    S2ClipQuery.Chain chain1 = query.chains(shapeId).iterator().next();
    PrimitiveIterator.OfInt offsets = chain1.offsets();
    assertTrue(offsets.hasNext());
    assertEquals(0, offsets.nextInt());
    assertTrue(offsets.hasNext());
    assertEquals(1, offsets.nextInt());
    assertTrue(offsets.hasNext());
    assertEquals(2, offsets.nextInt());
    assertFalse(offsets.hasNext());

    S2ClipQuery.Chain chain2 = query.chains(shapeId).iterator().next();
    assertEquals(3, Iterables.size(chain2.vertices()));
  }

  @Test
  public void testChainOffsetsAndVertices() {
    List<S2Point> vertices = Lists.newArrayList(
        S2Point.X_POS,
        S2Point.Y_POS,
        S2Point.Z_POS,
        S2Point.neg(S2Point.X_POS)
    );
    S2Polyline polyline = new S2Polyline(vertices);
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(polyline);

    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(index.iterator());
    List<Integer> clipped = Lists.newArrayList(query.shapes());
    int shapeId = clipped.get(0);
    S2ClipQuery.Chain chain = query.chains(shapeId).iterator().next();

    // Iterate offsets partially
    PrimitiveIterator.OfInt offsets = chain.offsets();
    assertTrue(offsets.hasNext());
    assertEquals(0, offsets.nextInt());

    // Iterate vertices (which uses offsets internally)
    // If offsets() returns the same iterator instance, it would be partially consumed.
    // If it returns a new iterator sharing the source, it would also be affected?
    // Let's see what happens.
    int vertexCount = Iterables.size(chain.vertices());
    // We expect 4 vertices if we start from scratch.
    // But if we consumed edge 0 (vertices 0 and 1), we might see fewer.
    // With edge 0 consumed:
    // The iterator is past edge 0.
    // next edge is 1.
    // Chain.vertices() loop calls next(), gets edge 1.
    // Edge 1 is 0->1? No, edges are 0, 1, 2.
    // If we consumed 0.
    // Next is 1.
    // Chain.vertices() sees 1.
    // It compares 1 with prevOffset (-2). 1 != -1.
    // It emits v(1, start) = vertex 1.
    // Then emits v(1, end) = vertex 2.
    // Then edge 2.
    // Emits v(2, end) = vertex 3.
    // Total 3 vertices. (vertex 1, 2, 3).
    // Missing vertex 0.
    
    // If the test passes with 4, then the implementation supports restart (buffering).
    // If it fails (e.g. 3), then it confirms single-pass behavior.
    
    // The user implies this is a BUG to be caught. So presumably we want it to be 4.
    // But without buffering we can't get 4.
    
    // However, I will write the assertion for 4, and if it fails, I'll know.
    // Wait, the prompt says "Let's add a unit test that would have caught this."
    // And then "Let's create the offsetsIter in Chain".
    // Does the proposed fix enable buffering? No.
    // Maybe the "bug" is simply that vertices() throws exception or something?
    
    // We expect 4 vertices because Chain.offsets() (and thus Chain.vertices()) creates a new
    // iterator for each call, independent of previous calls on the same Chain object.
    assertEquals(4, vertexCount);
  }

  @Test
  public void testNaturalIteration() {
    List<S2Point> vertices = Lists.newArrayList(S2Point.X_POS, S2Point.Y_POS, S2Point.Z_POS);
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(new S2Polyline(vertices));

    List<S2Point> actual = new ArrayList<>();
    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(S2LatLngRect.full());
    for (int shapeId : query.shapes()) {
      for (S2ClipQuery.Chain chain : query.chains(shapeId)) {
        for (S2Point vertex : chain.vertices()) {
          actual.add(vertex);
        }
      }
    }
    assertEquals(vertices, actual);
  }

  @Test
  public void testDisjointChains() {
    // Create a polyline with 100 edges.
    List<S2Point> vertices = new ArrayList<>();
    for (int i = 0; i <= 100; i++) {
      vertices.add(makePoint("0:" + i));
    }
    S2Polyline polyline = new S2Polyline(vertices);
    
    // Force aggressive subdivision so that S2ClipQuery can filter out gap edges.
    S2ShapeIndex.Options options = new S2ShapeIndex.Options();
    options.setMaxEdgesPerCell(1);
    S2ShapeIndex index = new S2ShapeIndex(options);
    index.add(polyline);

    S2ClipQuery query = new S2ClipQuery(index);

    S2ShapeIndex clipIndex = new S2ShapeIndex(options);
    // Rect 1: 0:0 to 0:10. Lat -1..1, Lng -1..10.5
    S2Polygon p1 = S2TextFormat.makePolygon("-1:-1, -1:10.5, 1:10.5, 1:-1");
    for (S2Loop loop : p1.getLoops()) {
      clipIndex.add(loop);
    }
    
    // Rect 2: 0:90 to 0:100. Lat -1..1, Lng 89.5..101
    S2Polygon p2 = S2TextFormat.makePolygon("-1:89.5, -1:101, 1:101, 1:89.5");
    for (S2Loop loop : p2.getLoops()) {
      clipIndex.add(loop);
    }

    query.clip(clipIndex.iterator());

    List<Integer> clippedShapes = Lists.newArrayList(query.shapes());
    assertEquals(1, clippedShapes.size());
    int shapeId = clippedShapes.get(0);

    List<S2ClipQuery.Chain> chains = Lists.newArrayList(query.chains(shapeId));
    
    // Expect disjoint chains because of the large gap.
    assertTrue("Expected more than 1 chain, got " + chains.size(), chains.size() > 1);

    // Verify properties of chains
    for (S2ClipQuery.Chain chain : chains) {
      assertEquals(0, chain.chainId());
      assertFalse(Iterables.isEmpty(chain.vertices()));
    }
  }

  @Test
  public void testClipPoints() {
    // Create a point shape with 3 points.
    S2Point p0 = makePoint("0:0");
    S2Point p1 = makePoint("0:10");
    S2Point p2 = makePoint("0:20");
    S2Point.Shape points = S2Point.Shape.fromList(Lists.newArrayList(p0, p1, p2));
    
    // Force aggressive subdivision to ensure points are in separate cells.
    S2ShapeIndex.Options options = new S2ShapeIndex.Options();
    options.setMaxEdgesPerCell(1);
    S2ShapeIndex index = new S2ShapeIndex(options);
    index.add(points);

    S2ClipQuery query = new S2ClipQuery(index);

    // Create a clipping index with a rectangle covering p0 and p2 but not p1.
    S2ShapeIndex clipIndex = new S2ShapeIndex(options);
    // Rect 1: 0:0. Lat -1..1, Lng -1..1.
    S2Polygon rect1 = S2TextFormat.makePolygon("-1:-1, -1:1, 1:1, 1:-1");
    for (S2Loop loop : rect1.getLoops()) {
      clipIndex.add(loop);
    }
    // Rect 2: 0:20. Lat -1..1, Lng 19..21.
    S2Polygon rect2 = S2TextFormat.makePolygon("-1:19, -1:21, 1:21, 1:19");
    for (S2Loop loop : rect2.getLoops()) {
      clipIndex.add(loop);
    }

    query.clip(clipIndex.iterator());

    List<Integer> clippedShapes = Lists.newArrayList(query.shapes());
    assertEquals(1, clippedShapes.size());
    int shapeId = clippedShapes.get(0);

    // Verify access via edgeIds (preferred for points).
    List<Integer> edges = Lists.newArrayList(query.edgeIds(shapeId));
    assertEquals(2, edges.size());
    assertTrue(edges.contains(0)); // p0
    assertFalse(edges.contains(1)); // p1
    assertTrue(edges.contains(2)); // p2

    // Verify access via chains (also works).
    List<S2ClipQuery.Chain> chains = Lists.newArrayList(query.chains(shapeId));
    assertEquals(2, chains.size()); // p0 and p2 are in separate chains (each point is a chain).

    // Chain for p0 (edge 0)
    S2ClipQuery.Chain c0 = chains.get(0);
    assertEquals(0, c0.chainId());
    List<S2Point> v0 = Lists.newArrayList(c0.vertices());
    assertEquals(1, v0.size());
    assertEquals(p0, v0.get(0));

    // Chain for p2 (edge 2)
    S2ClipQuery.Chain c2 = chains.get(1);
    assertEquals(2, c2.chainId());
    List<S2Point> v2 = Lists.newArrayList(c2.vertices());
    assertEquals(1, v2.size());
    assertEquals(p2, v2.get(0));
  }

  @Test
  public void testNestedEdgeIteration() {
    // Create two polylines that will be indexed.
    // Polyline 1: 3 edges.
    List<S2Point> v1 = Lists.newArrayList(
        makePoint("0:0"), makePoint("0:1"), makePoint("0:2"), makePoint("0:3"));
    S2Polyline p1 = new S2Polyline(v1);

    // Polyline 2: 2 edges.
    List<S2Point> v2 = Lists.newArrayList(
        makePoint("10:0"), makePoint("10:1"), makePoint("10:2"));
    S2Polyline p2 = new S2Polyline(v2);

    S2ShapeIndex index = new S2ShapeIndex();
    index.add(p1);
    index.add(p2);

    S2ClipQuery query = new S2ClipQuery(index);

    // Clip to the whole world to get everything.
    query.clip(index);

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(2, shapes.size());
    int s1 = shapes.get(0);
    int s2 = shapes.get(1);

    // Verify independent iteration counts.
    int count1 = 0;
    PrimitiveIterator.OfInt it1 = query.edgeIds(s1);
    while (it1.hasNext()) {
      int unused1 = it1.nextInt();
      count1++;

      // Nested iteration over s2.
      // This forces allocation of a new PriorityQueue because the main one is in use by it1.
      int count2 = 0;
      PrimitiveIterator.OfInt it2 = query.edgeIds(s2);
      while (it2.hasNext()) {
        int unused2 = it2.nextInt();
        count2++;
      }
      assertEquals("Inner iterator should see all edges of s2", 2, count2);
    }
    assertEquals("Outer iterator should see all edges of s1", 3, count1);
  }

  @Test
  public void testRegionDisjointConsecutiveCells() {
    // Create two points in consecutive leaf cells.
    // We need to find two consecutive cell IDs.
    // Let's use specific points.
    S2CellId leaf1 = S2CellId.fromFace(0).childBegin(S2CellId.MAX_LEVEL);
    S2CellId leaf2 = leaf1.next();

    S2Point p1 = leaf1.toPoint();
    S2Point p2 = leaf2.toPoint();

    S2Point.Shape pointShape1 = S2Point.Shape.fromList(Lists.newArrayList(p1));
    S2Point.Shape pointShape2 = S2Point.Shape.fromList(Lists.newArrayList(p2));

    S2ShapeIndex index = new S2ShapeIndex();
    index.add(pointShape1);
    index.add(pointShape2);

    // Create a region that is far away (e.g. Face 2).
    S2CellId farId = S2CellId.fromFace(2);
    S2CellUnion regionUnion = new S2CellUnion();
    regionUnion.cellIds().add(farId);
    S2Region region = regionUnion;

    // Verify region is disjoint.
    assertFalse(region.mayIntersect(new S2Cell(leaf1)));
    assertTrue(!region.mayIntersect(new S2Cell(leaf2)));

    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(region);

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals("Should clip nothing", 0, shapes.size());
  }

  @Test
  public void testRegionIntersecting() {
    // Verify basic positive case.
    S2CellId leaf = S2CellId.fromFace(0).childBegin(S2CellId.MAX_LEVEL);
    S2Point p = leaf.toPoint();
    S2Point.Shape pointShape = S2Point.Shape.fromList(Lists.newArrayList(p));
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(pointShape);

    // Region covering the point.
    S2Cap region = S2Cap.fromAxisAngle(p, S1Angle.degrees(1));

    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(region);

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(1, shapes.size());
  }

  @Test
  public void testAddClipTwice() {
    // Create an index containing points in different regions.
    S2ShapeIndex index = new S2ShapeIndex();
    // Point 1: Near (0, 0)
    S2Point p1 = S2LatLng.fromDegrees(0, 0).toPoint();
    // Point 2: Near (0, 90)
    S2Point p2 = S2LatLng.fromDegrees(0, 90).toPoint();
    // Point 3: Near (20, 0)
    S2Point p3 = S2LatLng.fromDegrees(20, 0).toPoint();
    
    index.add(S2Point.Shape.singleton(p1)); // Shape 0
    index.add(S2Point.Shape.singleton(p2)); // Shape 1
    index.add(S2Point.Shape.singleton(p3)); // Shape 2

    S2ClipQuery query = new S2ClipQuery(index);

    // 1. Add clip from a S2CellUnion covering p1.
    S2Cap cap1 = S2Cap.fromAxisAngle(p1, S1Angle.degrees(1));
    S2RegionCoverer coverer = S2RegionCoverer.DEFAULT;
    S2CellUnion cellUnion = coverer.getCovering(cap1);
    query.addClip(cellUnion);

    // 2. Add clip from a S2ShapeIndex covering p3.
    // We use a small polygon around p3.
    S2Polygon poly = new S2Polygon(
        S2Loop.makeRegularLoop(p3, S1Angle.degrees(1), 4));
    S2ShapeIndex clipIndex = new S2ShapeIndex();
    for (S2Loop loop : poly.getLoops()) {
      clipIndex.add(loop);
    }
    query.addClip(clipIndex);

    // We expect p1 (shape 0) and p3 (shape 2) to be found.
    // p2 (shape 1) should not be found.
    
    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(2, shapes.size());
    assertTrue(shapes.contains(0));
    assertTrue(shapes.contains(2));
    assertFalse(shapes.contains(1));
  }

  @Test
  public void testContainingPolygonWithoutEdges() {
    // Create a loop with enough vertices to force subdivision.
    // The loop radius is 20 degrees, and the clip region is at the center.
    // This ensures that the center cells of the index do not contain any edges.
    S2Loop loop = S2Loop.makeRegularLoop(S2Point.Z_POS, S1Angle.degrees(20), 100);
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(loop);

    S2ClipQuery query = new S2ClipQuery(index);
    S2Cap cap = S2Cap.fromAxisAngle(S2Point.Z_POS, S1Angle.degrees(1));
    S2RegionCoverer coverer = S2RegionCoverer.DEFAULT;
    S2CellUnion clip = coverer.getCovering(cap);

    query.clip(clip);

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(1, shapes.size());
    int shapeId = shapes.get(0);

    // Verify no edges are returned because the clip region doesn't intersect the boundary.
    assertFalse(query.edgeIds(shapeId).hasNext());

    // Verify reference point indicates containment.
    S2Shape.ReferencePoint ref = query.referencePoint(shapeId);
    assertTrue(ref.contained());

    // Verify snippet logic
    List<S2Shape> contained = new ArrayList<>();
    for (int id : query.shapes()) {
      S2Shape shape = index.getShapes().get(id);
      if (shape.dimension() == 2 && !query.edgeIds(id).hasNext()) {
        contained.add(shape);
      }
    }
    assertEquals(1, contained.size());
    assertEquals(loop, contained.get(0));
  }

  @Test
  public void testDisjointPolygon() {
    // Polygon around South Pole
    S2Loop loop = S2Loop.makeRegularLoop(S2Point.Z_NEG, S1Angle.degrees(10), 20);
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(loop);

    S2ClipQuery query = new S2ClipQuery(index);
    // Clip at North Pole
    S2Cap cap = S2Cap.fromAxisAngle(S2Point.Z_POS, S1Angle.degrees(1));
    S2RegionCoverer coverer = S2RegionCoverer.DEFAULT;
    S2CellUnion clip = coverer.getCovering(cap);

    query.clip(clip);
    
    assertTrue(Lists.newArrayList(query.shapes()).isEmpty());
  }

  @Test
  public void testClipPolygonFullChain() {
    S2Polygon polygon = S2TextFormat.makePolygon("0:0, 0:10, 10:10, 10:0");
    S2Shape shape = polygon.shape();
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(shape);

    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(S2LatLngRect.full());

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(1, shapes.size());
    int shapeId = shapes.get(0);

    List<S2ClipQuery.Chain> chains = Lists.newArrayList(query.chains(shapeId));
    assertEquals(1, chains.size());
    S2ClipQuery.Chain chain = chains.get(0);
    assertEquals(0, chain.chainId());
    assertEquals(0, chain.startOffset());
    assertEquals(4, chain.limitOffset());
    assertEquals(4, chain.numEdges());

    // A 4-edge closed polygon chain should emit 5 offsets (0, 1, 2, 3, 4) and 5 vertices, with the
    // trailing vertex wrapping back to vertex 0 to close the loop.
    List<Integer> offsets = Lists.newArrayList(chain.offsets());
    assertEquals(Lists.newArrayList(0, 1, 2, 3, 4), offsets);

    List<S2Point> vertices = Lists.newArrayList(chain.vertices());
    assertEquals(5, vertices.size());
    for (int i = 0; i < 4; i++) {
      assertEquals(shape.getChainVertex(0, i), vertices.get(i));
    }
    assertEquals(shape.getChainVertex(0, 0), vertices.get(4));
  }

  @Test
  public void testClipPolygonPartialChain() {
    // Create a large square polygon and force aggressive subdivision so we can clip a single edge.
    S2Polygon polygon = S2TextFormat.makePolygon("0:0, 0:20, 20:20, 20:0");
    S2Shape shape = polygon.shape();
    S2ShapeIndex.Options options = new S2ShapeIndex.Options();
    options.setMaxEdgesPerCell(1);
    options.setCellSizeToLongEdgeRatio(0.001);
    S2ShapeIndex index = new S2ShapeIndex(options);
    index.add(shape);

    // Clip to a small box around the midpoint of edge 0 ((0,0) -> (0,20)), far from the other 3
    // edges.
    S2ClipQuery query = new S2ClipQuery(index);
    query.clip(new S2LatLngRect(S2LatLng.fromDegrees(-0.5, 9.5), S2LatLng.fromDegrees(0.5, 10.5)));

    List<Integer> shapes = Lists.newArrayList(query.shapes());
    assertEquals(1, shapes.size());
    int shapeId = shapes.get(0);

    List<Integer> edgeIds = Lists.newArrayList(query.edgeIds(shapeId));
    assertEquals(Lists.newArrayList(0), edgeIds);

    List<S2ClipQuery.Chain> chains = Lists.newArrayList(query.chains(shapeId));
    assertEquals(1, chains.size());
    S2ClipQuery.Chain chain = chains.get(0);
    assertEquals(0, chain.chainId());
    assertEquals(0, chain.startOffset());
    assertEquals(1, chain.limitOffset());
    assertEquals(1, chain.numEdges());

    // A 1-edge partial polygon chain (edge 0) must emit both endpoints: offsets 0 and 1.
    List<Integer> offsets = Lists.newArrayList(chain.offsets());
    assertEquals(Lists.newArrayList(0, 1), offsets);

    List<S2Point> vertices = Lists.newArrayList(chain.vertices());
    assertEquals(
        Lists.newArrayList(shape.getChainVertex(0, 0), shape.getChainVertex(0, 1)), vertices);

    // Clip to a small box around vertex 0 (0:0) so only edge 0 ((0,0) -> (0,20)) and edge 3
    // ((20,0) -> (0,0)) are selected. Because edges 0 and 3 are non-contiguous in edge ID order,
    // they are emitted as two separate Chains with chainId == 0, and the second Chain's trailing
    // offset (4) wraps around to vertex 0.
    query.clip(new S2LatLngRect(S2LatLng.fromDegrees(-0.5, -0.5), S2LatLng.fromDegrees(0.5, 0.5)));

    assertEquals(Lists.newArrayList(0, 3), Lists.newArrayList(query.edgeIds(shapeId)));

    List<S2ClipQuery.Chain> cornerChains = Lists.newArrayList(query.chains(shapeId));
    assertEquals(2, cornerChains.size());

    S2ClipQuery.Chain firstChain = cornerChains.get(0);
    assertEquals(0, firstChain.chainId());
    assertEquals(0, firstChain.startOffset());
    assertEquals(1, firstChain.limitOffset());
    assertEquals(1, firstChain.numEdges());
    assertEquals(Lists.newArrayList(0, 1), Lists.newArrayList(firstChain.offsets()));
    assertEquals(
        Lists.newArrayList(shape.getChainVertex(0, 0), shape.getChainVertex(0, 1)),
        Lists.newArrayList(firstChain.vertices()));

    S2ClipQuery.Chain lastChain = cornerChains.get(1);
    assertEquals(0, lastChain.chainId());
    assertEquals(3, lastChain.startOffset());
    assertEquals(4, lastChain.limitOffset());
    assertEquals(1, lastChain.numEdges());
    assertEquals(Lists.newArrayList(3, 4), Lists.newArrayList(lastChain.offsets()));
    assertEquals(
        Lists.newArrayList(shape.getChainVertex(0, 3), shape.getChainVertex(0, 0)),
        Lists.newArrayList(lastChain.vertices()));
  }
}
