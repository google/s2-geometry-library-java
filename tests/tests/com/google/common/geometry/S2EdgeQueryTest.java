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

import static com.google.common.geometry.S2Projections.PROJ;
import static com.google.common.geometry.TestDataGenerator.makePoint;
import static com.google.common.geometry.TestDataGenerator.makePolyline;
import static java.lang.Math.pow;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.geometry.S2EdgeQuery.Edges;
import com.google.common.geometry.S2EdgeUtil.FaceSegment;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public strictfp class S2EdgeQueryTest extends GeometryTestCase {
  public S2Point perturbAtDistance(S1Angle distance, S2Point a0, S2Point b0) {
    S2Point p = S2EdgeUtil.interpolateAtDistance(distance, a0, b0);
    if (data.oneIn(2)) {
      p =
          new S2Point(
              data.oneIn(2) ? p.x + Platform.ulp(p.x) : p.x - Platform.ulp(p.x),
              data.oneIn(2) ? p.y + Platform.ulp(p.y) : p.y - Platform.ulp(p.y),
              data.oneIn(2) ? p.z + Platform.ulp(p.z) : p.z - Platform.ulp(p.z));
    }
    return p.normalize();
  }

  /**
   * Generates sub-edges of some edge ({@code a0}, {@code b0}). The length of the sub-edges is
   * distributed exponentially over a large range, and the endpoints may be slightly perturbed to
   * one side of ({@code a0}, {@code b0}) or the other.
   */
  public S2ShapeUtil.S2EdgeVectorShape getPerturbedSubEdges(S2Point a0, S2Point b0, int count) {
    S2ShapeUtil.S2EdgeVectorShape edges = new S2ShapeUtil.S2EdgeVectorShape();
    a0 = a0.normalize();
    b0 = b0.normalize();
    double length0 = a0.angle(b0);
    for (int i = 0; i < count; ++i) {
      double length = length0 * pow(1e-15, data.nextDouble());
      double offset = (length0 - length) * data.nextDouble();
      edges.add(
          perturbAtDistance(S1Angle.radians(offset), a0, b0),
          perturbAtDistance(S1Angle.radians(offset + length), a0, b0));
    }
    return edges;
  }

  /**
   * Generates edges whose center is randomly chosen from the given S2Cap, and whose length is
   * randomly chosen up to {@code maxLength}.
   */
  public S2ShapeUtil.S2EdgeVectorShape getCapEdges(S2Cap centerCap, S1Angle maxLength, int count) {
    S2ShapeUtil.S2EdgeVectorShape edges = new S2ShapeUtil.S2EdgeVectorShape();
    for (int i = 0; i < count; ++i) {
      S2Point center = data.samplePoint(centerCap);
      S2Cap edgeCap = S2Cap.fromAxisAngle(center, S1Angle.radians(0.5 * maxLength.radians()));
      S2Point p = data.samplePoint(edgeCap);
      edges.add(p, center.mul(2 * p.dotProd(center)).sub(p).normalize());
    }
    return edges;
  }

  private void checkAllCrossings(S2ShapeUtil.S2EdgeVectorShape edges) {
    MutableEdge tempEdge = new MutableEdge();
    S2ShapeUtil.S2EdgeVectorShape shape = new S2ShapeUtil.S2EdgeVectorShape();
    for (int i = 0; i < edges.size(); ++i) {
      shape.add(edges.get(i).getStart(), edges.get(i).getEnd());
    }
    // Force more subdivision than usual to make the test more challenging.
    S2ShapeIndex.Options options = new S2ShapeIndex.Options();
    options.setMaxEdgesPerCell(1);
    S2ShapeIndex index = new S2ShapeIndex(options);
    index.add(shape);
    // To check that candidates are being filtered reasonably, we count the total number of
    // candidates and the total number of edge pairs that either intersect or are very close to
    // intersecting.
    int numCandidates = 0;
    int numNearbyPairs = 0;
    for (int i = 0; i < edges.size(); ++i) {
      S2Point a = edges.get(i).getStart();
      S2Point b = edges.get(i).getEnd();
      S2EdgeQuery query = new S2EdgeQuery(index);
      // Shape id has to be 0 because only one shape was inserted.
      S2EdgeQuery.Edges candidates = query.getCandidates(a, b, index.shapes.get(0));

      Map<S2Shape, S2EdgeQuery.Edges> edgeMap = query.getCandidates(a, b);
      assertEquals(1, edgeMap.size());
      assertTrue(edgeMap.containsKey(shape));
      List<Integer> candidatesList = edgesToList(candidates);
      List<Integer> edgeMapList = edgesToList(edgeMap.get(shape));
      assertEquals(edgeMapList, candidatesList);
      assertFalse(candidatesList.isEmpty());

      // Now check the actual candidates.
      // Assert that 'candidates' is sorted.
      assertTrue(Ordering.natural().isOrdered(candidatesList));
      assertTrue(candidatesList.get(0) >= 0);
      assertTrue(candidatesList.get(candidatesList.size() - 1) < shape.numEdges());
      numCandidates += candidatesList.size();
      List<Integer> missingCandidates = Lists.newArrayList();
      List<Integer> expectedCrossings = Lists.newArrayList();
      for (int j = 0; j < shape.numEdges(); ++j) {
        shape.getEdge(j, tempEdge);
        S2Point c = tempEdge.getStart();
        S2Point d = tempEdge.getEnd();
        if (S2EdgeUtil.robustCrossing(a, b, c, d) >= 0) {
          expectedCrossings.add(j);
          ++numNearbyPairs;
          if (!candidatesList.contains(i)) {
            missingCandidates.add(i);
          }
        } else {
          double maxDist = PROJ.maxDiag.getValue(S2CellId.MAX_LEVEL);
          if (S2EdgeUtil.getDistance(a, c, d).radians() < maxDist
              || S2EdgeUtil.getDistance(b, c, d).radians() < maxDist
              || S2EdgeUtil.getDistance(c, a, b).radians() < maxDist
              || S2EdgeUtil.getDistance(d, a, b).radians() < maxDist) {
            ++numNearbyPairs;
          }
        }
      }
      assertTrue(missingCandidates.isEmpty());

      // Test that GetCrossings() returns only the actual crossing edges.
      assertEquals(expectedCrossings, edgesToList(query.getCrossings(a, b, shape)));

      // Verify that the second version of GetCrossings returns the same result.
      Map.Entry<S2Shape, Edges> result =
          Iterables.getOnlyElement(query.getCrossings(a, b).entrySet());
      assertEquals(shape, result.getKey());
      assertEquals(expectedCrossings, edgesToList(result.getValue()));
    }

    // There is nothing magical about this particular ratio; this check exists to catch changes that
    // dramatically increase the number of candidates.
    assertTrue(numCandidates < 3 * numNearbyPairs);
  }

  /** Returns a list of all edges in {@code candidates}. */
  private static List<Integer> edgesToList(S2EdgeQuery.Edges candidates) {
    List<Integer> list = Lists.newArrayList();
    while (!candidates.isEmpty()) {
      list.add(candidates.nextEdge());
    }
    return list;
  }

  /**
   * Tests edges that lie in the plane of one of the S2 cube edges. Such edges may lie on the
   * boundary between two cube faces, or pass through a cube vertex, or follow a 45 degree diagonal
   * across a cube face toward its center.
   *
   * <p>This test is sufficient to demonstrate that padding the cell boundaries is necessary for
   * correctness. (It fails if S2ShapeIndex.CELL_PADDING is set to zero.)
   */
  public void testPerturbedCubeEdges() {
    for (int iter = 0; iter < 10; ++iter) {
      int face = data.nextInt(6);
      double scale = pow(1e-15, data.nextDouble());
      R2Vector uv = new R2Vector(2 * data.nextInt(2) - 1, 2 * data.nextInt(2) - 1);
      S2Point a0 = S2Projections.faceUvToXyz(face, R2Vector.mul(uv, scale));
      S2Point b0 = a0.sub(S2Projections.getNorm(face).mul(2));
      // This test is slow because *every* crossing test needs S2Predicates.expensiveSign().
      checkAllCrossings(getPerturbedSubEdges(a0, b0, 30));
    }
  }

  /**
   * Tests edges that lie in the plane of one of the S2 cube face axes. These edges are special
   * because one coordinate is zero, and they lie on the boundaries between the immediate child
   * cells of the cube face.
   */
  public void testPerturbedCubeFaceAxes() {
    for (int iter = 0; iter < 5; ++iter) {
      int face = data.nextInt(6);
      double scale = pow(1e-15, data.nextDouble());
      S2Point axis = S2Projections.getUVWAxis(face, data.nextInt(2));
      S2Point a0 = axis.mul(scale).add(S2Projections.getNorm(face));
      S2Point b0 = axis.mul(scale).sub(S2Projections.getNorm(face));
      checkAllCrossings(getPerturbedSubEdges(a0, b0, 30));
    }
  }

  /**
   * Tests a random collection of edges near the S2 cube vertex where the Hilbert curve starts and
   * ends.
   */
  public void testCapEdgesNearCubeVertex() {
    S2Point p = new S2Point(-1, -1, 1).normalize();
    checkAllCrossings(
        getCapEdges(S2Cap.fromAxisAngle(p, S1Angle.radians(1e-3)), S1Angle.radians(1e-4), 1000));
  }

  public void testCollinearEdgesOnCellBoundaries() {
    // 9 * 8 / 2 = 36 edges
    int numEdgeIntervals = 8;
    for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
      S2Cell cell = new S2Cell(data.getRandomCellId(level));
      int v1 = data.nextInt(4);
      int v2 = (v1 + 1) & 3;
      S2Point p1 = cell.getVertexRaw(v1);
      S2Point p2 = cell.getVertexRaw(v2);
      S2Point delta = p2.sub(p1).div(numEdgeIntervals);
      S2ShapeUtil.S2EdgeVectorShape edges = new S2ShapeUtil.S2EdgeVectorShape();
      for (int i = 0; i <= numEdgeIntervals; ++i) {
        for (int j = 0; j < i; ++j) {
          edges.add(p1.add(delta.mul(i)).normalize(), p1.add(delta.mul(j)).normalize());
        }
      }
      checkAllCrossings(edges);
    }
  }

  public void testGetCandidatesMultipleShapes() {
    S2Point a = data.getRandomPoint();
    S2Point b = data.getRandomPoint();
    S2ShapeIndex.Options options = new S2ShapeIndex.Options();
    options.setMaxEdgesPerCell(1);
    S2ShapeIndex index = new S2ShapeIndex(options);
    int numShapes = 10;
    int numEdgesPerShape = 10;
    for (int s = 0; s < numShapes; ++s) {
      S2ShapeUtil.S2EdgeVectorShape shape = new S2ShapeUtil.S2EdgeVectorShape();
      for (int i = 0; i < numEdgesPerShape; ++i) {
        shape.add(data.getRandomPoint(), data.getRandomPoint());
      }
      index.add(shape);
    }

    S2EdgeQuery edgeQuery = new S2EdgeQuery(index);
    Map<S2Shape, S2EdgeQuery.Edges> edgeMap = edgeQuery.getCandidates(a, b);
    assertEquals(numShapes, edgeMap.size());
    for (int i = 0; i < index.shapes.size(); ++i) {
      assertTrue(edgeMap.containsKey(index.shapes.get(i)));
    }
  }

  public void testPolylineCrossings() {
    // Three zig-zag lines near the equator.
    List<S2Polyline> polylines =
        ImmutableList.of(
            makePolyline("0:0, 2:1, 0:2, 2:3, 0:4, 2:5, 0:6"),
            makePolyline("1:0, 3:1, 1:2, 3:3, 1:4, 3:5, 1:6"),
            makePolyline("2:0, 4:1, 2:2, 4:3, 2:4, 4:5, 2:6"));
    checkPolylineCrossings(polylines, makePoint("1:0"), makePoint("1:4"));
    checkPolylineCrossings(polylines, makePoint("5:5"), makePoint("6:6"));
  }

  /** Verifies the example in the class javadoc, with a few extras. */
  private void checkPolylineCrossings(List<S2Polyline> polylines, S2Point a0, S2Point a1) {
    S2ShapeIndex index = new S2ShapeIndex();
    for (S2Polyline line : polylines) {
      index.add(line);
    }
    S2EdgeQuery query = new S2EdgeQuery(index);
    for (Map.Entry<S2Shape, Edges> entry : query.getCrossings(a0, a1).entrySet()) {
      // Shapes with no crossings should be filtered out by this method.
      S2Polyline polyline = (S2Polyline) entry.getKey();
      assertFalse(edgesToList(query.getCrossings(a0, a1, polyline)).isEmpty());
      for (Edges edges = entry.getValue(); !edges.isEmpty(); ) {
        int j = edges.nextEdge();
        S2Point b0 = polyline.vertex(j);
        S2Point b1 = polyline.vertex(j + 1);
        assertTrue(S2EdgeUtil.robustCrossing(a0, a1, b0, b1) >= 0);
      }
    }
    // Also test that every edge with a crossing is in the results exactly once.
    for (S2Polyline line : polylines) {
      for (int i = 0; i < line.numVertices() - 1; i++) {
        if (S2EdgeUtil.robustCrossing(a0, a1, line.vertex(i), line.vertex(i + 1)) >= 0) {
          assertEquals(1, Collections.frequency(edgesToList(query.getCrossings(a0, a1, line)), i));
        }
      }
    }
  }

  /**
   * Tests that {@link S2EdgeQuery#getCells(S2Point, R2Vector, S2Point, R2Vector, S2PaddedCell,
   * List)} returns the correct value.
   */
  public void testGetCellsRegression() {
    int numTests = 10;
    for (int n = 0; n < numTests; ++n) {
      S2Point a = data.getRandomPoint();
      S2Cap cap = S2Cap.fromAxisAngle(a, S1Angle.degrees(2 * (n + 1)));
      S2Point b = data.samplePoint(cap);
      S2Loop loop = S2Loop.makeRegularLoop(a, S1Angle.radians(a.angle(b)), 10);
      S2ShapeIndex.Options options = new S2ShapeIndex.Options();
      options.setMaxEdgesPerCell(5);
      S2ShapeIndex index = new S2ShapeIndex(options);
      index.add(loop);
      S2EdgeQuery edgeQuery = new S2EdgeQuery(index);

      FaceSegment[] segments = FaceSegment.allFaces();
      int numSegments = S2EdgeUtil.getFaceSegments(a, b, segments);
      for (int i = 0; i < numSegments; ++i) {
        S2PaddedCell pCell = new S2PaddedCell(S2CellId.fromFace(segments[i].face), 0);
        List<S2ShapeIndex.Cell> cells = Lists.newArrayList();
        edgeQuery.getCells(a, segments[i].a, b, segments[i].b, pCell, cells);
        assertFalse(cells.isEmpty());
      }
    }
  }
}
