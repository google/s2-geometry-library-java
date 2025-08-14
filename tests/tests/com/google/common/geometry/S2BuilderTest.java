/*
 * Copyright 2023 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.makePolygonOrDie;
import static com.google.common.geometry.S2TextFormat.makePolylineOrDie;
import static com.google.common.geometry.S2TextFormat.parsePointsOrDie;
import static java.lang.Double.NaN;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Stopwatch;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.S2CellIdSnapFunction;
import com.google.common.geometry.S2BuilderUtil.GraphClone;
import com.google.common.geometry.primitives.IdSetLexicon.IdSet;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for S2Builder. */
@RunWith(JUnit4.class)
public final class S2BuilderTest extends GeometryTestCase {
  // Iteration multiplier for randomized tests.
  private static final int ITERATION_MULTIPLIER = 1;

  /** Verify that S2Builder.Builder.intersectionTolerance depends on splitCrossingEdges. */
  @Test
  public void testBuilderIntersectionTolerance() {
    S2Builder.Builder builder = new S2Builder.Builder();
    S1Angle nonSplitTolerance = builder.intersectionTolerance();
    builder.setSplitCrossingEdges(true);
    S1Angle splitTolerance = builder.intersectionTolerance();
    assertTrue(splitTolerance.radians() > nonSplitTolerance.radians());
  }

  /** Tests adding a polygon and getting the same polygon back again. */
  @Test
  public void testAddShape() {
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    S2Polygon input = makePolygonOrDie("0:0, 0:5, 5:5, 5:0; 1:1, 1:4, 4:4, 4:1");
    builder.addShape(input.linearSearchShape());

    assertBuildOk(builder);
    assertPolygonsEqual(input, layer.getPolygon());
  }

  /**
   * When IdentitySnapFunction is used (i.e., no special requirements on vertex locations), check
   * that vertices closer together than the snap radius are merged together.
   */
  @Test
  public void testSimpleVertexMerging() {
    S1Angle snapRadius = S1Angle.degrees(0.5);
    S2Builder.Builder options =
        new S2Builder.Builder(new S2BuilderSnapFunctions.IdentitySnapFunction(snapRadius));
    S2Builder builder = options.build();

    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    S2Polygon input =
        makePolygonOrDie("0:0, 0.2:0.2, 0.1:0.2, 0.1:0.9, 0:1, 0.1:1.1, 0.9:1, 1:1, 1:0.9");
    builder.addPolygon(input);

    assertBuildOk(builder);
    S2Polygon expected = makePolygonOrDie("0:0, 0:1, 1:0.9");
    assertPolygonsApproxEqual(expected, layer.getPolygon(), snapRadius);
  }

  /**
   * When S2CellIdSnapFunction is used, check that all output vertices are the centers of S2CellIds
   * at the specified level.
   */
  @Test
  public void testSimpleS2CellIdSnapping() {
    int level = S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.degrees(1));
    S2CellIdSnapFunction snapFunction = new S2CellIdSnapFunction(level);
    S2Builder.Builder options = new S2Builder.Builder(snapFunction);
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    S2Polygon input = makePolygonOrDie("2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3");
    builder.addPolygon(input);

    assertBuildOk(builder);
    S2Polygon output = layer.getPolygon();
    assertEquals(1, output.numLoops());
    S2Loop loop = output.loop(0);
    for (int i = 0; i < loop.numVertices(); ++i) {
      assertEquals(S2CellId.fromPoint(loop.vertex(i)).parent(level).toPoint(), loop.vertex(i));
    }
    assertPolygonsApproxEqual(input, output, snapFunction.snapRadius());
  }

  /**
   * When IntLatLngSnapFunction is used, check that all output vertices are snapped to lat/lngs
   * which are integers at the specified exponent.
   */
  @Test
  public void testSimpleIntLatLngSnapping() {
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0)); // E0 coords
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    S2Polygon input =
        makePolygonOrDie(
            "2.01:2.09, 3.24:4.49, 1.78:6.25, 3.51:5.49, 6.11:6.11, "
                + "5.22:3.88, 5.55:2.49, 4.49:2.51");
    S2Polygon expected = makePolygonOrDie("2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3");
    builder.addPolygon(input);

    assertBuildOk(builder);
    S2Polygon output = layer.getPolygon();
    assertEquals(1, output.numLoops());
    assertPolygonsEqual(expected, output);
  }

  /** Check that chains of closely spaced vertices do not collapse into a single vertex. */
  @Test
  public void testVerticesMoveLessThanSnapRadius() {
    S1Angle snapRadius = S1Angle.degrees(1);
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
    S2Builder builder = options.build();

    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    // The spacing between input vertices is about 2*pi*20/1000 = 0.125 degrees. The output
    // vertices are spaced between 1 and 2 degrees apart; the average spacing is about 1.33
    // degrees.
    S2Polygon input =
        new S2Polygon(S2Loop.makeRegularLoop(new S2Point(1, 0, 0), S1Angle.degrees(20), 1000));
    builder.addPolygon(input);

    assertBuildOk(builder);
    S2Polygon output = layer.getPolygon();
    assertEquals(1, output.numLoops());
    assertTrue(output.loop(0).numVertices() >= 90);
    assertTrue(output.loop(0).numVertices() <= 100);
    assertTrue(output.boundaryNear(input, snapRadius.radians()));
  }

  /**
   * Check that edges are separated from non-incident vertices by at least
   * minEdgeVertexSeparation(). This requires adding new vertices (not present in the input) in some
   * cases.
   */
  @Test
  public void testMinEdgeVertexSeparation() {
    // The input is a skinny right triangle with two legs of length 10 and 1, and whose diagonal
    // is subdivided into 10 short edges. Using a snap radius of 0.5, about half of the long leg
    // is snapped onto the diagonal (which causes that part of the polygon to be removed). But the
    // real problem is that the remaining part of the long leg gets too close to the remaining
    // vertices on the diagonal, i.e. it would violate the minimum edge-vertex separation
    // guarantee. S2Builder handles this by creating at least one vertex along the original long
    // leg, to keep the snapped edge far enough away from the diagonal.
    S2Polygon input =
        makePolygonOrDie("0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 6:.4, 7:.3, 8:.2, 9:.1, 10:0");
    S2Polygon expected =
        makePolygonOrDie("0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 4.00021862252687:0");
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(0.5)));
    S2Builder builder = options.build();

    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    builder.addPolygon(input);

    assertBuildOk(builder);
    assertPolygonsApproxEqual(expected, layer.getPolygon(), S1Angle.radians(1e-15));
  }

  /**
   * Test that when long edges are snapped to nearby sites, appropriate extra vertices are added so
   * that the snapped edge chain stays within options.maxEdgeDeviation() of the original edge.
   *
   * <p>We do this by constructing an edge AB with two nearly antipodal vertices, then adding
   * another vertex C slightly perturbed from A. Frequently AB will snap to C producing the edge
   * chain ACB, however edge CB diverges very far from the original edge AB. This requires S2Builder
   * to add a new vertex M near the middle of AB so that ACMB stays close everywhere to AB.
   */
  @Test
  public void testMaxEdgeDeviation() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSplitCrossingEdges(true); // Requires edges to be snapped.
    options.setIdempotent(false);
    S2Builder builder = options.build();

    // Even though we are using the default snap radius of zero, the edge snap radius is
    // S2EdgeUtil.INTERSECTION_ERROR because splitCrossingEdges() is true.
    assertExactly(options.edgeSnapRadius().radians(), S2EdgeUtil.INTERSECTION_ERROR);
    S1Angle maxDeviation = options.maxEdgeDeviation();
    int kIters = 50 * ITERATION_MULTIPLIER;

    // Test cases are constructed randomly, so not all tests are effective (i.e. AB might not snap
    // to the perturbed vertex C). Here we keep track of the number of effective tests.
    int numEffective = 0;
    for (int iter = 0; iter < kIters; ++iter) {
      data.setSeed(iter + 1); // Easier to reproduce a specific case.
      S2LaxPolylineLayer layer = new S2LaxPolylineLayer();
      builder.startLayer(layer);
      S2Point a = data.getRandomPoint();

      // B is slightly perturbed from -A, and C is slightly perturbed from A. The perturbation
      // amount is not critical except that it should not be too much larger than edgeSnapRadius()
      // (otherwise AB will rarely snap to C) and it should also not be so small that A and C are
      // likely to be the same point.
      S2Point b = a.neg().add(data.getRandomPoint().mul(5e-16)).normalize();
      S2Point c = a.add(data.getRandomPoint().mul(5e-16)).normalize();
      builder.addEdge(a, b);
      builder.forceVertex(c);

      assertBuildOk(builder);
      S2LaxPolylineShape output = layer.getPolyline();
      int n = output.numVertices();
      assertEquals(a, output.vertex(0));
      assertEquals(b, output.vertex(n - 1));
      for (int i = 0; i + 1 < n; ++i) {
        assertTrue(
            "Iteration " + iter
                + ": a=" + S2TextFormat.toString(a)
                + ", ab midpoint=" + S2TextFormat.toString(S2EdgeUtil.interpolate(0.5, a, b))
                + ", output polyline = " + S2TextFormat.toString(output)
                + " polyline edge from vertex " + i + " to " + (i + 1),
            S2EdgeUtil.isEdgeBNearEdgeA(
                a, b, output.vertex(i), output.vertex(i + 1), maxDeviation));
      }
      if (n > 2) {
        ++numEffective;
      }
    }
    // We require at least 20% of the test cases to be successful.
    assertGreaterOrEqual(numEffective * 5, kIters);
  }

  /**
   * Previously, S2Builder's chooseAllVerticesAsSites method could fail to deduplicate input
   * vertices that were within the same leaf S2Cell. This test checks that this no longer happens.
   */
  @Test
  public void testSnappingTinyLoopBugfix() {
    // Build a tiny loop as a convex hull around a single point.  All three points are within the
    // same level 30 S2Cell.
    var query = new S2ConvexHullQuery();
    query.addPoint(S2LatLng.fromDegrees(4.56, 1.23).toPoint());
    S2Loop loop = query.getConvexHull();

    // Use an identify snap function with a small snap radius. None of the points will be
    // moved by snapping, as although they're close together, they're not within the snap radius.
    S2Builder builder =
        new S2Builder.Builder(
            new S2BuilderSnapFunctions.IdentitySnapFunction(S1Angle.radians(1e-15)))
            .build();

    var pointLayer = new S2PointVectorLayer();
    var lineLayer =
        new S2PolylineVectorLayer(
            new S2PolylineVectorLayer.Options(S2Builder.EdgeType.DIRECTED).setValidate(true));
    var polygonLayer =
        new S2PolygonLayer(
            new S2PolygonLayer.Options(S2Builder.EdgeType.DIRECTED).setValidate(true));

    builder.startLayer(pointLayer);
    builder.startLayer(lineLayer);
    builder.startLayer(polygonLayer);
    builder.addLoop(loop);
    builder.addIsFullPolygonPredicate(unused -> false);

    var s2Error = new S2Error();
    boolean success = builder.build(s2Error);
    assertTrue(s2Error.toString(), success);

    assertEquals(0, pointLayer.getPointVector().size());
    assertEquals(0, lineLayer.getPolylines().size());
    assertEquals(1, polygonLayer.getPolygon().numLoops());
  }

  /**
   * This test checks that when vertices are closer together than minVertexSeparation(), they are
   * are snapped together, even when options.idempotent() is true.
   */
  @Test
  public void testIdempotencySnapsInadequatelySeparatedVertices() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(1.0)));
    S2Builder builder = options.build();
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:0, 0:0.9, 0:2"));

    assertBuildOk(builder);
    String expected = "0:0, 0:2";
    assertEquals(expected, S2TextFormat.toString(layer.getPolyline()));
  }

  /**
   * This test checks that even when the snap radius is zero, identical vertices are snapped
   * together.
   */
  @Test
  public void testIdempotencySnapsIdenticalVerticesWithZeroSnapRadius() {
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:1, 1:0"));
    builder.addPolyline(makePolylineOrDie("0:0, 0:1"));
    builder.addEdge(makePointOrDie("0:1"), makePointOrDie("0:1"));
    builder.addPolyline(makePolylineOrDie("1:0, 0:0"));

    assertBuildOk(builder);
    String expected = "0:0, 0:1, 1:0";

    S2Polygon output = layer.getPolygon();
    assertEquals(expected, S2TextFormat.toString(output));
  }

  /**
   * This test checks that identical vertices are snapped together even when the snap radius is zero
   * and options.splitCrossingEdges() is true.
   */
  @Test
  public void testIdempotencySnapsIdenticalVerticesWithZeroSnapRadiusEdgeSplitting() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSplitCrossingEdges(true);
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:1, 1:0"));
    builder.addPolyline(makePolylineOrDie("0:0, 0:1"));
    builder.addEdge(makePointOrDie("0:1"), makePointOrDie("0:1"));
    builder.addPolyline(makePolylineOrDie("1:0, 0:0"));

    assertBuildOk(builder);
    String expected = "0:0, 0:1, 1:0";

    S2Polygon output = layer.getPolygon();
    assertEquals(expected, S2TextFormat.toString(output));
  }

  /**
   * When idempotency is requested, no snapping is done unless S2Builder finds at least one vertex
   * or edge that could not be the output of a previous snapping operation. This test checks that
   * S2Builder detects vertices that are not at a valid location returned by the given snap
   * function.
   */
  @Test
  public void testIdempotencySnapsUnsnappedVertices() {
    // In this example we snap two vertices to integer lat/lng coordinates. The two vertices are
    // far enough apart (more than minVertexSeparation) so that they might be the result of a
    // previous snapping operation, but one of the two vertices does not have integer lat/lng
    // coordinates. We use internal knowledge of how snap sites are chosen (namely, that candidates
    // are considered in S2CellId order) to construct two different cases, one where the snapped
    // vertex is processed first and one where the unsnapped vertex is processed first. This
    // exercises two different code paths.
    IntLatLngSnapFunction snapFunction = new IntLatLngSnapFunction(0);
    assertGreaterOrEqual(snapFunction.snapRadius(), S1Angle.degrees(0.7));
    assertLessOrEqual(snapFunction.minVertexSeparation(), S1Angle.degrees(0.35));
    S2Builder.Builder options = new S2Builder.Builder(snapFunction);
    S2Builder builder = options.build();

    // In this example, the snapped vertex (0, 0) is processed first and is selected as a Voronoi
    // site (i.e., output vertex). The second vertex is closer than minVertexSeparation(), therefore
    // it is snapped to the first vertex and the polyline becomes degenerate.
    S2Point a = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Point b = S2LatLng.fromDegrees(0.01, 0.6).toPoint();
    assertTrue(S2CellId.fromPoint(a).compareTo(S2CellId.fromPoint(b)) < 0);
    S2Polyline input1 = new S2Polyline(ImmutableList.of(a, b));
    S2PolylineLayer layer1 = new S2PolylineLayer();
    builder.startLayer(layer1);
    builder.addPolyline(input1);

    assertBuildOk(builder);
    assertEquals("0:0, 0:1", S2TextFormat.toString(layer1.getPolyline()));

    // In this example the unsnapped vertex is processed first and is snapped to(0, 0). The second
    // vertex is further than snapRadius() away, so it is also snapped (which does nothing) and is
    // left at (0, 1).
    S2Point c = S2LatLng.fromDegrees(0.01, 0.4).toPoint();
    S2Point d = S2LatLng.fromDegrees(0, 1).toPoint();
    assertTrue(S2CellId.fromPoint(c).compareTo(S2CellId.fromPoint(d)) < 0);
    S2Polyline input2 = new S2Polyline(ImmutableList.of(c, d));
    S2PolylineLayer layer2 = new S2PolylineLayer();
    builder.startLayer(layer2);
    builder.addPolyline(input2);

    assertBuildOk(builder);
    assertEquals("0:0, 0:1", S2TextFormat.toString(layer2.getPolyline()));
  }

  /**
   * When idempotency is requested, no snapping is done unless S2Builder finds at least one vertex
   * or edge that could not be the output of a previous snapping operation. This test checks that
   * S2Builder detects edges that are too close to vertices even when the snap radius is very small
   * (e.g., S2EdgeUtil.INTERSECTION_ERROR).
   *
   * <p>Previously, the C++ implementation of S2Builder used a conservative approximation to decide
   * whether edges were too close to vertices; unfortunately this meant that when the snap radius
   * was very small then no snapping would be done at all, because even an edge/vertex distance of
   * zero was considered far enough apart.
   *
   * <p>This tests that the current code (which uses exact predicates) handles this situation
   * correctly (i.e., that an edge separated from a non-incident vertex by a distance of zero cannot
   * be the output of a previous snapping operation).
   */
  @Test
  public void testIdempotencySnapsEdgesWithTinySnapRadius() {
    S2Builder.Builder options =
        new S2Builder.Builder(
            new IdentitySnapFunction(S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR)));
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setDuplicateEdges(DuplicateEdges.MERGE);
    S2Builder builder = options.build();
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:0, 0:10"));
    builder.addPolyline(makePolylineOrDie("0:5, 0:7"));

    assertBuildOk(builder);
    List<S2Polyline> output = layer.getPolylines();
    assertEquals(1, output.size());
    assertEquals("0:0, 0:5, 0:7, 0:10", S2TextFormat.toString(output.get(0)));
  }

  /**
   * When idempotency is requested, no snapping is done unless S2Builder finds at least one vertex
   * or edge that could not be the output of a previous snapping operation. This test checks that
   * when an edge is further away than minEdgeVertexSeparation() then no snapping is done.
   */
  @Test
  public void testIdempotencyDoesNotSnapAdequatelySeparatedEdges() {
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0));
    options.setIdempotent(true); // Test fails if this is "false".
    S2Builder builder = options.build();

    S2PolygonLayer layer1 = new S2PolygonLayer();
    builder.startLayer(layer1);
    builder.addPolygon(makePolygonOrDie("1.49:0, 0:2, 0.49:3"));
    assertBuildOk(builder);
    S2Polygon output1 = layer1.getPolygon();
    String expected = "1:0, 0:2, 0:3";
    assertEquals(expected, S2TextFormat.toString(output1));

    S2PolygonLayer layer2 = new S2PolygonLayer();
    builder.startLayer(layer2);
    builder.addPolygon(output1);
    assertBuildOk(builder);
    S2Polygon output2 = layer2.getPolygon();
    assertEquals(expected, S2TextFormat.toString(output2));
  }

  /**
   * Verify that even when the splitCrossingEdges() option is used with a snap radius of zero, edges
   * are snapped to nearby vertices (those within a distance of S2EdgeUtil.INTERSECTION_ERROR). This
   * is necessary for correctness whenever new intersection vertices are created (since such
   * vertices are only within S2EdgeUtil.INTERSECTION_ERROR of the true intersection point), and it
   * is necessary for consistency even when they are not.
   */
  @Test
  public void testNearbyVerticesSnappedWithZeroSnapRadiusEdgeSplitting() {
    S2Builder.Builder builderOptions = new S2Builder.Builder();
    builderOptions.setSplitCrossingEdges(true);
    S2Builder builder = builderOptions.build();
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(S2BuilderGraph.PolylineType.WALK);
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:180, 0:3"));
    // The second value below produces an S2Point that is distinct from 0:180 and yet is so close
    // that 0:180 is the nearest representable value when it is converted back to an S2LatLng.
    builder.addPolyline(makePolylineOrDie("90:180, 0:179.9999999999999"));
    builder.addPolyline(makePolylineOrDie("10:10, 1e-15:10"));

    assertBuildOk(builder);

    List<S2Polyline> output = layer.getPolylines();
    assertEquals(3, output.size());

    // The first two points below are not duplicates (see above).
    assertEquals("0:180, 0:180, 1e-15:10, 0:3", S2TextFormat.toString(output.get(0)));
    assertEquals("90:180, 0:180", S2TextFormat.toString(output.get(1)));
    assertEquals("10:10, 1e-15:10", S2TextFormat.toString(output.get(2)));
  }

  /**
   * A simpler version of the test above that uses addIntersection() rather than
   * splitCrossingEdges().
   */
  @Test
  public void testNearbyIntersectionSnappedWithZeroSnapRadius() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setIntersectionTolerance(S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR));
    S2Builder builder = options.build();
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:0, 0:10"));
    builder.addIntersection(makePointOrDie("1e-16:5"));

    assertBuildOk(builder);
    assertEquals("0:0, 1e-16:5, 0:10", S2TextFormat.toString(layer.getPolyline()));
  }

  /**
   * Verify that even when the splitCrossingEdges() option is used with a snap radius of zero,
   * snapped edges do not cross vertices (i.e., the input topology is preserved up the creation of
   * degeneracies). In this situation the snap radius for vertices is zero (i.e. vertices are never
   * merged) but the snap radius for edges is S2EdgeUtil.INTERSECTION_ERROR (i.e., edges are snapped
   * to vertices up to this distance away).
   *
   * <p>For this test, we create an edge AB that is 47 degrees long. We then add two vertices X,Y
   * that are 1 degree away from the endpoints of A,B and offset from AB by a distance of 0.99 *
   * S2EdgeUtil.INTERSECTION_ERROR. This will cause AB to be snapped to AXYB, where the segment XY
   * is 45 degrees long. Because all edges are geodesics, the distance from the midpoint M of XY to
   * the original edge AB is approximately 1.07 * S2EdgeUtil.INTERSECTION_ERROR (see the comments in
   * {@link S2Builder.Builder#maxEdgeDeviation()}). Let vertex C be offset from the midpoint of AB
   * by 1.03 * S2EdgeUtil.INTERSECTION_ERROR. Then C is too far away from AB to be snapped to it,
   * but it is closer to AB than M. Therefore if nothing is done about it, the input topology would
   * change.
   *
   * <p>This test checks that the algorithm detects this situation and adds another vertex Z near
   * the projection of C onto AB. This causes AB to be snapped to AXZYB so that the snapped edge
   * passes on the correct side of C.
   *
   * <p>Vertex C is represented in the form of an edge CD perpendicular to AB so that we don't need
   * to use more than one output layer. We also need to turn off the idempotent() option to ensure
   * that AB is snapped to X and Y, since otherwise the algorithm correctly determines that there is
   * no need for snapping (since there are no crossing edges in this test case).
   */
  @Test
  public void testTopologyPreservedWithZeroSnapRadiusEdgeSplitting() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSplitCrossingEdges(true);
    options.setIdempotent(false);
    S2Builder builder = options.build();
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(S2BuilderGraph.PolylineType.WALK);
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);

    double kEdgeSnapRadDegrees = S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR).degrees();
    S2Point a = S2LatLng.fromDegrees(0, -1).toPoint();
    S2Point b = S2LatLng.fromDegrees(0, 46).toPoint();
    S2Point x = S2LatLng.fromDegrees(0.99 * kEdgeSnapRadDegrees, 0).toPoint();
    S2Point y = S2LatLng.fromDegrees(0.99 * kEdgeSnapRadDegrees, 45).toPoint();
    S2Point c = S2LatLng.fromDegrees(1.03 * kEdgeSnapRadDegrees, 22.5).toPoint();
    S2Point d = S2LatLng.fromDegrees(10, 22.5).toPoint();

    builder.addEdge(a, b);
    builder.forceVertex(x);
    builder.forceVertex(y);
    builder.addEdge(c, d);

    assertBuildOk(builder);
    List<S2Polyline> output = layer.getPolylines();

    assertEquals(2, output.size());
    // The snapped edge AXZYB described above. Note that the vertices are not exactly the same as
    // the C++ unit test, because INTERSECTION_ERROR in Java is not currently equal to
    // kIntersectionError in C++, so the kEdgeSnapRadDegrees used to construct the input geometry
    // is not the same.
    assertEquals(
        "0:-1, 1.00759972308764e-13:0, 0:22.5, 1.00759972308764e-13:45, 0:46",
        S2TextFormat.toString(output.get(0)));
    // The input edge CD. Also not exactly the same as C++, as above.
    assertEquals("1.04831082301038e-13:22.5, 10:22.5", S2TextFormat.toString(output.get(1)));
    assertLessThan(
        S2EdgeUtil.robustCrossing(
            output.get(0).vertex(1), output.get(0).vertex(2),
            output.get(1).vertex(0), output.get(1).vertex(1)),
        0);
  }

  /**
   * Verify that the input topology is preserved around vertices added using forceVertex (even
   * though there are no minimum separation guarantees for such vertices).
   *
   * <p>This test is the same as the one above except for the following:
   *
   * <ul>
   *   <li>splitCrossingEdges() is false
   *   <li>we use a snap radius of S2EdgeUtil.INTERSECTION_ERROR rather than zero
   *   <li>vertex C is added using forceVertex().
   * </ul>
   */
  @Test
  public void testTopologyPreservedWithForcedVertices() {
    S2Builder.Builder options =
        new S2Builder.Builder(
            new IdentitySnapFunction(S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR)));
    options.setIdempotent(false);
    S2Builder builder = options.build();
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(S2BuilderGraph.PolylineType.WALK);
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);
    double kEdgeSnapRadDegrees = S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR).degrees();
    S2Point a = S2LatLng.fromDegrees(0, -1).toPoint();
    S2Point b = S2LatLng.fromDegrees(0, 46).toPoint();
    S2Point x = S2LatLng.fromDegrees(0.99 * kEdgeSnapRadDegrees, 0).toPoint();
    S2Point y = S2LatLng.fromDegrees(0.99 * kEdgeSnapRadDegrees, 45).toPoint();
    S2Point c = S2LatLng.fromDegrees(1.03 * kEdgeSnapRadDegrees, 22.5).toPoint();
    S2Point d = S2LatLng.fromDegrees(10, 22.5).toPoint();
    builder.addEdge(a, b);
    builder.forceVertex(x);
    builder.forceVertex(y);
    builder.forceVertex(c);
    builder.addEdge(c, d);

    assertBuildOk(builder);
    List<S2Polyline> output = layer.getPolylines();
    assertEquals(2, output.size());

    // The snapped edge AXZYB described above. Note that the vertices are not exactly the same as
    // the C++ unit test, because INTERSECTION_ERROR in Java is not currently equal to
    // kIntersectionError in C++, so the kEdgeSnapRadDegrees used to construct the input geometry
    // is not the same.
    assertEquals(
        "0:-1, 1.00759972308764e-13:0, 0:22.5, 1.00759972308764e-13:45, 0:46",
        S2TextFormat.toString(output.get(0)));
    // The input edge CD. Also not exactly the same as C++, as above.
    assertEquals("1.04831082301038e-13:22.5, 10:22.5", S2TextFormat.toString(output.get(1)));
    assertLessThan(
        S2EdgeUtil.robustCrossing(
            output.get(0).vertex(1),
            output.get(0).vertex(2),
            output.get(1).vertex(0),
            output.get(1).vertex(1)),
        0);
  }

  /** Verify that maxSnapRadius will allow snapping at S2CellId level 0. */
  @Test
  public void testkMaxSnapRadiusCanSnapAtLevel0() {
    assertLessOrEqual(
        S2CellIdSnapFunction.minSnapRadiusForLevel(0), S2BuilderSnapFunctions.maxSnapRadius());
  }

  /** Tests the S2CellIdSnapFunction at every S2Cell level. */
  @Test
  public void testS2CellIdSnappingAtAllLevels() {
    S2Polygon input = makePolygonOrDie("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0");
    for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
      S2CellIdSnapFunction snapFunction = new S2CellIdSnapFunction(level);
      S2Builder.Builder options = new S2Builder.Builder(snapFunction);
      S2Builder builder = options.build();

      S2PolygonLayer layer = new S2PolygonLayer();
      builder.startLayer(layer);
      builder.addPolygon(input);

      assertBuildOk(builder);
      S2Polygon output = layer.getPolygon();
      assertTrue(output.isValid());

      // The approxContains calls below are not guaranteed to succeed in general because
      // approxContains works by snapping both polygons together using the given tolerance and then
      // checking for containment. Since approxContains snaps to an arbitrary subset of the input
      // vertices rather than to S2CellId centers at the current level, this means that
      // corresponding vertices in "input" and "output" can snap to different sites, which causes
      // the containment test to fail. Nevertheless, by using a larger tolerance of 2 * snapRadius,
      // all calls in this test succeed (and would be likely to succeed in other similar tests).
      // (To guarantee correctness we would need to use S2CellIdSnapFunction within the
      // approxContains implementation.)
      S1Angle tolerance =
          S1Angle.min(snapFunction.snapRadius().mul(2), S2BuilderSnapFunctions.maxSnapRadius());
      assertTrue(output.approxContains(input, tolerance));
      assertTrue(input.approxContains(output, tolerance));
    }
  }

  /** Tests that snapping does not rotate the vertices of a polygon. */
  @Test
  public void testSnappingDoesNotRotateVertices() {
    // This is already tested extensively elsewhere.
    S2Polygon input =
        makePolygonOrDie(
            "49.9305505:-124.8345463, 49.9307448:-124.8299657, "
                + "49.9332101:-124.8301996, 49.9331224:-124.8341368; "
                + "49.9311087:-124.8327042, 49.9318176:-124.8312621, "
                + "49.9318866:-124.8334451");
    S2Builder.Builder options = new S2Builder.Builder(new S2CellIdSnapFunction());
    S2Builder builder = options.build();

    S2PolygonLayer layer1 = new S2PolygonLayer();
    builder.startLayer(layer1);
    builder.addPolygon(input);

    assertBuildOk(builder);
    S2Polygon output1 = layer1.getPolygon();

    // This checks that the vertices are in the same cyclic order, and that vertices have not
    // moved by more than "snapRadius".
    assertPolygonsApproxEqual(input, output1, options.snapFunction().snapRadius());

    // Check that snapping twice doesn't rotate the vertices. This also verifies that S2Builder
    // can be used again after build() is called.
    S2PolygonLayer layer2 = new S2PolygonLayer();
    builder.startLayer(layer2);
    builder.addPolygon(output1);
    assertBuildOk(builder);
    S2Polygon output2 = layer2.getPolygon();
    assertPolygonsEqual(output1, output2);
  }

  /**
   * Check that when two edges of a polyline cross, the intersection point is added to both edges.
   */
  @Test
  public void testSelfIntersectingPolyline() {
    S2Builder.Builder options = new S2Builder.Builder();
    IntLatLngSnapFunction snapFunction = new IntLatLngSnapFunction(1); // Snap to E1 coordinates
    options.setSnapFunction(snapFunction);
    options.setSplitCrossingEdges(true);
    S2Builder builder = options.build();
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    S2Polyline input = makePolylineOrDie("3:1, 1:3, 1:1, 3:3");
    S2Polyline expected = makePolylineOrDie("3:1, 2:2, 1:3, 1:1, 2:2, 3:3");
    builder.addPolyline(input);

    assertBuildOk(builder);
    assertPolylinesEqual(expected, layer.getPolyline());
  }

  /**
   * Check that when two edge of a polygon cross, the intersection point is added to both edges, and
   * that the resulting (undirected) edges can be assembled into a valid polygon.
   */
  @Test
  public void testSelfIntersectingPolygon() {
    IntLatLngSnapFunction snapFunction = new IntLatLngSnapFunction(1); // Snap to E1 coordinates
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSnapFunction(snapFunction);
    options.setSplitCrossingEdges(true);
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer(new S2PolygonLayer.Options(EdgeType.UNDIRECTED));
    builder.startLayer(layer);
    S2Polyline input = makePolylineOrDie("3:1, 1:3, 1:1, 3:3, 3:1");
    S2Polygon expected = makePolygonOrDie("1:1, 1:3, 2:2; 3:3, 3:1, 2:2");
    builder.addPolyline(input);

    assertBuildOk(builder);
    S2Polygon output = layer.getPolygon();
    assertPolygonsEqual(expected, output);
  }

  /**
   * Check that when an edge passes between two equally distant vertices, that the choice of which
   * one to snap to does not depend on the edge direction.
   */
  @Test
  public void testTieBreakingIsConsistent() {
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(2)));
    options.setIdempotent(false);
    S2Builder builder = options.build();
    builder.forceVertex(S2LatLng.fromDegrees(1, 0).toPoint());
    builder.forceVertex(S2LatLng.fromDegrees(-1, 0).toPoint());
    S2PolylineLayer layer1 = new S2PolylineLayer();
    S2PolylineLayer layer2 = new S2PolylineLayer();
    builder.startLayer(layer1);
    builder.addPolyline(makePolylineOrDie("0:-5, 0:5"));
    builder.startLayer(layer2);
    builder.addPolyline(makePolylineOrDie("0:5, 0:-5"));

    assertBuildOk(builder);
    S2Polyline output1 = layer1.getPolyline();
    S2Polyline output2 = layer2.getPolyline();
    assertEquals(3, output1.numVertices());
    assertEquals(3, output2.numVertices());
    for (int i = 0; i < 3; ++i) {
      assertEquals(output1.vertex(i), output2.vertex(2 - i));
    }
  }

  /**
   * Ensure that the Graph objects passed to S2BuilderLayer.build() methods remain valid until all
   * layers have been built.
   */
  @Test
  public void testGraphPersistence() {
    ArrayList<S2BuilderGraph> graphs = new ArrayList<>();
    ArrayList<GraphClone> clones = new ArrayList<>();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    for (int i = 0; i < 20; ++i) {
      builder.startLayer(new GraphPersistenceLayer(new GraphOptions(), graphs, clones));
      for (int n = data.uniform(10); n > 0; --n) {
        builder.addEdge(data.getRandomPoint(), data.getRandomPoint());
      }
    }

    assertBuildOk(builder);
  }

  /** Simplify a perturbed edge chain into a single edge. */
  @Test
  public void testSimplifyOneEdge() {
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(1)));
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("0:0, 1:0.5, 2:-0.5, 3:0.5, 4:-0.5, 5:0"),
        ImmutableList.of("0:0, 5:0"),
        new S2PolylineLayer.Options(),
        options);
  }

  /** Verify that nothing goes wrong when attempting to simplify a nearly antipodal edge. */
  @Test
  public void testSimplifyNearlyAntipodal() {
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(1)));
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("0:180, 0:1e-09, 32:32"),
        ImmutableList.of("0:180, 0:1e-09, 32:32"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * Construct two layers, each containing a polyline that could be simplified to a single edge on
   * its own. However, the two polylines actually cross, so make sure that the output still contains
   * the intersection vertex.
   */
  @Test
  public void testSimplifyTwoLayers() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(0.5)));
    options.setSplitCrossingEdges(true);
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("-2:-1, -1:0, 1:0, 2:1", "1:-2, 0:-1, 0:1, -1:2"),
        ImmutableList.of("-2:-1, 0:0, 2:1", "1:-2, 0:0, -1:2"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * Simplify a regular loop with 1000 vertices and a radius of 20 degrees. Turning on edge chain
   * simplification yields a dramatically smaller number of vertices than snapping alone (10
   * vertices vs 95 vertices using a snap radius of 1 degree). This is because snapping alone yields
   * vertices that stay within 1 degree of the input *vertices*, while simplifying edge chains
   * yields edges that stay within 1 degree of the input *edges*.
   */
  @Test
  public void testSimplifyOneLoop() {
    for (int i = 0; i < 2; ++i) {
      EdgeType edgeType = i == 0 ? EdgeType.DIRECTED : EdgeType.UNDIRECTED;
      S1Angle snapRadius = S1Angle.degrees(1);
      S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
      options.setSimplifyEdgeChains(true);
      S2Builder builder = options.build();
      S2PolygonLayer layer = new S2PolygonLayer(new S2PolygonLayer.Options(edgeType));
      builder.startLayer(layer);
      // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
      S2Polygon input =
          new S2Polygon(S2Loop.makeRegularLoop(new S2Point(1, 0, 0), S1Angle.degrees(20), 1000));
      builder.addPolygon(input);

      assertBuildOk(builder);

      S2Polygon output = layer.getPolygon();
      assertEquals(1, output.numLoops());
      assertGreaterOrEqual(output.loop(0).numVertices(), 10);
      assertLessOrEqual(output.loop(0).numVertices(), 12);
      assertTrue(output.boundaryNear(input, snapRadius.radians()));
    }
  }

  /**
   * We build two layers with two polylines that follow the same circular arc in opposite
   * directions, and verify that they are snapped identically. (The snap radius is adjusted so that
   * the arc is simplified into a long edge and a short edge, and therefore we would get a different
   * result if the two layers followed the edge chain in different directions.)
   */
  @Test
  public void testSimplifyOppositeDirections() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(0.5)));
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of(
            "-4:0.83, -3:0.46, -2:0.2, -1:0.05, 0:0, 1:0.5, 2:0.2, 3:0.46, 4:0.83",
            "4:.83, 3:.46, 2:.2, 1:.05, 0:0, -1:.5, -2:.2, -3:.46, -4:.83"),
        ImmutableList.of(
            "-4:0.83, -2:0.2, 4:0.83", //
            "4:0.83, -2:0.2, -4:0.83"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * We build two layers each containing a polyline, such that the polyline in the first layer could
   * be simplified to a straight line except that then it would approach the second polyline too
   * closely.
   */
  @Test
  public void testSimplifyKeepsEdgeVertexSeparation() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(1.0)));
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"),
        ImmutableList.of("0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * Test simplifying an edge chain that backtracks on itself. (This should prevent simplification,
   * since edge chains are approximated parametrically rather than geometrically.)
   */
  @Test
  public void testSimplifyBacktrackingEdgeChain() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(0.5)));
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 4:0, 3:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0"),
        ImmutableList.of("0:0, 2:0, 5:0, 2:0, 5:0, 7:0"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * As noted above, edge chains must proceed monotonically away from the source vertex in order to
   * be simplified. However in rare cases, adding a new vertex to the chain (e.g. extending ABC with
   * a vertex D) can require us to avoid a nearby vertex (E) that is closer than the previous
   * endpoint of the chain (C). A previous version of the algorithm did not handle this case
   * correctly.
   *
   * <p>The first polyline (ABC) below cannot be simplified to AC, as that edge would pass too close
   * to the first vertex of the second polyline DE. Vertex D is not avoided while processing the
   * edge AB because AD > AB; instead it should be processed when edge BC is added to the simplified
   * chain.
   */
  @Test
  public void testSimplifyAvoidsBacktrackingVertices() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(1.0)));
    options.setSimplifyEdgeChains(true);
    assertLessThan(
        S2EdgeUtil.getDistance(
            makePointOrDie("0:1.05"), makePointOrDie("0:0"), makePointOrDie("1:2")),
        options.snapFunction().minEdgeVertexSeparation());
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("0:0, 1:0.1, 1:2", "0:1.05, -10:1.05"),
        ImmutableList.of("0:0, 1:0.1, 1:2", "0:1.05, -10:1.05"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * Make sure that simplification does not create long edges such that the midpoint of the edge
   * might be further than maxEdgeDeviation() from an input edge. In the example below, vertices are
   * snapped to integer lat/lng coordinates, and the snap radius is approximately 0.707 degrees.
   * Snapping moves the input vertices perpendicular to the input edge by just slightly less than
   * the snap radius (0.693 degrees). Now the midpoint of the snapped edge is about 0.98 degrees
   * from the input edge, which causes an extra site to be added at the midpoint of the original
   * edge.
   *
   * <p>When simplifyEdgeChains() is enabled, then usually an extra site like this would be
   * simplified away (because the simplified edge would still be within snapRadius() of all the
   * input vertices) except that there is an explicit check in S2Builder that prevents this. (If the
   * check is removed then this test fails.)
   */
  @Test
  public void testSimplifyLimitsEdgeDeviation() {
    // E0 coordinates
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0));
    options.setSimplifyEdgeChains(true);
    checkPolylineLayerBothEdgeTypes(
        ImmutableList.of("-30.49:-29.51, 29.51:30.49"),
        ImmutableList.of("-30:-30, -1:1, 30:30"),
        new S2PolylineLayer.Options(),
        options);
  }

  /**
   * Create several nested concentric loops, and verify that the loops are still nested after
   * simplification.
   */
  @Test
  public void testSimplifyPreservesTopology() {
    int kNumLoops = 20;
    int kNumVerticesPerLoop = 1000;
    S1Angle kBaseRadius = S1Angle.degrees(5);
    S1Angle kSnapRadius = S1Angle.degrees(0.1);
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(kSnapRadius));
    options.setSimplifyEdgeChains(true);
    S2Builder builder = options.build();
    ArrayList<S2Polygon> inputs = new ArrayList<>();
    ArrayList<S2PolygonLayer> outputLayers = new ArrayList<>();
    for (int j = 0; j < kNumLoops; ++j) {
      // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
      S1Angle radius =
          S1Angle.radians(kBaseRadius.radians() + 0.7 * j * j / kNumLoops * kSnapRadius.radians());
      S2Polygon polygon =
          new S2Polygon(S2Loop.makeRegularLoop(new S2Point(1, 0, 0), radius, kNumVerticesPerLoop));
      inputs.add(polygon);

      S2PolygonLayer polygonLayer = new S2PolygonLayer();
      outputLayers.add(polygonLayer);
      builder.startLayer(polygonLayer);
      builder.addPolygon(polygon);
    }

    assertBuildOk(builder);
    S2Polygon previousOutputPolygon = null;
    for (int j = 0; j < kNumLoops; ++j) {
      S2Polygon outputPolygon = outputLayers.get(j).getPolygon();
      assertTrue(outputPolygon.boundaryNear(inputs.get(j), kSnapRadius.radians()));
      if (j > 0) {
        assertTrue(outputPolygon.contains(previousOutputPolygon));
      }
      previousOutputPolygon = outputPolygon;
    }
  }

  @Test
  public void testSimplifyRemovesSiblingPairs() {
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0)); // E0 coords
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setSiblingPairs(GraphOptions.SiblingPairs.DISCARD);

    // Check that there is no sibling pair without simplification.
    checkPolylineVectorLayer(
        ImmutableList.of("0:0, 0:10", "0:10, 0.6:5, 0:0"),
        ImmutableList.of("0:0, 0:10, 1:5, 0:0"),
        layerOptions,
        options);

    // Now check that (1) simplification produces a sibling pair, and (2) the sibling pair is
    // removed (since we requested it).
    options.setSimplifyEdgeChains(true);
    checkPolylineVectorLayer(
        ImmutableList.of("0:0, 0:10", "0:10, 0.6:5, 0:0"),
        ImmutableList.of(),
        layerOptions,
        options);
  }

  @Test
  public void testSimplifyMergesDuplicateEdges() {
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0)); // E0 coords
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setDuplicateEdges(GraphOptions.DuplicateEdges.MERGE);

    // Check that there are no duplicate edges without simplification.
    checkPolylineVectorLayer(
        ImmutableList.of("0:0, 0:10", "0:0, 0.6:5, 0:10"),
        ImmutableList.of("0:0, 0:10", "0:0, 1:5, 0:10"),
        layerOptions,
        options);

    // Now check that (1) simplification produces a duplicate edge pair, and (2) the duplicate
    // pair is merged (since we requested it).
    options.setSimplifyEdgeChains(true);
    checkPolylineVectorLayer(
        ImmutableList.of("0:0, 0:10", "0:0, 0.6:5, 0:10"),
        ImmutableList.of("0:0, 0:10"),
        layerOptions,
        options);
  }

  @Test
  public void testSimplifyKeepsForcedVertices() {
    S2Builder.Builder options =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.radians(1e-15)));
    options.setSimplifyEdgeChains(true);
    S2Builder builder = options.build();
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:0, 0:1, 0:2, 0:3"));
    builder.forceVertex(makePointOrDie("0:1"));

    assertBuildOk(builder);
    assertEquals("0:0, 0:1, 0:3", S2TextFormat.toString(layer.getPolyline()));
  }

  /** Check that input edge ids are assigned in order. */
  @Test
  public void testInputEdgeIdAssignment() {
    checkInputEdgeIds(
        // A polyline with two edges.
        ImmutableList.of("0:0, 0:1, 0:2"),
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:1", IntVector.of(0)),  // First polyline edge has id 0.
                Pair.of("0:1, 0:2", IntVector.of(1)))), // Second polyline edge has id 1.
        new GraphOptions(),
        new S2Builder.Builder());
  }

  /** Check that the siblings of undirected edges do not have InputEdgeIds. */
  @Test
  public void testUndirectedSiblingsDontHaveInputEdgeIds() {
    GraphOptions graphOptions = new GraphOptions();
    graphOptions.setEdgeType(EdgeType.UNDIRECTED);
    checkInputEdgeIds(
        ImmutableList.of("0:0, 0:1, 0:2"),
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:1", IntVector.of(0)),
                Pair.of("0:1, 0:2", IntVector.of(1)),
                Pair.of("0:1, 0:0", IntVector.empty()),
                Pair.of("0:2, 0:1", IntVector.empty()))),
        graphOptions,
        new S2Builder.Builder());
  }

  /** Check that edges created by SiblingPairs.CREATE do not have InputEdgeIds. */
  @Test
  public void testCreatedSiblingsDontHaveInputEdgeIds() {
    GraphOptions graphOptions = new GraphOptions();
    graphOptions.setSiblingPairs(GraphOptions.SiblingPairs.CREATE);
    checkInputEdgeIds(
        ImmutableList.of("0:0, 0:1, 0:2"),
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:1", IntVector.of(0)),
                Pair.of("0:1, 0:2", IntVector.of(1)),
                Pair.of("0:1, 0:0", IntVector.empty()),
                Pair.of("0:2, 0:1", IntVector.empty()))),
        graphOptions,
        new S2Builder.Builder());
  }

  /** Tests that input edge ids are merged when directed edges are merged. */
  @Test
  public void testEdgeMergingDirected() {
    GraphOptions graphOptions = new GraphOptions();
    graphOptions.setDuplicateEdges(GraphOptions.DuplicateEdges.MERGE);
    checkInputEdgeIds(
        ImmutableList.of("0:0, 0:1", "0:0, 0:1"),
        new EdgeInputEdgeIds(ImmutableList.of(Pair.of("0:0, 0:1", IntVector.of(0, 1)))),
        graphOptions,
        new S2Builder.Builder());
  }

  /** Tests that input edge ids are merged when undirected edges are merged. */
  @Test
  public void testEdgeMergingUndirected() {
    GraphOptions graphOptions = new GraphOptions();
    graphOptions.setDuplicateEdges(GraphOptions.DuplicateEdges.MERGE);
    graphOptions.setSiblingPairs(GraphOptions.SiblingPairs.KEEP);
    checkInputEdgeIds(
        ImmutableList.of("0:0, 0:1, 0:2", "0:0, 0:1", "0:2, 0:1"),
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:1", IntVector.of(0, 2)),
                Pair.of("0:1, 0:2", IntVector.of(1)),
                Pair.of("0:2, 0:1", IntVector.of(3)))),
        graphOptions,
        new S2Builder.Builder());
  }

  /**
   * Check that when an input edge is snapped to a chain that includes degenerate edges, and the
   * edge chain is simplified, that the InputEdgeIds attached to those degenerate edges are
   * transferred to the simplified edge. For example (using integers for vertices), an edge chain
   * 1.2, 2.2, 2.3 that is simplified to 1.3 should get the InputEdgeIds associated with all three
   * original edges. (This ensures that the labels attached to those edges are also transferred.)
   *
   * <p>This also tests that degenerate edges at the start and end of the simplified chain are *not*
   * merged. (It's up to the output layer to decide what to do with these edges. The only reason we
   * merge degenerate edges in the interior of the interior of the simplified edge is because those
   * edges are being removed from the graph.)
   */
  @Test
  public void testSimplifyDegenerateEdgeMergingEasy() {
    GraphOptions graphOptions = new GraphOptions();
    graphOptions.setDegenerateEdges(GraphOptions.DegenerateEdges.KEEP);
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0));
    options.setSimplifyEdgeChains(true);
    checkInputEdgeIds(
        ImmutableList.of("0:0, 0:0.1, 0:1.1, 0:1, 0:0.9, 0:2, 0:2.1"),
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:0", IntVector.of(0)),
                Pair.of("0:0, 0:2", IntVector.of(1, 2, 3, 4)),
                Pair.of("0:2, 0:2", IntVector.of(5)))),
        graphOptions,
        options);
  }

  /**
   * This is a harder version of the test above. Now there are several edge chains that overlap each
   * other in both directions, and several degenerate edges at that middle vertex. This tests that
   * if exactly one edge chain contains a degenerate edge in input edge order (e.g., the input order
   * was AB, BB, BC), then the degenerate edge is assigned to that edge chain. Otherwise the edge is
   * assigned to an arbitrary chain.
   */
  @Test
  public void testSimplifyDegenerateEdgeMergingHard() {
    GraphOptions graphOptions = new GraphOptions(); // Default options keep everything.
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0));
    options.setSimplifyEdgeChains(true);
    ImmutableList<String> input =
        ImmutableList.of(
            "0:1, 0:1.1",
            "0:0, 0:1, 0:2", // Degenerate edge defined before chain
            "0:0, 0:0.9, 0:1, 0:1.1, 0:2", // Degenerate edge defined in chain
            "0:2, 0:1, 0:0.9, 0:0", // Defined in chain, chain reversed
            "0:2, 0:1, 0:0",
            "0:1.1, 0:1",
            "0:1, 0:1.1"); // Defined after chain

    EdgeInputEdgeIds expected =
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:2", IntVector.of(0, 1, 2)),
                Pair.of("0:0, 0:2", IntVector.of(3, 4, 5, 6)),
                Pair.of("0:2, 0:0", IntVector.of(7, 8, 9)),
                Pair.of("0:2, 0:0", IntVector.of(10, 11, 12, 13))));

    checkInputEdgeIds(input, expected, graphOptions, options);

    // Now try the same test with undirected edges. This results in four more simplified edges
    // that are not labelled with any input edge ids.
    ImmutableList.Builder<Pair<String, IntVector>> b = ImmutableList.builder();
    b.addAll(expected);
    b.addAll(
        ImmutableList.of(
            Pair.of("0:0, 0:2", IntVector.empty()),
            Pair.of("0:0, 0:2", IntVector.empty()),
            Pair.of("0:2, 0:0", IntVector.empty()),
            Pair.of("0:2, 0:0", IntVector.empty())));

    graphOptions.setEdgeType(EdgeType.UNDIRECTED);
    checkInputEdgeIds(input, new EdgeInputEdgeIds(b.build()), graphOptions, options);
  }

  /**
   * Check that degenerate edges are assigned to an edge in the correct layer when multiple edge
   * chains in different layers are simplified in the same way (i.e., yielding a set of identical or
   * reversed edges in different layers).
   */
  @Test
  public void testSimplifyDegenerateEdgeMergingMultipleLayers() {
    GraphOptions graphOptions = new GraphOptions(); // Default options keep everything.
    S2Builder.Builder options = new S2Builder.Builder(new IntLatLngSnapFunction(0));
    options.setSimplifyEdgeChains(true);

    // Note below that the edge chains in different layers have different vertex locations,
    // different number of interior vertices, different degenerate edges, etc, and yet they can all
    // be simplified together.
    ArrayList<List<String>> input = new ArrayList<>();
    ArrayList<EdgeInputEdgeIds> expected = new ArrayList<>();

    input.add(
        ImmutableList.of(
            // Edge 0 is degenerate after snapping to 0:5,0:5.
            "0.1:5, 0:5.2",
            // Edge 1 will snap to 0:0,0:10.
            "0.1:0, 0:9.9",
            // Edge 2 will snap to 0:10,0:0.
            "0:10.1, 0:0.1",
            // Degenerate Edge 3 0:3 is defined last in layer.
            "0:3.1, 0:2.9"
            ));
    expected.add(
        new EdgeInputEdgeIds(
            ImmutableList.of(
                // Two output edges, each with one of the degenerate input edges assigned to it.
                Pair.of("0:0, 0:10", IntVector.of(0, 1)),
                Pair.of("0:10, 0:0", IntVector.of(2, 3)))));

    input.add(
        ImmutableList.of(
            // Edge 4 is degenerate after snapping to 0:3.
            "0.1:3, 0:3.2",
            // Chain of Edges 5 and 6 snap and simplify to 0:0, 0:10
            "-0.1:0, 0:4.1, 0:9.9",
            // Chain of three edges 7,8,9 snaps and simplifies to 0:10, 0:0. Edge 8 is degenerate,
            // as both endpoints snap to 0:7.
            "0.1:9.9, 0:7, 0.1:6.9, 0.1:0.2"
            ));
    expected.add(
        new EdgeInputEdgeIds(
            ImmutableList.of(
                // Two output edges, each simplified from two non-degenerate and one degenerate
                // input edges.
                Pair.of("0:0, 0:10", IntVector.of(4, 5, 6)),
                Pair.of("0:10, 0:0", IntVector.of(7, 8, 9)))));

    input.add(
        ImmutableList.of(
            // Chain of three edges 10,11,12 snaps and simplifies to 0:0,0:10 with Edge 11
            // becoming degenerate after snapping both endpoints to 0:6.
            "0.2:0.3, 0.1:6, 0:5.9, 0.1:10.2",
            // Edge 13 snaps to 0:0,0:10, duplicating the previous chain.
            "0.1:0.1, 0:9.8",
            // Edge 14 is degenerate after snapping both endpoints to 0:2. It is assigned to the
            // previous chain, rather than the first one.
            "0.1:2, 0:2.1"
            ));
    expected.add(
        new EdgeInputEdgeIds(
            ImmutableList.of(
                Pair.of("0:0, 0:10", IntVector.of(10, 11, 12)),
                Pair.of("0:0, 0:10", IntVector.of(13, 14)))));

    S2Builder builder = options.build();
    for (int i = 0; i < expected.size(); ++i) {
      S2BuilderLayer layer = new InputEdgeIdCheckingLayer(expected.get(i), graphOptions);
      builder.startLayer(layer);
      for (String inputStr : input.get(i)) {
        S2Polyline polyline = makePolylineOrDie(inputStr);
        builder.addPolyline(polyline);
      }
      System.out.println();
    }

    assertBuildOk(builder);
  }

  /**
   * To produce correct output in this example, the algorithm needs fall back to high precision
   * predicates when the output of the normal predicates is uncertain.
   */
  @SuppressWarnings("FloatingPointLiteralPrecision") // to ensure exactly the same as the C++ test.
  @Test
  public void testHighPrecisionPredicates() {
    ImmutableList<S2Point> vertices =
        ImmutableList.of(
            new S2Point(-0.1053119128423491, -0.80522217121852213, 0.58354661852470235),
            new S2Point(-0.10531192039134209, -0.80522217309706012, 0.58354661457019508),
            new S2Point(-0.10531192039116592, -0.80522217309701472, 0.58354661457028933));
    S2Polyline input = new S2Polyline(vertices);
    S1Angle snapRadius = S2EdgeUtil.INTERSECTION_MERGE_RADIUS;
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
    options.setIdempotent(false);
    S2Builder builder = options.build();
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    builder.forceVertex(
        new S2Point(-0.10531192039134191, -0.80522217309705857, 0.58354661457019719));
    builder.addPolyline(input);

    assertBuildOk(builder);
  }

  /**
   * This test constructs many small, random inputs such that the output is likely to be
   * inconsistent unless high-precision predicates are used.
   */
  @Test
  public void testHighPrecisionStressTest() {
    S1Angle snapRadius = S2EdgeUtil.INTERSECTION_MERGE_RADIUS;
    // Some S2Builder calculations use an upper bound that takes into account S1ChordAngle errors.
    // We sometimes try perturbing points by very close to that distance in an attempt to expose
    // errors.
    S1ChordAngle ca = S1ChordAngle.fromS1Angle(snapRadius);
    S1Angle snapRadiusWithError =
        ca.plusError(
                ca.getS1AngleConstructorMaxError() + S2EdgeUtil.getUpdateMinDistanceMaxError(ca))
            .toAngle();

    int nonDegenerate = 0;
    int kIters = 8000 * ITERATION_MULTIPLIER;
    for (int iter = 0; iter < kIters; ++iter) {
      data.setSeed(iter + 1); // Easier to reproduce a specific case.

      // We construct a nearly degenerate triangle where one of the edges is sometimes very short.
      // Then we add a forced vertex somewhere near the shortest edge. Then after snapping, we
      // check that (1) the edges still form a loop, and (2) if the loop is non-degenerate, then
      // it has the same orientation as the original triangle.
      //
      // v1 is located randomly. (v0,v1) is the longest of the three edges.
      // v2 is located along (v0,v1) but is perturbed by up to 2 * snapRadius.
      S2Point v1 = choosePoint();
      S2Point v0Dir = choosePoint();
      double d0 = pow(1e-16, data.nextDouble());
      S2Point v0 = S2EdgeUtil.getPointOnLine(v1, v0Dir, S1Angle.radians(d0));
      double d2 = 0.5 * d0 * pow(1e-16, pow(data.nextDouble(), 2));
      S2Point v2 =
          data.samplePoint(
              S2Cap.fromAxisAngle(
                  S2EdgeUtil.getPointOnLine(v1, v0Dir, S1Angle.radians(d2)), snapRadius.mul(2)));
      // Vary the edge directions by randomly swapping v0 and v2.
      if (data.oneIn(2)) {
        S2Point tmp = v0;
        v0 = v2;
        v2 = tmp;
      }

      // The forced vertex (v3) is either located near the (v1, v2) edge. We perturb it either in
      // a random direction from v1 or v2, or perpendicular to (v1, v2) starting from an interior
      // edge point.
      S1Angle d3 = data.oneIn(2) ? snapRadius : snapRadiusWithError;
      if (data.oneIn(3)) {
        d3 = d3.mul(1.5 * data.nextDouble());
      }
      S2Point v3;
      if (data.oneIn(5)) {
        v3 = data.oneIn(2) ? v1 : v2;
        v3 = S2EdgeUtil.getPointOnLine(v3, choosePoint(), d3);
      } else {
        v3 = S2EdgeUtil.interpolate(pow(1e-16, data.nextDouble()), v1, v2);
        v3 = S2EdgeUtil.getPointOnLine(v3, v1.crossProd(v2).normalize(), d3);
      }
      S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
      options.setIdempotent(false);
      S2Builder builder = options.build();

      S2PolygonLayer layer = new S2PolygonLayer();
      builder.startLayer(layer);
      builder.forceVertex(v3);
      builder.addEdge(v0, v1);
      builder.addEdge(v1, v2);
      builder.addEdge(v2, v0);
      S2Error error = new S2Error();
      if (!builder.build(error)) {
        System.err.println("d0=" + d0 + ", d2=" + d2 + ", d3=" + d3);
      }

      if (error.ok() && !layer.getPolygon().isEmpty()) {
        S2Polygon output = layer.getPolygon();
        assertEquals(1, output.numLoops());
        if (output.numLoops() == 1) {
          assertTrue(output.isValid());
          assertEquals(
              "d0=" + d0 + ", d2=" + d2 + ", d3=" + d3,
              S2Predicates.sign(v0, v1, v2) > 0,
              output.loop(0).isNormalized());
          ++nonDegenerate;
        }
      }
    }
    System.out.println(nonDegenerate + " non-degenerate out of " + kIters);
    assertGreaterOrEqual(nonDegenerate, kIters / 10);
  }

  @Test
  public void testSelfIntersectionStressTest() {
    // How many iterations to run in total.
    int kIters = 50 * ITERATION_MULTIPLIER;

    // How many iterations to print details for.
    int kPrintIters = 50;

    // Note that the number of intersections (and the running time) is quadratic in the number of
    // vertices. With 200 input vertices, the output consists of about 2300 loops and 9000
    // vertices.
    // TODO(torrey): Increase this to at least 50.
    int numVertices = 10;

    for (int iter = 0; iter < kIters; ++iter) {
      if (iter < kPrintIters) {
        System.out.println("Iteration: " + iter);
      }
      data.setSeed(iter + 1); // Easier to reproduce a specific case.
      Stopwatch timer = Stopwatch.createStarted();

      // Select a cap that points will be sampled from. The minimum radius is about 36cm on the
      // Earth's surface. The performance is reduced for radii much smaller than this because
      // S2ShapeIndex only indexes regions down to about 1 cm across.
      S2Cap cap = data.getRandomCap(1e-14, 1e-2);
      S1Angle radius = cap.angle();

      S2Builder.Builder options = new S2Builder.Builder();
      options.setSplitCrossingEdges(true);

      // For half the tests, set an IntLatLng snap function. This may cause most of the random
      // points sampled within the cap to snap to a small number of sites, resulting in many
      // degenerate edges.
      if (data.oneIn(2)) {
        // The exponent such that vertices will not move by more than the cap radius.
        int minExp = IntLatLngSnapFunction.exponentForMaxSnapRadius(radius);
        // Select an exponent between the minimum and up to five levels higher, or the max.
        int exponent = min(IntLatLngSnapFunction.MAX_EXPONENT, minExp + data.uniform(5));
        S2Builder.SnapFunction f = new IntLatLngSnapFunction(exponent);
        options.setSnapFunction(f);
        if (iter < kPrintIters) {
          System.out.println(
              "  Snapping to IntLatLng exponent " + exponent
              + " with min vertex separation " + f.minVertexSeparation()
              + " within a cap of radius " + radius);
        }
      } else {
        if (iter < kPrintIters) {
          System.out.println(
              "  Snapping with an identity snap function, splitting crossing edges,"
              + " with min vertex separation " + options.snapFunction().minVertexSeparation()
              + " within a cap of radius " + radius);
        }
      }

      S2Builder builder = options.build();

      S2PolygonLayer layer = new S2PolygonLayer(new S2PolygonLayer.Options(EdgeType.UNDIRECTED));
      builder.startLayer(layer);
      S2Point[] vertices = new S2Point[numVertices];
      for (int i = 0; i < vertices.length; ++i) {
        vertices[i] = data.samplePoint(cap);
      }
      vertices[vertices.length - 1] = vertices[0];
      S2Polyline input = new S2Polyline(vertices);
      builder.addPolyline(input);
      if (iter < kPrintIters) {
        // System.out.println("  Input S2Polyline: " + S2TextFormat.toString(input));
      }

      assertBuildOk(builder);

      S2Polygon output = layer.getPolygon();
      S2Error error = new S2Error();
      boolean validationError = output.findValidationError(error);
      assertFalse(error.toString(), validationError);

      // Print timing and size details for the first 50 iterations.
      if (iter < kPrintIters) {
        // System.out.println("    Output S2Polygon: " + S2TextFormat.toString(output));
        Platform.printf(
            System.out,
            "    Timing for iter=%4d: microseconds=%4d, radius=%8.3g, loops=%d, vertices=%d\n",
            iter,
            timer.elapsed(MICROSECONDS),
            cap.angle().radians(),
            output.numLoops(),
            output.getNumVertices());
      }
    }
  }

  @Test
  public void testFractalStressTest() {
    // TODO(torrey): This should be 1000 * ITERATION_MULTIPLIER, but it gets very slow
    int kIters = 100 * ITERATION_MULTIPLIER;
    for (int iter = 0; iter < kIters; ++iter) {
      data.setSeed(iter + 1); // Easier to reproduce a specific case.
      S2FractalBuilder fractalBuilder = new S2FractalBuilder(data.rand);
      fractalBuilder.setLevelForApproxMaxEdges(12800);
      fractalBuilder.setLevelForApproxMinEdges(12);
      fractalBuilder.setFractalDimension(1.5 + 0.5 * data.nextDouble());
      S2Polygon input =
          new S2Polygon(fractalBuilder.makeLoop(data.getRandomFrame(), S1Angle.degrees(20)));
      S2Builder.Builder options = new S2Builder.Builder();
      if (data.oneIn(3)) {
        int exponent = data.uniform(11);
        options.setSnapFunction(new IntLatLngSnapFunction(exponent));
      } else if (data.oneIn(2)) {
        int level = data.uniform(20);
        options.setSnapFunction(new S2CellIdSnapFunction(level));
      } else {
        options.setSnapFunction(
            new IdentitySnapFunction(S1Angle.degrees(10 * pow(1e-4, data.nextDouble()))));
      }
      S2Builder builder = options.build();
      S2PolygonLayer layer = new S2PolygonLayer();
      builder.startLayer(layer);
      builder.addPolygon(input);

      assertBuildOk(builder);
      S2Polygon output = layer.getPolygon();
      S2Error error = new S2Error();
      assertFalse(error.toString(), output.findValidationError(error));

      // For debugging problems with a specific iteration.
      if (iter == -1) {
        System.out.println("S2Polygon: " + S2TextFormat.toString(input));
        System.out.println("S2Polygon: " + S2TextFormat.toString(output));
      }

      // Print details for the first 50 iters.
      if (iter < 50) {
        Platform.printf(
            System.out,
            "iter=%4d: input vertices=%d, output vertices=%d\n",
            iter, input.getNumVertices(), output.getNumVertices());
      }
    }
  }

  /**
   * The test for whether one Voronoi site excludes another along a given input edge boils down to a
   * test of whether two angle intervals "a" and "b" overlap. Let "ra" and "rb" be the semi-widths
   * of the two intervals, and let "d" be the angle between their centers. Then "a" contains "b" if
   * (rb + d <= ra), and "b" contains "a" if (rb - d >= ra). However the actual code uses the sines
   * of the angles, e.g. sin(rb + d) <= sin(ra). This works fine most of the time, but the first
   * condition (rb + d <= ra) also needs to check that rb + d < 90 degrees. This test verifies that
   * case.
   */
  @Test
  public void testAdjacentCoverageIntervalsSpanMoreThan90Degrees() {
    // The following 3 tests have d < 90, d = 90, and d > 90 degrees, but in all 3 cases,
    // rb + d > 90 degrees.
    checkSnappingWithForcedVertices("0:0, 0:80", S1Angle.degrees(60), "0:0, 0:70", "0:0, 0:70");
    checkSnappingWithForcedVertices("0:0, 0:80", S1Angle.degrees(60), "0:0, 0:90", "0:0, 0:90");
    checkSnappingWithForcedVertices("0:0, 0:80", S1Angle.degrees(60), "0:0, 0:110", "0:0, 0:110");

    // This test has d = 180 degrees, i.e. the two sites project to points that are 180 degrees
    // apart along the input edge. The snapped edge doesn't stay within maxEdgeDeviation() of the
    // input edge, so an extra site is added and it is snapped again (yielding two edges). The case
    // we are testing here is the first call to SnapEdge() before adding the site.
    checkSnappingWithForcedVertices(
        "0:10, 0:170", S1Angle.degrees(50), "47:0, 49:180", "47:0, 0:90, 49:180");

    // This test has d = 220 degrees, i.e. when the input edge is snapped it goes the "wrong way"
    // around the sphere. Again, the snapped edge is too far from the input edge so an extra site is
    // added and it is resnapped.
    checkSnappingWithForcedVertices(
        "0:10, 0:170", S1Angle.degrees(70), "0:-20, 0:-160", "0:-20, 0:90, 0:-160");

    // Without using forced vertices, the maximum angle between the coverage interval centers is
    // d = 300 degrees. This would use an edge 180 degrees long, and then place two sites 60 degrees
    // past either endpoint. With forced vertices we can increase the snap radius to 70 degrees and
    // get an angle of up to d = 320 degrees, but the sites are only 40 degrees apart (which is why
    // it requires forced vertices). The test below is an approximation of this situation with
    // d = 319.6 degrees.
    checkSnappingWithForcedVertices(
        "0:0.1, 0:179.9", S1Angle.degrees(70), "0:-69.8, 0:-110.2", "0:-69.8, 0:90, 0:-110.2");
  }

  /** Test behavior of S2Builder when some vertices are NaN. */
  @Test
  public void testNaNVertices() {
      ImmutableList<ImmutableList<S2Point>> loops = ImmutableList.of(
          ImmutableList.of(
              new S2Point(NaN, NaN, NaN),
              new S2Point(NaN, NaN, NaN),
              new S2Point(NaN, NaN, NaN))
      );

      S2Builder builder = new S2Builder.Builder(
          new IdentitySnapFunction(S1Angle.radians(1e-15))).build();
      S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
      builder.startLayer(layer);
      builder.addShape(S2LaxPolygonShape.create(loops));

      // The point of this test is to ensure that NaN values in inputs don't cause a crash (in C++)
      // or thrown exception in Java. We aren't too concerned about the actual error code, but in
      // both C++ and Java, the operation fails with an S2Error of BUILDER_SNAP_RADIUS_TOO_SMALL.
      // In C++ this is because the distance between two NaN vertices happens to be computed as Pi.
      // In Java, an IllegalArgumentException is thrown from the S1ChordAngle constructor when
      // attempting to compute the distance between NaN vertices, but the exception is caught by
      // S2Builder.build(), which sets an S2Error using the same error code for consistency.
      S2Error error = new S2Error();
      boolean result = builder.build(error);
      assertFalse(result);
      assertEquals(S2Error.Code.BUILDER_SNAP_RADIUS_TOO_SMALL, error.code());
    }

  /**
   * Test with a polygon that caused the obsolete S2PolygonBuilder class to generate an invalid
   * output polygon (duplicate edges).
   */
  @Test
  public void testOldS2PolygonBuilderBug() {
    S2Polygon input =
        makePolygonOrDie(
            "32.2983095:72.3416582, 32.2986281:72.3423059, "
                + "32.2985238:72.3423743, 32.2987176:72.3427807, "
                + "32.2988174:72.3427056, 32.2991269:72.3433480, "
                + "32.2991881:72.3433077, 32.2990668:72.3430462, "
                + "32.2991745:72.3429778, 32.2995078:72.3436725, "
                + "32.2996075:72.3436269, 32.2985465:72.3413832, "
                + "32.2984558:72.3414530, 32.2988015:72.3421839, "
                + "32.2991552:72.3429416, 32.2990498:72.3430073, "
                + "32.2983764:72.3416059");
    assertTrue(input.isValid());

    S1Angle snapRadius = metersToAngle(20 / 0.866);
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    builder.addPolygon(input);

    assertBuildOk(builder);
    S2Polygon output = layer.getPolygon();
    assertTrue(output.isValid());
    S2Polygon expected =
        makePolygonOrDie(
            "32.2991552:72.3429416, 32.2991881:72.3433077, 32.2996075:72.3436269; "
                + "32.2988015:72.3421839, 32.2985465:72.3413832, 32.2983764:72.3416059, "
                + "32.2985238:72.3423743, 32.2987176:72.3427807");
    assertPolygonsEqual(expected, output);
  }

  /**
   * This reproduces a bug (in the C++ implementation) that used to attempt to create a separation
   * site where none was necessary.
   */
  @SuppressWarnings("FloatingPointLiteralPrecision") // so as to be identical to the C++ test.
  @Test
  public void testIncorrectSeparationSiteBug() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setIdempotent(false);
    options.setSplitCrossingEdges(true);
    S2Builder builder = options.build();
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    builder.addEdge(
        new S2Point(-0.500944389640767, -0.8654794731750945, 0),
        new S2Point(1, 1.7786363250284876e-322, 4.772992939485661E-65));
    builder.forceVertex(new S2Point(1, 0, -4.772992939485661e-65));
    builder.forceVertex(new S2Point(1, 2.2604e-320, 4.772992939485662e-65));

    assertBuildOk(builder);
  }

  @Test
  public void testPushPopLabel() {
    // TODO(user): Test more thoroughly.
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    builder.pushLabel(1);
    builder.popLabel();
  }

  /**
   * This layer makes both a shallow and a deep copy of the S2BuilderGraph object passed to its
   * build() method and appends them to two vectors. Furthermore, it verifies that the shallow and
   * deep copies of any graphs previously appended to those vectors are still equal.
   */
  static class GraphPersistenceLayer implements S2BuilderLayer {
    private final GraphOptions graphOptions;
    // Shallow copies.
    private final List<S2BuilderGraph> graphs;
    // Deep copies.
    private final List<GraphClone> clones;

    /** Constructs a GraphPersistenceLayer using the provided lists for storage. */
    public GraphPersistenceLayer(
        GraphOptions graphOptions, List<S2BuilderGraph> graphs, List<GraphClone> clones) {
      this.graphOptions = graphOptions;
      this.graphs = graphs;
      this.clones = clones;
    }

    @Override
    public S2Builder.GraphOptions graphOptions() {
      return graphOptions;
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      // Verify that all graphs built so far are unchanged.
      for (int i = 0; i < graphs.size(); ++i) {
        assertGraphsEqual(clones.get(i).graph(), graphs.get(i));
      }
      graphs.add(g);
      clones.add(new GraphClone(g));
      return true;
    }

    /** Verifies that two graphs have the same vertices and edges. */
    private void assertGraphsEqual(S2BuilderGraph expected, S2BuilderGraph actual) {
      // List.equals() is true if same-indexed list elements are equals().
      assertEquals(expected.vertices(), actual.vertices());
      // Check EdgeLists for equality with EdgeList.isEqualTo(), inherited from IntPairVector.
      // Verifies that the source and destination vertex ids are the same, in the same order.
      assertTrue(expected.edges().isEqualTo(actual.edges()));
      // Verifies that the IdSetId for each input edge id is the same.
      assertTrue(expected.inputEdgeIdSetIds().equals(actual.inputEdgeIdSetIds()));
    }
  }

  private void checkPolylineLayer(
      List<String> inputStrs,
      List<String> expectedStrs,
      S2PolylineLayer.Options layerOptions,
      S2Builder.Builder builderOptions) {
    if (builderOptions == null) {
      builderOptions = new S2Builder.Builder();
    }

    S2Builder builder = builderOptions.build();
    ArrayList<S2PolylineLayer> outputLayers = new ArrayList<>();
    for (String inputStr : inputStrs) {
      S2PolylineLayer polylineLayer = new S2PolylineLayer(layerOptions);
      outputLayers.add(polylineLayer);
      builder.startLayer(polylineLayer);
      builder.addPolyline(makePolylineOrDie(inputStr));
    }
    assertBuildOk(builder);
    ArrayList<String> outputStrs = new ArrayList<>();
    for (S2PolylineLayer polylineLayer : outputLayers) {
      outputStrs.add(S2TextFormat.toString(polylineLayer.getPolyline()));
    }
    assertEquals(expectedStrs, outputStrs);
  }

  private void checkPolylineVectorLayer(
      List<String> inputStrs,
      List<String> expectedStrs,
      S2PolylineVectorLayer.Options layerOptions,
      S2Builder.Builder builderOptions) {
    if (builderOptions == null) {
      builderOptions = new S2Builder.Builder();
    }

    S2Builder builder = builderOptions.build();
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);
    for (String inputStr : inputStrs) {
      builder.addPolyline(makePolylineOrDie(inputStr));
    }
    assertBuildOk(builder);
    ArrayList<String> outputStrs = new ArrayList<>();
    for (S2Polyline polyline : layer.getPolylines()) {
      outputStrs.add(S2TextFormat.toString(polyline));
    }
    assertEquals(expectedStrs, outputStrs);
  }

  private void checkPolylineLayerBothEdgeTypes(
      List<String> inputStrs,
      List<String> expectedStrs,
      S2PolylineLayer.Options layerOptions,
      S2Builder.Builder builderOptions) {
    if (builderOptions == null) {
      builderOptions = new S2Builder.Builder();
    }
    layerOptions.setEdgeType(EdgeType.DIRECTED);
    checkPolylineLayer(inputStrs, expectedStrs, layerOptions, builderOptions);
    layerOptions.setEdgeType(EdgeType.UNDIRECTED);
    checkPolylineLayer(inputStrs, expectedStrs, layerOptions, builderOptions);
  }

  /**
   * A set of (edge string, IntVector) pairs representing the InputEdgeIds attached to the edges of
   * a graph. Edges are in S2TextFormat.toString() format, such as "1:3, 4:5".
   */
  private static class EdgeInputEdgeIds extends AbstractList<Pair<String, IntVector>> {
    private final List<Pair<String, IntVector>> content;

    public EdgeInputEdgeIds() {
      content = new ArrayList<>();
    }

    /** Constructor that takes ownership of the provided list. */
    public EdgeInputEdgeIds(List<Pair<String, IntVector>> list) {
      content = list;
    }

    @Override
    public int size() {
      return content.size();
    }

    @CanIgnoreReturnValue
    @Override
    public boolean add(Pair<String, IntVector> entry) {
      return content.add(entry);
    }

    @Override
    public Pair<String, IntVector> get(int index) {
      return content.get(index);
    }

    @Override
    public String toString() {
      StringBuilder sb = new StringBuilder();
      sb.append("EdgeInputEdgeIds {\n");
      for (Pair<String, IntVector> pair : content) {
        sb.append("  Output Edge \"").append(pair.first)
            .append("\" : with InputEdgeIds ").append(pair.second).append("\n");
      }
      return sb.append("}").toString();
    }

    public boolean contains(Pair<String, IntVector> entry) {
      // Cannot simply use "content.contains(entry)" because that relies on Objects.equal(), and
      // IntVector doesn't override Object.equal.
      for (Pair<String, IntVector> p : content) {
        if (p.first.contentEquals(entry.first) && p.second.isEqualTo(entry.second)) {
          return true;
        }
      }
      return false;
    }
  }

  private static final S2Error.Code INPUT_EDGE_ID_MISMATCH = S2Error.Code.USER_DEFINED_START;

  /**
   * A layer that checks actual input edge ids on the S2BuilderGraph vs. expected input edge ids.
   */
  static class InputEdgeIdCheckingLayer implements S2BuilderLayer {
    private final EdgeInputEdgeIds expected;
    private final GraphOptions graphOptions;

    public InputEdgeIdCheckingLayer(EdgeInputEdgeIds expected, GraphOptions graphOptions) {
      this.expected = expected;
      this.graphOptions = graphOptions;
    }

    @Override
    public GraphOptions graphOptions() {
      return graphOptions;
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      EdgeInputEdgeIds actual = new EdgeInputEdgeIds();
      for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
        S2Point src = g.vertex(g.edgeSrcId(edgeId));
        S2Point dst = g.vertex(g.edgeDstId(edgeId));
        String edge = S2TextFormat.s2PointsToString(ImmutableList.of(src, dst));
        IdSet ids = g.inputEdgeIds(edgeId);
        actual.add(Pair.of(edge, IntVector.copyOf(ids)));
      }

      // This comparison doesn't consider multiplicity, but that's fine.
      StringBuilder missing = new StringBuilder();
      StringBuilder extra = new StringBuilder();
      StringBuilder matched = new StringBuilder();
      for (Pair<String, IntVector> pair : expected) {
        if (actual.contains(pair)) {
          matched.append(toString(pair.first, pair.second));
          continue;
        }
        missing.append(toString(pair.first, pair.second));
      }
      for (Pair<String, IntVector> pair : actual) {
        if (expected.contains(pair)) {
          matched.append(toString(pair.first, pair.second));
          continue;
        }
        extra.append(toString(pair.first, pair.second));
      }

      if (missing.length() > 0 || extra.length() > 0) {
        error.init(INPUT_EDGE_ID_MISMATCH, "\nMatched:\n%sMissing:\n%sExtra:\n%s\n",
            matched, missing, extra);
      }
      return error.ok();
    }

    private String toString(String first, IntVector second) {
      StringBuilder sb = new StringBuilder();
      sb.append("  (").append(first).append(")={");
      if (!second.isEmpty()) {
        second.forEach(id -> sb.append(id).append(", "));
        sb.delete(sb.length() - 2, sb.length()); // remove trailing ", "
      }
      sb.append("}\n");
      return sb.toString();
    }
  }

  /**
   * Adds polylines in the given inputStrs to to an S2Builder created from the given options,
   * and checks that the edges and input edge ids attached to those edges in the resulting graph
   * match the given expected edges and input edge ids.
  */
  private void checkInputEdgeIds(
      List<String> inputStrs,
      EdgeInputEdgeIds expected,
      GraphOptions graphOptions,
      S2Builder.Builder options) {
    S2Builder builder = options.build();
    builder.startLayer(new InputEdgeIdCheckingLayer(expected, graphOptions));
    for (String inputStr : inputStrs) {
      builder.addPolyline(makePolylineOrDie(inputStr));
    }
    assertBuildOk(builder);
  }

  /**
   * Runs build() on the given builder and asserts that there was no error, printing the error
   * message otherwise.
   */
  static void assertBuildOk(S2Builder builder) {
    S2Error error = new S2Error();
    boolean built = builder.build(error);
    assertTrue(error.toString(), built);
  }

  /**
   * Chooses a random S2Point that is often near the intersection of one of the coordinate planes or
   * coordinate axes with the unit sphere. (It is possible to represent very small perturbations
   * near such points.)
   */
  private S2Point choosePoint() {
    double[] x = new double[3];
    data.getRandomPoint().fill(x);
    for (int i = 0; i < 3; ++i) {
      if (data.oneIn(3)) {
        x[i] *= pow(1e-50, data.nextDouble());
      }
    }
    return new S2Point(x).normalize();
  }

  private void checkSnappingWithForcedVertices(
      String inputStr, S1Angle snapRadius, String verticesStr, String expectedStr) {
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
    S2Builder builder = options.build();
    List<S2Point> vertices = parsePointsOrDie(verticesStr);
    for (S2Point vertex : vertices) {
      builder.forceVertex(vertex);
    }
    S2PolylineLayer layer = new S2PolylineLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie(inputStr));

    assertBuildOk(builder);
    S2Polyline output = layer.getPolyline();
    assertEquals(expectedStr, S2TextFormat.toString(output));
  }

  private void assertPolygonsEqual(S2Polygon expected, S2Polygon actual) {
    assertEquals(
        "\nExpected:\n"
            + S2TextFormat.toString(expected)
            + "\nActual:\n"
            + S2TextFormat.toString(actual),
        expected,
        actual);
  }

  private void assertPolygonsApproxEqual(S2Polygon expected, S2Polygon actual, S1Angle tolerance) {
    assertTrue(
        "\nExpected:  "
            + S2TextFormat.toString(expected)
            + "\nActual:    "
            + S2TextFormat.toString(actual)
            + "\nTolerance: "
            + tolerance.degrees(),
        expected.boundaryApproxEquals(actual, tolerance.radians()));
  }

  private void assertPolylinesEqual(S2Polyline expected, S2Polyline actual) {
    assertEquals(
        "\nExpected:\n"
            + S2TextFormat.toString(expected)
            + "\nActual:\n"
            + S2TextFormat.toString(actual),
        expected,
        actual);
  }
}
