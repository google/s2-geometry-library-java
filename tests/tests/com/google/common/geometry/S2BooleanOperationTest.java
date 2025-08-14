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

import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.makePolygonOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2BooleanOperation.OpType;
import com.google.common.geometry.S2BooleanOperation.PolygonModel;
import com.google.common.geometry.S2BooleanOperation.PolylineModel;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.S2BuilderUtil.IndexMatchingLayer;
import com.google.common.geometry.S2BuilderUtil.IndexedLayer;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for S2BooleanOperation. */
@RunWith(JUnit4.class)
public final class S2BooleanOperationTest extends GeometryTestCase {

  /** Tests constructing a S2BooleanOperation with a Builder and Options. */
  @Test
  public void testConstructors() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    // Verify that two of the options have the correct default values.
    S2BooleanOperation.Options options = builder.options();
    Preconditions.checkState(
        options.polylineModel() == S2BooleanOperation.Options.DEFAULT_POLYLINE_MODEL);
    Preconditions.checkState(
        options.polygonModel() == S2BooleanOperation.Options.DEFAULT_POLYGON_MODEL);

    // Modify the Builder state, verify the change.
    builder.setPolygonModel(PolygonModel.OPEN);
    builder.setPolylineModel(PolylineModel.OPEN);
    options = builder.options();
    Preconditions.checkState(options.polylineModel() == PolylineModel.OPEN);
    Preconditions.checkState(options.polygonModel() == PolygonModel.OPEN);

    // Create a new Builder from the options, verify it has the same values.
    S2BooleanOperation.Builder copiedBuilder = new S2BooleanOperation.Builder(options);
    options = copiedBuilder.options();
    Preconditions.checkState(options.polylineModel() == PolylineModel.OPEN);
    Preconditions.checkState(options.polygonModel() == PolygonModel.OPEN);

    // Modify the state of copiedBuilder and verify the modification.
    copiedBuilder.setPolygonModel(PolygonModel.CLOSED);
    copiedBuilder.setPolylineModel(PolylineModel.CLOSED);
    options = copiedBuilder.options();
    Preconditions.checkState(options.polylineModel() == PolylineModel.CLOSED);
    Preconditions.checkState(options.polygonModel() == PolygonModel.CLOSED);

    // Verify that the original builder is unchanged, independent of copiedBuilder.
    options = builder.options();
    Preconditions.checkState(options.polylineModel() == PolylineModel.OPEN);
    Preconditions.checkState(options.polygonModel() == PolygonModel.OPEN);

    // Construct a new S2BooleanOperation.
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
    S2BooleanOperation unused = builder.build(OpType.UNION, layer);
  }

  // TODO(ericv): Clean up or remove these notes.
  //
  // Options to test:
  //   polygonModel:                   OPEN, SEMI_OPEN, CLOSED
  //   polylineModel:                  OPEN, SEMI_OPEN, CLOSED
  //   polylineLoopsHaveBoundaries:  true, false
  //   conservative:                    true, false
  //
  // Geometry combinations to test:
  //
  // Point/point:
  //  - disjoint, coincident
  // Point/polyline:
  //  - Start vertex, end vertex, interior vertex, degenerate polyline
  //  - With polylineLoopsHaveBoundaries: start/end vertex, degenerate polyline
  // Point/polygon:
  //  - Polygon interior, exterior, vertex
  //  - Vertex of degenerate sibling pair shell, hole
  //  - Vertex of degenerate single point shell, hole
  // Polyline/polyline:
  //  - Vertex intersection:
  //    - Start, end, interior, degenerate, loop start/end, degenerate loop
  //    - Test cases where vertex is not emitted because an incident edge is.
  //  - Edge/edge: interior crossing, duplicate, reversed, degenerate
  //  - Test that degenerate edges are ignored unless polyline has a single edge.
  //    (For example, AA has one edge but AAA has no edges.)
  // Polyline/polygon:
  //  - Vertex intersection: polyline vertex cases already covered, but test
  //    polygon normal vertex, sibling pair shell/hole, single vertex shell/hole
  //    - Also test cases where vertex is not emitted because an edge is.
  //  - Edge/edge: interior crossing, duplicate, reversed
  //  - Edge/interior: polyline edge in polygon interior, exterior
  // Polygon/polygon:
  //  - Vertex intersection:
  //    - normal vertex, sibling pair shell/hole, single vertex shell/hole
  //    - Also test cases where vertex is not emitted because an edge is.
  //    - Test that polygons take priority when there is a polygon vertex and
  //      also isolated polyline vertices. (There should not be any points.)
  //  - Edge/edge: interior crossing, duplicate, reversed
  //  - Interior/interior: polygons in interior/exterior of other polygons

  /** Verify that degenerate polylines are preserved under all boundary models. */
  @Test
  public void testDegeneratePolylines() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    String a = "# 0:0, 0:0 #";
    String b = "# #";
    builder.setPolylineModel(PolylineModel.OPEN);
    expectResult(OpType.UNION, builder, a, b, a);
    builder.setPolylineModel(PolylineModel.SEMI_OPEN);
    expectResult(OpType.UNION, builder, a, b, a);
    builder.setPolylineModel(PolylineModel.CLOSED);
    expectResult(OpType.UNION, builder, a, b, a);
  }

  /**
   * Verify that degenerate polygon features (single-vertex and sibling pair shells and holes) are
   * preserved under all boundary models.
   */
  @Test
  public void testDegeneratePolygons() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    String a = "# # 0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 3:3; 6:6; 7:7, 8:8";
    String b = "# #";
    builder.setPolygonModel(PolygonModel.OPEN);
    expectResult(OpType.UNION, builder, a, b, a);
    builder.setPolygonModel(PolygonModel.SEMI_OPEN);
    expectResult(OpType.UNION, builder, a, b, a);
    builder.setPolygonModel(PolygonModel.CLOSED);
    expectResult(OpType.UNION, builder, a, b, a);
  }

  @Test
  public void testPointPoint() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    String a = "0:0 | 1:0 # #";
    String b = "0:0 | 2:0 # #";
    // Note that these results have duplicates, which is correct. Clients can eliminate the
    // duplicates with the appropriate GraphOptions.
    expectResult(OpType.UNION, builder, a, b, "0:0 | 0:0 | 1:0 | 2:0 # #");
    expectResult(OpType.INTERSECTION, builder, a, b, "0:0 | 0:0 # #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "1:0 # #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "1:0 | 2:0 # #");
  }

  /** Tests operations between an open polyline and its vertices. */
  @Test
  public void testPointOpenPolyline() {
    // The polyline "3:0, 3:0" consists of a single degenerate edge and contains no points (since
    // polylineModel() is OPEN). Since S2BooleanOperation preserves degeneracies, this means that
    // the union includes *both* the point 3:0 and the degenerate polyline 3:0, since they do not
    // intersect.
    //
    // This test uses Options.polylineLoopsHaveBoundaries() == true, which means that the loop
    // "4:0, 5:0, 4:0" does not contain the vertex "4:0".
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.OPEN);
    String a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
    String b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "0:0 | 2:0 | 3:0 | 4:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
    expectResult(OpType.INTERSECTION, builder, a, b, "1:0 | 5:0 # #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "0:0 | 2:0 | 3:0 | 4:0 # #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "0:0 | 2:0 | 3:0 | 4:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
  }

  @Test
  public void testPointOpenPolylineLoopBoundariesFalse() {
    // With Options.polylineLoopsHaveBoundaries() == false, the loop "4:0, 5:0, 4:0" has two
    // vertices, both of which are contained.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.OPEN);
    builder.setPolylineLoopsHaveBoundaries(false);
    String a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
    String b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "0:0 | 2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
    expectResult(OpType.INTERSECTION, builder, a, b, "1:0 | 4:0 | 5:0 # #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "0:0 | 2:0 | 3:0 # #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "0:0 | 2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
  }

  @Test
  public void testPointSemiOpenPolyline() {
    // Degenerate polylines are defined not contain any points under the SEMI_OPEN model either, so
    // again the point 3:0 and the degenerate polyline "3:0, 3:0" do not intersect.
    //
    // The result does not depend on Options.polylineLoopsHaveBoundaries().
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.SEMI_OPEN);
    for (boolean haveBoundaries : ImmutableList.of(false, true)) {
      builder.setPolylineLoopsHaveBoundaries(haveBoundaries);
      String a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
      String b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
      expectResult(
          OpType.UNION, builder, a, b, "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
      expectResult(OpType.INTERSECTION, builder, a, b, "0:0 | 1:0 | 4:0 | 5:0 # #");
      expectResult(OpType.DIFFERENCE, builder, a, b, "2:0 | 3:0 # #");
      expectResult(
          OpType.SYMMETRIC_DIFFERENCE,
          builder,
          a,
          b,
          "2:0 | 3:0 # 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
    }
  }

  @Test
  public void testPointClosedPolyline() {
    // Under the CLOSED model, the degenerate polyline 3:0 does contain its vertex. Since polylines
    // take precedence over points, the union of the point 3:0 and the polyline 3:0 is the polyline
    // only. Similarly, since subtracting a point from a polyline has no effect, the symmetric
    // difference includes only the polyline objects.
    //
    // The result does not depend on Options.polylineLoopsHaveBoundaries().
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.CLOSED);
    for (boolean haveBoundaries : ImmutableList.of(false, true)) {
      builder.setPolylineLoopsHaveBoundaries(haveBoundaries);
      String a = "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #";
      String b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #";
      expectResult(OpType.UNION, builder, a, b, "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
      expectResult(OpType.INTERSECTION, builder, a, b, "0:0 | 1:0 | 2:0 | 3:0 | 4:0 | 5:0 # #");
      expectResult(OpType.DIFFERENCE, builder, a, b, "# #");
      expectResult(
          OpType.SYMMETRIC_DIFFERENCE,
          builder,
          a,
          b,
          "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0, 4:0 #");
    }
  }

  @Test
  public void testLineTouchingPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setSnapFunction(new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));

    // A line "parallel to" and 1e-15 degrees north of the equator. (Not technically parallel)
    String a = "# 0.000000000000001:-1, 0.000000000000001:1 #";
    // A triangle with the southern point at 0:0. The line above crosses the tip of the southern
    // point of the triangle, _very slightly_ north of 1e-15 because the geodesic curves away from
    // the equator.
    String b = "# # 1:-1, 0:0, 1:1";

    // After snapping to the INTERSECTION_MERGE_RADIUS, the tiny bit of line intersecting the tip of
    // the triangle ends up as a degenerate line / single point:
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "#     1.000152328043908e-15:-1.0000000000000005e-15, "
            + "1.000152328043908e-15:-1.0000000000000005e-15#");

    // Before the IndexMatchingLayer had a matching tolerance, this test failed: the actual lat/lng
    // was so close to the expected one that it was printed the same, but was not identical. A more
    // precise expected result cannot be specified using S2TextFormat.
  }

  @Test
  public void testPointPolygonInterior() {
    // PolygonModel is irrelevant.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    // One interior point and one exterior point.
    String a = "1:1 | 4:4 # #";
    String b = "# # 0:0, 0:3, 3:0";
    expectResult(OpType.UNION, builder, a, b, "4:4 # # 0:0, 0:3, 3:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "1:1 # #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "4:4 # #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "4:4 # # 0:0, 0:3, 3:0");
  }

  @Test
  public void testPointOpenPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.OPEN);
    // See notes about the two vertices below.
    String a = "0:1 | 1:0 # #";
    String b = "# # 0:0, 0:1, 1:0";
    expectResult(OpType.UNION, builder, a, b, "0:1 | 1:0 # # 0:0, 0:1, 1:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "0:1 | 1:0 # #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "0:1 | 1:0 # # 0:0, 0:1, 1:0");
  }

  @Test
  public void testPointSemiOpenPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.SEMI_OPEN);
    // The two vertices are chosen such that the polygon contains one vertex but not the other under
    // PolygonModel.SEMI_OPEN. (The same vertices are used for all three PolygonModel builder.)
    S2Polygon polygon = makePolygonOrDie("0:0, 0:1, 1:0");
    assertTrue(polygon.contains(makePointOrDie("0:1")));
    assertFalse(polygon.contains(makePointOrDie("1:0")));
    String a = "0:1 | 1:0 # #";
    String b = "# # 0:0, 0:1, 1:0";
    expectResult(OpType.UNION, builder, a, b, "1:0 # # 0:0, 0:1, 1:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "0:1 # #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "1:0 # #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "1:0 # # 0:0, 0:1, 1:0");
  }

  @Test
  public void testPointClosedPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.CLOSED);
    // See notes about the two vertices above.
    String a = "0:1 | 1:0 # #";
    String b = "# # 0:0, 0:1, 1:0";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:1, 1:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "0:1 | 1:0 # #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# # 0:0, 0:1, 1:0");
  }

  @Test
  public void testPolylineVertexOpenPolylineVertex() {
    // Test first, last, and middle vertices of both polylines. Also test first/last and middle
    // vertices of two polyline loops.
    //
    // Degenerate polylines are tested in PolylineEdgePolylineEdgeOverlap below.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.OPEN);
    String a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
    String b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

    // The output consists of the portion of each input polyline that intersects the opposite
    // region, so the intersection vertex is present twice. This allows reassembling the individual
    // polylines that intersect, if desired. (Otherwise duplicates can be removed using
    // DuplicateEdges.MERGE.)
    expectResult(OpType.INTERSECTION, builder, a, b, "# 0:1, 0:1 | 0:1, 0:1 #");

    // Note that all operations are defined such that subtracting a lower-dimensional subset of an
    // object has no effect. In this case, subtracting the middle vertex of a polyline has no
    // effect.
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
  }

  @Test
  public void testPolylineVertexOpenPolylineVertexLoopBoundariesFalse() {
    // With Options.polylineLoopsHaveBoundaries() == false, the 3 polyline loops each have two
    // vertices, both of which are contained.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.OPEN);
    builder.setPolylineLoopsHaveBoundaries(false);
    String a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
    String b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

    // Note that the polyline "0:3, 0:4, 0:3" only has two vertices, not three. This means that 0:3
    // is emitted only once for that polyline, plus once for the other polyline, for a total of
    // twice.
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 0:1, 0:1 | 0:1, 0:1 | 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #");

    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
  }

  @Test
  public void testPolylineVertexSemiOpenPolylineVertex() {
    // The result does not depend on Options.polylineLoopsHaveBoundaries().
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.SEMI_OPEN);
    for (boolean haveBoundaries : ImmutableList.of(false, true)) {
      builder.setPolylineLoopsHaveBoundaries(haveBoundaries);
      String a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
      String b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
      expectResult(
          OpType.UNION,
          builder,
          a,
          b,
          "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
              + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
      expectResult(
          OpType.INTERSECTION,
          builder,
          a,
          b,
          "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 "
              + "| 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #");
      expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
      expectResult(
          OpType.SYMMETRIC_DIFFERENCE,
          builder,
          a,
          b,
          "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
              + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
    }
  }

  @Test
  public void testPolylineVertexClosedPolylineVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.CLOSED);
    String a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
    String b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

    // Since Options.polylineLoopsHaveBoundaries() == true, the polyline "0:3, 0:4, 0:3" has three
    // vertices. Therefore 0:3 is emitted twice for that polyline, plus once for the other polyline,
    // for a total of thrice.
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 "
            + "| 0:2, 0:2 | 0:2, 0:2 "
            + "| 0:3, 0:3 | 0:3, 0:3 | 0:3, 0:3 "
            + "| 0:4, 0:4 | 0:4, 0:4 | 0:4, 0:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
  }

  @Test
  public void testPolylineVertexClosedPolylineVertexLoopBoundariesFalse() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.CLOSED);
    builder.setPolylineLoopsHaveBoundaries(false);
    String a = "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #";
    String b = "# 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");

    // Since Options.polylineLoopsHaveBoundaries() == false, the polyline "0:3, 0:4, 0:3" has two
    // vertices. Therefore 0:3 is emitted once for that polyline, plus once for the other polyline,
    // for a total of twice.
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 0:0, 0:0 | 0:0, 0:0 | 0:1, 0:1 | 0:1, 0:1 "
            + "| 0:2, 0:2 | 0:2, 0:2 "
            + "| 0:3, 0:3 | 0:3, 0:3 | 0:4, 0:4 | 0:4, 0:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 0:1, 0:2 | 0:3, 0:4, 0:3 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:0, 0:1, 0:2 | 0:0, 1:0 | -1:1, 0:1, 1:1 | -1:2, 0:2 "
            + "| 0:3, 0:4, 0:3 | 1:3, 0:3, 1:3 | 0:4, 1:4, 0:4 #");
  }

  // The polygon used in the polyline/polygon vertex tests below.
  private static final String VERTEX_TEST_POLYGON_STR =
      "0:0, 0:1, 0:2, 0:3, 0:4, 0:5, 5:5, 5:4, 5:3, 5:2, 5:1, 5:0";

  @Test
  public void testTestSemiOpenPolygonVerticesContained() {
    // Verify whether certain vertices of the test polygon are contained under the semi-open
    // boundary model (for use in the tests below).
    S2Polygon polygon = makePolygonOrDie(VERTEX_TEST_POLYGON_STR);
    assertTrue(polygon.contains(makePointOrDie("0:1")));
    assertTrue(polygon.contains(makePointOrDie("0:2")));
    assertTrue(polygon.contains(makePointOrDie("0:3")));
    assertTrue(polygon.contains(makePointOrDie("0:4")));
    assertFalse(polygon.contains(makePointOrDie("5:1")));
    assertFalse(polygon.contains(makePointOrDie("5:2")));
    assertFalse(polygon.contains(makePointOrDie("5:3")));
    assertFalse(polygon.contains(makePointOrDie("5:4")));
  }

  // Don't bother testing every PolylineModel with every PolygonModel for vertex intersection, since
  // we have already tested the PolylineModels individually above. It is sufficient to use
  // PolylineModel.CLOSED with the various PolygonModel builder.
  @Test
  public void testPolylineVertexOpenPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.OPEN);

    // Define some constants to reduce code duplication. Test all combinations of polylines that
    // start or end on a polygon vertex, where the polygon vertex is open or closed using semi-open
    // boundaries, and where the incident edge is inside or outside the polygon.
    String a =
        "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
            + "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
    String b = "# # " + VERTEX_TEST_POLYGON_STR;

    String differenceResult =
        "# 0:1, 0:1 | 0:2, 0:2 | -1:3, 0:3 | 0:4, -1:4"
            + "| 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #";
    expectResult(OpType.UNION, builder, a, b, differenceResult + VERTEX_TEST_POLYGON_STR);
    expectResult(
        OpType.INTERSECTION, builder, a, b, "# 1:1, 0:1 | 0:2, 1:2 | 4:3, 5:3 | 5:4, 4:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, differenceResult);
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, differenceResult + VERTEX_TEST_POLYGON_STR);
  }

  // Like the test above, except that every polygon vertex is also incident to a closed polyline
  // vertex. This tests that when an open vertex and a closed vertex coincide with each other, the
  // result is considered closed.
  @Test
  public void testPolylineVertexOpenPolygonClosedPolylineVertex() {
    String testGeometrySuffix =
        "-2:0, 0:1 | -2:1, 0:2 | -2:2, 0:3 | -2:3, 0:4 | "
            + "7:0, 5:1 | 7:1, 5:2 | 7:2, 5:3 | 7:3, 5:4 # "
            + VERTEX_TEST_POLYGON_STR;

    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.OPEN);
    String a =
        "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
            + "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
    String b = "# " + testGeometrySuffix;

    String differencePrefix = "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        differencePrefix + " | 0:1, 0:1 | 0:2, 0:2 | 5:3, 5:3 | 5:4, 5:4 | " + testGeometrySuffix);
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4"
            + "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4"
            + "| 0:1, 0:1 | 0:2, 0:2 | 0:3, 0:3 | 0:4, 0:4"
            + "| 5:1, 5:1 | 5:2, 5:2 | 5:3, 5:3 | 5:4, 5:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, differencePrefix + " #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, differencePrefix + " | " + testGeometrySuffix);
  }

  @Test
  public void testPolylineVertexSemiOpenPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.SEMI_OPEN);
    // Test all combinations of polylines that start or end on a polygon vertex, where the polygon
    // vertex is open or closed using semi-open boundaries, and where the incident edge is inside or
    // outside the polygon.
    //
    // The vertices at latitude 0 used below are all closed while the vertices at latitude 5 are all
    // open (see TestSemiOpenPolygonVerticesContained).
    String a =
        "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
            + "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
    String b = "# # " + VERTEX_TEST_POLYGON_STR;
    String differenceResult =
        "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 | 5:3, 5:3 | 5:4, 5:4 #";
    expectResult(OpType.UNION, builder, a, b, differenceResult + VERTEX_TEST_POLYGON_STR);
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4 | 4:3, 5:3 | 5:4, 4:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, differenceResult);
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, differenceResult + VERTEX_TEST_POLYGON_STR);
  }

  @Test
  public void testPolylineVertexClosedPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.CLOSED);
    // Test all combinations of polylines that start or end on a polygon vertex, where the polygon
    // vertex is open or closed using semi-open boundaries, and where the incident edge is inside or
    // outside the polygon.
    String a =
        "# 1:1, 0:1 | 0:2, 1:2 | -1:3, 0:3 | 0:4, -1:4 "
            + "| 6:1, 5:1 | 5:2, 6:2 | 4:3, 5:3 | 5:4, 4:4 #";
    String b = "# # " + VERTEX_TEST_POLYGON_STR;
    String differenceResult = "# -1:3, 0:3 | 0:4, -1:4 | 6:1, 5:1 | 5:2, 6:2 #";
    expectResult(OpType.UNION, builder, a, b, differenceResult + VERTEX_TEST_POLYGON_STR);
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 1:1, 0:1 | 0:2, 1:2 | 0:3, 0:3 | 0:4, 0:4"
            + "| 5:1, 5:1 | 5:2, 5:2 | 4:3, 5:3 | 5:4, 4:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, differenceResult);
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, differenceResult + VERTEX_TEST_POLYGON_STR);
  }

  @Test
  public void testPolylineEdgePolylineEdgeCrossing() {
    // Two polyline edges that cross at a point interior to both edges.
    S2BooleanOperation.Builder builder = roundToE(1);
    String a = "# 0:0, 2:2 #";
    String b = "# 2:0, 0:2 #";
    expectResult(OpType.UNION, builder, a, b, "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 1:1, 1:1 | 1:1, 1:1 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 1:1, 2:2 #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 0:0, 1:1, 2:2 | 2:0, 1:1, 0:2 #");
  }

  @Test
  public void testPolylineEdgePolylineEdgeOverlap() {
    // The PolylineModel does not affect this calculation. In particular the intersection of a
    // degenerate polyline edge with itself is non-empty, even though the edge contains no points in
    // the OPEN and SEMI_OPEN models.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.OPEN);
    // Test edges in the same and reverse directions, and degenerate edges.
    String a = "# 0:0, 1:0, 2:0, 2:5 | 3:0, 3:0 | 6:0, 5:0, 4:0 #";
    String b = "# 0:0, 1:0, 2:0 | 3:0, 3:0 | 4:0, 5:0 #";
    // As usual, the expected output includes the relevant portions of *both* input polylines.
    // Duplicates can be removed using GraphOptions.
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 1:0, 2:0, 2:5 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 "
            + "| 6:0, 5:0, 4:0 | 4:0, 5:0 #");
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 0:0, 1:0, 2:0 | 0:0, 1:0, 2:0 | 3:0, 3:0 | 3:0, 3:0 | 5:0, 4:0 | 4:0, 5:0 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 2:0, 2:5 | 6:0, 5:0 #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 2:0, 2:5 | 6:0, 5:0 #");
  }

  @Test
  public void testPolylineLoopMultipleOpenPolylineEdge() {
    // Here we test a polyline loop ABCA with the pairs {AA, AB} and {AA, AC}. This tests not only
    // what happens when degenerate polylines intersect loop endpoints, but also what happens when
    // polylines intersect a degenerate and non-degenerate edge that overlap each other.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.OPEN);
    String a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
    String b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 "
            + "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
    expectResult(
        OpType.INTERSECTION, builder, a, b, "# 0:0, 0:1 | 0:0, 0:1 | 2:2, 3:2 | 3:2, 2:2 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:1, 1:0, 0:0 | 0:0, 0:0 | 2:2, 2:3, 3:2 | 2:2, 2:2 #");
  }

  @Test
  public void testPolylineLoopMultipleSemiOpenPolylineEdge() {
    // Like the test above but with SEMI_OPEN boundaries. In this case ABCA intersected with {AA,
    // AB} is {AA, AB, AB} but ABCA intersected with {AA, AC} is {AA, AA, AC, CA} since the chain
    // ABCA contains its start vertex but not its end vertex.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.SEMI_OPEN);
    String a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
    String b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 "
            + "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 0:0, 0:0 | 0:0, 0:1 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 | 3:2, 2:2 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
  }

  @Test
  public void testPolylineLoopMultipleClosedPolylineEdge() {
    // Like the test above but with CLOSED boundaries. In this case ABCA intersected with {AA, AB}
    // is {AA, AA, AB, AB} since the chain ABCA contains both its start vertex and end vertex.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.CLOSED);
    String a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
    String b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 "
            + "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# 0:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 | 0:0, 0:1 "
            + "| 2:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 | 3:2, 2:2 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
  }

  @Test
  public void testPolylineLoopMultiplePolylineEdgeLoopBoundariesFalse() {
    // Like the tests above but with polylineLoopsHaveBoundaries() == false. In this case the result
    // does not depend on the polyline model. The polyline AA intersects ABCA exactly once, and the
    // intersection of ABCA with {AA, AB} is {AA, AB, AB}.
    for (PolylineModel polylineModel :
        ImmutableList.of(PolylineModel.OPEN, PolylineModel.SEMI_OPEN, PolylineModel.CLOSED)) {
      S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
      builder.setPolylineModel(polylineModel);
      builder.setPolylineLoopsHaveBoundaries(false);
      String a = "# 0:0, 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2, 2:2 #";
      String b = "# 0:0, 0:0 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 #";
      expectResult(
          OpType.UNION,
          builder,
          a,
          b,
          "# 0:0, 0:1, 1:0, 0:0 | 0:0, 0:0 | 0:0, 0:1 "
              + "| 2:2, 2:3, 3:2, 2:2 | 2:2, 2:2 | 2:2, 3:2 #");
      expectResult(
          OpType.INTERSECTION,
          builder,
          a,
          b,
          "# 0:0, 0:0 | 0:0, 0:1 | 0:0, 0:1 | 2:2, 2:2 | 2:2, 3:2 | 3:2, 2:2 #");
      expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
      expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 0:1, 1:0, 0:0 | 2:2, 2:3, 3:2 #");
    }
  }

  @Test
  public void testPolylineEdgeOpenPolygonEdgeOverlap() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.OPEN);
    // A polygon and two polyline edges that coincide with the polygon boundary, one in the same
    // direction and one in the reverse direction.
    String a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
    String b = "# # 1:1, 1:3, 3:3, 3:1";
    expectResult(OpType.UNION, builder, a, b, "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
    expectResult(OpType.INTERSECTION, builder, a, b, "# #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 1:1, 1:3, 3:3 | 3:3, 1:3 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 1:1, 1:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
  }

  @Test
  public void testPolylineEdgeSemiOpenPolygonEdgeOverlap() {
    S2Polygon polygon = makePolygonOrDie("1:1, 1:3, 3:3, 3:1");
    assertFalse(polygon.contains(makePointOrDie("1:1")));
    assertTrue(polygon.contains(makePointOrDie("1:3")));
    assertFalse(polygon.contains(makePointOrDie("3:3")));
    assertFalse(polygon.contains(makePointOrDie("3:1")));
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.SEMI_OPEN);
    String a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
    String b = "# # 1:1, 1:3, 3:3, 3:1";
    expectResult(
        OpType.UNION, builder, a, b, "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 1:3, 1:3 | 1:1, 1:3, 3:3 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 1:1, 1:1 | 3:3, 3:3 | 3:3, 1:3 # 1:1, 1:3, 3:3, 3:1");
  }

  @Test
  public void testPolylineEdgeClosedPolygonEdgeOverlap() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.CLOSED);
    String a = "# 1:1, 1:3, 3:3 | 3:3, 1:3 # ";
    String b = "# # 1:1, 1:3, 3:3, 3:1";
    expectResult(OpType.UNION, builder, a, b, "# # 1:1, 1:3, 3:3, 3:1");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 1:1, 1:3, 3:3 | 3:3, 1:3 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# # 1:1, 1:3, 3:3, 3:1");
  }

  @Test
  public void testPolygonVertexMatching() {
    // This test shows that CrossingProcessor.processEdgeCrossings() must set a0MatchesPolygon and
    // a1MatchesPolygon correctly even when (a0, a1) itself is a polygon edge (or its sibling). (It
    // requires degenerate polygon geometry to demonstrate this.)
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolylineModel(PolylineModel.CLOSED);
    builder.setPolygonModel(PolygonModel.CLOSED);
    String a = "# 0:0, 1:1 # ";
    String b = "# # 0:0, 1:1";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 1:1");
  }

  @Test
  public void testPolylineEdgePolygonInterior() {
    // PolygonModel is irrelevant.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    // One normal and one degenerate polyline edge in the polygon interior, and similarly for the
    // polygon exterior.
    String a = "# 1:1, 2:2 | 3:3, 3:3 | 6:6, 7:7 | 8:8, 8:8 # ";
    String b = "# # 0:0, 0:5, 5:5, 5:0";
    expectResult(OpType.UNION, builder, a, b, "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 1:1, 2:2 | 3:3, 3:3 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 6:6, 7:7 | 8:8, 8:8 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 6:6, 7:7 | 8:8, 8:8 # 0:0, 0:5, 5:5, 5:0");
  }

  @Test
  public void testPolylineEdgeIsolatedStartVertexPlusInteriorCrossing() {
    // Tests a polyline XYZ that when intersected with a polygon results in an isolated vertex X
    // plus a clipped portion UYV. This case is unusual because the isolated vertex is handled by
    // creating a separate degenerate S2Builder input edge XX which is added before the actual edge
    // XY, and the crossing edge information needs to be associated with XY rather than XX in order
    // for GraphEdgeClipper to be able to do its work properly. The test is constructed such that if
    // the crossings are incorrectly associated with the degenerate edge XX then not only will the
    // output be incorrect, it will also trigger an internal assertion failure.
    S2BooleanOperation.Builder builder = roundToE(1);
    String a = "# 0:0, 0:10, 0:4 # "; // Polyline XYZ
    String b = "# # 0:0, -5:5, 5:5";
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 0:0 | 0:5, 0:10, 0:5 #");
  }

  @Test
  public void testPolygonEdgeIsolatedStartVertexPlusInteriorCrossing() {
    // Similar to the case above, but tests a polygon XYZ rather than a polyline. This requires
    // using the CLOSED polygon model and computing the intersection with a clockwise loop rather
    // than subtracting a CCW loop. The test is constructed such that if the crossings for the edge
    // 0:0, 0:8 are incorrectly associated with the degenerate edge 0:0, then not only will the
    // output be incorrect, it will also trigger an internal assertion failure.
    S2BooleanOperation.Builder builder = roundToE(1);
    builder.setPolygonModel(PolygonModel.CLOSED);
    String a = "# # 0:0, 5:5, -5:5";
    String b = "# # 1:4, 0:0, 0:8";
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 0:0; 0:5, 0:8, 0.8:5");
  }

  @Test
  public void testPolygonVertexOpenPolygonVertex() {
    // Two polygons that share a vertex.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.OPEN);
    String a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
    String b = "# # 0:0, 5:3, 5:2";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
    expectResult(OpType.INTERSECTION, builder, a, b, "# #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
  }

  @Test
  public void testPolygonVertexSemiOpenPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.SEMI_OPEN);
    String a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
    String b = "# # 0:0, 5:3, 5:2";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
    expectResult(OpType.INTERSECTION, builder, a, b, "# #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
  }

  @Test
  public void testPolygonVertexClosedPolygonVertex() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.CLOSED);
    String a = "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5";
    String b = "# # 0:0, 5:3, 5:2";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 0:0");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5");
    expectResult(OpType.DIFFERENCE, builder, b, a, "# # 0:0, 5:3, 5:2");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:5, 1:5, 0:0, 2:5, 3:5, 0:0, 5:3, 5:2");
  }

  @Test
  public void testPolygonEdgePolygonEdgeCrossing() {
    // Two polygons whose edges cross at points interior to both edges.
    S2BooleanOperation.Builder builder = roundToE(2);
    String a = "# # 0:0, 0:2, 2:2, 2:0";
    String b = "# # 1:1, 1:3, 3:3, 3:1";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:2, 1:2, 1:3, 3:3, 3:1, 2:1, 2:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 1:1, 1:2, 2:2, 2:1");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:0; 1:2, 1:3, 3:3, 3:1, 2:1, 2:2");
  }

  @Test
  public void testPolygonEdgeOpenPolygonEdgeOverlap() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    // One shape is a rectangle, the other consists of one triangle inside the rectangle and one
    // triangle outside the rectangle, where each triangle shares one edge with the rectangle. This
    // implies that the edges are in the same direction in one case and opposite directions in the
    // other case.
    builder.setPolygonModel(PolygonModel.OPEN);
    String a = "# # 0:0, 0:4, 2:4, 2:0";
    String b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0; 0:4, 1:5, 2:4");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 0:0, 1:1, 2:0");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0, 1:1");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
  }

  @Test
  public void testPolygonEdgeSemiOpenPolygonEdgeOverlap() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.SEMI_OPEN);
    String a = "# # 0:0, 0:4, 2:4, 2:0";
    String b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:4, 1:5, 2:4, 2:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 0:0, 1:1, 2:0");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0, 1:1");
    // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are normalized, i.e. the
    // output could contain siblings pairs (which can be discarded using S2Builder.GraphOptions).
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
  }

  @Test
  public void testPolygonEdgeClosedPolygonEdgeOverlap() {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.CLOSED);
    String a = "# # 0:0, 0:4, 2:4, 2:0";
    String b = "# # 0:0, 1:1, 2:0; 0:4, 1:5, 2:4";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:4, 1:5, 2:4, 2:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 0:0, 1:1, 2:0; 0:4, 2:4");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0, 1:1");
    // Note that SYMMETRIC_DIFFERENCE does not guarantee that results are normalized, i.e. the
    // output could contain siblings pairs.
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 2:4, 2:0, 1:1; 0:4, 1:5, 2:4");
  }

  @Test
  public void testPolygonPolygonInterior() {
    // PolygonModel is irrelevant.
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    // One loop in the interior of another polygon and one loop in the exterior.
    String a = "# # 0:0, 0:4, 4:4, 4:0";
    String b = "# # 1:1, 1:2, 2:2, 2:1; 5:5, 5:6, 6:6, 6:5";
    expectResult(OpType.UNION, builder, a, b, "# # 0:0, 0:4, 4:4, 4:0; 5:5, 5:6, 6:6, 6:5");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 1:1, 1:2, 2:2, 2:1");
    // Difference cuts a hole out of the big polygon, reversing the order of vertices.
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:4, 4:4, 4:0; 2:1, 2:2, 1:2, 1:1; 5:5, 5:6, 6:6, 6:5");
  }

  @Test
  public void testPolygonEdgesDegenerateAfterSnapping() {
    S2BooleanOperation.Builder builder = roundToE(0);
    // Two narrow rectangles forming a plus sign.
    String a = "# # 0:-1, 0:1, 0.1:1, 0.1:-1";
    String b = "# # -1:0.1, 1:0.1, 1:0, -1:0";
    // When snapping causes an output edge to become degenerate, it is still emitted (since
    // otherwise loops that contract to a single point would be lost). If the output layer doesn't
    // want such edges, they can be removed via DegenerateEdges.DISCARD.
    expectResult(OpType.UNION, builder, a, b, "# # 0:-1, 0:0, 0:1, 0:0 | -1:0, 0:0, 1:0, 0:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 0:0");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# # 0:-1, 0:0, 0:1, 0:0");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:-1, 0:0, 0:1, 0:0 | -1:0, 0:0, 1:0, 0:0");
  }

  // These are the characters representing possible input regions, in the order they are used in
  // the rules listed further below.
  private static final String INPUT_CHARS = ".pPudsSUDHh*";

  /** Returns a String containing the single char 'c' */
  private String str(char c) {
    char[] a = {c};
    return new String(a);
  }

  /**
   * Verifies that the S2BooleanOperation results for the given OpType and PolygonModel match the
   * given set of rules (encoded as described below).
   *
   * <p>The unit tests below based on runDegeneracyTest comprehensively test the handling of
   * degenerate geometry. They check over 4,000 individual cases encoded as a series of textual
   * tables.
   *
   * <p>The inputs to the test cases are intended to span all possible types of degenerate and
   * non-degenerate edge configurations in the vicinity of an individual input edge. Essentially we
   * are trying to cover all possible behaviors of
   * S2BooleanOperation.Impl.CrossingProcessor.ProcessEdge(), which is responsible for determining
   * which parts of each input edge should be emitted to the output.
   *
   * <p>For this purpose it is sufficient to consider two points A and B and the possible types of
   * degeneracies that involve just these two points. The actual locations of the points are
   * immaterial, although for descriptive purposes we suppose that B is located above A so that the
   * edge AB is considered "up" while the edge BA is considered "down". Only one point needs to be
   * considered for point degeneracies, so we use A for this purpose. This means that all inputs
   * consist of some subset of AA, AB, and BA represented as edges of dimensions 0, 1, or 2. For
   * example, a polyline edge intersected with a point located at one of its vertices would be
   * represented as the edge AB of dimension 1 in region 0 and the edge AA of dimension 0 in region
   * 1.
   *
   * <p>Each possibility is represented as a single letter as follows (with the corresponding edges
   * in braces):
   *
   * <pre>
   * Special:        . = empty {}
   *                 * = full sphere {}
   *
   * Dimension 0:    p = point {AA}
   *
   * Dimension 1:    P = point polyline {AA}
   *                 u = up edge {AB}
   *                 d = down edge {BA}
   *
   * Dimension 2:    s = point shell {AA}
   *                 S = sibling pair shell {AB, BA}
   *                 U = up edge {AB}
   *                 D = down edge {BA}
   *                 H = sibling pair hole {AB, BA}
   *                 h = point hole {AA}
   * </pre>
   *
   * <p>Using this encoding, the test case described above would be represented as 'u' for region 0
   * and 'p' for region 1. Note that while the test cases focus only on what happens to the edges
   * AB, BA, and AA, additional edges are needed in order to construct valid inputs. For example,
   * the test case U (the "up" polygon edge AB) is represented as a triangle ABC, while the test
   * case D (the "down" polygon edge BA) is represented as the triangle AB(-C).
   *
   * <p>The expected result of a given operation is then represented as a sequence of the characters
   * above. In some cases additional symbols are needed to define the expected results. For example,
   * in the closed polygon model the union of U and D is the quadrilateral A(-C)BC which does not
   * appear in the list of symbols above. So the results use the following additional symbols:
   *
   * <ul>
   *   <li>Q : the union of the "U" and "D" shapes (a quadrilateral)
   *   <li>B : a point polyline consisting only of the vertex B {BB}
   * </ul>
   *
   * <p>Finally, in some cases the expected result depends on the polyline model, or on whether one
   * of the two regions contains the vertex A. For example, the intersection of 'p' and 'U' in the
   * semi-open model depends on whether the representation of U (the triangle ABC mentioned above)
   * contains its vertex A. Conditional results of this sort are encoded using the following
   * operators:
   *
   * <ul>
   *   <li>~X : the complement of X (where X must be either U or D)
   *   <li>X<Y : the result is X if region 0 contains point A, otherwise Y.
   *   <li>X>Y : the result is X if region 1 contains point A, otherwise Y.
   *   <li>X|Y|Z : the result is X, Y, or Z depending on whether the polyline model is OPEN,
   *       SEMI_OPEN, or CLOSED.
   * </ul>
   *
   * <p>The operators above may be combined, e.g. X<>Y means the result is X if both regions contain
   * point A and Y otherwise.
   */
  private void runDegeneracyTest(OpType opType, PolygonModel polygonModel, List<String> rules) {
    final boolean verbose = false;

    Preconditions.checkArgument(rules.size() == INPUT_CHARS.length());

    // For the symmetric operators (i.e., all except difference) we only need to define half of the
    // expected results.
    final boolean symmetric = (opType != OpType.DIFFERENCE);

    // Possible values for Options.polylineModel().
    final ImmutableList<PolylineModel> polylineModels =
        ImmutableList.of(PolylineModel.OPEN, PolylineModel.SEMI_OPEN, PolylineModel.CLOSED);

    // Possible values for Options.polylineLoopsHaveBoundaries().
    final ImmutableList<Boolean> polylineLoopOptions = ImmutableList.of(true, false);

    // The set of characters representing polyline inputs.
    final String lineChars = "Pud";

    // The following nested loops iterate over all combinations of:
    //  - Input character 0
    //  - Input character 1
    //  - polyline model (if either input is a polyline)
    //  - whether a closed polyline loop is considered to have a boundary
    //    (if either input is a degenerate polyline)
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(polygonModel);
    for (int i = 0; i < INPUT_CHARS.length(); ++i) {
      char ch0 = INPUT_CHARS.charAt(i);
      List<String> row = strSplitSkipEmpty(" ", rules.get(i));
      // Verify and remove the row label.
      assertEquals(row.get(0), str(ch0));
      assertEquals("|", row.get(1));
      row = row.subList(2, row.size());
      int limit = symmetric ? (i + 1) : INPUT_CHARS.length();
      assertEquals(limit, row.size());
      for (int j = 0; j < limit; ++j) {
        char ch1 = INPUT_CHARS.charAt(j);
        // Only iterate over polyline models if at least one input is a polyline.
        int numLineModels = (lineChars.contains(str(ch0)) || lineChars.contains(str(ch1))) ? 3 : 1;
        for (int k = 0; k < numLineModels; ++k) {
          builder.setPolylineModel(polylineModels.get(k));
          // Only iterate over polyline loop boundary builder if at least one input is a
          // degenerate polyline.
          int numPolylineLoopOptions = (ch0 == 'P' || ch1 == 'P') ? 2 : 1;
          for (int m = 0; m < numPolylineLoopOptions; ++m) {
            builder.setPolylineLoopsHaveBoundaries(polylineLoopOptions.get(m));
            if (verbose) {
              System.out.println(Platform.formatString(
                  "!! RunDegeneracyTest[%s, '%c', '%c', polylineModel=%s/%d, "
                      + "polylineLoopsHaveBoundaries=%s/%d]:%n",
                  opType,
                  ch0,
                  ch1,
                  polylineModels.get(k),
                  numLineModels,
                  polylineLoopOptions.get(m),
                  numPolylineLoopOptions));
            }

            String code = row.get(j);
            // Process any '<' or '>' operators in the result.
            List<String> choices = strSplitSkipEmpty("[<>]", code);
            String result = choices.get(0);
            if (choices.size() > 1) {
              assertEquals(2, choices.size());
              // Test whether each input contains the point A. Note that we can't use
              // S2ContainsPointQuery because the containment test must be done using the given
              // S2BooleanOperation builder.
              boolean in0 =
                  S2BooleanOperation.contains(
                      makeIndex(str(ch0)), makeIndex("p"), builder.options());
              boolean in1 =
                  S2BooleanOperation.contains(
                      makeIndex(str(ch1)), makeIndex("p"), builder.options());
              // If the point containment conditions specified by the operators '<' and '>' are
              // satisfied then the result is choices[0], otherwise it is choices[1].
              if ((code.contains("<") && !in0) || (code.contains(">") && !in1)) {
                result = choices.get(1);
              }
            }
            // Next process any '|' operators in the result.
            choices = strSplit("\\|", result);
            if (choices.size() > 1) {
              assertEquals(
                  Platform.formatString(
                      "split result %s to get Choices.length=%d:(%s, %s). Expected numLineModels=3"
                          + " but was %d. ch0=%s, ch1=%s",
                      result,
                      choices.size(),
                      choices.get(0),
                      choices.get(1),
                      numLineModels,
                      ch0,
                      ch1),
                  3,
                  numLineModels);
              assertEquals(3, choices.size());
              result = choices.get(k);
            }
            checkResult(opType, builder, ch0, ch1, result);
            if (symmetric && j != i) {
              checkResult(opType, builder, ch1, ch0, result);
            }
            if (verbose) {
              System.out.println("!! Done.");
            }
          }
        }
      }
    }
  }

  private List<String> strSplitSkipEmpty(String separator, String input) {
    String[] split = input.split(separator, -1);
    List<String> output = new ArrayList<>();
    for (int j = 0; j < split.length; ++j) {
      if (split[j].isEmpty()) {
        continue;
      }
      output.add(split[j]);
    }
    return output;
  }

  private List<String> strSplit(String separator, String input) {
    String[] split = input.split(separator, -1);
    List<String> output = new ArrayList<>();
    for (int j = 0; j < split.length; ++j) {
      output.add(split[j]);
    }
    return output;
  }

  /** Concise conversion of a variable number of points into an immutable list. */
  private static ImmutableList<S2Point> loop(S2Point... points) {
    return ImmutableList.copyOf(points);
  }

  /** Concise conversion of a variable number of lists of points into an immutable list of lists. */
  private static ImmutableList<ImmutableList<S2Point>> loops(ImmutableList<S2Point>... loops) {
    return ImmutableList.copyOf(loops);
  }

  /**
   * Verifies that S2BooleanOperation returns the given result for the inputs represented by the
   * characters 'ch0' and 'ch1'.
   */
  private void checkResult(
      OpType opType, S2BooleanOperation.Builder builder, char ch0, char ch1, String result) {
    S2ShapeIndex index0 = makeIndex(str(ch0));
    S2ShapeIndex index1 = makeIndex(str(ch1));
    S2ShapeIndex expected = makeIndex(result);
    S2BooleanOperation.Options options = builder.options();

    StringBuilder sb = new StringBuilder();
    sb.append("\n  Operation: ").append(opType.toString());
    sb.append("\n  Options: ").append(options.toString());
    sb.append("\n  PolygonModel: ").append(options.polygonModel().toString());
    sb.append("\n  PolylineModel: ").append(options.polylineModel().toString());
    sb.append("\n  PolylineLoopsHaveBoundaries: ").append(options.polylineLoopsHaveBoundaries());
    sb.append("\n  Inputs: ").append(ch0).append(", ").append(ch1);
    sb.append("\n  Expected: ").append(result);

    expectResult(opType, builder, index0, index1, expected, sb.toString());
  }

  /** Returns true if the given shape contains the given point. */
  private boolean containsBruteForce(S2Shape shape, S2Point point) {
    if (shape.dimension() < 2) {
      return false;
    }

    S2Shape.ReferencePoint refPoint = shape.getReferencePoint();
    if (refPoint.equalsPoint(point)) {
      return refPoint.contained();
    }

    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(refPoint.point(), point);
    boolean inside = refPoint.contained();
    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
    for (int e = 0; e < shape.numEdges(); ++e) {
      shape.getEdge(e, edge);
      inside ^= crosser.edgeOrVertexCrossing(edge.getStart(), edge.getEnd());
    }
    return inside;
  }

  // Returns an S2ShapeIndex corresponding to the given string of characters.
  private S2ShapeIndex makeIndex(String chars) {
    // The locations of A, B, C are arbitrary, however some tests are sensitive as to whether
    // certain polygons contain the points {A, B} or not. If these points are moved then some
    // test results may need to be changed.
    S2Point a = new S2Point(1, 0, 0);
    S2Point b = new S2Point(0, 0, 1);
    S2Point c = new S2Point(0, 1, 0);

    S2ShapeIndex index = new S2ShapeIndex();
    for (int i = 0; i < chars.length(); ++i) {
      char ch = chars.charAt(i);
      if (ch == '.') { // Empty
      } else if (ch == 'p') { // Point
        index.add(S2Point.Shape.fromList(loop(a)));
      } else if (ch == 'P') { // Polyline consisting only of the point A
        index.add(S2LaxPolylineShape.create(loop(a, a)));
      } else if (ch == 'B') { // Polyline consisting only of the point B
        index.add(S2LaxPolylineShape.create(loop(b, b)));
      } else if (ch == 'u') { // Upwards polyline edge
        index.add(S2LaxPolylineShape.create(loop(a, b)));
      } else if (ch == 'd') { // Downwards polyline edge
        index.add(S2LaxPolylineShape.create(loop(b, a)));
      } else if (ch == 's') { // Point shell
        index.add(S2LaxPolygonShape.fromLoop(loop(a)));
      } else if (ch == 'S') { // Sibling pair shell
        index.add(S2LaxPolygonShape.create(loops(loop(a, b))));
      } else if (ch == 'U') { // Upwards polygon edge
        S2LaxPolygonShape laxPolygon = S2LaxPolygonShape.create(loops(loop(a, b, c.neg())));
        index.add(laxPolygon);
        // Some test results require that the U polygon contains A but not B.
        Preconditions.checkState(containsBruteForce(laxPolygon, a));
        Preconditions.checkState(!containsBruteForce(laxPolygon, b));
      } else if (ch == 'D') { // Downwards polygon edge
        S2LaxPolygonShape laxPolygon = S2LaxPolygonShape.create(loops(loop(b, a, c)));
        index.add(laxPolygon);
        // Some test cases require that the D polygon excludes both A and B.
        Preconditions.checkState(!containsBruteForce(laxPolygon, a));
        Preconditions.checkState(!containsBruteForce(laxPolygon, b));
      } else if (ch == '~') { // Complement of following region (U or D)
        ch = chars.charAt(++i);
        if (ch == 'U') {
          index.add(S2LaxPolygonShape.create(loops(loop(c.neg(), b, a))));
        } else if (ch == 'D') {
          index.add(S2LaxPolygonShape.create(loops(loop(c, a, b))));
        } else {
          throw new IllegalStateException("Unsupported character for ~ operator: " + str(ch));
        }
      } else if (ch == 'Q') { // Union of 'U' and 'D' shapes
        index.add(S2LaxPolygonShape.create(loops(loop(a, c, b, c.neg()))));
      } else if (ch == 'H') { // Sibling pair hole
        index.add(S2LaxPolygonShape.create(loops(loop(a, b), loop())));
      } else if (ch == 'h') { // Point hole
        index.add(S2LaxPolygonShape.create(loops(loop(a), loop())));
      } else if (ch == '*') { // Full sphere, as an S2LaxPolygonShape with one loop with no points.
        index.add(S2LaxPolygonShape.create(loops(loop())));
      } else {
        throw new IllegalStateException("Unknown degeneracy coverage symbol: " + str(ch));
      }
    }
    return index;
  }

  @Test
  public void testDegeneracyOpenIntersection() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  .    pp",
            "P |  .    p<.   PP",
            "u |  .    p<.   PP<.  uu",
            "d |  .    p<.   PP<.  ud    dd",
            "s |  .     .     .     .     .     s",
            "S |  .     .     .     .     .     .     S",
            "U |  .     .     .     .     .     .     .     U",
            "D |  .     .     .     .     .     .     .     .     D",
            "H |  .     .     .     .     .     .     .     U     D     H",
            "h |  .     .     .     u     d     .     S     U     D     H     h",
            "* |  .     p     P     u     d     s     S     U     D     H     h     *");
    runDegeneracyTest(OpType.INTERSECTION, PolygonModel.OPEN, rules);
  }

  @Test
  public void testDegeneracySemiOpenIntersection() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  .    pp",
            "P |  .    p<.   PP",
            "u |  .    p<.   PP<.  uu",
            "d |  .    p<.   PP<.  ud    dd",
            "s |  .     .     .     .     .     s",
            "S |  .    p<.   P<.    .     .    s<.    S",
            "U |  .    p<.   P<.    u    P<>.  s<.    .     U",
            "D |  .    p<.   P<.   P<>.   d    s<.    .     .     D",
            "H |  .    p<.   P<.    u     d    s<.    .     U     D     H",
            "h |  .     p     P     u     d     .     S     U     D     H     h",
            "* |  .     p     P     u     d     s     S     U     D     H     h     *");
    runDegeneracyTest(OpType.INTERSECTION, PolygonModel.SEMI_OPEN, rules);
  }

  @Test
  public void testDegeneracyClosedIntersection() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  .    pp",
            "P |  .    p<.   PP",
            "u |  .    p<.   PP<.  uu",
            "d |  .    p<.   PP<.  ud    dd",
            "s |  .     p     P    P>.   P>.    s",
            "S |  .     p     P     u     d     s     S",
            "U |  .     p     P     u     d     s     S     U",
            "D |  .     p     P     u     d     s     S     S     D",
            "H |  .     p     P     u     d     s     S     U     D     H",
            "h |  .     p     P     u     d     s     S     U     D     H     h",
            "* |  .     p     P     u     d     s     S     U     D     H     h     *");
    runDegeneracyTest(OpType.INTERSECTION, PolygonModel.CLOSED, rules);
  }

  @Test
  public void testDegeneracyOpenUnion() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  p    pp",
            "P |  P   P<Pp   PP",
            "u |  u   u<up   Pu    uu",
            "d |  d   d<dp   Pd    ud    dd",
            "s |  s    ps    Ps    us    ds     s",
            "S |  S    pS    PS    uS    dS     S     S",
            "U |  U    pU    PU    uU    dU     U     U     U",
            "D |  D    pD    PD    uD    dD     D     D    UD     D",
            "H |  H    pH    PH    uH    dH     H     H     H     H     H",
            "h |  h    ph    Ph   Ph>h  Ph>h    h     h     h     h     h     h",
            "* |  *     *     *     *     *     *     *     *     *     *     *     *");
    runDegeneracyTest(OpType.UNION, PolygonModel.OPEN, rules);
  }

  @Test
  public void testDegeneracySemiOpenUnion() {
    // CAVEAT: The results for (U,u) and (D,d) require the U polygon to contain vertex A but not
    // vertex B, and the D polygon to contain neither vertex. This differs from most of the other
    // tests, which encode the results conditionally using the '<' and '>' operators. That was not
    // practical in this case because (1) no conditional operators are defined for the 'B' vertex
    // and (2) encoding the full set of possibilities for all 12 cases (i.e., the 3 polyline models
    // and whether U contains A and/or B) would be unwieldy.
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  p    pp",
            "P |  P   P<Pp   PP",
            "u |  u   u<up   Pu    uu",
            "d |  d   d<dp   Pd    ud    dd",
            "s |  s    ps    Ps    us    ds     s",
            "S |  S   S<pS  S<PS   uS    dS     S     S",
            "U |  U   U<pU  U<PU U|U|BU  dU     U     U     U",
            "D |  D   D<pD  D<PD   uD  D|BD|PBD D     D     Q     D",
            "H |  H   H<pH  H<PH    H     H     H     *     *     *     H",
            "h |  h     h     h     h     h     *    *>h   *>h   *>h   *>h    h",
            "* |  *     *     *     *     *     *     *     *     *     *     *     *");
    runDegeneracyTest(OpType.UNION, PolygonModel.SEMI_OPEN, rules);
  }

  @Test
  public void testDegeneracyClosedUnion() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  p    pp",
            "P |  P   P<Pp   PP",
            "u |  u   u<up   Pu    uu",
            "d |  d   d<dp   Pd    ud    dd",
            "s |  s     s     s    us    ds     s",
            "S |  S     S     S     S     S     S     S",
            "U |  U     U     U     U     U     U     U     U",
            "D |  D     D     D     D     D     D     D     Q     D",
            "H |  H     H     H     H     H     H     *     *     *     H",
            "h |  h     h     h     h     h     *     *     *     *     *     h",
            "* |  *     *     *     *     *     *     *     *     *     *     *     *");
    runDegeneracyTest(OpType.UNION, PolygonModel.CLOSED, rules);
  }

  @Test
  public void testDegeneracyOpenDifference() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .     .     .     .     .     .     .     .     .     .     .     .",
            "p |  p     .    .>p   .>p   .>p    p     p     p     p     p     p     .",
            "P |  P     P     .    .>P   .>P    P     P     P     P     P     P     .",
            "u |  u     u     u     .   .|P|.   u     u     u     u     u    P<.    .",
            "d |  d     d     d   .|B|.   .     d     d     d     d     d    P<.    .",
            "s |  s     s     s     s     s     .     s     s     s     s     s     .",
            "S |  S     S     S     S     S     S     .     S     S     S     .     .",
            "U |  U     U     U     U     U     U     U     .     U     .     .     .",
            "D |  D     D     D     D     D     D     D     D     .     .     .     .",
            "H |  H     H     H     H     H     H     H    ~U    ~D     .     .     .",
            "h |  h     h     h     h     h     h     H    ~U    ~D     S     .     .",
            "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .");
    runDegeneracyTest(OpType.DIFFERENCE, PolygonModel.OPEN, rules);
  }

  @Test
  public void testDegeneracySemiOpenDifference() {
    // See SemiOpenUnion notes regarding (u,U) and (d,D).
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .     .     .     .     .     .     .     .     .     .     .     .",
            "p |  p     .    .>p   .>p   .>p    p     p    .>p   .>p    .     .     .",
            "P |  P     P     .    .>P   .>P    P     P    .>P   .>P    .     .     .",
            "u |  u     u     u     .   .|P|.   u     u   .|.|B   u     .     .     .",
            "d |  d     d     d   .|B|.   .     d     d     d   .|B|PB  .     .     .",
            "s |  s     s     s     s     s     .    .>s   .>s   .>s   .>s    s     .",
            "S |  S     S     S     S     S     S     .     .     .     S    s<.    .",
            "U |  U     U     U     U     U     U     U     .     U     .    s<.    .",
            "D |  D     D     D     D     D     D     D     D     .     .    s<.    .",
            "H |  H     H     H     H     H     H     H    ~U    ~D     .    s<.    .",
            "h |  h     h     h     h     h     h     H    ~U    ~D     S     .     .",
            "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .");
    runDegeneracyTest(OpType.DIFFERENCE, PolygonModel.SEMI_OPEN, rules);
  }

  @Test
  public void testDegeneracyClosedDifference() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .     .     .     .     .     .     .     .     .     .     .     .",
            "p |  p     .    .>p   .>p   .>p    .     .     .     .     .     .     .",
            "P |  P     P     .    .>P   .>P    .     .     .     .     .     .     .",
            "u |  u     u     u     .   .|P|.   u     .     .     .     .     .     .",
            "d |  d     d     d   .|B|.   .     d     .     .     .     .     .     .",
            "s |  s     s     s     s     s     .     .     .     .     .     s     .",
            "S |  S     S     S     S     S     S     .     .     .     S     .     .",
            "U |  U     U     U     U     U     U     U     .     U     .     .     .",
            "D |  D     D     D     D     D     D     D     D     .     .     .     .",
            "H |  H     H     H     H     H     H     H    ~U    ~D     .     .     .",
            "h |  h     h     h     h     h     h     H    ~U    ~D     S     .     .",
            "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .");
    runDegeneracyTest(OpType.DIFFERENCE, PolygonModel.CLOSED, rules);
  }

  @Test
  public void testDegeneracyOpenSymmetricDifference() {
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  p     .",
            "P |  P   P<Pp    .",
            "u |  u   u<up  u<uP    .",
            "d |  d   d<dp  d<dP .|PB|.  .",
            "s |  s    sp    sP    su    sd     .",
            "S |  S    Sp    SP    Su    Sd     S     .",
            "U |  U    Up    UP    Uu    Ud     U     U     .",
            "D |  D    Dp    DP    Du    Dd     D     D    UD     .",
            "H |  H    Hp    HP    Hu    Hd     H     H    ~U    ~D     .",
            "h |  h    hp    hP   hP>h  hP>h    h     H    ~U    ~D     S     .",
            "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .");
    runDegeneracyTest(OpType.SYMMETRIC_DIFFERENCE, PolygonModel.OPEN, rules);
  }

  @Test
  public void testDegeneracySemiOpenSymmetricDifference() {
    // See SemiOpenUnion notes regarding (U,u) and (D,d).
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  p     .",
            "P |  P   P<Pp    .",
            "u |  u   u<up  u<uP    .",
            "d |  d   d<dp  d<dP .|PB|.  .",
            "s |  s    sp    sP    su    sd     .",
            "S |  S    Sp    SP    Su    Sd     S     .",
            "U |  U   U<Up  U<UP U|U|UB  Ud     U     U     .",
            "D |  D   D<Dp  D<DP   Du  D|BD|PBD D     D    UD     .",
            "H |  H     H     H     H     H     H     H    ~U    ~D     .",
            "h |  h     h     h     h     h     h     H    ~U    ~D     S     .",
            "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .");
    runDegeneracyTest(OpType.SYMMETRIC_DIFFERENCE, PolygonModel.SEMI_OPEN, rules);
  }

  @Test
  public void testDegeneracyClosedSymmetricDifference() {
    // Note that (H,S)->H, (h,s)->h and (U,D)->UD. In all three cases the shared boundary is present
    // on both sides and therefore these edges should not be contained by the result, however this
    // is not possible under the CLOSED model. The indicated results are the best approximation.
    ImmutableList<String> rules =
        ImmutableList.of(
            //    .    p     P     u     d     s     S     U     D     H     h     *
            // |-----------------------------------------------------------------------
            ". |  .",
            "p |  p     .",
            "P |  P   P<Pp    .",
            "u |  u   u<up  u<uP    .",
            "d |  d   d<dp  d<dP .|PB|.  .",
            "s |  s     s     s    su    sd     .",
            "S |  S     S     S     S     S     S     .",
            "U |  U     U     U     U     U     U     U     .",
            "D |  D     D     D     D     D     D     D    UD     .",
            "H |  H     H     H     H     H     H     H    ~U    ~D     .",
            "h |  h     h     h     h     h     h     H    ~U    ~D     S     .",
            "* |  *     *     *     *     *     h     H    ~U    ~D     S     s     .");
    runDegeneracyTest(OpType.SYMMETRIC_DIFFERENCE, PolygonModel.CLOSED, rules);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // The remaining tests are intended to cover combinations of features or interesting special
  // cases.

  @Test
  public void testThreeOverlappingBars() {
    // Two vertical bars and a horizontal bar that overlaps both of the other bars and connects
    // them.

    // Round intersection points to E2 precision because the expected results were computed in
    // lat/lng space rather than using geodesics.
    S2BooleanOperation.Builder builder = roundToE(2);
    String a = "# # 0:0, 0:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3";
    String b = "# # 1:1, 1:4, 2:4, 2:1";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# # 0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0");
    expectResult(OpType.INTERSECTION, builder, a, b, "# # 1:1, 1:2, 2:2, 2:1; 1:3, 1:4, 2:4, 2:3");
    expectResult(
        OpType.DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; 0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0; "
            + "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3; "
            + "1:2, 1:3, 2:3, 2:2");
  }

  @Test
  public void testFourOverlappingBars() {
    // Two vertical bars and two horizontal bars.

    // Round intersection points to E2 precision because the expected results were computed in
    // lat/lng space rather than using geodesics.
    S2BooleanOperation.Builder builder = roundToE(2);
    String a = "# # 1:88, 1:93, 2:93, 2:88; -1:88, -1:93, 0:93, 0:88";
    String b = "# # -2:89, -2:90, 3:90, 3:89; -2:91, -2:92, 3:92, 3:91";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# # -1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, "
            + "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, "
            + "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88; "
            + "0:90, 1:90, 1:91, 0:91" /*CW*/);
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# # 1:89, 1:90, 2:90, 2:89; 1:91, 1:92, 2:92, 2:91; "
            + "-1:89, -1:90, 0:90, 0:89; -1:91, -1:92, 0:92, 0:91");
    expectResult(
        OpType.DIFFERENCE,
        builder,
        a,
        b,
        "# # 1:88, 1:89, 2:89, 2:88; 1:90, 1:91, 2:91, 2:90; "
            + "1:92, 1:93, 2:93, 2:92; -1:88, -1:89, 0:89, 0:88; "
            + "-1:90, -1:91, 0:91, 0:90; -1:92, -1:93, 0:93, 0:92");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # 1:88, 1:89, 2:89, 2:88; -1:88, -1:89, 0:89, 0:88; "
            + "1:90, 1:91, 2:91, 2:90; -1:90, -1:91, 0:91, 0:90; "
            + "1:92, 1:93, 2:93, 2:92; -1:92, -1:93, 0:93, 0:92; "
            + "-2:89, -2:90, -1:90, -1:89; -2:91, -2:92, -1:92, -1:91; "
            + "0:89, 0:90, 1:90, 1:89; 0:91, 0:92, 1:92, 1:91; "
            + "2:89, 2:90, 3:90, 3:89; 2:91, 2:92, 3:92, 3:91");
  }

  @Test
  public void testOverlappingDoughnuts() {
    // Two overlapping square doughnuts whose holes do not overlap. This means that the union
    // polygon has only two holes rather than three.

    // Round intersection points to E2 precision because the expected results were computed in
    // lat/lng space rather than using geodesics.
    S2BooleanOperation.Builder builder = roundToE(1);
    String a = "# # -1:-93, -1:-89, 3:-89, 3:-93; 0:-92, 2:-92, 2:-90, 0:-90" /*CW*/;
    String b = "# # -3:-91, -3:-87, 1:-87, 1:-91; -2:-90, 0:-90, 0:-88, -2:-88" /*CW*/;
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# # -1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93; "
            + "0:-92, 2:-92, 2:-90, 1:-90, 1:-91, 0:-91; " /*CW*/
            + "-2:-90, -1:-90, -1:-89, 0:-89, 0:-88, -2:-88" /*CW*/);
    expectResult(
        OpType.INTERSECTION,
        builder,
        a,
        b,
        "# # -1:-91, -1:-90, 0:-90, 0:-91; 0:-90, 0:-89, 1:-89, 1:-90");
    expectResult(
        OpType.DIFFERENCE,
        builder,
        a,
        b,
        "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, "
            + "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; "
            + "-1:-90, -1:-89, 0:-89, 0:-90");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# # -1:-93, -1:-91, 0:-91, 0:-92, 2:-92, "
            + "2:-90, 1:-90, 1:-89, 3:-89, 3:-93; "
            + "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88,-2:-88,-2:-90,-1:-90,-1:-91; "
            + "-1:-90, -1:-89, 0:-89, 0:-90; "
            + "1:-91, 0:-91, 0:-90, 1:-90");
  }

  @Test
  public void testPolylineEnteringRectangle() {
    // A polyline that enters a rectangle very close to one of its vertices.
    S2BooleanOperation.Builder builder = roundToE(1);
    String a = "# 0:0, 2:2 #";
    String b = "# # 1:1, 1:3, 3:3, 3:1";
    expectResult(OpType.UNION, builder, a, b, "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 1:1, 2:2 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:0, 1:1 #");
    expectResult(OpType.SYMMETRIC_DIFFERENCE, builder, a, b, "# 0:0, 1:1 # 1:1, 1:3, 3:3, 3:1");
  }

  @Test
  public void testPolylineCrossingRectangleTwice() {
    // A polyline that crosses a rectangle in one direction, then moves to a different side and
    // crosses the rectangle in the other direction. Note that the input polyline has a
    // self-intersection and that an extra vertex is *not* added at the intersection point. This is
    // an important feature since it allows polylines with many self-intersections (such as GPS
    // tracks) to be manipulated without possible quadratic size increases.
    S2BooleanOperation.Builder builder = roundToE(1);
    String a = "# 0:-5, 0:5, 5:0, -5:0 #";
    String b = "# # 1:1, 1:-1, -1:-1, -1:1";
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 "
            + "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 0:-1, 0:1 | 1:0, -1:0 #");
    expectResult(
        OpType.DIFFERENCE, builder, a, b, "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:-5, 0:-1 | 0:1, 0:5, 5:0, 1:0 | -1:0, -5:0 "
            + "# 1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1");
  }

  @Test
  public void testSelfIntersectingPolylines() {
    // Two polylines that intersect at the point 2:4, and that also have self-intersections at the
    // points 2:2 and 3:4 respectively. The intersection point should always be created, but the
    // self-intersection points should be created iff splitAllCrossingPolylineEdges() is true.

    S2BooleanOperation.Builder builder = roundToE(1);
    String a = "# 0:2, 4:2, 2:0, 2:5 #";
    String b = "# 0:4, 5:4, 3:6, 3:3 #";
    expectResult(
        OpType.UNION, builder, a, b, "# 0:2, 4:2, 2:0, 2:4, 2:5 | 0:4, 2:4, 5:4, 3:6, 3:3 #");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 2:4, 2:4 | 2:4, 2:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:2, 4:2, 2:0, 2:4, 2:5 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:2, 4:2, 2:0, 2:4, 2:5 | 0:4, 2:4, 5:4, 3:6, 3:3 #");

    builder.setSplitAllCrossingPolylineEdges(true);
    expectResult(
        OpType.UNION,
        builder,
        a,
        b,
        "# 0:2, 2:2, 4:2, 2:0, 2:2, 2:4, 2:5 | 0:4, 2:4, 3:4, 5:4, 3:6, 3:4, 3:3 #");
    expectResult(OpType.INTERSECTION, builder, a, b, "# 2:4, 2:4 | 2:4, 2:4 #");
    expectResult(OpType.DIFFERENCE, builder, a, b, "# 0:2, 2:2, 4:2, 2:0, 2:2, 2:4, 2:5 #");
    expectResult(
        OpType.SYMMETRIC_DIFFERENCE,
        builder,
        a,
        b,
        "# 0:2, 2:2, 4:2, 2:0, 2:2, 2:4, 2:5 | 0:4, 2:4, 3:4, 5:4, 3:6, 3:4, 3:3 #");
  }

  /**
   * Subtracts a degenerate loop along the 180 degree meridian from the given input geometry, and
   * compares the result to "expectedStr". The inputs should be in the format expected by
   * makeIndexOrDie().
   */
  private static void checkMeridianSplitting(
      String inputStr, String expectedStr, String j2clExpectedStr) {
    S2ShapeIndex input = makeIndexOrDie(inputStr);
    S2ShapeIndex meridian = new S2ShapeIndex();
    ImmutableList<ImmutableList<S2Point>> loops =
        loops(
            loop(
                new S2Point(0, 0, -1),
                new S2Point(-1, 0, 0),
                new S2Point(0, 0, 1),
                new S2Point(-1, 0, 0)));

    meridian.add(S2LaxPolygonShape.create(loops));
    S2ShapeIndex output = new S2ShapeIndex();

    S2BuilderLayer l0 = new IndexedLayer<>(output, new S2PointVectorLayer());
    S2BuilderLayer l1 = new IndexedLayer<>(output, new S2PolylineVectorLayer());
    S2BuilderLayer l2 = new IndexedLayer<>(output, new S2LaxPolygonLayer());
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    S2BooleanOperation op = builder.build(OpType.DIFFERENCE, l0, l1, l2);
    S2Error error = new S2Error();
    assertTrue(error.toString(), op.build(input, meridian, error));
    String actual = S2TextFormat.toString(output);
    if (j2clExpectedStr != null) {
      assertTrue(
          "Actual string: "
              + actual
              + "\nJava expected: '"
              + expectedStr
              + "'\nJ2CL expected: '"
              + j2clExpectedStr
              + "'",
          expectedStr.equals(actual) || j2clExpectedStr.equals(actual));
    } else {
      assertEquals(
          "Actual string: '" + actual + "'\nExpected: '" + expectedStr + "'",
          expectedStr, actual);
    }
  }

  private static void checkMeridianSplitting(String inputStr, String expectedStr) {
    checkMeridianSplitting(inputStr, expectedStr, null);
  }

  // This test demonstrates that S2 geometry can easily be transformed such that no edge crosses the
  // 180 degree meridian, as required by formats such as GeoJSON, by simply subtracting a degenerate
  // loop that follows the 180 degree meridian. This not only splits polylines along the meridian,
  // it also inserts the necessary extra vertices at the north/south poles. (The only extra step is
  // that the vertices along the 180 degree meridian or at the poles may need to be "doubled" into
  // two vertices, one at longitude 180 and one at longitude -180, in order to match the longitudes
  // of the adjacent vertices.)
  @Test
  public void testMeridianSplitting() {
    // A line along the equator crossing the 180 degree meridian.
    checkMeridianSplitting("# 0:-160, 0:170 #", "# 0:-160, 0:180, 0:170 #");

    // The northern hemisphere.
    checkMeridianSplitting("# # 0:0, 0:120, 0:-120", "# # 90:0, 0:180, 0:-120, 0:0, 0:120, 0:180");

    // A small square that crosses the 180th meridian. Notice that one input loop is split into two
    // output loops. S2TextFormat differences cause a slightly different result on J2CL, with one
    // fewer digit past the decimal point on some values.
    checkMeridianSplitting(
        "# # 9:179, 9:-179, 10:-179, 10:179",
        // Java expected output
        "# # 9.00134850712993:180, 9:-179, 10:-179, 10.0014925269841:180; "
            + "10.0014925269841:180, 10:179, 9:179, 9.00134850712993:180",
        // J2CL expected output
        "# # 9.0013485071299:180, 9:-179, 10:-179, 10.001492526984:180; "
            + "10.001492526984:180, 10:179, 9:179, 9.0013485071299:180");

    // An annulus that crosses the 180th meridian. This turns into two shells.
    checkMeridianSplitting(
        "# # 8:178, 8:-178, 11:-178, 11:178; 9:179, 10:179, 10:-179, 9:-179",
        // Java
        "# # 10.0014925269841:180, 10:-179, 9:-179, 9.00134850712993:180, "
            + "8.00481316618607:180, 8:-178, 11:-178, 11.00654129428:180; "
            + "9.00134850712993:180, 9:179, 10:179, 10.0014925269841:180, "
            + "11.00654129428:180, 11:178, 8:178, 8.00481316618607:180",
        // J2CL
        "# # 10.001492526984:180, 10:-179, 9:-179, 9.0013485071299:180, "
            + "8.004813166186:180, 8:-178, 11:-178, 11.00654129428:180; "
            + "9.0013485071299:180, 9:179, 10:179, 10.001492526984:180, "
            + "11.00654129428:180, 11:178, 8:178, 8.004813166186:180");

    // An annulus that crosses the 180th meridian. This turns into two shells.
    checkMeridianSplitting(
        "# # 8:178, 8:-178, 11:-178, 11:178; 9:179, 10:179, 10:-179, 9:-179",
        // Java expected
        "# # 10.0014925269841:180, 10:-179, 9:-179, 9.00134850712993:180, "
            + "8.00481316618607:180, 8:-178, 11:-178, 11.00654129428:180; "
            + "9.00134850712993:180, 9:179, 10:179, 10.0014925269841:180, "
            + "11.00654129428:180, 11:178, 8:178, 8.00481316618607:180",
        // J2CL expected
        "# # 10.001492526984:180, 10:-179, 9:-179, 9.0013485071299:180, "
            + "8.004813166186:180, 8:-178, 11:-178, 11.00654129428:180; "
            + "9.0013485071299:180, 9:179, 10:179, 10.001492526984:180, "
            + "11.00654129428:180, 11:178, 8:178, 8.004813166186:180");

    // The whole world except for a small square that crosses the 180th meridian. This is a single
    // loop that visits both poles. The result is correct except that (1) +180 or -180 needs to be
    // chosen consistently with the adjacent points, and (2) each pole needs to be duplicated (once
    // with longitude -180 and once with longitude 180).
    checkMeridianSplitting(
        "# # 9:-179, 9:179, 10:179, 10:-179",
        // Java expected
        "# # 0:180, 9.00134850712993:180, 9:179, 10:179, 10.0014925269841:180, "
            + "90:0, 10.0014925269841:180, 10:-179, 9:-179, 9.00134850712993:180, "
            + "0:180, -90:0",
        // J2CL expected);
        "# # 0:180, 9.0013485071299:180, 9:179, 10:179, 10.001492526984:180, "
            + "90:0, 10.001492526984:180, 10:-179, 9:-179, 9.0013485071299:180, "
            + "0:180, -90:0");
  }

  @CanIgnoreReturnValue
  private S2LaxPolygonShape computeTestUnion(
      ImmutableList<ImmutableList<S2Point>> aLoops,
      ImmutableList<ImmutableList<S2Point>> bLoops,
      S1Angle snapRadius) {
    S2ShapeIndex a = new S2ShapeIndex();
    S2ShapeIndex b = new S2ShapeIndex();
    a.add(S2LaxPolygonShape.create(aLoops));
    b.add(S2LaxPolygonShape.create(bLoops));

    S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
    S2BooleanOperation.Builder builder =
        new S2BooleanOperation.Builder(new IdentitySnapFunction(snapRadius));
    S2BooleanOperation op = builder.build(OpType.UNION, layer);
    S2Error error = new S2Error();
    assertTrue(error.toString(), op.build(a, b, error));
    S2LaxPolygonShape result = layer.getPolygon();
    assertFalse(
        "S2Polygon A: " + S2TextFormat.toString(a) + "\nS2Polygon B: " + S2TextFormat.toString(b),
        result.isEmpty());
    return result;
  }

  @Test
  public void testGetCrossedVertexIndexBug1() {
    // This test exercises a rare special case in getCrossedVertexIndex where two crossing edge
    // chains snap to a different permutation of the same vertices. In this example one input edge
    // crosses another edge from right to left, the first edge snaps to BCD and the second snaps to
    // ABDC, and triangle BCD is CCW. Since BCD is to the right of BD, this means that the first
    // edge has not yet crossed the second at vertex B, leaving C or D as the possible crossing
    // vertices.
    ImmutableList<ImmutableList<S2Point>> aLoops =
        loops(
            loop(
                new S2Point(-0.3830643798538849, -0.7492195533420621, 0.5403070809984629),
                new S2Point(-0.3830643798552798, -0.7492195533413425, 0.5403070809984718),
                new S2Point(-0.38306437985529124, -0.7492195533413641, 0.5403070809984336),
                new S2Point(-0.38306437985389635, -0.7492195533420838, 0.5403070809984247)));
    ImmutableList<ImmutableList<S2Point>> bLoops =
        loops(
            loop(
                new S2Point(-0.3830643798539096, -0.7492195533421059, 0.5403070809983846),
                new S2Point(-0.38306437985527797, -0.749219553341342, 0.5403070809984737),
                new S2Point(-0.3830643798552794, -0.749219553341344, 0.5403070809984701),
                new S2Point(-0.38306437985391095, -0.7492195533421078, 0.540307080998381)));
    S2LaxPolygonShape unused =
        computeTestUnion(aLoops, bLoops, S2EdgeUtil.INTERSECTION_MERGE_RADIUS);
  }

  @Test
  public void testGetCrossedVertexIndexBug2() {
    // This test exercises another rare case where the crossing vertices chosen by
    // getCrossedVertexIndex() are not ordered correctly along the edge being crossed. This is
    // handled by adding extra edges to the output in order to link up the crossings in the correct
    // order.
    ImmutableList<ImmutableList<S2Point>> aLoops =
        loops(
            loop(
                new S2Point(-0.3837392878495085, -0.7477800800281974, 0.5418201831546835),
                new S2Point(-0.38373928785696076, -0.7477800800212292, 0.5418201831590226),
                new S2Point(-0.3837392878570128, -0.7477800800212469, 0.5418201831589613),
                new S2Point(-0.38373928785703426, -0.7477800800212544, 0.5418201831589358),
                new S2Point(-0.38373947205489456, -0.747780142277955, 0.5418199667802881),
                new S2Point(-0.3837394720443441, -0.74778014228782, 0.5418199667741451),
                new S2Point(-0.38373947205872994, -0.7477801422818535, 0.5418199667721912),
                new S2Point(-0.38373947218468357, -0.7477801428893031, 0.5418199658446279),
                new S2Point(-0.3837396702525171, -0.7477802104436154, 0.5418197323311432),
                new S2Point(-0.38373967023137123, -0.7477802104633304, 0.5418197323189107),
                new S2Point(-0.38373947216030285, -0.7477801429079148, 0.541819965836209),
                new S2Point(-0.3837394721708758, -0.7477801428980574, 0.5418199658423253),
                new S2Point(-0.38373947215649007, -0.747780142904024, 0.5418199658442793),
                new S2Point(-0.3837394720305386, -0.7477801422965848, 0.5418199667718262),
                new S2Point(-0.38373928783586, -0.7477800800409594, 0.5418201831467369),
                new S2Point(-0.38373928784641037, -0.7477800800310942, 0.5418201831528797),
                new S2Point(-0.3837392878357865, -0.7477800800409342, 0.5418201831468237),
                new S2Point(-0.383739287835765, -0.7477800800409267, 0.5418201831468492)));
    ImmutableList<ImmutableList<S2Point>> bLoops =
        loops(
            loop(
                new S2Point(-0.38373923813692823, -0.7477800632164362, 0.5418202415655146),
                new S2Point(-0.3837392878569364, -0.7477800800212209, 0.5418201831590512),
                new S2Point(-0.38373928784640354, -0.7477800800310694, 0.5418201831529186),
                new S2Point(-0.3837392878463879, -0.7477800800310864, 0.5418201831529065),
                new S2Point(-0.38373928784638023, -0.7477800800310945, 0.5418201831529005),
                new S2Point(-0.383739287836921, -0.7477800800412459, 0.54182018314559),
                new S2Point(-0.38373928783691913, -0.7477800800412454, 0.5418201831455919),
                new S2Point(-0.3837392878463657, -0.7477800800311077, 0.5418201831528927),
                new S2Point(-0.3837392878463733, -0.7477800800310995, 0.5418201831528985),
                new S2Point(-0.3837392878358356, -0.7477800800409511, 0.5418201831467655),
                new S2Point(-0.38373923811582744, -0.7477800632361664, 0.5418202415532288),
                new S2Point(-0.3837385765031284, -0.7477798396184077, 0.5418210187539991),
                new S2Point(-0.3837385765242292, -0.7477798395986774, 0.5418210187662849)));
    S2LaxPolygonShape unused =
        computeTestUnion(aLoops, bLoops, S2EdgeUtil.INTERSECTION_MERGE_RADIUS);
  }

  @Test
  public void testGetCrossedVertexIndexBug3() {
    // This test exercises the special case in getCrossedVertexIndex() that requires checking the
    // orientation of a loop. This is done by adding up the turning angles at each vertex, which in
    // turn involves computing the edge normals and measuring the angles between them. However in
    // this test, some of the edge normals returned by S2.robustCrossProd() used to be so small that
    // there were floating-point underflows when computing the angles between them. This was fixed
    // by implementing the long-standing TODO of making S2.robustCrossProd() actually robust.
    ImmutableList<ImmutableList<S2Point>> aLoops =
        loops(
            loop(
                new S2Point(1, 0, 2.467823483526174E-72),
                new S2Point(0.9998476951563913, 0.01745240643728351, 1.8530922845942552e-27),
                new S2Point(0.9974025970361131, 0.06988184982643786, 0.01745240643728351)));
    ImmutableList<ImmutableList<S2Point>> bLoops =
        loops(
            loop(
                new S2Point(0.9999999999999999, 2.4674476220564615e-72, 2.467823483526174E-72),
                new S2Point(0.9999999999999999, 2.883798140665744E-169, 2.467823483526174E-72),
                new S2Point(1, 2.883798140665743E-169, 2.467823483526174E-72)));
    S2LaxPolygonShape unused = computeTestUnion(aLoops, bLoops, S1Angle.ZERO);
  }

  @Test
  public void testGetCrossedVertexIndexBug4() {
    // This example tests the "special case" in getCrossedVertexIndex() in situations where two
    // edges snap to the same sequence of vertices in different orders. The first two edges (a0, a1)
    // and (b0, b1) of the following polygons cross
    // such that after snapping, the corresponding edge chains are:
    //
    //   a0 a1 -> a0 b0 b1 x a1
    //   b0 b1 -> b0 x b1
    //
    // where "x" is the computed intersection point of (a0, a1) and (b0, b1). Previously there was a
    // bug such that the two edge chains did not choose the same vertex to represent the point where
    // the two chains cross: the (a0, a1) chain chose "x" as the crossing point while the (b0, b1)
    // chain chose "b0". This has been fixed such that both chains now choose "x". (Both "x" and
    // "b1" happen to be valid choices in this example, but it is essential that both subchains make
    // the same choice.)

    // S2LatLng coordinates are not accurate enough to reproduce this example.
    ImmutableList<ImmutableList<S2Point>> aLoops =
        loops(
            loop(
                // 51.5131559470858:-0.130381523356724
                new S2Point(0.622333310659119, -0.0014161759526823048, 0.7827510746653316),
                // 51.5131892038956:-0.130404244210776
                new S2Point(0.6223328557578689, -0.0014164217071954736, 0.7827514358937983),
                makePointOrDie("51.51317:-0.1306")));
    ImmutableList<ImmutableList<S2Point>> bLoops =
        loops(
            loop(
                // 51.5131559705551:-0.13038153939079
                new S2Point(0.6223333103380959, -0.001416176126110953, 0.78275107492025),
                // 51.5131559705551:-0.130381539390786
                new S2Point(0.6223333103380959, -0.0014161761261109063, 0.7827510749202501),
                makePointOrDie("51.52:-0.12"),
                makePointOrDie("51.52:-0.14")));
    computeTestUnion(aLoops, bLoops, S1Angle.ZERO);
  }

  @Test
  public void testGetCrossedVertexIndexBug5() {
    // Yet another bizarre situation where two crossing edges snap (correctly) to a sequence of
    // vertices in different orders. Using the internal vertex numbers assigned by S2Builder, input
    // edges 3 and 12 snap to the following vertex sequences:
    //
    //   Input edge  3:  14, 8, 4, 9, 2, 5
    //   Input edge 12:   2, 7, 8, 9
    //
    // Furthermore input edge 3 crosses input edge 12 from left to right. Schematically, here is
    // what edge 12 crossing edge 3 looks like:
    //
    //   14-->--8-->--4-->--9-->--2-->--5
    //          |\         /     /
    //          \ \--->---/     /
    //           \             /
    //            \--<--7--<--/
    //
    // And here is what edge 3 crossing edge 12 looks like:
    //
    //             14-->--\   /---4->-\
    //                     \ /         \
    //          2-->--7-->--8----->-----9
    //         / \                     /
    //  5--<--/   \---------<---------/
    //
    // In both cases, the only possible choice of crossing vertex consistent with the fact that edge
    // 3 crosses edge 12 from left to right is vertex 9. Determining this requires knowing that loop
    // (9, 2, 7, 8) is clockwise (the "special case" in GetCrossedVertexIndex). The code previously
    // didn't have quite the correct test to decide when this was necessary.
    ImmutableList<ImmutableList<S2Point>> aLoops =
        loops(
            loop(
                new S2Point(0.9998476951563913, 0, 0.01745240643728351),
                new S2Point(0.9992386149554826, 0.017441774902830158, 0.03489949670250097),
                new S2Point(0.9984774386394599, 0.05232798522331314, 0.01745240643728351),
                new S2Point(0.9980211966240684, 0.034851668155187324, 0.052335956242943835)));
    ImmutableList<ImmutableList<S2Point>> bLoops =
        loops(
            loop(
                new S2Point(0.9980211966240684, 0.034851668155187324, 0.052335956242943835),
                new S2Point(0.9961969233988566, 0.052208468483931986, 0.0697564737441253),
                new S2Point(0.9980209868161543, 0.03483971497214896, 0.05234791433446786),
                new S2Point(0.9974120827677868, 0.017411821260589495, 0.0697564737441253),
                new S2Point(0.9974121921010651, 0.01741134053876882, 0.06975503041925263),
                new S2Point(0.9974121164231596, 0.01740989325235717, 0.0697564737441253),
                new S2Point(0.9998476951563912, 4.950042464556023E-16, 0.017452406437284993),
                new S2Point(0.9998476951563913, 3.7368529835165677e-16, 0.017452406437284632),
                new S2Point(0.9998476951563912, 3.3065924905014365e-16, 0.017452406437284504),
                new S2Point(0.9998476951563913, 9.906003593224202E-16, 0.017452406437284504),
                new S2Point(0.9996954135095479, 0.017449748351250485, 0.01745240643728351)),
            loop(
                new S2Point(0.9998476951563912, 3.3065924905014365e-16, 0.017452406437284504),
                new S2Point(0.9998476951563912, 3.3006856770496304e-16, 0.017452406437284504),
                new S2Point(0.9998476951563913, 0, 0.017452406437284504),
                new S2Point(0.9998476951563913, 0, 0.01745240643728351)));
    S2LaxPolygonShape unused = computeTestUnion(aLoops, bLoops, S1Angle.ZERO);
  }

  @Test
  public void testGetCrossedVertexIndexBug6() {
    // This is another test of the code in getCrossedVertexIndex() that checks whether the B
    // subchain contains an interior vertex of the A edge.
    ImmutableList<ImmutableList<S2Point>> aLoops =
        loops(
            loop(
                new S2Point(0.9987048882355846, 0.026138065586168355, 0.04365028913720582),
                new S2Point(0.9987625943414924, 0.030513215246694664, 0.0392711578586665),
                new S2Point(0.9998476951563913, 0.01745240643728351, 0),
                new S2Point(0.998782023517925, 0.03486228668443791, 0.03491547600379121),
                new S2Point(0.9987820251299122, 0.03487823687206265, 0.03489949670250097),
                new S2Point(0.9975640502598242, 0.0697564737441253, 0),
                new S2Point(0.998779795837143, 0.034883478425067296, 0.034958008531414335),
                new S2Point(0.9961969233988566, 0.052208468483931986, 0.0697564737441253),
                new S2Point(0.9984758123481388, 0.017465633646566288, 0.05235459671364581),
                new S2Point(0.9975640502598242, 0, 0.0697564737441253),
                new S2Point(0.9984767425041021, 0.017444393356200013, 0.05234393774670617),
                new S2Point(0.9984774386394599, 0.017428488520812163, 0.052335956242943835),
                new S2Point(0.9998476951563913, 0, 0.01745240643728351)),
            loop(
                new S2Point(0.9961969233988566, 0.052208468483931986, 0.0697564737441253),
                new S2Point(0.9980211966196957, 0.0348516682804046, 0.052335956242943835),
                new S2Point(0.9987605225894034, 0.030527121154938986, 0.03931301808477241),
                new S2Point(0.9987032179652688, 0.0261619324398966, 0.04367419967013944)),
            loop(
                new S2Point(0.9961969233988566, 0.052208468483931986, 0.0697564737441253),
                new S2Point(0.9961969233988566, 0.06966087492121549, 0.052335956242943835),
                new S2Point(0.9951340343707851, 0.06958655048003272, 0.0697564737441253)));
    ImmutableList<ImmutableList<S2Point>> bLoops =
        loops(
            loop(
                new S2Point(0.998022004299885, 0.034828499898458924, 0.0523359773775543),
                new S2Point(0.9986295347545738, 0, 0.052335956242943835),
                new S2Point(0.9992379306151222, 0.017455729388178846, 0.03491211153074132),
                new S2Point(0.9992385908584587, 0.017443155365764275, 0.03489949670250097),
                new S2Point(0.9992379307614709, 0.01745573778081081, 0.034912103145779166),
                new S2Point(0.9992865072388355, 0.020934110218524152, 0.0314362764933699),
                new S2Point(1, 0, 0),
                new S2Point(0.9992998780878941, 0.022418034384064717, 0.029953053064335624),
                new S2Point(0.9993140623243144, 0.02616995393092059, 0.026201876881811362),
                new S2Point(0.9998476951563913, 0.01745240643728351, 0),
                new S2Point(0.9993057332020093, 0.029072747464899757, 0.023298646837028814),
                new S2Point(0.9986295347545738, 0.052335956242943835, 1.700986599320836e-73),
                new S2Point(0.9983851827700422, 0.03834718875939572, 0.04191085705972318),
                new S2Point(0.9961969233988567, 0.05220846848393198, 0.06975647374412529)),
            loop(
                new S2Point(0.9980211966240684, 0.05230407459247085, 0.03489949670250097),
                new S2Point(0.998477438346863, 0.05232799080639758, 0.01745240643728351),
                new S2Point(0.9961964528150565, 0.05220844382168006, 0.06976321231435134),
                new S2Point(0.9961969233988566, 0.052208468483932, 0.06975647374412532),
                new S2Point(0.9961969233988566, 0.052208468483931986, 0.0697564737441253),
                new S2Point(0.9961969233988568, 0.05220846848393199, 0.06975647374412532),
                new S2Point(0.9961969233988568, 0.052208468483931986, 0.0697564737441253),
                new S2Point(0.9961969233988567, 0.05220846848393198, 0.06975647374412529)));
    S2LaxPolygonShape unused = computeTestUnion(aLoops, bLoops, S1Angle.ZERO);
  }

  /**
   * Performs the given operation and compares the result to "expectedStr". All arguments are in
   * makeLaxPolygonOrDie() format.
   */
  private void expectPolygon(
      S2BooleanOperation.OpType opType, String aStr, String bStr, String expectedStr) {
    S2ShapeIndex a = makeIndexOrDie("# # " + aStr);
    S2ShapeIndex b = makeIndexOrDie("# # " + bStr);
    S2LaxPolygonLayer.Options polygonOptions = new S2LaxPolygonLayer.Options();
    polygonOptions.setDegenerateBoundaries(S2LaxPolygonLayer.DegenerateBoundaries.DISCARD);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(polygonOptions);
    S2BooleanOperation.Builder builder =
        new S2BooleanOperation.Builder(new IdentitySnapFunction(S1Angle.degrees(1.1)));
    S2BooleanOperation op = builder.build(opType, layer);
    S2Error error = new S2Error();
    assertTrue(error.toString(), op.build(a, b, error));
    assertEquals(expectedStr, S2TextFormat.toString(layer.getPolygon()));
  }

  @Test
  public void testFullAndEmptyResults() {
    // The following constants are all in makeLaxPolygonOrDie() format.
    String kEmpty = "";
    String kFull = "full";

    // Two complementary shell/hole pairs, together with alternative shells that are slightly
    // smaller or larger than the original.
    String kShell1 = "10:0, 10:10, 20:10";
    String kHole1 = "10:0, 20:10, 10:10";
    String kShell1Minus = "11:2, 11:9, 18:9";
    String kShell1Plus = "9:-2, 9:11, 22:11";
    String kShell2 = "10:20, 10:30, 20:30";
    String kHole2 = "10:20, 20:30, 10:30";

    // The northern and southern hemispheres.
    String kNorthHemi = "0:0, 0:120, 0:-120";
    String kSouthHemi = "0:0, 0:-120, 0:120";
    // These edges deviate from kSouthHemi by slightly more than 1 degree.
    String kSouthHemiPlus = "0.5:0, 0.5:-120, 0.5:120";

    // A shell and hole that cover complementary hemispheres, such that each hemisphere intersects
    // all six S2 cube faces. There are also alternative shells that are slightly smaller or larger
    // than the original.
    String k6FaceShell1 = "0:-45, 45:0, 45:90, 0:135, -45:180, -45:-90";
    String k6FaceHole1 = "0:-45, -45:-90, -45:180, 0:135, 45:90, 45:0";
    String k6FaceShell1Minus = "-1:-45, 44:0, 44:90, -1:135, -46:180, -46:-90";
    String k6FaceShell1Plus = "1:-45, 46:0, 46:90, 1:135, -44:180, -44:-90";

    // Two complementary shell/hole pairs that are small enough so that they will disappear when the
    // snap radius chosen above is used.
    String kAlmostEmpty1 = "2:0, 2:10, 3:0";
    String kAlmostFull1 = "2:0, 3:0, 2:10";
    String kAlmostEmpty2 = "4:0, 4:10, 5:0";
    String kAlmostFull2 = "4:0, 5:0, 4:10";

    // A polygon that intersects all 6 faces such but snaps to an empty polygon.
    String k6FaceAlmostEmpty1 = k6FaceShell1Minus + "; " + k6FaceHole1;

    // Test empty UNION results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.UNION, kEmpty, kEmpty, kEmpty);
    //  - Empty due to snapping, union does not intersect all 6 cube faces.
    expectPolygon(OpType.UNION, kAlmostEmpty1, kAlmostEmpty2, kEmpty);
    //  - Empty due to snapping, union intersects all 6 cube faces.
    expectPolygon(OpType.UNION, k6FaceAlmostEmpty1, k6FaceAlmostEmpty1, kEmpty);

    // Test full UNION results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.UNION, kEmpty, kFull, kFull);
    expectPolygon(OpType.UNION, kEmpty, kFull, kFull);
    expectPolygon(OpType.UNION, kFull, kFull, kFull);
    //  - Exact result, some input edges.
    expectPolygon(OpType.UNION, kFull, kShell1, kFull);
    expectPolygon(OpType.UNION, kHole1, kHole2, kFull);
    expectPolygon(OpType.UNION, kHole1, kShell1, kFull);
    //  - Full due to snapping, almost complementary polygons.
    expectPolygon(OpType.UNION, kHole1, kShell1Minus, kFull);
    expectPolygon(OpType.UNION, k6FaceHole1, k6FaceShell1Minus, kFull);

    // Test empty INTERSECTION results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.INTERSECTION, kEmpty, kEmpty, kEmpty);
    expectPolygon(OpType.INTERSECTION, kEmpty, kFull, kEmpty);
    expectPolygon(OpType.INTERSECTION, kFull, kEmpty, kEmpty);
    //  - Exact result, inputs do not both intersect all 6 cube faces.
    expectPolygon(OpType.INTERSECTION, kEmpty, kHole1, kEmpty);
    expectPolygon(OpType.INTERSECTION, kShell1, kShell2, kEmpty);
    expectPolygon(OpType.INTERSECTION, kShell1, kHole1, kEmpty);
    //  - Exact result, inputs both intersect all 6 cube faces.
    expectPolygon(OpType.INTERSECTION, k6FaceShell1, k6FaceHole1, kEmpty);
    //  - Empty due to snapping, inputs do not both intersect all 6 cube faces.
    expectPolygon(OpType.INTERSECTION, kShell1Plus, kHole1, kEmpty);
    //  - Empty due to snapping, inputs both intersect all 6 cube faces.
    expectPolygon(OpType.INTERSECTION, k6FaceShell1Plus, k6FaceHole1, kEmpty);

    // Test full INTERSECTION results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.INTERSECTION, kFull, kFull, kFull);
    //  - Full due to snapping, almost full input polygons.
    expectPolygon(OpType.INTERSECTION, kAlmostFull1, kAlmostFull2, kFull);

    // Test empty DIFFERENCE results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.DIFFERENCE, kEmpty, kEmpty, kEmpty);
    expectPolygon(OpType.DIFFERENCE, kEmpty, kFull, kEmpty);
    expectPolygon(OpType.DIFFERENCE, kFull, kFull, kEmpty);
    //  - Exact result, first input does not intersect all 6 cube faces.
    expectPolygon(OpType.DIFFERENCE, kEmpty, kShell1, kEmpty);
    expectPolygon(OpType.DIFFERENCE, kShell1, kFull, kEmpty);
    expectPolygon(OpType.DIFFERENCE, kShell1, kShell1, kEmpty);
    expectPolygon(OpType.DIFFERENCE, kShell1, kHole2, kEmpty);
    //  - Exact result, first input intersects all 6 cube faces.
    expectPolygon(OpType.DIFFERENCE, k6FaceShell1, k6FaceShell1Plus, kEmpty);
    //  - Empty due to snapping, first input does not intersect all 6 cube faces.
    expectPolygon(OpType.DIFFERENCE, kShell1Plus, kShell1, kEmpty);
    //  - Empty due to snapping, first input intersect all 6 cube faces.
    expectPolygon(OpType.DIFFERENCE, k6FaceShell1Plus, k6FaceShell1, kEmpty);

    // Test full DIFFERENCE results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.DIFFERENCE, kFull, kEmpty, kFull);
    //  - Full due to snapping, almost full/empty input polygons.
    expectPolygon(OpType.DIFFERENCE, kAlmostFull1, kAlmostEmpty2, kFull);

    // Test empty SYMMETRIC_DIFFERENCE results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kEmpty, kEmpty, kEmpty);
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kFull, kFull, kEmpty);
    //  - Exact result, union does not intersect all 6 cube faces.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1, kShell1, kEmpty);
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kNorthHemi, kEmpty);
    //  - Exact result, union intersects all 6 cube faces. This case is only handled correctly due
    // to the kBiasTowardsEmpty heuristic.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1, k6FaceShell1, kEmpty);
    //  - Empty due to snapping, union does not intersect all 6 cube faces.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1Plus, kShell1, kEmpty);
    //  - Empty due to snapping, union intersects all 6 cube faces. This case is only handled
    // correctly due to the kBiasTowardsEmpty heuristic.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Plus, k6FaceShell1, kEmpty);
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Minus, k6FaceShell1, kEmpty);

    // Test full SYMMETRIC_DIFFERENCE results.
    //  - Exact result, no input edges.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kFull, kEmpty, kFull);
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kEmpty, kFull, kFull);
    //  - Exact result, complementary input polygons.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1, kHole1, kFull);
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kAlmostEmpty1, kAlmostFull1, kFull);
    //  - Full due to snapping, almost complementary input polygons.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kShell1Plus, kHole1, kFull);
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kAlmostFull1, kAlmostEmpty2, kFull);
    //  - Exact result, complementary hemispheres, at least one input does not intersect all 6 cube
    // faces.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kSouthHemi, kFull);
    //  - Exact result, almost complementary hemispheres, at least one input does not intersect all
    // 6 cube faces.
    expectPolygon(OpType.SYMMETRIC_DIFFERENCE, kNorthHemi, kSouthHemiPlus, kFull);

    // TODO(ericv): The following case is not currently implemented.
    //  - Full result, complementary (to within the snap radius) input polygons each with an area of
    //    approximately 2*Pi, and both polygons intersect all 6 cube faces.
    /*
      expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1, k6FaceHole1, kFull);
      expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Plus, k6FaceHole1, kFull);
      expectPolygon(OpType.SYMMETRIC_DIFFERENCE, k6FaceShell1Minus, k6FaceHole1, kFull);
    */
  }

  /**
   * Tests S2BooleanOperation.equals, which computes the symmetric difference between two geometries
   * and tests whether the result is empty.
   *
   * <p>This also indirectly tests isEmpty(), which is used to implement contains() and
   * intersects().
   */
  @Test
  public void testEquals() {
    assertTrue(testEqual("# #", "# #"));
    assertTrue(testEqual("# # full", "# # full"));

    assertFalse(testEqual("# #", "# # full"));
    assertFalse(testEqual("0:0 # #", "# #"));
    assertFalse(testEqual("0:0 # #", "# # full"));
    assertFalse(testEqual("# 0:0, 1:1 #", "# #"));
    assertFalse(testEqual("# 0:0, 1:1 #", "# # full"));
    assertFalse(testEqual("# # 0:0, 0:1, 1:0 ", "# #"));
    assertFalse(testEqual("# # 0:0, 0:1, 1:0 ", "# # full"));
  }

  /** Tests contains() on empty and full geometries. */
  @Test
  public void testContainsEmptyAndFull() {
    S2ShapeIndex empty = makeIndexOrDie("# #");
    S2ShapeIndex full = makeIndexOrDie("# # full");
    assertTrue(S2BooleanOperation.contains(empty, empty));
    assertFalse(S2BooleanOperation.contains(empty, full));
    assertTrue(S2BooleanOperation.contains(full, empty));
    assertTrue(S2BooleanOperation.contains(full, full));
  }

  /** Tests intersects() on empty and full geometries. */
  @Test
  public void testIntersectsEmptyAndFull() {
    S2ShapeIndex empty = makeIndexOrDie("# #");
    S2ShapeIndex full = makeIndexOrDie("# # full");
    assertFalse(S2BooleanOperation.intersects(empty, empty));
    assertFalse(S2BooleanOperation.intersects(empty, full));
    assertFalse(S2BooleanOperation.intersects(full, empty));
    assertTrue(S2BooleanOperation.intersects(full, full));
  }

  /**
   * Returns S2BooleanOperation.Builder that round intersection points to a fixed precision in
   * degrees, (e.g., 2 decimals). This is used to compensate for differences between intersections
   * in the "expected" data below, which were computed in lat-lng space (i.e., the rectangular
   * projection), while the actual intersections are computed using geodesics.
   */
  private static S2BooleanOperation.Builder roundToE(int exp) {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setSnapFunction(new IntLatLngSnapFunction(exp));
    return builder;
  }

  /** Tests whether the two S2ShapeIndexes are equal according to S2BooleanOperation.equals(). */
  private static boolean testEqual(String aStr, String bStr) {
    S2ShapeIndex a = makeIndexOrDie(aStr);
    S2ShapeIndex b = makeIndexOrDie(bStr);
    return S2BooleanOperation.equals(a, b);
  }

  /**
   * Runs a boolean operation of type 'opType' with the given 'builder' on two input indexes 'a' and
   * 'b'. Checks that the resulting index matches 'expected'.
   */
  private void expectResult(
      S2BooleanOperation.OpType opType,
      S2BooleanOperation.Builder builder,
      S2ShapeIndex a,
      S2ShapeIndex b,
      S2ShapeIndex expected) {
    expectResult(opType, builder, a, b, expected, "");
  }

  /**
   * Runs a boolean operation of type 'opType' with the given 'builder' on two input indexes 'a' and
   * 'b'. Checks that the resulting index matches 'expected'. On error, includes the given
   * 'detailMessage' in the report.
   */
  private void expectResult(
      S2BooleanOperation.OpType opType,
      S2BooleanOperation.Builder builder,
      S2ShapeIndex a,
      S2ShapeIndex b,
      S2ShapeIndex expected,
      String detailMessage) {
    boolean verbose = false;

    // For boolean output, the S2Builder and these Layers are never built.
    ArrayList<S2BuilderLayer> layers = new ArrayList<>();
    for (int dim = 0; dim < 3; ++dim) {
      // Since all S2Builder polygon layers require DISCARD or DISCARD_EXCESS for degenerate edges,
      // we intentionally do not require any specific multiplicity for degenerate edges and sibling
      // pairs of dimension 2.
      GraphOptions graphOptions =
          new GraphOptions(
              EdgeType.DIRECTED,
              (dim == 2) ? DegenerateEdges.DISCARD_EXCESS : DegenerateEdges.KEEP,
              DuplicateEdges.KEEP,
              (dim == 2) ? SiblingPairs.DISCARD_EXCESS : SiblingPairs.KEEP);
      layers.add(new IndexMatchingLayer(graphOptions, expected, dim));
    }

    if (verbose) {
      System.err.println("ExpectResult for OpType " + opType + " with 3 expected layers: ");
      for (S2BuilderLayer layer : layers) {
        System.err.println("  " + layer);
      }
    }

    // This boolean operation should build the output layers: it is not boolean output.
    S2BooleanOperation op = builder.build(opType, layers.get(0), layers.get(1), layers.get(2));
    S2Error error = new S2Error();
    boolean result = op.build(a, b, error);

    assertTrue(
        opType
            + " failed:\nExpected result: '"
            + S2TextFormat.toString(expected)
            + "'\nS2Error: "
            + error
            + "\nDetail: "
            + detailMessage,
        result);

    if (verbose) {
      System.err.println("ExpectResult for OpType " + opType + " with 3 output layers passed.");
    }

    // Now try the same thing with boolean output.
    boolean expectEmpty = expected.getShapes().isEmpty();
    assertEquals(
        opType
            + " failed for boolean output:\nExpected result: '"
            + S2TextFormat.toString(expected)
            + "' ("
            + (expectEmpty ? "empty)\nS2Error: " : "not empty)\nS2Error: ")
            + error
            + "\nDetail: "
            + detailMessage,
        expectEmpty,
        S2BooleanOperation.isEmpty(opType, a, b, builder.options()));
    if (verbose) {
      System.err.println("ExpectResult for OpType " + opType + " with boolean output passed.");
    }
  }

  /**
   * Runs a boolean operation of type 'opType' with the given 'builder' on two input indexes defined
   * by the provided 'aStr' and 'bStr'. Checks that the resulting index matches 'expectedStr', where
   * all strings are in the format required by makeIndexOrDie().
   */
  private void expectResult(
      S2BooleanOperation.OpType opType,
      S2BooleanOperation.Builder builder,
      String aStr,
      String bStr,
      String expectedStr) {
    S2ShapeIndex a = makeIndexOrDie(aStr);
    S2ShapeIndex b = makeIndexOrDie(bStr);
    S2ShapeIndex expected = makeIndexOrDie(expectedStr);
    expectResult(opType, builder, a, b, expected);
  }

  // This used to crash, as described in comment #3 on b/407842119
  @Test
  public void testDissolveRegression() {
    // A 13 edge single line scribble across the polygon loops.
    S2ShapeIndex pointsAndLinesIndex =
        S2TextFormat.makeIndexOrDie(
            "# "
                + "  -39.0329729033509:135.029055987371, -39.006555061068:135.025411886111,"
                + " -39.0335186080345:135.008328395693, -39.0010516583349:135.033361416367,"
                + " -39.040013198064:135.034062216659, -39.0347377044487:135.048781590026,"
                + " -39.0387048770585:135.037934012991, -39.0130769779637:135.048055978205,"
                + " -39.0057599459958:134.998902372611, -39.0401417638199:135.043644033191,"
                + " -39.0475824239056:135.037502945469, -39.0344501374733:135.018973438999,"
                + " -38.9976152401773:135.006639145402, -39.0406103224914:135.025118895254"
                + " #");
    // A single polygon with three nested loops of three edges each.
    S2ShapeIndex polygonIndex =
        S2TextFormat.makeIndexOrDie(
            "# # -39.0093951082975:134.978394915587, -39.059653460008:135.013108406395,"
                + " -39.0111604245516:135.05176055349; -39.0131085425476:135.047094046527,"
                + " -39.0555392996262:135.013272584011, -39.0115638458088:134.982897245322;"
                + " -39.0137324104193:134.987399851501, -39.051425138717:135.01343674251,"
                + " -39.0150564747747:135.042427282141");

    // The resulting polygon layer will contain at most one polygon.
    S2BooleanOperation.Builder opBuilder = new S2BooleanOperation.Builder();
    opBuilder.setSnapFunction(new IdentitySnapFunction(S1Angle.ZERO));
    opBuilder.setPolygonModel(S2BooleanOperation.PolygonModel.CLOSED);
    opBuilder.setPolylineModel(S2BooleanOperation.PolylineModel.CLOSED);
    opBuilder.setPolylineLoopsHaveBoundaries(false);
    opBuilder.setSplitAllCrossingPolylineEdges(true);

    S2PolylineVectorLayer.Options polylineOptions =
        new S2PolylineVectorLayer.Options()
            .setEdgeType(S2Builder.EdgeType.UNDIRECTED)
            .setPolylineType(S2BuilderGraph.PolylineType.WALK)
            .setDuplicateEdges(S2Builder.GraphOptions.DuplicateEdges.MERGE);
    S2ShapeIndex outputIndex = new S2ShapeIndex();

    S2ClosedSetNormalizer normalizer =
        new S2ClosedSetNormalizer(
            new IndexedLayer<>(outputIndex, new S2PointVectorLayer()),
            new IndexedLayer<>(outputIndex, new S2PolylineVectorLayer(polylineOptions)),
            new IndexedLayer<>(outputIndex, new S2PolygonLayer()));

    S2BooleanOperation op =
        opBuilder.build(
            S2BooleanOperation.OpType.UNION,
            normalizer.pointLayer(),
            normalizer.lineLayer(),
            normalizer.polygonLayer());

    S2Error error = new S2Error();
    boolean ok = op.build(pointsAndLinesIndex, polygonIndex, error);
    if (!ok || !error.ok()) {
      Assert.fail("S2BooleanOperation failed: " + error);
    }
    System.err.println("Dissolved result\n" + S2TextFormat.toString(outputIndex));
  }
}
