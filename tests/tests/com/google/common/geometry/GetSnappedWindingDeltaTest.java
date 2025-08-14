/*
 * Copyright 2025 Google Inc.
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

import static com.google.common.geometry.GetSnappedWindingDelta.findFirstVertexId;
import static com.google.common.geometry.GetSnappedWindingDelta.getSnappedWindingDelta;
import static java.lang.Math.PI;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.GetSnappedWindingDelta.InputEdgeFilter;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2Shape.MutableEdge;
import it.unimi.dsi.fastutil.ints.Int2IntAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class GetSnappedWindingDeltaTest {

  /**
   * This S2Builder layer simply calls getSnappedWindingDelta() with the given "refInputEdgeId" and
   * compares the result to "expectedWindingDelta".
   */
  private static class WindingNumberComparingLayer implements S2BuilderLayer {

    private final int refInputEdgeId;
    private final S2Builder builder;
    private final int expectedWindingDelta;
    private final InputEdgeFilter inputEdgeFilter;

    /**
     * Constructs a WindingNumberComparingLayer to compare the expectedWindingDelta to the result of
     * getSnappedWindingDelta() on the given refInputEdgeId. Applies the given input edge filter.
     */
    public WindingNumberComparingLayer(
        int refInputEdgeId,
        S2Builder builder,
        int expectedWindingDelta,
        InputEdgeFilter inputEdgeFilter) {
      this.refInputEdgeId = refInputEdgeId;
      this.builder = builder;
      this.expectedWindingDelta = expectedWindingDelta;
      this.inputEdgeFilter = inputEdgeFilter;
    }

    @Override
    public GraphOptions graphOptions() {
      return new GraphOptions(
          EdgeType.DIRECTED, DegenerateEdges.KEEP, DuplicateEdges.MERGE, SiblingPairs.KEEP);
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      // Find the reference vertex before and after snapping.
      S2Point refIn = builder.inputEdgeSrcVertex(refInputEdgeId);
      int refVertexId = GetSnappedWindingDelta.findFirstVertexId(refInputEdgeId, g);
      assertTrue(refVertexId >= 0);
      int windingDelta =
          getSnappedWindingDelta(refIn, refVertexId, inputEdgeFilter, builder, g, error);
      // If a failure already occurred, don't overwrite the error.
      if (!error.ok()) {
        return false;
      }
      if (expectedWindingDelta != windingDelta) {
        error.init(
            S2Error.Code.USER_DEFINED_START,
            "Expected winding delta of %d, but got %d",
            expectedWindingDelta, windingDelta);
      }
      return error.ok();
    }
  }

  /**
   * Given a set of loops, a set of forced vertices, and a snap radius in degrees, verifies that the
   * change in winding number computed by getSnappedWindingDelta() for the degenerate edge
   * "refInputEdgeId" is "expectedWindingDelta", and that no errors occur due to graph invalidity
   * such as edges not forming loops, or other problems.
   */
  private void expectWindingDelta(
      String loopsStr,
      String forcedVerticesStr,
      double snapRadiusDegrees,
      int refInputEdgeId,
      int expectedWindingDelta) {
    expectWindingDelta(
        loopsStr,
        forcedVerticesStr,
        snapRadiusDegrees,
        refInputEdgeId,
        expectedWindingDelta,
        InputEdgeFilter.noFilter());
  }

  /**
   * Given a set of loops, a set of forced vertices, and a snap radius in degrees, verifies that the
   * change in winding number computed by getSnappedWindingDelta() for the degenerate edge
   * "refInputEdgeId" is "expectedWindingDelta", and that no errors occur due to graph invalidity
   * such as edges not forming loops, or other problems.
   */
  private void expectWindingDelta(
      String loopsStr,
      String forcedVerticesStr,
      double snapRadiusDegrees,
      int refInputEdgeId,
      int expectedWindingDelta,
      InputEdgeFilter inputEdgeFilter) {
    S2Builder builder =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(snapRadiusDegrees))).build();
    builder.startLayer(
        new WindingNumberComparingLayer(
            refInputEdgeId, builder, expectedWindingDelta, inputEdgeFilter));
    List<S2Point> forcedVertices = S2TextFormat.parsePointsOrDie(forcedVerticesStr);
    for (S2Point v : forcedVertices) {
      builder.forceVertex(v);
    }

    S2LaxPolygonShape loops = S2TextFormat.makeLaxPolygonOrDie(loopsStr);
    builder.addShape(loops);

    MutableEdge refEdge = new MutableEdge();
    builder.getInputEdge(refInputEdgeId, refEdge);
    assertTrue(
        "Reference edge must be degenerate", refEdge.getStart().equalsPoint(refEdge.getEnd()));

    S2Error error = new S2Error();
    boolean built = builder.build(error);
    assertTrue(error.text(), built);
  }

  // Since the GetSnappedWindingDelta algorithm doesn't depend on absolute vertex positions, we
  // use "0:0" as the snapped vertex position in the tests below (using forceVertex() to ensure
  // this). Most tests use a snap radius of 10 degrees since this makes it convenient to construct
  // suitable tests.
  //
  // FYI the S2.refDir() direction for "0:0" (used to determine whether a loop contains one of its
  // vertices) is just slightly north of due west, i.e. approximately 179.7 degrees CCW from the
  // positive longitude axis.

  // No edges except the degenerate edges that defines the reference vertex.
  @Test
  public void testNoOtherEdges() {
    expectWindingDelta("0:0", "0:0", 10.0, 0, 0);
    // Also, check that an incorrect winding delta does indeed cause a test failure.
    assertThrows(AssertionError.class, () -> expectWindingDelta("0:0", "0:0", 10.0, 0, 1));
  }

  // Degenerate input loops.
  @Test
  public void testDegenerateInputLoops() {
    expectWindingDelta("0:0; 1:1; 2:2", "0:0", 10.0, 0, 0);
  }

  // Duplicate degenerate input loops.
  @Test
  public void testDuplicateDegenerateInputLoops() {
    expectWindingDelta("0:0; 0:0; 1:1; 1:1", "0:0", 10.0, 0, 0);
  }

  // A triangular shell around the reference vertex that collapses to a single point.
  @Test
  public void testCollapsingShell() {
    // The original winding number of the reference vertex is one. Snapping collapses the shell onto
    // the reference vertex, so it does not enclose the reference vertex any more, making the
    // snapped winding number zero.
    expectWindingDelta("0:0; 1:1, 1:-2, -2:1", "0:0", 10.0, 0, -1);
    // As above, but with two collapsing shells, changing the winding number from two to zero.
    expectWindingDelta(
        "0:0; 1:1, 1:-2, -2:1; 2:2, 2:-3, -3:2",
        "0:0",
        10.0,
        0,
        -2);

    // Like the previous case with two collapsing shells, but now the second one is ignored by the
    // supplied edge filter, so the change in winding number is -1 again.
    expectWindingDelta(
        "0:0; 1:1, 1:-2, -2:1; 2:2, 2:-3, -3:2",
        "0:0",
        10.0,
        0,
        -1,
        // Ignore edges from the second shell, which has edge ids 4, 5, 6.
        edgeId -> edgeId > 3);
  }

  // A hole around the reference vertex that collapses to a single point.
  @Test
  public void testCollapsingHole() {
    expectWindingDelta("0:0; 1:1, -2:1, 1:-2", "0:0", 10.0, 0, +1);
  }

  // A single "shell" that winds around the reference vertex twice.
  @Test
  public void testCollapsingDoubleShell() {
    expectWindingDelta("0:0; 1:1, 1:-2, -2:1, 2:2, 2:-3, -3:2", "0:0", 10.0, 0, -2);
  }

  // A loop that enters the Voronoi region of the snapped reference vertex and then leaves again,
  // where the reference vertex is not contained by the loop and does not move during snapping.
  @Test
  public void testExternalLoopRefVertexStaysOutside() {
    expectWindingDelta("0:0; 20:0, 0:0, 0:20", "0:0", 10.0, 0, 0);
  }

  // Like the above, except that the reference vertex is contained by the loop. (See S2.refDir()
  // comments above.)
  @Test
  public void testExternalLoopRefVertexStaysInside() {
    expectWindingDelta("0:0; 0:-20, 0:0, 20:0", "0:0", 10.0, 0, 0);
  }

  // The reference vertex moves from outside to inside an external loop during snapping.
  @Test
  public void testExternalLoopRefVertexMovesInside() {
    expectWindingDelta("1:1; 0:-20, 1:-1, 20:0", "0:0", 10.0, 0, +1);
  }

  // A single loop edge crosses the Voronoi region of the reference vertex and the reference vertex
  // stays outside the loop during snapping.
  @Test
  public void testCrossingEdgeRefVertexStaysOutside() {
    expectWindingDelta("-1:-1; 20:-20, -20:20, 20:20", "0:0", 10.0, 0, 0);
  }

  // A single loop edge crosses the Voronoi region of the reference vertex and the reference vertex
  // moves outside the loop during snapping.
  @Test
  public void testCrossingEdgeRefVertexMovesOutside() {
    expectWindingDelta("1:1; 20:-20, -20:20, 20:20", "0:0", 10.0, 0, -1);
  }

  // An external loop that winds CW around the reference vertex twice, where the reference vertex
  // moves during snapping, and where the reference vertex is outside the loop after snapping (so
  // that its winding number only increases by 1).
  @Test
  public void testExternalLoopDoubleHoleToSingleHole() {
    expectWindingDelta("4:4; 0:20, 3:3, 6:3, 2:7, 2:2, 2:20", "0:0", 10.0, 0, +1);
  }

  // An external loop that winds CW around the reference vertex twice, where the reference vertex
  // moves during snapping, and where the reference vertex is inside the loop after snapping (so
  // that its winding number increases by 3).
  @Test
  public void testExternalLoopDoubleHoleToSingleShell() {
    expectWindingDelta("4:4; 0:-20, 6:2, 2:6, 2:2, 6:2, 2:6, 2:2, 20:0", "0:0", 10.0, 0, +3);
  }

  // This and the following tests verify that the partial loops formed by the local input and output
  // edges are closed consistently with each other (such that the hypothetical connecting edges can
  // deform from one to the other without passing through the reference vertex).
  //
  // An external loop where the input edges that enter/exit the Voronoi region cross, but the
  // snapped edges do not. (This can happen when the input edges snap to multiple edges and the
  // crossing occurs outside the Voronoi region of the reference vertex.)  In this particular test,
  // the entering/exiting edges snap to the same adjacent Voronoi site so that the snapped edges
  // form a loop with one external vertex.
  @Test
  public void testExternalEdgesCrossSnapToSameVertex() {
    expectWindingDelta("1:1; -5:30, 7:-3, -7:-3, 5:30", "0:0, 0:15", 10.0, 0, -1);
  }

  // This test is similar except that the entering/exiting edges snap to two different external
  // Voronoi sites. Again, the input edges cross but the snapped edges do not.
  @Test
  public void testExternalEdgesCrossSnapToDifferentVertices() {
    expectWindingDelta("1:1; -5:40, 7:-3, -7:-3, 5:40", "0:0, 6:10, -6:10", 10.0, 0, -1);
  }

  // Test cases where the winding numbers of the reference points Za and Zb in the algorithm
  // description change due to snapping. (The points Za and Zb are the centers of the great circles
  // defined by the first and last input edges.)  For their winding numbers to change, the input
  // loop needs to cross these points as it deforms during snapping.
  //
  // In all the test below we use perturbations of 70:-180, 5:0 as the first input edge (which
  // yields points close to 0:90 as Za) and perturbations of 0:5, 0:110 as the last input edge
  // (which yields points close to 90:0 as Zb). The first/last vertex can be adjusted slightly to
  // control which side of each edge Za/Zb is on.
  @Test
  public void testReferencePointWindingNumbersChange() {
    // Winding number of Za ~= 0.01:90 changes.
    expectWindingDelta("1:1; 70:-179.99, 5:0, 0:5, -0.01:110", "0:0, 1:90", 10.0, 0, 0);

    // Winding number of Zb ~= 89.99:90 changes.
    expectWindingDelta("1:1; 70:-179.99, 5:0, 0:5, -0.01:110", "0:0, 89:90", 10.0, 0, 0);

    // Winding numbers of Za and Zb both change.
    expectWindingDelta("1:1; 70:-179.99, 5:0, 0:5, -0.01:110", "0:0, 1:90, 89:90", 10.0, 0, 0);

    // Winding number of Za ~= -0.01:90 changes in the opposite direction.
    expectWindingDelta("1:1; 70:179.99, 5:0, 0:5, 0:110", "0:0, -1:20, 1:90", 10.0, 0, 0);
  }

  // This test demonstrates that a connecting vertex may be necessary in order to ensure that the
  // loops L and L' used to compute the change in winding number for reference points Za and Zb do
  // not pass through those points during snapping. (This can only happen when the edges A0A1 or
  // B0B1 snap to an edge chain longer than 180 degrees, i.e. where the shortest edge between their
  // new endpoints goes the wrong way around the sphere.)
  @Test
  public void testReferenceLoopsTopologicallyConsistent() {
    // A0A1 follows the equator, Za is at the north pole. A0A1 snaps to the polyline (0:148, 0:74,
    // 44:-39, -31:-48), where the last two vertices are A0' and R' respectively. (The perpendicular
    // bisector of A0' and R' just barely intersects A0A1 and therefore it snaps to both vertices.)
    // A connecting vertex is needed between A0 = 0:148 and A0' = 44:-39 to ensure that this edge
    // stays within the snap radius of A0A1.
    expectWindingDelta(
        "-45:24; 0:148, 0:0, -31:-48, 44:-39, -59:0", "-31:-48, 44:-39", 60.0, 0, -1);

    // This tests the symmetric case where a connecting vertex is needed between B1 and B1' to
    // ensure that the edge stays within the snap radius of B0B1.
    expectWindingDelta(
        "-45:24;  -59:0, 44:-39, -31:-48, 0:0, 0:148", "-31:-48, 44:-39", 60.0, 0, 1);
  }

  // A complex example with multiple loops that combines many of the situations tested individually
  // above.
  @Test
  public void testComplexExample() {
    expectWindingDelta(
        "1:1; "
            + "70:179.99, 5:0, 0:5, 0:110; "
            + "70:179.99, 0:0, 0:3, 3:0, 0:-1, 0:110; "
            + "10:-10, -10:10, 10:10; "
            + "2:2, 1:-2, -1:2, 2:2, 1:-2, -1:2 ",
        "0:0, -1:90, 1:90, 45:-5",
        10.0,
        0,
        -5);
  }

  // This test demonstrates the necessity of the algorithm step that reverses the sign of Za, Zb if
  // necessary to point away from the Voronoi site R'. Further examples can be generated by running
  // the "randomLoops" test below for enough iterations.
  @Test
  public void testEnsureZaZbNotInVoronoiRegion() {
    expectWindingDelta("30:42, 30:42; -27:52, 66:131, 30:-93", "", 67.0, 0, -1);
  }

  // This test demonstrates the necessity of closing the "chainDiff" loop used by the algorithm.
  // Further examples can be found by running RandomLoops.
  @Test
  public void testEnsureChainDiffLoopIsClosed() {
    expectWindingDelta("8:26, 8:26; -36:70, -64:-35, -41:48", "", 66, 0, 0);
  }

  // This test previously failed due to a bug in getVoronoiSiteExclusion() involving long edges
  // (near 180 degrees) and large snap radii.
  @Test
  public void testVoronoiExclusionBug() {
    expectWindingDelta(
        "24.97:102.02, 24.97:102.02; " + "25.84:131.46, -29.23:-166.58, 29.40:173.03, -18.02:-5.83",
        "",
        64.83,
        0,
        -1);
  }

  /** Used to build a histogram of winding numbers. */
  private static class WindingTally extends Int2IntAVLTreeMap {
    // absl::btree_map<int, int> in C++
  }

  /**
   * This S2BuilderLayer checks that the change in winding number due to snapping computed by
   * getSnappedWindingDelta() is correct for the given configuration of input edges.
   *
   * <p>"refInputEdgeId" should be a degenerate edge "SS" that specifies the reference vertex R
   * whose change in winding number is verified.
   *
   * <p>"isolatedInputEdgeId" should be a degenerate edge "II" that is not expected to snap together
   * with any other edges. (This ensures that the winding number of its vertex I does not change due
   * to snapping.) "I" should be chosen to be as far away as possible from other vertices and edges
   * used for testing purposes. If more than one edge happens to snap to this vertex,
   * S2Error.Code.FAILED_PRECONDITION is returned.
   */
  private static class WindingNumberCheckingLayer implements S2BuilderLayer {
    private final TestDataGenerator data;
    private final int refInputEdgeId;
    private final int isolatedInputEdgeId;
    private final S2Builder builder;
    private final WindingTally windingTally;

    /** Constructs a WindingNumberCheckingLayer with the given parameters. */
    public WindingNumberCheckingLayer(
        TestDataGenerator data,
        int refInputEdgeId,
        int isolatedInputEdgeId,
        S2Builder builder,
        WindingTally windingTally) {
      this.data = data;
      this.refInputEdgeId = refInputEdgeId;
      this.isolatedInputEdgeId = isolatedInputEdgeId;
      this.builder = builder;
      this.windingTally = windingTally;
    }

    @Override
    public GraphOptions graphOptions() {
      // Some of the graph options are chosen randomly.
      return new GraphOptions(
          EdgeType.DIRECTED,
          DegenerateEdges.KEEP,
          data.nextBoolean() ? DuplicateEdges.KEEP : DuplicateEdges.MERGE,
          data.nextBoolean() ? SiblingPairs.KEEP : SiblingPairs.CREATE);
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      // First we locate the vertices R, R', I, I'.
      S2Point refIn = builder.inputEdgeSrcVertex(refInputEdgeId);
      int refVertexId = findFirstVertexId(refInputEdgeId, g);
      S2Point refOut = g.vertex(refVertexId);

      S2Point isoIn = builder.inputEdgeSrcVertex(isolatedInputEdgeId);
      int isoVertexId = findFirstVertexId(isolatedInputEdgeId, g);
      S2Point isoOut = g.vertex(isoVertexId);

      // If more than one edge snapped to the isolated vertex I, we skip this test since we have no
      // way to independently verify correctness of the results.
      S2BuilderGraph.VertexOutMap outMap = new S2BuilderGraph.VertexOutMap(g);
      if (outMap.outDegree(isoVertexId) != 1
          || g.inputEdgeIds(outMap.edgeIds(isoVertexId).begin()).size() != 1) {
        error.init(S2Error.Code.FAILED_PRECONDITION, "Isolated vertex not isolated");
        return false;
      }

      // Next we compute the winding number of R relative to I by counting signed edge crossings of
      // the input edges, and the winding number of R' related to I' by counting signed edge
      // crossings of the output edges. (In order to support DuplicateEdges.MERGE and
      // SiblingEdges.CREATE, we also need to take into account the number of input edges that
      // snapped to each output edge.)
      int windingIn = 0;
      EdgeCrosser crosser = new EdgeCrosser(isoIn, refIn);
      for (int e = 0; e < builder.numInputEdges(); ++e) {
        windingIn +=
            crosser.signedEdgeOrVertexCrossing(
                builder.inputEdgeSrcVertex(e), builder.inputEdgeDstVertex(e));
      }

      int windingOut = 0;
      crosser.init(isoOut, refOut);
      for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
        windingOut +=
            g.inputEdgeIds(edgeId).size()
                * crosser.signedEdgeOrVertexCrossing(g.edgeSrc(edgeId), g.edgeDst(edgeId));
      }

      // Finally, check that getSnappedWindingDelta() computes the difference between the two (using
      // only local snapping information).
      int windingDelta =
          getSnappedWindingDelta(refIn, refVertexId, InputEdgeFilter.noFilter(), builder, g, error);
      assertTrue(error.text(), error.ok());

      if (windingOut - windingIn != windingDelta) {
        error.init(
            S2Error.Code.USER_DEFINED_START,
            "Expected windingDelta(%d) to equal windingOut(%d) - windingIn(%d) but it was not.",
            windingDelta,
            windingOut,
            windingIn);
        return false;
      }

      windingTally.compute(windingDelta, (k, v) -> v == null ? 1 : v + 1);
      return true;
    }
  }

  @Test
  public void testGetSnappedWindingDeltaRandomLoops() {
    // Count the number of tests for each winding number result and also the number of tests where
    // the isolated vertex was not isolated, to verify that the test is working as intended.

    // This test passes with at least 10 million iterations, but that takes too long to run as
    // part of normal presubmits, etc. so only 1000 iterations are run normally.
    int numIters = 1000;
    int numNotIsolated = 0;
    WindingTally windingTally = new WindingTally();
    TestDataGenerator data = new TestDataGenerator();

    for (int iter = 0; iter < numIters; ++iter) {
      StringBuilder sb = new StringBuilder();
      sb.append("Iteration ").append(iter);

      // Choose a random snap radius up to the allowable maximum.
      S1Angle snapRadius = S2BuilderSnapFunctions.maxSnapRadius().mul(data.uniform(0.0, 1.0));
      sb.append("\n  snapRadius = ").append(snapRadius.degrees()).append(" degrees");
      S2Builder builder = new S2Builder.Builder(new IdentitySnapFunction(snapRadius)).build();
      builder.startLayer(
          new WindingNumberCheckingLayer(
              data, /* refInputEdgeId= */ 0, /* isolatedInputEdgeId= */ 1, builder, windingTally));

      // Choose a random reference vertex, and an isolated vertex that is as far away as possible.
      // (The small amount of perturbation reduces the number of calls to
      // S2Predicates.Sign.expensive() and is not necessary for correctness.)
      S2Point ref = data.getRandomPoint();
      sb.append("\n  ref = ").append(ref.toDegreesString());
      S2Point isolated = ref.neg().add(S2.ortho(ref).mul(1e-12)).normalize();
      sb.append("\n  isolated = ").append(isolated.toDegreesString());
      builder.addEdge(ref, ref); // Reference vertex edge.
      builder.addEdge(isolated, isolated); // Isolated vertex edge.

      // Now repeatedly build and add loops. Loops consist of 1 or more random vertices where
      // approximately 2/3 are expected to be within snapRadius of the reference vertex. Some
      // vertices are duplicates of previous vertices. Vertices are kept at least snapRadius away
      // from the isolated vertex to reduce the chance that edges will snap to it. (This can still
      // happen with long edges, or because the reference vertex snapped to a new location far away
      // from its original location.)
      List<S2Point> verticesUsed = new ArrayList<>();
      List<S2Point> loop = new ArrayList<>();
      for (int numLoops = data.uniformInt(1, 6); --numLoops >= 0; ) {
        for (int numVertices = data.uniformInt(1, 10); --numVertices >= 0; ) {
          if (!verticesUsed.isEmpty() && data.oneIn(4)) {
            loop.add(verticesUsed.get(data.uniformInt(0, verticesUsed.size())));
          } else if (data.oneIn(3)) {
            S2Point v =
                data.samplePoint(S2Cap.fromAxisAngle(ref, S1Angle.radians(PI).sub(snapRadius)));
            loop.add(v);
            verticesUsed.add(v);
          } else {
            S2Point v = data.samplePoint(S2Cap.fromAxisAngle(ref, snapRadius));
            loop.add(v);
            verticesUsed.add(v);
          }
        }
        S2Loop shape = new S2Loop(loop);
        builder.addShape(shape);
        sb.append("\n  loop ").append(S2TextFormat.toString(shape));
        loop.clear();
      }
      S2Error error = new S2Error();
      if (!builder.build(error)) {
        if (error.code() == S2Error.Code.INVALID_ARGUMENT) {
          throw new AssertionError(error.text() + "\n" + sb);
        }
        if (error.code() == S2Error.Code.USER_DEFINED_START) {
          throw new AssertionError("FAILURE: " + error.text() + "\n" + sb);
        }
        ++numNotIsolated;
      }
    }

    // We expect that at most 20% of tests will result in an isolated vertex.
    assertTrue(numNotIsolated <= 0.2 * numIters);
    System.err.println("Histogram of winding number deltas tested:");
    for (Int2IntMap.Entry entry : windingTally.int2IntEntrySet()) {
      System.err.println(entry.getIntKey() + " : " + entry.getIntValue() + "\n");
    }
  }

  /**
   * Given a set of edges, a set of forced vertices, and a snap radius in degrees, verifies that the
   * change in winding number computed by getSnappedWindingDelta() for the degenerate edge
   * "refInputEdgeId" is "expectedWindingDelta", if it gets that far. Returns an S2Error indicating
   * if the test passed or failed.
   */
  public S2Error expectWindingDeltaErrorForEdges(
      String edgesStr,
      String forcedVerticesStr,
      double snapRadiusDegrees,
      int refInputEdgeId,
      int expectedWindingDelta) {
    S2Builder builder =
        new S2Builder.Builder(new IdentitySnapFunction(S1Angle.degrees(snapRadiusDegrees))).build();
    builder.startLayer(
        new WindingNumberComparingLayer(
            refInputEdgeId, builder, expectedWindingDelta, InputEdgeFilter.noFilter()));
    List<S2Point> forcedVertices = S2TextFormat.parsePointsOrDie(forcedVerticesStr);
    for (S2Point v : forcedVertices) {
      builder.forceVertex(v);
    }

    List<S2Edge> edges = S2TextFormat.makeEdgesOrDie(edgesStr);
    for (S2Edge edge : edges) {
      builder.addEdge(edge.getStart(), edge.getEnd());
    }

    MutableEdge refEdge = new MutableEdge();
    builder.getInputEdge(refInputEdgeId, refEdge);
    assertTrue(
        "Reference edge must be degenerate", refEdge.getStart().equalsPoint(refEdge.getEnd()));

    S2Error error = new S2Error();
    boolean built = builder.build(error);
    // If build fails, the error must be set.
    assertEquals(error.ok(), built);
    return error;
  }

  @Test
  public void testEdgesDoNotFormLoops() {
    // Besides the reference vertex, a single edge where both ends snap to the reference vertex.
    S2Error error = expectWindingDeltaErrorForEdges("0:0, 0:0; 1:1, 2:2;", "0:0", 10.0, 0, 0);
    assertEquals(S2Error.Code.INVALID_ARGUMENT, error.code());
    assertTrue("Actual text is '" + error.text() + "'", error.text().contains("do not form loops"));
  }
}
