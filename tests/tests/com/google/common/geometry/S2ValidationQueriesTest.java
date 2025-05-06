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

import static com.google.common.geometry.S2ValidationQueries.sortEdgesCcw;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2ValidationQueries.S2LegacyValidQuery;
import com.google.common.geometry.S2ValidationQueries.S2ValidQuery;
import com.google.common.geometry.ValidationQueriesTestUtils.ValidationQueryTest;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit Tests for S2ValidationQueries. */
@RunWith(JUnit4.class)
public class S2ValidationQueriesTest extends GeometryTestCase {
  private static final S2Point ORIGIN = S2LatLng.fromRadians(0, 0).toPoint();

  /** Tests all (both) validation queries defined by S2ValidationQueries. */
  private final ValidationQueryTest avq =
      new ValidationQueryTest(new S2ValidQuery(), new S2LegacyValidQuery());

  /** Tests validation queries that allow multidimensional geometry. */
  private final ValidationQueryTest mdq = new ValidationQueryTest(new S2ValidQuery());

  /** Tests S2ValidQuery only. */
  private final ValidationQueryTest s2vq = new ValidationQueryTest(new S2ValidQuery());

  /** Tests S2LegacyValidQuery only. */
  private final ValidationQueryTest s2lvq = new ValidationQueryTest(new S2LegacyValidQuery());

  @Test
  public void testSortEdgesCcwSSortsEdges() {
    int kNumEdges = 10;

    // Generate edges in order around the origin.
    List<S2Edge> sorted = TestDataGenerator.ccwEdgesAbout(ORIGIN, kNumEdges);

    for (int i = 0; i < kNumEdges; ++i) {
      // Shift sorted edges left by one.
      rotate(sorted, 1);

      // Copy the list of edges, then randomize the copy.
      List<S2Edge> shuffled = new ArrayList<>(sorted);
      shuffle(shuffled);

      // Sorting about the first edge should always give back the sorted array.
      sortEdgesCcw(ORIGIN, sorted.get(0), shuffled);
      assertEquals(shuffled, sorted);
    }
  }

  S2Edge reverseOf(S2Edge edge) {
    return new S2Edge(edge.getEnd(), edge.getStart());
  }

  @Test
  public void testSortEdgesCcwSortsEdgesFlipped() {
    int kNumEdges = 10;

    // Generate edges in order around the origin.
    List<S2Edge> sorted = TestDataGenerator.ccwEdgesAbout(ORIGIN, kNumEdges);

    // Flip the orientation of some of the edges. This only changes their direction, not their
    // ordering around the vertex.
    sorted.set(3, reverseOf(sorted.get(3)));
    sorted.set(8, reverseOf(sorted.get(8)));

    for (int i = 0; i < kNumEdges; ++i) {
      // Shift sorted edges left by one.
      rotate(sorted, 1);

      // Randomize the edges.
      List<S2Edge> shuffled = new ArrayList<>(sorted);
      shuffle(shuffled);

      // Sorting about the first edge should always give back the sorted array.
      sortEdgesCcw(ORIGIN, sorted.get(0), shuffled);
      assertEquals(shuffled, sorted);
    }
  }

  @Test
  public void testSortEdgesCcwStartEdgeAlwaysFirst() {
    int kNumEdges = 10;

    // Generate edges in order around the origin.
    List<S2Edge> sorted = TestDataGenerator.ccwEdgesAbout(ORIGIN, kNumEdges);

    for (int i = 0; i < kNumEdges; ++i) {
      // Randomize the edges.
      List<S2Edge> shuffled = new ArrayList<>(sorted);
      shuffle(shuffled);

      // The reference edge we sort around should always come out first.
      sortEdgesCcw(ORIGIN, sorted.get(i), shuffled);
      // Reference equality is actually what we want here.
      assertSame(shuffled.get(0), sorted.get(i));
    }
  }

  @Test
  public void testSortEdgesCcwReverseDuplicatesOrdered() {
    int kNumEdges = 10;

    // Generate edges in order around the origin.
    List<S2Edge> sorted = TestDataGenerator.ccwEdgesAbout(ORIGIN, kNumEdges);

    // Insert two reverse duplicate edges.
    sorted.add(8, reverseOf(sorted.get(8)));
    sorted.add(3, reverseOf(sorted.get(3)));

    // Randomize the edges.
    List<S2Edge> shuffled = new ArrayList<>(sorted);
    shuffle(shuffled);

    // After sorting, the reverse duplicates should be right after their sibling.
    sortEdgesCcw(ORIGIN, sorted.get(4), shuffled);

    S2Point common = sorted.get(4).getStart();
    assertEquals(shuffled.get(0), reverseOf(shuffled.get(1)));
    assertEquals(shuffled.get(0).getStart(), common);
    assertEquals(shuffled.get(6), reverseOf(shuffled.get(7)));
    assertEquals(shuffled.get(6).getStart(), common);
  }

  //////////////////   All Validation Query Tests   ////////////////////

  @Test
  public void testAllValidationQueriesBasicGeometryOk() {
    ValidationQueriesTestUtils.expectAllValidationQueriesBasicGeometryOk(avq);
  }

  @Test
  public void testAllValidationQueriesEmptyGeometryOk() {
    ValidationQueriesTestUtils.expectAllValidationQueriesEmptyGeometryOk(avq);
  }

  @Test
  public void testAllValidationQueriesFullGeometryOk() {
    ValidationQueriesTestUtils.expectAllValidationQueriesFullGeometryOk(avq);
  }

  @Test
  public void testAllValidationQueriesInteriorOnRightRegression() {
    ValidationQueriesTestUtils.expectAllValidationQueriesInteriorOnRightRegression(avq);
  }

  @Test
  public void testAllValidationQueriesTangentPolygonsOk() {
    ValidationQueriesTestUtils.expectAllValidationQueriesTangentPolygonsOk(avq);
  }

  @Test
  public void testAllValidationQueriesAntipodalEdgeFails() {
    ValidationQueriesTestUtils.expectAllValidationQueriesAntipodalEdgeFails(avq);
  }

  @Test
  public void testAllValidationQueriesBadlyDimensionedFails() {
    ValidationQueriesTestUtils.expectAllValidationQueriesBadlyDimensionedFails(avq);
  }

  @Test
  public void testAllValidationQueriesOpenChainFails() {
    ValidationQueriesTestUtils.expectAllValidationQueriesOpenChainFails(avq);
  }

  @Test
  public void testAllValidationQueriesDuplicatePolygonEdgesFail() {
    ValidationQueriesTestUtils.expectAllValidationQueriesDuplicatePolygonEdgesFail(avq);
  }

  @Test
  public void testAllValidationQueriesChainsTouchingOk() {
    ValidationQueriesTestUtils.expectAllValidationQueriesChainsTouchingOk(avq);
  }

  @Test
  public void testAllValidationQueriesNestedShellsFail() {
    ValidationQueriesTestUtils.expectAllValidationQueriesNestedShellsFail(avq);
  }

  @Test
  public void testAllValidationQueriesChainsCannotCross() {
    ValidationQueriesTestUtils.expectAllValidationQueriesChainsCannotCross(avq);
  }

  @Test
  public void testAllValidationQueriesShellInHoleFails() {
    ValidationQueriesTestUtils.expectAllValidationQueriesShellInHoleFails(avq);
  }

  @Test
  public void testAllValidationQueriesLoopsCrossing() {
    ValidationQueriesTestUtils.expectAllValidationQueriesLoopsCrossing(data, avq);
  }

  //////////////////   Multidimensional Query Tests   ////////////////////

  // Basic multi-dimensional geometry should be valid.
  @Test
  public void testMultiDimensionalQueriesBasicGeometryOk() {
    ValidationQueriesTestUtils.expectMultiDimensionalQueriesBasicGeometryOk(mdq);
  }

  @Test
  public void testMultiDimensionalQueriesContainedGeometryFails() {
    ValidationQueriesTestUtils.expectMultiDimensionalQueriesContainedGeometryFails(mdq);
  }

  //////////////////   S2ValidQuery Tests   ////////////////////

  @Test
  public void testS2ValidTestQuiltIsValid() {
    s2vq.expectShapeValid(ValidationQueriesTestUtils.makeQuilt());
  }

  @Test
  public void testS2ValidTestDegenerateRingsAllowed() {
    // Point loops with a single edge should be allowed.
    s2vq.expectGeometryValid("## 0:0");

    // Degenerate loops with a sibling edge pair should be allowed.
    s2vq.expectGeometryValid("## 0:0, 1:1");
  }

  @Test
  public void testS2ValidTestSplitInteriorsOk() {
    // Chains touching at two points (splitting interior), should be fine.
    s2vq.expectGeometryValid("## 3:0, 0:-3, -3:0, 0:+3;" + "   3:0, 0:+1, -3:0, 0:-1");
  }

  @Test
  public void testS2ValidTestPolylineEdgesCrossSemanticsOk() {
    // Interior crossings between polylines are ok.
    s2vq.expectGeometryValid("# 0:0, 1:1, 0:2, 1:3, 0:4 " + "| 1:0, 0:1, 1:2, 0:3, 1:4 #");

    // Vertex crossings between polylines are ok.
    s2vq.expectGeometryValid(
        "# 0:0, 1:1, 2:2, 1:3, 0:4, 1:5, 2:6, 1:7, 0:8"
            + "| 2:0, 1:1, 0:2, 1:3, 2:4, 1:5, 0:6, 1:7, 2:8 #");

    // Interior crossings within a polyline are ok.
    s2vq.expectGeometryValid("# 0:0, 1:1, 0:2, 1:3, 0:4, 1:4, 0:3, 1:2, 0:1, 1:0 #");

    // Vertex crossings within a polyline are ok.
    s2vq.expectGeometryValid(
        "# 0:0, 1:1, 2:2, 1:3, 0:4, 1:5, 2:6, 1:7, 0:8,"
            + "  2:0, 1:1, 0:2, 1:3, 2:4, 1:5, 0:6, 1:7, 2:8 #");

    // A closed loop that touches at the endpoints is ok.
    s2vq.expectGeometryValid("# 2:1, 1:0, 0:1, 1:2, 2:1 #");

    // Two polylines touching at endpoints is also OK;
    s2vq.expectGeometryValid("# 0:0, 1:1, 0:2" + "| 1:3, 0:4, 1:5 #");
  }

  @Test
  public void testS2ValidTestReverseDuplicateOnCenterWorks() {
    // A pair of reverse duplicate edges touching the cell center should work.
    s2vq.expectGeometryValid("## 2:0, 0:-2, -2:0, 0:2;" + "   0:0, 1:1");
  }

  @Test
  public void testS2ValidTestPolygonOnCentersWorks() {
    // Draw a nested-diamond polygon using 8 cells straddling the equator/prime
    // meridian. The chain orientation check will have to fall back to a brute
    // force check but these should still work.
    ImmutableList<ImmutableList<S2Point>> kPoints =
        ImmutableList.of(
            ImmutableList.of(
                new S2Cell(S2CellId.fromToken("0ec")).getCenter(),
                new S2Cell(S2CellId.fromToken("044")).getCenter(),
                new S2Cell(S2CellId.fromToken("1bc")).getCenter(),
                new S2Cell(S2CellId.fromToken("114")).getCenter()),
            ImmutableList.of(
                new S2Cell(S2CellId.fromToken("104")).getCenter(),
                new S2Cell(S2CellId.fromToken("1ac")).getCenter(),
                new S2Cell(S2CellId.fromToken("054")).getCenter(),
                new S2Cell(S2CellId.fromToken("0fc")).getCenter()));

    s2vq.expectShapeValid(S2LaxPolygonShape.create(kPoints));
  }

  @Test
  public void testS2ValidTestDegeneratePolygonOnCentersworks() {
    // Now to be really rude and test a polygon consisting of nothing but
    // reverse-duplicate pairs between cell centers.
    ImmutableList<S2Point> kPoints =
        ImmutableList.of(
            new S2Cell(S2CellId.fromToken("0ec")).getCenter(),
            new S2Cell(S2CellId.fromToken("044")).getCenter(),
            new S2Cell(S2CellId.fromToken("1bc")).getCenter(),
            new S2Cell(S2CellId.fromToken("114")).getCenter(),
            new S2Cell(S2CellId.fromToken("1bc")).getCenter(),
            new S2Cell(S2CellId.fromToken("044")).getCenter());

    s2vq.expectShapeValid(S2LaxPolygonShape.fromLoop(kPoints));

    // Same thing but the polygon just goes out at a diagonal and back.
    String[] tokens = {"1004", "1014", "1044", "1054", "1104", "1114"};

    List<List<S2Point>> loops = new ArrayList<>();
    List<S2Point> loop = new ArrayList<>();
    for (int i = 0; i < 6; ++i) {
      loop.add(new S2Cell(S2CellId.fromToken(tokens[i])).getCenter());
    }

    for (int i = 4; i > 0; --i) {
      loop.add(new S2Cell(S2CellId.fromToken(tokens[i])).getCenter());
    }
    loops.add(loop);

    s2vq.expectShapeValid(S2LaxPolygonShape.create(loops));
  }

  //////////////////   S2LegacyValid Tests   ////////////////////

  @Test
  public void testS2LegacyValidTestQuiltIsNotValid() {
    // The quilt has reverse duplicate edges near the poles.
    s2lvq.expectShapeInvalid(
        ValidationQueriesTestUtils.makeQuilt(), S2Error.Code.OVERLAPPING_GEOMETRY);
  }

  @Test
  public void testS2LegacyValidTestMultiDimensionalFails() {
    // Multi-dimensional geometry shouldn't be allowed.
    s2lvq.expectGeometryInvalid(
        "  3:0| 0:-3| -3:0| 0:3" + "# 2:0, 0:-2, -2:0, 0:2" + "# 1:0, 0:-1, -1:0, 0:1",
        S2Error.Code.INVALID_DIMENSION);
  }

  @Test
  public void testS2LegacyValidTestSplitInteriorsOk() {
    // Chains touching at two points (splitting interior), should be fine.
    s2lvq.expectGeometryValid("## 3:0, 0:-3, -3:0, 0:+3;" + "   3:0, 0:+1, -3:0, 0:-1");
  }

  @Test
  public void testS2LegacyValidTestSelfTouchingLoopFails() {
    // A loop that touches itself should fail.
    s2lvq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, -1:1, 0:-2, 1:1", S2Error.Code.DUPLICATE_VERTICES);
  }

  @Test
  public void testS2LegacyValidTestDegenerateEdgesFail() {
    // A degenerate polygon edge should fail.
    s2lvq.expectGeometryInvalid("## 2:0, 2:0, 0:-2, -2:0, 0:-2", S2Error.Code.DUPLICATE_VERTICES);

    // A degenerate polyline edge should fail.
    s2lvq.expectGeometryInvalid("# 0:0, 0:0, 1:1, 2:2 #", S2Error.Code.DUPLICATE_VERTICES);
  }

  @Test
  public void testS2LegacyValidTestShortChainsFail() {
    S2Error.Code kCode = S2Error.Code.LOOP_NOT_ENOUGH_VERTICES;

    // A polygon chain with one or two edges should fail.
    s2lvq.expectGeometryInvalid("## 0:0", kCode);
    s2lvq.expectGeometryInvalid("## 0:0, 1:1", kCode);
  }
}
