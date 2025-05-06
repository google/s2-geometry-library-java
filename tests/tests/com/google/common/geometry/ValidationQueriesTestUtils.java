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

import static com.google.common.geometry.S2.M_SQRT1_2;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Shape.ChainPosition;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ValidationQueries.S2ValidationQueryBase;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Utility classes and methods for testing validation queries. Used by S2ValidationQueriesTest as
 * well as unit tests for other validation queries outside S2 that extend the queries defined in
 * {@link S2ValidationQueries}.
 */
public final class ValidationQueriesTestUtils {

  //////////////////   All Validation Query Tests   ////////////////////

  /** Tests that basic geometry is valid, which should be true for all validation queries. */
  public static void expectAllValidationQueriesBasicGeometryOk(ValidationQueryTest avq) {
    // Basic polygon should be valid.
    avq.expectGeometryValid("## 1:0, 0:-1, -1:0, 0:1");

    // Basic polyline should be valid.
    avq.expectGeometryValid("# 0:0, 1:0, 0:-1, -1:0, 0:1 #");

    // Basic multipoint should be valid.
    avq.expectGeometryValid("0:0 | 1:0 | 0:-1 | -1:0 | 0:1 ##");

    // Basic polygon with a hole should be valid.
    avq.expectGeometryValid("## 2:0, 0:-2, -2:0, 0:2;" + "   0:1, -1:0, 0:-1, 1:0;");

    // Polygon with improperly oriented hole should fail.
    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 1:0, 0:-1, -1:0, 0:1;",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
  }

  /** Tests that empty geometry is valid, which should be true for all validation queries. */
  public static void expectAllValidationQueriesEmptyGeometryOk(ValidationQueryTest avq) {
    avq.expectGeometryValid("##");
  }

  /** Tests that full geometry is valid, which should be true for all validation queries. */
  public static void expectAllValidationQueriesFullGeometryOk(ValidationQueryTest avq) {
    avq.expectGeometryValid("## full");
  }

  /** Tests geometry which should be valid for all validation queries. */
  public static void expectAllValidationQueriesInteriorOnRightRegression(ValidationQueryTest avq) {
    // Found via fuzzing in C++, this should not test as invalid. The root cause was that we
    // weren't clearing the incident edge list before gathering edges, causing a duplicate edge to
    // appear and give us the wrong interior state.
    avq.expectGeometryValid("## 0:4, 3:128, 4:2, 0:0");
  }

  /** Two polygons touching at one vertex should be valid for all validation queries. */
  public static void expectAllValidationQueriesTangentPolygonsOk(ValidationQueryTest avq) {
    avq.expectGeometryValid("## 1:0, 0:-1, -1:0, 0:1" + "|  0:1, -1:2,  0:3, 1:2");
  }

  /**
   * Tests that an edge between antipodal points is invalid, which should be true for all validation
   * queries.
   */
  public static void expectAllValidationQueriesAntipodalEdgeFails(ValidationQueryTest avq) {
    // We need to specify the points instead of using makeIndexOrDie to ensure they're exactly
    // opposite sign to be antipodal.
    ImmutableList<S2Point> kPoints =
        ImmutableList.of(
            new S2Point(M_SQRT1_2, M_SQRT1_2, 0),
            new S2Point(0, 1, 0),
            new S2Point(-1, 0, 0),
            new S2Point(1, 0, 0));

    avq.expectShapeInvalid(S2LaxPolygonShape.fromLoop(kPoints), S2Error.Code.ANTIPODAL_VERTICES);
  }

  /** Tests geometry which should be invalid for all validation queries. */
  public static void expectAllValidationQueriesBadlyDimensionedFails(ValidationQueryTest avq) {
    avq.expectShapeInvalid(
        new ValidationQueriesTestUtils.BadDimensionShape(), S2Error.Code.INVALID_DIMENSION);
  }

  /** Tests that a loop that isn't closed should be invalid for all validation queries. */
  public static void expectAllValidationQueriesOpenChainFails(ValidationQueryTest avq) {
    avq.expectShapeInvalid(
        new ValidationQueriesTestUtils.OpenShape(), S2Error.Code.LOOP_NOT_ENOUGH_VERTICES);
  }

  /** Two polygons sharing an edge should be invalid for all validation queries. */
  public static void expectAllValidationQueriesDuplicatePolygonEdgesFail(ValidationQueryTest avq) {
    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2" + " | 2:0, 0:-2,  0:0", S2Error.Code.OVERLAPPING_GEOMETRY);
  }

  /**
   * Polygon with a hole with chains touching at a point should be valid for all validation queries.
   */
  public static void expectAllValidationQueriesChainsTouchingOk(ValidationQueryTest avq) {
    avq.expectGeometryValid("## 2:0, 0:-2, -2:0, 0:2;" + " 0:2, -1:0, 0:-1, 1:0;");
    avq.expectGeometryValid("## 2:0, 0:-2, -2:0, 0:2;" + " 0:1, -2:0, 0:-1, 1:0;");

    avq.expectGeometryInvalid(
        "## 2:0,  0:-2, -2:0, 0:2;" + " 1:0,  0:-2, -1:0, 0:2;",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
  }

  /** Polygon with improperly oriented hole should fail for all validation queries. */
  public static void expectAllValidationQueriesNestedShellsFail(ValidationQueryTest avq) {
    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 1:0, 0:-1, -1:0, 0:1",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

    // Even if they're touching.
    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 2:0, 0:-1, -1:0, 0:1",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 2:0, 0:-1, -2:0, 0:1",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 1:0, 0:-2, -1:0, 0:1",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 1:0, 0:-1, -2:0, 0:1",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

    avq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2;" + " 1:0, 0:-1, -1:0, 0:2",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
  }

  /** Polygon chains crossing is an error, regardless of orientation. */
  public static void expectAllValidationQueriesChainsCannotCross(ValidationQueryTest avq) {
    avq.expectGeometryInvalid(
        "## 3:0, 0:-3, -3:0, 0:+3;" + " 3:2, 0:-1, -3:2, 0:+5",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

    avq.expectGeometryInvalid(
        "## 0:3, 3:0,   0:-3, -3:0;" + " 3:2, 0:+5, -3:2,  0:-1",
        S2Error.Code.OVERLAPPING_GEOMETRY);

    // Crossing at a vertex isn't allowed either. This test vector is carefully constructed to
    // avoid triggering the interior-on-the-right check.
    avq.expectGeometryInvalid(
        "## 0:-6, -6:0, 0:6, 6:0 ;" + "  0:0, 3:0, 6:0, 6:3, 6:6, 3:6, 0:6, 0:3",
        S2Error.Code.OVERLAPPING_GEOMETRY);
  }

  public static void expectAllValidationQueriesShellInHoleFails(ValidationQueryTest avq) {
    avq.expectGeometryInvalid(
        "## 0:0, 10:10, 10:0; 5:21, 8:21, 6:23",
        S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
  }

  /** None of the validation semantics allow loops to cross. */
  public static void expectAllValidationQueriesLoopsCrossing(
      TestDataGenerator data, ValidationQueryTest avq) {
    int kIters = 100;

    for (int iter = 0; iter < kIters; ++iter) {
      data.addConcentricLoops(avq.loops, 2, 4 /*min_vertices*/);

      // Grab the two loops that were just added.
      List<S2Point> loop0 = avq.loops.get(avq.loops.size() - 2);
      List<S2Point> loop1 = avq.loops.get(avq.loops.size() - 1);

      // Both loops have the same number of vertices, and vertices at the same index position are
      // collinear with the center point, so we can create a crossing by simply exchanging two
      // vertices at the same index position.
      int n = loop0.size();
      int i = data.uniform(n);

      // Swap loop 0 index i with loop 1 index i.
      S2Point tmp = loop0.get(i);
      loop0.set(i, loop1.get(i));
      loop1.set(i, tmp);

      if (data.oneIn(2)) {
        // By copying the two vertices adjacent to the previously swapped vertex from one loop to
        // the other, we can ensure that the crossings happen at vertices rather than edges.
        loop0.set((i + 1) % n, loop1.get((i + 1) % n));
        loop0.set((i + n - 1) % n, loop1.get((i + n - 1) % n));
      }

      // Disable internal assertions to allow invalid polygons to be built.
      S2Polygon polygon = GeometryTestCase.uncheckedCreate(avq::makePolygon);

      avq.expectShapeInvalid(
          polygon.shape(),
          S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
          S2Error.Code.OVERLAPPING_GEOMETRY);
    }
  }

  //////////////////   Multidimensional Query Tests   ////////////////////

  /** Expects that basic multi-dimensional geometry is valid. */
  public static void expectMultiDimensionalQueriesBasicGeometryOk(ValidationQueryTest mdq) {
    mdq.expectGeometryValid(
        "  3:0| 0:-3| -3:0| 0:3" + "# 2:0, 0:-2, -2:0, 0:2" + "# 1:0, 0:-1, -1:0, 0:1");
  }

  /** Expects that multi-dimensional geometry with contained geometry is invalid. */
  public static void expectMultiDimensionalQueriesContainedGeometryFails(ValidationQueryTest mdq) {
    // Point contained in polygon is invalid.
    mdq.expectGeometryInvalid("0:0 ## 2:0, 0:-2, -2:0, 0:2", S2Error.Code.OVERLAPPING_GEOMETRY);

    // Polyline contained in polygon is invalid.
    mdq.expectGeometryInvalid(
        "# 0:-1, 0:1 # 2:0, 0:-2, -2:0, 0:2", S2Error.Code.OVERLAPPING_GEOMETRY);

    // Polygon contained in polygon is invalid.
    mdq.expectGeometryInvalid(
        "## 2:0, 0:-2, -2:0, 0:2" + " | 1:0, 0:-1, -1:0, 0:1", S2Error.Code.OVERLAPPING_GEOMETRY);

    // Either end of a polyline contained in a polygon is invalid.
    mdq.expectGeometryInvalid(
        "# 0:-3, 0:1 # 2:0, 0:-2, -2:0, 0:2", S2Error.Code.OVERLAPPING_GEOMETRY);
    mdq.expectGeometryInvalid(
        "# 0:-1, 0:3 # 2:0, 0:-2, -2:0, 0:2", S2Error.Code.OVERLAPPING_GEOMETRY);

    // Polyline edges crossing are fine though.
    mdq.expectGeometryValid("# 0:-1, 0:1 | 1:0, -1:0 #");
  }

  //////////////////   End of Query Tests   ////////////////////

  // Diamond shaped pattern of points around (0,0).
  public static final S2Point[] diamondPoints = {
    S2LatLng.fromDegrees(+1.0, +0.0).toPoint(),
    S2LatLng.fromDegrees(-1.0, +0.0).toPoint(),
    S2LatLng.fromDegrees(+0.0, -1.0).toPoint(),
    S2LatLng.fromDegrees(+0.0, +1.0).toPoint()
  };

  /** Verifies that one or more validation queries have expected validation results. */
  public static class ValidationQueryTest {
    public final List<S2ValidationQueryBase> queries = new ArrayList<>();
    public final List<List<S2Point>> loops = new ArrayList<>();

    /** Constructs a new ValidationQueryTest for the given S2ValidationQueries. */
    public ValidationQueryTest(S2ValidationQueryBase... queries) {
      Collections.addAll(this.queries, queries);
    }

    /** Checks that a particular shape is valid in each of the configured queries. */
    public void expectShapeValid(S2Shape shape) {
      S2ShapeIndex index = new S2ShapeIndex();
      index.add(shape);

      S2Error error = new S2Error();
      for (S2ValidationQueryBase query : queries) {
        boolean result = query.validate(index, error);
        String message =
            "Expected query "
                + query.getClass().getSimpleName()
                + " to confirm this geometry as valid '"
                + S2TextFormat.toString(index)
                + "' but it returned error '"
                + error
                + "'";
        assertTrue(message, result);
        assertTrue(message, error.ok());
        error.clear();
      }
    }

    /**
     * Checks that a particular shape fails to validate in each of the configured queries, and that
     * the validation error matches one of the given codes.
     */
    public void expectShapeInvalid(S2Shape shape, S2Error.Code... codes) {
      List<S2Error.Code> codesList = ImmutableList.copyOf(codes);
      S2ShapeIndex index = new S2ShapeIndex();
      index.add(shape);

      S2Error error = new S2Error();
      for (S2ValidationQueryBase query : queries) {
        boolean result = query.validate(index, error);
        String message =
            "Expected query "
                + query.getClass().getSimpleName()
                + " to find an error in '"
                + Joiner.on(",").join(codes)
                + "' in shape '"
                + S2TextFormat.toString(index)
                + "' but it did not.";
        assertFalse(message, result);
        assertTrue(message, codesList.contains(error.code()));
        error.clear();
      }
    }

    /** Checks that the given geometry string validates in each of the configured queries. */
    public void expectGeometryValid(String geometry) {
      S2ShapeIndex index = S2TextFormat.makeIndexOrDie(geometry);
      S2Error error = new S2Error();
      for (S2ValidationQueryBase query : queries) {
        boolean result = query.validate(index, error);
        String message =
            "Expected query "
                + query.getClass().getSimpleName()
                + " to confirm this geometry as valid '"
                + geometry
                + "' but it returned error '"
                + error
                + "'";
        assertTrue(message, result);
        assertTrue(message, error.ok());
        error.clear();
      }
    }

    /**
     * Checks that the given geometry string does not validate in each of the configured queries,
     * and that the given error code is returned.
     */
    public void expectGeometryInvalid(String geometry, S2Error.Code code) {
      S2ShapeIndex index = S2TextFormat.makeIndexOrDie(geometry);
      S2Error error = new S2Error();
      for (S2ValidationQueryBase query : queries) {
        boolean result = query.validate(index, error);
        String message =
            "Expected query "
                + query.getClass().getSimpleName()
                + " to find error '"
                + code
                + "' in geometry '"
                + geometry
                + "' but it did not.";
        assertFalse(message, result);
        assertEquals(message, code, error.code());
        error.clear();
      }
    }

    /**
     * Creates and returns a new polygon from the current loops after shuffling their order. The
     * polygon may be invalid, so internal assertions are disabled.
     */
    S2Polygon makePolygon() {
      List<S2Loop> shuffledLoops = new ArrayList<>();
      for (List<S2Point> loop : loops) {
        shuffledLoops.add(new S2Loop(loop));
      }

      Collections.shuffle(loops);
      return GeometryTestCase.uncheckedCreate(() -> new S2Polygon(shuffledLoops));
    }
  }

  // An S2Shape implementation with one chain that's open.
  public static class OpenShape implements S2Shape {
    @Override
    public int dimension() {
      return 2;
    }

    @Override
    public int numEdges() {
      return 3;
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge edge) {
      Preconditions.checkArgument(chainId == 0);
      getEdge(offset, edge);
    }

    @Override
    public void getEdge(int id, MutableEdge edge) {
      Preconditions.checkArgument(id >= 0);
      Preconditions.checkArgument(id < 3);

      switch (id) {
        case 0:
          edge.set(diamondPoints[0], diamondPoints[1]);
          break;
        case 1:
          edge.set(diamondPoints[1], diamondPoints[2]);
          break;
        case 2:
          edge.set(diamondPoints[2], diamondPoints[3]);
          break;
        default:
          throw new IllegalArgumentException();
      }
    }

    @Override
    public S2Point getChainVertex(int chainId, int offset) {
      Preconditions.checkArgument(chainId == 0);
      switch (offset) {
        case 0:
          return diamondPoints[0];
        case 1:
          return diamondPoints[1];
        case 2:
          return diamondPoints[2];
        case 3:
          return diamondPoints[3];
        default:
          Preconditions.checkArgument(false);
      }
      return null;
    }

    @Override
    public int numChains() {
      return 1;
    }

    @Override
    public int getChainStart(int chainId) {
      return 0;
    }

    @Override
    public int getChainLength(int chainId) {
      Preconditions.checkArgument(chainId == 0);
      return numEdges();
    }

    @Override
    public ReferencePoint getReferencePoint() {
      return ReferencePoint.create(new S2Point(1, 0, 0), true);
    }

    @Override
    public boolean containsOrigin() {
      return false;
    }

    @Override
    public void getChainPosition(int id, ChainPosition pos) {
      pos.set(0, id);
    }

    @Override
    public boolean hasInterior() {
      return false;
    }
  }

  /**
   * An S2Shape subclass returning an invalid dimension value. It has one chain with no edges, which
   * would typically represent a full shape.
   */
  public static class BadDimensionShape implements S2Shape {
    @Override
    public int dimension() {
      return 42;
    }

    @Override
    public int numEdges() {
      return 0;
    }

    @Override
    public void getEdge(int id, MutableEdge edge) {
      throw new IllegalArgumentException();
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge edge) {
      throw new IllegalArgumentException();
    }

    @Override
    public int getChainLength(int chainId) {
      Preconditions.checkArgument(chainId == 0);
      return 0;
    }

    @Override
    public int getChainStart(int chainId) {
      Preconditions.checkArgument(chainId == 0);
      return 0; // This doesn't actually make sense. S2Shape doesn't define what getChainStart
      // should return when a chain has no edges.
    }

    @Override
    public int numChains() {
      return 1;
    }

    @Override
    public void getChainPosition(int edgeId, ChainPosition pos) {
      throw new IllegalArgumentException();
    }

    @Override
    public S2Point getChainVertex(int chainId, int offset) {
      Preconditions.checkArgument(false);
      return null;
    }

    @Override
    public ReferencePoint getReferencePoint() {
      return ReferencePoint.create(new S2Point(1, 0, 0), true);
    }

    @Override
    public boolean hasInterior() {
      return true;
    }

    @Override
    public boolean containsOrigin() {
      return true;
    }
  }

  /**
   * An S2Shape subclass returning a chain that's too short: i.e. numEdges is 2, numChains is 1,
   * getChainLength(0) is 2, but actually there is only one edge. Attempting to get chain edges or
   * vertices fails.
   */
  public static class BadChainLengthShape implements S2Shape {
    @Override
    public int dimension() {
      return 2;
    }

    @Override
    public int numEdges() {
      return 2;
    }

    @Override
    public void getEdge(int id, MutableEdge edge) {
      Preconditions.checkArgument(id >= 0);
      Preconditions.checkArgument(id < 2);

      switch (id) {
        case 0:
          edge.set(diamondPoints[0], diamondPoints[1]);
          break;
        case 1:
          edge.set(diamondPoints[1], diamondPoints[2]);
          break;
        default:
          throw new IllegalArgumentException();
      }
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge edge) {}

    @Override
    public int numChains() {
      return 1;
    }

    @Override
    public int getChainStart(int chainId) {
      return 0;
    }

    @Override
    public int getChainLength(int chainId) {
      return numEdges();
    }

    @Override
    public void getChainPosition(int edgeId, ChainPosition pos) {
      // Bad shape, does nothing
    }

    @Override
    public S2Point getChainVertex(int chainId, int offset) {
      // Bad shape, does nothing
      return null;
    }

    // Claims to contain the reference point, but that's wrong and impossible for a line.
    @Override
    public ReferencePoint getReferencePoint() {
      return ReferencePoint.create(new S2Point(1, 0, 0), true);
    }

    @Override
    public boolean containsOrigin() {
      return false;
    }

    @Override
    public boolean hasInterior() {
      return false;
    }
  }

  /**
   * Defines a grid where points are separated by 15 degrees in longitude and 30 degrees in
   * latitude. The integral indices for x thus range from [0, 24) and y ranges from [0, 13]. X
   * values are automatically wrapped into range.
   */
  public static S2Point gridPoint(int x, int y) {
    double kLatDelta = 15;  // degrees;
    double kLonDelta = 15;  // degrees;

    Preconditions.checkArgument(0 <= y && y <= 12);
    Preconditions.checkArgument(x >= 0);
    x %= 24;

    // Return exact polar points to avoid numerical issues.
    if (y == 0) {
      return new S2Point(0, 0, -1);
    }
    if (y == 12) {
      return new S2Point(0, 0, +1);
    }
    return S2LatLng.fromDegrees(-90 + kLatDelta * y, -180 + kLonDelta * x).toPoint();
  }

  /**
   * Builds a "quilt" test shape which is a series of rings that stretch from the south to north
   * pole and form a grid where every vertex has at least two chains incident on it.
   */
  public static S2LaxPolygonShape makeQuilt() {
    // Build loops that touch at every-other point. This will give us a diamond quilt pattern where
    // every vertex has two chains incident on it.
    List<List<S2Point>> loops = new ArrayList<>();

    for (int x = 0; x < 24; x += 2) {
      for (int y = 0; y < 12; y += 2) {
        List<S2Point> loop = new ArrayList<>();
        loop.add(gridPoint(x + 0, y + 1));
        loop.add(gridPoint(x + 1, y + 2));
        loop.add(gridPoint(x + 2, y + 1));
        loop.add(gridPoint(x + 1, y + 0));
        loops.add(loop);
      }
    }
    return S2LaxPolygonShape.create(loops);
  }

  private ValidationQueriesTestUtils() {}
}
