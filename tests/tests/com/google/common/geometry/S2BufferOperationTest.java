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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2EdgeUtil.getUpdateMinDistanceMaxError;
import static com.google.common.geometry.S2ShapeUtil.s2PointLoopList;
import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static com.google.common.geometry.S2TextFormat.makeLaxPolylineOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.parsePointsOrDie;
import static java.lang.Math.PI;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2BooleanOperation.PolygonModel;
import com.google.common.geometry.S2BooleanOperation.PolylineModel;
import com.google.common.geometry.S2BufferOperation.EndCapStyle;
import com.google.common.geometry.S2BufferOperation.PolylineSide;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class S2BufferOperationTest extends GeometryTestCase {

  /**
   * Convenience function that calls the given lambda to add input to an S2BufferOperation and
   * returns the buffered result.
   */
  private static S2LaxPolygonShape doBuffer(
      Consumer<S2BufferOperation> inputCallback, S2BufferOperation.Options options) {
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
    S2BufferOperation op = new S2BufferOperation(layer, options);
    inputCallback.accept(op);
    S2Error error = new S2Error();
    boolean unused = op.build(error);
    assertTrue(error.text(), error.ok());

    S2LaxPolygonShape output = layer.getPolygon();
    // if (ABSL_VLOG_IS_ON(1) && output.numVertices() < 1000) {
    System.err.println("\nS2Polygon: " + S2TextFormat.toString(output) + "\n");
    // }
    return output;
  }

  /** Simpler version that accepts a buffer radius and error fraction. */
  private static S2LaxPolygonShape doBuffer(
      Consumer<S2BufferOperation> inputCallback, S1Angle bufferRadius, double errorFraction) {
    S2BufferOperation.Options options = new S2BufferOperation.Options();
    options.setBufferRadius(bufferRadius);
    options.setErrorFraction(errorFraction);
    return doBuffer(inputCallback, options);
  }

  /**
   * Given a callback that adds empty geometry to the S2BufferOperation, verifies that the result is
   * empty after buffering.
   */
  private static void verifyBufferEmpty(Consumer<S2BufferOperation> input) {
    // Test code paths that involve buffering by negative, zero, and positive
    // values, and also values where the result is usually empty or full.
    assertTrue(doBuffer(input, S1Angle.degrees(-200), 0.1).isEmpty());
    assertTrue(doBuffer(input, S1Angle.degrees(-1), 0.1).isEmpty());
    assertTrue(doBuffer(input, S1Angle.degrees(0), 0.1).isEmpty());
    assertTrue(doBuffer(input, S1Angle.degrees(1), 0.1).isEmpty());
    assertTrue(doBuffer(input, S1Angle.degrees(200), 0.1).isEmpty());
  }

  @Test
  public void testOptions() {
    // Provide test coverage for `options()`.
    S2BufferOperation.Options options = new S2BufferOperation.Options();
    assertEquals(EndCapStyle.ROUND, options.endCapStyle());
    options.setEndCapStyle(EndCapStyle.FLAT);
    assertEquals(EndCapStyle.FLAT, options.endCapStyle());

    options = new S2BufferOperation.Options(S1Angle.radians(1e-12));
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
    S2BufferOperation op = new S2BufferOperation(layer, options);
    assertEquals(options.bufferRadius(), op.options().bufferRadius());
  }

  @Test
  public void testNoInput() {
    verifyBufferEmpty(op -> {});
  }

  @Test
  public void testEmptyPolyline() {
    // Note that polylines with 1 vertex are defined to have no edges.
    verifyBufferEmpty(
        op -> {
          op.addPolyline(ImmutableList.of(new S2Point(1, 0, 0)));
        });
  }

  @Test
  public void testEmptyLoop() {
    verifyBufferEmpty(
        op -> {
          op.addLoop(ImmutableList.of());
        });
  }

  @Test
  public void testEmptyPointShape() {
    verifyBufferEmpty(
        op -> {
          op.addShape(S2Point.Shape.fromList(ImmutableList.of()));
        });
  }

  @Test
  public void testEmptyPolylineShape() {
    verifyBufferEmpty(
        op -> {
          op.addShape(S2TextFormat.makeLaxPolylineOrDie(""));
        });
  }

  @Test
  public void testEmptyPolygonShape() {
    verifyBufferEmpty(
        op -> {
          op.addShape(S2TextFormat.makeLaxPolygonOrDie(""));
        });
  }

  @Test
  public void testEmptyShapeIndex() {
    verifyBufferEmpty(
        op -> {
          op.addShapeIndex(S2TextFormat.makeIndexOrDie("# #"));
        });
  }

  @Test
  public void testPoorlyNormalizedPoint() {
    // Verify that debugging assertions are not triggered when an input point is not unit length
    // (but within the limits guaranteed by S2Point.normalize).
    //
    // The purpose of this test is not to check that the result is correct (which is done
    // elsewhere), simply that no assertions occur.
    S2LaxPolygonShape unused = doBuffer(
        op -> {
          S2Point p = new S2Point(1 - 2 * DBL_EPSILON, 0, 0); // Maximum error allowed.
          assertTrue(S2.isUnitLength(p));
          op.addPoint(p);
        },
        S1Angle.degrees(1),
        0.01);
  }

  /**
   * Given a callback that adds the full polygon to the S2BufferOperation, verifies that the result
   * is full after buffering.
   */
  void verifyBufferFull(Consumer<S2BufferOperation> input) {
    // Test code paths that involve buffering by negative, zero, and positive
    // values, and also values where the result is usually empty or full.
    assertTrue(doBuffer(input, S1Angle.degrees(-200), 0.1).isFull());
    assertTrue(doBuffer(input, S1Angle.degrees(-1), 0.1).isFull());
    assertTrue(doBuffer(input, S1Angle.degrees(0), 0.1).isFull());
    assertTrue(doBuffer(input, S1Angle.degrees(1), 0.1).isFull());
    assertTrue(doBuffer(input, S1Angle.degrees(200), 0.1).isFull());
  }

  @Test
  public void testFullPolygonShape() {
    verifyBufferFull(
        op -> {
          op.addShape(S2TextFormat.makeLaxPolygonOrDie("full"));
        });
  }

  @Test
  public void testFullShapeIndex() {
    verifyBufferFull(
        op -> {
          op.addShapeIndex(S2TextFormat.makeIndexOrDie("# # full"));
        });
  }

  @Test
  public void testPointsAndPolylinesAreRemoved() {
    // Test that points and polylines are removed with a negative buffer radius.
    S2LaxPolygonShape output =
        doBuffer(
            op -> {
              op.addShapeIndex(S2TextFormat.makeIndexOrDie("0:0 # 2:2, 2:3#"));
            },
            S1Angle.degrees(-1),
            0.1);
    assertTrue(output.isEmpty());
  }

  @Test
  public void testBufferedPointsAreSymmetric() {
    // Test that points are buffered into regular polygons. (This is not
    // guaranteed by the API but makes the output nicer to look at. :)
    S2LaxPolygonShape output =
        doBuffer(
            op -> {
              op.addPoint(new S2Point(1, 0, 0));
            },
            S1Angle.degrees(5),
            0.001234567);

    // We use the length of the last edge as our reference length.
    int n = output.numVertices();
    S1Angle edgeLength = new S1Angle(output.getChainVertex(0, 0), output.getChainVertex(0, n - 1));
    for (int i = 1; i < n; ++i) {
      assertLessOrEqual(
          edgeLength
              .sub(new S1Angle(output.getChainVertex(0, i - 1), output.getChainVertex(0, i)))
              .abs(),
          S1Angle.radians(1e-14));
    }
  }

  @Test
  public void testSetCircleSegments() {
    // Test that when a point is buffered with a small radius the number of edges matches
    // options.circleSegments(). (This is not true for large radii because large circles on the
    // sphere curve less than 360 degrees.) Using a tiny radius helps to catch rounding problems.
    S2BufferOperation.Options options = new S2BufferOperation.Options(S1Angle.radians(1e-12));
    for (int circleSegments = 3; circleSegments <= 20; ++circleSegments) {
      options.setCircleSegments(circleSegments);
      // The error on round-trip of cos() and acos() is fairly large. C++ checks that these are
      // almost equal only as 32-bit floats.
      assertEquals((double) circleSegments, options.circleSegments(), 1e-8);
      S2LaxPolygonShape output =
          doBuffer(
              op -> {
                op.addPoint(new S2Point(1, 0, 0));
              },
              options);
      assertEquals(circleSegments, output.numVertices());
    }
  }

  @Test
  public void testSetSnapFunction() {
    // Verify that the snap function is passed through to S2Builder.
    // We use the default buffer radius of zero to make the test simpler.
    S2BufferOperation.Options options = new S2BufferOperation.Options();
    options.setSnapFunction(new IntLatLngSnapFunction(0));
    S2LaxPolygonShape output =
        doBuffer(
            op -> op.addPoint(makePointOrDie("0.1:-0.4")),
            options);
    assertEquals(1, output.numVertices());
    assertEquals(makePointOrDie("0:0"), output.getChainVertex(0, 0));
  }

  @Test
  public void testNegativeBufferRadiusMultipleLayers() {
    // Verify that with a negative buffer radius, at most one polygon layer is allowed.
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
    S2BufferOperation op =
        new S2BufferOperation(layer, new S2BufferOperation.Options(S1Angle.radians(-1)));
    op.addLoop(s2PointLoopList(parsePointsOrDie("0:0, 0:1, 1:0")));
    op.addShapeIndex(makeIndexOrDie("# # 2:2, 2:3, 3:2"));
    S2Error error = new S2Error();
    assertFalse(op.build(error));
    assertEquals(S2Error.Code.FAILED_PRECONDITION, error.code());
  }

  /**
   * If bufferRadius > maxError, tests that "output" contains "input". If bufferRadius < -maxError
   * tests that "input" contains "output". Otherwise does nothing.
   */
  private static void verifyContainment(
      S2ShapeIndex input, S2ShapeIndex output, S1Angle bufferRadius, S1Angle maxError) {
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    builder.setPolygonModel(PolygonModel.CLOSED);
    builder.setPolylineModel(PolylineModel.CLOSED);

    if (bufferRadius.greaterThan(maxError)) {
      // For positive buffer radii, the output should contain the input.
      assertTrue(S2BooleanOperation.contains(output, input, builder.options()));
    } else if (bufferRadius.lessThan(maxError.neg())) {
      // For negative buffer radii, the input should contain the output.
      assertTrue(S2BooleanOperation.contains(input, output, builder.options()));
    }
  }

  /**
   * Tests that the minimum distance from the boundary of "output" to the boundary of "input" is at
   * least "minDist" using exact predicates.
   */
  void verifyMinimumDistance(S2ShapeIndex input, S2ShapeIndex output, S1ChordAngle minDist) {
    if (minDist.isZero()) {
      return;
    }

    // We do one query to find the edges of "input" that might be too close to "output", then for
    // each such edge we do another query to find the edges of "output" that might be too close to
    // it. Then we check the distance between each edge pair (A, B) using exact predicates.

    // We make the distance limit big enough to find all edges whose true distance might be less
    // than "minDist".
    S2ClosestEdgeQuery.Builder queryBuilder = new S2ClosestEdgeQuery.Builder();
    queryBuilder.setIncludeInteriors(false);
    queryBuilder.setMaxDistance(minDist.plusError(getUpdateMinDistanceMaxError(minDist)));

    S2ClosestEdgeQuery<S1ChordAngle> inQuery = queryBuilder.build(input);
    S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> outTarget =
        S2ClosestEdgeQuery.createShapeIndexTarget(output);
    outTarget.setIncludeInteriors(false);
    S2ClosestEdgeQuery<S1ChordAngle> outQuery = queryBuilder.build(output);
    MutableEdge a = new MutableEdge();
    MutableEdge b = new MutableEdge();
    inQuery.findClosestEdges(
        outTarget,
        (S1ChordAngle aUnused, int aShapeId, int aEdgeId) -> {
          input.getShapes().get(aShapeId).getEdge(aEdgeId, a);
          S2ClosestEdgeQuery.EdgeTarget<S1ChordAngle> inTarget =
              new S2ClosestEdgeQuery.EdgeTarget<>(a.getStart(), a.getEnd());
          outQuery.findClosestEdges(
              inTarget,
              (S1ChordAngle bUnused, int bShapeId, int bEdgeId) -> {
                output.getShapes().get(bShapeId).getEdge(bEdgeId, b);
                // The distance from edge A to edge B must be greater than or equal to minDist.
                assertFalse(
                    S2EdgeUtil.isEdgePairDistanceLess(
                        a.getStart(), a.getEnd(), b.getStart(), b.getEnd(), minDist));
                return true;
              });
          return true;
        });
  }

  /**
   * Tests that the Hausdorff distance from the boundary of "output" to the boundary of "input" is
   * at most (1 + errorFraction) * bufferRadius. The implementation approximates this by measuring
   * the distance at a set of points along the boundary of "output".
   */
  private static void verifyHausdorffDistance(
      S2ShapeIndex input, S2ShapeIndex output, S1ChordAngle maxDist) {
    S2ClosestEdgeQuery.Builder queryBuilder = new S2ClosestEdgeQuery.Builder();
    queryBuilder.setIncludeInteriors(false);
    queryBuilder.setMaxDistance(maxDist.plusError(getUpdateMinDistanceMaxError(maxDist)));

    S2ClosestEdgeQuery<S1ChordAngle> inQuery = queryBuilder.build(input);
    MutableEdge e = new MutableEdge();
    MutableEdge a = new MutableEdge();
    for (S2Shape outShape : output.getShapes()) {
      for (int i = 0; i < outShape.numEdges(); ++i) {
        outShape.getEdge(i, e);
        // Measure the distance at 5 points along the edge.
        for (double t = 0; t <= 1.0; t += 0.25) {
          S2Point b = S2EdgeUtil.interpolate(e.getStart(), e.getEnd(), t);
          S2ClosestEdgeQuery.PointTarget<S1ChordAngle> outTarget =
              new S2ClosestEdgeQuery.PointTarget<>(b);
          // We check the distance bound using exact predicates.
          inQuery.findClosestEdges(
              outTarget,
              (dist, shapeId, edgeId) -> {
                input.getShapes().get(shapeId).getEdge(edgeId, a);
                int cmp =
                    S2Predicates.compareEdgeDistance(
                        b, a.getStart(), a.getEnd(), maxDist.getLength2());
                assertLessOrEqual(cmp, 0);
                return true;
              });
        }
      }
    }
  }

  /**
   * Buffers the given input with the given bufferRadius and errorFraction and verifies that the
   * output is correct.
   */
  void verifyBuffer(S2ShapeIndex input, S1Angle bufferRadius, double errorFraction) {
    // Ideally we would verify the correctness of buffering as follows. Suppose
    // that B = Buffer(A, r) and let ~X denote the complement of region X. Then
    // if r > 0, we would verify:
    //
    //   1a. Minimum distance between ~B and A >= r_min
    //   1b. Directed Hausdorff distance from B to A <= r_max
    //
    // Buffering A by r < 0 is equivalent to buffering ~A by |r|, so instead we
    // would verify the following (where r_min and r_max are based on |r|):
    //
    //   2a. Minimum distance between B and ~A >= r_min
    //   2b. Directed Hausdorff distance from ~B to ~A <= r_max
    //
    // Conditions 1a and 2a can be implemented as follows:
    //
    //   1a*: B.Contains(A) && minimum distance between @B and @A >= r_min
    //   2a*: A.Contains(B) && minimum distance between @B and @A >= r_min
    //
    // Note that if r_min <= 0 then there is nothing to be tested, since the
    // containment condition may not hold. (Note that even when the specified
    // buffer radius is zero, edges can move slightly when crossing edges are
    // split during the snapping step.)  The correct approach would be to test
    // instead that the directed Hausdorff distance from A to ~B is at most
    // -r_min, but Hausdorff distance is not yet implemented.
    //
    // Similarly, conditions 1b and 2b need to be approximated because Hausdorff
    // distance is not yet implemented. We do this by measuring the distance at
    // a set of points on the boundary of B:
    //
    //   1b*: Minimum distance from P to @A <= r_max for a set of points P on @B
    //   2b*: Minimum distance from P to @A <= r_max for a set of points P on @B
    //
    // This is not perfect (e.g., it won't detect cases where an entire boundary
    // loop of B is missing, such as returning a disc in the place of an
    // annulus) but it is sufficient to detect many types of errors.
    S2BufferOperation.Options options = new S2BufferOperation.Options();
    options.setBufferRadius(bufferRadius);
    options.setErrorFraction(errorFraction);
    S2ShapeIndex output = new S2ShapeIndex();
    output.add(doBuffer(op -> op.addShapeIndex(input), options));

    //String trace = Platform.formatString(
    //        "\nradius = %.17g, errorFraction = %.17g\ninput = %s\noutput = %s",
    //        bufferRadius.radians(),
    //        errorFraction,
    //        S2TextFormat.toString(input),
    //        S2TextFormat.toString(output));
    //System.err.println(trace);

    // Check the 1a*/1b* condition above.
    S1Angle maxError = options.maxError();
    verifyContainment(input, output, bufferRadius, maxError);

    S1ChordAngle minDist = S1ChordAngle.fromS1Angle(S1Angle.max(S1Angle.ZERO, bufferRadius.abs().sub(maxError)));
    verifyMinimumDistance(input, output, minDist);

    // Check the 2a*/2b* condition (i.e., directed Hausdorff distance).
    S1ChordAngle maxDist = S1ChordAngle.fromS1Angle(bufferRadius.abs().add(maxError));
    verifyHausdorffDistance(input, output, maxDist);
  }

  /** Convenience function that takes an S2ShapeIndex in S2TextFormat format. */
  void verifyBuffer(String indexStr, S1Angle bufferRadius, double errorFraction) {
    verifyBuffer(S2TextFormat.makeIndexOrDie(indexStr), bufferRadius, errorFraction);
  }

  /** Convenience function that tests buffering using +/- the given radius. */
  void verifySignedBuffer(String indexStr, S1Angle bufferRadius, double errorFraction) {
    verifyBuffer(indexStr, bufferRadius, errorFraction);
    verifyBuffer(indexStr, bufferRadius.neg(), errorFraction);
  }

  @Test
  public void testPointShell() {
    verifySignedBuffer("# # 0:0", S1Angle.radians(M_PI_2), 0.01);
  }

  @Test
  public void testSiblingPairShell() {
    verifySignedBuffer("# # 0:0, 0:5", S1Angle.radians(M_PI_2), 0.01);
  }

  @Test
  public void testSiblingPairHole() {
    verifySignedBuffer("# # 0:0, 0:10, 7:7; 3:4, 3:6", S1Angle.degrees(1), 0.01);
  }

  @Test
  public void testSquare() {
    verifySignedBuffer("# # -3:-3, -3:3, 3:3, 3:-3", S1Angle.degrees(1), 0.01);
    verifySignedBuffer("# # -3:-3, -3:3, 3:3, 3:-3", S1Angle.degrees(170), 1e-4);
  }

  @Test
  public void testHollowSquare() {
    verifySignedBuffer(
        "# # -3:-3, -3:3, 3:3, 3:-3; 2:2, -2:2, -2:-2, 2:-2", S1Angle.degrees(1), 0.01);
  }

  @Test
  public void testZigZagLoop() {
    verifySignedBuffer("# # 0:0, 0:7, 5:3, 5:10, 6:10, 6:1, 1:5, 1:0", S1Angle.degrees(0.2), 0.01);
  }

  @Test
  public void testFractals() {
    TestDataGenerator data = new TestDataGenerator();
    for (double dimension : new double[] {1.02, 1.8}) {
      S2FractalBuilder fractal = new S2FractalBuilder(data.rand);
      fractal.setLevelForApproxMaxEdges(3 * 64);
      fractal.setFractalDimension(dimension);
      S2Loop loop = fractal.makeLoop(S2.getFrame(new S2Point(1, 0, 0)), S1Angle.degrees(10));
      S2ShapeIndex input = new S2ShapeIndex();
      input.add(loop);
      verifyBuffer(input, S1Angle.degrees(0.4), 0.01);
    }
  }

  @Test
  public void testS2Curve() {
    // Tests buffering the S2 curve by an amount that yields the full polygon.
    int kLevel = 2; // Number of input edges == 6 * (4 ** kLevel)
    List<S2Point> points = new ArrayList<>();
    for (S2CellId id = S2CellId.begin(kLevel); !id.equals(S2CellId.end(kLevel)); id = id.next()) {
      points.add(id.toPoint());
    }
    // Duplicate the beginning point at the end to close the curve.
    points.add(points.get(0));

    // Buffering by this amount or more is guaranteed to yield the full polygon.
    // (Note that the bound is not tight for S2CellIds at low levels.)
    S1Angle fullRadius = S1Angle.radians(0.5 * S2Projections.MAX_DIAG.getValue(kLevel));
    assertTrue(
        doBuffer(
                op -> op.addShape(S2LaxPolylineShape.create(points)),
                fullRadius,
                0.1)
            .isFull());
  }

  /**
   * Tests buffering the given S2ShapeIndex with a variety of radii and error fractions. This method
   * is intended to be used with relatively simple shapes since calling it is quite expensive.
   */
  void verifyRadiiAndErrorFractions(String indexStr) {
    // Try the full range of radii with a representative error fraction.
    double kFrac = 0.01;
    double[] kTestRadiiRadians =
        new double[] {
          0,
          1e-300,
          1e-15,
          2e-15,
          3e-15,
          1e-5,
          0.01,
          0.1,
          1.0,
          (1 - kFrac) * M_PI_2,
          M_PI_2 - 1e-15,
          M_PI_2,
          M_PI_2 + 1e-15,
          (1 - kFrac) * PI,
          PI - 1e-6,
          PI,
          1e300
        };
    for (double radius : kTestRadiiRadians) {
      verifySignedBuffer(indexStr, S1Angle.radians(radius), kFrac);
    }

    // Now try the full range of error fractions with a few selected radii.
    double[] kTestErrorFractions =
        new double[] {S2BufferOperation.Options.MIN_ERROR_FRACTION, 0.001, 0.01, 0.1, 1.0};
    for (double errorFraction : kTestErrorFractions) {
      verifyBuffer(indexStr, S1Angle.radians(-1e-6), errorFraction);
      verifyBuffer(indexStr, S1Angle.radians(1e-14), errorFraction);
      verifyBuffer(indexStr, S1Angle.radians(1e-2), errorFraction);
      verifyBuffer(indexStr, S1Angle.radians(PI - 1e-3), errorFraction);
    }
  }

  @Test
  public void testRadiiAndErrorFractionCoverage() {
    // Test buffering simple shapes with a wide range of different buffer radii and error fractions.

    // A single point.
    verifyRadiiAndErrorFractions("1:1 # #");

    // A zig-zag polyline.
    verifyRadiiAndErrorFractions("# 0:0, 0:30, 30:30, 30:60 #");

    // A triangular polygon with a triangular hole. (The hole is clockwise.)
    verifyRadiiAndErrorFractions("# # 0:0, 0:100, 70:50; 10:20, 50:50, 10:80");

    // A triangle with one very short and two very long edges.
    verifyRadiiAndErrorFractions("# # 0:0, 0:179.99999999999, 1e-300:0");
  }

  /** A helper class for verifying that polylines are buffered correctly. */
  static class VerifyBufferPolyline {

    private static final double ARC_LO = 0.001;
    private static final double ARC_HI = 0.999;
    private static final int ARC_SAMPLES = 7;

    private final List<S2Point> polyline;
    private final S2ShapeIndex output;
    private final S1Angle bufferRadius;
    private final S1Angle maxError;
    private final S1ChordAngle minDist;
    private final S1ChordAngle maxDist;
    private final boolean round;
    private final boolean twoSided;

    /**
     * Tests buffering a polyline with the given options. This method is intended only for testing
     * Options.EndCapStyle and Options.PolylineSide; if these options have their default values then
     * verifyBuffer() should be used instead. Similarly verifyBuffer should be used to test negative
     * buffer radii and polylines with 0 or 1 vertices.
     */
    @CanIgnoreReturnValue
    public VerifyBufferPolyline(String inputStr, S2BufferOperation.Options options) {
      List<S2Point> parsedPoints = parsePointsOrDie(inputStr);

      // Left-sided buffering is tested by reversing the polyline and then testing whether it has
      // been buffered correctly on the right.
      this.polyline =
          (options.polylineSide() == PolylineSide.LEFT)
              ? Lists.reverse(parsedPoints)
              : parsedPoints;
      this.output = new S2ShapeIndex();
      this.bufferRadius = options.bufferRadius();
      this.maxError = options.maxError();
      this.minDist =
          S1ChordAngle.fromS1Angle(S1Angle.max(S1Angle.ZERO, bufferRadius.sub(maxError)));
      this.maxDist = S1ChordAngle.fromS1Angle(bufferRadius.add(maxError));
      this.round = (options.endCapStyle() == EndCapStyle.ROUND);
      this.twoSided = (options.polylineSide() == PolylineSide.BOTH);

      assertGreaterOrEqual(polyline.size(), 2);
      assertTrue(bufferRadius.greaterThan(S1Angle.ZERO));

      S2ShapeIndex input = new S2ShapeIndex();
      input.add(makeLaxPolylineOrDie(inputStr));
      output.add(doBuffer(op -> op.addShapeIndex(input), options));

      // Even with one-sided buffering and flat end caps the Hausdorff distance criterion should
      // still be true. (This checks that the buffered result is never further than (bufferRadius +
      // maxError) from the input.)
      verifyHausdorffDistance(input, output, maxDist);

      // However the minimum distance criterion is different; it only applies to the portions of the
      // boundary that are buffered using the given radius. We check this approximately by walking
      // along the polyline and checking that (1) on portions of the polyline that should be
      // buffered, the output contains the offset point at distance (bufferRadius - maxError) and
      // (2) on portions of the polyline that should not be buffered, the output does not contain
      // the offset point at distance maxError. The tricky part is that both of these conditions
      // have exceptions: (1) may not hold if the test point is closest to the non-buffered side of
      // the polyline (see the last caveat in the documentation for Options.polylineSide), and (2)
      // may not hold if the test point is within (bufferRadius + maxError) of the buffered side of
      // any portion of the polyline.
      if (minDist.isZero()) {
        return;
      }

      int n = polyline.size();
      S2Point start0 = polyline.get(0);
      S2Point start1 = polyline.get(1);
      S2Point startBegin = getEdgeAxis(start0, start1);
      S2Point startMid = start0.crossProd(startBegin);
      verifyVertexArc(start0, startBegin, startMid, round && twoSided);
      verifyVertexArc(start0, startMid, startBegin.neg(), round);
      for (int i = 0; i < n - 2; ++i) {
        verifyEdgeAndVertex(polyline.get(i), polyline.get(i + 1), polyline.get(i + 2), true);
      }
      S2Point end0 = polyline.get(n - 1);
      S2Point end1 = polyline.get(n - 2);
      S2Point endBegin = getEdgeAxis(end0, end1);
      S2Point endMid = end0.crossProd(endBegin);
      verifyEdgeArc(endBegin, end1, end0, true);
      verifyVertexArc(end0, endBegin, endMid, round);
      verifyVertexArc(end0, endMid, endBegin.neg(), round && twoSided);
      for (int i = n - 3; i >= 0; --i) {
        verifyEdgeAndVertex(polyline.get(i + 2), polyline.get(i + 1), polyline.get(i), twoSided);
      }
      verifyEdgeArc(startBegin, start1, start0, twoSided);
    }

    private S2Point getEdgeAxis(S2Point a, S2Point b) {
      return S2RobustCrossProd.robustCrossProd(a, b).normalize();
    }

    private boolean pointBufferingUncertain(S2Point p, boolean expectContained) {
      // The only case where a point might be excluded from the buffered output is if it is on the
      // unbuffered side of the polyline.
      if (expectContained && twoSided) {
        return false;
      }

      int n = polyline.size();
      for (int i = 0; i < n - 1; ++i) {
        S2Point a = polyline.get(i);
        S2Point b = polyline.get(i + 1);
        if (!twoSided) {
          // Ignore points on the buffered side if expectContained is true, and on the unbuffered
          // side if expectContained is false.
          if ((S2Predicates.sign(a, b, p) < 0) == expectContained) {
            continue;
          }
        }
        // TODO(ericv): Depending on how the erasing optimization is implemented, it might be
        // possible to add "&& expectContained" to the test below.
        if (round) {
          if (S2EdgeUtil.isDistanceLess(p, a, b, maxDist)) {
            return true;
          }
        } else {
          if (S2EdgeUtil.isInteriorDistanceLess(p, a, b, maxDist)) {
            return true;
          }
          if (i > 0 && new S1ChordAngle(p, a).lessThan(maxDist)) {
            return true;
          }
          if (i == n - 2 && new S1ChordAngle(p, b).lessThan(maxDist)) {
            return true;
          }
        }
      }
      return false;
    }

    private void verifyPoint(S2Point p, S2Point dir, boolean expectContained) {
      S2Point x = S2EdgeUtil.getPointOnRay(p, dir, expectContained ? bufferRadius.sub(maxError) : maxError);
      if (!pointBufferingUncertain(x, expectContained)) {
        assertEquals(expectContained, new S2ContainsPointQuery(output).contains(x));
      }
    }

    private void verifyVertexArc(S2Point p, S2Point start, S2Point end, boolean expectContained) {
      for (double t = ARC_LO; t < 1; t += (ARC_HI - ARC_LO) / ARC_SAMPLES) {
        S2Point dir = S2EdgeUtil.interpolate(start, end, t);
        verifyPoint(p, dir, expectContained);
      }
    }

    private void verifyEdgeArc(S2Point baAxis, S2Point a, S2Point b, boolean expectContained) {
      for (double t = ARC_LO; t < 1; t += (ARC_HI - ARC_LO) / ARC_SAMPLES) {
        S2Point p = S2EdgeUtil.interpolate(a, b, t);
        verifyPoint(p, baAxis, expectContained);
      }
    }

    private void verifyEdgeAndVertex(S2Point a, S2Point b, S2Point c, boolean expectContained) {
      S2Point baAxis = getEdgeAxis(b, a);
      S2Point cbAxis = getEdgeAxis(c, b);
      verifyEdgeArc(baAxis, a, b, expectContained);
      verifyVertexArc(b, baAxis, cbAxis, expectContained);
    }
  }

  @Test
  public void testZigZagPolyline() {
    S2BufferOperation.Options options = new S2BufferOperation.Options(S1Angle.degrees(1));
    for (PolylineSide polylineSide :
        new PolylineSide[] {PolylineSide.LEFT, PolylineSide.RIGHT, PolylineSide.BOTH}) {
      for (EndCapStyle endCapStyle : new EndCapStyle[] {EndCapStyle.ROUND, EndCapStyle.FLAT}) {
        String trace =
            Platform.formatString("polylineSide = %s, endCapStyle = %s", polylineSide, endCapStyle);
        System.err.println(trace);
        options.setPolylineSide(polylineSide);
        options.setEndCapStyle(endCapStyle);
        new VerifyBufferPolyline("0:0, 0:7, 5:3, 5:10", options);
        new VerifyBufferPolyline("10:0, 0:0, 5:1", options);
      }
    }
  }
}
