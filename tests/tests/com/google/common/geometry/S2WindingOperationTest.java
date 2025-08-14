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

import static com.google.common.geometry.S2TextFormat.makeLaxPolygonOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.parsePointsOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.S2BuilderUtil.IndexedLayer;
import com.google.common.geometry.S2WindingOperation.WindingRule;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class S2WindingOperationTest {

  /**
   * Verifies that the S2WindingOperation with the given arguments produces the given result. In
   * order to ensure that results are not affected by the cyclic order of the loop vertices or the
   * S2ShapeIndex loop ordering, we compute the symmetric difference using S2BooleanOperation. (We
   * don't use IndexMatchingLayer because that class does not distinguish empty from full polygons,
   * and we don't need its ability to match edge multiplicities here.)
   */
  private void expectWindingResult(
      S2WindingOperation.Options options,
      List<String> loopStrs,
      String refPointStr,
      int refWinding,
      WindingRule rule,
      String expectedStr) {
    S2ShapeIndex expected = new S2ShapeIndex();
    expected.add(makeLaxPolygonOrDie(expectedStr));

    S2ShapeIndex actual = new S2ShapeIndex();
    IndexedLayer<S2LaxPolygonLayer> indexedLayer = new IndexedLayer<>(actual, new S2LaxPolygonLayer());
    S2WindingOperation windingOp = new S2WindingOperation(indexedLayer, options);
    for (String loopStr : loopStrs) {
      windingOp.addLoop(parsePointsOrDie(loopStr));
    }

    S2Error error = new S2Error();
    boolean unused = windingOp.build(makePointOrDie(refPointStr), refWinding, rule, error);
    assertTrue(error.text(), error.ok());

    S2LaxPolygonLayer differenceLayer = new S2LaxPolygonLayer();
    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
    S2BooleanOperation diffOp =
        builder.build(S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, differenceLayer);

    unused = diffOp.build(actual, expected, error);
    assertTrue(error.text(), error.ok());
    S2LaxPolygonShape difference = differenceLayer.getPolygon();
    if (!difference.isEmpty()) {
      String report ="Expected no difference, got:\n" + S2TextFormat.toString(difference)
          + "\nActual=" + S2TextFormat.toString(actual)
          + "\nExpected=" + S2TextFormat.toString(expected)
          + "\n";
      assertTrue(report, difference.isEmpty());
    }
  }

  /**
   * Like expectWindingResult(), but with two different expected results depending on whether
   * options.includeDegeneracies() is false or true.
   */
  private void expectDegenerateWindingResult(
      S2WindingOperation.Options options,
      List<String> loopStrs,
      String refPointStr,
      int refWinding,
      WindingRule rule,
      String expectedStrFalse,
      String expectedStrTrue) {
    options.setIncludeDegeneracies(false);
    expectWindingResult(options, loopStrs, refPointStr, refWinding, rule, expectedStrFalse);
    options.setIncludeDegeneracies(true);
    expectWindingResult(options, loopStrs, refPointStr, refWinding, rule, expectedStrTrue);
  }

  @Test
  public void testEmpty() {
    expectWindingResult(
        new S2WindingOperation.Options(), ImmutableList.of(""), "5:5", 0, WindingRule.POSITIVE, "");
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of(""),
        "5:5",
        1,
        WindingRule.POSITIVE,
        "full");
  }

  @Test
  public void testNoDegeneraciesByDefault() {
    // By default, a degenerate input is ignored.
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("2:2"),
        "5:5",
        0,
        WindingRule.POSITIVE,
        "");
  }

  @Test
  public void testPointLoop() {
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("2:2"),
        "5:5",
        0,
        WindingRule.POSITIVE,
        "",
        "2:2");
  }

  @Test
  public void testSiblingPairLoop() {
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("2:2, 3:3"),
        "5:5",
        0,
        WindingRule.POSITIVE,
        "",
        "2:2, 3:3");
  }

  @Test
  public void testRectangle() {
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("0:0, 0:10, 10:10, 10:0"),
        "5:5",
        1,
        WindingRule.POSITIVE,
        "0:0, 0:10, 10:10, 10:0");
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("0:0, 0:10, 10:10, 10:0"),
        "5:5",
        1,
        WindingRule.NEGATIVE,
        "");
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("0:0, 0:10, 10:10, 10:0"),
        "5:5",
        1,
        WindingRule.NON_ZERO,
        "0:0, 0:10, 10:10, 10:0");
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("0:0, 0:10, 10:10, 10:0"),
        "5:5",
        1,
        WindingRule.ODD,
        "0:0, 0:10, 10:10, 10:0");
  }

  @Test
  public void testBowTie() {
    // Note that NEGATIVE, NON_ZERO, and ODD effectively reverse the orientation
    // of one of the two triangles that form the bow tie.
    expectWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(1))),
        ImmutableList.of("5:-5, -5:5, 5:5, -5:-5"),
        "10:0",
        0,
        WindingRule.POSITIVE,
        "0:0, -5:5, 5:5");
    expectWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(1))),
        ImmutableList.of("5:-5, -5:5, 5:5, -5:-5"),
        "10:0",
        0,
        WindingRule.NEGATIVE,
        "-5:-5, 0:0, 5:-5");
    expectWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(1))),
        ImmutableList.of("5:-5, -5:5, 5:5, -5:-5"),
        "10:0",
        0,
        WindingRule.NON_ZERO,
        "0:0, -5:5, 5:5; -5:-5, 0:0, 5:-5");
    expectWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(1))),
        ImmutableList.of("5:-5, -5:5, 5:5, -5:-5"),
        "10:0",
        0,
        WindingRule.ODD,
        "0:0, -5:5, 5:5; -5:-5, 0:0, 5:-5");
  }

  @Test
  public void testCollapsingShell() {
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(5))),
        ImmutableList.of("0:0, 0:3, 3:3"),
        "10:0",
        0,
        WindingRule.POSITIVE,
        "",
        "0:0");
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(5))),
        ImmutableList.of("0:0, 0:3, 3:3"),
        "1:1",
        1,
        WindingRule.POSITIVE,
        "",
        "0:0");
    expectWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(5))),
        ImmutableList.of("0:0, 3:3, 0:3"),
        "10:0",
        1,
        WindingRule.POSITIVE,
        "full");
    expectWindingResult(
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.degrees(5))),
        ImmutableList.of("0:0, 3:3, 0:3"),
        "1:1",
        0,
        WindingRule.POSITIVE,
        "full");
  }

  // Two triangles that touch along a common boundary.
  @Test
  public void testTouchingTriangles() {
    // The touch edges are considered to form a degenerate hole. Such holes are
    // removed by WindingRule.POSITIVE since they are not needed for computing
    // N-way unions. They are kept by WindingRule.ODD since they are needed in
    // order to compute N-way symmetric differences.
    expectWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("0:0, 0:8, 8:8", "0:0, 8:8, 8:0"),
        "1:1",
        1,
        WindingRule.POSITIVE,
        "0:0, 0:8, 8:8, 8:0");
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(),
        ImmutableList.of("0:0, 0:8, 8:8", "0:0, 8:8, 8:0"),
        "2:2",
        1,
        WindingRule.ODD,
        "0:0, 0:8, 8:8, 8:0",
        "0:0, 0:8, 8:8; 0:0, 8:8, 8:0");
  }

  // Like the test above, but the triangles only touch after snapping.
  @Test
  public void testTouchingTrianglesAfterSnapping() {
    // The snap function below rounds coordinates to the nearest degree.
    expectWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(0)),
        ImmutableList.of("0.1:0.2, 0:7.8, 7.6:8.2", "0.3:0.2, 8.1:7.8, 7.6:0.4"),
        "6:2",
        1,
        WindingRule.POSITIVE,
        "0:0, 0:8, 8:8, 8:0");
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(0)),
        ImmutableList.of("0.1:0.2, 0:7.8, 7.6:8.2", "0.3:0.2, 8.1:7.8, 7.6:0.4"),
        "2:6",
        1,
        WindingRule.ODD,
        "0:0, 0:8, 8:8, 8:0",
        "0:0, 0:8, 8:8; 0:0, 8:8, 8:0");
  }

  // This tests an N-way union of 5 overlapping squares forming a "staircase".
  @Test
  public void testUnionOfSquares() {
    // Five overlapping squares.
    ImmutableList<String> squares = ImmutableList.of(
            "0:0, 0:4, 4:4, 4:0",
            "1:1, 1:5, 5:5, 5:1",
            "2:2, 2:6, 6:6, 6:2",
            "3:3, 3:7, 7:7, 7:3",
            "4:4, 4:8, 8:8, 8:4");

    expectWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(1)),
        squares,
        "0.5:0.5",
        1,
        WindingRule.POSITIVE,
        "7:4, 7:3, 6:3, 6:2, 5:2, 5:1, 4:1, 4:0, 0:0, 0:4, "
            + "1:4, 1:5, 2:5, 2:6, 3:6, 3:7, 4:7, 4:8, 8:8, 8:4");

    // This computes the region overlapped by at least two squares.
    expectWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(1)),
        squares,
        "0.5:0.5",
        0,
        WindingRule.POSITIVE,
        "6:4, 6:3, 5:3, 5:2, 4:2, 4:1, 1:1, 1:4, 2:4, 2:5, 3:5, 3:6, 4:6, 4:7, 7:7, 7:4");

    // This computes the region overlapped by at least three squares.
    expectWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(1)),
        squares,
        "0.5:0.5",
        -1,
        WindingRule.POSITIVE,
        "5:4, 5:3, 4:3, 4:2, 2:2, 2:4, 3:4, 3:5, 4:5, 4:6, 6:6, 6:4");

    // This computes the region overlapped by at least four squares.
    expectWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(1)),
        squares,
        "0.5:0.5",
        -2,
        WindingRule.POSITIVE,
        "3:3, 3:4, 4:4, 4:3; 4:4, 4:5, 5:5, 5:4");

    // WindingRule.ODD yields a pattern reminiscent of a checkerboard.
    expectWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(1)),
        squares,
        "0.5:0.5",
        1,
        WindingRule.ODD,
        "4:1, 4:0, 0:0, 0:4, 1:4, 1:1; "
            + "4:3, 4:2, 2:2, 2:4, 3:4, 3:3; "
            + "1:4, 1:5, 2:5, 2:4; "
            + "5:4, 5:3, 4:3, 4:4; "
            + "5:2, 5:1, 4:1, 4:2; "
            + "2:5, 2:6, 3:6, 3:5; "
            + "6:3, 6:2, 5:2, 5:3; "
            + "3:6, 3:7, 4:7, 4:6; "
            + "3:4, 3:5, 4:5, 4:4; "
            + "7:4, 7:3, 6:3, 6:4; "
            + "4:7, 4:8, 8:8, 8:4, 7:4, 7:7; "
            + "4:5, 4:6, 6:6, 6:4, 5:4, 5:5");
  }

  // This tests that WindingRule.ODD can be used to compute the symmetric
  // difference even for input geometries with degeneracies, e.g. one geometry
  // has a degenerate hole or degenerate shell that the other does not.
  @Test
  public void testSymmetricDifferenceDegeneracies() {
    expectDegenerateWindingResult(
        new S2WindingOperation.Options(new IntLatLngSnapFunction(1)),
        ImmutableList.of(
            "0:0, 0:3, 3:3, 3:0",
            "1:1",
            "2:2",
            "4:4", // Geometry 1
            "0:0, 0:3, 3:3, 3:0",
            "1:1",
            "4:4",
            "5:5"), // Geometry 2
        "10:10",
        0,
        WindingRule.ODD,
        "",
        "2:2; 5:5");
  }

  @Test
  public void testSetGetSnapFunction() {
    S2WindingOperation.Options opts = new S2WindingOperation.Options();
    opts.setSnapFunction(new IdentitySnapFunction(S1Angle.degrees(1)));
    assertEquals(opts.snapFunction().snapRadius(), S1Angle.degrees(1));
  }
}
