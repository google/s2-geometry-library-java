/*
 * Copyright 2022 Google Inc.
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

import com.google.common.collect.ImmutableList;

/** Unit tests for {@link S2PolylineSimplifier}. */
public final class S2PolylineSimplifierTest extends GeometryTestCase {

  /**
   * Tests if extending an edge from 'srcPoint' to 'dstPoint' is possible given a set of
   * constraints: The edge must intersect discs of radius 'radiusDegrees' centered at
   * 'targetPoints', and must not intersect discs of radius 'radiusDegrees' centered at
   * 'avoidPoints'. The edge must also be on the left or right of each avoided disc as specified by
   * the 'discOnLeft' list, which must have an entry for each avoidPoint.
   *
   * <p>'expectedResult' indicates if extending an edge subject to these constraints is supposed to
   * be possible or not.
   */
  private void checkSimplify(
      String srcPoint,
      String dstPoint,
      String targetPoints,
      String avoidPoints,
      ImmutableList<Boolean> discOnLeft,
      double radiusDegrees,
      boolean expectedResult) {
    S1ChordAngle radius = S1ChordAngle.fromDegrees(radiusDegrees);
    S2PolylineSimplifier simplifier = new S2PolylineSimplifier();
    simplifier.init(S2TextFormat.makePointOrDie(srcPoint));
    for (S2Point p : S2TextFormat.parsePointsOrDie(targetPoints)) {
      simplifier.targetDisc(p, radius);
    }
    int i = 0;
    for (S2Point p : S2TextFormat.parsePointsOrDie(avoidPoints)) {
      simplifier.avoidDisc(p, radius, discOnLeft.get(i++));
    }
    assertEquals(
        "\nsrc = "
            + srcPoint
            + "\ndst = "
            + dstPoint
            + "\ntargetPoints = "
            + targetPoints
            + "\navoid = "
            + avoidPoints,
        expectedResult,
        simplifier.extend(S2TextFormat.makePointOrDie(dstPoint)));
  }

  /** Check that Init() can be called more than once. */
  public void testReuse() {
    S2PolylineSimplifier s = new S2PolylineSimplifier();
    S1ChordAngle radius = S1ChordAngle.fromDegrees(10);
    s.init(new S2Point(1, 0, 0));
    assertTrue(s.targetDisc(new S2Point(1, 1, 0).normalize(), radius));
    assertTrue(s.targetDisc(new S2Point(1, 1, 0.1).normalize(), radius));
    assertFalse(s.extend(new S2Point(1, 1, 0.4).normalize()));

    s.init(new S2Point(0, 1, 0));
    assertTrue(s.targetDisc(new S2Point(1, 1, 0.3).normalize(), radius));
    assertTrue(s.targetDisc(new S2Point(1, 1, 0.2).normalize(), radius));
    assertFalse(s.extend(new S2Point(1, 1, 0).normalize()));
  }

  /** Tests simplification with no constraints. */
  public void testNoConstraints() {
    // No constraints, dst == src.
    checkSimplify("0:1", "0:1", "", "", ImmutableList.of(), 0, true);

    // No constraints, dst != src.
    checkSimplify("0:1", "1:0", "", "", ImmutableList.of(), 0, true);

    // No constraints, (src, dst) longer than 90 degrees (not supported).
    checkSimplify("0:0", "0:91", "", "", ImmutableList.of(), 0, false);
  }

  /** Tests simplification of lines from 0:0 to 0:2 with a single target in various cases. */
  public void testTargetOnePoint() {
    // Three points on a straight line where the target is exactly between the source and
    // destination. In theory zero tolerance for the radius should work, but in practice there are
    // floating point errors.
    checkSimplify("0:0", "0:2", "0:1", "", ImmutableList.of(), 1e-10, true);

    // Three points where the middle point is too far away.
    checkSimplify("0:0", "0:2", "1:1", "", ImmutableList.of(), 0.9, false);

    // A target disc that contains the source vertex.
    checkSimplify("0:0", "0:2", "0:0.1", "", ImmutableList.of(), 1.0, true);

    // A target disc that contains the destination vertex.
    checkSimplify("0:0", "0:2", "0:2.1", "", ImmutableList.of(), 1.0, true);
  }

  /** Tests simplification of lines from 0:0 to 0:2 with a single "avoid" point in various cases. */
  public void testAvoidOnePoint() {
    // Three points on a straight line, attempting to avoid the middle point.
    checkSimplify("0:0", "0:2", "", "0:1", ImmutableList.of(true), 1e-10, false);

    // Three points where the middle point can be successfully avoided.
    checkSimplify("0:0", "0:2", "", "1:1", ImmutableList.of(true), 0.9, true);

    // Three points where the middle point is on the left, but where the client requires the point
    // to be on the right of the edge.
    checkSimplify("0:0", "0:2", "", "1:1", ImmutableList.of(false), 1e-10, false);

    // Check cases where the point to be avoided is behind the source vertex. In this situation
    // "discOnLeft" should not affect the result.
    checkSimplify("0:0", "0:2", "", "1:-1", ImmutableList.of(false), 1.4, true);
    checkSimplify("0:0", "0:2", "", "1:-1", ImmutableList.of(true), 1.4, true);
    checkSimplify("0:0", "0:2", "", "-1:-1", ImmutableList.of(false), 1.4, true);
    checkSimplify("0:0", "0:2", "", "-1:-1", ImmutableList.of(true), 1.4, true);
  }

  /**
   * Tests that when several discs are avoided but none are targeted, the range of acceptable edge
   * directions is not represented a single interval. (The output edge can proceed through any gap
   * between the discs as long as the "discOnLeft" criteria are satisfied.)
   */
  public void testAvoidSeveralPoints() {
    // This test involves 3 very small discs spaced 120 degrees apart around the source vertex,
    // where "discOnLeft" is true for all discs. This means that each disc blocks the 90 degrees of
    // potential output directions just to its left, leave 3 gaps measuring about 30 degrees each.
    for (String dst : ImmutableList.of("0:2", "1.732:-1", "-1.732:-1")) {
      checkSimplify(
          "0:0",
          dst,
          "",
          "0.01:2, 1.732:-1.01, -1.732:-0.99",
          ImmutableList.of(true, true, true),
          0.00001,
          true);

      // Also test that directions prohibited by "discOnLeft" are avoided.
      checkSimplify(
          "0:0",
          dst,
          "",
          "0.01:2, 1.732:-1.01, -1.732:-0.99",
          ImmutableList.of(false, false, false),
          0.00001,
          false);
    }
  }

  /** Tests constraints on both targeting and avoiding discs. */
  public void testTargetAndAvoid() {
    // Target several points that are separated from the proposed edge by about 0.7 degrees, and
    // avoid several points that are separated from the proposed edge by about 1.4 degrees.
    checkSimplify(
        "0:0",
        "10:10",
        "2:3, 4:3, 7:8",
        "4:2, 7:5, 7:9",
        ImmutableList.of(true, true, false),
        1.0,
        true);

    // The same example, but one point to be targeted is 1.4 degrees away.
    checkSimplify(
        "0:0",
        "10:10",
        "2:3, 4:6, 7:8",
        "4:2, 7:5, 7:9",
        ImmutableList.of(true, true, false),
        1.0,
        false);

    // The same example, but one point to be avoided is 0.7 degrees away.
    checkSimplify(
        "0:0",
        "10:10",
        "2:3, 4:3, 7:8",
        "4:2, 6:5, 7:9",
        ImmutableList.of(true, true, false),
        1.0,
        false);
  }

  /** Tests that targeting and avoiding discs has high enough precision. */
  public void testPrecision() {
    // This is a rough upper bound on both the error in constructing the disc locations (i.e.,
    // S2.getPointOnLine, etc.), and also on the padding that S2PolylineSimplifier
    // uses to ensure that its results are conservative (i.e., the error calculated by
    // getSemiwidth).
    final S1Angle kMaxError = S1Angle.radians(25 * DBL_EPSILON);

    // We repeatedly generate a random edge. We then target several discs that barely overlap the
    // edge, and avoid several discs that barely miss the edge. About half the time, we choose one
    // disc and make it slightly too large or too small so
    // that targeting fails.
    final int kIters = 1000; // Also passes with 1 million iterations on Java.
    S2PolylineSimplifier simplifier = new S2PolylineSimplifier();
    for (int iter = 0; iter < kIters; ++iter) {
      data.setSeed(iter + 1); // Easier to reproduce a specific case.
      S2Point src = data.getRandomPoint();
      simplifier.init(src);
      S2Point dst =
          S2EdgeUtil.getPointOnLine(
              src, data.getRandomPoint(), S1Angle.radians(data.uniform(0, 1)));
      S2Point n = S2.robustCrossProd(src, dst).normalize();

      // If badDisc >= 0, then we make targeting fail for that disc.
      final int kNumDiscs = 5;
      int badDisc = data.uniform(2 * kNumDiscs) - kNumDiscs;
      for (int i = 0; i < kNumDiscs; ++i) {
        System.err.println("testPrecision iteration " + iter + " disc " + i);
        // The center of the disc projects to a point that is the given fraction "f" along the edge
        // (src, dst). If f < 0, the center is located behind "src" (in order to test this case).
        double f = data.uniform(-0.5, 1.0);
        S2Point a = src.mul(1 - f).add(dst.mul(f)).normalize();
        S1Angle r = S1Angle.radians(data.uniform(0, 1));
        boolean onLeft = data.oneIn(2);
        S2Point x = S2EdgeUtil.getPointOnLine(a, onLeft ? n : n.neg(), r);
        // If the disc is behind "src", adjust its radius so that it just touches "src" rather than
        // just touching the line through (src, dst).
        if (f < 0) {
          r = new S1Angle(src, x);
        }
        // We grow the radius slightly if we want to target the disc and shrink it otherwise,
        // *unless* we want targeting to fail for this disc, in which case these actions are
        // reversed.
        boolean avoid = data.oneIn(2);
        boolean growRadius = (avoid == (i == badDisc));
        S1ChordAngle radius =
            S1ChordAngle.fromS1Angle(growRadius ? r.add(kMaxError) : r.sub(kMaxError));
        if (avoid) {
          simplifier.avoidDisc(x, radius, onLeft);
        } else {
          simplifier.targetDisc(x, radius);
        }
      }
      // The result is true iff all the disc constraints were satisfiable.
      assertEquals(badDisc < 0, simplifier.extend(dst));
    }
  }
}
