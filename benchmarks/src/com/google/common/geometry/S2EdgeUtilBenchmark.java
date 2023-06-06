/*
 * Copyright 2021 Google Inc.
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
package com.google.common.geometry.benchmarks;

import static java.lang.Math.PI;
import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S1ChordAngle;
import com.google.common.geometry.S2EdgeUtil;
import com.google.common.geometry.S2Point;
import com.google.errorprone.annotations.CheckReturnValue;
import java.io.IOException;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/** Benchmarks for {@link S2EdgeUtil}. */
@CheckReturnValue
public class S2EdgeUtilBenchmark {
  private S2EdgeUtilBenchmark() {}

  /** Benchmarks the interpolation and crossing functions on a list of randomly generated points. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 10, time = 10, timeUnit = SECONDS)
  public static class BenchmarkCrossingState extends S2BenchmarkBaseState {
    // We want to avoid cache effects, so these arrays should be small enough to fit in L1 cache.
    // An S2Point is at most 40 bytes, so 250 will take about ~10 KiB.
    protected static final int NUM_POINTS = 250;
    protected static final int NUM_ANGLES = 120;

    // Randomly selected values between 0 and 1.
    protected double[] fractions = new double[NUM_ANGLES];
    // Randomly selected angles between 0 and PI/4.
    protected S1Angle[] angles = new S1Angle[NUM_ANGLES];
    protected S1ChordAngle[] chordAngles = new S1ChordAngle[NUM_ANGLES];

    // Randomly selected points.
    protected S2Point[] points = new S2Point[NUM_POINTS];
    // Precomputed S1Angle distances between successive points.
    protected S1Angle[] precomputedS1Angle = new S1Angle[NUM_POINTS];

    protected int pointIndex;
    protected int distanceIndex;
    protected S2Point a;
    protected S2Point b;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      for (int i = 0; i < NUM_POINTS; ++i) {
        points[i] = data.getRandomPoint();
      }
      for (int i = 0; i < NUM_POINTS; ++i) {
        S2Point a = points[i];
        S2Point b = points[(i + 1) % NUM_POINTS];
        precomputedS1Angle[i] = new S1Angle(a, b);
      }
      pointIndex = 0;

      for (int i = 0; i < NUM_ANGLES; ++i) {
        fractions[i] = data.uniform(0, 1);
        angles[i] = S1Angle.radians(data.uniform(0, PI / 4));
        chordAngles[i] = S1ChordAngle.fromS1Angle(angles[i]);
      }
      distanceIndex = 0;

      // Approximately 1/4th of points will cross the edge 'ab'.
      a = data.getRandomPoint();
      b = S2Point.neg(a).add(new S2Point(0.1, 0.1, 0.1)).normalize();
    }

    /**
     * Measures overhead due to iterating through and obtaining four S2Points with this benchmark
     * State. Subtract the time measured by this benchmark from the Crossing benchmarks below to get
     * a more accurate measure of their real cost.
     *
     * <p>For example, one test run found the following:
     * <ul>
     * <li> 37.787 +/- 0.499 ns/op for benchmarkFourPointOverhead.
     * <li>117.601 +/- 15.786 ns/op for edgeOrVertexCrossing, so true cost is 79 +/- 16 ns.
     * <li>113.949 +/- 1.507 ns/op for robustCrossing, so true cost is 76 +/- 2 ns.
     * </ul>
     */
    @Benchmark
    public int benchmarkFourPointOverhead(Blackhole bh) {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S2Point c = points[(pointIndex + 2) % NUM_POINTS];
      S2Point d = points[(pointIndex + 3) % NUM_POINTS];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      bh.consume(a);
      bh.consume(b);
      bh.consume(c);
      bh.consume(d);
      return pointIndex;
    }

    /** Benchmarks a single call to edgeOrVertexCrossing() with random points. */
    @Benchmark
    public boolean edgeOrVertexCrossing() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S2Point c = points[(pointIndex + 2) % NUM_POINTS];
      S2Point d = points[(pointIndex + 3) % NUM_POINTS];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      return S2EdgeUtil.edgeOrVertexCrossing(a, b, c, d);
    }

    /** Benchmarks a single call to robustCrossing() with random points. */
    @Benchmark
    public int robustCrossing() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S2Point c = points[(pointIndex + 2) % NUM_POINTS];
      S2Point d = points[(pointIndex + 3) % NUM_POINTS];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      return S2EdgeUtil.robustCrossing(a, b, c, d);
    }

    /**
     * Benchmarks setup of an EdgeCrosser for three points followed by 100
     * EdgeCrosser.robustCrossing() calls with random forth points, where approximately 1/4 of the
     * edges cross.
     */
    @Benchmark
    public void edgeCrosser100RobustCrossings(Blackhole bh) {
      S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a, b, points[0]);
      for (int r = 100; r > 0; --r) {
        S2Point d = points[r % NUM_POINTS];
        bh.consume(crosser.robustCrossing(d));
      }
    }

    /**
     * Measures overhead due to iterating through and obtaining two S2Points with this benchmark
     * State. Subtract the time measured by this benchmark from the interpolation benchmarks below
     * to get a more accurate measure of their real cost.
     */
    @Benchmark
    public int benchmarkTwoPointOverhead(Blackhole bh) {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      bh.consume(a);
      bh.consume(b);
      return pointIndex;
    }

    /** Benchmarks a single call to interpolate() with random points and distances. */
    @Benchmark
    public S2Point interpolate() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      double t = fractions[distanceIndex % NUM_ANGLES];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      distanceIndex = (distanceIndex + 1) % NUM_ANGLES;
      return S2EdgeUtil.interpolate(t, a, b);
    }

    /**
     * Benchmarks {@link S2EdgeUtil#getPointOnLine(S2Point, S2Point, S1Angle)} with random
     * points and distances.
     */
    @Benchmark
    public S2Point getPointOnLineS1Angle() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S1Angle angle = angles[distanceIndex % NUM_ANGLES];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      distanceIndex = (distanceIndex + 1) % NUM_ANGLES;
      return S2EdgeUtil.getPointOnLine(a, b, angle);
    }

    /**
     * Benchmarks {@link S2EdgeUtil#getPointOnLine(S2Point, S2Point, S1Angle)} with random
     * points and distances.
     */
    @Benchmark
    public S2Point getPointOnLineS1ChordAngle() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S1ChordAngle chordAngle = chordAngles[distanceIndex % NUM_ANGLES];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      distanceIndex = (distanceIndex + 1) % NUM_ANGLES;
      return S2EdgeUtil.getPointOnLine(a, b, chordAngle);
    }

    /**
     * Benchmarks {@link S2EdgeUtil#interpolateAtDistance(S1Angle, S2Point, S2Point)} with random
     * points and distances.
     */
    @Benchmark
    public S2Point interpolateAtDistance() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S1Angle angle = angles[distanceIndex % NUM_ANGLES];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      distanceIndex = (distanceIndex + 1) % NUM_ANGLES;
      return S2EdgeUtil.interpolateAtDistance(angle, a, b);
    }

    /**
     * Benchmarks {@link S2EdgeUtil#interpolateAtDistance(S1Angle, S2Point, S2Point, S1Angle)} with
     * precomputed S1Angle distance AB for random points and distances.
     */
    @Benchmark
    public S2Point interpolateAtDistanceWithPrecomputedS1AngleAB() {
      S2Point a = points[(pointIndex + 0) % NUM_POINTS];
      S2Point b = points[(pointIndex + 1) % NUM_POINTS];
      S1Angle ab = precomputedS1Angle[pointIndex % NUM_POINTS];
      S1Angle angle = angles[distanceIndex % NUM_ANGLES];
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      distanceIndex = (distanceIndex + 1) % NUM_ANGLES;
      return S2EdgeUtil.interpolateAtDistance(angle, a, b, ab);
    }
  }
}
