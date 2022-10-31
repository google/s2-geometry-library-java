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

import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S2EdgeUtil;
import com.google.common.geometry.S2Point;
import com.google.errorprone.annotations.CheckReturnValue;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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

  /** Benchmarks the crossing functions on a list of randomly generated points. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 15, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 15, timeUnit = SECONDS)
  public static class BenchmarkCrossingState extends S2BenchmarkBaseState {
    protected static final int NUM_POINTS = 400;

    protected List<S2Point> points;
    protected int pointIndex;
    protected S2Point a;
    protected S2Point b;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      // We want to avoid cache effects, so numPoints should be small enough so that the points can
      // be in L1 cache.  The size of an S2Point is at most 40 bytes, so 400 will only take at most
      // ~16 KiB of 64 KiB of L1 cache.
      points = new ArrayList<>(NUM_POINTS);
      for (int i = 0; i < NUM_POINTS; ++i) {
        points.add(data.getRandomPoint());
      }
      pointIndex = 0;

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
     * <li> 37.787 +/- 0.499 ns/op for benchmarkOverhead.
     * <li>117.601 +/- 15.786 ns/op for edgeOrVertexCrossing, so true cost is 79 +/- 16 ns.
     * <li>113.949 +/- 1.507 ns/op for robustCrossing, so true cost is 76 +/- 2 ns.
     * </ul>
     */
    @Benchmark
    public int benchmarkOverhead(Blackhole bh) {
      S2Point a = points.get((pointIndex + 0) % NUM_POINTS);
      S2Point b = points.get((pointIndex + 1) % NUM_POINTS);
      S2Point c = points.get((pointIndex + 2) % NUM_POINTS);
      S2Point d = points.get((pointIndex + 3) % NUM_POINTS);
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
      S2Point a = points.get((pointIndex + 0) % NUM_POINTS);
      S2Point b = points.get((pointIndex + 1) % NUM_POINTS);
      S2Point c = points.get((pointIndex + 2) % NUM_POINTS);
      S2Point d = points.get((pointIndex + 3) % NUM_POINTS);
      pointIndex = (pointIndex + 1) % NUM_POINTS;
      return S2EdgeUtil.edgeOrVertexCrossing(a, b, c, d);
    }

    /** Benchmarks a single call to robustCrossing() with random points. */
    @Benchmark
    public int robustCrossing() {
      S2Point a = points.get((pointIndex + 0) % NUM_POINTS);
      S2Point b = points.get((pointIndex + 1) % NUM_POINTS);
      S2Point c = points.get((pointIndex + 2) % NUM_POINTS);
      S2Point d = points.get((pointIndex + 3) % NUM_POINTS);
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
      S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a, b, points.get(0));
      for (int r = 100; r > 0; --r) {
        S2Point d = points.get(r % NUM_POINTS);
        bh.consume(crosser.robustCrossing(d));
      }
    }
  }
}
