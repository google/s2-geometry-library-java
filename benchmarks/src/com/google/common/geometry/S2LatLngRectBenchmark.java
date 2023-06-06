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

import com.google.common.geometry.S2LatLng;
import com.google.common.geometry.S2LatLngRect;
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

/** Benchmarks for S2LatLngRect and its Builder. */
@CheckReturnValue
public class S2LatLngRectBenchmark {
  public S2LatLngRectBenchmark() {}

  /**
   * Benchmark state that provides an array of S2LatLngs and a corresponding array of S2Points,
   * placed "randomly" in a range of -30 degrees to +45 degrees for both lat and lng. This range is
   * arbitrary -- the exact range does not affect the runtime of the benchmarks. The random points
   * are the same for every iteration, as setup() resets the random number generator.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class S2ArrayOfPointsState extends S2BenchmarkBaseState {
    static final int ARRAY_SIZE = 128;
    private S2LatLng[] latlngs;
    private S2Point[] points;
    private int index;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      latlngs = new S2LatLng[ARRAY_SIZE];
      points = new S2Point[ARRAY_SIZE];
      for (index = 0; index < ARRAY_SIZE; index++) {
        S2LatLng ll = S2LatLng.fromDegrees(data.uniform(-30, 45), data.uniform(-30, 45));
        latlngs[index] = ll;
        points[index] = ll.toPoint();
      }
      index = 0;
    }

    /**
     * Measures the time to create an empty S2LatLngRect.Builder, add 8 S2LatLngs to it, and then
     * build and return the resulting S2LatLngRect.
     */
    @Benchmark
    public S2LatLngRect addS2LatLngsX8() {
      S2LatLngRect.Builder builder = new S2LatLngRect.Builder();
      for (int i = 0; i < 8; i++) {
        builder.addPoint(latlngs[index + i]);
      }
      index = (index + 8) % ARRAY_SIZE;
      return builder.build();
    }

    /**
     * Measures the time to create an empty S2LatLngRect.Builder, add 8 S2LatLngs to it, first
     * converting each to an S2Point, and then build and return the resulting S2LatLngRect.
     * Note that converting S2Latlngs to S2Points before adding them to the builder is a bad idea.
     * The addPoint(S2Point) call immediately converts them back to S2LatLngs, and these two
     * conversions make it about 8 times slower than adding S2LatLngs directly.
     */
    @Benchmark
    public S2LatLngRect addS2LatLngsAsPointsX8() {
      S2LatLngRect.Builder builder = new S2LatLngRect.Builder();
      for (int i = 0; i < 8; i++) {
        builder.addPoint(latlngs[index + i].toPoint());
      }
      index = (index + 8) % ARRAY_SIZE;
      return builder.build();
    }

    /**
     * Measures the time to create an empty S2LatLngRect.Builder, add 8 S2Points to it, and then
     * build and return the resulting S2LatLngRect.
     */
    @Benchmark
    public S2LatLngRect addS2PointsX8() {
      S2LatLngRect.Builder builder = new S2LatLngRect.Builder();
      for (int i = 0; i < 8; i++) {
        builder.addPoint(points[index + i]);
      }
      index = (index + 8) % ARRAY_SIZE;
      return builder.build();
    }
  }
}
