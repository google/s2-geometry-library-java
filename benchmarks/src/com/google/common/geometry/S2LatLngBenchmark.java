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

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2LatLng;
import com.google.common.geometry.S2Point;
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

/** Benchmarks for S2LatLng. */
public class S2LatLngBenchmark {
  public S2LatLngBenchmark() {}

  /**
   * Benchmark state that supports iteration through an array of S2LatLngs in a line at 56 degrees
   * longitude, spaced 1e-7 of a degree apart, which is less than one centimeter at that longitude.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class S2LatLngState extends S2BenchmarkBaseState {
    static final int ARRAY_SIZE = 128;
    private S2LatLng ll;
    private S2LatLng x;
    private S2LatLng[] array;
    private int index;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      ll = new S2LatLng(S1Angle.e7(0x150bc888), S1Angle.e7(0x5099d63f));
      x = S2LatLng.fromDegrees(25.0, -78.0);
      array = new S2LatLng[ARRAY_SIZE];
      for (index = 0; index < ARRAY_SIZE; index++) {
        array[index] = new S2LatLng(S1Angle.e7(index), S1Angle.degrees(56.0));
      }
      index = 0;
    }

    /**
     * Measures the time to convert an S2LatLng to an S2Point. Uses the same S2LatLng for each
     * invocation. Compare to toPointMany below, which uses a different S2LatLng for each
     * invocation.
     */
    @Benchmark
    public double toPointSingle() {
      return ll.toPoint().getX();
    }

    /**
     * Measures the time to convert an S2LatLng to an S2Point. Iterates through the array to use a
     * different S2LatLng for each invocation. The idea here is to compare to the results of
     * toPointSingle just above to determine how much of a difference it makes to use a range of
     * different values.
     */
    @Benchmark
    public double toPointMany() {
      index = (index + 1) % ARRAY_SIZE;
      return array[index].toPoint().getX();
    }

    /**
     * Measures the time to get the distance between two S2LatLngs. Repeatedly iterates through the
     * array of S2LatLngs, so each invocation is different.
     */
    @Benchmark
    public S1Angle getDistance() {
      index = (index + 1) % ARRAY_SIZE;
      return x.getDistance(array[index]);
    }
  }

  /** Benchmark state that provides a single S2Point which is always the same. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class S2PointState extends S2BenchmarkBaseState {
    protected S2Point point;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      point = data.getRandomPoint();
    }

    /**
     * Measures the time to compute the latitude of a random S2Point, which is the same for each
     * invocation.
     */
    @Benchmark
    public S1Angle staticLatitudeFromPoint() {
      return S2LatLng.latitude(point);
    }

    /**
     * Measures the time to compute the longitude of a random S2Point, which is the same for each
     * invocation.
     */
    @Benchmark
    public S1Angle staticLongitudeFromPoint() {
      return S2LatLng.longitude(point);
    }
  }
}
