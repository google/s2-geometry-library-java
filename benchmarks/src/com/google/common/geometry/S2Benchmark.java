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
package com.google.common.geometry.benchmarks;

import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.TestDataGenerator;
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

/** Benchmarks for S2. */
public class S2Benchmark {
  public S2Benchmark() {}

  /**
   * Benchmark state that supports iteration through arrays of S2Points which are nearly parallel or
   * nearly antipodal to a single fixed but random point.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class AngleState extends S2BenchmarkBaseState {
    static final int ARRAY_SIZE = 128;

    // Size of the sample caps as the angle between the cap axis and cap boundary.
    S1Angle capSize = TestDataGenerator.metersToAngle(10);

    // A random point.
    private S2Point a;
    // A list of points sampled from a small cap around 'a'.
    private final S2Point[] nearParallel = new S2Point[ARRAY_SIZE];
    // A list of points sampled from a small cap antipodal to 'a'.
    private final S2Point[] nearAntipodal = new S2Point[ARRAY_SIZE];
    // An index into the lists above.
    private int index;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      a = data.getRandomPoint();
      S2Cap sampleCapNearParallel = S2Cap.fromAxisAngle(a, capSize);
      S2Cap sampleCapNearAntipodal = S2Cap.fromAxisAngle(a.neg(), capSize);

      for (index = 0; index < ARRAY_SIZE; index++) {
        nearParallel[index] = data.samplePoint(sampleCapNearParallel);
        nearAntipodal[index] = data.samplePoint(sampleCapNearAntipodal);
      }
      index = 0;
    }

    /**
     * Measures the time to compute the angle between 'a' and a point nearly parallel to 'a' using
     * {@link S2Point#angle(S2Point)}.
     */
    @Benchmark
    public double s2PointAngleNearParallel() {
      index = (index + 1) % ARRAY_SIZE;
      return a.angle(nearParallel[index]);
    }

    /**
     * Measures the time to compute the angle between 'a' and a point nearly parallel to 'a' using
     * {@link S2#stableAngle(S2Point, S2Point)}.
     */
    @Benchmark
    public double stableAngleNearParallel() {
      index = (index + 1) % ARRAY_SIZE;
      return S2.stableAngle(a, nearParallel[index]);
    }

    /**
     * Measures the time to compute the angle between 'a' and a point nearly antipodal to 'a' using
     * {@link S2Point#angle(S2Point)}.
     */
    @Benchmark
    public double s2PointAngleNearAntipodal() {
      index = (index + 1) % ARRAY_SIZE;
      return a.angle(nearAntipodal[index]);
    }

    /**
     * Measures the time to compute the angle between 'a' and a point nearly antipodal to 'a' using
     * {@link S2#stableAngle(S2Point, S2Point)}.
     */
    @Benchmark
    public double stableAngleNearAntipodal() {
      index = (index + 1) % ARRAY_SIZE;
      return S2.stableAngle(a, nearAntipodal[index]);
    }
  }
}
