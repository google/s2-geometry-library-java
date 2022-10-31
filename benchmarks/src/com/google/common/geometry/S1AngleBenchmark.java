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

/**
 * Benchmarks for S1Angle. These verify that the conversion to E5/E6/E7 is not much slower than the
 * conversion from E5/E6/E7. (Float-to-integer conversions can be quite slow on some platforms.) We
 * only check the times for E6; the times for E5/E7 should be similar.
 *
 * <p>Inspect the results to make sure that the To/From E6 times are not much different. The
 * difference factor slightly less than 2 on an x86_64.
 */
public final class S1AngleBenchmark {
  private S1AngleBenchmark() {}

  /** Benchmark state that supports iterating through lists of angles in radians and e6. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class RadiansToAndFromE6State {
    private static final int NUM_ANGLES = 1 << 5;
    private static final int MASK_ANGLES = NUM_ANGLES - 1;

    // Angles in radians, for converting to e6.
    private double[] angles;

    // Angles in e6 (microdegrees), for converting to radians.
    private int[] e6s;

    protected int index;

    /** Fills the angles and e6s arrays with NUM_ANGLES evenly spaced angles, from -PI to PI. */
    @Setup(Level.Trial)
    public void setup() {
      angles = new double[NUM_ANGLES];
      e6s = new int[NUM_ANGLES];
      double delta = (2 * PI) / NUM_ANGLES;
      double angle = -PI;
      for (int i = 0; i < NUM_ANGLES; i++) {
        angles[i] = angle;
        e6s[i] = S1Angle.radians(angle).e6();
        angle += delta;
      }
      index = 0;
    }

    /**
     * Measures the amount of time it takes to convert e6 to radians. To reduce the impact of
     * benchmark overhead, we do 8 operations per invocation. Equivalent to BM_E6ToRadians in the
     * C++ implementation.
     */
    @Benchmark
    public void e6ToRadiansX8(Blackhole bh) {
      bh.consume(S1Angle.e6(e6s[index + 0]).radians());
      bh.consume(S1Angle.e6(e6s[index + 1]).radians());
      bh.consume(S1Angle.e6(e6s[index + 2]).radians());
      bh.consume(S1Angle.e6(e6s[index + 3]).radians());
      bh.consume(S1Angle.e6(e6s[index + 4]).radians());
      bh.consume(S1Angle.e6(e6s[index + 5]).radians());
      bh.consume(S1Angle.e6(e6s[index + 6]).radians());
      bh.consume(S1Angle.e6(e6s[index + 7]).radians());

      index = (index + 8) & MASK_ANGLES;
    }

    /**
     * Measures the amount of time it takes to convert from radians to e6. To reduce the impact of
     * benchmark overhead, we do 8 operations per invocation. Equivalent to BM_RadiansToE6 in the
     * C++ implementation.
     */
    @Benchmark
    public void radiansToE6X8(Blackhole bh) {
      bh.consume(S1Angle.radians(angles[index + 0]).e6());
      bh.consume(S1Angle.radians(angles[index + 1]).e6());
      bh.consume(S1Angle.radians(angles[index + 2]).e6());
      bh.consume(S1Angle.radians(angles[index + 3]).e6());
      bh.consume(S1Angle.radians(angles[index + 4]).e6());
      bh.consume(S1Angle.radians(angles[index + 5]).e6());
      bh.consume(S1Angle.radians(angles[index + 6]).e6());
      bh.consume(S1Angle.radians(angles[index + 7]).e6());

      index = (index + 8) & MASK_ANGLES;
    }
  }
}
