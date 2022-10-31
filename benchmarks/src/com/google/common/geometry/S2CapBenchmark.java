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

import static com.google.common.geometry.S2.M_PI_4;
import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2Cell;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/** Benchmarks for S2Cap. */
public final class S2CapBenchmark {
  private S2CapBenchmark() {}

  // About 9 times the double-precision roundoff relative error.
  private static final double K_EPS = 1e-15;

  /** Benchmark state that provides constant S2Cells and S2Caps. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class ConstantState {
    private final S2Cell rootCell = S2Cell.fromFace(0);
    private final S2Cap emptyCap = S2Cap.empty();
    private final S2Cap fullCap = S2Cap.full();

    // This "bulging" cap has the same center as the root cell, intersects the cell's edges, but
    // does not contain any of its vertices.
    private final S2Cap bulgingCap =
        S2Cap.fromAxisAngle(rootCell.getCenter(), S1Angle.radians(M_PI_4 + K_EPS));

    /**
     * Measures the amount of time it takes to check mayIntersect for an empty cap and a root cell.
     * Equivalent to BM_MayIntersectWholeFaceEmptyCap in the C++ implementation.
     */
    @Benchmark
    public void mayIntersectWholeFaceEmptyCap(Blackhole bh) {
      bh.consume(emptyCap.mayIntersect(rootCell));
    }

    /**
     * Measures the amount of time it takes to check mayIntersect for a full cap and a root cell.
     * Equivalent to BM_MayIntersectWholeFaceFullCap in the C++ implementation.
     */
    @Benchmark
    public void mayIntersectWholeFaceFullCap(Blackhole bh) {
      bh.consume(fullCap.mayIntersect(rootCell));
    }

    /**
     * Measures the amount of time it takes to check mayIntersect for the "bulging" cap and the root
     * cell that intersects it. Equivalent to BM_MayIntersectBulging in the C++ implementation.
     */
    @Benchmark
    public void mayIntersectWholeFaceBulgingCap(Blackhole bh) {
      bh.consume(bulgingCap.mayIntersect(rootCell));
    }
  }
}
