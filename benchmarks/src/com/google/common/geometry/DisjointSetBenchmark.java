/*
 * Copyright 2024 Google Inc.
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

import com.google.common.geometry.primitives.DisjointSet;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

/**
 * Benchmark for {@link DisjointSet}.
 */
public final class DisjointSetBenchmark {

  private DisjointSetBenchmark() { }

  /** Benchmark state that sets up a disjoint set with consecutive Integer values. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class ConsecutiveIntegerState extends S2BenchmarkBaseState {
    @Param({"64", "512", "4096", "32768", "262144"})
    int numElements;

    private final DisjointSet<Integer> set = new DisjointSet<>(numElements);

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
    }

    /**
     * Measures the amount of time it takes to add numElements, union consecutive sets together, and
     * then find the root of the first element. Equivalent to BM_FindRoot in the C++ implementation.
     */
    @Benchmark
    public int findRoot() {
      set.clear();
      for (int i = 0; i < numElements; i++) {
        set.add(i);
      }
      for (int i = 0; i < numElements; i++) {
        set.union(i, i + 1);
      }

      return set.findRoot(1);
    }
  }
}
