/*
 * Copyright 2014 Google Inc.
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

import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S2Iterator;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2ShapeIndex;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.ArrayList;
import java.util.List;
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

/** Benchmarks for S2ShapeIndex. */
@CheckReturnValue
public class S2ShapeIndexBenchmark {
  private S2ShapeIndexBenchmark() {}

  /**
   * Benchmark state that provides a list of regular S2Loops, and also an S2Polygon made from them.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class ManyLoopsState extends S2BenchmarkBaseState {
    @Param({"1", "32", "64"})
    int numLoops;

    @Param({"4", "16", "256", "4096"})
    int numVertices;

    private final List<S2Loop> loops = new ArrayList<>();

    @Setup(Level.Iteration)
    @Override
    public void setup() {
      super.setup();
      loops.clear();
      for (int k = 0; k < numLoops; k++) {
        loops.add(S2Loop.makeRegularLoop(data.getRandomPoint(), kmToAngle(5), numVertices));
      }
    }

    /** Measures the time to construct an S2ShapeIndex with the provided S2Loops. */
    @Benchmark
    public S2Iterator<Cell> s2ShapeIndexInitialization() {
      S2ShapeIndex idx = new S2ShapeIndex();
      for (S2Loop loop : loops) {
        idx.add(loop);
      }
      return idx.iterator();
    }
  }
}
