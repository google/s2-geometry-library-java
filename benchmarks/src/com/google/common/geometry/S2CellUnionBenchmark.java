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

import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2CellUnion;
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

/** Benchmarks for {@link S2CellUnion}. */
public class S2CellUnionBenchmark {
  private S2CellUnionBenchmark() {}

  /**
   * Benchmark state that provides two lists of random leaf cell ids, aCells and bCells, such that
   * for each index position i, aCells[i] < bCells[i].
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class TwoListsOfLeafCellIdsState extends S2BenchmarkBaseState {
    private static final int NUM_TEST_CELLS = 100;

    protected S2CellUnion cellUnion;
    protected List<S2CellId> aCells;
    protected List<S2CellId> bCells;
    protected int rep;

    @Setup(Level.Trial)
    public void setUp() throws IOException {
      super.setup();
      cellUnion = new S2CellUnion();
      // Populate two lists of cell ids to use for cell union initializing.
      aCells = new ArrayList<>(NUM_TEST_CELLS);
      bCells = new ArrayList<>(NUM_TEST_CELLS);
      for (int i = 0; i < NUM_TEST_CELLS; ++i) {
        S2CellId a = data.getRandomCellId(S2CellId.MAX_LEVEL);
        S2CellId b = data.getRandomCellId(S2CellId.MAX_LEVEL);
        if (a.greaterThan(b)) {
          aCells.add(b);
          bCells.add(a);
        } else {
          aCells.add(a);
          bCells.add(b);
        }
      }
      rep = 0;
    }

    /** Measures the time required to construct a cell union from a range of leaf cell ids. */
    @Benchmark
    public int initFromBeginEnd() {
      cellUnion.initFromBeginEnd(
          aCells.get(rep % NUM_TEST_CELLS),
          bCells.get(rep % NUM_TEST_CELLS));
      rep++;
      return cellUnion.size();
    }
  }
}
