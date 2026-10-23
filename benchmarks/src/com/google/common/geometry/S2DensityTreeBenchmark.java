/*
 * Copyright 2025 Google Inc.
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

import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2DensityTree;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder;
import com.google.common.geometry.S2DensityTree.Cell;
import com.google.common.geometry.S2DensityTree.CellVisitor;
import com.google.common.geometry.S2DensityTree.FeatureDensityWeigher;
import com.google.common.geometry.S2Shape;
import com.google.common.geometry.S2ShapeIndex;
import com.google.common.geometry.TestDataGenerator;
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

/** Benchmarks for S2DensityTree. */
public class S2DensityTreeBenchmark {

  private S2DensityTreeBenchmark() {}

  /**
   * Benchmark state for S2DensityTree benchmarks. Creates an S2ShapeIndex with a parameterized
   * number of regular loop shapes.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class DensityTreeState extends S2BenchmarkBaseState {

    // These sizes match the C++ BM_DecodeCellsByVisitor benchmark parameters.
    @Param({"10", "100", "1000", "10000", "100000"})
    int numShapes;

    @Param({"100", "1000", "10000", "100000"})
    int treeSizeBytes;

    int maxTreeLevel = 30;

    private S2ShapeIndex index;
    private S2DensityTree tree;
    private FeatureDensityWeigher weigher;

    @Override
    @Setup(Level.Trial)
    public void setup() {
      super.setup();
      index = new S2ShapeIndex();
      TestDataGenerator.RegularLoopShapeFactory factory =
          data.new RegularLoopShapeFactory();
      for (int i = 0; i < numShapes; i++) {
        S2Cap cap = S2Cap.fromAxisAngle(data.getRandomPoint(), TestDataGenerator.kmToAngle(1.0));
        // We use 5 vertices per loop to match the C++ benchmark.
        S2ShapeIndex shape = factory.getShape(cap, 5);
        index.add(shape.getShapes().get(0));
      }

      weigher = new FeatureDensityWeigher();
      for (S2Shape shape : index.getShapes()) {
        weigher.add(shape, 2);
      }

      tree = new BreadthFirstTreeBuilder(treeSizeBytes, maxTreeLevel).build(weigher);
    }

    /**
     * Measures the time to visit all cells in a density tree. Corresponds to
     * BM_DecodeCellsByVisitor in C++.
     */
    @Benchmark
    public int decodeCellsByVisitor() {
      CounterVisitor visitor = new CounterVisitor();
      tree.visitCells(visitor);
      return visitor.count;
    }

    /**
     * Measures the time to build a density tree from an index of polygons. Corresponds to
     * BM_InitToFeatureDensity in C++.
     */
    @Benchmark
    public S2DensityTree buildFeatureDensity() {
      // The weigher has an index associated with it after we called 'add', and that index is used
      // in 'applyAsLong'. We need to build a new tree based on that weigher.
      return new BreadthFirstTreeBuilder(treeSizeBytes, maxTreeLevel).build(weigher);
    }

    /** A simple visitor that counts the number of cells visited. */
    private static class CounterVisitor implements CellVisitor {
      int count = 0;
      @Override
      public Action visit(S2CellId cell, Cell node) {
        count++;
        return Action.ENTER_CELL;
      }
    }
  }
}
