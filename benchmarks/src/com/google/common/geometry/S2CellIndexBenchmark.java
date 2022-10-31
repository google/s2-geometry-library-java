/*
 * Copyright 2019 Google Inc.
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

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2CellIndex;
import com.google.common.geometry.S2CellUnion;
import com.google.common.geometry.S2RegionCoverer;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
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

/** Benchmarks for S2CellIndex. */
public class S2CellIndexBenchmark {
  private S2CellIndexBenchmark() {}

  /**
   * Benchmark state which creates NUM_TARGET_CAPS target S2Caps in a query area, and stores
   * approximate coverings for each target cap in a list. Each cap has an approximate covering
   * consisting of numCellsPerCap cells.
   *
   * <p>The state also creates an S2CellIndex, containing many S2Cap coverings generated in the same
   * way as the target cap coverings. This is set up such that a random target cap will intersect
   * approximately (numCaps / numCapsPerArea) caps in the index.
   *
   * <p>The parameters cover a variety of values for numCaps, numCapsPerArea, and numCellsPerCap.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class IndexAndIntersectingTargetsState extends S2BenchmarkBaseState {

    public static enum SizeParams {
      CASE1(100, 1000, 16), // ~0.1 labels/query
      CASE2(100, 1000, 128), // ~0.1 labels/query
      CASE3(10000, 100000, 16), // ~0.1 labels/query
      CASE4(10000, 10000, 16), // ~1 label/query
      CASE5(10000, 1000, 16); // ~10 labels/query

      private final int numCaps;
      private final int numCapsPerArea;
      private final int numCellsPerCap;

      SizeParams(int numCaps, int numCapsPerArea, int numCellsPerCap) {
        this.numCaps = numCaps;
        this.numCapsPerArea = numCapsPerArea;
        this.numCellsPerCap = numCellsPerCap;
      }
    }

    private static final int NUM_TARGET_CAPS = 100;
    /** The approximate radius of S2Cap within which test caps are constructed. */
    private static final S1Angle kQueryCapRadius = kmToAngle(10);

    @Param
    SizeParams sizes;

    protected final S2CellIndex cellIndex = new S2CellIndex();
    protected final List<S2CellUnion> targets = new ArrayList<>();
    protected int targetIndex;

    /** Setup only needs to be done once, and can be reused for all iterations. */
    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      // Build a query cap.
      S2Cap queryCap = S2Cap.fromAxisAngle(data.getRandomPoint(), kQueryCapRadius);

      // Rebuild the cellIndex with the current SizeParams.
      cellIndex.clear();
      for (int i = 0; i < sizes.numCaps; ++i) {
        double areaFraction = data.nextDouble() / sizes.numCapsPerArea;
        cellIndex.add(getCapCovering(queryCap, areaFraction, sizes.numCellsPerCap), i);
      }
      cellIndex.build();

      // Rebuild the query targets with the current value of the sizes parameter.
      targets.clear();
      for (int i = 0; i < NUM_TARGET_CAPS; i++) {
        double areaFraction = data.nextDouble() / sizes.numCapsPerArea;
        targets.add(getCapCovering(queryCap, areaFraction, sizes.numCellsPerCap));
      }
      targetIndex = 0;
    }

    private S2CellUnion getCapCovering(S2Cap queryCap, double areaFraction, int numCellsPerCap) {
      S2Cap cap = S2Cap.fromAxisArea(data.samplePoint(queryCap), areaFraction * queryCap.area());
      return S2RegionCoverer.builder().setMaxCells(numCellsPerCap).build().getCovering(cap);
    }

    /**
     * (Method 1) Benchmark how long it takes to find the set of coverings in the index that
     * intersect a random query cap constructed in the same way as the indexed coverings. Iterates
     * through the targets, calling getIntersectingLabels() for each. The benchmark returns the
     * size() of the returned Labels object to ensure that results are computed.
     */
    @Benchmark
    public int getIntersectingLabelsForCapCovering() {
      int intersectingCount = cellIndex.getIntersectingLabels(targets.get(targetIndex)).size();
      targetIndex = (targetIndex + 1) % targets.size();
      return intersectingCount;
    }

    /**
     * (Method 2) Benchmark how long it takes to find the set of coverings in the index that
     * intersect a random query cap constructed in the same way as the indexed coverings. Iterates
     * through the targets, calling visitIntersectingCells() for each. The CellVisitor will visit
     * each intersecting (cellId, label) pair in the index, and tracks the set of unique labels seen
     * along the way. The benchmark returns the size of the resulting set to ensure that results are
     * computed.
     */
    @Benchmark
    public int getIntersectingLabelsByVisitingIntersectingCells() {
      HashSet<Integer> labelSet = new HashSet<>();
      cellIndex.visitIntersectingCells(targets.get(targetIndex), (cellId, label) -> {
        labelSet.add(label);
        return true;
      });
      targetIndex = (targetIndex + 1) % targets.size();
      return labelSet.size();
    }
  }
}
