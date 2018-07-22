package com.google.common.geometry;

import com.google.caliper.BeforeExperiment;
import com.google.caliper.Benchmark;

import java.util.ArrayList;
import java.util.List;

/** 
 * A benchmark for {@link S2CellUnion}. Note that return values of benchmark functions are
 * ignored by Caliper.
 */
public class S2CellUnionBenchmark {
  private int numTestCells = 100;
  private S2CellUnion cellUnion;
  List<S2CellId> aCells, bCells;
  
  @BeforeExperiment void setUp() {
    GeometryTestCase testUtils = new GeometryTestCase();
    testUtils.setUp();
    cellUnion = new S2CellUnion();
    // Populate a list of cell ids to use for cell union initializing.
    aCells = new ArrayList<>(numTestCells);
    bCells = new ArrayList<>(numTestCells);
    for (int i = 0; i < numTestCells; ++i) {
      S2CellId a = testUtils.getRandomCellId(S2CellId.MAX_LEVEL);
      S2CellId b = testUtils.getRandomCellId(S2CellId.MAX_LEVEL);
      if (a.greaterThan(b)) {
        aCells.add(b);
        bCells.add(a);
      } else {
        aCells.add(a);
        bCells.add(b);
      }
    }
  }
  
  /** Benchmarks the speed of constructing a cell union from a range of cell ids. */
  @Benchmark void initFromBeginEnd(int reps) {
    for (int r = reps; r > 0; --r) {
      cellUnion.initFromBeginEnd(aCells.get(r % numTestCells), bCells.get(r % numTestCells));
    }
  }
}