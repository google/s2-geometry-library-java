package com.google.common.geometry;

import com.google.caliper.Benchmark;
import com.google.caliper.runner.CaliperMain;

/** A benchmark for {@link S2CellId#toPoint()}. */

public class S2CellIdBenchmark {
  /**
   * Test the amount of time it takes to convert from leaf cells to points. Caliper will choose 
   * values for reps.
   */
  @Benchmark double toPoint(int reps) {
    if (reps == 0) {
      return 0.0; // Avoids division by 0.
    }
    S2CellId begin = S2CellId.begin(S2CellId.MAX_LEVEL);
    S2CellId end = S2CellId.end(S2CellId.MAX_LEVEL);
    long delta = (end.id() - begin.id()) / reps;
    delta &= ~0x1L; // Make delta's last bit equal to 0, so that all id's are leaf cells.

    double xSum = 0;
    S2CellId id = begin;
    for (int i = reps; i > 0; --i) {
      xSum += id.toPoint().x;
      id = new S2CellId(id.id() + delta);
    }

    // Return value is ignored by Caliper, but we compute and return the sum of the x coordinates
    // so that the for loop doesn't get optimized away.
    return xSum;
  }

  public static void main(String[] args) {
    CaliperMain.main(S2CellIdBenchmark.class, args);
  }
}