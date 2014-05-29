package com.google.common.geometry;

import com.google.caliper.BeforeExperiment;
import com.google.caliper.Benchmark;
import com.google.common.collect.Lists;

import java.util.List;
import java.util.Random;

/** 
 * A benchmark for {@link S2CellId}. Caliper ignores returned values such as idSum and levelSum,
 * but these are computed and returned so that corresponding for loops don't get optimized away.
 */
public class S2CellIdBenchmark {
  private Random random;

  /**
   * Ensure that each benchmark is always run with the same random numbers.
   */
  @BeforeExperiment void setUp() {
    random = new Random(0);
  }

  /**
   * Test the amount of time it takes to convert from leaf cells to points.
   */
  @Benchmark double toPoint(int reps) {
    if (reps == 0) {
      // Avoids division by 0.
      return 0.0;
    }
    S2CellId begin = S2CellId.begin(S2CellId.MAX_LEVEL);
    S2CellId end = S2CellId.end(S2CellId.MAX_LEVEL);
    long delta = (end.id() - begin.id()) / reps;
    // Make delta's last bit equal to 0, so that all id's are leaf cells.
    delta &= ~0x1L;

    double xSum = 0;
    S2CellId id = begin;
    for (int r = reps; r > 0; --r) {
      xSum += id.toPoint().x;
      id = new S2CellId(id.id() + delta);
    }

    return xSum;
  }

  /**
   * Measure the cost of calling {@link S2CellId#fromPoint(S2Point)} on 20 points in a helical
   * spiral around the z-axis.
   */
  @Benchmark long fromPoint(int reps) {
    // Sample points used for benchmarking fromPoint().
    List<S2Point> points = Lists.newArrayList();
    // The sample points follow a spiral curve that completes one revolution around the z-axis
    // every 1/dt samples.  The z-coordinate increases from -4 to +4 over 'count' samples.
    int count = 20;
    S2Point p = new S2Point(1, 0, -4);
    double dz = (-2 * p.z) / count;
    double dt = 0.05;
    for (int r = count; r > 0; --r) {
      points.add(p);
      // Cheap revolution around the z-axis
      p = S2Point.add(p, new S2Point(-dt * p.y, dt * p.x, dz));
    }
    
    long idSum = 0;
    for (int r = reps; r > 0; --r) {
      for (S2Point pt : points) {
        idSum += S2CellId.fromPoint(pt).id();
      }
    }

    return idSum;
  }

  /**
   * Test the amount of time it takes to return a leaf cell given its cube face (range 0..5)
   * and i- and j-coordinates.
   */
  @Benchmark long fromFaceIJ(int reps) {
    // These can be the same for each iteration because fromFaceIJ's performance is not input
    // dependant.
    int face = random.nextInt(6);
    // i, j are random integers in [0, 2^30 - 1].
    int i = random.nextInt(1 << 30);
    int j = random.nextInt(1 << 30);

    long idSum = 0;
    for (int r = reps; r > 0; --r) {
      idSum += S2CellId.fromFaceIJ(face, i, j).id();
    }

    return idSum;
  }
    
  /**
   * Test the amount of time it takes to compute the level of a level 10 cell.
   */
  @Benchmark int level10(int reps) {
    return level(reps, 10);
  }
  
  /**
   * Test the amount of time it takes to compute the level of a level 20 cell.
   */
  @Benchmark int level20(int reps) {
    return level(reps, 20);
  }
  
  /**
   * Test the amount of time it takes to compute the level of a level 30 cell.
   */
  @Benchmark int level30(int reps) {
    return level(reps, 30);
  }

  /**
   * Compute the sum of the levels of 'reps' cells of a given level.
   */
  int level(int reps, int level) {
    // level() only depends on how many zero bits are at the end of the cell id, so we do not need
    // to test it on points with different degrees.
    S2CellId id = S2CellId.fromLatLng(S2LatLng.fromDegrees(5, 6)).parent(level);
    int levelSum = 0;
    for (int r = reps; r > 0; --r) {
      levelSum += id.level();
    }

    return levelSum;
  }
}