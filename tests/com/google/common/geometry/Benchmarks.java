package com.google.common.geometry;

import com.google.caliper.config.InvalidConfigurationException;
import com.google.caliper.runner.CaliperMain;
import com.google.caliper.runner.InvalidBenchmarkException;
import com.google.caliper.util.InvalidCommandException;

import java.io.PrintWriter;

/**
 * A class that allows us to call other Benchmark classes.  This is necessary in order for Caliper
 * to see static inner classes.  Currently we cannot do this with java_benchmarks since the
 * java_benchmarks rule only takes in files, not classes.
 */
public class Benchmarks {  
  public static void main(String[] args) throws Exception {
    run(S2CellIdBenchmark.class);
    run(S2CellIdBenchmark.Level.class);
    run(S2LoopBenchmark.class);
    run(S2LoopBenchmark.LoopConstructor.class);
    run(S2LoopBenchmark.IsValid.class);
    run(S2LoopBenchmark.LoopContains.class);
    run(S2LoopBenchmark.ConstructAndContains.class);
    run(S2LoopBenchmark.CompareBoundaryBenchmarks.class);
    run(S2LoopBenchmark.ContainsVsRadiusMeters.class);
    run(S2LoopBenchmark.ContainsVsGapEdgeMultiple.class);
    run(S2LoopBenchmark.IntersectsCrossesVsLogRadiusRatio.class);
    run(S2PolygonBenchmark.ConstructorSingleLoop.class);
    run(S2PolygonBenchmark.ConstructorLoopGrid.class);
    run(S2PolygonBenchmark.IsValidConcentricLoops.class);
    run(S2PolygonBenchmark.ContainsPointLoopGrid.class);
    run(S2PolygonBenchmark.ContainsPointNestedFractals.class);
    run(S2PolygonBenchmark.Covering.class);
    run(S2PolygonBenchmark.ContainsSelfLoopGrid.class);
    run(S2PolygonBenchmark.ContainsSelfNestedFractals.class);
    run(S2PolygonBenchmark.IntersectsCmplLoopGrid.class);
    run(S2PolygonBenchmark.IntersectsCmplNestedFractals.class);
    run(S2PolygonBenchmark.IntersectFractalWithCovering.class);
    run(S2PolygonBenchmark.UnionNestedFractalWithSelf.class);
    run(S2PolygonBenchmark.UnionLoopGridWithBound.class);
  }
  
  private static void run(Class<?> benchmarkClass)
      throws InvalidCommandException, InvalidBenchmarkException, InvalidConfigurationException {
    // Invoke exitlessMain so we can run additional experiments in this binary.
    CaliperMain.exitlessMain(new String[] {
        // Run only the runtime instrument.  The allocation instrument
        // disables Java optimizations, which makes this take too long.
        "-iruntime",
        // Pass in the name of the benchmark class.
        benchmarkClass.getName()},
        new PrintWriter(System.out, true), new PrintWriter(System.err, true));
  }
}