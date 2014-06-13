package com.google.common.geometry;

import com.google.caliper.runner.CaliperMain;

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
    run(S2CellUnionBenchmark.class);
    run(S2EdgeUtilBenchmark.class);
    run(S2LatLngBenchmark.class);
    run(S2LoopBenchmark.LoopConstructor.class);
    run(S2LoopBenchmark.IsValid.class);
    run(S2LoopBenchmark.LoopContains.class);
    run(S2LoopBenchmark.ConstructAndContains.class);
    run(S2LoopBenchmark.CompareBoundaryContains.class);
    run(S2LoopBenchmark.CompareBoundaryCrosses.class);
    run(S2LoopBenchmark.CompareBoundaryDisjoint.class);
    run(S2LoopBenchmark.ContainsContains.class);
    run(S2LoopBenchmark.ContainsCrosses.class);
    run(S2LoopBenchmark.ContainsDisjoint.class);
    run(S2LoopBenchmark.IntersectsContains.class);
    run(S2LoopBenchmark.IntersectsCrosses.class);
    run(S2LoopBenchmark.IntersectsDisjoint.class);
    run(S2LoopBenchmark.ConstructAndContains.class);
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
    run(S2PolygonBenchmark.IntersectLoopGridWithCovering.class);
    run(S2PolygonBenchmark.UnionNestedFractalWithSelf.class);
    run(S2PolygonBenchmark.UnionLoopGridWithBound.class);
    run(S2ShapeIndexBenchmark.class);
  }

  static void run(Class<?> benchmarkClass) throws Exception {
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
