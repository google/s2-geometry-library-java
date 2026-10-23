// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package com.google.common.geometry.benchmarks;

import static com.google.common.geometry.S2Loop.makeRegularLoop;
import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S1ChordAngle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2Cell;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2CellIndex;
import com.google.common.geometry.S2ClosestCellQuery;
import com.google.common.geometry.S2LaxPolygonShape;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2RegionCoverer;
import com.google.common.geometry.S2ShapeIndex;
import com.google.common.geometry.TestDataGenerator;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.ArrayList;
import java.util.List;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

@State(Scope.Thread)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(MICROSECONDS)
@Warmup(iterations = 1, time = 1, timeUnit = SECONDS)
@Measurement(iterations = 1, time = 1, timeUnit = SECONDS)
@Fork(1)
@CheckReturnValue
public class S2ClosestCellQueryBenchmark {
  /** The type of index to use. */
  public enum IndexFactory {
    POINT_CLOUD, COARSE_CAPS
  }

  /** The type of target to use. */
  public enum TargetType {
    POINT, EDGE, CELL, COARSE_CELL_UNION, FINE_CELL_UNION, COARSE_SHAPE_INDEX, FINE_SHAPE_INDEX
  }

  @State(Scope.Thread)
  public static class BaseState {
    public IndexFactory indexFactory;
    public TargetType targetType;
    public int numCells;
    public double maxDistFraction;

    public S2ClosestCellQuery query;
    public List<S2ClosestCellQuery.Target> targets = new ArrayList<>();
    public int targetIndex = 0;
    
    public void init(TestDataGenerator data) {
      S1Angle indexRadius = kmToAngle(100);
      S2Cap indexCap = S2Cap.fromAxisAngle(data.getRandomPoint(), indexRadius);
      S2CellIndex index = new S2CellIndex();
      
      if (indexFactory == IndexFactory.POINT_CLOUD) {
        for (int i = 0; i < numCells; i++) {
          index.add(S2CellId.fromPoint(data.samplePoint(indexCap)), i);
        }
      } else {
        int maxCellsPerCap = 8;
        double capDensity = 0.1;
        int numCaps = (numCells - 1) / maxCellsPerCap + 1;
        double maxArea = indexCap.area() * capDensity / numCaps;
        for (int i = 0; i < numCaps; ++i) {
          S2Cap cap = S2Cap.fromAxisArea(data.samplePoint(indexCap), data.uniform(0.0, maxArea));
          S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(maxCellsPerCap).build();
          index.add(coverer.getCovering(cap), i);
        }
      }
      index.build();
      query = new S2ClosestCellQuery(index);

      for (int i = 0; i < 100; i++) {
        targets.add(createTarget(data, indexCap, indexRadius));
      }
    }

    private S2ClosestCellQuery.Target createTarget(TestDataGenerator data, S2Cap indexCap, S1Angle indexRadius) {
      S1Angle targetDist = S1Angle.radians(data.uniform(0.0, 2.0 * indexRadius.radians()));
      S2Point targetCenter = data.samplePoint(S2Cap.fromAxisAngle(indexCap.axis(), targetDist));

      switch (targetType) {
        case POINT -> {
          return new S2ClosestCellQuery.PointTarget(targetCenter);
        }
        case EDGE -> {
          S1Angle edgeLen = S1Angle.radians(data.uniform(0.0, 0.1 * indexRadius.radians()));
          S2Point a = data.samplePoint(S2Cap.fromAxisAngle(targetCenter, edgeLen));
          S2Point b = data.samplePoint(S2Cap.fromAxisAngle(targetCenter, edgeLen));
          return new S2ClosestCellQuery.EdgeTarget(a, b);
        }
        case CELL -> {
          S1Angle cellRadius = S1Angle.radians(data.uniform(0.0, 0.01 * indexRadius.radians()));
          // Approximate the level.
          int level = Math.min(30, Math.max(0, (int) -Math.log(cellRadius.radians())));
          S2CellId cellId = S2CellId.fromPoint(targetCenter).parent(level);
          return new S2ClosestCellQuery.CellTarget(new S2Cell(cellId));
        }
        case COARSE_CELL_UNION, FINE_CELL_UNION -> {
          int numTargetCells = targetType == TargetType.COARSE_CELL_UNION ? 8 : 100;
          S1Angle capRadius = S1Angle.radians(data.uniform(0.0, 0.01 * indexRadius.radians()));
          S2Cap cap = S2Cap.fromAxisAngle(targetCenter, capRadius);
          S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(numTargetCells).build();
          return new S2ClosestCellQuery.CellUnionTarget(coverer.getCovering(cap));
        }
        case COARSE_SHAPE_INDEX, FINE_SHAPE_INDEX -> {
          int numTargetEdges = targetType == TargetType.COARSE_SHAPE_INDEX ? 12 : 768;
          S1Angle capRadius = S1Angle.radians(data.uniform(0.0, 0.01 * indexRadius.radians()));
          S2Cap cap = S2Cap.fromAxisAngle(targetCenter, capRadius);
          S2ShapeIndex targetIndex = new S2ShapeIndex();
          targetIndex.add(
              S2LaxPolygonShape.fromLoop(
                  makeRegularLoop(cap.axis(), cap.angle(), numTargetEdges).vertices()));
          return new S2ClosestCellQuery.ShapeIndexTarget(targetIndex);
        }
      }
      throw new IllegalArgumentException();
    }
  }

  @State(Scope.Thread)
  public static class FindClosestState extends BaseState {
    @Param({"12", "768", "12288", "49152"})
    public int numCellsParam;

    @Param({"POINT_CLOUD", "COARSE_CAPS"})
    public IndexFactory indexFactoryParam;

    @Param({"POINT", "EDGE", "CELL", "COARSE_CELL_UNION", "FINE_CELL_UNION", "COARSE_SHAPE_INDEX", "FINE_SHAPE_INDEX"})
    public TargetType targetTypeParam;

    @Setup(Level.Trial)
    public void setup() {
      this.numCells = numCellsParam;
      this.indexFactory = indexFactoryParam;
      this.targetType = targetTypeParam;
      this.maxDistFraction = -1.0;
      init(new TestDataGenerator());
    }
  }

  @State(Scope.Thread)
  public static class DistanceLessState extends BaseState {
    @Param({"12,-3", "768,-3", "12288,-3", "49152,-6", "49152,-3", "49152,0", "49152,3"})
    public String numCellsAndLog2Results;

    @Param({"POINT_CLOUD", "COARSE_CAPS"})
    public IndexFactory indexFactoryParam;

    @Param({"POINT", "EDGE", "CELL", "COARSE_CELL_UNION", "FINE_CELL_UNION", "COARSE_SHAPE_INDEX", "FINE_SHAPE_INDEX"})
    public TargetType targetTypeParam;

    public int log2NumResults;

    @Setup(Level.Trial)
    public void setup() {
      String[] parts = numCellsAndLog2Results.split(",");
      this.numCells = Integer.parseInt(parts[0]);
      this.log2NumResults = Integer.parseInt(parts[1]);
      this.indexFactory = indexFactoryParam;
      this.targetType = targetTypeParam;
      this.maxDistFraction = Math.sqrt(2 * Math.pow(2.0, log2NumResults) / numCells);
      init(new TestDataGenerator());
    }
  }

  @Benchmark
  public int findClosest(FindClosestState state) {
    S2ClosestCellQuery.Target target = state.targets.get(state.targetIndex++ % state.targets.size());
    state.query.setMaxResults(1);
    state.query.setMaxDistance(S1Angle.INFINITY);
    S2ClosestCellQuery.Result result = state.query.findClosestCell(target);
    return result.label();
  }

  @Benchmark
  public boolean isDistanceLess(DistanceLessState state) {
    S2ClosestCellQuery.Target target = state.targets.get(state.targetIndex++ % state.targets.size());
    S1Angle distLimit = S1Angle.radians(state.maxDistFraction * kmToAngle(100).radians());
    return state.query.isDistanceLess(target, S1ChordAngle.fromS1Angle(distLimit));
  }

  @Benchmark
  public int findCellsNearby(DistanceLessState state) {
    S2ClosestCellQuery.Target target = state.targets.get(state.targetIndex++ % state.targets.size());
    S1Angle distLimit = S1Angle.radians(state.maxDistFraction * kmToAngle(100).radians());
    state.query.setMaxResults(Integer.MAX_VALUE);
    state.query.setMaxDistance(distLimit);
    List<S2ClosestCellQuery.Result> results = state.query.findClosestCells(target);
    return results.size();
  }
}