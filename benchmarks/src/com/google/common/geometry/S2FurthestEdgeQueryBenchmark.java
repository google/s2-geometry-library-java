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

import com.google.common.base.Preconditions;
import com.google.common.geometry.S1ChordAngle;
import com.google.common.geometry.S2BestEdgesQueryBase.Result;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2FurthestEdgeQuery;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.benchmarks.BestEdgesBenchmarkHarness.BestEdgeQueryWithTargets;
import com.google.common.geometry.benchmarks.BestEdgesBenchmarkHarness.BestEdgesBenchmarkBaseState;
import com.google.common.geometry.benchmarks.BestEdgesBenchmarkHarness.FactoryType;
import com.google.common.geometry.benchmarks.BestEdgesBenchmarkHarness.TargetType;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.infra.BenchmarkParams;

/** Benchmarks for {@link S2FurthestEdgeQuery}. */
public class S2FurthestEdgeQueryBenchmark {
  private S2FurthestEdgeQueryBenchmark() {}

  /**
   * A container for an index, a furthest edge query on that index, and a set of targets to test
   * with.
   */
  public static class FurthestEdgeQueryWithTargets extends BestEdgeQueryWithTargets {
    public S2FurthestEdgeQuery<S1ChordAngle> query;
    public ArrayList<S2FurthestEdgeQuery.Target<S1ChordAngle>> targets;

  }

  public static class FurthestEdgeBenchmarkBaseState extends BestEdgesBenchmarkBaseState {
    protected ArrayList<FurthestEdgeQueryWithTargets> furthestEdgeQueriesAndTargets;
    protected S2FurthestEdgeQuery<S1ChordAngle> currentQuery;
    protected S2FurthestEdgeQuery.Target<S1ChordAngle> currentTarget;

    @Override
    void setCurrentQuery() {
      currentQuery = furthestEdgeQueriesAndTargets.get(queryIndex).query;
    }

    @Override
    void setCurrentTarget() {
      currentTarget = furthestEdgeQueriesAndTargets.get(queryIndex).targets.get(targetIndex);
    }

    @Override
    protected void createQueriesAndTargets() throws IOException {
      // Build the furthest edge query options.
      S2FurthestEdgeQuery.Builder builder =
          S2FurthestEdgeQuery.builder()
              .setUseBruteForce(useBruteForce)
              .setIncludeInteriors(includeInteriors)
              .setMaxResults(maxResults);

      if (distanceFraction > 0) {
        builder.setMinDistance(indexRadius.mul(distanceFraction));
      }
      if (maxErrorFraction > 0) {
        builder.setMaxError(indexRadius.mul(maxErrorFraction));
      }

      /*
       * Create and fill the list of ClosestEdgeQueryWithTargets, each containing an index, the
       * S2ClosestEdgeQuery for that index, and query targets, with geometry as specified by the
       * current parameters. If encodeIndex is true, the index and the shapes in it will encoded
       * before the query is constructed on the index.
       */
      furthestEdgeQueriesAndTargets = new ArrayList<>();
      for (int i = 0; i < numIndexSamples; i++) {
        furthestEdgeQueriesAndTargets.add(
            generateFurthestEdgeQueryWithTargets(numTargetSamples, encodeIndex, builder));
      }
    }

    /**
     * Generates a single FurthestEdgeQueryWithTargets containing an index, query, and targets,
     * using the given S2FurthestEdgeQuery.Builder to produce the query.
     */
    public FurthestEdgeQueryWithTargets generateFurthestEdgeQueryWithTargets(
        int numTargets, boolean encodeIndex, S2FurthestEdgeQuery.Builder builder)
        throws IOException {
      // The empty container to fill.
      FurthestEdgeQueryWithTargets output = new FurthestEdgeQueryWithTargets();

      // Generate and set output.index. If the index is encoded, also fills output.shapes.
      // The index is built at a random location which is returned in indexCap.
      S2Cap indexCap = fillOutputIndex(output, encodeIndex);

      // Build the output.query on the new index, using the provided Builder.
      output.query = builder.build(output.index);

      // Generate and set output.targets for the query.
      output.targets = new ArrayList<>();
      for (int i = 0; i < numTargets; ++i) {
        S2Cap targetCap = getTargetCap(indexCap);

        if (targetType == TargetType.POINT) {
          output.targets.add(
              new S2FurthestEdgeQuery.PointTarget<>(getTargetPoint(output, targetCap)));
        } else if (targetType == TargetType.EDGE) {
          MutableEdge edge = getTargetEdge(output, targetCap);
          output.targets.add(new S2FurthestEdgeQuery.EdgeTarget<>(edge.getStart(), edge.getEnd()));
        } else if (targetType == TargetType.CELL) {
          output.targets.add(
              new S2FurthestEdgeQuery.CellTarget<>(getTargetCell(output, targetCap)));
        } else {
          Preconditions.checkState(targetType == TargetType.INDEX);
          S2FurthestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
              S2FurthestEdgeQuery.createShapeIndexTarget(getTargetIndex(output, targetCap));
          target.setIncludeInteriors(includeInteriors);
          output.targets.add(target);
        }
      }
      return output;
    }

    /**
     * Run the findFurthestEdge() on the current query and target, and then go to the next ones. The
     * value for maxResults that was set on the query is ignored.
     */
    protected Optional<Result<S1ChordAngle>> findFurthestEdgeAndAdvance() {
      Optional<Result<S1ChordAngle>> result = currentQuery.findFurthestEdge(currentTarget);
      next();
      return result;
    }

    /**
     * Run the findFurthestEdges() on the current query and target, and then go to the next ones.
     */
    protected List<Result<S1ChordAngle>> findFurthestEdgesAndAdvance() {
      List<Result<S1ChordAngle>> result = currentQuery.findFurthestEdges(currentTarget);
      next();
      return result;
    }
  }

  /**
   * Benchmark for determining good values for maxBruteForceIndexSize() for the different targets.
   * For usage, see the comments on
   * {@link S2ClosestEdgesQueryBenchmark.DetermineGoodBruteForceIndexSizes}.
   *
   * <p>As of July 20, 2022, testing on a Cloudtop instance, the crossover point on index sizes was
   * for Point: between 64 and 96, for edge: between 40 and 42, for cell: Over 96, and for Index
   * about 32.
   */
  public static class DetermineGoodBruteForceIndexSizes extends FurthestEdgeBenchmarkBaseState {
    @Param({"POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"20", "24", "32", "40", "42", "48", "56", "64", "96"})
    int numIndexEdgesParam;

    @Param({"POINT", "EDGE", "CELL", "INDEX"})
    TargetType targetTypeParam;

    @Param({"true", "false"})
    boolean useBruteForceParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;
      useBruteForce = useBruteForceParam;

      targetType = targetTypeParam;
      includeInteriors = false;
      centerSeparationFraction = 1.0;

      // numTargetEdges only matters for the INDEX target.
      numTargetEdges = 16;
      chooseTargetFromIndex = false;
      targetRadiusFraction = 1.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to get the single furthest edge of various indexed geometries to a point
   * target in the general vicinity. Compare to BM_FindFurthest and BM_FindFurthestInterior in the
   * C++ tests.
   */
  public static class FindFurthestEdgeToPointTarget extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "48", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Param({"true", "false"})
    boolean includeInteriorsParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.POINT;
      includeInteriors = includeInteriorsParam;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Like the benchmark above, but benchmarks the time to get multiple furthest edges of various
   * indexed geometries to a point target in the general vicinity.
   */
  public static class FindManyFurthestEdgesToPointTarget extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"48", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Param({"1", "8", "16", "32"})
    int maxResultsParam;

    @Param({"true", "false"})
    boolean includeInteriorsParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;
      maxResults = maxResultsParam;

      targetType = TargetType.POINT;
      includeInteriors = includeInteriorsParam;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public List<Result<S1ChordAngle>> findFurthestEdges() {
      List<Result<S1ChordAngle>> furthestEdges = findFurthestEdgesAndAdvance();
      // In a unit test, numIndexEdges is forced set to 4, ignoring the @Param.
      Preconditions.checkState(isUnitTest() || furthestEdges.size() == maxResults,
          "Expected furthestEdges.size(%s) to equal maxResults (%s)",
          furthestEdges.size(),
          maxResults);
      return furthestEdges;
    }
  }

  /**
   * Benchmarks the time to get the single furthest edge of various indexed geometries to a point
   * target in the general vicinity allowing 1% error, which should make searches faster. Compare to
   * BM_FindFurthestMaxErrorPct.
   */
  public static class FindFurthestEdgeToPointTargetMaxErrorPct
      extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"49152"})
    int numIndexEdgesParam;

    @Param({"1", "10"})
    float maxErrorPercent;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.POINT;
      maxErrorFraction = 0.01 * maxErrorPercent;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to get furthest edges of various indexed geometries to a single edge or
   * cell target in the general vicinity.
   *
   * <p>Compare to BM_FindFurthestToEdge, BM_FindFurthestToCell, BM_FindFurthestToEdgeInterior,
   * and BM_FindFurthestToCellInterior.
   */
  public static class FindFurthestToEdgeOrCell extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "768", "49152"})
    int numIndexEdgesParam;

    @Param({"EDGE", "CELL"})
    TargetType targetTypeParam;

    @Param({"true", "false"})
    boolean includeInteriorsParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = targetTypeParam;
      includeInteriors = includeInteriorsParam;
      targetRadiusFraction = -1.0;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to measure the distance from an S2ShapeIndex with "numIndexEdges" edges to
   * a small (12 edge) S2ShapeIndex target such that their bounding S2Caps touch.
   */
  public static class FindFurthestToSmallAbuttingIndex extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"48", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = 12;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Like the benchmark above, except that the two indexes are reversed.
   *
   * <p>TODO(torrey): In C++ it is apparently somewhat faster to use the larger of the two
   * S2ShapeIndexes as the target. Compare.
   */
  public static class FindFurthestFromSmallAbuttingIndex extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"768", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = 12;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /** Like the benchmark above, except that the two indexes have the same size. */
  public static class FindFurthestFromSameSizeAbuttingIndex extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Measure the distance between two indexes with the same number of edges and whose bounding
   * S2Caps are distant relative to their radius.
   */
  public static class FindFurthestFromSameSizeDistantIndex extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 10.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Like the test above, but instead just compare the maximum distance to a threshold that is
   * slightly greater than the maximum possible value. (The centers of the two geometries are 10
   * radii apart, so the distance is at most 12 radii. We set the threshold slightly higher than
   * this.)
   */
  public static class IsDistanceGreaterSameSizeDistantIndexFalse
      extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      distanceFraction = 12.4;
      maxErrorFraction = 12.4;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 10.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }

  /**
   * Like the test above, but instead just compare the maximum distance to a threshold that is
   * usually slightly less than the true value. (The centers of the two geometries are 10 radii
   * apart, and by experimentation the maximum distance between the two fractals is at most about 12
   * radii. We set the threshold slightly lower than this.)
   */
  public static class IsDistanceGreaterSameSizeDistantIndexTrue
      extends FurthestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      distanceFraction = 11.6;
      maxErrorFraction = 11.6;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 10.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findFurthestEdge() {
      return findFurthestEdgeAndAdvance();
    }
  }
}
