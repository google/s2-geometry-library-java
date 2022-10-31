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
import com.google.common.geometry.S2BestEdgesQueryBase;
import com.google.common.geometry.S2BestEdgesQueryBase.Result;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2ClosestEdgeQuery;
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

/** Benchmarks for {@link S2ClosestEdgeQuery}. */
public class S2ClosestEdgeQueryBenchmark {
  private S2ClosestEdgeQueryBenchmark() {}

  /**
   * A container for an index, a closest edge query on that index, and a set of targets to test
   * with.
   */
  public static class ClosestEdgeQueryWithTargets extends BestEdgeQueryWithTargets {
    public S2ClosestEdgeQuery<S1ChordAngle> query;
    public ArrayList<S2ClosestEdgeQuery.Target<S1ChordAngle>> targets;

    @Override
    public S2BestEdgesQueryBase<S1ChordAngle> query() {
      return query;
    }

    @Override
    public ArrayList<S2ClosestEdgeQuery.Target<S1ChordAngle>> targets() {
      return targets;
    }
  }

  public static class ClosestEdgeBenchmarkBaseState extends BestEdgesBenchmarkBaseState {
    protected ArrayList<ClosestEdgeQueryWithTargets> closestEdgeQueriesAndTargets;
    protected S2ClosestEdgeQuery<S1ChordAngle> currentQuery;
    protected S2ClosestEdgeQuery.Target<S1ChordAngle> currentTarget;

    @Override
    void setCurrentQuery() {
      currentQuery = closestEdgeQueriesAndTargets.get(queryIndex).query;
    }

    @Override
    void setCurrentTarget() {
      currentTarget = closestEdgeQueriesAndTargets.get(queryIndex).targets.get(targetIndex);
    }

    @Override
    protected void createQueriesAndTargets() throws IOException {
      // Build the closest edge query options.
      S2ClosestEdgeQuery.Builder builder =
          S2ClosestEdgeQuery.builder()
              .setUseBruteForce(useBruteForce)
              .setIncludeInteriors(includeInteriors)
              .setMaxResults(maxResults);

      if (distanceFraction > 0) {
        builder.setMaxDistance(indexRadius.mul(distanceFraction));
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
      closestEdgeQueriesAndTargets = new ArrayList<>();
      for (int i = 0; i < numIndexSamples; i++) {
        closestEdgeQueriesAndTargets.add(
            generateClosestEdgeQueryWithTargets(numTargetSamples, encodeIndex, builder));
      }
    }

    /**
     * Generates a single ClosestEdgeQueryWithTargets containing an index, query, and targets, using
     * the given S2ClosestEdgeQueryBuilder to produce the query.
     */
    public ClosestEdgeQueryWithTargets generateClosestEdgeQueryWithTargets(
        int numTargets, boolean encodeIndex, S2ClosestEdgeQuery.Builder builder)
        throws IOException {
      // The empty container to fill.
      ClosestEdgeQueryWithTargets output = new ClosestEdgeQueryWithTargets();

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
              new S2ClosestEdgeQuery.PointTarget<>(getTargetPoint(output, targetCap)));
        } else if (targetType == TargetType.EDGE) {
          MutableEdge edge = getTargetEdge(output, targetCap);
          output.targets.add(new S2ClosestEdgeQuery.EdgeTarget<>(edge.getStart(), edge.getEnd()));
        } else if (targetType == TargetType.CELL) {
          output.targets.add(new S2ClosestEdgeQuery.CellTarget<>(getTargetCell(output, targetCap)));
        } else {
          Preconditions.checkState(targetType == TargetType.INDEX);
          S2ClosestEdgeQuery.ShapeIndexTarget<S1ChordAngle> target =
              S2ClosestEdgeQuery.createShapeIndexTarget(getTargetIndex(output, targetCap));
          target.setIncludeInteriors(includeInteriors);
          output.targets.add(target);
        }
      }
      return output;
    }

    /**
     * Run findClosestEdge() on the current query and target, and then go to the next ones. The
     * value for maxResults that was set on the query is ignored.
     */
    protected Optional<Result<S1ChordAngle>> findClosestEdgeAndAdvance() {
      Optional<Result<S1ChordAngle>> result = currentQuery.findClosestEdge(currentTarget);
      next();
      return result;
    }

    /** Run findClosestEdges() on the current query and target, and then go to the next ones. */
    protected List<Result<S1ChordAngle>> findClosestEdgesAndAdvance() {
      List<Result<S1ChordAngle>> result = currentQuery.findClosestEdges(currentTarget);
      next();
      return result;
    }
  }

  /**
   * Benchmark for determining good values for maxBruteForceIndexSize() for the different targets.
   * These are the *index* sizes at which point a the ClosestEdgeQuery should switch to the
   * optimized algorithm. See {@link S2BestEdgesQueryBase#findBestEdgesInternal()}.
   *
   * <p>To use this benchmark:
   * <ol>
   * <li> Set all the maxBruteForceIndexSize methods to return 1 in the S2ClosestEdgeQuery.Target
   *      implementations, so it uses the optimized algorithm unless forced not to.
   * <li> Run this benchmark as follows:
   * <li> export B=javatests/com/google/common/geometry/benchmarks
   * <li> blaze run -c opt $B:S2ClosestEdgeQueryBenchmark_jmh DetermineGoodBruteForceIndexSizes
   * <li> In the results, find the index size crossover points where the brute force is no longer
   *      faster than the optimized approach.  Adjust the @Param for numIndexEdges as needed to
   *      narrow the range, but it is not important to get an exact number as the best value likely
   *      depends on many additional factors not considered here.
   * <li> Update the code comments below about where the crossover occurs and the benchmark
   *      environment used.
   * <li> Update the maxBruteForceIndexSize implementations.
   * </ol>
   *
   * <p>As of July 20, 2022, testing on a Cloudtop instance: For POINT, EDGE, and INDEX (with 16
   * edges) the crossover point is between 40 and 42. For CELL the crossover point is between 16
   * and 20.
   */
  public static class DetermineGoodBruteForceIndexSizes extends ClosestEdgeBenchmarkBaseState {
    @Param({"POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"16", "20", "32", "40", "42", "48", "56"})
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
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to get the single closest edge of various indexed geometries to a point
   * target in the general vicinity, or (when includeInteriors is true) an edge or interior of those
   * geometries.
   */
  public static class FindClosestEdgeToPointTarget extends ClosestEdgeBenchmarkBaseState {
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
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Like the benchmark above, but benchmarks the time to get multiple closest edges of various
   * indexed geometries to a point target in the general vicinity, or (when includeInteriors is
   * true) an edge or interior of those geometries.
   */
  public static class FindManyClosestEdgesToPointTarget extends ClosestEdgeBenchmarkBaseState {
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
    public List<Result<S1ChordAngle>> findClosestEdges() {
      List<Result<S1ChordAngle>> closestEdges = findClosestEdgesAndAdvance();
      // In a unit test, numIndexEdges is forced set to 4, ignoring the @Param.
      Preconditions.checkState(isUnitTest() || closestEdges.size() == maxResults,
          "Expected closestEdges.size(%s) to equal maxResults (%s)",
          closestEdges.size(),
          maxResults);
      return closestEdges;
    }
  }

  /**
   * Benchmarks the time to get closest edges of various indexed geometries to a point target in the
   * general vicinity, using the brute force implementation.
   */
  public static class FindClosestBruteForce extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "48", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.POINT;
      centerSeparationFraction = -2.0;
      useBruteForce = true;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to get closest edges of various geometries to a point target in the general
   * vicinity, allowing error tolerance.
   *
   * <p>Compare to C++ benchmark BM_FindClosestMaxErrorPct.
   */
  public static class FindClosestMaxErrorPct extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"49152"})
    int numIndexEdgesParam;

    @Param({"1", "10"})
    int maxErrorPercent;

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
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks searching with a small distance limit. This case is important for applications where
   * the distance limit is small enough that there may not be any edges in the result.
   *
   * <p>Compare to C++ benchmark BM_FindClosestMaxDistPow10
   */
  public static class FindClosestMaxDistPow10 extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"49152"})
    int numIndexEdgesParam;

    // For indexRadiusRadians = 1000 km (the default), this benchmarks with maxDistance = 10km (-2)
    // and one meter (-6).
    @Param({"-2", "-6"})
    int maxDistanceFractionPow10;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.POINT;
      distanceFraction = Math.pow(10, maxDistanceFractionPow10);
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmark searching where the query point is a vertex in the index. This measures performance
   * in the case where an edge is guaranteed to be nearby.
   *
   * <p>Compare to C++ benchmark BM_FindClosestNearVertex.
   */
  public static class FindClosestNearVertex extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"768", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.POINT;
      chooseTargetFromIndex = true;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmark searching where the query point is a vertex and the distance limit is small. This
   * case is important for clustering nearby points.
   *
   * <p>Compare to C++ benchmark BM_FindClosestNearVertexMaxDistPow10.
   */
  public static class FindClosestNearVertexMaxDistPow10 extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"49152"})
    int numIndexEdgesParam;

    // For indexRadiusRadians = 1000 km (the default), this sets maxDistance to one meter.
    @Param({"-6"})
    int maxDistanceFractionPow10;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.POINT;
      distanceFraction = Math.pow(10, maxDistanceFractionPow10);
      chooseTargetFromIndex = true;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to get closest edges of various indexed geometries to a single edge or cell
   * target in the general vicinity.
   *
   * <p>Compare to C++ benchmarks BM_FindClosestToEdge, BM_FindClosestToEdgeInterior, and
   * BM_FindClosestToCell.
   */
  public static class FindClosestToEdgeOrCell extends ClosestEdgeBenchmarkBaseState {
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
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to get the closest edge of a fractal loop to either a single edge target
   * chosen from one of the index edges, or a cell target chosen from one of the index cells.
   *
   * <p>Compare to C++ benchmarks BM_FindClosestToEdgeNearEdge and BM_FindClosestToCellInIndex.
   */
  public static class FindClosestToEdgeOrCellInIndex extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"768", "49152"})
    int numIndexEdgesParam;

    @Param({"EDGE", "CELL"})
    TargetType targetTypeParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = targetTypeParam;
      chooseTargetFromIndex = true;
      targetRadiusFraction = -1.0;
      centerSeparationFraction = -2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time to measure the distance from an S2ShapeIndex with "numIndexEdges" edges to
   * a small (12 edge) S2ShapeIndex target such that their bounding S2Caps touch.
   *
   * <p>Compare to C++ benchmark BM_FindClosestToSmallAbuttingIndex.
   */
  public static class FindClosestToSmallAbuttingIndex extends ClosestEdgeBenchmarkBaseState {
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
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Like the benchmark above, except that the two indexes are reversed.
   *
   * <p>Compare to C++ benchmark BM_FindClosestFromSmallAbuttingIndex.
   */
  public static class FindClosestFromSmallAbuttingIndex extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"48", "768", "12288", "49152"})
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
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time required to measure the distance between two indexes with the same number
   * of edges and whose bounding S2Caps touch each other.
   */
  public static class FindClosestToSameSizeAbuttingIndex extends ClosestEdgeBenchmarkBaseState {
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
      chooseTargetFromIndex = true;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 2.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time required to measure the distance from an S2ShapeIndex to a contained
   * S2ShapeIndex target with the same number of edges, *including* the interiors of both indexes.
   * This operation is very fast for polygonal geometry because the index interior contains the
   * target.
   */
  public static class FindClosestToSameSizeContainedIndex extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"48", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      includeInteriors = true;
      chooseTargetFromIndex = true;
      targetRadiusFraction = 0.25;
      centerSeparationFraction = 0.25;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time required to measure the distance from an S2ShapeIndex to a containing
   * S2ShapeIndex target with the same number of edges, *including* the interiors of both indexes.
   * This is slower for polygonal geometry than the benchmark above because the index interior does
   * not contain the target.
   */
  public static class FindClosestToSameSizeContainingIndex extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"48", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      chooseTargetFromIndex = true;
      targetRadiusFraction = 4.0;
      centerSeparationFraction = 1.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Benchmarks the time required to measure the distance between two indexes with the same number
   * of edges and whose bounding S2Caps are distant relative to their radius.
   */
  public static class FindClosestToSameSizeDistantIndex extends ClosestEdgeBenchmarkBaseState {
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
      chooseTargetFromIndex = true;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 10.0;

      super.initialize(params);
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }
  }

  /**
   * Like the test above, but instead just benchmarks the time required to compare the minimum
   * distance to thresholds that are slightly less or more than the true value. The centers of the
   * two geometries are 10 radii apart, so the distance is at least 8 radii. For the false case, we
   * set the threshold slightly lower than this, at 7.6 radii, and for the true case, to slightly
   * higher, at 8.4 radii.
   */
  public static class IsDistanceLessSameSizeDistantIndex extends ClosestEdgeBenchmarkBaseState {
    @Param({"FRACTAL_LOOP", "REGULAR_LOOP", "POINT_CLOUD"})
    FactoryType factoryTypeParam;

    @Param({"12", "768", "12288", "49152"})
    int numIndexEdgesParam;

    @Param({"true", "false"})
    boolean isLess;

    private S1ChordAngle threshold;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      factoryType = factoryTypeParam;
      numIndexEdges = numIndexEdgesParam;

      targetType = TargetType.INDEX;
      numTargetEdges = numIndexEdgesParam;
      distanceFraction = 7.6;
      maxErrorFraction = 7.6;
      chooseTargetFromIndex = true;
      targetRadiusFraction = 1.0;
      centerSeparationFraction = 10.0;

      super.initialize(params);

      double radii = isLess ? 7.6 : 8.4;
      threshold = S1ChordAngle.fromS1Angle(fractionToRadiusOrRange(radii, indexRadius));
    }

    @Benchmark
    public Optional<Result<S1ChordAngle>> findClosestEdge() {
      return findClosestEdgeAndAdvance();
    }

    @Benchmark
    public boolean isDistanceLess() {
      boolean result = currentQuery.isDistanceLess(currentTarget, threshold);
      next();
      return result;
    }
  }
}
