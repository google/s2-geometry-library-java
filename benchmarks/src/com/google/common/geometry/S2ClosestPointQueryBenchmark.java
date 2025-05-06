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

import static com.google.common.geometry.TestDataGenerator.addPoints;
import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2ClosestPointQuery;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2PointIndex;
import com.google.common.geometry.TestDataGenerator;
import com.google.common.geometry.TestDataGenerator.PointFactory;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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

/** Benchmarks for {@link S2ClosestPointQuery}. */
@CheckReturnValue
public class S2ClosestPointQueryBenchmark {
  private S2ClosestPointQueryBenchmark() {}

  private static final S1Angle QUERY_RADIUS = kmToAngle(10.0);

  /**
   * Runs either findClosestPoints or findClosestPointsToEdge for one of the preconfigured queries
   * in the given {@link BMData} and its corresponding preset sample points or edges.
   */
  public enum TargetType {
    POINT_TARGET {
      @Override
      List<S2ClosestPointQuery.Result<Integer>> measure(BMData benchmarkData, int rep) {
        int size = benchmarkData.queries.size();
        int cap = rep % size;
        S2Point[] samples = benchmarkData.samples[cap];
        S2Point p = samples[rep % samples.length];
        return benchmarkData.queries.get(cap).findClosestPoints(p);
      }
    },
    EDGE_TARGET {
      @Override
      List<S2ClosestPointQuery.Result<Integer>> measure(BMData benchmarkData, int rep) {
        int size = benchmarkData.queries.size();
        int cap = rep % size;
        S2ClosestPointQuery<Integer> query = benchmarkData.queries.get(cap);
        S2Point[] samples = benchmarkData.samples[cap];
        S2Point a = samples[rep % samples.length];
        S2Point b = samples[(rep + 1) % samples.length];
        return query.findClosestPointsToEdge(a, b);
      }
    };

    abstract List<S2ClosestPointQuery.Result<Integer>> measure(BMData benchmarkData, int reps);
  }

  /**
   * Provides BENCHMARK_CAP_COUNT preconfigured S2ClosestPointQueries, each constructed over a
   * different S2PointIndex containing a set of indexed points.
   *
   * <p>In addition, provides multiple arrays of samples, one for each preconfigured query, which
   * are used by the {@link TargetType} above to repeatedly run the findClosestPoints or
   * findClosestPointsToEdge queries.
   *
   * <p>The queries build up internal state over multiple operations, so they are periodically
   * reset by calling maybeResetQuery(), which resets one query every RESET_FREQUENCY calls, so each
   * query is used around 50 times before being reset.
   */
  private static class BMData {
    /** Number of preconfigured queries and their corresponding caps to create. */
    private static final int BENCHMARK_CAP_COUNT = 8;
    /** Number of query points to sample in each cap. */
    private static final int BENCHMARK_SAMPLE_COUNT = 512;
    /** The preconfigured queries. They should be reset periodically. */
    private final List<S2ClosestPointQuery<Integer>> queries = new ArrayList<>();
    /** Samples for all iterations; the inner arrays are parallel to {@link #queries}. */
    private final S2Point[][] samples = new S2Point[BENCHMARK_CAP_COUNT][BENCHMARK_SAMPLE_COUNT];

    /** How often to reset one of the preconfigured queries. */
    private static final int RESET_FREQUENCY = 50 * BMData.BENCHMARK_CAP_COUNT;
    /** Which preconfigured query to reset next. */
    private int resetQuery;

    /**
     * Constructs the preconfigured queries, their query caps, and sample points for each query to
     * index, using the given PointFactory to sample points for indexing. The queries are configured
     * using the provided maxDistanceFraction to set the query maxDistance, with maxPoints to return
     * always equal to 1. If the provided useBruteForce is true, the queries are set to
     * useBruteForce, but otherwise the queries will determine internally to use brute force or not.
     * If chooseTargetFromIndex is true, the sample points are selected from the indexed points.
     * Otherwise, the sample points are selected at random as appropriate for the given targetType
     * of points or edges.
     */
    BMData(
        TestDataGenerator data,
        PointFactory factory,
        int numPoints,
        double maxDistanceFraction,
        TargetType targetType,
        boolean chooseTargetFromIndex,
        boolean useBruteForce) {
      // When creating edge endpoints, from what size cap around the first endpoint should the
      // second be sampled?
      S1Angle edgeLengthMax = S1Angle.radians(QUERY_RADIUS.radians() / 2.0);

      // Create BENCHMARK_CAP_COUNT queries over samples of points from various random caps.
      // Performance is affected by how the points are positioned relative to the S2Cell hierarchy.
      for (int i = 0; i < BENCHMARK_CAP_COUNT; i++) {
        S2PointIndex<Integer> index = new S2PointIndex<>();
        S2Cap queryCap = S2Cap.fromAxisAngle(data.getRandomPoint(), QUERY_RADIUS);
        List<S2Point> indexPoints = factory.createPoints(data, queryCap, numPoints);
        addPoints(index, indexPoints);
        S2ClosestPointQuery<Integer> query = new S2ClosestPointQuery<>(index);
        query.setMaxPoints(1);
        if (maxDistanceFraction > 0) {
          query.setMaxDistance(S1Angle.radians(maxDistanceFraction * QUERY_RADIUS.radians()));
        }
        // Don't call useBruteForce(false) since this would disable the automatic use of brute force
        // for small numbers of edges.
        if (useBruteForce) {
          query.useBruteForce(true);
        }
        queries.add(query);
        if (chooseTargetFromIndex) {
          for (int j = 0; j < BENCHMARK_SAMPLE_COUNT; j++) {
            samples[i][j] = indexPoints.get(j % indexPoints.size());
          }
          indexPoints.toArray(samples[i]);
          Collections.shuffle(Arrays.asList(samples[i]));
        } else if (targetType == TargetType.POINT_TARGET) {
          for (int j = 0; j < BENCHMARK_SAMPLE_COUNT; j++) {
            samples[i][j] = data.samplePoint(queryCap);
          }
        } else if (targetType == TargetType.EDGE_TARGET) {
          for (int j = 0; j < BENCHMARK_SAMPLE_COUNT; j += 2) {
            // Create random line endpoints by sampling one end from anywhere in the query cap, and
            // the other end within a cap around the first endpoint. That end might be outside the
            // query cap.
            S2Point p = data.samplePoint(queryCap);
            samples[i][j] = p;
            samples[i][j + 1] = data.samplePoint(S2Cap.fromAxisAngle(p, edgeLengthMax));
          }
        } else {
          throw new RuntimeException("Unimplemented target type");
        }
      }
    }

    // Instead of periodically resetting all queries at once, which would cause a large difference
    // in performance for the benchmark rep where that happens, periodically reset just one query.
    public void maybeResetQuery(int rep) {
      if (rep % RESET_FREQUENCY == 0) {
        queries.get(resetQuery).reset();
        resetQuery = (resetQuery + 1) % queries.size();
      }
    }
  }

  /**
   * Benchmark state for finding the closest point or edge within the general vicinity of the index
   * points, with no limit on maximum search distance. Targets are sampled from the query cap
   * independently of the index points.
   */
  @State(Scope.Thread)
  public static class FindClosestState extends S2BenchmarkBaseState {
    @Param // PointFactories: CIRCLE, FRACTAL, and GRID.
    PointFactory factory;

    @Param // POINT_TARGET and EDGE_TARGET.
    TargetType target;

    @Param({"12", "48", "768", "12288", "196608"})
    int numPoints;

    protected BMData benchmarkData;
    protected int rep;

    // Reset all the data for each iteration, as the queries build up internal state.
    @Setup(Level.Iteration)
    @Override
    public void setup() {
      super.setup();
      benchmarkData = new BMData(data, factory, numPoints, -1.0, target, false, false);
      rep = 0;
    }

    /**
     * Measures the time to run findClosestPoints() or findClosestPointsToEdge() with the
     * parameterized benchmark state configurations defined above.
     */
    @Benchmark
    @BenchmarkMode(Mode.AverageTime)
    @OutputTimeUnit(MICROSECONDS)
    @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
    @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
    public List<S2ClosestPointQuery.Result<Integer>> bmFindClosest() {
      List<S2ClosestPointQuery.Result<Integer>> results = target.measure(benchmarkData, rep);
      benchmarkData.maybeResetQuery(rep);
      rep++;
      return results;
    }
  }

  /**
   * Benchmark state for finding the closest point with a distance limit. This case is important for
   * applications where the distance limit is small enough that there may not be any points in the
   * result. We don't benchmark searching with an edge target for small distance limits because it
   * is only beneficial when the edge is small compared to the spacing of the index points, which is
   * not the case here. This benchmark state also compares brute force vs. using the S2PointIndex.
   */
  @State(Scope.Thread)
  public static class FindClosestMaxDistPow10State extends S2BenchmarkBaseState {
    @Param PointFactory factory;

    @Param({"-2", "-6"})
    int maxDistLog10;

    protected final TargetType target = TargetType.POINT_TARGET;

    @Param({"true", "false"})
    boolean useBruteForce;

    protected BMData benchmarkData;
    protected int rep;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      benchmarkData =
          new BMData(
              data, factory, 196608, Math.pow(10, maxDistLog10), target, false, useBruteForce);
      rep = 0;
    }

    /**
     * Measures the time to run findClosestPoints() or findClosestPointsToEdge() with the
     * parameterized benchmark state configurations defined above.
     */
    @Benchmark
    @BenchmarkMode(Mode.AverageTime)
    @OutputTimeUnit(NANOSECONDS)
    @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
    @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
    public List<S2ClosestPointQuery.Result<Integer>>  bmFindClosest() {
      List<S2ClosestPointQuery.Result<Integer>> results = target.measure(benchmarkData, rep);
      benchmarkData.maybeResetQuery(rep);
      rep++;
      return results;
    }
  }

  /**
   * Benchmark state for finding the closest point or edge where the query points (or edge
   * endpoints) are selected from the indexed points ("sites").
   */
  @State(Scope.Thread)
  public static class FindClosestNearSiteState extends S2BenchmarkBaseState {
    @Param PointFactory factory;

    @Param TargetType target;

    @Param({"768", "196608"})
    int numPoints;

    protected BMData benchmarkData;
    protected int rep;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      benchmarkData = new BMData(data, factory, numPoints, -1, target, true, false);
    }

    /**
     * Measures the time to run findClosestPoints() or findClosestPointsToEdge() with the
     * parameterized benchmark state configurations defined above.
     */
    @Benchmark
    @BenchmarkMode(Mode.AverageTime)
    @OutputTimeUnit(NANOSECONDS)
    @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
    @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
    public List<S2ClosestPointQuery.Result<Integer>>  bmFindClosest() {
      List<S2ClosestPointQuery.Result<Integer>> results = target.measure(benchmarkData, rep);
      benchmarkData.maybeResetQuery(rep);
      rep++;
      return results;
    }
  }
}
