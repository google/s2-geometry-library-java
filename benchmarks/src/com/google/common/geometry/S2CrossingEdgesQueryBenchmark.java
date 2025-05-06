/*
 * Copyright 2023 Google Inc.
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

import static java.lang.Math.PI;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.geometry.Matrix;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2CrossingEdgesQuery;
import com.google.common.geometry.S2CrossingEdgesQuery.CrossingType;
import com.google.common.geometry.S2CrossingEdgesQuery.EdgePairVisitor;
import com.google.common.geometry.S2LatLng;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polyline;
import com.google.common.geometry.S2ShapeIndex;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.IOException;
import java.util.ArrayList;
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
import org.openjdk.jmh.infra.BenchmarkParams;

/** Benchmarks for S2CrossingEdgesQuery. */
public class S2CrossingEdgesQueryBenchmark {

  private S2CrossingEdgesQueryBenchmark() {}

  private static class IndexPair {
    public final S2ShapeIndex first;
    public final S2ShapeIndex second;

    public IndexPair() {
      first = new S2ShapeIndex();
      second = new S2ShapeIndex();
    }
  }

  /**
   * Benchmark state that benchmarks finding crossing edges between on two indexes. Subclasses must
   * set up indexPair and opType and should call init().
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 2, time = 3, timeUnit = SECONDS)
  @Measurement(iterations = 3, time = 3, timeUnit = SECONDS)
  public static class TwoIndexState extends S2BenchmarkBaseState {
    protected CrossingType crossingType;
    protected IndexPair indexPair;

    protected void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
    }

    /**
     * Creates an S2CrossingEdgesQuery with the current crossingType, and calls
     * visitCrossingEdgePairs() for the two indexes in 'indexPair' .
     */
    public int visitCrossingEdgePairs() {
      S2CrossingEdgesQuery query = new S2CrossingEdgesQuery(crossingType);
      int[] numCrossings = {0};
      query.visitCrossingEdgePairs(
          indexPair.first,
          indexPair.second,
          new EdgePairVisitor() {
            @Override
            public boolean visit(
                int shapeIdA,
                int edgeIdA,
                S2Point edgeASrc,
                S2Point edgeADst,
                int shapeIdB,
                int edgeIdB,
                S2Point edgeBSrc,
                S2Point edgeBDst,
                boolean isInterior) {
              numCrossings[0]++;
              return true; // Continue visiting.
            }
          });
      return numCrossings[0];
    }

    /** Applies updates to the given IndexPair and returns it. */
    @CanIgnoreReturnValue
    public static IndexPair buildIndexPair(IndexPair indexes) {
      indexes.first.applyUpdates();
      indexes.second.applyUpdates();
      return indexes;
    }
  }

  /** Benchmarks many intersections. */
  public static class PerpendicularZigZagsState extends TwoIndexState {
    @Param({"INTERIOR", "ALL"})
    CrossingType crossingTypeParam;
    @Param({"25", "50", "100"})
    int numIndexEdgesParam;

    @Override
    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup(params);
      crossingType = crossingTypeParam;
      indexPair = makePerpendicularZigZags(numIndexEdgesParam);
    }

    @Benchmark
    public int run() {
      return super.visitCrossingEdgePairs();
    }

    /**
     * Creates two zig-zag polylines that intersect in O(n**2) points. One polyline with O(n)
     * vertices is added to each index. The indexes are returned in an IndexPair.
     */
    public IndexPair makePerpendicularZigZags(int n) {
      ArrayList<S2LatLng> aVertices = new ArrayList<>();
      ArrayList<S2LatLng> bVertices = new ArrayList<>();
      double kSizeDegrees = 50;
      for (int i = 0; i < n; ++i) {
        double x = i * kSizeDegrees / (n - 1);
        double y = (i & 1) * kSizeDegrees;
        aVertices.add(S2LatLng.fromDegrees(x, y));
        bVertices.add(S2LatLng.fromDegrees(y, x));
      }
      return makePolylineIndexPair(aVertices, bVertices);
    }

    /**
     * Converts the two given lists of LatLngs into two polylines, returning them in an IndexPair.
     */
    private static IndexPair makePolylineIndexPair(
        List<S2LatLng> aVertices, List<S2LatLng> bVertices) {
      IndexPair indexes = new IndexPair();
      List<S2Point> aPoints = Lists.transform(aVertices, S2LatLng::toPoint);
      List<S2Point> bPoints = Lists.transform(bVertices, S2LatLng::toPoint);
      indexes.first.add(new S2Polyline(aPoints));
      indexes.second.add(new S2Polyline(bPoints));
      return buildIndexPair(indexes);
    }
  }

  /** Tests many edges with intersections. */
  public static class OffsetLoopsState extends TwoIndexState {
    @Param({"INTERIOR", "ALL"})
    CrossingType crossingTypeParam;
    @Param({"25", "50", "100"})
    int numIndexEdgesParam;

    @Override
    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup(params);
      crossingType = crossingTypeParam;
      indexPair = makeOffsetLoops(numIndexEdgesParam);
    }

    @Benchmark
    public int run() {
      return super.visitCrossingEdgePairs();
    }

    /**
     * Creates two loops such that every edge of one loop intersects exactly one edge of the other
     * loop. This results in many edges that have an intersection (as opposed to many
     * intersections). The two loops are returned in an IndexPair.
     */
    public IndexPair makeOffsetLoops(int n) {
      S1Angle kRadius = S1Angle.degrees(50);
      IndexPair indexes = new IndexPair();

      // Construct two regular loops around the same center point C where one loop is rotated by
      // half the angle subtended by each edge.
      S2Point a = new S2Point(0, 1, 0);
      S2Point b = new S2Point(0, 0, 1);
      S2Point c = new S2Point(1, 0, 0);
      indexes.first.add(new S2Loop(S2Loop.makeRegularLoop(Matrix.fromCols(a, b, c), kRadius, n)));
      a = a.rotate(c, PI / n);
      b = b.rotate(c, PI / n);
      indexes.second.add(new S2Loop(S2Loop.makeRegularLoop(Matrix.fromCols(a, b, c), kRadius, n)));
      return buildIndexPair(indexes);
    }
  }

  /** Tests finding crossing edges between two point clouds. There will be none. */
  public static class PointCloudState extends TwoIndexState {
    @Param({"INTERIOR", "ALL"})
    CrossingType crossingTypeParam;
    @Param({"500", "1000", "200"})
    int numIndexEdgesParam;

    @Override
    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup(params);
      crossingType = crossingTypeParam;
      indexPair = makeDisjointPointClouds(numIndexEdgesParam);
    }

    @Benchmark
    public void run() {
      Preconditions.checkState(super.visitCrossingEdgePairs() == 0);
    }

    /** Creates two point clouds with 'n' points each, returning them in an IndexPair. */
    public IndexPair makeDisjointPointClouds(int n) {
      IndexPair indexes = new IndexPair();
      indexes.first.add(makePointCloud(n));
      indexes.second.add(makePointCloud(n));
      return buildIndexPair(indexes);
    }

    /** Returns a cloud of 'n' random points. */
    public S2Point.Shape makePointCloud(int n) {
      ArrayList<S2Point> vertices = new ArrayList<>();
      for (int i = 0; i < n; ++i) {
        vertices.add(data.getRandomPoint());
      }
      return S2Point.Shape.fromList(vertices);
    }
  }
}
