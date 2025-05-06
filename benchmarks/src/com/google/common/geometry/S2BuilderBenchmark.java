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

import static com.google.common.geometry.TestDataGenerator.metersToAngle;
import static java.util.concurrent.TimeUnit.MILLISECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.base.Preconditions;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2Builder;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2EdgeUtil;
import com.google.common.geometry.S2Error;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polyline;
import com.google.common.geometry.S2PolylineLayer;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.ArrayList;
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

/** Benchmarks for {@link S2Builder}. */
@CheckReturnValue
public class S2BuilderBenchmark {
  private S2BuilderBenchmark() {}

  /**
   * Benchmark state for S2Builder. For each trial, constructs a polyline with 'numEdges' edges that
   * goes back and forth repeatedly between two disc-shaped areas spaced 100 meters apart and 100
   * microns in radius. The vertices are randomly located within the two discs which creates a large
   * number of self-intersections (approximately 0.25 * p0 ** 2). This benchmark uses
   * splitCrossingEdges() which adds a new vertex at each intersection, and since the lines are very
   * close together many lines snap to nearby intersection points in addition to their own.
   * S2Builder also adds extra vertices to ensure that edges are well-separated from non-incident
   * vertices (see {@link S2Builder.SnapFunction#minEdgeVertexSeparation()}).
   *
   * <p>With 1000 input edges there are ~250,000 intersection points. Using a 1 micron snap radius
   * each edge snaps to ~4,500 of these intersection points and so the output has 4.5 million edges.
   * Using a 10 micron snap radius each edge snaps to ~45,000 intersection points and the output has
   * about 45 million edges.
   */
  @State(Scope.Thread)
  public static class NearlyParallelCrossingEdgesState extends S2BenchmarkBaseState {
    @Param({"10", "25", "50", "100", "200"})
    int numEdgesParam;

    @Param({"0", "1"})
    int snapRadiusMicronsParam;

    protected S2Polyline polyline;
    protected S2Builder.Builder options;
    protected S2Error error = new S2Error();

    @Setup(Level.Iteration)
    @Override
    public void setup() {
      super.setup();
      options = new S2Builder.Builder(
          new IdentitySnapFunction(metersToAngle(snapRadiusMicronsParam / 1e6)));
      options.setSplitCrossingEdges(true);

      // Construct the input polyline.
      S1Angle edgeLength = metersToAngle(100);
      S1Angle perturbRadius = metersToAngle(0.0001); // 100 microns

      S2Point a = data.getRandomPoint();
      S2Point b = S2EdgeUtil.getPointOnLine(a, data.getRandomPoint(), edgeLength);
      S2Cap capA = S2Cap.fromAxisAngle(a, perturbRadius);
      S2Cap capB = S2Cap.fromAxisAngle(b, perturbRadius);

      ArrayList<S2Point> vertices = new ArrayList<>();
      for (int i = 0; i < numEdgesParam; i++) {
        vertices.add(data.samplePoint(capA));
        vertices.add(data.samplePoint(capB));
      }
      polyline = new S2Polyline(vertices);
    }

    /**
     * Measures the time to create an S2Builder, create and add a layer, and run S2Builder.build()
     * on the input polyline constructed above.
     */
    @Benchmark
    @BenchmarkMode(Mode.AverageTime)
    @OutputTimeUnit(MILLISECONDS)
    @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
    @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
    public S2Polyline bmBuilderConstructAndBuild() {
      S2Builder builder = options.build();
      S2PolylineLayer layer = new S2PolylineLayer();
      builder.startLayer(layer);
      error.clear();
      builder.addPolyline(polyline);
      boolean ok = builder.build(error);
      Preconditions.checkState(ok, error);
      return layer.getPolyline();
    }
  }
}
