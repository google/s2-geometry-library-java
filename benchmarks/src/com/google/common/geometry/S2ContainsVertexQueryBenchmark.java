/*
 * Copyright 2025 Google Inc.
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

import static com.google.common.geometry.TestDataGenerator.DEFAULT_LOOP_RADIUS;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S2ContainsVertexQuery;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polygon;
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

/**
 * A simple benchmark for S2ContainsVertexQuery.
 */
public class S2ContainsVertexQueryBenchmark {

  private S2ContainsVertexQueryBenchmark() {}

  /**
   * Benchmark state consisting of a polygon of many triangles, all around a common center point,
   * forming a fan or windmill shape.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class PolygonState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768"})
    int numTriangles;

    private S2Polygon polygon;
    private S2Point center;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      center = data.getRandomPoint();
      List<S2Point> regularVertices = S2Loop.makeRegularVertices(
          center, DEFAULT_LOOP_RADIUS, numTriangles * 2);

      List<S2Loop> loops = new ArrayList<>();
      for (int tri = 0; tri < numTriangles; tri++) {
        List<S2Point> loopVertices = new ArrayList<>();
        loopVertices.add(center);
        loopVertices.add(regularVertices.get(tri * 2 + 0));
        loopVertices.add(regularVertices.get(tri * 2 + 1));
        loops.add(new S2Loop(loopVertices));
      }

      polygon = new S2Polygon(loops);
    }

    /**
     * Measures the average time to create an S2ContainsVertexQuery, add all the incoming and
     * outgoing edges of the polygon (w.r.t. the center point), and determine if the center point is
     * contained or not.
     */
    @Benchmark
    public int query() {
      S2ContainsVertexQuery query = new S2ContainsVertexQuery(center);

      for (S2Loop loop : polygon.getLoops()) {
        // Loop vertex 0 is the center point.
        query.addOutgoing(loop.vertex(1)); // The first loop edge is outgoing from the center.
        // The edge from vertex 1 to vertex 2 isn't incident on the center.
        query.addIncoming(loop.vertex(2)); // The third loop edge is incoming to the center.
      }

      return query.containsSign();
    }
  }


}
