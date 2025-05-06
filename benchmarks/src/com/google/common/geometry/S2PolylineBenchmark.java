/*
 * Copyright 2022 Google Inc.
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

import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polyline;
import com.google.common.geometry.TestDataGenerator.PointFactory;
import com.google.errorprone.annotations.CheckReturnValue;
import java.io.IOException;
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

/**
 * Benchmarks for {@link S2Polyline}.
 *
 * <p>TODO(user): Improve benchmarks for S2Polyline. Only a few methods are benchmarked.
 */
@CheckReturnValue
public class S2PolylineBenchmark {
  private S2PolylineBenchmark() {}

  /**
   * Benchmark state which provides a list of points with a fractal shape approximating an S2Cap
   * boundary, for constructing an S2Polyline. Also provides the polyline constructed from that
   * list, for testing other methods.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class SinglePolylineState extends S2BenchmarkBaseState {
    @Param({"8", "1024", "1048576"})
    int totalNumVertices;

    // A list of points which can be used to benchmark polyline construction.
    List<S2Point> testVertices;
    // A polyline constructed from testVertices which can be used to benchmark other methods.
    S2Polyline testPolyline;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        totalNumVertices = 8;
      }
      super.setup();

      S1Angle radius = S1Angle.degrees(1);
      S2Point center = data.getRandomPoint();
      testVertices =
          PointFactory.FRACTAL.createPoints(
              data, S2Cap.fromAxisAngle(center, radius), totalNumVertices);
      testPolyline = new S2Polyline(testVertices);
    }

    /** Measures the time to construct an S2Polyline. */
    @Benchmark
    public S2Polyline constructorFromList() {
      return new S2Polyline(testVertices);
    }

    /** Measures the time to compute getArclengthAngle. */
    @Benchmark
    public S1Angle getArclengthAngle() {
      return testPolyline.getArclengthAngle();
    }

    /** Measures the time to compute the polyline centroid. */
    @Benchmark
    public S2Point getCentroid() {
      return testPolyline.getCentroid();
    }
  }
}
