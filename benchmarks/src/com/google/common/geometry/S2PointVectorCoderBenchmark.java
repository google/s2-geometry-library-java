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

import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static java.util.concurrent.TimeUnit.MILLISECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2PointVectorCoder;
import com.google.common.geometry.S2Polygon;
import java.io.ByteArrayOutputStream;
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
import org.openjdk.jmh.infra.Blackhole;

/** Benchmarks for {@link S2PointVectorCoder}. */
public class S2PointVectorCoderBenchmark {
  private S2PointVectorCoderBenchmark() {}

  /**
   * An enum for which implementation of S2PointVectorCoder we are benchmarking.
   */
  public enum CoderType {
    FAST,
    COMPACT;

    public S2PointVectorCoder getCoder() {
      if (this == FAST) {
        return S2PointVectorCoder.FAST;
      } else {
        return S2PointVectorCoder.COMPACT;
      }
    }
  };

  /**
   * Benchmark the S2PointVectorCoder with a {@code List<S2Point>s}, created by encoding and then
   * decoding a randomly located regular loop. A CoderClass enum parameter switches between using
   * S2PointVectorCoder.FAST and S2PointVectorCoder.COMPACT. For both of those coders, decoding
   * actually happens only when the list contents are accessed, so the benchmarks here test the
   * speed of decoding by accessing the List elements.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MILLISECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class EncodedVectorState extends S2BenchmarkBaseState {
    @Param({"100000", "1000000"})
    int numVertices;

    @Param
    CoderType coderType;

    private List<S2Point> vector;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numVertices = 10000; // Unit tests just need to exercise the code paths, quickly.
      }
      super.setup();

      S2Loop loop = S2Loop.makeRegularLoop(data.getRandomPoint(), kmToAngle(5), numVertices);
      S2Polygon polygon = new S2Polygon();
      polygon.initToSnapped(new S2Polygon(loop), S2CellId.MAX_LEVEL);
      ByteArrayOutputStream bytesOutput = new ByteArrayOutputStream();
      coderType.getCoder().encode(polygon.loop(0).vertices(), bytesOutput);
      Bytes bytes = Bytes.fromByteArray(bytesOutput.toByteArray());
      vector = coderType.getCoder().decode(bytes, bytes.cursor());
    }

    @Benchmark
    public void getConsumer(Blackhole bh) {
      // Use a Blackhole to consume the data, otherwise the JVM might optimize it away. Normally,
      // loops should be avoided in benchmarks, but in this case it's realistic, as decoding will
      // typically be called in a tight loop that does something with the decoded points.
      // Note that this approach to obtaining the content of the vector benchmarked as slightly
      // faster than using a for loop.
      vector.forEach(bh::consume);
    }
  }
}
