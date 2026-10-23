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

import static com.google.common.geometry.S2TextFormat.makeLoop;
import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static com.google.common.geometry.TestDataGenerator.makeFractal;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2Coder;
import com.google.common.geometry.S2LatLng;
import com.google.common.geometry.S2LaxPolygonShape;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polygon;
import com.google.common.geometry.S2Shape;
import com.google.common.geometry.S2ShapeIndex;
import com.google.errorprone.annotations.CheckReturnValue;
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

/** Benchmarks for S2LaxPolygon. */
@CheckReturnValue
public class S2LaxPolygonBenchmark {
  private S2LaxPolygonBenchmark() {}

  private static final S2Point A = S2LatLng.fromDegrees(0, 0).toPoint();
  private static final S2Point B = S2LatLng.fromDegrees(0, 1).toPoint();
  private static final S2Point C = S2LatLng.fromDegrees(0, 2).toPoint();
  private static final S2Point D = S2LatLng.fromDegrees(1, 0).toPoint();
  private static final S2Point E = S2LatLng.fromDegrees(1, 1).toPoint();
  private static final S2Point F = S2LatLng.fromDegrees(1, 2).toPoint();
  private static final S2Point G = S2LatLng.fromDegrees(2, 0).toPoint();

  /** Several test geometries. */
  public static enum Geo {
    TRIANGLE(makeLoop(A, B, C)),
    SIMPLE(makeLoop(A, B, C, F, G, D)),
    COVERING_TRIANGLES(
        makeLoop(A, B, D),
        makeLoop(E, D, B),
        makeLoop(B, C, E),
        makeLoop(F, E, C),
        makeLoop(G, D, E)),
    FRACTAL_LOOP(makeFractal(20, 2000)),
    OVERLAPPING_FRACTALS(
        makeFractal(1, 100),
        makeFractal(2, 200),
        makeFractal(3, 300),
        makeFractal(4, 400),
        makeFractal(5, 500),
        makeFractal(6, 600));

    final ImmutableList<ImmutableList<S2Point>> vertices;
    final ImmutableList<ImmutableList<S2CellId>> cells;

    Geo(S2Loop... loops) {
      ImmutableList.Builder<ImmutableList<S2Point>> verts = ImmutableList.builder();
      ImmutableList.Builder<ImmutableList<S2CellId>> cells = ImmutableList.builder();

      for (S2Loop loop : loops) {
        List<S2Point> points = new ArrayList<>();
        List<S2CellId> ids = new ArrayList<>();
        for (S2Point v : loop.vertices()) {
          points.add(v);
          ids.add(S2CellId.fromPoint(v));
        }
        verts.add(ImmutableList.copyOf(points));
        cells.add(ImmutableList.copyOf(ids));
      }

      this.vertices = verts.build();
      this.cells = cells.build();
    }

    S2Polygon poly() {
      return new S2Polygon(new ArrayList<>(Lists.transform(vertices, S2Loop::new)));
    }
  }

  /** Several factories for S2Polygon.shape() and S2LaxPolygonShapes. */
  public static enum ShapeFactory {
    S2POLYGON() {
      @Override
      S2Polygon.Shape create(Geo geo) {
        return geo.poly().shape();
      }
    },
    LAX_ARRAY() {
      @Override
      S2Shape create(Geo geo) {
        return S2LaxPolygonShape.create(geo.vertices);
      }
    },
    LAX_PACKED() {
      @Override
      S2Shape create(Geo geo) {
        return S2LaxPolygonShape.createPacked(geo.vertices);
      }
    },
    LAX_SNAPPED() {
      @Override
      S2Shape create(Geo geo) {
        return S2LaxPolygonShape.createSnapped(geo.cells);
      }
    };

    abstract S2Shape create(Geo shape);
  }

  /**
   * Benchmark state which provides an S2Polygon and an S2Shape for each @Param Geometries, but only
   * the S2POLYGON ShapeFactory.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class PolygonShapeStates extends S2BenchmarkBaseState {
    @Param Geo geo;

    S2Polygon poly;
    S2Shape shape;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      ShapeFactory shapeFactory = ShapeFactory.S2POLYGON;
      shape = shapeFactory.create(geo);
      poly = geo.poly();
    }

    /** Measures the time to construct the building of the S2Polygon index. */
    @Benchmark
    public S2ShapeIndex indexS2Polygon() {
      S2ShapeIndex index = new S2ShapeIndex();
      poly.getLoops().forEach(index::add);
      Preconditions.checkState(index.iterator().atBegin());
      return index;
    }
  }

  /** Benchmark state which covers all combinations of the @Param Geometries and ShapeFactories. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class AllShapes extends S2BenchmarkBaseState {
    @Param Geo geo;
    @Param ShapeFactory shapeFactory;

    private S2Shape shape;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      shape = shapeFactory.create(geo);
    }

    /** Measures the time to construct each kind of shape. */
    @Benchmark
    public void create() {
      Preconditions.checkState(shapeFactory.create(geo).dimension() == 2);
    }

    /** Measures the time to build and index an S2ShapeIndex containing each shape. */
    @Benchmark
    public S2ShapeIndex indexShapes() {
      S2ShapeIndex index = new S2ShapeIndex();
      index.add(shape);
      Preconditions.checkState(index.iterator().atBegin());
      return index;
    }
  }

  /**
   * Constructs a polygon with the given number of concentric loops, each having the given number of
   * vertices. Concentric radii are scaled so all loops remain localized near the center point.
   */
  private static S2LaxPolygonShape makeMultiLoopPolygon(int numLoops, int verticesPerLoop) {
    S2Point center = S2LatLng.fromDegrees(23, 16).toPoint();
    List<List<S2Point>> loops = new ArrayList<>(numLoops);
    for (int i = 0; i < numLoops; ++i) {
      // Radii range from ~1 km (first loop) to 11 km (outermost loop).
      S1Angle radius = kmToAngle(1.0 + 10.0 * (i + 1.0) / numLoops);
      loops.add(S2Loop.makeRegularVertices(center, radius, verticesPerLoop));
    }
    return S2LaxPolygonShape.create(loops);
  }

  /**
   * Benchmarks {@link S2LaxPolygonShape} decoding with varying numbers of loops to measure the
   * overhead of loop offset decoding and validation.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class DecodeMultiLoop extends S2BenchmarkBaseState {
    /** An enum for which implementation of S2Coder we are benchmarking. */
    public enum CoderType {
      FAST,
      COMPACT;

      public S2Coder<S2LaxPolygonShape> getCoder() {
        return this == FAST ? S2LaxPolygonShape.FAST_CODER : S2LaxPolygonShape.COMPACT_CODER;
      }
    }

    @Param CoderType coderType;

    @Param({"4", "32"})
    int verticesPerLoop;

    @Param({"1", "10", "100", "1000", "10000"})
    int numLoops;

    private Bytes encodedData;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = Math.min(numLoops, 10);
      }
      super.setup();

      S2LaxPolygonShape polygon = makeMultiLoopPolygon(numLoops, verticesPerLoop);
      encodedData = Bytes.fromByteArray(coderType.getCoder().encode(polygon));
    }

    @Benchmark
    public S2LaxPolygonShape decode() throws IOException {
      return coderType.getCoder().decode(encodedData);
    }
  }
}
