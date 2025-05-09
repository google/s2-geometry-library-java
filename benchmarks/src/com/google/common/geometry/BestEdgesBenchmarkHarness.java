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

import static com.google.common.geometry.S2Projections.MAX_DIAG;
import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2Cell;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Shape;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex;
import com.google.common.geometry.S2ShapeIndexCoder;
import com.google.common.geometry.S2ShapeUtil;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;
import com.google.common.geometry.TestDataGenerator;
import com.google.common.geometry.TestDataGenerator.ShapeFactory;
import com.google.common.geometry.VectorCoder;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.BenchmarkParams;

/**
 * Infrastructure that is common to S2FurthestEdgeQueryBenchmark and S2ClosestEdgeQueryBenchmark.
 */
public final class BestEdgesBenchmarkHarness {
  private BestEdgesBenchmarkHarness() {}

  /** Targets may be of the following types. */
  public enum TargetType {
    POINT,
    EDGE,
    CELL,
    INDEX;
  }

  /** Factories may be of the following types. */
  public enum FactoryType {
    FRACTAL_LOOP,
    REGULAR_LOOP,
    POINT_CLOUD;
  }

  /**
   * A container for an index, a closest or furthest edge query on that index, and a set of targets
   * to test with.
   */
  public abstract static class BestEdgeQueryWithTargets {
    // The index may be replaced with an EncodedS2ShapeIndex.
    public S2ShapeIndex index = new S2ShapeIndex();

    // The number of edges in 'index'. Stored here to avoid repeated counting.
    public int numIndexEdges;

    // If the index is encoded, this is the list of shapes in it, which are encoded separately.
    public List<S2Shape> shapes;
  }

  /**
   * A base benchmark @State for measuring the performance of S2ClosestEdgeQuery.findClosestEdges()
   * and S2FurthestEdgeQuery.findFurthestEdges() under a variety of scenarios.
   *
   * <p>The setup() method builds NUM_INDEX_SAMPLES different S2ShapeIndexes, each with
   * approximately "numIndexEdges" edges generated by a specified factory. Each S2ShapeIndex
   * geometry is generated within an S2Cap of the radius specified by DEFAULT_INDEX_RADIUS_KM, the
   * "index radius".
   *
   * <p>A S2ClosestEdgeQuery or S2FurthestEdgeQuery is constructed for each of those S2ShapeIndexes,
   * with options set according to the parameters to setup(). If maxDistanceFraction > 0, then
   * maxDistance() is set to the given fraction of the index radius. If maxErrorFraction > 0, then
   * maxError() is set to the given fraction of the index radius.
   *
   * <p>For each query, NUM_TARGET_SAMPLES different targets will also be built, of the given {@code
   * TargetType}. If "targetType" is INDEX, then the target will have approximately "numTargetEdges"
   * edges. The query and its targets are stored in either a ClosestEdgeQueryWithTargets or
   * FurthestEdgeQueryWithTargets, built by the {@code S2EdgeQueryGenerator}, which generates the
   * geometry for the index and targets according to the parameters passed to it, some of which are
   * described above. See that class for explanation of the other parameters.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public abstract static class BestEdgesBenchmarkBaseState extends S2BenchmarkBaseState {
    // Default radius for indexed geometry, 1000 km.
    protected static final double DEFAULT_INDEX_RADIUS_KM = 1000;

    // Several different random S2ShapeIndexes are built to run queries against. (Performance is
    // affected by how the shapes are positioned relative to the S2Cell hierarchy.)
    protected static final int DEFAULT_NUM_INDEX_SAMPLES = 8;

    // How many targets will be tested for each query index?
    protected static final int DEFAULT_NUM_TARGET_SAMPLES = 16;

    /**
     * Must be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Specifies the {@link
     * TestDataGenerator.ShapeIndexFactory} to use to generate the index geometry. The index
     * geometry is always bounded by a cap of size {@link indexRadiusKm}. No default is provided.
     */
    protected FactoryType factoryType = null;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. True if the index
     * should be encoded with S2ShapeIndexCoder. Defaults to false.
     */
    protected boolean encodeIndex = false;

    /**
     * Should be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. The approximate
     * number of edges to be added to the index. The actual value may depend on the FactoryType.
     */
    protected int numIndexEdges = 1024;

    /**
     * Should be set by benchmarks inheriting from BestEdgesBenchmarkBaseState unless only
     * findClosestEdge() / findFurthestEdge() is used in the benchmark. The number of results to
     * obtain.
     */
    protected int maxResults = 1;

    /**
     * Should be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Determines if
     * interiors of indexed polygons be considered when measuring distance? Defaults to false.
     */
    protected boolean includeInteriors = false;

    /**
     * Must be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. The type of targets to
     * create.
     */
    protected TargetType targetType = null;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Specifies the number of
     * edges to use if targetType is INDEX, and ignored otherwise. Defaults to 4.
     */
    protected int numTargetEdges = 4;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. If positive, determines
     * the maxDistance to be used for S2ClosestEdgeQuery.Options, or minDistance used for the
     * S2FurthestEdgeQuery.Options as a fraction of {@link indexRadiusKm}.
     *
     * <p>A positive value indicates a fixed value, i.e. 0.5 for 500 km, while a negative value
     * indicates that no maximum distance (for closest edge queries) or minimum distance (for
     * furthest edge queries) will be set.
     *
     * <p>Defaults to -1, for no limits on distance.
     */
    protected double distanceFraction = -1;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Determines the maxError
     * to be used for the query Options as a fraction of {@link indexRadiusKm}.
     *
     * <p>A positive value indicates a fixed value, i.e. 0.5 for 500 km. A negative value, which is
     * the default, specifies that maxError will not be set, so will be zero.
     */
    protected double maxErrorFraction = -1;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. If true, the target
     * will be created by sampling edges from the index. Otherwise (the default) the target will be
     * created by constructing a shape within the target cap.
     */
    protected boolean chooseTargetFromIndex = false;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Determines the size of
     * the target cap, as a fraction of the {@link indexRadiusKm}.
     *
     * <p>A positive value indicates a fixed value, i.e. 0.5 for 500 km, while a negative value
     * indicates the maximum for a uniform distribution of values starting at zero. For example,
     * -0.5 indicates a uniformly distributed value between 0 and 500 km.
     *
     * <p>Defaults to -1, for uniformly distributed values between 0 and {@link indexRadiusKm}.
     */
    protected double targetRadiusFraction = -1;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Determines the distance
     * between the center of the S2Cap bounding the index and the center of the S2Cap bounding the
     * target.
     *
     * <p>A positive value indicates a fixed value, i.e. 0.5 for 500 km, while a negative value
     * indicates the maximum for a uniform distribution of values starting at zero. For example,
     * -0.5 indicates a uniformly distributed value between 0 and 500 km. Note that if
     * centerSeparationFraction is less than (0.5 * (1 + targetRadiusFraction)) then the S2Cap
     * bounding the target and the S2Cap bounding the index will overlap, and edges of the target
     * and index might intersect, or a polygon in the index might completely contain a connected
     * component of the target, or vice versa.
     *
     * <p>Defaults to -1, for uniformly distributed values between 0 and {@link indexRadiusKm}.
     */
    protected double centerSeparationFraction = -1;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. Determines the size of
     * generated indexes and other parameters above which are specified as fractions of this value.
     * Defaults to DEFAULT_INDEX_RADIUS_KM.
     */
    protected double indexRadiusKm = DEFAULT_INDEX_RADIUS_KM;

    /** The indexRadiusKm converted to an angle, which happens in initialize(). */
    protected S1Angle indexRadius;

    /**
     * May be set by benchmarks inheriting from BestEdgesBenchmarkBaseState. If true,
     * Options.setUseBruteForce(true) will be set. Defaults to false.
     */
    protected boolean useBruteForce = false;

    // An instance of ShapeFactory specified by the factoryType.
    protected ShapeFactory factory;

    // Current position as the benchmark iterates through the targets and queries.
    protected int queryIndex = 0;
    protected int targetIndex = 0;

    // The actual number of index samples to use.
    protected int numIndexSamples = DEFAULT_NUM_INDEX_SAMPLES;
    // The actual number of target samples to use.
    protected int numTargetSamples = DEFAULT_NUM_TARGET_SAMPLES;

    // Subclasses implement this to update the current query from the value of queryIndex.
    abstract void setCurrentQuery();

    // Subclasses implement this to update the current target from the values of queryIndex and
    // targetIndex.
    abstract void setCurrentTarget();

    // Subclasses implement this to create the list of queries with their indexes and targets.
    abstract void createQueriesAndTargets() throws IOException;

    /**
     * If "fraction" is positive, returns the given "radius" multiplied by "fraction". If "fraction"
     * is negative, returns a random S1Angle uniformly distributed between zero and a positive
     * maximum of "radius" * abs("fraction").
     */
    protected S1Angle fractionToRadiusOrRange(double fraction, S1Angle radius) {
      if (fraction < 0) {
        fraction = -(data.nextDouble() * fraction);
      }
      return radius.mul(fraction);
    }

    /**
     * Must be called by benchmarks inheriting from BestEdgesBenchmarkBaseState in their JMH @Setup.
     * Initializes the edgeQueriesAndTargets according to the current values that have been set for
     * the index and target generation and prepares to iterate through query indexes and targets.
     */
    protected void initialize(BenchmarkParams params) throws IOException {
      // Determine if we are running as part of a unit test.
      setUnitTestModeFromParams(params);

      // For unit tests, override setup parameters for speed.
      if (isUnitTest()) {
        numIndexEdges = 4;
        numTargetEdges = 4;
        numIndexSamples = 1;
        numTargetSamples = 1;
      }
      super.setup();
      indexRadius = kmToAngle(indexRadiusKm);

      // Create the factory from the factoryType.
      switch (factoryType) {
        case FRACTAL_LOOP:
          this.factory = data.new FractalLoopShapeFactory();
          break;
        case REGULAR_LOOP:
          this.factory = data.new RegularLoopShapeFactory();
          break;
        case POINT_CLOUD:
          this.factory = data.new PointCloudShapeFactory();
          break;
      }

      // Create indexes, with closest or furthest edge queries on them, and collections of closest
      // or furthest targets. This method is provided by subclasses.
      createQueriesAndTargets();

      // Reset iteration over targets and queries.
      queryIndex = 0;
      targetIndex = 0;
      setCurrentQuery();
      setCurrentTarget();
    }

    /** Generates an index according to the current specification, and returns the bounding cap. */
    protected S2Cap fillOutputIndex(BestEdgeQueryWithTargets output, boolean encodeIndex)
        throws IOException {
      // Use the factory to add edges to the index.
      S2Cap indexCap = S2Cap.fromAxisAngle(data.getRandomPoint(), kmToAngle(indexRadiusKm));
      output.index = factory.getShape(indexCap, numIndexEdges);

      // Count the number of edges added when done.
      output.numIndexEdges = S2ShapeUtil.countEdges(output.index);

      // Encode the index if required.
      if (encodeIndex) {
        ByteArrayOutputStream stream = new ByteArrayOutputStream();
        VectorCoder.FAST_SHAPE.encode(output.index.getShapes(), stream);
        S2ShapeIndexCoder.INSTANCE.encode(output.index, stream);

        Bytes data = Bytes.fromByteArray(stream.toByteArray());
        Cursor cursor = data.cursor();

        // Store the encoded index with its shapes in the output, replacing the un-encoded index.
        output.shapes = VectorCoder.FAST_SHAPE.decode(data, cursor);
        output.index = new S2ShapeIndexCoder(output.shapes).decode(data, cursor);
      }
      return indexCap;
    }

    protected S2Cap getTargetCap(S2Cap indexCap) {
      S1Angle targetDist = fractionToRadiusOrRange(centerSeparationFraction, indexRadius);
      return S2Cap.fromAxisAngle(
          data.sampleBoundary(S2Cap.fromAxisAngle(indexCap.axis(), targetDist)),
          fractionToRadiusOrRange(targetRadiusFraction, indexRadius));
    }

    protected S2Point getTargetPoint(BestEdgeQueryWithTargets output, S2Cap targetCap) {
      return chooseTargetFromIndex
          ? data.sampleEdge(output.index, output.numIndexEdges).getStart()
          : targetCap.axis();
    }

    protected MutableEdge getTargetEdge(BestEdgeQueryWithTargets output, S2Cap targetCap) {
      if (chooseTargetFromIndex) {
        return data.sampleEdge(output.index, output.numIndexEdges);
      }
      MutableEdge edge = new MutableEdge();
      edge.set(data.sampleBoundary(targetCap), data.sampleBoundary(targetCap));
      return edge;
    }

    protected S2Cell getTargetCell(BestEdgeQueryWithTargets output, S2Cap targetCap) {
      S2CellId cellId =
          chooseTargetFromIndex
              ? data.sampleCell(output.index)
              : S2CellId.fromPoint(targetCap.axis())
                  .parent(MAX_DIAG.getClosestLevel(targetCap.radius().radians()));
      return new S2Cell(cellId);
    }

    protected S2ShapeIndex getTargetIndex(BestEdgeQueryWithTargets output, S2Cap targetCap) {
      S2ShapeIndex targetIndex = new S2ShapeIndex();
      if (chooseTargetFromIndex) {
        // If the index factoryType is POINT_CLOUD the edges are points, and can't be added to
        // a S2EdgeVectorShape, because it checks that the edges are not degenerate. In that case,
        // use an S2Point.Shape.
        S2Shape shape;
        if (factoryType == FactoryType.POINT_CLOUD) {
          List<S2Point> points = new ArrayList<>();
          for (int p = 0; p < numTargetEdges; ++p) {
            MutableEdge edge = data.sampleEdge(output.index, output.numIndexEdges);
            points.add(edge.getStart());
          }
          shape = S2Point.Shape.fromList(points);
        } else {
          S2EdgeVectorShape edgeVectorshape = new S2EdgeVectorShape();
          for (int e = 0; e < numTargetEdges; ++e) {
            MutableEdge edge = data.sampleEdge(output.index, output.numIndexEdges);
            edgeVectorshape.add(edge.getStart(), edge.getEnd());
          }
          shape = edgeVectorshape;
        }
        targetIndex.add(shape);
      } else {
        targetIndex = factory.getShape(targetCap, numTargetEdges);
      }

      targetIndex.applyUpdates();
      return targetIndex;
    }

    /**
     * Iterate through the targets. When all the targets for the current query are done, go to the
     * next query. Loop back to the start at the end of the queries.
     */
    public void next() {
      if (++targetIndex == numTargetSamples) {
        targetIndex = 0;
        if (++queryIndex == numIndexSamples) {
          queryIndex = 0;
        }
        setCurrentQuery();
      }
      setCurrentTarget();
    }
  }
}
