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

import static com.google.common.geometry.TestDataGenerator.getBoundingPolygon;
import static com.google.common.geometry.TestDataGenerator.index;

import com.google.common.geometry.S2Cell;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polygon;
import com.google.common.geometry.S2RegionCoverer;
import com.google.common.geometry.TestDataGenerator;
import com.google.common.geometry.TestDataGenerator.PolygonFactory;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.infra.BenchmarkParams;

/**
 * A base class for JMH @State in S2 benchmarks which provides utilities for building complex
 * benchmark states and test data.
 */
public class S2BenchmarkBaseState {
  protected TestDataGenerator data;

  /** Is this benchmark run actually a unit test? */
  private boolean isUnitTest = false;

  /**
   * Reinitializes the TestDataGenerator, which reinitializes the contained Random to a known state.
   * This is equivalent to calling data.resetSeed().
   */
  public void setup() {
    data = new TestDataGenerator();
  }

  /**
   * Sets the "unit test mode" for this benchmark state if the provided BenchmarkParams indicate
   * that the current benchmark run is actually a unit test. This can be used to disable expensive
   * benchmarks, or adjust parameters as required so that dry run unit tests complete quickly.
   *
   * <p>If there are no warmup iterations, no forks, and the benchmark mode is single-shot, this is
   * a unit test, as real benchmarking does not use that combination of parameters.
   */
  public void setUnitTestModeFromParams(BenchmarkParams params) {
    isUnitTest =
        (params.getWarmup().getCount() == 0)
            && (params.getForks() == 0)
            && (params.getMode() == Mode.SingleShotTime);
  }

  /**
   * Returns true iff this benchmark run is actually a unit test, as set in a previous call to
   * {@link setUnitTestModeFromParams }.
   */
  public boolean isUnitTest() {
    return isUnitTest;
  }

  /**
   * Abstract benchmark state that cycles through a list of polygons from a provided factory, and a
   * list of query points for each of those polygons which are sampled from the bound of the
   * polygon. This is intended for benchmarking containment of points in polygons.
   */
  public abstract static class ContainsPointBaseState extends S2BenchmarkBaseState {
    protected static final int NUM_QUERY_SAMPLES = 10;
    protected static final int NUM_POLYGONS = 5;

    protected List<S2Polygon> polygons = new ArrayList<>(NUM_POLYGONS);
    protected List<List<S2Point>> queries = new ArrayList<>(NUM_POLYGONS);
    protected int polygonIndex;
    protected int queryIndex;

    // The current polygon for checking containment against.
    protected S2Polygon currentPolygon;

    protected void setupPolygonsAndQueries(
        TestDataGenerator data,
        PolygonFactory factory,
        int numLoops,
        int totalNumVertices,
        boolean preindexOption)
        throws IOException {
      polygons.clear();
      queries.clear();

      for (int i = 0; i < NUM_POLYGONS; ++i) {
        S2Polygon polygon = factory.newPolygon(data, numLoops, totalNumVertices);
        if (preindexOption) {
          polygons.add(index(polygon));
        } else {
          polygons.add(polygon);
        }

        List<S2Point> polygonSpecificQueries = new ArrayList<>(NUM_QUERY_SAMPLES);
        for (int j = 0; j < NUM_QUERY_SAMPLES; ++j) {
          polygonSpecificQueries.add(data.samplePoint(polygon.getRectBound()));
        }
        queries.add(polygonSpecificQueries);
      }

      polygonIndex = 0;
      queryIndex = 0;
      currentPolygon = polygons.get(polygonIndex);
    }

    /**
     * Benchmarks should call advanceQuery(), which will be included in the benchmark cost. When
     * queries for the current polygon are done, automatically advances to the next polygon.
     */
    protected void advanceQuery() {
      queryIndex = (1 + queryIndex) % NUM_QUERY_SAMPLES;

      // When queries loop around, advance to the next polygon.
      if (queryIndex == 0) {
        polygonIndex = (1 + polygonIndex) % NUM_POLYGONS;
        currentPolygon = polygons.get(polygonIndex);
      }
    }

    /** Return the current query for the current polygon. */
    protected S2Point currentQuery() {
      return queries.get(polygonIndex).get(queryIndex);
    }
  }

  /**
   * Abstract benchmark state that cycles through a list of polygons from a provided polygon
   * factory. The state also provides a copy of each polygon, and the bound of each polygon for
   * checking intersection, containment, etc. Optionally, the polygons can be pre-indexed.
   */
  public abstract static class PolygonListState extends S2BenchmarkBaseState {
    protected static final int NUM_POLYGONS = 10;
    protected ArrayList<S2Polygon> polygons = new ArrayList<>(NUM_POLYGONS);
    protected ArrayList<S2Polygon> copies = new ArrayList<>(NUM_POLYGONS);
    protected ArrayList<S2Polygon> bounds = new ArrayList<>(NUM_POLYGONS);
    protected int polygonIndex;

    protected void setupPolygons(
        TestDataGenerator data,
        PolygonFactory factory,
        int numLoops,
        int totalNumVertices,
        boolean preindexOption)
        throws IOException {
      super.setup();

      polygons.clear();
      copies.clear();
      bounds.clear();

      for (int i = 0; i < NUM_POLYGONS; ++i) {
        S2Polygon testPolygon = factory.newPolygon(data, numLoops, totalNumVertices);
        if (preindexOption) {
          polygons.add(index(testPolygon));
          copies.add(index(testPolygon));
          bounds.add(index(getBoundingPolygon(testPolygon)));
        } else {
          polygons.add(testPolygon);
          copies.add(testPolygon);
          bounds.add(getBoundingPolygon(testPolygon));
        }
      }
      polygonIndex = 0;
      currentPolygon = polygons.get(polygonIndex);
      currentPolygonCopy = copies.get(polygonIndex);
      currentPolygonBound = bounds.get(polygonIndex);
    }

    // Benchmarks can use the following.
    protected S2Polygon currentPolygon;
    protected S2Polygon currentPolygonCopy;
    protected S2Polygon currentPolygonBound;

    // Benchmarks should proceed around the cycle of test polygons by calling advancePolygon. The
    // cost of advancePolygon() is therefore included in benchmark times.
    protected void advancePolygon() {
      polygonIndex = (1 + polygonIndex) % NUM_POLYGONS;
      currentPolygon = polygons.get(polygonIndex);
      currentPolygonCopy = copies.get(polygonIndex);
      currentPolygonBound = bounds.get(polygonIndex);
    }
  }

  /**
   * Abstract benchmark state that extends PolygonCircularListState to add the polygon complements,
   * where the complements are obtained by subtracting the polygon from the bound of the polygon.
   */
  public abstract static class PolygonAndComplementListState extends PolygonListState {
    protected List<S2Polygon> complements = new ArrayList<>(NUM_POLYGONS);
    ;
    protected S2Polygon currentPolygonComplement;

    @Override
    protected void setupPolygons(
        TestDataGenerator data,
        PolygonFactory factory,
        int numLoops,
        int totalNumVertices,
        boolean preindexOption)
        throws IOException {
      super.setupPolygons(data, factory, numLoops, totalNumVertices, preindexOption);

      // Build the list of complement polygons to match the test polygons.
      complements.clear();
      for (int i = 0; i < NUM_POLYGONS; i++) {
        S2Polygon complement = new S2Polygon();
        complement.initToDifference(bounds.get(i), polygons.get(i));
        if (preindexOption) {
          complements.add(index(complement));
        } else {
          complements.add(complement);
        }
      }
      currentPolygonComplement = complements.get(polygonIndex);
    }

    @Override
    protected void advancePolygon() {
      super.advancePolygon();
      currentPolygonComplement = complements.get(polygonIndex);
    }
  }

  /**
   * Abstract benchmark state that cycles through a list of polygons and for each polygon, a list of
   * covering cell ids.
   */
  public abstract static class PolygonAndCoveringListState extends PolygonListState {
    protected ArrayList<ArrayList<S2CellId>> coverings;
    // The current cell of the current polygon's covering.
    protected int cellIndex;

    protected void setupPolygonsAndCoverings(
        TestDataGenerator data,
        PolygonFactory factory,
        int numLoops,
        int totalNumVertices,
        boolean preindexOption,
        int maxCells)
        throws IOException {
      super.setupPolygons(data, factory, numLoops, totalNumVertices, preindexOption);

      // Build the list of coverings to match the test polygons.
      S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(maxCells).build();
      coverings = new ArrayList<>(NUM_POLYGONS);
      for (int i = 0; i < NUM_POLYGONS; i++) {
        S2Polygon testPolygon = polygons.get(i);
        ArrayList<S2CellId> covering = new ArrayList<>(16);
        coverer.getCovering(testPolygon, covering);
        coverings.add(covering);
      }
      cellIndex = 0;
    }

    /** Returns the current S2CellId. */
    protected S2CellId currentCellId() {
      return coverings.get(polygonIndex).get(cellIndex);
    }

    /** Construct and index a new S2Polygon for the current cell of the current covering. */
    protected S2Polygon constructIndexedPolygonForCurrentCell() {
      return index(new S2Polygon(new S2Cell(currentCellId())));
    }

    /*
     * Advance to the next cell in the current covering. When the current covering is complete,
     * automatically advances to the next polygon.
     */
    protected void advanceCovering() {
      cellIndex = (1 + cellIndex) % coverings.get(polygonIndex).size();
      if (cellIndex == 0) {
        advancePolygon();
      }
    }
  }
}
