/*
 * Copyright 2014 Google Inc.
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
package com.google.common.geometry;

import com.google.caliper.BeforeExperiment;
import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Polygon.S2PolygonIndex;

import java.util.List;

/**
 * Simple comparison of old and new indexes with a variety of loop and vertex counts.
 */
public class S2ShapeIndexBenchmark {
  @Param({"1", "32", "64"}) int numLoops;
  @Param({"4", "16", "256", "4096"}) int numVertices;

  private List<S2Loop> loops = Lists.newArrayList();
  private S2Polygon poly;

  @BeforeExperiment void setUp() {
    GeometryTestCase testUtils = new GeometryTestCase();
    testUtils.setUp();
    for (int k = 0; k < numLoops; k++) {
      loops.add(S2Loop.makeRegularLoop(
          testUtils.randomPoint(),
          GeometryTestCase.kmToAngle(5),
          numVertices));
    }
    poly = new S2Polygon(Lists.newArrayList(loops));
  }

  @Benchmark S2PolygonIndex oldIndex(int reps) {
    S2PolygonIndex idx = null;
    for (int i = 0; i < reps; i++) {
      idx = new S2PolygonIndex(poly);
      idx.computeIndex();
    }
    return idx;
  }

  @Benchmark S2ShapeIndex newIndex(int reps) {
    S2ShapeIndex idx = null;
    for (int i = 0; i < reps; i++) {
      idx = new S2ShapeIndex();
      for (S2Loop loop : loops) {
        idx.add(loop);
      }
      idx.iterator();
    }
    return idx;
  }
}
