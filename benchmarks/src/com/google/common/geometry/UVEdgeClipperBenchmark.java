/*
 * Copyright 2024 Google Inc.
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

import com.google.common.collect.Iterators;
import com.google.common.geometry.R2Rect;
import com.google.common.geometry.S2Cell;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2Edge;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Projections;
import com.google.common.geometry.UVEdgeClipper;
import java.util.Iterator;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

@State(Scope.Thread)
public class UVEdgeClipperBenchmark extends S2BenchmarkBaseState {
  private final S2CellId cell = S2CellId.fromToken("89c25c1");
  private final UVEdgeClipper clipper = new UVEdgeClipper();
  {
    clipper.init(new S2Cell(cell));
  }

  private Iterator<S2Edge> it;

  /** The level of cell above 'cell' to generate edges in; hit rate should be ~4^-parentLevel. */
  @Param({"0", "2", "12"})
  int level;

  /** Whether the endpoints of generated edges should be connected. */
  @Param({"false", "true"})
  boolean connected;

  @Override
  @Setup(Level.Trial)
  public void setup() {
    super.setup();
    // Build edges on one face some level above 'cell', connected iff requested.
    int face = cell.face();
    R2Rect bound = cell.parent(level).getBoundUV();
    S2Edge[] edges = new S2Edge[100_000];
    if (connected) {
      S2Point last = randomEndpoint(face, bound);
      for (int i = 0; i < edges.length; i++) {
        S2Point next = randomEndpoint(face, bound);
        edges[i] = new S2Edge(last, next);
        last = next;
      }
    } else {
      for (int i = 0; i < edges.length; i++) {
        edges[i] = new S2Edge(randomEndpoint(face, bound), randomEndpoint(face, bound));
      }
    }
    it = Iterators.cycle(edges);
  }

  private S2Point randomEndpoint(int face, R2Rect bound) {
    return S2Projections.faceUvToXyz(face,
            data.uniform(bound.x().lo(), bound.x().hi()),
            data.uniform(bound.y().lo(), bound.y().hi())).normalize();
  }

  /** Benchmarks clipping edges in UV space that may be connected in sequence. */
  @Benchmark
  @Warmup(iterations = 1000, time = 16, timeUnit = MICROSECONDS)
  @Measurement(iterations = 1000, time = 16, timeUnit = MICROSECONDS)
  public void clipEdge(Blackhole bh) {
    S2Edge edge = it.next();
    bh.consume(clipper.clipEdge(edge.getStart(), edge.getEnd(), connected));
  }
}
