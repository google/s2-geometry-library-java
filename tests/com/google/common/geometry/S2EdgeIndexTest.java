/*
 * Copyright 2011 Google Inc.
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

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import java.util.HashSet;
import java.util.List;
import java.util.logging.Logger;

/**
 * Tests for {@link S2EdgeIndex}.
 *
 * @author andriy@google.com (Andriy Bihun) ported from util/geometry
 * @author pilloff@google.com (Mark Pilloff) original author
 */
public strictfp class S2EdgeIndexTest extends GeometryTestCase {
  private static final Logger log = Logger.getLogger(S2EdgeIndexTest.class.getCanonicalName());

  public static class EdgeVectorIndex extends S2EdgeIndex {
    private List<S2Edge> edges;

    public EdgeVectorIndex(List<S2Edge> edges) {
      this.edges = edges;
    }

    @Override
    protected int getNumEdges() {
      return edges.size();
    }

    @Override
    protected S2Point edgeFrom(int index) {
      return edges.get(index).getStart();
    }

    @Override
    protected S2Point edgeTo(int index) {
      return edges.get(index).getEnd();
    }
  }

  /**
   * Generates a random edge whose center is in the given cap.
   */
  private S2Edge randomEdgeCrossingCap(double maxLengthMeters, S2Cap cap) {
    // Pick the edge center at random.
    S2Point edgeCenter = samplePoint(cap);
    // Pick two random points in a suitably sized cap about the edge center.
    S2Cap edgeCap = S2Cap.fromAxisAngle(
        edgeCenter, S1Angle.radians(maxLengthMeters / S2LatLng.EARTH_RADIUS_METERS / 2));
    S2Point p1 = samplePoint(edgeCap);
    S2Point p2 = samplePoint(edgeCap);
    return new S2Edge(p1, p2);
  }

  /*
   * Generates "numEdges" random edges, of length at most "edgeLengthMetersMax"
   * and each of whose center is in a randomly located cap with radius
   * "capSpanMeters", and puts results into "edges".
   */
  private void generateRandomEarthEdges(
      double edgeLengthMetersMax, double capSpanMeters, int numEdges, List<S2Edge> edges) {
    S2Cap cap = S2Cap.fromAxisAngle(
        randomPoint(), S1Angle.radians(capSpanMeters / S2LatLng.EARTH_RADIUS_METERS));
    for (int i = 0; i < numEdges; ++i) {
      edges.add(randomEdgeCrossingCap(edgeLengthMetersMax, cap));
    }
  }

  private void checkAllCrossings(
      List<S2Edge> allEdges, int minCrossings, int maxChecksCrossingsRatio) {
    EdgeVectorIndex index = new EdgeVectorIndex(allEdges);
    index.computeIndex();
    EdgeVectorIndex.DataEdgeIterator it = new EdgeVectorIndex.DataEdgeIterator(index);
    double totalCrossings = 0;
    double totalIndexChecks = 0;

    for (int in = 0; in < allEdges.size(); ++in) {
      S2Edge e = allEdges.get(in);

      HashSet<Integer> candidateSet = Sets.newHashSet();

      StringBuilder sb = new StringBuilder();
      for (it.getCandidates(e.getStart(), e.getEnd()); it.hasNext(); it.next()) {
        candidateSet.add(it.index());
        sb.append(it.index()).append("/");
        ++totalIndexChecks;
      }

      for (int i = 0; i < allEdges.size(); ++i) {
        int crossing = S2EdgeUtil.robustCrossing(
            e.getStart(), e.getEnd(), allEdges.get(i).getStart(), allEdges.get(i).getEnd());
        if (crossing >= 0) {
          StringBuilder sbError = new StringBuilder();
          sbError
              .append("\n==CHECK_ERROR===================================\n")
              .append("CandidateSet: ")
              .append(sb)
              .append("\nin=")
              .append(in)
              .append(" i=")
              .append(i)
              .append(" robustCrossing=")
              .append(crossing)
              .append("\nfrom:\n")
              .append(e)
              .append("\nto:\n")
              .append(allEdges.get(i))
              .append("\n==================================================");
          assertTrue(sbError.toString(), candidateSet.contains(i));
          ++totalCrossings;
        }
      }
    }

    log.info(
        "Pairs/num crossings/check crossing ratio: "
            + Integer.toString(allEdges.size() * allEdges.size()) + "/"
            + Double.toString(totalCrossings) + "/"
            + Double.toString(totalIndexChecks / totalCrossings));
    assertTrue(minCrossings <= totalCrossings);
    assertTrue(totalCrossings * maxChecksCrossingsRatio >= totalIndexChecks);
  }

  /*
   * Generates random edges and tests, for each edge, that all those that cross
   * are candidates.
   */
  private void tryCrossingsRandomInCap(int numEdges, double edgeLengthMax, double capSpanMeters,
      int minCrossings, int maxChecksCrossingsRatio) {
    List<S2Edge> allEdges = Lists.newArrayList();
    generateRandomEarthEdges(edgeLengthMax, capSpanMeters, numEdges, allEdges);
    checkAllCrossings(allEdges, minCrossings, maxChecksCrossingsRatio);
  }

  public void testSpecificEdges() {
    List<S2Point> ps = Lists.newArrayList();
    ps.add(new S2Point(0.8088625416501157, -0.40633615485481134, 0.4250086092929434));
    ps.add(new S2Point(0.8088939911085784, -0.40631384442755236, 0.4249700824469155));
    ps.add(new S2Point(0.8088088971141814, -0.40642839367135375, 0.425022503835579));
    ps.add(new S2Point(0.8088643962606756, -0.406333410696549, 0.4250077032402616));
    List<S2Edge> allEdges = Lists.newArrayList();
    allEdges.add(new S2Edge(ps.get(0), ps.get(1)));
    allEdges.add(new S2Edge(ps.get(2), ps.get(3)));
    checkAllCrossings(allEdges, 0, 16);
  }

  public void testLoopCandidateOfItself() {
    List<S2Point> ps = Lists.newArrayList(); // A diamond loop around 0,180.
    ps.add(makePoint("0:178"));
    ps.add(makePoint("-1:180"));
    ps.add(makePoint("0:-179"));
    ps.add(makePoint("1:-180"));
    List<S2Edge> allEdges = Lists.newArrayList();
    for (int i = 0; i < 4; ++i) {
      allEdges.add(new S2Edge(ps.get(i), ps.get((i + 1) % 4)));
    }
    checkAllCrossings(allEdges, 0, 16);
  }

  public void testRandomEdgeCrossings() {
    tryCrossingsRandomInCap(2000, 30, 5000, 500, 2);
    tryCrossingsRandomInCap(1000, 100, 5000, 500, 3);
    tryCrossingsRandomInCap(1000, 1000, 5000, 1000, 40);
    tryCrossingsRandomInCap(500, 5000, 5000, 5000, 20);
  }

  public void testRandomEdgeCrossingsSparse() {
    for (int i = 0; i < 5; ++i) {
      tryCrossingsRandomInCap(2000, 100, 5000, 500, 8);
      tryCrossingsRandomInCap(2000, 300, 50000, 1000, 10);
    }
  }
}
