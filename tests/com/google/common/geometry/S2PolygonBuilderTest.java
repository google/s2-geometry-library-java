/*
 * Copyright 2006 Google Inc.
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

import java.util.List;
import java.util.logging.Logger;

/**
 * Tests for {@link S2Loop}.
 *
 */
public strictfp class S2PolygonBuilderTest extends GeometryTestCase {
  private static final Logger log = Logger.getLogger(S2PolygonBuilderTest.class.getCanonicalName());

  // A chain represents either a polyline or a loop, depending
  // on whether "closed" is true.
  private class Chain {
    String str;
    boolean closed;

    public Chain(String str, boolean closed) {
      this.str = str;
      this.closed = closed;
    }
  }

  private class TestCase {
    // +1 = undirected, -1 = directed, 0 = either one
    int undirectedEdges;

    // +1 = XOR, -1 = don't XOR, 0 = either one
    int xorEdges;

    // Minimum and maximum merge distances for this test case in degrees.
    double minMerge;
    double maxMerge;

    // Each test case consists of a set of input loops and polylines.
    Chain[] chainsIn;

    // The expected set of output loops, directed appropriately.
    String[] loopsOut;

    // The expected number of unused edges.
    int numUnusedEdges;

    public TestCase(int undirectedEdges,
        int xorEdges,
        double minMerge,
        double maxMerge,
        Chain[] chainsIn,
        String[] loopsOut,
        int numUnusedEdges) {
      this.undirectedEdges = undirectedEdges;
      this.xorEdges = xorEdges;
      this.minMerge = minMerge;
      this.maxMerge = maxMerge;
      this.chainsIn = chainsIn;
      this.loopsOut = loopsOut;
      this.numUnusedEdges = numUnusedEdges;
    }
  }

  TestCase[] testCases = new TestCase[] {
  // 0: No loops.
  new TestCase(0, 0, 0.0, 10.0, new Chain[] {new Chain(null, false)}, new String[] {}, 0),

      // 1: One loop with some extra edges.
      new TestCase(0,
          0,
          0.0,
          4.0,
          new Chain[] {new Chain("0:0, 0:10, 10:5", true), new Chain("0:0, 5:5", false),
              new Chain("10:5, 20:7, 30:10, 40:15, 50:3, 60:-20", false)},
          new String[] {"0:0, 0:10, 10:5"},
          6),

      // 2: One loop that has an edge removed by XORing, plus lots of
      // extra edges.
      new TestCase(0, 1, 0.0, 1.0, // XOR
          new Chain[] {new Chain("0:0, 0:10, 5:15, 10:10, 10:0", true),
              new Chain("10:10, 12:12, 14:14, 16:16, 18:18", false),
              new Chain("14:14, 14:16, 14:18, 14:20", false),
              new Chain("14:18, 16:20, 18:22", false),
              new Chain("18:12, 16:12, 14:12, 12:12", false),
              new Chain("20:18, 18:16, 16:14, 14:12", false),
              new Chain("20:14, 18:14, 16:14", false),
              new Chain("5:15, 0:10", false)},
          new String[] {},
          21),

      // 3: Three loops (two shells and one hole) that combine into one.
      new TestCase(0, 1, 0.0, 4.0, // XOR
          new Chain[] {new Chain("0:0, 0:10, 5:10, 10:10, 10:5, 10:0", true),
              new Chain("0:10, 0:15, 5:15, 5:10", true),
              new Chain("10:10, 5:10, 5:5, 10:5", true), },
          new String[] {"0:0, 0:10, 0:15, 5:15, 5:10, 5:5, 10:5, 10:0"},
          0),

      // 4: A big CCW triangle contained 3 CW triangular holes. The whole thing
      // looks like a pyramid of nine small triangles (with two extra edges).
      new TestCase(-1, 0, 0.0, 0.9, // Directed edges required for unique result.
          new Chain[] {new Chain("0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1", true),
              new Chain("0:2, 1:1, 1:3", true),
              new Chain("0:4, 1:3, 1:5", true),
              new Chain("1:3, 2:2, 2:4", true),
              new Chain("0:0, 0:1", false),
              new Chain("1:3, 5:7", false)},
          new String[] {"0:0, 0:2, 1:1",
              "0:2, 0:4, 1:3",
              "0:4, 0:6, 1:5",
              "1:1, 1:3, 2:2",
              "1:3, 1:5, 2:4",
              "2:2, 2:4, 3:3"},
          2),

      // 5: A square divided into four subsquares. In this case we want
      // to extract the four loops rather than taking their union.
      // There are four extra edges as well.
      new TestCase(0, -1, 0.0, 4.0, // Don't XOR
          new Chain[] {new Chain("0:0, 0:5, 5:5, 5:0", true),
              new Chain("0:5, 0:10, 5:10, 5:5", true),
              new Chain("5:0, 5:5, 10:5, 10:0", true),
              new Chain("5:5, 5:10, 10:10, 10:5", true),
              new Chain("0:10, 0:15, 0:20", false),
              new Chain("20:0, 15:0, 10:0", false)},
          new String[] {"0:0, 0:5, 5:5, 5:0", "0:5, 0:10, 5:10, 5:5", "5:0, 5:5, 10:5, 10:0",
              "5:5, 5:10, 10:10, 10:5"},
          4),

      // 6: Five nested loops that touch at a point.
      new TestCase(0,
          0,
          0.0,
          0.8,
          new Chain[] {new Chain("0:0, 0:10, 10:10, 10:0", true),
              new Chain("0:0, 1:9, 9:9, 9:1", true), new Chain("0:0, 2:8, 8:8, 8:2", true),
              new Chain("0:0, 3:7, 7:7, 7:3", true), new Chain("0:0, 4:6, 6:6, 6:4", true)},
          new String[] {"0:0, 0:10, 10:10, 10:0", "0:0, 1:9, 9:9, 9:1", "0:0, 2:8, 8:8, 8:2",
              "0:0, 3:7, 7:7, 7:3", "0:0, 4:6, 6:6, 6:4"},
          0),


      // 7: Four diamonds nested within each other touching at two points.
      new TestCase(-1, 0, 0.0, 4.0, // Directed edges required for unique result.
          new Chain[] {new Chain("0:-20, -10:0, 0:20, 10:0", true),
              new Chain("0:10, -10:0, 0:-10, 10:0", true),
              new Chain("0:-10, -5:0, 0:10, 5:0", true), new Chain("0:5, -5:0, 0:-5, 5:0", true)},
          new String[] {"0:-20, -10:0, 0:-10, 10:0", "0:-10, -5:0, 0:-5, 5:0",
              "0:5, -5:0, 0:10, 5:0", "0:10, -10:0, 0:20, 10:0"},
          0),

      // 8: Seven diamonds nested within each other touching at one
      // point between each nested pair.
      new TestCase(0,
          0,
          0.0,
          9.0,
          new Chain[] {new Chain("0:-70, -70:0, 0:70, 70:0", true),
              new Chain("0:-70, -60:0, 0:60, 60:0", true),
              new Chain("0:-50, -60:0, 0:50, 50:0", true),
              new Chain("0:-40, -40:0, 0:50, 40:0", true),
              new Chain("0:-30, -30:0, 0:30, 40:0", true),
              new Chain("0:-20, -20:0, 0:30, 20:0", true),
              new Chain("0:-10, -20:0, 0:10, 10:0", true)},
          new String[] {"0:-70, -70:0, 0:70, 70:0",
              "0:-70, -60:0, 0:60, 60:0",
              "0:-50, -60:0, 0:50, 50:0",
              "0:-40, -40:0, 0:50, 40:0",
              "0:-30, -30:0, 0:30, 40:0",
              "0:-20, -20:0, 0:30, 20:0",
              "0:-10, -20:0, 0:10, 10:0"},
          0),

      // 9: A triangle and a self-intersecting bowtie.
      new TestCase(0,
          0,
          0.0,
          4.0,
          new Chain[] {new Chain("0:0, 0:10, 5:5", true), new Chain("0:20, 0:30, 10:20", false),
              new Chain("10:20, 10:30, 0:20", false)},
          new String[] {"0:0, 0:10, 5:5"},
          4),

      // 10: Two triangles that intersect each other.
      new TestCase(0,
          0,
          0.0,
          2.0,
          new Chain[] {new Chain("0:0, 0:10, 5:5", true), new Chain("2:2, 2:12, 7:7", true)},
          new String[] {},
          6),

      // 11: Four squares that combine to make a big square. The nominal
      // edges of the square are at +/-8.5 degrees in latitude and longitude.
      // All vertices except the center vertex are perturbed by up to 0.5
      // degrees in latitude and/or longitude. The various copies of the
      // center vertex are misaligned by more than this (i.e. they are
      // structured as a tree where adjacent vertices are separated by at
      // most 1 degree in latitude and/or longitude) so that the clustering
      // algorithm needs more than one iteration to find them all. Note that
      // the merged position of this vertex doesn't matter because it is XORed
      // away in the output.
      new TestCase(0, 1, 1.5, 5.8, // XOR, min_merge > sqrt(2), max_merge < 6.
          new Chain[] {new Chain("-8:-8, -8:0", false),
              new Chain("-8:1, -8:8", false),
              new Chain("0:-9, -2:0", false),
              new Chain("-1:1, 1:9", false),
              new Chain("0:8, 2:2", false),
              new Chain("0:-2, 1:-8", false),
              new Chain("8:9, 9:1", false),
              new Chain("9:0, 8:-9", false),
              new Chain("9:-9, 0:-8", false),
              new Chain("1:-9, -9:-9", false),
              new Chain("8:0, 1:0", false),
              new Chain("1:2, -8:0", false),
              new Chain("-8:1, 1:-1", false),
              new Chain("0:1, 8:1", false),
              new Chain("-9:8, 1:8", false),
              new Chain("0:9, 8:8", false)},
          new String[] {"8.5:8.5, 8.5:0.5, 8.5:-8.5, 0.5:-8.5, "
              + "-8.5:-8.5, -8.5:0.5, -8.5:8.5, 0.5:8.5"},
          0)};

  private void getVertices(String str,
      S2Point x,
      S2Point y,
      S2Point z,
      double maxPerturbation,
      List<S2Point> vertices) {

    // Parse the vertices, perturb them if desired, and transform them into the
    // given frame.
    S2Polyline line = makePolyline(str);

    for (int i = 0; i < line.numVertices(); ++i) {
      S2Point p = line.vertex(i);
      // (p[0]*x + p[1]*y + p[2]*z).Normalize()
      S2Point axis = S2Point.normalize(
          S2Point.add(S2Point.add(S2Point.mul(x, p.x), S2Point.mul(y, p.y)), S2Point.mul(z, p.z)));
      S2Cap cap = S2Cap.fromAxisAngle(axis, S1Angle.radians(maxPerturbation));
      vertices.add(samplePoint(cap));
    }
  }

  private boolean loopsEqual(S2Loop a, S2Loop b, double maxError) {
    // Return true if two loops have the same cyclic vertex sequence.

    if (a.numVertices() != b.numVertices()) {
      return false;
    }
    for (int offset = 0; offset < a.numVertices(); ++offset) {
      if (S2.approxEquals(a.vertex(offset), b.vertex(0), maxError)) {
        boolean success = true;
        for (int i = 0; i < a.numVertices(); ++i) {
          if (!S2.approxEquals(a.vertex(i + offset), b.vertex(i), maxError)) {
            success = false;
            break;
          }
        }
        if (success) {
          return true;
        }
        // Otherwise continue looping. There may be more than one candidate
        // starting offset since vertices are only matched approximately.
      }
    }
    return false;
  }

  private boolean findLoop(S2Loop loop, List<S2Loop> candidates, double maxError) {
    for (int i = 0; i < candidates.size(); ++i) {
      if (loopsEqual(loop, candidates.get(i), maxError)) {
        return true;
      }
    }
    return false;
  }

  boolean findMissingLoops(
      List<S2Loop> actual, List<S2Loop> expected, double maxError, String label) {
    // Dump any loops from "actual" that are not present in "expected".
    boolean found = false;
    for (int i = 0; i < actual.size(); ++i) {
      if (findLoop(actual.get(i), expected, maxError)) {
        continue;
      }
      System.err.print(label + " loop " + i + ":\n");
      S2Loop loop = actual.get(i);
      for (int j = 0; j < loop.numVertices(); ++j) {
        S2Point p = loop.vertex(j);
        System.err.print("   [" + p.x + ", " + p.y + ", " + p.z + "]\n");
      }
      found = true;
    }
    return found;
  }

  void addChain(Chain chain,
      S2Point x,
      S2Point y,
      S2Point z,
      double maxPerturbation,
      S2PolygonBuilder builder) {

    // Transform the given edge chain to the frame (x,y,z), perturb each vertex
    // up to the given distance, and add it to the builder.

    List<S2Point> vertices = Lists.newArrayList();
    getVertices(chain.str, x, y, z, maxPerturbation, vertices);
    if (chain.closed) {
      vertices.add(vertices.get(0));
    }
    for (int i = 1; i < vertices.size(); ++i) {
      builder.addEdge(vertices.get(i - 1), vertices.get(i));
    }
  }

  boolean evalTristate(int state) {
    return (state > 0) ? true : (state < 0) ? false : (rand.nextDouble() > 0.5);
  }

  boolean testBuilder(TestCase test) {
    for (int iter = 0; iter < 200; ++iter) {
      // Initialize to the default options, which are changed below
      S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR;

      options.setUndirectedEdges(evalTristate(test.undirectedEdges));
      options.setXorEdges(evalTristate(test.xorEdges));

      // Each test has a minimum and a maximum merge distance. The merge
      // distance must be at least the given minimum to ensure that all expected
      // merging will take place, and it must be at most the given maximum to
      // ensure that no unexpected merging takes place.
      //
      // If the minimum and maximum values are different, we have some latitude
      // to perturb the vertices as long as the merge distance is adjusted
      // appropriately. If "p" is the maximum perturbation distance, "min" and
      // "max" are the min/max merge distances, and "m" is the actual merge
      // distance for this test, we require that
      //
      // x >= min + 2*p and x <= max - 2*p .
      //
      // This implies that p <= 0.25 * (max - min). We choose "p" so that it is
      // zero half of the time, and otherwise chosen randomly up to this limit.

      double minMerge = S1Angle.degrees(test.minMerge).radians();
      double maxMerge = S1Angle.degrees(test.maxMerge).radians();
      double r = Math.max(0.0, 2 * rand.nextDouble() - 1);
      double maxPerturbation = r * 0.25 * (maxMerge - minMerge);

      // Now we set the merge distance chosen randomly within the limits above
      // (min + 2*p and max - 2*p). Half of the time we set the merge distance
      // to the minimum value.

      r = Math.max(0.0, 2 * rand.nextDouble() - 1);
      options.setMergeDistance(S1Angle.radians(
          minMerge + 2 * maxPerturbation + r * (maxMerge - minMerge - 4 * maxPerturbation)));

      options.setValidate(true);
      S2PolygonBuilder builder = new S2PolygonBuilder(options);

      // On each iteration we randomly rotate the test case around the sphere.
      // This causes the S2PolygonBuilder to choose different first edges when
      // trying to build loops.
      S2Point x = randomPoint();
      S2Point y = S2Point.normalize(S2Point.crossProd(x, randomPoint()));
      S2Point z = S2Point.normalize(S2Point.crossProd(x, y));

      for (Chain chain : test.chainsIn) {
        addChain(chain, x, y, z, maxPerturbation, builder);
      }
      List<S2Loop> loops = Lists.newArrayList();
      List<S2Edge> unusedEdges = Lists.newArrayList();
      if (test.xorEdges < 0) {
        builder.assembleLoops(loops, unusedEdges);
      } else {
        S2Polygon polygon = new S2Polygon();
        builder.assemblePolygon(polygon, unusedEdges);
        polygon.release(loops);
      }
      List<S2Loop> expected = Lists.newArrayList();
      for (String loop : test.loopsOut) {
        List<S2Point> vertices = Lists.newArrayList();
        getVertices(loop, x, y, z, 0, vertices);
        expected.add(new S2Loop(vertices));
      }
      // We assume that the vertex locations in the expected output polygon
      // are separated from the corresponding vertex locations in the input
      // edges by at most half of the minimum merge distance. Essentially
      // this means that the expected output vertices should be near the
      // centroid of the various input vertices.
      double maxError = 0.5 * minMerge + maxPerturbation;

      // Note single "|" below so that we print both sets of loops.
      if (findMissingLoops(loops, expected, maxError, "Actual")
          | findMissingLoops(expected, loops, maxError, "Expected")) {
        System.err.print(
            "During iteration " + iter + ", undirected: " + options.getUndirectedEdges() + ", xor: "
                + options.getXorEdges() + "\n\n");
        return false;
      }
      if (unusedEdges.size() != test.numUnusedEdges) {
        System.err.print("Wrong number of unused edges: " + unusedEdges.size() + "%d (should be "
            + test.numUnusedEdges + ")\n");
        return false;
      }
    }
    return true;
  }

  public void testAssembleLoops() {
    boolean success = true;
    for (int i = 0; i < testCases.length; ++i) {
      log.info("Starting test case " + i);

      boolean caseSuccess = testBuilder(testCases[i]);

      log.info("Test case " + i + " finished: " + ((caseSuccess) ? "SUCCESS" : "FAILED"));

      success &= caseSuccess;
    }
    assertTrue(success);
  }
}
