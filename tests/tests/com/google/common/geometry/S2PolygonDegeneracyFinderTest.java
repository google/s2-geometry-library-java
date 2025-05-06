/*
 * Copyright 2023 Google Inc.
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

import static com.google.common.geometry.S2PolygonDegeneracyFinder.findPolygonDegeneracies;
import static com.google.common.geometry.S2PolygonDegeneracyFinder.isFullyDegenerate;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2Builder.IsFullPolygonPredicate;
import com.google.common.geometry.S2PolygonDegeneracyFinder.PolygonDegeneracyList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2PolygonDegeneracyFinder}. */
@RunWith(JUnit4.class)
public final class S2PolygonDegeneracyFinderTest extends GeometryTestCase {

  @Test
  public void testEmptyPolygon() {
    expectDegeneracies("", ds());
  }

  @Test
  public void testNoDegeneracies() {
    expectDegeneracies("0:0, 0:1, 1:0", ds());
  }

  @Test
  public void testPointShell() {
    expectDegeneracies("0:0", ds(d("0:0, 0:0", false)));
  }

  @Test
  public void testSiblingPairShells() {
    expectDegeneracies(
        "0:0, 0:1, 1:0; 1:0, 0:1, 0:0",
        ds(
            d("0:0, 0:1", false),
            d("0:1, 0:0", false),
            d("0:1, 1:0", false),
            d("1:0, 0:1", false),
            d("0:0, 1:0", false),
            d("1:0, 0:0", false)));
  }

  @Test
  public void testAttachedSiblingPairShells() {
    expectDegeneracies( //
        "0:0, 0:1, 1:0; 1:0, 2:0", //
        ds(d("1:0, 2:0", false), d("2:0, 1:0", false)));
  }

  @Test
  public void testAttachedSiblingPairHoles() {
    expectDegeneracies( //
        "0:0, 0:3, 3:0; 0:0, 1:1", //
        ds(d("0:0, 1:1", true), d("1:1, 0:0", true)));
  }

  @Test
  public void testAttachedSiblingPairShellsAndHoles() {
    expectDegeneracies(
        "0:0, 0:3, 3:0; 3:0, 1:1; 3:0, 5:5",
        ds(
            d("3:0, 1:1", true), d("1:1, 3:0", true),
            d("3:0, 5:5", false), d("5:5, 3:0", false)));
  }

  @Test
  public void testDegenerateShellsOutsideLoop() {
    expectDegeneracies(
        "0:0, 0:3, 3:3, 3:0; 4:4, 5:5; 6:6",
        ds(d("4:4, 5:5", false), d("5:5, 4:4", false), d("6:6, 6:6", false)));
  }

  @Test
  public void testDegenerateHolesWithinLoop() {
    expectDegeneracies(
        "0:0, 0:5, 5:5, 5:0; 1:1, 2:2; 3:3",
        ds(d("1:1, 2:2", true), d("2:2, 1:1", true), d("3:3, 3:3", true)));
  }

  @Test
  public void testPointHoleWithinFull() {
    expectDegeneracies("full; 0:0", ds(d("0:0, 0:0", true)));
  }

  @Test
  public void testSiblingPairHolesWithinFull() {
    expectDegeneracies(
        "full; 0:0, 0:1, 1:0; 1:0, 0:1, 0:0",
        ds(
            d("0:0, 0:1", true),
            d("0:1, 0:0", true),
            d("0:1, 1:0", true),
            d("1:0, 0:1", true),
            d("0:0, 1:0", true),
            d("1:0, 0:0", true)));
  }

  /** A container holdering a description of a degeneracies' edges, and if it is a hole. */
  private static class TestDegeneracy implements Comparable<TestDegeneracy> {
    public final String edgeStr;
    public final boolean isHole;

    /** Constructs taking a string describing the edge, and if it is a hole. */
    public TestDegeneracy(String edgeStr, boolean isHole) {
      this.edgeStr = edgeStr;
      this.isHole = isHole;
    }

    @Override
    public int compareTo(TestDegeneracy that) {
      if (this.edgeStr.equals(that.edgeStr)) {
        return Boolean.compare(this.isHole, that.isHole);
      }
      return this.edgeStr.compareTo(that.edgeStr);
    }

    @Override
    public boolean equals(Object o) {
      if (o instanceof TestDegeneracy) {
        TestDegeneracy that = (TestDegeneracy) o;
        return this.edgeStr.equals(that.edgeStr) && this.isHole == that.isHole;
      }
      return false;
    }

    @Override
    public int hashCode() {
      return edgeStr.hashCode() + (isHole ? 1 : 0);
    }

    @Override
    public String toString() {
      return edgeStr + " " + (isHole ? "Hole" : "Shell");
    }
  }

  private static class DegeneracyCheckingLayer implements S2BuilderLayer {
    private final ArrayList<TestDegeneracy> expected;

    public DegeneracyCheckingLayer(List<TestDegeneracy> expected) {
      // Copy the expected degeneracies into an ArrayList so they may be sorted.
      this.expected = new ArrayList<>(expected);
    }

    @Override
    public GraphOptions graphOptions() {
      return new GraphOptions(
          EdgeType.DIRECTED,
          DegenerateEdges.DISCARD_EXCESS,
          DuplicateEdges.KEEP,
          SiblingPairs.DISCARD_EXCESS);
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      PolygonDegeneracyList degeneracies = findPolygonDegeneracies(g);
      // Convert the output into a human-readable format.
      ArrayList<TestDegeneracy> actual = new ArrayList<>();
      for (int i = 0; i < degeneracies.size(); i++) {
        int edgeId = degeneracies.edgeId(i);
        ImmutableList<S2Point> p = ImmutableList.of(
            g.vertex(g.edges().getSrcId(edgeId)), g.vertex(g.edges().getDstId(edgeId)));
        actual.add(new TestDegeneracy(S2TextFormat.s2PointsToString(p), degeneracies.isHole(i)));
      }

      Collections.sort(actual);
      Collections.sort(expected);

      StringBuilder sb = new StringBuilder();
      sb.append("Expected Degeneracies:\n");
      append(sb, expected);
      sb.append("\nActual Degeneracies:\n");
      append(sb, actual);
      sb.append("\n");

      assertEquals(sb.toString(), expected.size(), actual.size());

      for (int i = 0; i < actual.size(); i++) {
        TestDegeneracy expectedDegeneracy = expected.get(i);
        TestDegeneracy actualDegeneracy = actual.get(i);
        assertEquals(sb.toString(), expectedDegeneracy, actualDegeneracy);
      }

      assertEquals(sb.toString(), isFullyDegenerate(g), degeneracies.size() == g.numEdges());
      return true;
    }

    private static void append(StringBuilder sb, List<TestDegeneracy> degeneracies) {
      sb.append("[");
      for (TestDegeneracy degeneracy : degeneracies) {
        sb.append(degeneracy).append("; ");
      }
      sb.append("]");
    }
  }

  private void expectDegeneracies(String polygonStr, List<TestDegeneracy> expected) {
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    DegeneracyCheckingLayer layer = new DegeneracyCheckingLayer(expected);
    builder.startLayer(layer);
    S2LaxPolygonShape polygon = S2TextFormat.makeLaxPolygonOrDie(polygonStr);
    builder.addIsFullPolygonPredicate(
        new IsFullPolygonPredicate() {
          @Override
          public boolean test(S2BuilderGraph graph) {
            return polygon.getReferencePoint().contained();
          }
        });
    builder.addShape(polygon);
    S2Error error = new S2Error();
    boolean built = builder.build(error);
    assertTrue(error.toString(), built);
  }

  private ImmutableList<TestDegeneracy> ds(TestDegeneracy... degeneracies) {
    return ImmutableList.copyOf(degeneracies);
  }

  private TestDegeneracy d(String edgeStr, boolean isHole) {
    return new TestDegeneracy(edgeStr, isHole);
  }
}
