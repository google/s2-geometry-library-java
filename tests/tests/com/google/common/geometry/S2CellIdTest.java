/*
 * Copyright 2005 Google Inc.
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

import static com.google.common.geometry.S2CellId.MAX_LEVEL;
import static com.google.common.geometry.S2Projections.PROJ;
import static java.lang.Math.min;
import static java.lang.Math.round;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Unit tests for S2CellId.
 * @author ericv@google.com (Eric Veach)
 */
public strictfp class S2CellIdTest extends GeometryTestCase {
  private static S2CellId getCellId(double latDegrees, double lngDegrees) {
    return S2CellId.fromLatLng(S2LatLng.fromDegrees(latDegrees, lngDegrees));
  }

  public void testDefaultConstructor() {
    S2CellId id = new S2CellId();
    assertEquals(0, id.id());
    assertFalse(id.isValid());
  }

  public void testFaceDefinitions() {
    assertEquals(0, getCellId(0, 0).face());
    assertEquals(1, getCellId(0, 90).face());
    assertEquals(2, getCellId(90, 0).face());
    assertEquals(3, getCellId(0, 180).face());
    assertEquals(4, getCellId(0, -90).face());
    assertEquals(5, getCellId(-90, 0).face());
  }

  public void testFromFace() {
    for (int face = 0; face < 6; ++face) {
      assertEquals(S2CellId.fromFace(face), S2CellId.fromFacePosLevel(face, 0, 0));
    }
  }

  public void testParentChildRelationships() {
    S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, MAX_LEVEL - 4);
    assertTrue(id.isValid());
    assertEquals(3, id.face());
    assertEquals(0x12345700, id.pos());
    assertEquals(MAX_LEVEL - 4, id.level());
    assertFalse(id.isLeaf());

    assertEquals(0x12345610, id.childBegin(id.level() + 2).pos());
    assertEquals(0x12345640, id.childBegin().pos());
    assertEquals(0x12345400, id.parent().pos());
    assertEquals(0x12345000, id.parent(id.level() - 2).pos());

    // Check ordering of children relative to parents.
    assertTrue(id.childBegin().lessThan(id));
    assertTrue(id.childEnd().greaterThan(id));
    assertEquals(id.childEnd(), id.childBegin().next().next().next().next());
    assertEquals(id.rangeMin(), id.childBegin(MAX_LEVEL));
    assertEquals(id.rangeMax().next(), id.childEnd(MAX_LEVEL));

    // Check that cells are represented by the position of their center along the Hilbert curve.
    assertEquals(2 * id.id(), id.rangeMin().id() + id.rangeMax().id());
  }

  public void testCenterSiTi() {
    S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, MAX_LEVEL);
    // Check that the (si, ti) coordinates of the center end in a 1 followed by (30 - level) 0s.

    // Leaf level, 30.
    long center = id.getCenterSiTi();
    assertEquals(1 << 0, S2CellId.getSi(center) & 1);
    assertEquals(1 << 0, S2CellId.getTi(center) & 1);

    // Level 29.
    center = id.parent(MAX_LEVEL - 1).getCenterSiTi();
    assertEquals(1 << 1, S2CellId.getSi(center) & 3);
    assertEquals(1 << 1, S2CellId.getTi(center) & 3);

    // Level 28.
    center = id.parent(MAX_LEVEL - 2).getCenterSiTi();
    assertEquals(1 << 2, S2CellId.getSi(center) & 7);
    assertEquals(1 << 2, S2CellId.getTi(center) & 7);

    // Level 20.
    center = id.parent(MAX_LEVEL - 10).getCenterSiTi();
    assertEquals(1 << 10, S2CellId.getSi(center) & ((1 << 11) - 1));
    assertEquals(1 << 10, S2CellId.getTi(center) & ((1 << 11) - 1));

    // Level 10.
    center = id.parent(MAX_LEVEL - 20).getCenterSiTi();
    assertEquals(1 << 20, S2CellId.getSi(center) & ((1 << 21) - 1));
    assertEquals(1 << 20, S2CellId.getTi(center) & ((1 << 21) - 1));

    // Level 0.
    center = id.parent(0).getCenterSiTi();
    assertEquals(1L << 30, S2CellId.getSi(center) & ((1L << 31) - 1));
    assertEquals(1L << 30, S2CellId.getTi(center) & ((1L << 31) - 1));
  }

  /** Check wrapping from beginning of Hilbert curve to end and vice versa. */
  public void testWrapping() {
    assertEquals(S2CellId.end(0).prev(), S2CellId.begin(0).prevWrap());
    assertEquals(
        S2CellId.fromFacePosLevel(5, ~0L >>> S2CellId.FACE_BITS, MAX_LEVEL),
        S2CellId.begin(MAX_LEVEL).prevWrap());
    assertEquals(
        S2CellId.fromFacePosLevel(5, ~0L >>> S2CellId.FACE_BITS, MAX_LEVEL),
        S2CellId.begin(MAX_LEVEL).advanceWrap(-1));
    assertEquals(S2CellId.begin(4), S2CellId.end(4).prev().nextWrap());
    assertEquals(S2CellId.begin(4), S2CellId.end(4).advance(-1).advanceWrap(1));
    assertEquals(
        S2CellId.fromFacePosLevel(0, 0, MAX_LEVEL), S2CellId.end(MAX_LEVEL).prev().nextWrap());
    assertEquals(
        S2CellId.fromFacePosLevel(0, 0, MAX_LEVEL),
        S2CellId.end(MAX_LEVEL).advance(-1).advanceWrap(1));
    assertEquals(1, S2CellId.fromFacePosLevel(1, 0x9e5bd23c2e4b694L, 29).prevWrap().face());
  }

  public void testAdvance() {
    S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, MAX_LEVEL - 4);
    // Check basic properties of advance().
    assertEquals(S2CellId.end(0), S2CellId.begin(0).advance(7));
    assertEquals(S2CellId.end(0), S2CellId.begin(0).advance(12));
    assertEquals(S2CellId.begin(0), S2CellId.end(0).advance(-7));
    assertEquals(S2CellId.begin(0), S2CellId.end(0).advance(-12000000));
    int numLevel5Cells = 6 << (2 * 5);
    assertEquals(S2CellId.end(5).advance(500 - numLevel5Cells), S2CellId.begin(5).advance(500));
    assertEquals(id.next().childBegin(MAX_LEVEL), id.childBegin(MAX_LEVEL).advance(256));
    assertEquals(
        S2CellId.fromFacePosLevel(5, 0, MAX_LEVEL),
        S2CellId.fromFacePosLevel(1, 0, MAX_LEVEL).advance(4L << (2 * MAX_LEVEL)));

    // Check basic properties of advanceWrap().
    assertEquals(S2CellId.fromFace(1), S2CellId.begin(0).advanceWrap(7));
    assertEquals(S2CellId.begin(0), S2CellId.begin(0).advanceWrap(12));
    assertEquals(S2CellId.fromFace(4), S2CellId.fromFace(5).advanceWrap(-7));
    assertEquals(S2CellId.begin(0), S2CellId.begin(0).advanceWrap(-12000000));
    assertEquals(S2CellId.begin(5).advanceWrap(6644), S2CellId.begin(5).advanceWrap(-11788));
    assertEquals(id.next().childBegin(MAX_LEVEL), id.childBegin(MAX_LEVEL).advanceWrap(256));
    assertEquals(
        S2CellId.fromFacePosLevel(1, 0, MAX_LEVEL),
        S2CellId.fromFacePosLevel(5, 0, MAX_LEVEL).advanceWrap(2L << (2 * MAX_LEVEL)));
  }

  public void testDistanceFromBegin() {
    assertEquals(6, S2CellId.end(0).distanceFromBegin());
    assertEquals(6 * (1L << (2 * MAX_LEVEL)), S2CellId.end(MAX_LEVEL).distanceFromBegin());

    assertEquals(0, S2CellId.begin(0).distanceFromBegin());
    assertEquals(0, S2CellId.begin(MAX_LEVEL).distanceFromBegin());

    S2CellId id = S2CellId.fromFacePosLevel(3, 0x12345678, MAX_LEVEL - 4);
    assertEquals(id, S2CellId.begin(id.level()).advance(id.distanceFromBegin()));
  }

  public void testInverses() {
    // Check the conversion of random leaf cells to S2LatLngs and back.
    for (int i = 0; i < 200000; ++i) {
      S2CellId id = data.getRandomCellId(MAX_LEVEL);
      assertTrue(id.isLeaf());
      assertEquals(MAX_LEVEL, id.level());
      S2LatLng center = id.toLatLng();
      assertEquals(id.id(), S2CellId.fromLatLng(center).id());
    }
  }

  public void testGetCommonAncestorLevel() {
    // Two identical cell ids.
    S2CellId face0 = S2CellId.fromFace(0);
    S2CellId leaf0 = face0.childBegin(30);
    assertEquals(0, face0.getCommonAncestorLevel(face0));
    assertEquals(30, leaf0.getCommonAncestorLevel(leaf0));

    // One cell id is a descendant of the other.
    S2CellId face5 = S2CellId.fromFace(5);
    assertEquals(0, leaf0.getCommonAncestorLevel(face0));
    assertEquals(0, face5.getCommonAncestorLevel(face5.childEnd(30).prev()));

    // Two cells that have no common ancestor.
    assertEquals(-1, face0.getCommonAncestorLevel(face5));
    assertEquals(
        -1,
        S2CellId.fromFace(2)
            .childBegin(30)
            .getCommonAncestorLevel(S2CellId.fromFace(3).childEnd(20)));

    // Two cells that have a common ancestor distinct from both of them.
    S2CellId face5child9 = face5.childBegin(9);
    S2CellId face0child2 = face0.childBegin(2);
    assertEquals(
        8, face5child9.next().childBegin(15).getCommonAncestorLevel(face5child9.childBegin(20)));
    assertEquals(
        1, face0child2.childBegin(30).getCommonAncestorLevel(face0child2.next().childBegin(5)));
  }

  public void testTokens() {
    // Test random cell ids at all levels.
    for (int i = 0; i < 10000; ++i) {
      S2CellId id = data.getRandomCellId();
      String token = id.toToken();
      assertTrue(token.length() <= 16);
      assertEquals(id, S2CellId.fromToken(token));
    }

    // Check that invalid cell ids can be encoded.
    String token = S2CellId.none().toToken();
    assertEquals(S2CellId.none(), S2CellId.fromToken(token));
  }

  private static final int MAX_EXPAND_LEVEL = 3;

  private static void expandCell(
      S2CellId parent, List<S2CellId> cells, Map<S2CellId, S2CellId> parentMap) {
    cells.add(parent);
    if (parent.level() == MAX_EXPAND_LEVEL) {
      return;
    }

    for (int pos = 0; pos < 4; pos++) {
      S2CellId child = parent.child(pos);
      // Do some basic checks on the children
      assertEquals(child, parent.child(pos));
      assertEquals(parent.level() + 1, child.level());
      assertFalse(child.isLeaf());
      assertEquals(parent.face(), child.face());
      assertEquals(parent.getOrientation() ^ S2.posToOrientation(pos), child.getOrientation());
      parentMap.put(child, parent);
      expandCell(child, cells, parentMap);
    }
  }

  /** Test contains() and intersects(). */
  public void testContainment() {
    Map<S2CellId, S2CellId> parentMap = Maps.newHashMap();
    List<S2CellId> cells = Lists.newArrayList();
    for (int face = 0; face < 6; ++face) {
      expandCell(S2CellId.fromFace(face), cells, parentMap);
    }
    for (int i = 0; i < cells.size(); ++i) {
      S2CellId ci = cells.get(i);
      for (int j = 0; j < cells.size(); ++j) {
        S2CellId cj = cells.get(j);
        boolean contained = true;
        for (S2CellId id = cj; !id.equals(ci); id = parentMap.get(id)) {
          if (!parentMap.containsKey(id)) {
            contained = false;
            break;
          }
        }
        assertEquals(contained, ci.contains(cj));
        assertEquals(
            contained, cj.greaterOrEquals(ci.rangeMin()) && cj.lessOrEquals(ci.rangeMax()));
        assertEquals(ci.contains(cj) || cj.contains(ci), ci.intersects(cj));
      }
    }
  }

  private static final int MAX_WALK_LEVEL = 8;

  /**
   * Verifies that sequentially increasing cell ids form a continuous path over the surface of the
   * sphere, i.e. there are no discontinuous jumps from one region to another.
   */
  public void testContinuity() {
    double maxDist = PROJ.maxEdge.getValue(MAX_WALK_LEVEL);
    S2CellId end = S2CellId.end(MAX_WALK_LEVEL);
    for (S2CellId id = S2CellId.begin(MAX_WALK_LEVEL); !id.equals(end); id = id.next()) {
      assertTrue(id.toPointRaw().angle(id.nextWrap().toPointRaw()) <= maxDist);
      assertEquals(id.nextWrap(), id.advanceWrap(1));
      assertEquals(id, id.nextWrap().advanceWrap(-1));

      // Check that the toPointRaw() returns the center of each cell in (s,t) coordinates.
      S2Point p = id.toPointRaw();
      int face = S2Projections.xyzToFace(p);
      R2Vector uv = S2Projections.validFaceXyzToUv(face, p);
      final double cellSize = 1.0 / (1 << MAX_WALK_LEVEL);
      assertEquals(0.0, drem(PROJ.uvToST(uv.x()), 0.5 * cellSize), 1e-15);
      assertEquals(0.0, drem(PROJ.uvToST(uv.y()), 0.5 * cellSize), 1e-15);
    }
  }

  private static double drem(double num, double dem) {
    return num - round(num / dem) * dem;
  }

  /**
   * Verifies that random points on the sphere can be represented to the expected level of accuracy,
   * which in the worst case is sqrt(2/3) times the maximum arc length between the points on the
   * sphere associated with adjacent values of "i" or "j". (It is sqrt(2/3) rather than 1/2 because
   * the cells at the corners of each face are stretched -- they have 60 and 120 degree angles.)
   */
  public void testCoverage() {
    double maxDist = 0.5 * PROJ.maxDiag.getValue(MAX_LEVEL);
    for (int i = 0; i < 1000000; ++i) {
      S2Point p = data.getRandomPoint();
      S2Point q = S2CellId.fromPoint(p).toPointRaw();
      assertTrue(p.angle(q) <= maxDist);
    }
  }

  private static void checkAllNeighbors(S2CellId id, int level) {
    assertTrue(level >= id.level());
    assertTrue(level < MAX_LEVEL);

    // We compute appendAllNeighbors, and then add in all the children of "id" at the given level.
    // We then compare this against the result of finding all the vertex neighbors of all the
    // vertices of children of "id" at the given level. These should give the same result.
    List<S2CellId> all = Lists.newArrayList();
    List<S2CellId> expected = Lists.newArrayList();
    id.getAllNeighbors(level, all);
    S2CellId end = id.childEnd(level + 1);
    for (S2CellId c = id.childBegin(level + 1); !c.equals(end); c = c.next()) {
      all.add(c.parent());
      c.getVertexNeighbors(level, expected);
    }

    // Verify the results contain the same sets of unique cells.
    assertEquals(Sets.newTreeSet(all), Sets.newTreeSet(expected));
  }

  /** Checks the edge neighbors of face 1. */
  public void testNeighbors() {
    final int[] outFaces = {5, 3, 2, 0};
    S2CellId[] faceNbrs = new S2CellId[4];
    S2CellId.fromFace(1).getEdgeNeighbors(faceNbrs);
    for (int i = 0; i < 4; ++i) {
      assertTrue(faceNbrs[i].isFace());
      assertEquals(outFaces[i], faceNbrs[i].face());
    }

    // Check the edge neighbors of the corner cells at all levels. This case is trickier because it
    // requires projecting onto adjacent faces.
    final int maxIJ = S2CellId.MAX_SIZE - 1;
    for (int level = 1; level <= MAX_LEVEL; ++level) {
      S2CellId id = S2CellId.fromFaceIJ(1, 0, 0).parent(level);
      S2CellId[] nbrs = new S2CellId[4];
      id.getEdgeNeighbors(nbrs);
      // These neighbors were determined manually using the face and axis relationships defined in
      // s2.cc.
      int sizeIJ = S2CellId.getSizeIJ(level);
      assertEquals(nbrs[0], S2CellId.fromFaceIJ(5, maxIJ, maxIJ).parent(level));
      assertEquals(nbrs[1], S2CellId.fromFaceIJ(1, sizeIJ, 0).parent(level));
      assertEquals(nbrs[2], S2CellId.fromFaceIJ(1, 0, sizeIJ).parent(level));
      assertEquals(nbrs[3], S2CellId.fromFaceIJ(0, maxIJ, 0).parent(level));
    }

    // Check the vertex neighbors of the center of face 2 at level 5.
    List<S2CellId> nbrs = Lists.newArrayList();
    S2CellId.fromPoint(S2Point.Z_POS).getVertexNeighbors(5, nbrs);
    Collections.sort(nbrs);
    for (int child = 0; child < 4; child++) {
      int i = (1 << 29) - (child < 2 ? 1 : 0);
      int j = (1 << 29) - ((child == 0 || child == 3) ? 1 : 0);
      assertEquals(S2CellId.fromFaceIJ(2, i, j).parent(5), nbrs.get(child));
    }
    nbrs.clear();

    // Check the vertex neighbors of the corner of faces 0, 4, and 5.
    S2CellId id = S2CellId.fromFacePosLevel(0, 0, MAX_LEVEL);
    id.getVertexNeighbors(0, nbrs);
    Collections.sort(nbrs);
    assertEquals(
        Lists.newArrayList(S2CellId.fromFace(0), S2CellId.fromFace(4), S2CellId.fromFace(5)), nbrs);

    // Check that appendAllNeighbors produces results that are consistent with getVertexNeighbors
    // for a bunch of random cells.
    for (int i = 0; i < TestPlatform.S2_CELL_ID_TEST_RANDOM_NEIGHBORS_CHECKS; ++i) {
      do {
        id = data.getRandomCellId();
      } while (id.level() <= 2 || id.level() >= 29);

      // testAllNeighbors computes approximately 2**(2*(diff+1)) cell ids, so it's not reasonable to
      // use large values of "diff".
      int maxDiff = min(6, MAX_LEVEL - id.level() - 1);
      int level = id.level() + data.nextInt(maxDiff);
      checkAllNeighbors(id, level);
    }
  }

  public void testToLoop() {
    try {
      S2CellId.FACE_CELLS[0].child(0).toLoop(0);
      fail();
    } catch (Exception e) {
      // pass.
    }

    for (S2CellId id = S2CellId.FACE_CELLS[0]; id.level() < 30 - 3; id = id.child(id.level() % 4)) {
      for (int level = id.level(); level < id.level() + 3; level++) {
        S2Loop actual = id.toLoop(level);
        S2PolygonBuilder builder = new S2PolygonBuilder();
        for (S2CellId child : id.childrenAtLevel(level)) {
          builder.addLoop(new S2Loop(new S2Cell(child)));
        }
        S2Polygon expected = builder.assemblePolygon();
        assertEquals(1, expected.numLoops());
        assertEquals(actual.numVertices(), expected.loop(0).numVertices());
        int offset = -1;
        for (int i = 0; i < expected.loop(0).numVertices(); i++) {
          if (actual.vertex(0).equals(expected.loop(0).vertex(i))) {
            offset = i;
            break;
          }
        }
        assertTrue(offset >= 0);
        for (int i = 1; i < actual.numVertices(); i++) {
          assertEquals(expected.loop(0).vertex(i + offset), actual.vertex(i));
        }
      }
    }
  }

  /**
   * Tests that invalid face arguments for {@link S2Projections#validFaceXyzToUv(int, S2Point)}
   * result in index out of bounds exceptions.
   */
  public void testInvalidFaceValuesThrow() {
    S2Point arbitrary = S2Point.X_POS;

    // Test arbitrary negative face value.
    try {
      S2Projections.validFaceXyzToUv(-1, arbitrary);
      fail("Expected index out of bounds exception.");
    } catch (IndexOutOfBoundsException | NullPointerException expected) {
      // Pass because the exception was expected.
      //
      // NPE is okay because translation to Javascript results in using Javascript arrays, which
      // have no bounds checks.
    }

    // Test invalid positive face value.
    try {
      S2Projections.validFaceXyzToUv(6, arbitrary);
      fail("Expected index out of bounds exception.");
    } catch (IndexOutOfBoundsException | NullPointerException expected) {
      // Pass because the exception was expected.
      //
      // NPE is okay because translation to Javascript results in using Javascript arrays, which
      // have no bounds checks.
    }
  }

  public void testCoder()  {
    List<S2Coder<S2CellId>> coders = ImmutableList.of(S2CellId.CODER, S2CellId.TOKEN_CODER);
    for (S2Coder<S2CellId> coder : coders) {
      assertEquals(S2CellId.FACE_CELLS[5], roundtrip(coder, S2CellId.FACE_CELLS[5]));
    }
  }
}
