/*
 * Copyright 2022 Google Inc.
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

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2CrossingEdgesQuery.CrossingType;
import com.google.common.geometry.S2CrossingEdgesQuery.EdgePairVisitor;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;
import java.util.HashSet;
import java.util.List;

/** Tests S2CrossingEdgesQuery. */
public class S2CrossingEdgesQueryTest extends GeometryTestCase {

  /** Check that an empty index has no crossing edge pairs. */
  public void testGetCrossingEdgePairsEmptyIndex() {
    S2ShapeIndex index = new S2ShapeIndex();
    checkGetCrossingEdgePairs(index, CrossingType.ALL, 0);
    checkGetCrossingEdgePairs(index, CrossingType.INTERIOR, 0);
  }

  /** Check that two empty indexes have no crossing edge pairs. */
  public void testGetCrossingEdgePairsEmptyIndexTwoIndexes() {
    S2ShapeIndex indexA = new S2ShapeIndex();
    S2ShapeIndex indexB = new S2ShapeIndex();
    checkGetCrossingEdgePairs(indexA, indexB, CrossingType.ALL, 0);
    checkGetCrossingEdgePairs(indexA, indexB, CrossingType.INTERIOR, 0);
  }

  /** Test crossing edges in a grid of 11 vertical and 11 horizontal lines. */
  public void testGetCrossingEdgePairsEdgeGrid() {
    int kGridSize = 10; // grid is zero to kGridSize inclusive
    double epsilon = 1e-10;

    S2ShapeIndex index = new S2ShapeIndex();
    S2EdgeVectorShape shape = new S2EdgeVectorShape();
    for (int i = 0; i <= kGridSize; ++i) {
      double e = (i == 0 || i == kGridSize) ? 0 : epsilon;
      shape.add(
          S2LatLng.fromDegrees(-e, i).toPoint(), S2LatLng.fromDegrees(kGridSize + e, i).toPoint());
      shape.add(
          S2LatLng.fromDegrees(i, -e).toPoint(), S2LatLng.fromDegrees(i, kGridSize + e).toPoint());
    }
    index.add(shape);

    // Compare brute force and optimized approaches to getting all crossing edge pairs. There are
    // 11 horizontal and 11 vertical lines. The expected number of interior crossings is 9x9, plus 9
    // "touching" intersections along each of the left, right, and bottom edges. "epsilon" is used
    // to make the interior lines slightly longer so the "touches" actually cross, otherwise ~3 of
    // the 27 touches are not considered to intersect. However, the vertical lines do not reach the
    // top line as it curves on the surface of the sphere: despite epsilon those 9 are not even very
    // close to intersecting. Thus 9*12 = 108 interior and four more at the corners when
    // CrossingType.ALL is used.
    checkGetCrossingEdgePairs(index, CrossingType.ALL, 112);
    checkGetCrossingEdgePairs(index, CrossingType.INTERIOR, 108);
  }

  /**
   * Test crossing edges between two indexes, using the same setup as the previous test but with
   * horizontal lines in one index and vertical lines in the other.
   */
  public void testGetCrossingEdgePairsEdgeGridTwoIndexes() {
    int kGridSize = 10; // grid is zero to kGridSize inclusive
    double epsilon = 1e-10;

    S2ShapeIndex indexA = new S2ShapeIndex();
    S2ShapeIndex indexB = new S2ShapeIndex();
    S2EdgeVectorShape shapeA = new S2EdgeVectorShape();
    S2EdgeVectorShape shapeB = new S2EdgeVectorShape();
    for (int i = 0; i <= kGridSize; ++i) {
      double e = (i == 0 || i == kGridSize) ? 0 : epsilon;
      shapeA.add(
          S2LatLng.fromDegrees(-e, i).toPoint(), S2LatLng.fromDegrees(kGridSize + e, i).toPoint());
      shapeB.add(
          S2LatLng.fromDegrees(i, -e).toPoint(), S2LatLng.fromDegrees(i, kGridSize + e).toPoint());
    }
    indexA.add(shapeA);
    indexB.add(shapeB);

    // See comments on the previous test regarding the number of crossings.
    checkGetCrossingEdgePairs(indexA, indexB, CrossingType.ALL, 112);
    checkGetCrossingEdgePairs(indexA, indexB, CrossingType.INTERIOR, 108);
  }

  /**
   * Gets the crossing edge pairs in the given index, using {@link
   * S2VisitCrossingEdgePairs#visitCrossingEdgePairs()}.
   */
  private HashSet<ShapeEdgeIdPair> getCrossingEdgePairsOptimized(
      S2ShapeIndex index, CrossingType type) {
    HashSet<ShapeEdgeIdPair> edgePairs = new HashSet<>();
    S2CrossingEdgesQuery query = new S2CrossingEdgesQuery(type);
    query.visitCrossingEdgePairs(index,
        new EdgePairVisitor() {
          @Override
          public boolean visit(
              S2Shape shapeA, int edgeIdA, S2Point edgeASrc, S2Point edgeADst,
              S2Shape shapeB, int edgeIdB, S2Point edgeBSrc, S2Point edgeBDst,
              boolean isInterior) {
            edgePairs.add(new ShapeEdgeIdPair(shapeA, edgeIdA, shapeB, edgeIdB));
            return true; // Continue visiting.
          }
        });

    return edgePairs;
  }

  /**
   * Gets the crossing edge pairs between the two given shape indexes, using {@link
   * S2VisitCrossingEdgePairs#visitCrossingEdgePairs()}.
   */
  private HashSet<ShapeEdgeIdPair> getCrossingEdgePairsOptimized(
      S2ShapeIndex indexA, S2ShapeIndex indexB, CrossingType type) {
    HashSet<ShapeEdgeIdPair> edgePairs = new HashSet<>();
    S2CrossingEdgesQuery query = new S2CrossingEdgesQuery(type);
    query.visitCrossingEdgePairs(indexA, indexB,
        new EdgePairVisitor() {
          @Override
          public boolean visit(
              S2Shape shapeA, int edgeIdA, S2Point edgeASrc, S2Point edgeADst,
              S2Shape shapeB, int edgeIdB, S2Point edgeBSrc, S2Point edgeBDst,
              boolean isInterior) {
            edgePairs.add(new ShapeEdgeIdPair(shapeA, edgeIdA, shapeB, edgeIdB));
            return true; // Continue visiting.
          }
        });

    return edgePairs;
  }

  /** A simple container for a shape id and edge id. */
  private static final class ShapeIdEdgeId {
    public final int shapeId;
    public final int edgeId;

    public ShapeIdEdgeId(int shapeId, int edgeId) {
      this.shapeId = shapeId;
      this.edgeId = edgeId;
    }
  }

  /**
   * Returns the next (shapeId, edgeId) in the index after the given shapeId and edgeId, advancing
   * to the next edge of the given shape, or if the given edge is the last edge of the shape, to the
   * first edge of the next shape that has at least one edge (if any). If there are no more edges,
   * the returned shapeId will be equal to the number of shapes in the index.
   */
  private ShapeIdEdgeId next(S2ShapeIndex index, int shapeId, int edgeId) {
    Preconditions.checkArgument(shapeId < index.getShapes().size());
    int numEdges = index.getShapes().get(shapeId).numEdges();
    edgeId++;
    while (edgeId >= numEdges) {
      shapeId++;
      if (shapeId >= index.getShapes().size()) {
        break;
      }
      numEdges = index.getShapes().get(shapeId).numEdges();
      edgeId = 0;
    }
    return new ShapeIdEdgeId(shapeId, edgeId);
  }

  /**
   * Gets the crossing edge pairs in the given index using a brute force algorithm that tests all
   * pairs of edges.
   */
  private HashSet<ShapeEdgeIdPair> getCrossingEdgePairsBruteForce(
      S2ShapeIndex index, CrossingType type) {
    HashSet<ShapeEdgeIdPair> edgePairs = new HashSet<>();
    int minSign = (type == CrossingType.ALL) ? 0 : 1;

    List<S2Shape> shapes = index.getShapes();
    S2Shape.MutableEdge edgeA = new S2Shape.MutableEdge();
    S2Shape.MutableEdge edgeB = new S2Shape.MutableEdge();

    // Edge A iterates from the first edge of the first shape to the last edge of the last shape.
    for (int shapeIdA = 0; shapeIdA < shapes.size(); ++shapeIdA) {
      S2Shape shapeA = shapes.get(shapeIdA);
      for (int edgeIdA = 0; edgeIdA < shapeA.numEdges(); ++edgeIdA) {
        shapeA.getEdge(edgeIdA, edgeA);

        // Find the next edge in the index after (shapeIdA, edgeIdA).
        ShapeIdEdgeId nextShapeEdgeId = next(index, shapeIdA, edgeIdA);

        // Edge B starts on the next edge after Edge A, if there is one, and goes to the last edge
        // of the last shape.
        for (int shapeIdB = nextShapeEdgeId.shapeId; shapeIdB < shapes.size(); ++shapeIdB) {
          S2Shape shapeB = shapes.get(shapeIdB);
          for (int edgeIdB = nextShapeEdgeId.edgeId; edgeIdB < shapeB.numEdges(); ++edgeIdB) {
            shapeB.getEdge(edgeIdB, edgeB);

            // Do Edge A and Edge B have a crossing of the type we care about?
            int crossing = S2EdgeUtil.robustCrossing(edgeA.a, edgeA.b, edgeB.a, edgeB.b);
            if (crossing >= minSign) {
              edgePairs.add(new ShapeEdgeIdPair(shapeA, edgeIdA, shapeB, edgeIdB));
            }
          }
        }
      }
    }
    return edgePairs;
  }

  /**
   * Gets the crossing edge pairs between the two given shape indexes using a brute force algorithm
   * that tests all pairs of edges.
   */
  private HashSet<ShapeEdgeIdPair> getCrossingEdgePairsBruteForce(
      S2ShapeIndex indexA, S2ShapeIndex indexB, CrossingType type) {
    HashSet<ShapeEdgeIdPair> edgePairs = new HashSet<>();
    int minSign = (type == CrossingType.ALL) ? 0 : 1;

    List<S2Shape> shapesA = indexA.getShapes();
    List<S2Shape> shapesB = indexB.getShapes();
    S2Shape.MutableEdge edgeA = new S2Shape.MutableEdge();
    S2Shape.MutableEdge edgeB = new S2Shape.MutableEdge();

    // Edge A iterates over all the indexA shapes and edges.
    for (int shapeIdA = 0; shapeIdA < shapesA.size(); ++shapeIdA) {
      S2Shape shapeA = shapesA.get(shapeIdA);
      for (int edgeIdA = 0; edgeIdA < shapeA.numEdges(); ++edgeIdA) {
        shapeA.getEdge(edgeIdA, edgeA);

        // Edge B iterates over all the indexB shapes and edges.
        for (int shapeIdB = 0; shapeIdB < shapesB.size(); ++shapeIdB) {
          S2Shape shapeB = shapesB.get(shapeIdB);
          for (int edgeIdB = 0; edgeIdB < shapeB.numEdges(); ++edgeIdB) {
            shapeB.getEdge(edgeIdB, edgeB);

            // Do Edge A and Edge B have a crossing of the type we care about?
            int crossing = S2EdgeUtil.robustCrossing(edgeA.a, edgeA.b, edgeB.a, edgeB.b);
            if (crossing >= minSign) {
              edgePairs.add(new ShapeEdgeIdPair(shapeA, edgeIdA, shapeB, edgeIdB));
            }
          }
        }
      }
    }
    return edgePairs;
  }

  /**
   * Gets crossing edge pairs using both the brute force implementation above and the implementation
   * in {@link S2VisitCrossingEdgePairs#visitCrossingEdgePairs()}. Compares the results to ensure
   * they are the same, and checks that the number of crossings found was correct for both the brute
   * force and visitCrossingEdgePairs implementations.
   */
  private void checkGetCrossingEdgePairs(
      S2ShapeIndex index, CrossingType type, int expectedCrossingCount) {
    HashSet<ShapeEdgeIdPair> bruteForcePairs = getCrossingEdgePairsBruteForce(index, type);
    HashSet<ShapeEdgeIdPair> visitCrossingPairs = getCrossingEdgePairsOptimized(index, type);
    checkCrossingEdgePairs(expectedCrossingCount, bruteForcePairs, visitCrossingPairs);
  }

  /**
   * Gets crossing edge pairs between two shape indexes using both the brute force implementation
   * and the implementation in {@link S2VisitCrossingEdgePairs#visitCrossingEdgePairs()}. Compares
   * the results to ensure they are the same, and checks that the number of crossings found was
   * correct for both the brute force and visitCrossingEdgePairs implementations.
   */
  private void checkGetCrossingEdgePairs(
      S2ShapeIndex indexA, S2ShapeIndex indexB, CrossingType type, int expectedCrossingCount) {
    HashSet<ShapeEdgeIdPair> bruteForcePairs = getCrossingEdgePairsBruteForce(indexA, indexB, type);
    HashSet<ShapeEdgeIdPair> visitCrossingPairs =
        getCrossingEdgePairsOptimized(indexA, indexB, type);
    checkCrossingEdgePairs(expectedCrossingCount, bruteForcePairs, visitCrossingPairs);
  }

  /**
   * Compares the two sets of results to ensure they are the same, and checks that the number of
   * crossings found was correct for both.
   */
  private void checkCrossingEdgePairs(int expectedCrossingCount,
      HashSet<ShapeEdgeIdPair> bruteForcePairs,
      HashSet<ShapeEdgeIdPair> visitCrossingPairs) {
    boolean equal = true;
    for (ShapeEdgeIdPair edgePair : bruteForcePairs) {
      if (!visitCrossingPairs.contains(edgePair)) {
        System.err.println("Missing edge pair (in brute force, not in optimized): " + edgePair);
        equal = false;
      }
    }
    for (ShapeEdgeIdPair edgePair : visitCrossingPairs) {
      if (!bruteForcePairs.contains(edgePair)) {
        System.err.println("Unexpected edge pair (in optimized, not in brute force): " + edgePair);
        equal = false;
      }
    }

    assertEquals(
        "visitCrossingEdgePairs crossing count was not equal to expectedCrossingCount.",
        expectedCrossingCount,
        visitCrossingPairs.size());
    assertEquals(
        "Brute force crossing count was not equal to expectedCrossingCount.",
        expectedCrossingCount,
        bruteForcePairs.size());

    assertTrue(
        "visitCrossingPairs were not the same as brute force pairs; see details above. Expected "
            + bruteForcePairs.size()
            + " edge pairs, got "
            + visitCrossingPairs.size()
            + " edge pairs.",
        equal);
  }

  /**
   * Stores a (S2Shape, edgeId) pair, and defines equality and hashcode so that the order of the two
   * edges does not matter.
   */
  private static class ShapeEdgeIdPair {
    public final S2Shape shapeA;
    public final int edgeIdA;
    public final S2Shape shapeB;
    public final int edgeIdB;

    public ShapeEdgeIdPair(S2Shape shapeA, int edgeIdA, S2Shape shapeB, int edgeIdB) {
      this.shapeA = shapeA;
      this.edgeIdA = edgeIdA;
      this.shapeB = shapeB;
      this.edgeIdB = edgeIdB;
    }

    private String shapeEdgeToString(S2Shape shape, int edgeId) {
      S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
      shape.getEdge(edgeId, edge);
      return edge.getStart().toDegreesString() + " - " + edge.getEnd().toDegreesString();
    }

    @Override
    public String toString() {
      return "EdgeId A="
          + edgeIdA
          + ", EdgeId B="
          + edgeIdB
          + ", Edge A=("
          + shapeEdgeToString(shapeA, edgeIdA)
          + "), Edge B=("
          + shapeEdgeToString(shapeB, edgeIdB)
          + ")\n";
    }

    @Override
    public int hashCode() {
      return (shapeA.hashCode() + shapeB.hashCode()) * 7 + (edgeIdA + edgeIdB);
    }

    @Override
    public boolean equals(Object other) {
      if (this == other) {
        return true;
      }
      if (!(other instanceof ShapeEdgeIdPair)) {
        return false;
      }
      ShapeEdgeIdPair that = (ShapeEdgeIdPair) other;
      return (this.shapeA == that.shapeA
              && this.edgeIdA == that.edgeIdA
              && this.shapeB == that.shapeB
              && this.edgeIdB == that.edgeIdB)
          || (this.shapeA == that.shapeB
              && this.edgeIdA == that.edgeIdB
              && this.shapeB == that.shapeA
              && this.edgeIdB == that.edgeIdA);
    }
  }
}
