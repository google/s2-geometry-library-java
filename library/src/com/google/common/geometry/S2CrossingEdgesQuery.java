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
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.S2ShapeUtil.RangeIterator;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A query for visiting pairs of crossing edges in one S2ShapeIndex, or pairs of crossing edges
 * between two S2ShapeIndexes where one edge is from each index.
 *
 * <p>See also {@link S2EdgeQuery}, which finds edges or shapes that are crossed by a single edge at
 * a time. If you have one or just a few edges that aren't in an index, and want to find their
 * crossings of indexed edges or shapes, S2EdgeQuery may be faster than S2CrossingEdgesQuery.
 *
 * <p>Otherwise, for finding crossings in bulk between all edges in an index or between two indexes,
 * this S2CrossingEdgesQuery will be much faster.
 */
public class S2CrossingEdgesQuery {
  /**
   * What crossings should be returned? INTERIOR or ALL. For INTERIOR, only crossings between
   * endpoints are returned. For ALL, edges that merely touch at endpoints are also considered
   * crossings.
   */
  private final CrossingType crossingType;
  /** minCrossingSign is 1 if crossingType is INTERIOR, 0 if it is ALL. */
  private final int minCrossingSign;

  /**
   * If "needAdjacent" is false, then edge pairs of the form (AB, BC) may optionally be ignored,
   * even if the two edges belong to different edge chains.
   */
  private final boolean needAdjacent;

  private final IndexCrosser abCrosser = new IndexCrosser();
  private final IndexCrosser baCrosser = new IndexCrosser();
  private final EdgeCrosser crosser = new EdgeCrosser();

  /** The S2Shape, plus the source, destination, and ids of clipped edges in clipped shape 'A'. */
  private final LoadedShape a = new LoadedShape();
  /** The S2Shape, plus the source, destination, and ids of clipped edges in clipped shape 'B'. */
  private final LoadedShape b = new LoadedShape();

  /**
   * Constructs a new S2CrossingEdgesQuery for the given CrossingType.
   *
   * @param type indicates whether all crossings should be visited, or only interior crossings.
   */
  public S2CrossingEdgesQuery(CrossingType type) {
    crossingType = type;
    needAdjacent = (type == CrossingType.ALL);
    minCrossingSign = (type == CrossingType.INTERIOR) ? 1 : 0;
  }

  /**
   * Visits all pairs of crossing edges in the given S2ShapeIndex, terminating early if the given
   * {@link EdgePairVisitor} returns false, in which case visitCrossings returns false as well.
   *
   * <p>CAVEAT: Crossings may be visited more than once.
   */
  @CanIgnoreReturnValue
  public boolean visitCrossingEdgePairs(S2ShapeIndex index, EdgePairVisitor visitor) {
    return visitCrossings(index, visitor);
  }

  /**
   * Visits all pairs of crossing edges where one edge comes from each of the given S2ShapeIndexes,
   * terminating early if the given {@link EdgePairVisitor} returns false, in which case
   * visitCrossingEdgePairs returns false as well.
   *
   * <p>CAVEAT: Crossings may be visited more than once.
   */
  @CanIgnoreReturnValue
  public boolean visitCrossingEdgePairs(
      S2ShapeIndex aIndex, S2ShapeIndex bIndex, EdgePairVisitor visitor) {
    // We look for S2CellId ranges where the indexes of A and B overlap, and then test those edges
    // for crossings.
    // TODO(user): Use brute force if the total number of edges is small enough (using a
    // larger threshold if the S2ShapeIndex is not constructed yet).
    RangeIterator<Cell> ai = new RangeIterator<>(aIndex.iterator());
    RangeIterator<Cell> bi = new RangeIterator<>(bIndex.iterator());
    // Tests A against B
    abCrosser.init(aIndex, bIndex, crossingType, visitor, false);
    // Tests B against A
    baCrosser.init(bIndex, aIndex, crossingType, visitor, true);

    // TODO(user): Replace the following with S2CellIteratorJoin when that is available.
    while (!ai.done() || !bi.done()) {
      if (ai.rangeMax().lessThan(bi.rangeMin())) {
        // The A and B cells don't overlap, and A precedes B.
        ai.seekTo(bi);
      } else if (bi.rangeMax().lessThan(ai.rangeMin())) {
        // The A and B cells don't overlap, and B precedes A.
        bi.seekTo(ai);
      } else {
        // One cell contains the other. Determine which cell is larger.
        long abRelation = ai.id().lowestOnBit() - bi.id().lowestOnBit();
        if (abRelation > 0) {
          // A's index cell is larger.
          if (!abCrosser.visitCrossings(ai, bi)) {
            return false;
          }
        } else if (abRelation < 0) {
          // B's index cell is larger.
          if (!baCrosser.visitCrossings(bi, ai)) {
            return false;
          }
        } else {
          // The A and B cells are the same.
          if (ai.iterator().entry().numEdges() > 0 && bi.iterator().entry().numEdges() > 0) {
            if (!abCrosser.visitCellCellCrossings(ai.iterator().entry(), bi.iterator().entry())) {
              return false;
            }
          }
          ai.next();
          bi.next();
        }
      }
    }
    return true;
  }

  /**
   * Visits all pairs of crossing edges in the given S2ShapeIndex, terminating early if the given
   * {@link EdgePairVisitor} returns false (in which case visitCrossings returns false as well).
   */
  private boolean visitCrossings(S2ShapeIndex index, EdgePairVisitor visitor) {
    // TODO(user): Use brute force if the total number of edges is small enough (using a
    // larger threshold if the S2ShapeIndex is not constructed yet).
    for (S2Iterator<Cell> it = index.iterator(); !it.done(); it.next()) {
      if (!visitCellCrossings(it.entry(), visitor)) {
        return false;
      }
    }
    return true;
  }

  /** Given a S2ShapeIndex.Cell, visit all pairs of crossing edges in the cell. */
  private boolean visitCellCrossings(Cell cell, EdgePairVisitor visitor) {
    // In the outer loop, iterate over and load every clipped shape in the cell.
    for (int aShapeNum = 0; aShapeNum < cell.numShapes(); aShapeNum++) {
      a.load(cell.clipped(aShapeNum));

      // Iterate over and load clipped shapes starting at the current shape in the outer loop.
      for (int bShapeNum = aShapeNum; bShapeNum < cell.numShapes(); bShapeNum++) {
        b.load(cell.clipped(bShapeNum));

        // Cross edges of 'a' and 'b', noting that 'a' and 'b' are the same shape when starting.
        for (int ai = 0; ai < a.size(); ai++) {
          S2Point aSrc = a.srcs.get(ai);
          S2Point aDst = a.dsts.get(ai);
          int aId = a.ids[ai];

          crosser.init(aSrc, aDst);
          // If 'a' and 'b' are the same, start iterating B on the edge after A.
          for (int bi = (aShapeNum == bShapeNum) ? (ai + 1) : 0; bi < b.size(); bi++) {
            S2Point bSrc = b.srcs.get(bi);
            S2Point bDst = b.dsts.get(bi);
            int bId = b.ids[bi];

            // A common situation is that an edge AB is followed by an edge BC. We only need to
            // visit such crossings if "needAdjacent" is true (even if AB and BC belong to different
            // edge chains). Otherwise, skip that edge.
            if (!needAdjacent && aDst.equalsPoint(bSrc)) {
              continue;
            }

            // Does edge A cross edge B?
            int sign = crosser.robustCrossing(bSrc, bDst);
            if (sign >= minCrossingSign) {
              if (!visitor.visit(
                  a.shape, aId, aSrc, aDst,
                  b.shape, bId, bSrc, bDst,
                  sign == 1)) {
                return false;
              }
            }
          }
        }
      }
    }
    return true;
  }

  /** A parameter that controls the reporting of edge intersections. */
  public enum CrossingType {
    /**
     * CrossingType.INTERIOR reports intersections that occur at a point interior to both edges
     * (i.e., not at a vertex).
     */
    INTERIOR,
    /**
     * CrossingType.ALL reports all intersections, even those where two edges intersect only because
     * they share a common vertex.
     */
    ALL
  };

  /**
   * A visitor that is called with pairs of crossing edges. The predicate may return false in order
   * to request that the algorithm should be terminated, i.e. no further crossings are needed.
   */
  public interface EdgePairVisitor {
    /**
     * visit() is called for each pair of crossing edges until it returns false or no more crossings
     * are found.
     *
     * @param isInterior indicates that the crossing is at a point interior to both edges (i.e., not
     *     at a vertex). The calling function already has this information and it is moderately
     *     expensive to recompute.
     */
    boolean visit(
        S2Shape aShape, int aEdgeId, S2Point aEdgeSrc, S2Point aEdgeDst,
        S2Shape bShape, int bEdgeId, S2Point bEdgeSrc, S2Point bEdgeDst,
        boolean isInterior);
  }

  /** A reusable container for loading and visiting the clipped edges of a clipped shape. */
  static class LoadedShape {
    /** The source points of the clipped shape edges. */
    public final ArrayList<S2Point> srcs = new ArrayList<>();

    /** The destination points of the clipped shape edges. */
    public final ArrayList<S2Point> dsts = new ArrayList<>();

    /**
     * The edge ids of the clipped shape edges. To avoid reallocation, the array may be larger than
     * required. The number of elements actually used is equal to srcs.size() (and dsts.size()).
     */
    public int[] ids;

    /** A MutableEdge used to load edge end points. */
    private final MutableEdge edge = new MutableEdge();

    /** The S2Shape underlying the currently clipped shape. */
    public S2Shape shape;

    /** Constructs a LoadedShape. */
    public LoadedShape() {}

    /** Loads the edges of the given clipped shape, replacing the current content. */
    public void load(S2ClippedShape clipped) {
      shape = clipped.shape();
      srcs.clear();
      dsts.clear();

      // Only (re)allocate the ids array if required.
      if (ids == null) {
        ids = new int[clipped.numEdges()];
      } else if (ids.length < clipped.numEdges()) {
        ids = Arrays.copyOf(ids, clipped.numEdges());
      }

      for (int i = 0; i < clipped.numEdges(); i++) {
        int id = clipped.edge(i);
        clipped.shape().getEdge(id, edge);
        ids[i] = id;
        srcs.add(edge.getStart());
        dsts.add(edge.getEnd());
      }
    }

    /** How many clipped edges does this LoadedShape currently contain? */
    public int size() {
      return srcs.size();
    }

    /**
     * Visits the loaded clipped edges until the visitor returns false, in which case false is
     * returned, or all edges have been visited, in which case true is returned.
     */
    public boolean visit(Visitor visitor) {
      for (int i = 0; i < srcs.size(); i++) {
        if (!visitor.visit(shape, ids[i], srcs.get(i), dsts.get(i))) {
          return false;
        }
      }
      return true;
    }

    /**
     * A visitor for the contents of a LoadedShape. May return false to indicate that no more
     * clipped edges should be visited.
     */
    public interface Visitor {
      boolean visit(S2Shape shape, int edgeId, S2Point src, S2Point dst);
    }
  }


  /**
   * IndexCrosser is a helper class for finding edge crossings between a pair of indexes with {@link
   * #visitCrossingEdgePairs(S2ShapeIndex, S2ShapeIndex, EdgePairVisitor)}. It is reusable by
   * calling init() repeatedly. Two instances are used by visitCrossingEdgePairs(), one for the
   * index pair (A, B) and one for the index pair (B, A), in order to be able to test edge crossings
   * in the most efficient order.
   */
  private static class IndexCrosser {
    // These fields are set during init().
    private EdgePairVisitor visitor;
    private int minCrossingSign;
    private boolean swapped;
    private S2EdgeQuery bQuery;

    // Reused collections and objects.
    private final EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser();
    private final ArrayList<Cell> bCells = new ArrayList<>();

    /** The source, destination, and ids of clipped edges in clipped shape 'A'. */
    private final LoadedShape a = new LoadedShape();
    /** The source, destination, and ids of clipped edges in clipped shape 'B'. */
    private final LoadedShape b = new LoadedShape();

    /** Constructs an IndexCrosser. */
    public IndexCrosser() {}

    /**
     * Initializes this IndexCrosser to find edge crossings between "aIndex" and "bIndex". If
     * "swapped" is true, the shape indexes A and B have been swapped. This affects how arguments
     * are passed to the visitor, since for example A.contains(B) is not the same as B.contains(A).
     */
    public void init(
        S2ShapeIndex aIndex,
        S2ShapeIndex bIndex,
        CrossingType type,
        EdgePairVisitor visitor,
        boolean swapped) {
      this.visitor = visitor;
      minCrossingSign = (type == CrossingType.INTERIOR) ? 1 : 0;
      this.swapped = swapped;
      bQuery = new S2EdgeQuery(bIndex);
    }

    /** Calls the visitor with the given edge information, swapping the order if needed. */
    private boolean visitEdgePair(
        S2Shape aShape, int aEdgeId, S2Point aEdgeSrc, S2Point aEdgeDst,
        S2Shape bShape, int bEdgeId, S2Point bEdgeSrc, S2Point bEdgeDst,
        boolean isInterior) {
      if (swapped) {
        return visitor.visit(
            bShape, bEdgeId, bEdgeSrc, bEdgeDst, aShape, aEdgeId, aEdgeSrc, aEdgeDst, isInterior);
      } else {
        return visitor.visit(
            aShape, aEdgeId, aEdgeSrc, aEdgeDst, bShape, bEdgeId, bEdgeSrc, bEdgeDst, isInterior);
      }
    }

    /**
     * Visits the crossings of all the edges of the current LoadedShape 'a' vs. all the edges of
     * LoadedShape 'b'.
     */
    private boolean visitLoadedShapeCrossings() {
      for (int ai = 0; ai < a.srcs.size(); ai++) {
        int aId = a.ids[ai];
        S2Point aSrc = a.srcs.get(ai);
        S2Point aDst = a.dsts.get(ai);
        crosser.init(aSrc, aDst);

        for (int bi = 0; bi < b.srcs.size(); bi++) {
          int bId = b.ids[bi];
          S2Point bSrc = b.srcs.get(bi);
          S2Point bDst = b.dsts.get(bi);

          int sign = crosser.robustCrossing(bSrc, bDst);
          if (sign >= minCrossingSign) {
            if (!visitEdgePair(
                a.shape, aId, aSrc, aDst,
                b.shape, bId, bSrc, bDst,
                sign == 1)) {
              return false;
            }
          }
        }
      }
      return true;
    }

    /**
     * Given two iterators positioned such that ai.id().contains(bi.id()), visits all crossings
     * between edges of A and B that intersect a.id(). Terminates early and returns false if the
     * "visitor" provided to the constructor returns false.
     *
     * <p>Advances both iterators past ai.id().
     */
    @CanIgnoreReturnValue
    public boolean visitCrossings(RangeIterator<Cell> ai, RangeIterator<Cell> bi) {
      Preconditions.checkArgument(ai.id().contains(bi.id()));

      // If the cell at 'ai' has no edges, simply advance 'bi' past it.
      if (ai.entry().numEdges() == 0) {
        // Skip over the cells of B using binary search.
        bi.seekBeyond(ai);
      } else {
        // If ai.id() intersects many edges of B, then it is faster to use S2EdgeQuery to narrow
        // down the candidates. But if it intersects only a few edges, it is faster to check all
        // the crossings directly. We handle this by advancing "bi" and keeping track of how many
        // edges we would need to test.
        // TODO(torrey): Set this appropriately for Java with a benchmark. "23" is what C++ uses.
        final int kEdgeQueryMinEdges = 23;
        int bEdges = 0;
        bCells.clear();
        do {
          int cellEdges = bi.entry().numEdges();
          if (cellEdges > 0) {
            bEdges += cellEdges;
            if (bEdges >= kEdgeQueryMinEdges) {
              // There are too many edges, so use an S2EdgeQuery.
              if (!visitSubcellCrossings(ai.entry(), ai.id())) {
                return false;
              }
              bi.seekBeyond(ai);
              return true;
            }
            bCells.add(bi.entry());
          }
          bi.next();
        } while (bi.id().lessOrEquals(ai.rangeMax()));

        if (!bCells.isEmpty()) {
          // Test all the edge crossings directly.
          if (!visitCellCellsCrossings(ai.entry(), bCells)) {
            return false;
          }
        }
      }
      ai.next();
      return true;
    }

    /**
     * Visits all crossings of any edge in "aCell" with any index cell of B that is a descendant of
     * "bId". Terminates early and returns false if the "visitor" supplied to the constructor
     * returns false.
     */
    private boolean visitSubcellCrossings(S2ShapeIndex.Cell aCell, S2CellId bId) {
      final S2PaddedCell bRoot = new S2PaddedCell(bId, 0);
      // Test all shape edges in "aCell" against the edges contained in B index cells that are
      // descendants of "bId".
      for (int aShapeNum = 0; aShapeNum < aCell.numShapes(); aShapeNum++) {
        S2ClippedShape aClippedShape = aCell.clipped(aShapeNum);
        a.load(aClippedShape);
        if (!a.visit(
            (shape, edgeId, src, dst) -> {
              bCells.clear();
              bQuery.getCells(src, dst, bRoot, bCells);
              for (Cell bCell : bCells) {
                if (!visitEdgeCellCrossings(shape, edgeId, src, dst, bCell)) {
                  return false;
                }
              }
              return true;
            })) {
          return false;
        }
      }
      return true;
    }

    /** Visit edge crossings of edge A with any edge in bCell. */
    private boolean visitEdgeCellCrossings(
        S2Shape aShape, int aId, S2Point aSrc, S2Point aDst, S2ShapeIndex.Cell bCell) {
      crosser.init(aSrc, aDst);
      for (int bShapeNum = 0; bShapeNum < bCell.numShapes(); bShapeNum++) {
        S2ClippedShape bClippedShape = bCell.clipped(bShapeNum);
        b.load(bClippedShape);
        if (!b.visit(
            (bShape, bId, bSrc, bDst) -> {
              int sign = crosser.robustCrossing(bSrc, bDst);
              if (sign >= minCrossingSign) {
                if (!visitEdgePair(aShape, aId, aSrc, aDst, bShape, bId, bSrc, bDst, sign == 1)) {
                  return false;
                }
              }
              return true;
            })) {
          return false;
        }
      }
      return true;
    }

    /**
     * Given a S2ShapeIndex.Cell "aCell" and a list of S2ShapeIndex.Cells "bCells", visit all
     * crossings between edges of the aCell and the bCells. Terminates early and returns false if
     * the "visitor" supplied to the constructor returns false. Note that edge pairs may be visited
     * more than once.
     */
    private boolean visitCellCellsCrossings(
        S2ShapeIndex.Cell aCell, List<S2ShapeIndex.Cell> bCells) {
      for (int aShapeNum = 0; aShapeNum < aCell.numShapes(); aShapeNum++) {
        a.load(aCell.clipped(aShapeNum));
        for (Cell bCell : bCells) {
          for (int bShapeNum = 0; bShapeNum < bCell.numShapes(); bShapeNum++) {
            b.load(bCell.clipped(bShapeNum));
            if (!visitLoadedShapeCrossings()) {
              return false;
            }
          }
        }
      }
      return true;
    }

    /**
     * Given two index cells, visits all crossings between edges of those cells. Terminates early
     * and returns false if the "visitor" supplied to the constructor returns false.
     */
    private boolean visitCellCellCrossings(S2ShapeIndex.Cell aCell, S2ShapeIndex.Cell bCell) {
      for (int aShapeNum = 0; aShapeNum < aCell.numShapes(); aShapeNum++) {
        a.load(aCell.clipped(aShapeNum));
        for (int bShapeNum = 0; bShapeNum < bCell.numShapes(); bShapeNum++) {
          b.load(bCell.clipped(bShapeNum));
          if (!visitLoadedShapeCrossings()) {
            return false;
          }
        }
      }
      return true;
    }
  }
}
