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
package com.google.common.geometry;

import static com.google.common.geometry.S2CellRangeIterator.makeS2CellRangeIterator;

import com.google.common.geometry.S2ShapeIndex.CellRelation;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.List;
import java.util.function.BiPredicate;
import java.util.function.Predicate;

/**
 * An inner join operation on any two S2Iterators. S2CellIteratorJoin takes an optional distance
 * value which specifies a "buffer" around cells in the iterator. If cells are within that distance,
 * then they're considered to overlap. This allows us to support "tolerant" versions of queries.
 * Currently the tolerance must not be negative.
 *
 * <p>The iterator will find each pair of cells that are less than or equal to the tolerance
 * distance from each other and call a visitor with a reference to each of the properly positioned
 * iterators. The visitor can then process the overlapping iterators however it wishes, returning
 * true if it wishes to continue iterating, and false otherwise.
 *
 * <p>Example usage:
 *
 * {@snippet :
 * // Join two S2ShapeIndex iterators.
 * // Process cell pairs until we find two cells with different edges.
 *
 * S2ShapeIndex indexA = ...;
 * S2ShapeIndex indexB = ...;
 * makeS2CellIteratorJoin(indexA, indexB, distance).join(
 *     (iterA, iterB) -> allEdgesAreEqual(iterA.entry(), iterB.entry()));
 * }
 *
 * But we're not limited to joining iterators of the same type. Any implementations of the {@code
 * S2Iterator<>} API will work. For example, we can join an S2ShapeIndex with an S2PointIndex:
 *
 * {@snippet :
 * boolean processCellPoints(S2Iterator<Cell> iterA, S2Iterator<Entry<String>> iterB) {
 *   // Process shapes and points from overlapping cells.
 *   return true;
 * }
 *
 * S2ShapeIndex indexA = ...;
 * S2PointIndex<String> indexB = ...;
 *
 * makeS2CellIteratorJoin(indexA, indexB).join(processCellPoints);
 * }
 *
 * <p>We use a visitor pattern because it allows us to avoid materializing intermediate join results
 * before processing, which is important for joins because the number of results can be very large.
 * Also, as opposed to something like the next()/entry()/done() pattern that S2Iterator itself uses,
 * we can pass the visitor to a different class or function to actually implement the join,
 * invisibly to the user.
 *
 * <p>In addition to allowing us to delegate iteration nicely, the visitor pattern lets us separate
 * concerns when writing spatial algorithms. The actual processing of indexed data is separated from
 * the logic of positioning the iterators, which allows for much cleaner implementations.
 */
public class S2CellIteratorJoin<E1 extends S2Iterator.Entry, E2 extends S2Iterator.Entry> {
  /** Wraps the S2Iterator "iteratorA" given to the constructor. */
  private final S2CellRangeIterator<E1> iterA;

  /** Wraps the S2Iterator "iteratorB" given to the constructor. */
  private final S2CellRangeIterator<E2> iterB;

  /** Cells within this distance are considered to overlap. */
  private final S1ChordAngle tolerance;

  /** Reusable storage for S2Cells on one side of the join. */
  private final List<S2Cell> matchedCells = new ArrayList<>();

  /** Reused in processCellPairs. */
  int idx;

  /** Reusable array. New array entries are created each time the array is used. */
  S2Cell[] reusableChildCells = new S2Cell[4];

  /** Constructor that takes two S2Iterators and uses a tolerance of zero. */
  public S2CellIteratorJoin(S2Iterator<E1> iteratorA, S2Iterator<E2> iteratorB) {
    this(iteratorA, iteratorB, S1ChordAngle.ZERO);
  }

  /** Constructor that takes two S2Iterators and a specified tolerance. */
  public S2CellIteratorJoin(
      S2Iterator<E1> iteratorA, S2Iterator<E2> iteratorB, S1ChordAngle tolerance) {
    this.iterA = makeS2CellRangeIterator(iteratorA);
    this.iterB = makeS2CellRangeIterator(iteratorB);
    this.tolerance = tolerance;
  }

  /**
   * Executes the join. Returns false if the visitor ever does, true otherwise. The provided visitor
   * is called with each pair of cells that are less than or equal to the tolerance distance from
   * each other. It should return true if the iteration should continue, and false otherwise.
   */
  @CanIgnoreReturnValue
  public <V extends BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>> boolean join(
      V visitor) {
    if (tolerance.isZero()) {
      return exactJoin(visitor);
    } else {
      return tolerantJoin(visitor);
    }
  }

  /** Conducts a simple join, sending left results to the given consumer while it returns true. */
  @CanIgnoreReturnValue
  public boolean simpleJoin(Predicate<S2CellRangeIterator<E1>> leftResults) {
    return join(
        new BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>() {
          S2CellId lastId = S2CellId.SENTINEL;

          @Override
          public boolean test(S2CellRangeIterator<E1> left, S2CellRangeIterator<E2> right) {
            if (lastId.equals(left.id())) {
              return true;
            }
            lastId = left.id();
            return leftResults.test(left);
          }
        });
  }

  /** Performs an exact inner join (when the tolerance is zero). */
  private <V extends BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>>
      boolean exactJoin(V visitor) {
    iterA.begin();
    iterB.begin();

    // Iterate until we hit the end of an iterator or visitor tells us to stop.
    while (!iterA.done() && !iterB.done()) {
      int order = iterA.relation(iterB);
      switch (order) {
        case -1:
          // A precedes B, seek A.
          iterA.seekTo(iterB);
          break;

        case +1:
          // B precedes A, seek B.
          iterB.seekTo(iterA);
          break;

        case 0:
          // Iterators overlap.
          // If visitor rejects pair, then we're done.
          if (!visitor.test(iterA, iterB)) {
            return false;
          }

          // Move the smaller of the cells forward.
          long lsbA = iterA.id().lowestOnBit();
          long lsbB = iterB.id().lowestOnBit();
          int cmp = UnsignedLongs.compare(lsbA, lsbB);
          switch (cmp) {
            case -1: // lsbA < lsbB
              iterA.next();
              break;
            case 1: // lsbA > lsbB
              iterB.next();
              break;
            case 0: //  Cells are the same size and overlap, they must be the same. Advance both.
              iterA.next();
              iterB.next();
              break;
            default:
              throw new IllegalStateException("Unexpected compareUnsigned result: " + cmp);
          }
          break;

        default:
          throw new IllegalStateException("Unexpected order: " + order);
      }
    }
    return true;
  }

  // ---- Tolerant join related code.

  // Maximum number of cross-terms before we recurse.
  static final int MAX_CROSS_PRODUCT = 25;
  static final int COVER_LIMIT = MAX_CROSS_PRODUCT / 2;

  /**
   * Initialize the given covering with cells that cover the range of the given
   * S2CellRangeIterator's current position.
   */
  private void coverCurrentPosition(S2CellRangeIterator<?> iter, S2CellUnion covering) {
    iter.restart();
    covering.clear();
    if (iter.done()) {
      return;
    }

    S2CellId min = iter.rangeMin();
    iter.finish();
    if (!iter.prev()) {
      return; // Empty iterator.
    }
    S2CellId max = iter.rangeMax();
    covering.initFromMinMax(min, max);
  }

  /** Does a tolerant join (when the tolerance is non-zero). */
  private <V extends BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>>
      boolean tolerantJoin(V visitor) {
    // Seed the recursion with a coarse covering of each input iterator.
    S2CellUnion coveringA = new S2CellUnion();
    S2CellUnion coveringB = new S2CellUnion();
    coverCurrentPosition(iterA, coveringA);
    coverCurrentPosition(iterB, coveringB);

    List<S2Cell> nearbyCells = new ArrayList<>();
    for (S2CellId cellIdA : coveringA.cellIds()) {
      nearbyCells.clear();
      S2Cell cellA = new S2Cell(cellIdA);

      for (S2CellId cellIdB : coveringB.cellIds()) {
        S2Cell cellB = new S2Cell(cellIdB);
        if (cellA.isDistanceLessOrEqual(cellB, tolerance)) {
          nearbyCells.add(cellB);
        }
      }

      if (!nearbyCells.isEmpty()) {
        if (!processCellPairs(cellA, nearbyCells, visitor)) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Given an S2Cell 'cellA' and a list of S2Cells 'cellsB', filters the cells from B that are
   * within the tolerance distance distance of cellA A and passes them to processCellPairs().
   * Returns false if the visitor ever does, true otherwise.
   */
  private <V extends BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>>
      boolean processNearby(S2Cell cellA, List<S2Cell> cellsB, V visitor) {
    // TODO(torrey): Avoid allocating a new list here every time this is recursively called.
    List<S2Cell> nearbyCells = new ArrayList<>();

    for (S2Cell cellB : cellsB) {
      if (cellA.isDistanceLessOrEqual(cellB, tolerance)) {
        nearbyCells.add(cellB);
      }
    }

    if (!nearbyCells.isEmpty()) {
      if (!processCellPairs(cellA, nearbyCells, visitor)) {
        return false;
      }
    }
    return true;
  }

  /** As above but subdivides the given 'cellA' before processing it. */
  private <V extends BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>>
      boolean processNearbySubdivided(S2Cell cellA, List<S2Cell> cellsB, V visitor) {
    // TODO(torrey): Avoid allocating a new array here every time this is recursively called.
    S2Cell[] childCells = new S2Cell[4];
    for (int i = 0; i < 4; ++i) {
      childCells[i] = new S2Cell();
    }
    cellA.subdivide(childCells);
    for (S2Cell childA : childCells) {
      if (!processNearby(childA, cellsB, visitor)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Processes a pair of cells known to be within tolerance of each other. If the portion of each
   * index that's covered by the cells is small enough, then we report pairs to the visitor,
   * otherwise the cells are subdivided and we recurse.
   *
   * <p>Since there are only thirty levels to the cell hierarchy, this recursion is safe as we'll
   * never go more than 30 stack frames deep.
   */
  private <V extends BiPredicate<S2CellRangeIterator<E1>, S2CellRangeIterator<E2>>>
      boolean processCellPairs(S2Cell cellA, List<S2Cell> cellsB, V visitor) {
    // Estimate how many index cells the A cell covers.
    int numCoveredA = estimateCoveredCells(iterA, cellA.id());
    if (numCoveredA == 0) {
      return true;
    }

    // Scan the cells of the B union. Prune any cells that don't cover any of the index, and
    // subdivide any that cover too much.
    boolean subdivided = false;
    List<S2Cell> subdividedB = new ArrayList<>();
    for (S2Cell cellB : cellsB) {
      int numCoveredB = estimateCoveredCells(iterB, cellB.id());
      if (numCoveredB == 0) {
        continue;
      } else if (numCoveredB < COVER_LIMIT) {
        subdividedB.add(cellB);
      } else {
        // Subdivide cellB and add its children to the list of cells to process.
        appendChildren(cellB, subdividedB);
        subdivided = true;
      }
    }

    // If A covers too many index cells or we had to subdivide B, then continue the recursion by
    // pairing the A cells and B cells that are within the tolerance.
    if (numCoveredA >= COVER_LIMIT || subdivided) {
      // If the A cell covers too many index cells, subdivide it.
      return (numCoveredA >= COVER_LIMIT)
          ? processNearbySubdivided(cellA, subdividedB, visitor)
          : processNearby(cellA, subdividedB, visitor);
    }

    // Otherwise A and the B union are small enough we can pair them up and report nearby pairs to
    // the visitor now.

    // Pre-compute the B S2Cells to avoid replicating work in the inner loop.
    matchedCells.clear();
    for (S2Cell cellB : cellsB) {
      scanCellRange(
          iterB,
          cellB.id(),
          iter -> {
            matchedCells.add(new S2Cell(iter.id()));
            return true;
          });
    }

    return scanCellRange(
        iterA,
        cellA.id(),
        iterA -> {
          // If the index cell doesn't have it's endpoint in this cell then ignore it.
          // This makes it so that we only see each index cell in a single A cell.
          // TODO(user): This conditional is almost certainly not needed any more, as the
          // cell iterator join is driven from side A, and only ever subdivides the B side. If
          // unit test coverage of S2CellIteratorJoin clients confirms this, remove it.
          if (!cellA.id().intersects(iterA.id().rangeMin())) {
            return true;
          }

          S2Cell subCellA = new S2Cell(iterA.id());

          // For each sub cell in the cellA range, scan each B cell in the cell union and send the
          // resulting (A,B) pairs to the visitor. If it returns false at any point then stop the
          // join and return false.

          idx = 0;
          boolean success = true;
          for (int i = 0; i < cellsB.size() && success; ++i) {
            S2CellId id = cellsB.get(i).id();
            success &=
                scanCellRange(
                    iterB,
                    id,
                    iterB -> {
                      if (subCellA.isDistanceLessOrEqual(matchedCells.get(idx++), tolerance)) {
                        if (!visitor.test(iterA, iterB)) {
                          return false;
                        }
                      }
                      return true;
                    });
          }
          return success;
        });
  }

  /** Subdivides the given cell and appends the children to the given list of cells. */
  private void appendChildren(S2Cell cell, List<S2Cell> cells) {
    for (int i = 0; i < 4; ++i) {
      reusableChildCells[i] = new S2Cell();
    }
    cell.subdivide(reusableChildCells);
    for (int i = 0; i < 4; ++i) {
      cells.add(reusableChildCells[i]);
    }
  }

  /**
   * Positions an S2CellRangeIterator and visits every position that intersects a given S2CellId.
   * The iterator is passed to the visitor at each position. Returns false if the visitor ever does,
   * otherwise true.
   */
  @CanIgnoreReturnValue
  private static <E extends S2Iterator.Entry, V extends Predicate<S2CellRangeIterator<E>>>
      boolean scanCellRange(S2CellRangeIterator<E> iter, S2CellId id, V visitor) {
    CellRelation unused = iter.locate(id);
    for (; !iter.done() && iter.id().intersects(id); iter.next()) {
      if (!visitor.test(iter)) {
        return false;
      }
    }
    return true;
  }

  /** Estimates how many values are covered by a given S2CellId. */
  private static int estimateCoveredCells(S2Iterator<?> iter, S2CellId cell) {
    CellRelation relation = iter.locate(cell);
    switch (relation) {
      case DISJOINT:
        return 0;

      case INDEXED:
        return 1;

      case SUBDIVIDED:
        // Note: If we develop a Distance method for Iterators that can tell us how many cells are
        // within a range efficiently, then we can avoid the linear scan here that just has us
        // repeatedly aborting.
        S2CellId end = cell.rangeMax();
        int matches = 0;
        for (; !iter.done() && iter.compareTo(end) <= 0; iter.next()) {
          matches++;

          // There's too many matches, so give up. This is a heuristic that will help us decide
          // whether to recurse on the left or right.
          if (matches > COVER_LIMIT) {
            return COVER_LIMIT;
          }
        }
        return matches;
    }
    throw new IllegalStateException("Unexpected CellRelation: " + relation);
  }
}
