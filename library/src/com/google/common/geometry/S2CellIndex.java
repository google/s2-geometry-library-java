/*
 * Copyright 2019 Google Inc.
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
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * S2CellIndex stores a collection of (cellId, label) pairs. The S2CellIds may be overlapping or
 * contain duplicate values. For example, an S2CellIndex could store a collection of S2CellUnions,
 * where each S2CellUnion has its own label. Labels are 32-bit non-negative integers, and are
 * typically used to map the results of queries back to client data structures.
 *
 * <p>To build an S2CellIndex, call add() for each (cellId, label) pair, and then call the build()
 * method. For example:
 *
 * {@snippet :
 *  List<S2CellId> contents = ...;
 *  for (int i = 0; i < contents.size(); ++i) {
 *    index.add(contents.get(i), i);  // i is used as the label.
 *  }
 *  index.build();
 * }
 *
 * <p>There is also a convenience method to add all the cells in a union, associated with one label.
 *
 * <p>Note that the index is not dynamic; the contents of the index cannot be changed once it has
 * been built. However, the same memory space can be reused by {@link #clear()} clearing the index.
 *
 * <p>There are several options for retrieving data from the index. The simplest is to use a
 * built-in method such as getIntersectingLabels, which returns the labels of all cells that
 * intersect a given target S2CellUnion. Alternatively, you can access the index contents directly.
 *
 * <p>Internally, the index consists of a set of non-overlapping leaf cell ranges that subdivide the
 * sphere and such that each range intersects a particular set of (cellId, label) pairs. Data is
 * accessed using the following iterator types:
 *
 * <ul>
 *   <li>{@link RangeIterator}: used to seek and iterate over the non-overlapping leaf cell ranges.
 *   <li>{@link NonEmptyRangeIterator}: like RangeIterator, but skips ranges whose contents are
 *       empty.
 *   <li>{@link ContentsIterator}: iterates over the (cellId, label) pairs that intersect a given
 *       range.
 *   <li>{@link CellIterator}: iterates over the entire set of (cellId, label) pairs.
 * </ul>
 *
 * <p>Note that these are low-level, efficient types intended mainly for implementing new query
 * classes. Most clients should use either the built-in methods such as {@link
 * #visitIntersectingCells(S2CellUnion, CellVisitor)} and {@link
 * #getIntersectingLabels(S2CellUnion)}.
 */
public class S2CellIndex {
  /**
   * A tree of (cellId, label) pairs such that if X is an ancestor of Y, then X.cellId contains
   * Y.cellId. The contents of a given range of leaf cells can be represented by pointing to a node
   * of this tree.
   */
  private final List<CellNode> cellNodes = new ArrayList<>();

  /**
   * The last element of rangeNodes is a sentinel value, which is necessary in order to represent
   * the range covered by the previous element.
   */
  private final ArrayList<RangeNode> rangeNodes = new ArrayList<>();

  /** Returns the number of (cellId, label) pairs in the index. */
  public int numCells() {
    return cellNodes.size();
  }

  /**
   * Adds the given (cellId, label) pair to the index. Note that the index is not valid until
   * {@link #build} is called.
   *
   * <p>The S2CellIds in the index may overlap (including duplicate values). Duplicate (cellId,
   * label) pairs are also allowed, although query tools often remove duplicates.
   *
   * <p>Results are undefined unless all cells are {@link S2CellId#isValid()} valid.
   */
  public void add(S2CellId cellId, int label) {
    assert cellId.isValid();
    assert label >= 0;
    cellNodes.add(new CellNode(cellId, label, -1));
  }

  /** Convenience function that adds a collection of cells with the same label. */
  public void add(Iterable<S2CellId> cellIds, int label) {
    for (S2CellId cellId : cellIds) {
      add(cellId, label);
    }
  }

  /**
   * Builds the index. This method may only be called once. No iterators may be used until the index
   * is built.
   */
  public void build() {
    // Create two deltas for each (cellId, label) pair: one to add the pair to the stack (at the
    // start of its leaf cell range), and one to remove it from the stack (at the end of its leaf
    // cell range).
    Delta[] deltas = new Delta[2 * cellNodes.size() + 2];
    int size = 0;
    for (CellNode node : cellNodes) {
      deltas[size++] = new Delta(node.cellId.rangeMin(), node.cellId, node.label);
      deltas[size++] = new Delta(node.cellId.rangeMax().next(), S2CellId.sentinel(), -1);
    }

    // We also create two special deltas to ensure that a RangeNode is emitted at the beginning and
    // end of the S2CellId range.
    deltas[size++] = new Delta(S2CellId.begin(S2CellId.MAX_LEVEL), S2CellId.none(), -1);
    deltas[size++] = new Delta(S2CellId.end(S2CellId.MAX_LEVEL), S2CellId.none(), -1);
    Arrays.sort(deltas, Delta.BY_START_CELL_NEG_LABEL);

    // Now walk through the deltas to build the leaf cell ranges and cell tree (which is essentially
    // a permanent form of the "stack" described above).
    cellNodes.clear();
    rangeNodes.ensureCapacity(deltas.length);
    int contents = -1;
    for (int i = 0; i < deltas.length; ) {
      // Process all the deltas associated with the current startId.
      S2CellId startId = deltas[i].startId;
      for (; i < deltas.length && deltas[i].startId.equals(startId); i++) {
        if (deltas[i].label >= 0) {
          cellNodes.add(new CellNode(deltas[i].cellId, deltas[i].label, contents));
          contents = cellNodes.size() - 1;
        } else if (deltas[i].cellId.equals(S2CellId.sentinel())) {
          contents = cellNodes.get(contents).parent;
        }
      }
      rangeNodes.add(new RangeNode(startId, contents));
    }
  }

  /** Returns an iterator over the cells of this index. */
  public CellIterator cells() {
    Preconditions.checkState(!rangeNodes.isEmpty(), "Call build() first.");
    return new CellIterator();
  }

  /** Returns an iterator over the ranges of this index. */
  public RangeIterator ranges() {
    Preconditions.checkState(!rangeNodes.isEmpty(), "Call build() first.");
    return new RangeIterator();
  }

  /** Returns an iterator over the non-empty ranges of this index. */
  public NonEmptyRangeIterator nonEmptyRanges() {
    return new NonEmptyRangeIterator();
  }

  /** Returns an iterator over the contents of this index. */
  public ContentsIterator contents() {
    Preconditions.checkState(!rangeNodes.isEmpty(), "Call build() first.");
    return new ContentsIterator();
  }

  /**
   * To build the cell tree and leaf cell ranges, we maintain a stack of (cellId, label) pairs that
   * contain the current leaf cell. This class represents an instruction to push or pop a (cellId,
   * label) pair.
   *
   * <p>If label >= 0, the (cellId, label) pair is pushed on the stack. If cellId ==
   * S2CellId.SENTINEL, a pair is popped from the stack. Otherwise the stack is unchanged but a
   * RangeNode is still emitted.
   */
  private static final class Delta {
    /**
     * Deltas are sorted first by startId, then in reverse order by cellId, and then by label. This
     * is necessary to ensure that (1) larger cells are pushed on the stack before smaller cells,
     * and (2) cells are popped off the stack before any new cells are added.
     */
    public static final Comparator<Delta> BY_START_CELL_NEG_LABEL =
        (a, b) -> {
          int result = a.startId.compareTo(b.startId);
          if (result != 0) {
            return result;
          }
          result = -a.cellId.compareTo(b.cellId);
          if (result != 0) {
            return result;
          }
          return Integer.compare(a.label, b.label);
        };

    private final S2CellId startId;
    private final S2CellId cellId;
    private final int label;

    private Delta(S2CellId startId, S2CellId cellId, int label) {
      this.startId = startId;
      this.cellId = cellId;
      this.label = label;
    }
  }

  /** Clears the index so that it can be re-used. */
  public void clear() {
    cellNodes.clear();
    rangeNodes.clear();
  }

  /**
   * Visits all (cellId, label) pairs in the given index that intersect the given S2CellUnion
   * "target" and returns true, or terminates early and returns false if {@link
   * CellVisitor#visit(S2CellId, int)} ever returns false.
   *
   * <p>Each (cellId, label) pair in the index is visited at most once. If the index contains
   * duplicates, then each copy is visited.
   */
  @CanIgnoreReturnValue
  public boolean visitIntersectingCells(S2CellUnion target, CellVisitor visitor) {
    if (target.size() == 0) {
      return true;
    }

    ContentsIterator contents = contents();
    RangeIterator range = ranges();
    for (int i = 0; i < target.size(); ) {
      S2CellId id = target.cellId(i);

      // Only seek the range to this target cell when necessary.
      if (range.limitId().lessOrEquals(id.rangeMin())) {
        range.seek(id.rangeMin());
      }

      // Visit contents of this range that intersect this cell.
      for (; range.startId().lessOrEquals(id.rangeMax()); range.next()) {
        for (contents.startUnion(range); !contents.done(); contents.next()) {
          if (!visitor.visit(contents.cellId(), contents.label())) {
            return false;
          }
        }
      }

      // Check whether the next target cell is also contained by the leaf cell range that we just
      // processed. If so, we can skip over all such cells using binary search. This speeds up
      // benchmarks by 2-10x when the average number of intersecting cells is small (< 1).
      i++;
      if (i != target.size() && target.cellId(i).rangeMax().lessThan(range.startId())) {
        // Skip to the first target cell that extends past the previous range.
        i = S2ShapeUtil.lowerBound(i + 1, target.size(),
            j -> range.startId().greaterThan(target.cellId(j)));
        if (target.cellId(i - 1).rangeMax().greaterOrEquals(range.startId())) {
          i--;
        }
      }
    }

    return true;
  }

  /** Returns the distinct sorted labels that intersect the given target. */
  public Labels getIntersectingLabels(S2CellUnion target) {
    Labels result = new Labels();
    getIntersectingLabels(target, result);
    result.normalize();
    return result;
  }

  /** Appends labels intersecting 'target', in unspecified order, with possible duplicates. */
  public void getIntersectingLabels(S2CellUnion target, Labels results) {
    visitIntersectingCells(target, (cellId, label) -> results.add(label));
  }

  /** A function that is called with each (cellId, label) pair to be visited. */
  public interface CellVisitor {
    /** Provides a (cellId, label) pair to this visitor, which may return true to keep searching. */
    boolean visit(S2CellId cellId, int label);
  }

  /**
   * A set of labels that can be grown by {@link #getIntersectingLabels(S2CellUnion, Labels)} and
   * shrunk via {@link #clear()} or {@link #normalize()}. May contain duplicates or be unsorted
   * unless {@link #normalize()} is called.
   */
  public static class Labels extends AbstractList<Integer> {
    private int[] labels = new int[8];
    private int size;

    @Override
    public void clear() {
      size = 0;
    }

    private boolean add(int label) {
      if (labels.length == size) {
        labels = Arrays.copyOf(labels, size * 2);
      }
      labels[size++] = label;
      return true;
    }

    @Override
    public int size() {
      return size;
    }

    @Override
    @SuppressWarnings("unusable-by-js")
    public Integer get(int index) {
      return getInt(index);
    }

    /** As {@link #get(int)} but without the overhead of boxing. */
    public int getInt(int index) {
      return labels[index];
    }

    /** Sorts the labels and removes duplicates. */
    public void normalize() {
      if (size == 0) {
        return;
      }
      Arrays.sort(labels, 0, size);
      int lastIndex = 0;
      for (int i = 1; i < size; i++) {
        if (labels[lastIndex] != labels[i]) {
          labels[++lastIndex] = labels[i];
        }
      }
      size = lastIndex + 1;
    }
  }

  /**
   * Represents a node in the (cellId, label) tree. Cells are organized in a tree such that the
   * ancestors of a given node contain that node.
   */
  private static final class CellNode {
    private S2CellId cellId;
    private int label;
    private int parent;

    private CellNode(S2CellId cellId, int label, int parent) {
      this.cellId = cellId;
      this.label = label;
      this.parent = parent;
    }

    private void setFrom(CellNode node) {
      this.cellId = node.cellId;
      this.label = node.label;
      this.parent = node.parent;
    }
  }

  /** An iterator over all (cellId, label) pairs in an unspecified order. */
  public final class CellIterator {
    /** Offset into {@link S2CellIndex#cellNodes}. */
    private int offset;
    /** Current node pointed to by 'offset', or null if {@link #done()}. */
    private CellNode cell;

    // Initializes a CellIterator for the S2CellIndex, positioned at the first cell (if any).
    private CellIterator() {
      Preconditions.checkState(!rangeNodes.isEmpty(), "Call build() first.");
      seek(0);
    }

    /** Returns the S2CellId of the current (cellId, label) pair. */
    public S2CellId cellId() {
      assert !done();
      return cell.cellId;
    }

    /** Returns the label of the current (cellId, label) pair. */
    public int label() {
      assert !done();
      return cell.label;
    }

    /** Returns true if all (cellId, label) pairs have been visited. */
    public boolean done() {
      return offset == cellNodes.size();
    }

    /** Advances this iterator to the next (cellId, label) pair. */
    public void next() {
      assert !done();
      seek(offset + 1);
    }

    /** Sets the offset and sets 'cell' accordingly. */
    private void seek(int offset) {
      this.offset = offset;
      this.cell = done() ? null : cellNodes.get(offset);
    }
  }

  /**
   * A RangeNode represents a range of leaf S2CellIds. The range starts at startId (a leaf cell) and
   * ends at the startId field of the next RangeNode. "contents" points to the node of cellNodes
   * representing the cells that overlap this range.
   */
  private static class RangeNode {
    /** First leaf cell contained by this range. */
    private final S2CellId startId;
    /** Index in {@link S2CellIndex#cellNodes} for the cells that overlap this range. */
    private final int contents;

    private RangeNode(S2CellId startId, int contents) {
      this.startId = startId;
      this.contents = contents;
    }
  }

  /**
   * An iterator that seeks and iterates over a set of non-overlapping leaf cell ranges that cover
   * the entire sphere. The indexed (S2CellId, Label) pairs that intersect the current leaf cell
   * range can be visited using ContentsIterator (see below).
   */
  public class RangeIterator {
    /** Offset into {@link S2CellIndex#rangeNodes}. */
    private int offset;
    /** Current node pointed to by 'offset'. */
    private RangeNode node = rangeNodes.get(offset);

    /**
     * Returns the start of the current range of leaf S2CellIds. When {@link #done()}, this returns
     * the {@link S2CellId#end(int)} of the {@link S2CellId#MAX_LEVEL} max level of cells, so that
     * most loops may test this method instead of done().
     */
    public S2CellId startId() {
      return node.startId;
    }

    /** The (non-inclusive) end of the current range of leaf S2CellIds. */
    public S2CellId limitId() {
      assert (!done());
      return rangeNodes.get(offset + 1).startId;
    }

    /** Returns true if the iterator is positioned beyond the last valid range. */
    public boolean done() {
      // Note that the last element of rangeNodes is a sentinel value.
      return offset >= rangeNodes.size() - 1;
    }

    /** Positions this iterator at the first range of leaf cells (if any). */
    public void begin() {
      seekAndLoad(0);
    }

    /** Positions the iterator so that done() is true. */
    public void finish() {
      // Note that the last element of rangeNodes is a sentinel value.
      seekAndLoad(rangeNodes.size() - 1);
    }

    /** Advances the iterator to the next range of leaf cells. */
    public void next() {
      assert (!done());
      seekAndLoad(offset + 1);
    }

    /**
     * Returns false if the iterator was already positioned at the beginning, otherwise positions
     * the iterator at the previous entry and returns true.
     */
    public boolean prev() {
      if (offset == 0) {
        return false;
      }
      seekAndLoad(offset - 1);
      return true;
    }

    /**
     * Positions the iterator at the range containing "target". Such a range exists as long as the
     * target is a valid leaf cell.
     *
     * @param target a valid leaf (level 30) cell to seek to
     */
    public void seek(S2CellId target) {
      assert target.isLeaf();
      seekAndLoad(S2ShapeUtil.upperBound(0, rangeNodes.size(),
          i -> target.lessThan(rangeNodes.get(i).startId)) - 1);
    }

    /** Returns true if no (S2CellId, Label) pairs intersect this range, or if {@link #done()}. */
    public boolean isEmpty() {
      return node.contents == ContentsIterator.DONE;
    }

    /**
     * Advances this iterator 'n' times and returns true, or if doing so would advance this iterator
     * past the end, leaves the iterator unmodified and returns false.
     */
    public boolean advance(int n) {
      // Note that the last element of rangeNodes is a sentinel value.
      if (n >= rangeNodes.size() - 1 - offset) {
        return false;
      }
      seekAndLoad(offset + n);
      return true;
    }

    private void seekAndLoad(int offset) {
      this.offset = offset;
      this.node = rangeNodes.get(offset);
    }
  }

  /** As {@link RangeIterator} but only visits range nodes that overlap (cellId, label) pairs. */
  public class NonEmptyRangeIterator extends RangeIterator {

    /** Positions the iterator at the first non-empty range of leaf cells. */
    @Override
    public void begin() {
      super.begin();
      while (isEmpty() && !done()) {
        super.next();
      }
    }

    /** Advances the iterator to the next non-empty range of leaf cells. */
    @Override
    public void next() {
      do {
        super.next();
      } while (isEmpty() && !done());
    }

    /**
     * If the iterator is already positioned at the beginning, returns false. Otherwise positions
     * the iterator at the previous non-empty entry and returns true.
     */
    @Override
    public boolean prev() {
      while (super.prev()) {
        if (!isEmpty()) {
          return true;
        }
      }
      // Return the iterator to its original position.
      if (isEmpty() && !done()) {
        next();
      }
      return false;
    }

    // Positions the iterator at the range that contains or follows "target", or at the end if no
    // such range exists. (Note that start_id() may still be called in the latter case.)
    @Override
    public void seek(S2CellId target) {
      super.seek(target);
      while (isEmpty() && !done()) {
        super.next();
      }
    }
  }

  /**
   * An iterator that visits the (cellId, label) pairs that cover a set of leaf cell ranges (see
   * RangeIterator). To use it, construct an instance or {@link #clear()} an existing instance, and
   * {@link #startUnion(RangeIterator)} to visit the contents of each desired leaf cell range.
   *
   * <p>Note that when multiple leaf cell ranges are visited, this class only guarantees that each
   * result will be reported at least once, i.e. duplicate values may be suppressed. If you want
   * duplicate values to be reported again, be sure to call {@link #clear()} first.
   *
   * <p>In particular, the implementation guarantees that when multiple leaf cell ranges are visited
   * in monotonically increasing order, then each (cellId, label) pair is reported exactly once.
   */
  public class ContentsIterator {
    /** A special label indicating that {@link #done} is true. */
    private static final int DONE = -1;

    /**
     * The value of it.startId() from the previous call to startUnion(). This is used to check
     * whether these values are monotonically increasing.
     */
    private S2CellId prevStartId;

    /**
     * The maximum index within {@link #cellNodes} visited during the previous call to startUnion().
     * This is used to eliminate duplicate values when startUnion() is called multiple times.
     */
    private int nodeCutoff;

    /**
     * The maximum index within {@link #cellNodes} visited during the current call to startUnion().
     * This is used to update nodeCutoff.
     */
    private int nextNodeCutoff;

    /** A copy of the current node in the cell tree. */
    private final CellNode node = new CellNode(null, DONE, -1);

    /** Creates a new iterator. Call {@link #startUnion(RangeIterator)} next. */
    private ContentsIterator() {
      clear();
    }

    /** Clears all state with respect to which range(s) have been visited. */
    public void clear() {
      prevStartId = S2CellId.none();
      nodeCutoff = -1;
      nextNodeCutoff = -1;
      setDone();
    }

    /**
     * Positions the ContentsIterator at the first (cellId, label) pair that covers the given leaf
     * cell range. Note that when multiple leaf cell ranges are visited using the same
     * ContentsIterator, duplicate values may be suppressed. If you don't want this behavior, call
     * clear() first.
     */
    public void startUnion(RangeIterator range) {
      if (range.startId().lessThan(prevStartId)) {
        // Can't automatically eliminate duplicates.
        nodeCutoff = -1;
      }
      prevStartId = range.startId();
      int contents = range.node.contents;
      if (contents <= nodeCutoff) {
        setDone();
      } else {
        node.setFrom(cellNodes.get(contents));
      }

      // When visiting ancestors, we can stop as soon as the node index is smaller than any
      // previously visited node index. Because indexes are assigned using a preorder traversal,
      // such nodes are guaranteed to have already been reported.
      nextNodeCutoff = contents;
    }

    /** Returns the S2CellId of the current (cellId, label) pair. */
    public S2CellId cellId() {
      assert !done();
      return node.cellId;
    }

    /** Returns the label of the current (cellId, label) pair. */
    public int label() {
      assert !done();
      return node.label;
    }

    /** Returns true if all (cellId, label) pairs have been visited. */
    public boolean done() {
      return node.label == DONE;
    }

    /**
     * Advances the iterator to the next (cellId, label) pair covered by the current leaf cell
     * range.
     */
    public void next() {
      assert !done();
      if (node.parent <= nodeCutoff) {
        // We have already processed this node and its ancestors.
        nodeCutoff = nextNodeCutoff;
        setDone();
      } else {
        node.setFrom(cellNodes.get(node.parent));
      }
    }

    /** Sets the current node label to DONE to indicate that iteration has finished. */
    private void setDone() {
      node.label = DONE;
    }
  }
}
