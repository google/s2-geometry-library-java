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

import static com.google.common.geometry.S2Projections.PROJ;

import com.google.common.annotations.GwtCompatible;
import com.google.common.collect.Lists;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * An S2CellUnion is a region consisting of cells of various sizes. Typically a cell union is used
 * to approximate some other shape. There is a tradeoff between the accuracy of the approximation
 * and how many cells are used. Unlike polygons, cells have a fixed hierarchical structure. This
 * makes them more suitable for optimizations based on preprocessing.
 *
 * <p>An S2CellUnion is represented as a vector of sorted, non-overlapping S2CellIds. By default the
 * vector is also "normalized", meaning that groups of 4 child cells have been replaced by their
 * parent cell whenever possible. S2CellUnions are not required to be normalized, but certain
 * operations will return different results if they are not, e.g. {@link #contains(S2CellUnion)}.
 *
 */
@GwtCompatible(serializable = true)
public strictfp class S2CellUnion implements S2Region, Iterable<S2CellId>, Serializable {
  private static final long serialVersionUID = 1L;

  private static final byte LOSSLESS_ENCODING_VERSION = 1;

  /** The CellIds that form the Union */
  private ArrayList<S2CellId> cellIds = new ArrayList<S2CellId>();

  public S2CellUnion() {}

  /**
   * Populates a cell union with the given S2CellIds, and then calls normalize(). This directly uses
   * the input list, without copying it.
   */
  public void initFromCellIds(ArrayList<S2CellId> cellIds) {
    initRawCellIds(cellIds);
    normalize();
  }

  /** Populates a cell union with the given 64-bit cell ids, and then calls normalize(). */
  public void initFromIds(List<Long> cellIds) {
    initRawIds(cellIds);
    normalize();
  }

  /**
   * Populates a cell union with the given S2CellIds. The input list is copied, and then cleared.
   */
  public void initSwap(List<S2CellId> cellIds) {
    initRawSwap(cellIds);
    normalize();
  }

  /**
   * Populates a cell union with the given S2CellIds. This does not call normalize, see {@link
   * #initRawSwap} for details. This directly uses the input list, without copying it.
   */
  public void initRawCellIds(ArrayList<S2CellId> cellIds) {
    this.cellIds = cellIds;
  }

  /**
   * Populates a cell union with the given 64 bit cell ids. This does not call normalize, see {@link
   * #initRawSwap} for details. The input list is copied.
   */
  // TODO(eengle): Make a constructed S2CellUnion immutable, and port other init methods from C++.
  public void initRawIds(List<Long> cellIds) {
    int size = cellIds.size();
    this.cellIds = new ArrayList<S2CellId>(size);
    for (Long id : cellIds) {
      this.cellIds.add(new S2CellId(id));
    }
  }

  /**
   * Like the initFrom*() constructors, but does not call normalize(). The cell union *must* be
   * normalized before doing any calculations with it, so it is the caller's * responsibility to
   * make sure that the input is normalized. This method is useful when converting cell unions to
   * another representation and back.
   *
   * <p>The input list is copied, and then cleared.
   */
  public void initRawSwap(List<S2CellId> cellIds) {
    this.cellIds = new ArrayList<S2CellId>(cellIds);
    cellIds.clear();
  }

  /**
   * Create a cell union that corresponds to a continuous range of cell ids. The output is a
   * normalized collection of cell ids that covers the leaf cells between "minId" and "maxId"
   * inclusive.
   *
   * <p>Requires that minId.isLeaf(), maxId.isLeaf(), and minId <= maxId.
   */
  public void initFromMinMax(S2CellId minId, S2CellId maxId) {
    // assert minId.isLeaf();
    // assert maxId.isLeaf();
    // assert minId.compareTo(maxId) <= 0;
    // assert minId.isValid() && maxId.isValid();
    initFromBeginEnd(minId, maxId.next());
  }

  /**
   * As {@link #initFromMinMax(S2CellId, S2CellId)}, except that the union covers the range of leaf
   * cells from "begin" (inclusive) to "end" (exclusive.) If {@code begin.equals(end)}, the result
   * is empty.
   *
   * <p>Requires that begin.isLeaf(), end.isLeaf(), and begin <= end.
   */
  public void initFromBeginEnd(S2CellId begin, S2CellId end) {
    // assert (begin.isLeaf());
    // assert (end.isLeaf());
    // assert (begin.compareTo(end) <= 0);

    // We repeatedly add the largest cell we can, in sorted order.
    cellIds.clear();
    for (S2CellId nextBegin = begin; nextBegin.compareTo(end) < 0; ) {
      // assert(nextBegin.isLeaf());

      // Find the largest cell that starts at "next_begin" and ends before "end".
      S2CellId nextId = nextBegin;
      while (!nextId.isFace()
          && nextId.parent().rangeMin().equals(nextBegin)
          && nextId.parent().rangeMax().compareTo(end) < 0) {
        nextId = nextId.parent();
      }
      cellIds.add(nextId);
      nextBegin = nextId.rangeMax().next();
    }

    // The output should already be sorted and normalized.
    // assert(!normalize());
  }

  public int size() {
    return cellIds.size();
  }

  /** Convenience methods for accessing the individual cell ids. */
  public S2CellId cellId(int i) {
    return cellIds.get(i);
  }

  /** Enable iteration over the union's cells. */
  @Override
  public Iterator<S2CellId> iterator() {
    return cellIds.iterator();
  }

  /** Direct access to the underlying vector for iteration . */
  public ArrayList<S2CellId> cellIds() {
    return cellIds;
  }

  /**
   * Returns true if the cell union is valid, meaning that the S2CellIds are non-overlapping and
   * sorted in increasing order.
   */
  public boolean isValid() {
    for (int i = 1; i < cellIds.size(); i++) {
      if (cellIds.get(i - 1).rangeMax().compareTo(cellIds.get(i).rangeMin()) >= 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if the cell union is normalized, meaning that it {@link #isValid()} is true and
   * that no four cells have a common parent.
   *
   * <p>Certain operations such as {@link #contains(S2CellUnion)} may return a different result if
   * the cell union is not normalized.
   */
  public boolean isNormalized() {
    for (int i = 1; i < cellIds.size(); i++) {
      if (cellIds.get(i - 1).rangeMax().compareTo(cellIds.get(i).rangeMin()) >= 0) {
        return false;
      }
      if (i >= 3
          && areSiblings(
              cellIds.get(i - 3), cellIds.get(i - 2),
              cellIds.get(i - 1), cellIds.get(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if the given four cells have a common parent.
   *
   * <p>Requires the four cells are distinct.
   */
  private static boolean areSiblings(S2CellId a, S2CellId b, S2CellId c, S2CellId d) {
    // A necessary (but not sufficient) condition is that the XOR of the four cells must be zero.
    // This is also very fast to test.
    if ((a.id() ^ b.id() ^ c.id()) != d.id()) {
      return false;
    }

    // Now we do a slightly more expensive but exact test.  First, compute a mask that blocks out
    // the two bits that encode the child position of "id" with respect to its parent, then check
    // that the other three children all agree with "mask".
    long mask = d.lowestOnBit() << 1;
    mask = ~(mask + (mask << 1));
    long idMasked = d.id() & mask;
    return !d.isFace()
        && (a.id() & mask) == idMasked
        && (b.id() & mask) == idMasked
        && (c.id() & mask) == idMasked;
  }

  /**
   * Replaces "output" with an expanded version of the cell union where any cells whose level is
   * less than "min_level" or where (level - min_level) is not a multiple of "level_mod" are
   * replaced by their children, until either both of these conditions are satisfied or the maximum
   * level is reached.
   *
   * <p>This method allows a covering generated by S2RegionCoverer using min_level() or level_mod()
   * constraints to be stored as a normalized cell union (which allows various geometric
   * computations to be done) and then converted back to the original list of cell ids that
   * satisfies the desired constraints.
   */
  public void denormalize(int minLevel, int levelMod, ArrayList<S2CellId> output) {
    // assert (minLevel >= 0 && minLevel <= S2CellId.MAX_LEVEL);
    // assert (levelMod >= 1 && levelMod <= 3);

    output.clear();
    output.ensureCapacity(size());
    for (S2CellId id : this) {
      int level = id.level();
      int newLevel = Math.max(minLevel, level);
      if (levelMod > 1) {
        // Round up so that (new_level - min_level) is a multiple of level_mod.
        // (Note that S2CellId::kMaxLevel is a multiple of 1, 2, and 3.)
        newLevel += (S2CellId.MAX_LEVEL - (newLevel - minLevel)) % levelMod;
        newLevel = Math.min(S2CellId.MAX_LEVEL, newLevel);
      }
      if (newLevel == level) {
        output.add(id);
      } else {
        S2CellId end = id.childEnd(newLevel);
        for (id = id.childBegin(newLevel); !id.equals(end); id = id.next()) {
          output.add(id);
        }
      }
    }
  }

  /**
   * If there are more than "excess" elements of the cell_ids() vector that are allocated but
   * unused, reallocate the array to eliminate the excess space. This reduces memory usage when many
   * cell unions need to be held in memory at once.
   */
  public void pack() {
    cellIds.trimToSize();
  }

  /**
   * Return true if the cell union contains the given cell id. Containment is defined with respect
   * to regions, e.g. a cell contains its 4 children. This is a fast operation (logarithmic in the
   * size of the cell union).
   *
   * <p>CAVEAT: If you have constructed a valid but non-normalized S2CellUnion, note that groups of
   * 4 child cells are <em>not</em> considered to contain their parent cell. To get this behavior
   * you must construct a normalized cell union, or call {@link #normalize()} prior to this method.
   */
  public boolean contains(S2CellId id) {
    // This is an exact test. Each cell occupies a linear span of the S2
    // space-filling curve, and the cell id is simply the position at the center
    // of this span. The cell union ids are sorted in increasing order along
    // the space-filling curve. So we simply find the pair of cell ids that
    // surround the given cell id (using binary search). There is containment
    // if and only if one of these two cell ids contains this cell.

    int pos = Collections.binarySearch(cellIds, id);
    if (pos < 0) {
      pos = -pos - 1;
    }
    if (pos < cellIds.size() && cellIds.get(pos).rangeMin().lessOrEquals(id)) {
      return true;
    }
    return pos != 0 && cellIds.get(pos - 1).rangeMax().greaterOrEquals(id);
  }

  /**
   * Return true if the cell union intersects the given cell id. This is a fast operation
   * (logarithmic in the size of the cell union).
   */
  public boolean intersects(S2CellId id) {
    // This is an exact test; see the comments for Contains() above.
    int pos = Collections.binarySearch(cellIds, id);

    if (pos < 0) {
      pos = -pos - 1;
    }

    if (pos < cellIds.size() && cellIds.get(pos).rangeMin().lessOrEquals(id.rangeMax())) {
      return true;
    }
    return pos != 0 && cellIds.get(pos - 1).rangeMax().greaterOrEquals(id.rangeMin());
  }

  /**
   * Returns true if this cell union contains {@code that}.
   *
   * <p>CAVEAT: If you have constructed a valid but non-normalized S2CellUnion, note that groups of
   * 4 child cells are <em>not</em> considered to contain their parent cell. To get this behavior
   * you must construct a normalized cell union, or call {@link #normalize()} prior to this method.
   */
  public boolean contains(S2CellUnion that) {
    S2CellUnion result = new S2CellUnion();
    result.getIntersection(this, that);
    return result.cellIds.equals(that.cellIds);
  }

  /** This is a fast operation (logarithmic in the size of the cell union). */
  @Override
  public boolean contains(S2Cell cell) {
    return contains(cell.id());
  }

  /** Return true if this cell union intersects {@code union}. */
  public boolean intersects(S2CellUnion union) {
    S2CellUnion result = new S2CellUnion();
    result.getIntersection(this, union);
    return result.size() > 0;
  }

  /** Sets this cell union to the union of {@code x} and {@code y}. */
  public void getUnion(S2CellUnion x, S2CellUnion y) {
    // assert (x != this && y != this);
    cellIds.clear();
    cellIds.ensureCapacity(x.size() + y.size());
    cellIds.addAll(x.cellIds);
    cellIds.addAll(y.cellIds);
    normalize();
  }

  /**
   * Specialized version of GetIntersection() that gets the intersection of a cell union with the
   * given cell id. This can be useful for "splitting" a cell union into chunks.
   *
   * <p><b>Note:</b> {@code x} and {@code y} must be normalized.
   */
  public void getIntersection(S2CellUnion x, S2CellId id) {
    // assert (x != this);
    cellIds.clear();
    if (x.contains(id)) {
      cellIds.add(id);
    } else {
      int pos = Collections.binarySearch(x.cellIds, id.rangeMin());

      if (pos < 0) {
        pos = -pos - 1;
      }

      S2CellId idmax = id.rangeMax();
      int size = x.cellIds.size();
      while (pos < size && x.cellIds.get(pos).lessOrEquals(idmax)) {
        cellIds.add(x.cellIds.get(pos++));
      }
    }
    // assert isNormalized() || !x.isNormalized();
  }

  /**
   * Initializes this cell union to the intersection of the two given cell unions. Requires: x !=
   * this and y != this.
   *
   * <p><b>Note:</b> {@code x} and {@code y} must be normalized.
   */
  public void getIntersection(S2CellUnion x, S2CellUnion y) {
    getIntersection(x.cellIds, y.cellIds, cellIds);
    // The output is normalized as long as at least one input is normalized.
    // assert isNormalized() || (!x.isNormalized() && !y.isNormalized());
  }

  /**
   * Like {@code #getIntersection(S2CellUnion, S2CellUnion)}, but works directly with lists of
   * S2CellIds, and this method has slightly more relaxed normalization requirements: the input
   * vectors may contain groups of 4 child cells that all have the same parent. (In a normalized
   * S2CellUnion, such groups are always replaced by the parent cell.)
   *
   * <p><b>Note:</b> {@code x} and {@code y} must be sorted.
   */
  public static void getIntersection(List<S2CellId> x, List<S2CellId> y, List<S2CellId> results) {
    // assert (x != results && y != results);

    // This is a fairly efficient calculation that uses binary search to skip
    // over sections of both input vectors. It takes constant time if all the
    // cells of "x" come before or after all the cells of "y" in S2CellId order.
    results.clear();
    int i = 0;
    int j = 0;
    while (i < x.size() && j < y.size()) {
      S2CellId xCell = x.get(i);
      S2CellId xMin = xCell.rangeMin();
      S2CellId yCell = y.get(j);
      S2CellId yMin = yCell.rangeMin();
      if (xMin.greaterThan(yMin)) {
        // Either j->contains(xCell) or the two cells are disjoint.
        if (xCell.lessOrEquals(yCell.rangeMax())) {
          results.add(xCell);
          i++;
        } else {
          // Advance "j" to the first cell possibly contained by xCell.
          j = indexedBinarySearch(y, xMin, j + 1);
          // The previous cell *(j-1) may now contain xCell.
          if (xCell.lessOrEquals(y.get(j - 1).rangeMax())) {
            --j;
          }
        }
      } else if (yMin.greaterThan(xMin)) {
        // Identical to the code above with "i" and "j" reversed.
        if (yCell.lessOrEquals(xCell.rangeMax())) {
          results.add(yCell);
          j++;
        } else {
          i = indexedBinarySearch(x, yMin, i + 1);
          if (yCell.lessOrEquals(x.get(i - 1).rangeMax())) {
            --i;
          }
        }
      } else {
        // "i" and "j" have the same rangeMin(), so one contains the other.
        if (xCell.lessThan(yCell)) {
          results.add(xCell);
          i++;
        } else {
          results.add(yCell);
          j++;
        }
      }
    }
  }

  /** Initiaizes this cell union to the difference of the two given cell unions. */
  public void getDifference(S2CellUnion x, S2CellUnion y) {
    // TODO(user): this is approximately O(N*log(N)), but could probably use similar techniques as
    // getIntersection() to be more efficient.
    cellIds.clear();
    for (S2CellId id : x) {
      getDifferenceInternal(id, y);
    }
    // The output is normalized as long as the first argument is normalized.
    // assert isNormalized() || !x.isNormalized();
  }

  private void getDifferenceInternal(S2CellId cell, S2CellUnion y) {
    // Add the difference between cell and y to cellIds. If they intersect but the difference is
    // non-empty, divide and conquer.
    if (!y.intersects(cell)) {
      cellIds.add(cell);
    } else if (!y.contains(cell)) {
      for (int i = 0; i < 4; i++) {
        getDifferenceInternal(cell.child(i), y);
      }
    }
  }

  /**
   * Just as normal binary search, except that it allows specifying the starting value for the lower
   * bound.
   *
   * @return The position of the searched element in the list (if found), or the position where the
   *     element could be inserted without violating the order.
   */
  private static int indexedBinarySearch(List<S2CellId> l, S2CellId key, int low) {
    int high = l.size() - 1;

    while (low <= high) {
      int mid = (low + high) >> 1;
      S2CellId midVal = l.get(mid);
      int cmp = midVal.compareTo(key);

      if (cmp < 0) {
        low = mid + 1;
      } else if (cmp > 0) {
        high = mid - 1;
      } else {
        return mid; // key found
      }
    }
    return low; // key not found
  }

  /**
   * Expands the cell union such that it contains all cells of the given level that are adjacent to
   * any cell of the original union. Two cells are defined as adjacent if their boundaries have any
   * points in common, i.e. most cells have 8 adjacent cells (not counting the cell itself).
   *
   * <p>Note that the size of the output is exponential in "level". For example, if level == 20 and
   * the input has a cell at level 10, there will be on the order of 4000 adjacent cells in the
   * output. For most applications the Expand(min_fraction, min_distance) method below is easier to
   * use.
   */
  public void expand(int level) {
    ArrayList<S2CellId> output = new ArrayList<S2CellId>();
    long levelLsb = S2CellId.lowestOnBitForLevel(level);
    for (int i = size(); --i >= 0; ) {
      S2CellId id = cellId(i);
      if (id.lowestOnBit() < levelLsb) {
        id = id.parent(level);
        // Optimization: skip over any cells contained by this one. This is
        // especially important when very small regions are being expanded.
        while (i > 0 && id.contains(cellId(i - 1))) {
          --i;
        }
      }
      output.add(id);
      id.getAllNeighbors(level, output);
    }
    initSwap(output);
  }

  /**
   * Expand the cell union such that it contains all points whose distance to the cell union is at
   * most minRadius, but do not use cells that are more than maxLevelDiff levels higher than the
   * largest cell in the input. The second parameter controls the tradeoff between accuracy and
   * output size when a large region is being expanded by a small amount (e.g. expanding Canada by
   * 1km).
   *
   * <p>For example, if maxLevelDiff == 4, the region will always be expanded by approximately 1/16
   * the width of its largest cell. Note that in the worst case, the number of cells in the output
   * can be up to 4 * (1 + 2 ** maxLevelDiff) times larger than the number of cells in the input.
   */
  public void expand(S1Angle minRadius, int maxLevelDiff) {
    int minLevel = S2CellId.MAX_LEVEL;
    for (S2CellId id : this) {
      minLevel = Math.min(minLevel, id.level());
    }
    // Find the maximum level such that all cells are at least "min_radius"
    // wide.
    int radiusLevel = PROJ.minWidth.getMaxLevel(minRadius.radians());
    if (radiusLevel == 0 && minRadius.radians() > PROJ.minWidth.getValue(0)) {
      // The requested expansion is greater than the width of a face cell.
      // The easiest way to handle this is to expand twice.
      expand(0);
    }
    expand(Math.min(minLevel + maxLevelDiff, radiusLevel));
  }

  // NOTE: This should be marked as @Override, but clone() isn't present in GWT's version of
  // Object, so we can't mark it as such.
  @SuppressWarnings("MissingOverride")
  public S2Region clone() {
    S2CellUnion copy = new S2CellUnion();
    copy.initRawCellIds(Lists.newArrayList(cellIds));
    return copy;
  }

  @Override
  public S2Cap getCapBound() {
    // Compute the approximate centroid of the region. This won't produce the
    // bounding cap of minimal area, but it should be close enough.
    if (cellIds.isEmpty()) {
      return S2Cap.empty();
    }
    S2Point centroid = S2Point.ORIGIN;
    for (S2CellId id : this) {
      double area = S2Cell.averageArea(id.level());
      centroid = S2Point.add(centroid, S2Point.mul(id.toPoint(), area));
    }
    if (centroid.equalsPoint(S2Point.ORIGIN)) {
      centroid = S2Point.X_POS;
    } else {
      centroid = S2Point.normalize(centroid);
    }

    // Use the centroid as the cap axis, and expand the cap angle so that it
    // contains the bounding caps of all the individual cells. Note that it is
    // *not* sufficient to just bound all the cell vertices because the bounding
    // cap may be concave (i.e. cover more than one hemisphere).
    S2Cap cap = S2Cap.fromAxisChord(centroid, S1ChordAngle.ZERO);
    for (S2CellId id : this) {
      cap = cap.addCap(new S2Cell(id).getCapBound());
    }
    return cap;
  }

  @Override
  public S2LatLngRect getRectBound() {
    S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();
    for (S2CellId id : this) {
      builder.union(new S2Cell(id).getRectBound());
    }
    return builder.build();
  }

  /** This is a fast operation (logarithmic in the size of the cell union). */
  @Override
  public boolean mayIntersect(S2Cell cell) {
    return intersects(cell.id());
  }

  /**
   * The point 'p' does not need to be normalized. This is a fast operation (logarithmic in the size
   * of the cell union).
   */
  @Override
  public boolean contains(S2Point p) {
    return contains(S2CellId.fromPoint(p));
  }

  /**
   * The number of leaf cells covered by the union. This will be no more than 6*2^60 for the whole
   * sphere.
   *
   * @return the number of leaf cells covered by the union
   */
  public long leafCellsCovered() {
    long numLeaves = 0;
    for (S2CellId cellId : cellIds) {
      int invertedLevel = S2CellId.MAX_LEVEL - cellId.level();
      numLeaves += (1L << (invertedLevel << 1));
    }
    return numLeaves;
  }

  /**
   * Approximate this cell union's area by summing the average area of each contained cell's average
   * area, using {@link S2Cell#averageArea()}. This is equivalent to the number of leaves covered,
   * multiplied by the average area of a leaf.
   *
   * <p>Note that {@link S2Cell#averageArea()} does not take into account distortion of cell, and
   * thus may be off by up to a factor of 1.7. NOTE: Since this is proportional to
   * LeafCellsCovered(), it is always better to use the other function if all you care about is the
   * relative average area between objects.
   *
   * @return the sum of the average area of each contained cell's average area
   */
  public double averageBasedArea() {
    return S2Cell.averageArea(S2CellId.MAX_LEVEL) * leafCellsCovered();
  }

  /**
   * Calculates this cell union's area by summing the approximate area for each contained cell,
   * using {@link S2Cell#approxArea()}.
   *
   * @return approximate area of the cell union
   */
  public double approxArea() {
    double area = 0;
    for (S2CellId cellId : cellIds) {
      area += new S2Cell(cellId).approxArea();
    }
    return area;
  }

  /**
   * Calculates this cell union's area by summing the exact area for each contained cell, using the
   * {@link S2Cell#exactArea()}.
   *
   * @return the exact area of the cell union
   */
  public double exactArea() {
    double area = 0;
    for (S2CellId cellId : cellIds) {
      area += new S2Cell(cellId).exactArea();
    }
    return area;
  }

  /** Return true if two cell unions are identical. */
  @Override
  public boolean equals(Object that) {
    if (!(that instanceof S2CellUnion)) {
      return false;
    }
    S2CellUnion union = (S2CellUnion) that;
    return this.cellIds.equals(union.cellIds);
  }

  @Override
  public int hashCode() {
    int value = 17;
    for (S2CellId id : this) {
      value = 37 * value + id.hashCode();
    }
    return value;
  }

  /**
   * Normalizes the cell union by discarding cells that are contained by other cells, replacing
   * groups of 4 child cells by their parent cell whenever possible, and sorting all the cell ids in
   * increasing order. Returns true if the number of cells was reduced.
   *
   * <p>This method *must* be called before doing any calculations on the cell union, such as
   * Intersects() or Contains().
   *
   * @return true if the normalize operation had any effect on the cell union, false if the union
   *     was already normalized
   */
  public boolean normalize() {
    return normalize(cellIds);
  }

  /** Like {@link #normalize()}, but works directly with a vector of S2CellIds. */
  public static boolean normalize(List<S2CellId> ids) {
    // Optimize the representation by looking for cases where all subcells of a parent cell are
    // present.
    Collections.sort(ids);
    int out = 0;
    for (int i = 0; i < ids.size(); i++) {
      S2CellId id = ids.get(i);

      // Check whether this cell is contained by the previous cell.
      if (out > 0 && ids.get(out - 1).contains(id)) {
        continue;
      }

      // Discard any previous cells contained by this cell.
      while (out > 0 && id.contains(ids.get(out - 1))) {
        out--;
      }

      // Check whether the last 3 elements of "output" plus "id" can be collapsed into a single
      // parent cell.
      while (out >= 3) {
        // A necessary (but not sufficient) condition is that the XOR of the
        // four cells must be zero. This is also very fast to test.
        if ((ids.get(out - 3).id() ^ ids.get(out - 2).id() ^ ids.get(out - 1).id()) != id.id()) {
          break;
        }

        // Now we do a slightly more expensive but exact test. First, compute a
        // mask that blocks out the two bits that encode the child position of
        // "id" with respect to its parent, then check that the other three
        // children all agree with "mask.
        long mask = id.lowestOnBit() << 1;
        mask = ~(mask + (mask << 1));
        long idMasked = (id.id() & mask);
        if ((ids.get(out - 3).id() & mask) != idMasked
            || (ids.get(out - 2).id() & mask) != idMasked
            || (ids.get(out - 1).id() & mask) != idMasked
            || id.isFace()) {
          break;
        }

        // Replace four children by their parent cell.
        id = id.parent();
        out -= 3;
      }

      ids.set(out++, id);
    }

    int size = ids.size();
    boolean trimmed = out < size;
    while (out < size) {
      size--;
      ids.remove(size);
    }
    return trimmed;
  }

  /**
   * Writes a simple lossless encoding of this cell union to the given output stream. This encoding
   * uses 1 byte for a version number, and N+1 64-bit longs where the first is the number of longs
   * that follow.
   *
   * @throws IOException there is a problem writing to the underlying stream
   */
  public void encode(OutputStream output) throws IOException {
    encode(new LittleEndianOutput(output));
  }

  /**
   * As {@link #encode(OutputStream)}, but avoids creating a little endian output helper.
   *
   * <p>Use this method if a number of S2 objects will be encoded to the same underlying stream.
   */
  public void encode(LittleEndianOutput output) throws IOException {
    output.writeByte(LOSSLESS_ENCODING_VERSION);
    output.writeLong(cellIds.size());
    for (S2CellId cellId : this) {
      output.writeLong(cellId.id());
    }
  }

  /**
   * Decodes an S2CellUnion encoded with Encode(). Returns true on success.
   *
   * @throws IOException there is a problem reading from the underlying stream, the version number
   *     doesn't match, or the number of elements to read is not between 0 and 2^31-1.
   */
  public static S2CellUnion decode(InputStream input) throws IOException {
    return decode(new LittleEndianInput(input));
  }

  /**
   * As {@link #decode(InputStream)}, but avoids creating a little endian input helper.
   *
   * <p>Use this method if a number of S2 objects will be decoded from the same underlying stream.
   */
  public static S2CellUnion decode(LittleEndianInput input) throws IOException {
    // Should contain at least version and vector length.
    byte version = input.readByte();
    if (version != LOSSLESS_ENCODING_VERSION) {
      throw new IOException("Unrecognized version number " + version);
    }
    long numCells = input.readLong();
    if (numCells < 0 || numCells > Integer.MAX_VALUE) {
      throw new IOException("Unsupported number of cells encountered: " + numCells);
    }
    S2CellUnion result = new S2CellUnion();
    for (int i = 0; i < numCells; i++) {
      result.cellIds().add(new S2CellId(input.readLong()));
    }
    return result;
  }
}
