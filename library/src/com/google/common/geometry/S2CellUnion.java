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

import static com.google.common.geometry.S2CellId.isFace;
import static com.google.common.geometry.S2CellId.parentAsLong;
import static com.google.common.geometry.S2CellId.rangeMaxAsLong;
import static com.google.common.geometry.S2CellId.rangeMinAsLong;
import static com.google.common.geometry.S2CellId.unsignedLongLessThan;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

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
 * @author danieldanciu (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
@SuppressWarnings("Assertion")
public class S2CellUnion implements S2Region, Iterable<S2CellId>, Serializable {
  private static final long serialVersionUID = 1L;

  private static final byte LOSSLESS_ENCODING_VERSION = 1;

  /** An {@link S2Coder} of cell unions that uses {@link #encode} and {@link #decode}. */
  public static final S2Coder<S2CellUnion> FAST_CODER =
      new S2Coder<S2CellUnion>() {
        @JsIgnore // OutputStream is not available to J2CL.
        @Override
        public void encode(S2CellUnion value, OutputStream output) throws IOException {
          value.encode(output);
        }

        @Override
        public S2CellUnion decode(Bytes data, Cursor cursor) throws IOException {
          return S2CellUnion.decode(data.toInputStream(cursor));
        }

        @Override
        public boolean isLazy() {
          return false;
        }
      };

  /** A compact coder that compresses the given cells by around 4-5x in many cases. */
  public static final S2Coder<S2CellUnion> COMPACT_CODER =
      S2CellIdVectorCoder.INSTANCE.delegating(cells -> cells.cellIds, S2CellUnion::copyFrom);

  /** The CellIds that form the Union */
  private ArrayList<S2CellId> cellIds = new ArrayList<>();

  @JsConstructor
  public S2CellUnion() {}

  /** Clears the union contents, leaving it empty. */
  public void clear() {
    cellIds.clear();
  }

  /** Creates a new cell union as as copy of the given cell union 'other'. */
  public static S2CellUnion copyFrom(S2CellUnion other) {
    S2CellUnion copy = new S2CellUnion();
    copy.initRawCellIds(new ArrayList<>(other.cellIds));
    return copy;
  }

  /** Creates a new cell union from a copy of the given cells. */
  @JsIgnore // J2CL warning "Iterable<S2CellId> ... is not usable by JavaScript" but not clear why.
  public static S2CellUnion copyFrom(Iterable<S2CellId> cells) {
    S2CellUnion result = new S2CellUnion();
    // Note that if 'cells' are an AbstractList over lazily decoded data, addAll may throw
    // NoSuchElementException.
    Iterables.addAll(result.cellIds, cells);
    return result;
  }

  /** Constructs a cell union for the whole sphere. */
  public static S2CellUnion wholeSphere() {
    return S2CellUnion.copyFrom(Arrays.asList(S2CellId.FACE_CELLS));
  }

  /**
   * Populates this cell union with the given S2CellIds, and then calls normalize(). This directly
   * uses the input list, without copying it.
   */
  @CanIgnoreReturnValue
  public S2CellUnion initFromCellIds(ArrayList<S2CellId> cellIds) {
    initRawCellIds(cellIds);
    normalize();
    return this;
  }

  /** Populates this cell union with the given 64-bit cell ids, and then calls normalize(). */
  @CanIgnoreReturnValue
  public S2CellUnion initFromIds(List<Long> cellIds) {
    initRawIds(cellIds);
    normalize();
    return this;
  }

  /**
   * Populates this cell union with the given S2CellIds. The input list is copied, and then cleared.
   */
  @CanIgnoreReturnValue
  public S2CellUnion initSwap(List<S2CellId> cellIds) {
    initRawSwap(cellIds);
    normalize();
    return this;
  }

  /**
   * Populates this cell union with the given S2CellIds. This does not call normalize, see {@link
   * #initRawSwap} for details. This directly uses the input list, without copying it.
   */
  @CanIgnoreReturnValue
  public S2CellUnion initRawCellIds(ArrayList<S2CellId> cellIds) {
    this.cellIds = cellIds;
    return this;
  }

  /** Populates this cell union with the single given 64 bit cell id, which must be valid. */
  @CanIgnoreReturnValue
  public S2CellUnion initFromId(long cellId) {
    assert S2CellId.isValid(cellId);
    cellIds.clear();
    cellIds.add(new S2CellId(cellId));
    return this;
  }

  /** Populates this cell union with the single given S2CellId, which must be valid. */
  @CanIgnoreReturnValue
  public S2CellUnion initFromCellId(S2CellId cellId) {
    assert cellId.isValid();
    cellIds.clear();
    cellIds.add(cellId);
    return this;
  }

  /**
   * Populates this cell union with the given 64 bit cell ids. This does not call normalize, see
   * {@link #initRawSwap} for details. The input list is copied.
   */
  // TODO(user): Make a constructed S2CellUnion immutable, and port other init methods from
  // C++.
  @CanIgnoreReturnValue
  public S2CellUnion initRawIds(List<Long> cellIds) {
    int size = cellIds.size();
    this.cellIds = new ArrayList<>(size);
    for (Long id : cellIds) {
      this.cellIds.add(new S2CellId(id));
    }
    return this;
  }

  /**
   * Like the initFrom*() constructors, but does not call normalize(). The cell union *must* be
   * normalized before doing any calculations with it, so it is the caller's responsibility to make
   * sure that the input is normalized. This method is useful when converting cell unions to another
   * representation and back.
   *
   * <p>The input list is copied, and then cleared.
   */
  public void initRawSwap(List<S2CellId> cellIds) {
    this.cellIds = new ArrayList<>(cellIds);
    cellIds.clear();
  }

  /**
   * Create a cell union that corresponds to a continuous range of cell ids. The output is a
   * normalized collection of cell ids that covers the leaf cells between "minId" and "maxId"
   * inclusive.
   *
   * <p>Requires that {@code minId.isLeaf(), maxId.isLeaf()}, and {@code minId <= maxId}.
   */
  public void initFromMinMax(S2CellId minId, S2CellId maxId) {
    assert minId.isLeaf();
    assert maxId.isLeaf();
    assert minId.compareTo(maxId) <= 0;
    assert minId.isValid() && maxId.isValid();
    initFromBeginEnd(minId, maxId.next());
  }

  /**
   * As {@link #initFromMinMax(S2CellId, S2CellId)}, except that the union covers the range of leaf
   * cells from "begin" (inclusive) to "end" (exclusive.) If {@code begin.equals(end)}, the result
   * is empty.
   *
   * <p>Requires that {@code begin.isLeaf(), end.isLeaf()}, and {@code begin <= end}.
   */
  public void initFromBeginEnd(S2CellId begin, S2CellId end) {
    assert begin.isLeaf();
    assert end.isLeaf();
    assert begin.compareTo(end) <= 0;

    // We repeatedly add the largest cell we can, in sorted order.
    cellIds.clear();
    for (S2CellId nextBegin = begin; nextBegin.compareTo(end) < 0; ) {
      assert nextBegin.isLeaf();

      // Find the largest cell that starts at "nextBegin" and ends before "end". This loop uses
      // longs rather than S2CellIds to avoid many allocations of S2CellId objects.
      long nextId = nextBegin.id();
      while (!isFace(nextId)
          && rangeMinAsLong(parentAsLong(nextId)) == nextBegin.id()
          && unsignedLongLessThan(rangeMaxAsLong(parentAsLong(nextId)), end.id())) {
        nextId = parentAsLong(nextId);
      }
      S2CellId nextCellId = new S2CellId(nextId);
      cellIds.add(nextCellId);
      nextBegin = nextCellId.rangeMax().next();
    }

    // The output should already be sorted and normalized.
    assert !normalize();
  }

  public int size() {
    return cellIds.size();
  }

  /** Convenience methods for accessing the individual cell ids. */
  public S2CellId cellId(int i) {
    return cellIds.get(i);
  }

  /**
   * Provides an S2Iterator over this union's cells. The cell union must be normalized, or results
   * will be undefined.
   */
  public S2Iterator<S2CellId> s2Iterator() {
    assert isNormalized();
    return S2Iterator.fromList(cellIds);
  }

  /** Enable iteration over the union's cells. See also {@link #s2Iterator()}. */
  @Override
  @JsIgnore
  public Iterator<S2CellId> iterator() {
    return cellIds.iterator();
  }

  /** Direct access to the underlying vector for iteration . */
  public ArrayList<S2CellId> cellIds() {
    return cellIds;
  }

  /** Returns true if the cell union is empty. */
  public boolean isEmpty() {
    return cellIds.isEmpty();
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
   * that no four cells at the same level have a common parent.
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
   * Returns true if the given four cells are at the same level and have a common parent.
   *
   * <p>Requires the four cells are distinct.
   */
  private static boolean areSiblings(S2CellId a, S2CellId b, S2CellId c, S2CellId d) {
    // A necessary (but not sufficient) condition is that the XOR of the four cells must be zero.
    // This is also very fast to test.
    if ((a.id() ^ b.id() ^ c.id()) != d.id()) {
      return false;
    }

    // Now we do a slightly more expensive but exact test. First, compute a mask that blocks out the
    // two bits that encode the child position of "id" with respect to its parent, then check that
    // the other three children all agree with "mask".
    long mask = d.lowestOnBit() << 1;
    mask = ~(mask + (mask << 1));
    long idMasked = d.id() & mask;
    return !d.isFace()
        && (a.id() & mask) == idMasked
        && (b.id() & mask) == idMasked
        && (c.id() & mask) == idMasked;
  }

  /**
   * Returns a list of the cell ids in this cell union after replacing any cells whose level is less
   * than "minLevel" with their children, until the required minLevel is reached. The provided
   * minLevel must be in the range [0, S2CellId.MAX_LEVEL].
   */
  public List<S2CellId> denormalized(int minLevel) {
    ArrayList<S2CellId> output = new ArrayList<>();
    denormalize(minLevel, 1, output);
    return output;
  }

  /**
   * Replaces "output" with an expanded version of the cell union where any cells whose level is
   * less than "minLevel" are replaced by their children, until the required minLevel is reached.
   *
   * <p>The provided minLevel must be in the range [0, S2CellId.MAX_LEVEL].
   */
  @JsIgnore
  public void denormalize(int minLevel, ArrayList<S2CellId> output) {
    denormalize(minLevel, 1, output);
  }

  /**
   * Replaces "output" with an expanded version of the cell union where any cells whose level is
   * less than "minLevel" or where (level - minLevel) is not a multiple of "levelMod" are replaced
   * by their children, until either both of these conditions are satisfied or the maximum level is
   * reached.
   *
   * <p>This method allows a covering generated by S2RegionCoverer using minLevel() or levelMod()
   * constraints to be stored as a normalized cell union (which allows various geometric
   * computations to be done) and then converted back to the original list of cell ids that
   * satisfies the desired constraints.
   *
   * <p>The provided minLevel must be in the range [0, S2CellId.MAX_LEVEL]. The provided levelMod
   * must be in the range [1, 3].
   */
  public void denormalize(int minLevel, int levelMod, ArrayList<S2CellId> output) {
    assert minLevel >= 0;
    assert minLevel <= S2CellId.MAX_LEVEL;
    assert levelMod >= 1;
    assert levelMod <= 3;

    output.clear();
    output.ensureCapacity(size());
    for (S2CellId id : this) {
      int level = id.level();
      int newLevel = max(minLevel, level);
      if (levelMod > 1) {
        // Round up so that (newLevel - minLevel) is a multiple of levelMod.
        // (Note that S2CellId.MAX_LEVEL is a multiple of 1, 2, and 3.)
        newLevel += (S2CellId.MAX_LEVEL - (newLevel - minLevel)) % levelMod;
        newLevel = min(S2CellId.MAX_LEVEL, newLevel);
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
   * If there are more than "excess" elements of the cellIds() vector that are allocated but unused,
   * reallocate the array to eliminate the excess space. This reduces memory usage when many cell
   * unions need to be held in memory at once.
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
  @JsMethod(name = "containsCellId")
  public boolean contains(S2CellId id) {
    // This is an exact test. Each cell occupies a linear span of the S2 space-filling curve, and
    // the cell id is simply the position at the center of this span. The cell union ids are sorted
    // in increasing order along the space-filling curve. So we simply find the pair of cell ids
    // that surround the given cell id (using binary search). There is containment if and only if
    // one of these two cell ids contains this cell.
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
  @JsMethod(name = "intersectsCellId")
  public boolean intersects(S2CellId id) {
    // This is an exact test; see the comments for contains(S2CellId) above.
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
  @JsMethod(name = "containsCell")
  public boolean contains(S2Cell cell) {
    return contains(cell.id());
  }

  /** Return true if this cell union intersects {@code union}. */
  public boolean intersects(S2CellUnion union) {
    S2CellUnion result = new S2CellUnion();
    result.getIntersection(this, union);
    return result.size() > 0;
  }

  /** Returns the union of two S2CellUnions. */
  public static S2CellUnion union(S2CellUnion x, S2CellUnion y) {
    S2CellUnion result = new S2CellUnion();
    result.getUnion(x, y);
    return result;
  }

  /**
   * Sets this cell union to the union of {@code x} and {@code y}, which both must be different cell
   * unions than this one.
   */
  @SuppressWarnings("ReferenceEquality") // Precondition check is checking reference equality.
  public void getUnion(S2CellUnion x, S2CellUnion y) {
    Preconditions.checkArgument(x != this);
    Preconditions.checkArgument(y != this);
    cellIds.clear();
    cellIds.ensureCapacity(x.size() + y.size());
    cellIds.addAll(x.cellIds);
    cellIds.addAll(y.cellIds);
    normalize();
  }

  /**
   * Specialized version of getIntersection() that gets the intersection of a cell union with the
   * given cell id. This can be useful for "splitting" a cell union into chunks.
   *
   * <p><b>Note:</b> {@code x} must be normalized, and must be a different cell union than this one.
   */
  @SuppressWarnings("ReferenceEquality") // Precondition check is checking reference equality.
  public void getIntersection(S2CellUnion x, S2CellId id) {
    Preconditions.checkArgument(x != this);
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
    assert isNormalized() || !x.isNormalized();
  }

  /** Returns the intersection of two S2CellUnions. */
  public static S2CellUnion intersection(S2CellUnion x, S2CellUnion y) {
    S2CellUnion result = new S2CellUnion();
    result.getIntersection(x, y);
    return result;
  }

  /**
   * Initializes this cell union to the intersection of the two given cell unions. Requires: x !=
   * this and y != this.
   *
   * <p><b>Note:</b> {@code x} and {@code y} must both be normalized to ensure the output is
   * normalized.
   */
  @JsMethod(name = "getIntersectionCellUnion")
  @SuppressWarnings("ReferenceEquality") // Precondition check is checking reference equality.
  public void getIntersection(S2CellUnion x, S2CellUnion y) {
    // It's fine if 'this' and 'x' or 'y' are different cell unions with the same cells, but they
    // may not be the same object.
    Preconditions.checkArgument(x != this && y != this);
    getIntersection(x.cellIds, y.cellIds, cellIds);
    // The output is normalized as long as both inputs are normalized.
    assert isNormalized() || (!x.isNormalized() || !y.isNormalized());
  }

  /**
   * Like {@code #getIntersection(S2CellUnion, S2CellUnion)}, but works directly with lists of
   * S2CellIds, and this method has slightly more relaxed normalization requirements: the input
   * vectors may contain groups of 4 child cells that all have the same parent. (In a normalized
   * S2CellUnion, such groups are always replaced by the parent cell.)
   *
   * <p><b>Note:</b> {@code x} and {@code y} must be sorted.
   */
  @JsIgnore
  public static void getIntersection(List<S2CellId> x, List<S2CellId> y, List<S2CellId> results) {
    assert x != results;
    assert y != results;

    // This is a fairly efficient calculation that uses binary search to skip over sections of both
    // input vectors. It takes constant time if all the cells of "x" come before or after all the
    // cells of "y" in S2CellId order.
    results.clear();
    int i = 0;
    int j = 0;
    while (i < x.size() && j < y.size()) {
      S2CellId xCell = x.get(i);
      S2CellId xMin = xCell.rangeMin();
      S2CellId yCell = y.get(j);
      S2CellId yMin = yCell.rangeMin();
      if (xMin.greaterThan(yMin)) {
        // Either j.contains(xCell) or the two cells are disjoint.
        if (xCell.lessOrEquals(yCell.rangeMax())) {
          results.add(xCell);
          i++;
        } else {
          // Advance "j" to the first cell possibly contained by xCell.
          j = indexedBinarySearch(y, xMin, j + 1);
          // The previous cell (j-1) may now contain xCell.
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

  /** Initializes this cell union to the difference of the two given cell unions. */
  public void getDifference(S2CellUnion x, S2CellUnion y) {
    // TODO(user): this is approximately O(N*log(N)), but could probably use similar
    // techniques as getIntersection() to be more efficient.
    cellIds.clear();
    for (S2CellId id : x) {
      getDifferenceInternal(id, y);
    }
    // The output is normalized as long as the first argument is normalized.
    assert isNormalized() || !x.isNormalized();
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
   * Expands the cell union by adding a buffer of cells at "expandLevel" around the union boundary.
   *
   * <p>For each cell "c" in the union, we add all neighboring cells at level "expandLevel" that are
   * adjacent to "c". Note that there can be many such cells if "c" is large compared to
   * "expandLevel". If "c" is smaller than "expandLevel", we first add the parent of "c" at
   * "expandLevel", and then add all the neighbors of that cell. Note that this can cause the
   * expansion around such cells to be up to almost twice as large as expected, because the union
   * boundary is moved outward by up to the size difference between "c" and its parent at
   * "expandLevel", (depending on the position of "c" in that parent) *before* the buffering
   * neighbors are added.
   *
   * <p>Note that the size of the output is exponential in "level". For example, if level == 20 and
   * the input has a cell at level 10, there will be on the order of 4000 adjacent cells in the
   * output. For most applications the expand(minRadius, maxLevelDiff) method below is easier to
   * use.
   */
  @JsMethod(name = "expandAtLevel")
  public void expand(int expandLevel) {
    ArrayList<S2CellId> output = new ArrayList<>();
    long levelLsb = S2CellId.lowestOnBitForLevel(expandLevel);
    for (int i = size(); --i >= 0; ) {
      S2CellId id = cellId(i);
      if (id.lowestOnBit() < levelLsb) {
        id = id.parent(expandLevel);
        // Optimization: skip over any cells contained by this one. This is especially important
        // when very small regions are being expanded.
        while (i > 0 && id.contains(cellId(i - 1))) {
          --i;
        }
      }
      output.add(id);
      id.getAllNeighbors(expandLevel, output);
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
      minLevel = min(minLevel, id.level());
    }
    // Find the maximum level such that all cells are at least "minRadius" wide.
    int radiusLevel = S2Projections.MIN_WIDTH.getMaxLevel(minRadius.radians());
    if (radiusLevel == 0 && minRadius.radians() > S2Projections.MIN_WIDTH.getValue(0)) {
      // The requested expansion is greater than the width of a face cell. The easiest way to handle
      // this is to expand twice.
      expand(0);
    }
    expand(min(minLevel + maxLevelDiff, radiusLevel));
  }

  @Override
  public S2Cap getCapBound() {
    // Compute the approximate centroid of the region. This won't produce the bounding cap of
    // minimal area, but it should be close enough.
    if (cellIds.isEmpty()) {
      return S2Cap.empty();
    }
    S2Point centroid = S2Point.ZERO;
    for (S2CellId id : this) {
      double area = S2Cell.averageArea(id.level());
      centroid = centroid.add(id.toPoint().mul(area));
    }
    if (centroid.equalsPoint(S2Point.ZERO)) {
      centroid = S2Point.X_POS;
    } else {
      centroid = centroid.normalize();
    }

    // Use the centroid as the cap axis, and expand the cap angle so that it contains the bounding
    // caps of all the individual cells. Note that it is *not* sufficient to just bound all the
    // cell vertices because the bounding cap may be concave (i.e. cover more than one hemisphere).
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

  @Override
  public void getCellUnionBound(Collection<S2CellId> results) {
    results.clear();
    results.addAll(cellIds);
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
  @JsMethod(name = "containsPoint")
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

  @Override
  public String toString() {
    return cellIds.toString();
  }

  /**
   * Normalizes the cell union by discarding cells that are contained by other cells, replacing
   * groups of 4 child cells by their parent cell whenever possible, and sorting all the cell ids in
   * increasing order. Returns true if the number of cells was reduced.
   *
   * <p>This method *must* be called before doing any calculations on the cell union, such as
   * intersects() or contains().
   */
  @CanIgnoreReturnValue
  public boolean normalize() {
    return normalize(cellIds);
  }

  /** Like {@link #normalize()}, but works directly with a vector of S2CellIds. */
  @JsIgnore
  @CanIgnoreReturnValue
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
        // A necessary (but not sufficient) condition is that the XOR of the four cells must be
        // zero. This is also very fast to test.
        if ((ids.get(out - 3).id() ^ ids.get(out - 2).id() ^ ids.get(out - 1).id()) != id.id()) {
          break;
        }

        // Now we do a slightly more expensive but exact test. First, compute a mask that blocks out
        // the two bits that encode the child position of "id" with respect to its parent, then
        // check that the other three children all agree with "mask.
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
  @JsIgnore
  public void encode(OutputStream output) throws IOException {
    output.write(LOSSLESS_ENCODING_VERSION);
    LittleEndianOutput.writeLong(output, cellIds.size());
    for (S2CellId cellId : this) {
      LittleEndianOutput.writeLong(output, cellId.id());
    }
  }

  /**
   * Decodes an S2CellUnion encoded with encode(). Returns true on success. Decodes all the cell ids
   * immediately, i.e. is not lazy.
   *
   * <p>Use this method if a number of S2 objects will be decoded from the same underlying stream.
   *
   * @throws IOException there is a problem reading from the underlying stream, the version number
   *     doesn't match, or the number of elements to read is not between 0 and 2^31-1.
   */
  @JsIgnore
  public static S2CellUnion decode(InputStream input) throws IOException {
    // Should contain at least version and vector length.
    byte version = (byte) input.read();
    if (version != LOSSLESS_ENCODING_VERSION) {
      throw new IOException("Unrecognized version number " + version);
    }
    long numCells = LittleEndianInput.readLong(input);
    if (numCells < 0 || numCells > Integer.MAX_VALUE) {
      throw new IOException("Unsupported number of cells encountered: " + numCells);
    }
    S2CellUnion result = new S2CellUnion();
    for (int i = 0; i < numCells; i++) {
      result.cellIds().add(new S2CellId(LittleEndianInput.readLong(input)));
    }
    return result;
  }
}
