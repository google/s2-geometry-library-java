/*
 * Copyright 2014 Google Inc.
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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2Projections.faceSiTiToXyz;
import static com.google.common.geometry.S2Projections.siTiToSt;
import static com.google.common.geometry.S2Projections.stToIj;
import static com.google.common.geometry.S2Projections.stToUV;
import static com.google.common.geometry.S2Projections.uvToST;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.geometry.R1Interval.Endpoint;

/**
 * S2PaddedCell represents an S2Cell whose (u,v)-range has been expanded on all sides by a given
 * amount of "padding". Unlike S2Cell, its methods and representation are optimized for clipping
 * edges against S2Cell boundaries to determine which cells are intersected by a given set of edges.
 */
public class S2PaddedCell {
  /** The cell being padded. */
  private final S2CellId id;

  /** UV padding on all sides. */
  private final double padding;

  /** Bound in (u,v)-space. Includes padding. */
  private final R2Rect bound;

  /**
   * The rectangle in (u,v)-space that belongs to all four padded children. It is computed on demand
   * by the middle() accessor method.
   */
  private R2Rect middle;

  /** Minimum (i,j)-coordinates of this cell, before padding. */
  private final int iLo;

  private final int jLo;

  /** Hilbert curve orientation of this cell. */
  private final int orientation;

  /** Level of this cell. */
  private final int level;

  /** Construct an S2PaddedCell for the given cell id and padding. */
  public S2PaddedCell(S2CellId id, double padding) {
    this.id = id;
    this.padding = padding;
    if (id.isFace()) {
      // Fast path for constructing a top-level face (the most common case).
      double limit = 1 + padding;
      bound = new R2Rect(new R1Interval(-limit, limit), new R1Interval(-limit, limit));
      middle = new R2Rect(new R1Interval(-padding, padding), new R1Interval(-padding, padding));
      iLo = jLo = 0;
      orientation = id.face() & 1;
      level = 0;
    } else {
      long ijo = id.toIJOrientation();
      int i = S2CellId.getI(ijo);
      int j = S2CellId.getJ(ijo);
      orientation = S2CellId.getOrientation(ijo);
      level = id.level();
      bound = S2CellId.ijLevelToBoundUv(i, j, level).expanded(padding);
      int ijSize = S2CellId.getSizeIJ(level);
      iLo = i & -ijSize;
      jLo = j & -ijSize;
    }
  }

  /**
   * Construct the child of this cell with the given (i,j) index. The four child cells have indices
   * of (0,0), (0,1), (1,0), (1,1), where the i and j indices correspond to increasing u- and
   * v-values respectively.
   */
  public S2PaddedCell childAtIJ(int i, int j) {
    return new S2PaddedCell(this, S2.ijToPos(orientation, i * 2 + j), i, j);
  }

  /** Construct the child of this cell with the given Hilbert curve position, from 0 to 3. */
  public S2PaddedCell childAtPos(int pos) {
    int ij = S2.posToIJ(orientation, pos);
    return new S2PaddedCell(this, pos, ij >> 1, ij & 1);
  }

  /** Private constructor to create a new S2PaddedCell for the child at the given (i,j) position. */
  private S2PaddedCell(S2PaddedCell parent, int pos, int i, int j) {
    this.padding = parent.padding;
    this.bound = new R2Rect(parent.bound);
    this.level = parent.level + 1;

    // Compute the position and orientation of the child incrementally from the orientation of the
    // parent.
    id = parent.id.child(pos);
    int ijSize = S2CellId.getSizeIJ(level);
    iLo = parent.iLo + i * ijSize;
    jLo = parent.jLo + j * ijSize;
    orientation = parent.orientation ^ S2.posToOrientation(pos);

    // For each child, one corner of the bound is taken directly from the parent while the
    // diagonally opposite corner is taken from middle().
    R2Rect middle = parent.middle();
    Endpoint uEnd = i == 0 ? Endpoint.HI : Endpoint.LO;
    bound.x().setValue(uEnd, middle.x().getValue(uEnd));
    Endpoint vEnd = j == 0 ? Endpoint.HI : Endpoint.LO;
    bound.y().setValue(vEnd, middle.y().getValue(vEnd));
  }

  /** Returns the ID of this padded cell. */
  public S2CellId id() {
    return id;
  }

  /** Returns the padding around this cell. */
  public double padding() {
    return padding;
  }

  /** Returns the level of this cell. */
  public int level() {
    return level;
  }

  /** Returns the orientation of this cell. */
  public int orientation() {
    return orientation;
  }

  /** Returns the bound for this cell (including padding.) */
  public R2Rect bound() {
    return bound;
  }

  /**
   * Return the "middle" of the padded cell, defined as the rectangle that belongs to all four
   * children.
   *
   * <p>Note that this method is *not* thread-safe, because the return value is computed on demand
   * and cached. (It is expected that this class will be mainly useful in the context of single-
   * threaded recursive algorithms.)
   */
  public R2Rect middle() {
    // We compute this field lazily because it is not needed the majority of the time (i.e., for
    // cells where the recursion terminates.)
    if (middle == null) {
      int ijSize = S2CellId.getSizeIJ(level);
      double u = stToUV(siTiToSt(2L * iLo + ijSize));
      double v = stToUV(siTiToSt(2L * jLo + ijSize));
      middle =
          new R2Rect(
              new R1Interval(u - padding, u + padding), new R1Interval(v - padding, v + padding));
    }
    return middle;
  }

  /**
   * Returns the smallest cell that contains all descendants of this cell whose bounds intersect
   * "rect". For algorithms that use recursive subdivision to find the cells that intersect a
   * particular object, this method can be used to skip all the initial subdivision steps where only
   * one child needs to be expanded.
   *
   * <p>Note that this method is not the same as returning the smallest cell that contains the
   * intersection of this cell with "rect". Because of the padding, even if one child completely
   * contains "rect" it is still possible that a neighboring child also intersects "rect".
   *
   * <p>Results are undefined if {@link #bound()} does not intersect the given rectangle.
   */
  public S2CellId shrinkToFit(R2Rect rect) {
    // assert bound().intersects(rect);

    // Quick rejection test: if "rect" contains the center of this cell along either axis, then no
    // further shrinking is possible.
    int ijSize = S2CellId.getSizeIJ(level);
    if (level == 0) {
      // Fast path (most calls to this function start with a face cell).
      if (rect.x().contains(0) || rect.y().contains(0)) {
        return id();
      }
    } else {
      if (rect.x().contains(stToUV(siTiToSt(2L * iLo + ijSize)))
          || rect.y().contains(stToUV(siTiToSt(2L * jLo + ijSize)))) {
        return id();
      }
    }
    // Otherwise we expand "rect" by the given padding() on all sides and find the range of
    // coordinates that it spans along the i- and j-axes. We then compute the highest bit position
    // at which the min and max coordinates differ. This corresponds to the first cell level at
    // which at least two children intersect "rect".

    // Increase the padding to compensate for the error in uvToST().
    // (The constant below is a provable upper bound on the additional error.)
    R2Rect padded = rect.expanded(padding() + 1.5 * DBL_EPSILON);
    int iMin = max(iLo, stToIj(uvToST(padded.x().lo())));
    int jMin = max(jLo, stToIj(uvToST(padded.y().lo())));
    int iMax = min(iLo + ijSize - 1, stToIj(uvToST(padded.x().hi())));
    int jMax = min(jLo + ijSize - 1, stToIj(uvToST(padded.y().hi())));
    int iXor = iMin ^ iMax;
    int jXor = jMin ^ jMax;

    // Compute the highest bit position where the two i- or j-endpoints differ, and then choose the
    // cell level that includes both of these endpoints. So if both pairs of endpoints are equal we
    // choose MAX_LEVEL; if they differ only at bit 0, we choose (MAX_LEVEL - 1), and so on.
    int levelMsb = ((iXor | jXor) << 1) + 1;
    int level = S2CellId.MAX_LEVEL - floorLog2(levelMsb);
    if (level <= this.level) {
      return id();
    }
    return S2CellId.fromFaceIJ(id().face(), iMin, jMin).parent(level);
  }

  /** Returns the floor of the log2 of x, assuming x is positive. */
  private static final int floorLog2(long x) {
    return 63 - Long.numberOfLeadingZeros(x);
  }

  /** Returns the center of this cell. */
  public S2Point getCenter() {
    int ijSize = S2CellId.getSizeIJ(level);
    long si = 2L * iLo + ijSize;
    long ti = 2L * jLo + ijSize;
    return faceSiTiToXyz(id.face(), si, ti).normalize();
  }

  /** Returns the vertex where the S2 space-filling curve enters this cell. */
  public S2Point getEntryVertex() {
    // The curve enters at the (0,0) vertex unless the axis directions are reversed, in which case
    // it enters at the (1,1) vertex.
    int i = iLo;
    int j = jLo;
    if ((orientation & S2.INVERT_MASK) != 0) {
      int ijSize = S2CellId.getSizeIJ(level);
      i += ijSize;
      j += ijSize;
    }
    return faceSiTiToXyz(id.face(), 2L * i, 2L * j).normalize();
  }

  /** Returns the vertex where the S2 space-filling curve exits this cell. */
  public S2Point getExitVertex() {
    // The curve exits at the (1,0) vertex unless the axes are swapped or inverted but not both, in
    // which case it exits at the (0,1) vertex.
    int i = iLo;
    int j = jLo;
    int ijSize = S2CellId.getSizeIJ(level);
    if (orientation == 0 || orientation == S2.SWAP_MASK + S2.INVERT_MASK) {
      i += ijSize;
    } else {
      j += ijSize;
    }
    return faceSiTiToXyz(id.face(), 2L * i, 2L * j).normalize();
  }
}
