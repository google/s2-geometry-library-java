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
import static com.google.common.primitives.Ints.min;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.base.Objects;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import javax.annotation.Nullable;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * An S2RegionCoverer is a class that allows arbitrary regions to be approximated as unions of cells
 * (S2CellUnion). This is useful for implementing various sorts of search and precomputation
 * operations.
 *
 * <p>Typical usage: {@code S2RegionCoverer coverer =
 * S2RegionCoverer.builder().setMaxCells(5).build(); S2Cap cap = S2Cap.fromAxisAngle(...);
 * S2CellUnion covering; coverer.getCovering(cap, covering);}
 *
 * <p>This yields a cell union of at most 5 cells that is guaranteed to cover the given cap (a
 * disc-shaped region on the sphere).
 *
 * <p>The approximation algorithm is not optimal but does a pretty good job in practice. The output
 * does not always use the maximum number of cells allowed, both because this would not always yield
 * a better approximation, and because maxCells() is a limit on how much work is done exploring the
 * possible covering as well as a limit on the final output size.
 *
 * <p>One can also generate interior coverings, which are sets of cells which are entirely contained
 * within a region. Interior coverings can be empty, even for non-empty regions, if there are no
 * cells that satisfy the provided constraints and are contained by the region. Note that for
 * performance reasons, it is wise to specify a maxLevel when computing interior coverings -
 * otherwise for regions with small or zero area, the algorithm may spend a lot of time subdividing
 * cells all the way to leaf level to try to find contained cells.
 *
 * @author danieldanciu (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
public final strictfp class S2RegionCoverer implements Serializable {

  /**
   * A S2RegionCoverer configured with the default options. The min level, max level, and level mod
   * are unrestricted, and maxCells is {@link Builder#DEFAULT_MAX_CELLS}. See {@link Builder} for
   * details.
   */
  public static final S2RegionCoverer DEFAULT = builder().build();

  private static final ImmutableList<S2Cell> FACE_CELLS;

  static {
    ImmutableList.Builder<S2Cell> builder = ImmutableList.builder();
    for (int face = 0; face < 6; ++face) {
      builder.add(S2Cell.fromFace(face));
    }
    FACE_CELLS = builder.build();
  }

  private final int minLevel;
  private final int maxLevel;
  private final int levelMod;
  private final int maxCells;

  static class Candidate {
    private S2Cell cell;
    private boolean isTerminal; // Cell should not be expanded further.
    private int numChildren; // Number of children that intersect the region.
    private Candidate[] children; // Actual size may be 0, 4, 16, or 64 elements.
  }

  static class QueueEntry {
    private final int id;
    private final Candidate candidate;

    public QueueEntry(int id, Candidate candidate) {
      this.id = id;
      this.candidate = candidate;
    }
  }

  /**
   * We define our own comparison function on QueueEntries in order to make the results
   * deterministic.
   */
  static class QueueEntriesComparator implements Comparator<QueueEntry> {
    @Override
    public int compare(S2RegionCoverer.QueueEntry x, S2RegionCoverer.QueueEntry y) {
      return x.id < y.id ? 1 : (x.id > y.id ? -1 : 0);
    }
  }

  /**
   * Returns a new Builder with default values, which can be used to construct an S2RegionCoverer
   * instance.
   */
  public static Builder builder() {
    return new Builder();
  }

  /**
   * Construct from a {@link Builder}. Users should construct with
   * S2RegionCoverer.builder().build(), or use the DEFAULT instance.
   */
  private S2RegionCoverer(Builder builder) {
    minLevel = builder.getMinLevel();
    maxLevel = builder.getMaxLevel();
    levelMod = builder.getLevelMod();
    maxCells = builder.getMaxCells();
  }

  /** A Build to construct a {@link S2RegionCoverer} with options. */
  @JsType
  public static final class Builder {
    /**
     * By default, the covering uses at most 8 cells at any level. This gives a reasonable tradeoff
     * between the number of cells used and the accuracy of the approximation (see table below).
     */
    private static final int DEFAULT_MAX_CELLS = 8;

    private int minLevel = 0;
    private int maxLevel = S2CellId.MAX_LEVEL;
    private int levelMod = 1;
    private int maxCells = DEFAULT_MAX_CELLS;

    /** Users should create a Builder via the S2RegionCoverer.builder() method. */
    private Builder() {}

    // Set the minimum and maximum cell level to be used. The default is to use all cell levels.
    // Requires: maxLevel() >= minLevel().
    //
    // To find the cell level corresponding to a given physical distance, use the S2Cell metrics
    // defined in S2Projections.java. For example, to find the cell level that corresponds to an
    // average edge length of 10km, use:
    // {@code
    // import static com.google.common.geometry.S2Projections.PROJ;
    //
    // int level = PROJ.avgEdge.getClosestLevel(S2Earth.kmToRadians(10));
    // }
    // Note: minLevel() takes priority over maxCells(), i.e. cells below the given level will never
    // be used even if this causes a large number of cells to be returned.

    /**
     * Sets the minimum level to be used.
     *
     * <p>Default: 0
     */
    @CanIgnoreReturnValue
    public Builder setMinLevel(int minLevel) {
      // assert (minLevel >= 0 && minLevel <= S2CellId.MAX_LEVEL);
      this.minLevel = max(0, min(S2CellId.MAX_LEVEL, minLevel));
      return this;
    }

    /** Returns the minimum cell level to be used. */
    public int getMinLevel() {
      return minLevel;
    }

    /**
     * Sets the maximum level to be used.
     *
     * <p>Default: S2CellId.MAX_LEVEL
     */
    @CanIgnoreReturnValue
    public Builder setMaxLevel(int maxLevel) {
      // assert (maxLevel >= 0 && maxLevel <= S2CellId.MAX_LEVEL);
      this.maxLevel = max(0, min(S2CellId.MAX_LEVEL, maxLevel));
      return this;
    }

    /** Returns the maximum cell level to be used. */
    public int getMaxLevel() {
      return maxLevel;
    }

    /**
     * Only cells where (level - minLevel) is a multiple of "levelMod" will be used (default 1).
     * This effectively allows the branching factor of the S2CellId hierarchy to be increased.
     * Currently the only parameter values allowed are 1, 2, or 3, corresponding to branching
     * factors of 4, 16, and 64 respectively.
     *
     * <p>Default: 1
     */
    @CanIgnoreReturnValue
    public Builder setLevelMod(int levelMod) {
      // assert (levelMod >= 1 && levelMod <= 3);
      this.levelMod = max(1, min(3, levelMod));
      return this;
    }

    /** Returns the level mod. */
    public int getLevelMod() {
      return levelMod;
    }

    /**
     * Sets the maximum desired number of cells in the approximation (defaults to
     * DEFAULT_MAX_CELLS). Note the following:
     *
     * <ul>
     *   <li>For any setting of maxCells(), up to 6 cells may be returned if that is the minimum
     *       number of cells required (e.g. if the region intersects all six face cells). Up to 3
     *       cells may be returned even for very tiny convex regions if they happen to be located at
     *       the intersection of three cube faces.
     *   <li>For any setting of maxCells(), an arbitrary number of cells may be returned if
     *       minLevel() is too high for the region being approximated.
     *   <li>If maxCells() is less than 4, the area of the covering may be arbitrarily large
     *       compared to the area of the original region even if the region is convex (e.g. an S2Cap
     *       or S2LatLngRect).
     * </ul>
     *
     * <p>Accuracy is measured by dividing the area of the covering by the area of the original
     * region. The following table shows the median and worst case values for this area ratio on a
     * test case consisting of 100,000 spherical caps of random size (generated using
     * s2regioncoverer_unittest):
     *
     * <pre>
     * maxCells:          3      4      5      6      8     12     20    100   1000
     * median ratio:   5.33   3.32   2.73   2.34   1.98   1.66   1.42   1.11   1.01
     * worst case:   215518  14.41   9.72   5.26   3.91   2.75   1.92   1.20   1.02
     * </pre>
     *
     * <p>Default: 8
     */
    @CanIgnoreReturnValue
    public Builder setMaxCells(int maxCells) {
      this.maxCells = maxCells;
      return this;
    }

    /** Returns the maximum desired number of cells to be used. */
    public int getMaxCells() {
      return maxCells;
    }

    /** Constructs a {@link S2RegionCoverer} with this Builders options. */
    public S2RegionCoverer build() {
      return new S2RegionCoverer(this);
    }
  }

  @Override
  public boolean equals(Object obj) {
    if (obj instanceof S2RegionCoverer) {
      S2RegionCoverer that = (S2RegionCoverer) obj;
      return this.minLevel == that.minLevel
          && this.maxLevel == that.maxLevel
          && this.levelMod == that.levelMod
          && this.maxCells == that.maxCells;
    }
    return false;
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(minLevel, maxLevel, levelMod, maxCells);
  }

  public int minLevel() {
    return minLevel;
  }

  public int maxLevel() {
    return maxLevel;
  }

  public int maxCells() {
    return maxCells;
  }

  public int levelMod() {
    return levelMod;
  }

  /**
   * Computes a list of cell ids that covers the given region and satisfies the various restrictions
   * specified above.
   *
   * @param region The region to cover
   * @param covering The list filled in by this method
   */
  @JsMethod(name = "getCoveringList")
  public void getCovering(S2Region region, ArrayList<S2CellId> covering) {
    // Rather than just returning the raw list of cell ids generated by getCoveringInternal(), we
    // construct a cell union and then denormalize it. This has the effect of replacing four child
    // cells with their parent whenever this does not violate the covering parameters specified
    // (minLevel, levelMod, etc). This strategy significantly reduces the number of cells returned
    // in many cases, and it is cheap compared to computing the covering in the first place.
    S2CellUnion tmp = getCovering(region);
    tmp.denormalize(minLevel(), levelMod(), covering);
  }

  /**
   * Computes a list of cell ids that is contained within the given region and satisfies the various
   * restrictions specified above; note that if the max cell level is not specified very carefully
   * this method can try to create an enormous number of cells, wasting a lot of time and memory, so
   * care should be taken to set a max level suitable for the scale of the region being covered.
   *
   * @param region The region to fill
   * @param interior The list filled in by this method
   */
  @JsMethod(name = "getInteriorCoveringList")
  public void getInteriorCovering(S2Region region, ArrayList<S2CellId> interior) {
    S2CellUnion tmp = getInteriorCovering(region);
    tmp.denormalize(minLevel(), levelMod(), interior);
  }

  /**
   * Return a normalized cell union that covers the given region and satisfies the restrictions
   * *EXCEPT* for minLevel() and levelMod(). These criteria cannot be satisfied using a cell union
   * because cell unions are automatically normalized by replacing four child cells with their
   * parent whenever possible. (Note that the list of cell ids passed to the cell union constructor
   * does in fact satisfy all the given restrictions.)
   */
  public S2CellUnion getCovering(S2Region region) {
    S2CellUnion covering = new S2CellUnion();
    getCovering(region, covering);
    return covering;
  }

  @JsMethod(name = "getCoveringCellUnion")
  public void getCovering(S2Region region, S2CellUnion covering) {
    ActiveCovering state = new ActiveCovering(false, region);
    state.getCoveringInternal();
    covering.initSwap(state.result);
  }

  /**
   * Return a normalized cell union that is contained within the given region and satisfies the
   * restrictions *EXCEPT* for minLevel() and levelMod().
   */
  public S2CellUnion getInteriorCovering(S2Region region) {
    S2CellUnion covering = new S2CellUnion();
    getInteriorCovering(region, covering);
    return covering;
  }

  @JsMethod(name = "getInteriorCoveringCellUnion")
  public void getInteriorCovering(S2Region region, S2CellUnion covering) {
    ActiveCovering state = new ActiveCovering(true, region);
    state.getCoveringInternal();
    covering.initSwap(state.result);
  }

  /**
   * Given a connected region and a starting point, return a set of cells at the given level that
   * cover the region.
   */
  public static void getSimpleCovering(
      S2Region region, S2Point start, int level, ArrayList<S2CellId> output) {
    floodFill(region, S2CellId.fromPoint(start).parent(level), output);
  }

  /**
   * Like getCovering(), except that this method is much faster and the coverings are not as tight.
   *
   * <p>All of the usual parameters are respected (maxCells, minLevel, maxLevel, and levelMod),
   * except that the implementation makes no attempt to take advantage of large values of maxCells.
   * (A small number of cells will always be returned.)
   *
   * <p>This function is useful as a starting point for algorithms that recursively subdivide cells.
   */
  public void getFastCovering(S2Cap cap, ArrayList<S2CellId> results) {
    getRawFastCovering(cap, maxCells(), results);
    normalizeCovering(results);
  }

  /**
   * Compute a covering of the given cap. In general the covering consists of at most 4 cells
   * (except for very large caps, which may need up to 6 cells). The output is not sorted.
   *
   * <p>{@code maxCellsHint} can be used to request a more accurate covering (but is currently
   * ignored).
   */
  private static void getRawFastCovering(
      S2Cap cap, @SuppressWarnings("unused") int maxCellsHint, List<S2CellId> covering) {
    // TODO(user): The covering could be made quite a bit tighter by mapping the cap to a
    // rectangle in (i,j)-space and finding a covering for that.
    covering.clear();

    // Find the maximum level such that the cap contains at most one cell vertex and such that
    // S2CellId.appendVertexNeighbors() can be called.
    int level = S2Projections.PROJ.minWidth.getMaxLevel(2 * cap.angle().radians());
    level = min(level, S2CellId.MAX_LEVEL - 1);

    if (level == 0) {
      // Don't bother trying to optimize the level == 0 case, since more than four face cells may be
      // required.
      Collections.addAll(covering, S2CellId.FACE_CELLS);
    } else {
      // The covering consists of the 4 cells at the given level that share the cell vertex that is
      // closest to the cap center.
      S2CellId id = S2CellId.fromPoint(cap.axis());
      id.getVertexNeighbors(level, covering);
    }
  }

  /**
   * Normalize "covering" so that it conforms to the current covering parameters (maxCells,
   * minLevel, maxLevel, and levelMod).
   */
  public void normalizeCovering(ArrayList<S2CellId> covering) {
    // This method makes no attempt to be optimal. In particular, if minMevel() > 0 or levelMod()
    // > 1, then it may return more than the desired number of cells even when this isn't necessary.
    //
    // Note that when the covering parameters have their default values, almost all of the code in
    // this function is skipped.

    // If any cells are too small, or don't satisfy levelMod(), then replace them with ancestors.
    if (maxLevel() < S2CellId.MAX_LEVEL || levelMod() > 1) {
      for (int i = 0; i < covering.size(); i++) {
        S2CellId id = covering.get(i);
        int level = id.level();
        int newLevel = adjustLevel(min(level, maxLevel()));
        if (newLevel != level) {
          covering.set(i, id.parent(newLevel));
        }
      }
    }

    // Sort the cells and simplify them.
    S2CellUnion.normalize(covering);

    // If there are still too many cells, then repeatedly replace two adjacent cells in S2CellId
    // order by their lowest common ancestor.
    while (covering.size() > maxCells()) {
      int bestIndex = -1;
      int bestLevel = -1;
      for (int i = 0; i + 1 < covering.size(); i++) {
        int level = covering.get(i).getCommonAncestorLevel(covering.get(i + 1));
        level = adjustLevel(level);
        if (level > bestLevel) {
          bestLevel = level;
          bestIndex = i;
        }
      }
      if (bestLevel < minLevel()) {
        break;
      }
      covering.set(bestIndex, covering.get(bestIndex).parent(bestLevel));
      S2CellUnion.normalize(covering);
    }

    // Make sure that the covering satisfies minLevel() and levelMod(), possibly at the expense of
    // satisfying maxCells().
    if (minLevel() > 0 || levelMod() > 1) {
      S2CellUnion result = new S2CellUnion();
      result.initRawSwap(covering);
      result.denormalize(minLevel(), levelMod(), covering);
    }
  }

  /**
   * If level > minLevel(), then reduce {@code level} if necessary so that it also satisfies
   * levelMod(). Levels smaller than minLevel() are not affected (since cells at these levels are
   * eventually expanded).
   */
  private int adjustLevel(int level) {
    if (levelMod() > 1 && level > minLevel()) {
      level -= (level - minLevel()) % levelMod();
    }
    return level;
  }

  /** This class tracks the state of a covering while it is underway. */
  final class ActiveCovering {
    /** True if we're covering the interior. */
    final boolean interiorCovering;

    /** The region being covered. */
    final S2Region region;

    /** Counter of number of candidates created, for performance evaluation. */
    int candidatesCreatedCounter = 0;

    /** Cell ids that have been added to the covering so far. */
    final ArrayList<S2CellId> result = new ArrayList<>();

    /** Prioritized candidates to explore next. */
    final PriorityQueue<QueueEntry> candidateQueue =
        new PriorityQueue<>(10, new QueueEntriesComparator());

    ActiveCovering(boolean interior, S2Region region) {
      this.interiorCovering = interior;
      this.region = region;
    }

    /**
     * If the cell intersects the given region, return a new candidate with no children, otherwise
     * return null. Also marks the candidate as "terminal" if it should not be expanded further.
     */
    @Nullable
    private Candidate newCandidate(S2Cell cell) {
      if (!region.mayIntersect(cell)) {
        return null;
      }

      boolean isTerminal = false;
      if (cell.level() >= minLevel) {
        if (interiorCovering) {
          if (region.contains(cell)) {
            isTerminal = true;
          } else if (cell.level() + levelMod > maxLevel) {
            return null;
          }
        } else {
          if (cell.level() + levelMod > maxLevel || region.contains(cell)) {
            isTerminal = true;
          }
        }
      }
      Candidate candidate = new Candidate();
      candidate.cell = cell;
      candidate.isTerminal = isTerminal;
      if (!isTerminal) {
        candidate.children = new Candidate[1 << maxChildrenShift()];
      }
      candidatesCreatedCounter++;
      return candidate;
    }

    /** Return the log base 2 of the maximum number of children of a candidate. */
    private int maxChildrenShift() {
      return 2 * levelMod;
    }

    /**
     * Process a candidate by either adding it to the result list or expanding its children and
     * inserting it into the priority queue. Passing a null argument does nothing.
     */
    private void addCandidate(Candidate candidate) {
      if (candidate == null) {
        return;
      }

      if (candidate.isTerminal) {
        result.add(candidate.cell.id());
        return;
      }
      // assert (candidate.numChildren == 0);

      // Expand one level at a time until we hit minLevel to ensure that
      // we don't skip over it.
      int numLevels = (candidate.cell.level() < minLevel) ? 1 : levelMod;
      int numTerminals = expandChildren(candidate, candidate.cell, numLevels);

      if (candidate.numChildren == 0) {
        // Do nothing
      } else if (!interiorCovering
          && numTerminals == 1 << maxChildrenShift()
          && candidate.cell.level() >= minLevel) {
        // Optimization: add the parent cell rather than all of its children. We can't do this for
        // interior coverings, since the children just intersect the region, but may not be
        // contained by it - we need to subdivide them further.
        candidate.isTerminal = true;
        addCandidate(candidate);

      } else {
        // We negate the priority so that smaller absolute priorities are returned first. The
        // heuristic is designed to refine the largest cells first, since those are where we have
        // the largest potential gain. Among cells at the same level, we prefer the cells with the
        // smallest number of intersecting children. Finally, we prefer cells that have the smallest
        // number of children that cannot be refined any further.
        int priority =
            -((((candidate.cell.level() << maxChildrenShift()) + candidate.numChildren)
                    << maxChildrenShift())
                + numTerminals);
        candidateQueue.add(new QueueEntry(priority, candidate));
        // logger.info("Push: " + candidate.cell.id() + " (" + priority + ") ");
      }
    }

    /**
     * Populate the children of "candidate" by expanding the given number of levels from the given
     * cell. Returns the number of children that were marked "terminal".
     */
    private int expandChildren(Candidate candidate, S2Cell cell, int numLevels) {
      numLevels--;
      S2Cell[] childCells = new S2Cell[4];
      for (int i = 0; i < 4; ++i) {
        childCells[i] = new S2Cell();
      }
      cell.subdivide(childCells);
      int numTerminals = 0;
      for (int i = 0; i < 4; ++i) {
        if (numLevels > 0) {
          if (region.mayIntersect(childCells[i])) {
            numTerminals += expandChildren(candidate, childCells[i], numLevels);
          }
          continue;
        }
        Candidate child = newCandidate(childCells[i]);
        if (child != null) {
          candidate.children[candidate.numChildren++] = child;
          if (child.isTerminal) {
            ++numTerminals;
          }
        }
      }
      return numTerminals;
    }

    /** Computes a set of initial candidates that cover the given region. */
    private void getInitialCandidates() {
      // Optimization: if at least 4 cells are desired (the normal case), start with a 4-cell
      // covering of the region's bounding cap. This lets us skip quite a few levels of refinement
      // when the region to be covered is relatively small.
      if (maxCells >= 4) {
        // Find the maximum level such that the bounding cap contains at most one cell vertex at
        // that level.
        S2Cap cap = region.getCapBound();
        int level =
            min(
                PROJ.minWidth.getMaxLevel(2 * cap.angle().radians()),
                maxLevel(),
                S2CellId.MAX_LEVEL - 1);
        if (levelMod() > 1 && level > minLevel()) {
          level -= (level - minLevel()) % levelMod();
        }
        // We don't bother trying to optimize the level == 0 case, since more than four face cells
        // may be required.
        if (level > 0) {
          // Find the leaf cell containing the cap axis, and determine which subcell of the parent
          // cell contains it.
          ArrayList<S2CellId> base = new ArrayList<>(4);
          S2CellId id = S2CellId.fromPoint(cap.axis());
          id.getVertexNeighbors(level, base);
          for (int i = 0; i < base.size(); ++i) {
            addCandidate(newCandidate(new S2Cell(base.get(i))));
          }
          return;
        }
      }
      // Default: start with all six cube faces.
      for (int face = 0; face < 6; ++face) {
        addCandidate(newCandidate(FACE_CELLS.get(face)));
      }
    }

    /** Generates a covering and stores it in result. */
    private void getCoveringInternal() {
      // Strategy: Start with the 6 faces of the cube. Discard any that do not intersect the shape.
      // Then repeatedly choose the largest cell that intersects the shape and subdivide it.
      //
      // "result" contains the cells that will be part of the output, while the priority queue
      // contains cells that we may still subdivide further. Cells that are entirely contained
      // within the region are immediately added to the output, while cells that do not intersect
      // the region are immediately discarded.
      //
      // Therefore candidateQueue only contains cells that partially intersect the region.
      // Candidates are prioritized first according to cell size (larger cells first), then by the
      // number of intersecting children they have (fewest children first), and then by the number
      // of fully contained children (fewest children first).

      Preconditions.checkState(candidateQueue.isEmpty() && result.isEmpty());

      getInitialCandidates();
      while (!candidateQueue.isEmpty() && (!interiorCovering || result.size() < maxCells)) {
        Candidate candidate = candidateQueue.poll().candidate;
        // For interior covering we keep subdividing no matter how many children candidate has. If
        // we reach maxCells before expanding all children, we will just use some of them. For
        // exterior covering we cannot do this, because result has to cover the whole region, so all
        // children have to be used. The candidate.numChildren == 1 case takes care of the situation
        // when we already have more than maxCells in result (minLevel is too high).
        if (interiorCovering
            || candidate.cell.level() < minLevel
            || candidate.numChildren == 1
            || result.size() + candidateQueue.size() + candidate.numChildren <= maxCells) {
          // Expand this candidate into its children.
          for (int i = 0; i < candidate.numChildren; ++i) {
            if (!interiorCovering || result.size() < maxCells) {
              addCandidate(candidate.children[i]);
            }
          }
        } else {
          candidate.isTerminal = true;
          addCandidate(candidate);
        }
      }
    }
  }

  /**
   * Given a region and a starting cell, return the set of all the edge-connected cells at the same
   * level that intersect "region". The output cells are returned in arbitrary order.
   */
  private static void floodFill(S2Region region, S2CellId start, ArrayList<S2CellId> output) {
    HashSet<S2CellId> all = new HashSet<>();
    ArrayList<S2CellId> frontier = new ArrayList<>();
    output.clear();
    all.add(start);
    frontier.add(start);
    while (!frontier.isEmpty()) {
      S2CellId id = frontier.get(frontier.size() - 1);
      frontier.remove(frontier.size() - 1);
      if (!region.mayIntersect(new S2Cell(id))) {
        continue;
      }
      output.add(id);

      S2CellId[] neighbors = new S2CellId[4];
      id.getEdgeNeighbors(neighbors);
      for (int edge = 0; edge < 4; ++edge) {
        S2CellId nbr = neighbors[edge];
        if (!all.contains(nbr)) {
          frontier.add(nbr);
          all.add(nbr);
        }
      }
    }
  }
}
