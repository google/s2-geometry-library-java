/*
 * Copyright 2010 Google Inc.
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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

/**
 * Implementation of {@link S2IndexingHelper}. See the documentation of that class for a description
 * of the problem being addressed.
 *
 * <p>Since S2 cells are hierarchical, a cell intersects only with itself, its ancestors, and its
 * descendants. This fundamental principle is at the basis of this algorithm.
 *
 * <p>RULE 1: The conceptually simplest approach is as follows. Given the covering of a document
 * region:
 *
 * <ul>
 *   <li>we index each cell as a covering term; and
 *   <li>we index each ancestor of each cell as an ancestor term.
 * </ul>
 *
 * <p>To perform a search given the covering of a query region, we combine the following result
 * sets:
 *
 * <ul>
 *   <li>we find all the document regions which intersect because a query cell is the same as a
 *       document cell, by looking up each cell as a covering term;
 *   <li>we find all the document regions which intersect because a query cell is an ancestor of a
 *       document cell, by looking up each cell as an ancestor term; and
 *   <li>we find all the document regions which intersect because a query cell is a descendent of a
 *       document cell (i.e. because a document cell is an ancestor of a query cell), by looking up
 *       each ancestor of each cell as a covering term.
 * </ul>
 *
 * <p>We refer to this approach as 'optimize for space', by contrast with the following variant.
 *
 * <p>RULE 2: An alternative approach indexes each cell in the covering of a document region as an
 * ancestor term as well as indexing it as a covering term. This means that the first of the three
 * lookups above is not necessary, because it is subsumed by the second. That is, given the covering
 * of a document region:
 *
 * <ul>
 *   <li>we index each cell as both a covering term and an ancestor term; and
 *   <li>we index each ancestor of each cell as an ancestor term.
 * </ul>
 *
 * <p>To perform a search given the covering of a query region:
 *
 * <ul>
 *   <li>we find all the document regions which intersect because a query cell is the same as or an
 *       ancestor of a document cell, by looking up each cell as an ancestor term; and
 *   <li>we find all the document regions which intersect because a query cell is an descendent of a
 *       document cell, by looking up each ancestor of each cell as a covering term.
 * </ul>
 *
 * <p>This variant increases the size of the index and reduces the complexity of the queries. So we
 * to this approach as 'optimize for complexity', by contrast with the previous 'optimize for space'
 * approach.
 *
 * <p>RULE 3: A leaf cell (i.e. one at the maximum level we are dealing with) will never be the
 * ancestor of another cell. ` *
 *
 * <ul>
 *   <li>In 'optimize for complexity' mode, where we would normally index a covering cell in a
 *       document region as both a covering term and an ancestor term, we can safely index a leaf
 *       covering cell only as an ancestor term. This reduces the index size without causing any
 *       change in the query or results from the search.
 *   <li>In 'optimize for space' mode, where we would normally index a covering cell in a document
 *       region as a covering term, we can safely index a leaf covering cell as either a covering
 *       cell or an ancestor cell, without causing any change in the query or results. We choose to
 *       index it as an ancestor cell, to make subsequent optimizations possible.
 *   <li>In 'optimize for space' mode, where we would normally look up a covering cell in a query
 *       region as both a covering term and an ancestor term, we can safely look up a leaf covering
 *       cell only as an ancestor term. This reduces the query complexity without causing any change
 *       in the results.
 *   <li>In 'optimize for complexity' mode, we continue to look up a leaf cell as an ancestor term.
 * </ul>
 *
 * <p>RULE 4: A point is most accurately represented by the single leaf cell which covers it. We can
 * see from rules 1, 2, and 3 that we will index this leaf cell and all its ancestors as ancestor
 * terms. This is the case whether we are in 'optimize for space' or 'optimize for complexity' mode.
 * Note that we create no covering terms.
 *
 * <p>RULE 5: If we know that only points have been indexed (as their leaf cells), rule 4 shows that
 * the index will contain no covering terms. Therefore, we can safely omit all covering terms from
 * the lookups.
 *
 * <p>RULE 6: A point is most accurately represented by the single leaf cell which covers it. We can
 * see from rules 1, 2, and 3 that we will search for this by looking up the leaf cell as an
 * ancestor term and looking up all its ancestors as covering terms. This is the case whether we are
 * in 'optimize for space' or 'optimize for complexity' mode.
 *
 * <p>This class is thread-safe. However, this is achieved through synchronization of the
 * calculation of coverings, so if you need high throughput from concurrent usage of the methods
 * which take {@link S2Region}s then you should consider thread-local instances. The {@link Builder}
 * class is thread-unsafe.
 *
 * @author Pete Gillin (peteg@google.com)
 */
public class S2IndexingHelperImpl implements S2IndexingHelper {

  /** A builder for {@link S2IndexingHelper} instances. */
  public static class Builder {

    private S2RegionCoverer.Builder covererBuilder;
    private boolean onlyPointsIndexed = false;
    private boolean optimizeForSpace = false;

    // Private constructor, use builder().
    private Builder() {
      covererBuilder = S2RegionCoverer.builder();
    }

    /**
     * Build the {@link S2IndexingHelper} instance. Invalidates this builder.
     *
     * @return The instance
     */
    public S2IndexingHelperImpl build() {
      S2RegionCoverer coverer = covererBuilder.build();
      S2IndexingHelperImpl impl =
          new S2IndexingHelperImpl(coverer, onlyPointsIndexed, optimizeForSpace);
      covererBuilder = null;
      return impl;
    }

    /**
     * Sets the minimum level to be used (the one closest to the root of the hierarchy). See
     * {@link S2RegionCoverer#setMinLevel()}.
     *
     * @param level The minimum level
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setMinLevel(int level) {
      covererBuilder.setMinLevel(level);
      return this;
    }

    /**
     * Sets the maximum level to be used (the one closest to the leaves of the hierarchy). See
     * {@link S2RegionCoverer#setMaxLevel()}.
     *
     * @param level The minimum level
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setMaxLevel(int level) {
      covererBuilder.setMaxLevel(level);
      return this;
    }

    /**
     * Sets the spacing of the levels in use: only levels an integer multiple of this above the
     * minimum will be used. See {@link S2RegionCoverer#setLevelMod()}.
     *
     * @param spacing The level spacing
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setLevelMod(int spacing) {
      covererBuilder.setLevelMod(spacing);
      return this;
    }

    /**
     * Sets a limit on the number of cells to be used in a covering. This is a soft limit: it may be
     * exceeded if it is not possible to cover the region in this number of cells, or if it would
     * only be possible to use cells larger than the minimum level (i.e. {@link #setMinLevel} takes
     * precedence over {@link #setMaxCells}). See {@link S2RegionCoverer#setMaxCells()}.
     *
     * @param limit The limit on the number of cells
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setMaxCells(int limit) {
      covererBuilder.setMaxCells(limit);
      return this;
    }

    /**
     * Use 'optimize for complexity' mode (see above). This is the default.
     *
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setOptimizeForComplexity() {
      optimizeForSpace = false;
      return this;
    }

    /**
     * Use 'optimize for space' mode (see above).
     *
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setOptimizeForSpace() {
      optimizeForSpace = true;
      return this;
    }

    /**
     * Allow the indexing of both points and regions (see above). This is the default.
     *
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setPointAndRegionsIndexed() {
      onlyPointsIndexed = false;
      return this;
    }

    /**
     * Allow the indexing only of points. Causes {@link S2IndexingHelper#getIndexTerms(S2Region)}
     * and {@link S2IndexingHelper#getIndexTerms(Iterable)} of S2CellId to throw exceptions.
     *
     * @return This builder
     */
    @CanIgnoreReturnValue
    public Builder setOnlyPointIndexed() {
      onlyPointsIndexed = true;
      return this;
    }
  }

  /** @return A new {@link Builder} with default values */
  public static Builder builder() {
    return new Builder();
  }

  /** Simple pojo implementation of {@code Term}. */
  @VisibleForTesting
  public static class TermImpl implements Term {

    private final TermType type;
    private final long cellId;

    public TermImpl(TermType type, S2CellId cellId) {
      this.type = type;
      this.cellId = cellId.id();
    }

    @Override
    public TermType type() {
      return type;
    }

    @Override
    public long cellId() {
      return cellId;
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (!(o instanceof TermImpl)) {
        return false;
      }

      TermImpl term = (TermImpl) o;

      if (cellId != term.cellId) {
        return false;
      }
      if (type != term.type) {
        return false;
      }

      return true;
    }

    @Override
    public int hashCode() {
      int result;
      result = (type != null ? type.hashCode() : 0);
      result = 31 * result + (int) (cellId ^ (cellId >>> 32));
      return result;
    }

    @Override
    public String toString() {
      return String.format("%s:%016X", type.name(), cellId);
    }
  }

  private final S2RegionCoverer coverer;
  private final boolean onlyPointsIndexed;
  private final boolean optimizeForSpace;

  // Private constructor, use the Builder.
  private S2IndexingHelperImpl(
      S2RegionCoverer coverer, boolean onlyPointsIndexed, boolean optimizeForSpace) {
    this.coverer = coverer;
    this.onlyPointsIndexed = onlyPointsIndexed;
    this.optimizeForSpace = optimizeForSpace;
  }

  @Override
  public int minLevel() {
    return coverer.minLevel();
  }

  @Override
  public int maxLevel() {
    return coverer.maxLevel();
  }

  @Override
  public int levelMod() {
    return coverer.levelMod();
  }

  @Override
  public Collection<Term> getIndexTerms(S2LatLng latlng) {
    // See Rule 4 in the class javadoc, about indexing points.
    S2CellId cellId = S2CellId.fromLatLng(latlng).parent(maxLevel());
    ImmutableList.Builder<Term> termBuilder = ImmutableList.builder();
    for (int level = minLevel(); level <= maxLevel(); level += levelMod()) {
      termBuilder.add(new TermImpl(TermType.ANCESTOR, cellId.parent(level)));
    }
    return termBuilder.build();
  }

  @Override
  public Collection<Term> getIndexTerms(Iterable<S2CellId> cellIds) {
    // See Rules 1 and 2 in the class javadoc, about the basic approaches.
    Preconditions.checkState(
        !onlyPointsIndexed,
        "Attempting to index a region when helper has been told only to index points");
    ImmutableList.Builder<Term> termBuilder = ImmutableList.builder();
    Set<S2CellId> ancestors = Sets.newHashSet();
    for (S2CellId cellId : cellIds) {
      int level = cellId.level();
      checkValidLevel(level);
      if (canBeAncestor(level)) {
        termBuilder.add(new TermImpl(TermType.COVERING, cellId));
        if (!optimizeForSpace) {
          termBuilder.add(new TermImpl(TermType.ANCESTOR, cellId));
        }
      } else {
        // See Rule 3 in the class javadoc, about indexing leaf cells.
        termBuilder.add(new TermImpl(TermType.ANCESTOR, cellId));
      }
      buildTermsForAncestors(cellId, level, termBuilder, ancestors, TermType.ANCESTOR);
    }
    return termBuilder.build();
  }

  @Override
  public Collection<Term> getIndexTerms(S2Region region) {
    S2CellUnion cellUnion = regionToCellUnion(region);
    return getIndexTerms(getDenormalizedCellIds(cellUnion));
  }

  @Override
  public Collection<Term> getIndexTerms(List<S2Region> regions) {
    S2CellUnion cellUnion = regionListToCellUnion(regions);
    return getIndexTerms(getDenormalizedCellIds(cellUnion));
  }

  @Override
  public Collection<Term> getQueryTerms(S2LatLng latlng) {
    // See Rule 6 in the class javadoc, about queries for points.
    S2CellId cellId = S2CellId.fromLatLng(latlng).parent(maxLevel());
    ImmutableList.Builder<Term> termBuilder = ImmutableList.builder();
    for (int level = minLevel(); level <= maxLevel(); level += levelMod()) {
      if ((maxLevel() - level) < levelMod()) {
        termBuilder.add(new TermImpl(TermType.ANCESTOR, cellId.parent(level)));
      } else if (!onlyPointsIndexed) {
        // See Rule 5 in the class javadoc, about queries when the index contains only points.
        termBuilder.add(new TermImpl(TermType.COVERING, cellId.parent(level)));
      }
    }
    return termBuilder.build();
  }

  @Override
  public Collection<Term> getQueryTerms(Iterable<S2CellId> cellIds) {
    // See Rules 1 and 2 in the class javadoc, about the basic approaches.
    ImmutableList.Builder<Term> termBuilder = ImmutableList.builder();
    Set<S2CellId> ancestors = Sets.newHashSet();
    for (S2CellId cellId : cellIds) {
      int level = cellId.level();
      checkValidLevel(level);
      termBuilder.add(new TermImpl(TermType.ANCESTOR, cellId));
      if (!onlyPointsIndexed) {
        // See Rule 5 in the class javadoc, about queries when the index contains only points.
        if (optimizeForSpace && canBeAncestor(level)) {
          // See Rule 3 in the class javadoc, about queries for leaf cells.
          termBuilder.add(new TermImpl(TermType.COVERING, cellId));
        }
        buildTermsForAncestors(cellId, level, termBuilder, ancestors, TermType.COVERING);
      }
    }
    return termBuilder.build();
  }

  @Override
  public Collection<Term> getQueryTerms(S2Region region) {
    S2CellUnion cellUnion = regionToCellUnion(region);
    return getQueryTerms(getDenormalizedCellIds(cellUnion));
  }

  @Override
  public Collection<Term> getQueryTerms(List<S2Region> regions) {
    S2CellUnion cellUnion = regionListToCellUnion(regions);
    return getQueryTerms(getDenormalizedCellIds(cellUnion));
  }

  private void checkValidLevel(int level) {
    Preconditions.checkArgument(
        level >= minLevel(), "Found cell at level %s, below the minimum of %s", level, minLevel());
    Preconditions.checkArgument(
        level <= maxLevel(), "Found cell at level %s, above the maximum of %s", level, maxLevel());
    Preconditions.checkArgument(
        ((level - minLevel()) % levelMod()) == 0,
        "Found cell at level %s, not an integer number of %s-level steps over %s",
        level,
        levelMod(),
        minLevel());
  }

  private boolean canBeAncestor(int level) {
    return level < maxLevel();
  }

  private void buildTermsForAncestors(
      S2CellId cellId,
      int level,
      ImmutableList.Builder<Term> termBuilder,
      Set<S2CellId> ancestors,
      TermType type) {
    while (true) {
      level -= levelMod();
      if (level < minLevel()) {
        break;
      }
      cellId = cellId.parent(level);
      if (!ancestors.add(cellId)) {
        break;
      }
      termBuilder.add(new TermImpl(type, cellId));
    }
  }

  // This method is synchronized because S2RegionCoverer#getCovering is thread-unsafe.
  private synchronized S2CellUnion regionToCellUnion(S2Region region) {
    // As an optimization, do not compute the covering of a region which is already a cell union.
    if (region instanceof S2CellUnion) {
      return (S2CellUnion) region;
    } else {
      return coverer.getCovering(region);
    }
  }

  private S2CellUnion regionListToCellUnion(List<S2Region> regions) {
    Preconditions.checkArgument(!regions.isEmpty());
    ArrayList<S2CellId> cellIds =
        Lists.newArrayList(
            Iterables.concat(
                Lists.transform(regions, region -> regionToCellUnion(region).cellIds())));
    S2CellUnion cellUnion = new S2CellUnion();
    cellUnion.initFromCellIds(cellIds);
    return cellUnion;
  }

  private List<S2CellId> getDenormalizedCellIds(S2CellUnion cellUnion) {
    ArrayList<S2CellId> cellIds = Lists.newArrayList();
    cellUnion.denormalize(minLevel(), levelMod(), cellIds);
    return ImmutableList.copyOf(cellIds);
  }
}
