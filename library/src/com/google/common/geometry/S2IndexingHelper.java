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

import java.util.Collection;
import java.util.List;

/**
 * A helper for S2-based spatial indexing. The aim is to index a collection of document regions, and
 * then to find all document regions which intersect with a query region. To do this, each region is
 * converted into the S2 cells which cover it, and the intersection of the coverings is found.
 * Assuming an exterior covering, this may find some false positives, where the coverings intersect
 * but the regions do not. The rate of false positives is decreased by using smaller cells, but this
 * also increases the size of the index and the complexity of the queries.
 *
 * <p>The coverings may be calculated using a restricted subset of the level in the cell hierarchy
 * (see {@link S2RegionCoverer}. The helper may perform assumption based on this restriction.
 *
 * <p>This helper is designed to be independent of your search platform. You will need to convert
 * the {@link S2IndexingHelper.Term}s into the terms appropriate to your platform. To add a document
 * to an index, you call {@link S2IndexingHelper#getIndexTerms} and you index it with each of the
 * terms that returns. To perform a search, you call {@link S2IndexingHelper#getQueryTerms} and you
 * create a query which looks up documents matching any of those terms (i.e. you OR them together).
 *
 * @author Pete Gillin (peteg@google.com)
 */
public interface S2IndexingHelper {

  /**
   * The type of an index term. See the algorithm in the implementation for variants and details.
   */
  enum TermType {

    /** Term which indexes a cell which is part of the covering of the region. */
    COVERING,

    /** Term which indexes an ancestor of a cell which is part of the covering of the region. */
    ANCESTOR,
  }

  /** An index term. */
  interface Term {

    /**
     * @return The type of the term
     */
    TermType type();

    /**
     * @return The cell ID of the term
     */
    long cellId();
  }

  /**
   * @return The minimum level in use (the one closest to the root of the hierarchy)
   */
  public int minLevel();

  /**
   * @return The maximum level in use (the one closest to the leaves of the hierarchy)
   */
  public int maxLevel();

  /**
   * @return The spacing of the levels in use: only levels an integer multiple of this above the
   *     minimum are used
   */
  public int levelMod();

  /**
   * Gets the terms to be used to index a region consisting of a single cell around a point.
   *
   * @param latlng The coordinates of the point
   * @return The terms
   */
  Collection<Term> getIndexTerms(S2LatLng latlng);

  /**
   * Gets the terms to be used to index a region consisting of a collection of cells.
   *
   * @param cellIds The cell IDs which make up the region (which must all be at levels is use as
   *     defined by {@link #minLevel}, {@link #maxLevel}, and {@link #levelMod})
   * @return The terms
   */
  Collection<Term> getIndexTerms(Iterable<S2CellId> cellIds);

  /**
   * Gets the terms to be used to index a general region. An exterior cell covering of the region
   * will be created.
   *
   * @param region The region
   * @return The terms
   */
  Collection<Term> getIndexTerms(S2Region region);

  /**
   * Gets the terms to be used to index the union of several regions. Exterior cell coverings of the
   * regions will be created.
   *
   * @param regions The regions
   * @return The terms
   */
  Collection<Term> getIndexTerms(List<S2Region> regions);

  /**
   * Gets the terms to be used to search for intersections with a region consisting of a single cell
   * around a point.
   *
   * @param latlng The coordinates of the point
   * @return The terms
   */
  Collection<Term> getQueryTerms(S2LatLng latlng);

  /**
   * Gets the terms to be used to search for intersections with a region consisting of a collection
   * of cells.
   *
   * @param cellIds The cell IDs which make up the region (which must all be at levels is use as
   *     defined by {@link #minLevel}, {@link #maxLevel}, and {@link #levelMod})
   * @return The terms
   */
  Collection<Term> getQueryTerms(Iterable<S2CellId> cellIds);

  /**
   * Gets the terms to be used to search for intersections with a general region. An exterior cell
   * covering of the region will be created.
   *
   * @param region The region
   * @return The terms
   */
  Collection<Term> getQueryTerms(S2Region region);

  /**
   * Gets the terms to be used to search for intersections with the union of several general region.
   * Exterior cell coverings of the regions will be created.
   *
   * @param regions The regions
   * @return The terms
   */
  Collection<Term> getQueryTerms(List<S2Region> regions);
}
