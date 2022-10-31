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

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2IndexingHelper.Term;
import com.google.common.geometry.S2IndexingHelper.TermType;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import junit.framework.TestCase;

/** @author Pete Gillin (peteg@google.com) */

public class S2IndexingHelperImplTest extends TestCase {

  public void testGetIndexTermsLatLng() {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId = S2CellId.fromLatLng(latlng);
    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 4)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(16)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(10)),
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(latlng);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetIndexTermsLatLng_optimizeForSpace() {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId = S2CellId.fromLatLng(latlng);
    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 4)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(16)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(10)),
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOptimizeForSpace()
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(latlng);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetIndexTermsLatLng_onlyPointsIndexed() {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId = S2CellId.fromLatLng(latlng).parent(16);
    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 4)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(16)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOnlyPointIndexed()
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(latlng);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetIndexTermsIterableCellId() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));
    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(10)),
      // A non-leaf cell gets indexed as a covering term and an ancestor term
      // (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // A non-leaf cell's ancestors get indexed as ancestor terms (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(cellIds);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetIndexTermsIterableCellId_optimizeForSpace() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));
    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(10)),
      // A non-leaf cell gets indexed as a covering term (rule 1)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(14)),
      // A non-leaf cell's ancestors get indexed as ancestor terms (rule 1)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOptimizeForSpace()
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(cellIds);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetIndexTermsIterableCellId_onlyPointsIndexed() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOnlyPointIndexed()
            .build();
    try {
      helper.getIndexTerms(cellIds);
      fail("Expected exception");
    } catch (IllegalStateException e) {
      // expected exception
    }
  }

  public void testGetIndexTermsRegion() {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2CellId cellId = S2CellId.fromLatLng(latlng).parent(15);
    S2Cell cell = new S2Cell(cellId); // Region consisting of a level-15 cell
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;
    // With these options, the level-15 cell is denormalized into its level-16 children
    S2CellId childCellId1 = cellId.childBegin();
    S2CellId childCellId2 = childCellId1.next();
    S2CellId childCellId3 = childCellId2.next();
    S2CellId childCellId4 = childCellId3.next();

    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 4)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId1),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId2),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId3),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId4),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(10)),
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(cell);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetIndexTermsListRegion() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2Region> regions =
        ImmutableList.<S2Region>of(new S2Cell(cellId1.parent(16)), new S2Cell(cellId2.parent(14)));
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    Term[] expectedTerms = {
      // A leaf cell and its ancestors get indexed as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(10)),
      // A non-leaf cell gets indexed as a covering term and an ancestor term (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // A non-leaf cell's ancestors get indexed as ancestor terms (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getIndexTerms(regions);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testQueryTermsLatLng() throws Exception {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId = S2CellId.fromLatLng(latlng);
    Term[] expectedTerms = {
      // A leaf cell gets looked up as an ancestor term (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(16)),
      // A leaf cell's ancestors get looked up as covering terms (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(10)),
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(latlng);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testQueryTermsLatLng_optimizeForSpace() throws Exception {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId = S2CellId.fromLatLng(latlng);
    Term[] expectedTerms = {
      // A leaf cell gets looked up as an ancestor term (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(16)),
      // A leaf cell's ancestors get looked up as covering terms (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(10)),
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOptimizeForSpace()
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(latlng);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testQueryTermsLatLng_onlyPointsIndexed() throws Exception {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId = S2CellId.fromLatLng(latlng);
    Term[] expectedTerms = {
      // A leaf cell gets looked up as an ancestor term (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId.parent(16)),
      // All lookups of covering terms omitted (rule 5)
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOnlyPointIndexed()
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(latlng);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetQueryTermsIterableCellId() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));
    Term[] expectedTerms = {
      // A leaf cell get looked up as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      // A leaf cell's ancestors get looked up as covering terms (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(10)),
      // A non-leaf cell gets looked up as an ancestor term (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // A non-leaf cell's ancestors get looked up as covering terms (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(cellIds);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetQueryTermsIterableCellId_optimizeForSpace() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));
    Term[] expectedTerms = {
      // A leaf cell get looked up as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      // A leaf cell's ancestors get looked up as covering terms (rule 1)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(10)),
      // A non-leaf cell gets looked up as a covering term and as an ancestor term (rule 1)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // A non-leaf cell's ancestors get looked up as covering terms (rule 1)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOptimizeForSpace()
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(cellIds);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetQueryTermsIterableCellId_onlyPointsIndexed() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));
    Term[] expectedTerms = {
      // A leaf cell get looked up as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      // A non-leaf cell gets looked up as an ancestor term (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // All lookups of covering terms omitted (rule 5)
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOnlyPointIndexed()
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(cellIds);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testGetQueryTermsIterableCellId_optimizeForSpaceAndOnlyPointIndexed() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    List<S2CellId> cellIds = ImmutableList.of(cellId1.parent(16), cellId2.parent(14));
    Term[] expectedTerms = {
      // A leaf cell get looked up as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      // A non-leaf cell gets looked up as an ancestor term (rules 1 and 5)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // All lookups of covering terms omitted (rule 5)
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .setOptimizeForSpace()
            .setOnlyPointIndexed()
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(cellIds);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  public void testQueryTermsRegion() throws Exception {
    S2LatLng latlng = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2CellId cellId = S2CellId.fromLatLng(latlng).parent(15);
    S2Cell cell = new S2Cell(cellId); // Region consisting of a level-15 cell
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;
    // With these options, the level-15 cell is denormalized into its level-16 children
    S2CellId childCellId1 = cellId.childBegin();
    S2CellId childCellId2 = childCellId1.next();
    S2CellId childCellId3 = childCellId2.next();
    S2CellId childCellId4 = childCellId3.next();

    Term[] expectedTerms = {
      // A leaf cell gets looked up as an ancestor term (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId1),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId2),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId3),
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, childCellId4),
      // A leaf cell's ancestors get looked up as covering terms (rule 6)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId.parent(10)),
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(cell);

    assertEquals(new HashSet<Term>(Arrays.asList(expectedTerms)), new HashSet<Term>(actualTerms));
  }

  public void testGetQueryTermsListRegion() {
    S2LatLng latlng1 = S2LatLng.fromDegrees(51.49483, -0.14653);
    S2LatLng latlng2 = S2LatLng.fromDegrees(31.232867, 121.47676);
    S2CellId cellId1 = S2CellId.fromLatLng(latlng1);
    S2CellId cellId2 = S2CellId.fromLatLng(latlng2);
    ImmutableList<S2Region> regions =
        ImmutableList.<S2Region>of(new S2Cell(cellId1.parent(16)), new S2Cell(cellId2.parent(14)));
    int minLevel = 10;
    int maxLevel = 16;
    int levelMod = 2;

    Term[] expectedTerms = {
      // A leaf cell get looked up as ancestor terms (rule 3)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId1.parent(16)),
      // A leaf cell's ancestors get looked up as covering terms (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(14)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId1.parent(10)),
      // A non-leaf cell gets looked up as an ancestor term (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.ANCESTOR, cellId2.parent(14)),
      // A non-leaf cell's ancestors get looked up as covering terms (rule 2)
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(12)),
      new S2IndexingHelperImpl.TermImpl(TermType.COVERING, cellId2.parent(10))
    };

    S2IndexingHelperImpl helper =
        S2IndexingHelperImpl.builder()
            .setMinLevel(minLevel)
            .setMaxLevel(maxLevel)
            .setLevelMod(levelMod)
            .build();
    Collection<Term> actualTerms = helper.getQueryTerms(regions);

    checkExpectedVsActual(expectedTerms, actualTerms);
  }

  // Compare contents as sets, because ordering does not matter.
  private void checkExpectedVsActual(Term[] expected, Collection<Term> actual) {
    assertEquals(new HashSet<Term>(Arrays.asList(expected)), new HashSet<Term>(actual));
  }
}
