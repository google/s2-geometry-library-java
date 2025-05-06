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
package com.google.common.geometry.primitives;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.primitives.IdSetLexicon.IdSet;
import com.google.common.geometry.primitives.IdSetLexicon.TestIdSet;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.IntConsumer;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.NoSuchElementException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link IdSetLexicon}. */
@RunWith(JUnit4.class)
public final class IdSetLexiconTest {

  /** Verifies that the contents and size of the "actual" IdSet are the same as "expected". */
  public void assertIdSetsEqual(IdSet expected, IdSet actual) {
    assertEquals(expected.size(), actual.size());
    OfInt iter = expected.intIterator();
    actual.forEach((IntConsumer) actualId -> assertEquals(iter.nextInt(), actualId));
  }

  @Test
  public void testEmptyIdSet() {
    IdSetLexicon lexicon = new IdSetLexicon();
    assertEquals(new TestIdSet(), lexicon.idSet(lexicon.add(ImmutableList.of())));

    IdSet emptyIdSet = lexicon.idSet(IdSetLexicon.EMPTY_SET_ID);
    assertEquals(0, emptyIdSet.size());
    assertTrue(emptyIdSet.isEmpty());

    try {
      int unused = emptyIdSet.first();
      fail();
    } catch (NoSuchElementException expectedException) {
      // Expected.
    }

    OfInt emptyIter = emptyIdSet.intIterator();
    assertFalse(emptyIter.hasNext());
    try {
      int unused = emptyIter.nextInt();
      fail();
    } catch (NoSuchElementException expectedException) {
      // Expected.
    }

    emptyIdSet.forEach((int i) -> fail());
  }

  @Test
  public void testTestIdSet() {
    TestIdSet simple = new TestIdSet(1, 3, 5, 7);
    // Check sort and deduplicate.
    TestIdSet messy = new TestIdSet(1, 3, 1, 7, 3, 1, 5, 3, 1, 5);
    assertIdSetsEqual(simple, messy);

    // Construct from a Collection of Integers.
    TestIdSet fromCollection = new TestIdSet(ImmutableList.of(7, 3, 5, 3, 1, 1, 5));
    assertIdSetsEqual(simple, fromCollection);

    IdSetLexicon lexicon = new IdSetLexicon();

    // Empty TestIdSet.
    assertIdSetsEqual(new TestIdSet(), lexicon.idSet(IdSetLexicon.EMPTY_SET_ID));
  }

  @Test
  public void testSingletonSets() {
    IdSetLexicon lexicon = new IdSetLexicon();
    assertEquals(5, lexicon.add(ImmutableList.of(5)));
    // Deduplication produces singleton sets.
    assertEquals(20, lexicon.add(ImmutableIntSequence.of(20, 20)));
    assertEquals(0, lexicon.add(ImmutableList.of(0, 0, 0)));
    assertEquals(1, lexicon.addSingleton(1));

    assertIdSetsEqual(new TestIdSet(0), lexicon.idSet(0));
    assertIdSetsEqual(new TestIdSet(1), lexicon.idSet(1));
    assertIdSetsEqual(new TestIdSet(5), lexicon.idSet(5));

    // We don't actually have to add singleton sets before obtaining them.
    IdSet setOfEleven = lexicon.idSet(11);
    assertIdSetsEqual(new TestIdSet(11), setOfEleven);
    assertEquals(1, setOfEleven.size());
    assertFalse(setOfEleven.isEmpty());
    assertEquals(11, setOfEleven.first());

    OfInt elevenIter = setOfEleven.intIterator();
    assertTrue(elevenIter.hasNext());
    assertEquals(11, elevenIter.nextInt());
    assertFalse(elevenIter.hasNext());

    setOfEleven.forEach((int i) -> assertEquals(11, i));
  }

  @Test
  public void testSetsAreSorted() {
    IdSetLexicon lexicon = new IdSetLexicon();
    int set25Id = lexicon.add(ImmutableIntSequence.of(2, 5));
    int set325Id = lexicon.add(ImmutableList.of(3, 2, 5));

    // sets are sorted, duplicate sets are the same set
    assertEquals(set25Id, lexicon.add(ImmutableList.of(5, 2)));
    // and duplicate values are ignored.
    assertEquals(set325Id, lexicon.add(ImmutableIntSequence.of(5, 3, 2, 5)));

    // When iterated with forEach, the idSet elements are presented in order.
    assertIdSetsEqual(new TestIdSet(2, 5), lexicon.idSet(set25Id));
    assertIdSetsEqual(new TestIdSet(2, 3, 5), lexicon.idSet(set325Id));
  }

  @Test
  public void testClearAndCopyConstructor() {
    IdSetLexicon lexicon = new IdSetLexicon();
    int set12Id = lexicon.add(ImmutableList.of(1, 2));
    int set34Id = lexicon.add(ImmutableIntSequence.of(3, 4));

    // Make a copy of the lexicon using the copy constructor, and then clear it.
    IdSetLexicon copy = new IdSetLexicon(lexicon);
    lexicon.clear();

    // After clear(), set ids are reused, they don't depend on the set contents, just the order
    // that the sets are added to the lexicon.
    assertEquals(set12Id, lexicon.add(ImmutableList.of(3, 4)));
    assertEquals(set34Id, lexicon.add(ImmutableIntSequence.of(1, 2)));

    // Meanwhile the copy still has the original sets.
    assertIdSetsEqual(new TestIdSet(1, 2), copy.idSet(set12Id));
    assertIdSetsEqual(new TestIdSet(3, 4), copy.idSet(set34Id));
  }

  @Test
  public void testRetrieveNonexistent() {
    IdSetLexicon lexicon = new IdSetLexicon();

    // Try to get a non-singleton, non-empty set id that wasn't added.
    try {
      IdSetLexicon.IdSet unused = lexicon.idSet(~1);
      fail();
    } catch (IllegalArgumentException expectedException) {
      // Expected.
    }
  }

  @Test
  public void testNegativeValuesRejected() {
    IdSetLexicon lexicon = new IdSetLexicon();

    // Add a set with a negative element.
    try {
      int unused = lexicon.add(ImmutableList.of(1, -1, 2));
      fail();
    } catch (IllegalArgumentException expectedException) {
      // Expected.
    }
  }
}
