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

import static com.google.common.geometry.primitives.Ints.MIXING_HASHER;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.IntSequenceHasher;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Tests for {@link SequenceLexicon}, the sequences it contains, and the IntSequenceHasher it uses.
 */
@RunWith(JUnit4.class)
public final class SequenceLexiconTest {

  private static final int SEED = 123455;
  private final Random rand = new Random(SEED);

  /**
   * This implementation of IntSequenceHasher was designed for testing the behavior of
   * SequenceLexicon when hash collisions do occur. Every input sequence is hashed to its first
   * value, or zero if the input sequence is empty. The implementation of equals() is delegated to
   * MixingHasher.
   */
  public static class CollidingHasher implements IntSequenceHasher {
    @Override
    public int hashCode(IntSequence iterable) {
      OfInt iter = iterable.intIterator();
      if (!iter.hasNext()) {
        return 0;
      }
      return iter.nextInt();
    }

    @Override
    public boolean equals(IntSequence a, IntSequence b) {
      return MIXING_HASHER.equals(a, b);
    }
  }

  // Exercises the basic operations of SequenceLexicon.
  private void assertLexiconBehavior(SequenceLexicon lex) {
    // Check that the sequence is starting out empty.
    assertEquals(0, lex.numElements());
    assertEquals(0, lex.numBlocks());

    // Sequence ids are assigned sequentially.
    assertEquals(0, lex.add(ImmutableIntSequence.of()));
    assertEquals(1, lex.add(ImmutableIntSequence.of(5)));

    // Two sequences were added so far, but only one data element.
    assertEquals(1, lex.numElements());

    // Empty list is already present, either as a list or an array.
    assertEquals(0, lex.add(ImmutableIntSequence.of()));

    // Add two more sequences and get two more consecutive ids. After this, 5 elements are stored.
    assertEquals(2, lex.add(ImmutableIntSequence.of(5, 5)));
    assertEquals(3, lex.add(ImmutableIntSequence.of(5, 6)));
    assertEquals(5, lex.numElements());

    // When testing with CollidingHasher, sequences 1 through 4 all have the same hash code, but
    // that should not affect correct operation.
    assertEquals(4, lex.add(ImmutableIntSequence.of(5, 0, -3)));

    // singleton sequence '5' is already present.
    assertEquals(1, lex.add(ImmutableIntSequence.of(5)));
    // sequence 5,0,-3 is already present.
    assertEquals(4, lex.add(ImmutableIntSequence.of(5, 0, -3)));
    // sequence 5,6 is already present.
    assertEquals(3, lex.add(ImmutableIntSequence.of(5, 6)));

    // There are now five distinct sequences in the lexicon.
    assertEquals(5, lex.size());

    // Check that each of the five sequences has the expected contents.
    assertTrue(MIXING_HASHER.equals(ImmutableIntSequence.of(), lex.sequence(0)));
    assertTrue(MIXING_HASHER.equals(ImmutableIntSequence.of(5), lex.sequence(1)));
    assertTrue(MIXING_HASHER.equals(ImmutableIntSequence.of(5, 5), lex.sequence(2)));
    assertTrue(MIXING_HASHER.equals(ImmutableIntSequence.of(5, 6), lex.sequence(3)));
    assertTrue(MIXING_HASHER.equals(ImmutableIntSequence.of(5, 0, -3), lex.sequence(4)));
  }

  private void assertEqualsIterators(SequenceLexicon lex) {
    assertEquals(0, lex.add(ImmutableIntSequence.of(100, 101, 102, 103)));
    ImmutableIntSequence s = lex.sequence(0);
    assertEquals(4, s.size());

    OfInt intIterator = s.intIterator();
    assertTrue(intIterator.hasNext());
    assertEquals(100, intIterator.nextInt());

    assertTrue(intIterator.hasNext());
    assertEquals(101, intIterator.nextInt());

    assertTrue(intIterator.hasNext());
    assertEquals(102, intIterator.nextInt());

    assertTrue(intIterator.hasNext());
    assertEquals(103, intIterator.nextInt());

    assertFalse(intIterator.hasNext());
  }

  private void assertIteratorValues(SequenceLexicon lex) {
    assertEquals(0, lex.add(ImmutableIntSequence.of(100, 101, 102, 103, 104)));
    ImmutableIntSequence s = lex.sequence(0);
    assertEquals(5, s.size());

    OfInt i = s.intIterator();
    assertEquals(100, i.nextInt());
    assertEquals(101, i.nextInt());
    assertEquals(102, i.nextInt());
    assertEquals(103, i.nextInt());
    assertEquals(104, i.nextInt());
    assertFalse(i.hasNext());

    // forEach is checked by assertEquals
    assertTrue(MIXING_HASHER.equals(ImmutableIntSequence.of(100, 101, 102, 103, 104), s));
  }

  // Verify that clearing the lexicon causes sequence ids to be reused.
  private void assertLexiconClear(SequenceLexicon lex) {
    assertEquals(0, lex.numElements());
    assertEquals(0, lex.add(ImmutableIntSequence.of(1)));
    assertEquals(1, lex.add(ImmutableIntSequence.of(2)));
    assertEquals(2, lex.size());
    assertEquals(2, lex.numElements());
    lex.clear();
    assertEquals(0, lex.size());
    assertEquals(0, lex.numElements());
    assertEquals(0, lex.add(ImmutableIntSequence.of(2)));
    assertEquals(1, lex.add(ImmutableIntSequence.of(1)));
    assertEquals(2, lex.size());
  }

  // Test SequenceLexicon operations on a normal SequenceLexicon.
  @Test
  public void testSequenceLexicon() {
    SequenceLexicon lex = new SequenceLexicon();
    assertLexiconBehavior(lex);
    lex.clear();
    assertEqualsIterators(lex);
    lex.clear();
    assertIteratorValues(lex);
    lex.clear();
    assertLexiconClear(lex);
  }

  @Test
  public void testSameSequence() {
    SequenceLexicon lexOne = new SequenceLexicon();
    int s1113id = lexOne.add(ImmutableIntSequence.of(1, 1, 1, 3));
    int s1114id = lexOne.add(ImmutableIntSequence.of(1, 1, 1, 4));

    SequenceLexicon lexTwo = new SequenceLexicon();
    int s3113idB = lexTwo.add(ImmutableIntSequence.of(3, 1, 1, 3));
    int s1113idB = lexTwo.add(ImmutableIntSequence.of(1, 1, 1, 3));

    // Equal sequences in the same lexicon (will be referring to the same data)
    assertTrue(MIXING_HASHER.equals(lexOne.sequence(s1113id), lexOne.sequence(s1113id)));

    // Different sequences in the same lexicon
    assertFalse(MIXING_HASHER.equals(lexOne.sequence(s1113id), lexOne.sequence(s1114id)));

    // Equal sequences in different lexicons
    assertTrue(MIXING_HASHER.equals(lexOne.sequence(s1113id), lexTwo.sequence(s1113idB)));

    // Different sequences in different lexicons, but they have the same sequence id
    assertEquals(s3113idB, s1113id);
    assertFalse(MIXING_HASHER.equals(lexTwo.sequence(s3113idB), lexOne.sequence(s1113id)));

    // ImmutableIntSequence.of equals
    assertTrue(MIXING_HASHER.equals(lexOne.sequence(s1114id), ImmutableIntSequence.of(1, 1, 1, 4)));
    assertFalse(
        MIXING_HASHER.equals(lexOne.sequence(s1114id), ImmutableIntSequence.of(2, 1, 1, 4)));
    assertFalse(
        MIXING_HASHER.equals(lexOne.sequence(s1114id), ImmutableIntSequence.of(1, 1, 1, 7)));
  }

  private ImmutableList<Integer> generateRandomSequence(int length, int range) {
    ImmutableList.Builder<Integer> b = new ImmutableList.Builder<>();
    for (int i = 0; i < length; i++) {
      b.add(rand.nextInt(range));
    }
    return b.build();
  }

  // Generate 100 random sequences and store them in the sequence lexicon. Verify that they can be
  // retrieved by id. Make a deep copy of the sequence lexicon and check it as well.
  @Test
  public void testSequenceLexiconConstructFromRandom() {
    SequenceLexicon lex = new SequenceLexicon();
    HashMap<Integer, ImmutableList<Integer>> randomSequences = new HashMap<>();

    // Store 100 sequences in the lexicon, and also in a HashMap.
    for (int i = 0; i < 1; i++) {
      ImmutableList<Integer> seq = generateRandomSequence(rand.nextInt(32) + 1, 10000);
      int id = lex.add(seq);
      randomSequences.put(id, seq);
    }

    // Iterate over the HashMap and verify sequences looked up in the Lexicon match.
    for (Map.Entry<Integer, ImmutableList<Integer>> entry : randomSequences.entrySet()) {
      // Convert the ImmutableList<Integer> to a Sequence to call equals().
      ImmutableIntSequence expectedSequence = ImmutableIntSequence.of(entry.getValue());
      assertEquals(expectedSequence, lex.sequence(entry.getKey()));
    }

    // Make a deep copy of the SequenceLexicon, and then clear the original.
    SequenceLexicon copiedLex = new SequenceLexicon(lex);
    lex.clear();

    // Verify that the copy contains the same sequences.
    for (Map.Entry<Integer, ImmutableList<Integer>> entry : randomSequences.entrySet()) {
      ImmutableIntSequence expectedSequence = ImmutableIntSequence.of(entry.getValue());
      assertEquals(expectedSequence, copiedLex.sequence(entry.getKey()));
    }
  }

  @Test
  public void testSequenceLexiconReallocate() {
    SequenceLexicon lex = new SequenceLexicon();
    // Add exactly enough sequence data that the data array will be full. It has room for 32 blocks
    // of 2048 entries each, for a total of 65536 entries. Add 4K Sequences of size 16.
    for (int i = 0; i < 4096; i++) {
      lex.add(generateRandomSequence(16, 10000));
    }

    // 4096 sequences.
    assertEquals(4096, lex.size());
    // A total of 64K sequence data elements.
    assertEquals(65536, lex.numElements());
    // Stored in 32 blocks of size 2048, and there is no more room in data for additional blocks.
    assertEquals(2048, lex.blockSize());
    assertEquals(32, lex.dataLength());
    assertEquals(32, lex.numBlocks());

    // The 4K sequence begin points fit in two blocks, but begins are written in advance, so the
    // third block is already allocated.
    assertEquals(3, lex.numBeginsBlocks());
    // There are room for 32 blocks of begins.
    assertEquals(32, lex.beginsLength());

    // Now write one more sequence to the lexicon, which will trigger reallocation.
    lex.add(generateRandomSequence(16, 10000));

    // Check that everything updated as expected.
    assertEquals(4097, lex.size());
    assertEquals(65552, lex.numElements());
    // In particular, the data.length was increased to 64, and 33 blocks are now used.
    assertEquals(64, lex.dataLength());
    assertEquals(33, lex.numBlocks());
  }

  /**
   * Throwing an exception on hash collisions is only intended for testing. Check that when it is
   * enabled, it works.
   */
  @Test
  public void testCollisionDetectionThrowsIfEnabled() {
    SequenceLexicon collidingLexicon =
        new SequenceLexicon(new CollidingHasher(), /* throwOnCollision= */ true);

    // Even with CollidingHasher there should be no collisions for these first three sequences.
    collidingLexicon.add(ImmutableIntSequence.of(1, 1, 1, 4));
    collidingLexicon.add(ImmutableIntSequence.of(2, 1, 1, 4));
    collidingLexicon.add(ImmutableIntSequence.of(3, 1, 1, 4));

    try {
      // CollidingHasher hashes to the first value of a list, so 2214 will collide with 2114.
      collidingLexicon.add(ImmutableIntSequence.of(2, 2, 1, 4));
      fail();
    } catch (IllegalStateException expectedException) {
      // Expected.
    }

    // Same again but adding an Sequence which will collide with the sequence 1114.
    try {
      collidingLexicon.add(ImmutableIntSequence.of(1, 2));
      fail();
    } catch (IllegalStateException expectedException) {
      // Expected.
    }
  }

  /**
   * Generate and add 1000 4-element sequences of uniformly distributed random single digit numbers
   * to the lexicon. There are 10,000 possible such random sequences. The expected number of unique
   * values in a set of k values chosen with replacement from n possible values is:
   *
   * {@snippet :
   *    E = n * [ 1 - ( (n-1)/n )^k ]
   *
   * With n = 10000 and k = 1000, we expect
   *
   *    E = 10000 * [ 1 - 0.904832894 ]
   *     ~= 951
   * }
   *
   * So, there should be about 1000 - 951 = 49 duplicates. Check that the sequence lexicon detects
   * close to the expected number of duplicates.
   */
  private void assertManySequences(SequenceLexicon lex) {
    int duplicates = 0;
    int largestSequenceId = lex.size();
    for (int i = 0; i < 1000; i++) {
      // Sequence ids are consecutive unless the sequence is already present.
      int sequenceId = lex.add(generateRandomSequence(4, 10));
      if (sequenceId != largestSequenceId + 1) {
        duplicates++;
      } else {
        largestSequenceId = sequenceId;
      }
    }
    System.out.println(
        "Added 1000 random four-digit sequences and found "
            + duplicates
            + " duplicates, expected about 49.");
    assertTrue(duplicates > 44);
    assertTrue(duplicates < 54);
  }

  @Test
  public void testManySequences() {
    SequenceLexicon lex = new SequenceLexicon();
    assertManySequences(lex);
  }

  @Test
  public void testCorrectBehaviorWithHashCollisions() {
    CollidingHasher collidingHasher = new CollidingHasher();

    // First, check that if hash collisions are expected to throw exceptions, they do.
    SequenceLexicon collidingLexicon = new SequenceLexicon(collidingHasher, true);
    // There should be no collisions for these first three sequences.
    collidingLexicon.add(ImmutableIntSequence.of(1, 1, 1, 4));
    collidingLexicon.add(ImmutableIntSequence.of(2, 1, 1, 4));
    collidingLexicon.add(ImmutableIntSequence.of(3, 1, 1, 4));

    try {
      // CollidingHasher hashes to the first value of a list, so 2214 will collide with 2114.
      collidingLexicon.add(ImmutableIntSequence.of(2, 2, 1, 4));
      fail();
    } catch (IllegalStateException expectedException) {
      // Expected.
    }

    // Now test that lexicon behavior is correct even if there are collisions.
    collidingLexicon = new SequenceLexicon(collidingHasher, false);
    assertLexiconBehavior(collidingLexicon);
    collidingLexicon.clear();
    assertEqualsIterators(collidingLexicon);
    collidingLexicon.clear();
    assertIteratorValues(collidingLexicon);
    collidingLexicon.clear();
    assertLexiconClear(collidingLexicon);

    assertManySequences(collidingLexicon);
  }
}
