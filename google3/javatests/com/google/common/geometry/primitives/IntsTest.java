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
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.IntIterators;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.IntSequenceHasher;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.HashSet;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link Ints}. */
@RunWith(JUnit4.class)
public final class IntsTest {
  private static final int SEED = 123455;
  private final Random rand = new Random(SEED);

  private void assertNext(OfInt iter, int expected, boolean expectedHasNext) {
    assertTrue(iter.hasNext());
    assertEquals(expected, iter.nextInt());
    assertEquals(expectedHasNext, iter.hasNext());
  }

  // Tests the static factory methods in IntIterator.
  @Test
  public void testIntIterators() {
    assertFalse(IntIterators.EMPTY.hasNext());

    OfInt single = IntIterators.of(9);
    assertNext(single, 9, false);

    OfInt triple = IntIterators.of(4, 5, 6);
    assertNext(triple, 4, true);
    assertNext(triple, 5, true);
    assertNext(triple, 6, false);

    OfInt wrappingBoxed = IntIterators.of(ImmutableList.of(9, 10, 11));
    assertNext(wrappingBoxed, 9, true);
    assertNext(wrappingBoxed, 10, true);
    assertNext(wrappingBoxed, 11, false);
  }

  // Tests default methods on IntSequence.
  @Test
  public void testIntSequenceDefaultMethods() {
    // Empty sequence.
    IntSequence empty = ImmutableIntSequence.of();
    assertTrue(empty.isEmpty());
    assertEquals(0, empty.size());
    empty.forEach(i -> fail("empty.forEach() should not be called"));
    empty.forEach((k, v) -> fail("empty.forEach() should not be called"));
    assertFalse(empty.contains(0));
    assertFalse(empty.stream().max().isPresent());
    int[] emptyArray = empty.toArray();
    assertEquals(0, emptyArray.length);
    assertEquals("IntSequence of 0 elements.", empty.debugString());

    // Single element sequence.
    IntSequence single = ImmutableIntSequence.of(9);
    assertFalse(single.isEmpty());
    assertEquals(1, single.size());
    single.forEach(i -> assertEquals(9, i));
    single.forEach((k, v) -> assertTrue(k == 0 && v == 9));
    assertTrue(single.contains(9));
    assertFalse(single.contains(0));
    assertEquals(9, single.stream().max().getAsInt());
    assertEquals(9, single.stream().min().getAsInt());
    int[] singleArray = single.toArray();
    assertEquals(1, singleArray.length);
    assertEquals(9, singleArray[0]);
    assertEquals("IntSequence of 1 elements: [9]", single.debugString());

    // Three element sequence.
    IntSequence triple = ImmutableIntSequence.of(4, 5, 6);
    assertFalse(triple.isEmpty());
    assertEquals(3, triple.size());
    triple.forEach(i -> assertTrue(i == 4 || i == 5 || i == 6));
    triple.forEach(
        (k, v) -> assertTrue((k == 0 && v == 4) || (k == 1 && v == 5) || (k == 2 && v == 6)));
    assertTrue(triple.contains(4));
    assertTrue(triple.contains(5));
    assertTrue(triple.contains(6));
    assertFalse(triple.contains(0));
    assertEquals(6, triple.stream().max().getAsInt());
    assertEquals(4, triple.stream().min().getAsInt());
    int[] tripleArray = triple.toArray();
    assertEquals(3, tripleArray.length);
    assertEquals(4, tripleArray[0]);
    assertEquals(5, tripleArray[1]);
    assertEquals(6, tripleArray[2]);
    assertEquals("IntSequence of 3 elements: [4, 5, 6]", triple.debugString());
  }

  // Tests the factory methods producing ImmutableIntSequences, and equals() implementation.
  @Test
  public void testImmutableIntSequenceEquals() {
    assertTrue(ImmutableIntSequence.of().equals(ImmutableIntSequence.of(ImmutableList.of())));
    assertTrue(ImmutableIntSequence.of().equals(IntSequence.of(ImmutableList.of())));
    assertTrue(ImmutableIntSequence.EMPTY.equals(ImmutableIntSequence.of()));

    assertFalse(ImmutableIntSequence.of(5).equals(ImmutableIntSequence.of()));
    assertFalse(ImmutableIntSequence.of(5).equals(ImmutableIntSequence.of(4)));
    assertTrue(ImmutableIntSequence.of(5).equals(ImmutableIntSequence.of(5)));

    assertFalse(ImmutableIntSequence.of(5, 6).equals(ImmutableIntSequence.of(5)));
    assertTrue(
        ImmutableIntSequence.of(3, 5).equals(ImmutableIntSequence.of(ImmutableList.of(3, 5))));
    assertTrue(ImmutableIntSequence.of(3, 5).equals(IntSequence.of(ImmutableList.of(3, 5))));
    assertFalse(ImmutableIntSequence.of(3, 5).equals(IntSequence.of(ImmutableList.of(3, 5, 7))));
    assertFalse(ImmutableIntSequence.of(3, 5, 7).equals(IntSequence.of(ImmutableList.of(3, 5))));
  }

  @Test
  public void testSequenceDebugString() {
    assertEquals("IntSequence of 0 elements.", ImmutableIntSequence.of().debugString());
    assertEquals("IntSequence of 1 elements: [3]", ImmutableIntSequence.of(3).debugString());
    assertEquals("IntSequence of 2 elements: [9, 6]", ImmutableIntSequence.of(9, 6).debugString());
    assertEquals(
        "IntSequence of 16 elements: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]",
        ImmutableIntSequence.of(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
            .debugString());
    assertEquals(
        "IntSequence of 17 elements: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16...]",
        ImmutableIntSequence.of(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
            .debugString());
  }

  @Test
  public void testSequenceHasher() {
    IntSequenceHasher hasher = new Ints.MixingHasher();

    assertEquals(
        hasher.hashCode(ImmutableIntSequence.of()), hasher.hashCode(ImmutableIntSequence.of()));
    assertEquals(
        hasher.hashCode(ImmutableIntSequence.of()), hasher.hashCode(ImmutableIntSequence.of()));

    assertFalse(
        hasher.hashCode(ImmutableIntSequence.of(5)) == hasher.hashCode(ImmutableIntSequence.of()));
    assertFalse(
        hasher.hashCode(ImmutableIntSequence.of(5)) == hasher.hashCode(ImmutableIntSequence.of(4)));
    assertEquals(
        hasher.hashCode(ImmutableIntSequence.of(5)), hasher.hashCode(ImmutableIntSequence.of(5)));

    assertFalse(
        hasher.hashCode(ImmutableIntSequence.of(5, 6))
            == hasher.hashCode(ImmutableIntSequence.of(5)));
    assertEquals(
        hasher.hashCode(ImmutableIntSequence.of(3, 5)),
        hasher.hashCode(ImmutableIntSequence.of(ImmutableList.of(3, 5))));
  }

  @Test
  public void testMixingHasherCollision() {
    IntSequenceHasher hasher = new Ints.MixingHasher();

    long hc1512 = hasher.hashCode(ImmutableIntSequence.of(1, 5, 1, 2));
    System.out.println("hc1512 = " + hc1512);

    long hc9367 = hasher.hashCode(ImmutableIntSequence.of(9, 3, 6, 7));
    System.out.println("hc9367 = " + hc9367);
    // This fails if signed rather than unsigned right shift is used in MixingHasher.mix().
    assertFalse(hc9367 == hc1512);
  }

  private IntSequence generateRandomSequence(int length) {
    int[] a = new int[length];
    for (int i = 0; i < length; i++) {
      a[i] = rand.nextInt();
    }
    return ImmutableIntSequence.of(a);
  }

  /**
   * Generate 100,000 random sequences of 10 random positive integers and hashes them with the
   * MixingHasher. Check that exactly one hash collision occurs. This depends on the random seed and
   * the implementation of Random.
   */
  @Test
  public void testRarityOfMixingHasherHashCollisions() {
    rand.setSeed(SEED);
    HashSet<Integer> seen = new HashSet<>();
    int collisions = 0;
    for (int i = 0; i < 100000; i++) {
      int hash = Ints.MIXING_HASHER.hashCode(generateRandomSequence(10));
      collisions += seen.add(hash) ? 0 : 1;
    }
    assertEquals(1, collisions);
  }
}
