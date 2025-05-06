/*
 * Copyright 2024 Google Inc.
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

import junit.framework.TestCase;

/** Tests for {@link CountingBitset}. */
public class CountingBitsetTest extends TestCase {

  public void testCountingBitsetWorks() {
    int numBits = 113;

    CountingBitset bitmap = new CountingBitset(numBits);
    assertEquals(numBits, bitmap.bits());
    assertEquals(0, bitmap.numSet());

    for (int i = 0; i < numBits - 1; ++i) {
      bitmap.set(i);
      assertFalse(bitmap.full());
      assertEquals(i + 1, bitmap.numSet());
    }
    bitmap.set(numBits - 1);
    assertTrue(bitmap.full());

    bitmap.clear();
    assertEquals(0, bitmap.numSet());

    // Test that reset() sets size as expected.
    bitmap.reset(128);
    assertEquals(128, bitmap.bits());
  }

  public void testCountingBitsetIsIdempotent() {
    int numBits = 113;
    CountingBitset bitmap = new CountingBitset(numBits);
    assertEquals(numBits, bitmap.bits());
    assertEquals(0, bitmap.numSet());

    bitmap.set(37);
    assertTrue(bitmap.get(37));
    assertEquals(1, bitmap.numSet());

    // Setting again shouldn't change anything.
    bitmap.set(37);
    assertTrue(bitmap.get(37));
    assertEquals(1, bitmap.numSet());

    // Same with clearing.
    assertFalse(bitmap.get(52));
    assertEquals(1, bitmap.numSet());

    bitmap.set(52, false);
    assertFalse(bitmap.get(52));
    assertFalse(bitmap.empty());
    assertEquals(1, bitmap.numSet());

    // But clearing our one set bit should decrement
    bitmap.set(37, false);
    assertFalse(bitmap.get(37));
    assertTrue(bitmap.empty());
    assertEquals(0, bitmap.numSet());
  }
}
