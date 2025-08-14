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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Comparator;
import java.util.Iterator;
import java.util.NoSuchElementException;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link PooledList}. */
@RunWith(JUnit4.class)
public final class PooledListTest {

  /** A mutable String container for testing PooledList. Counts allocations. */
  @SuppressWarnings("NonFinalStaticField")
  private static class MutableString {
    public static int totalAllocations = 0;

    public String value;
    public MutableString() {
      totalAllocations++;
      value = "";
    }

    public void set(String value) {
      this.value = value;
    }

    public String get() {
      return value;
    }
  }

  @Test
  public void testBasics() {
    PooledList<MutableString> list = new PooledList<>(MutableString::new);
    assertEquals(0, list.size());
    assertTrue(list.isEmpty());

    MutableString m1 = list.add();
    m1.set("foo");
    assertEquals(1, list.size());
    assertFalse(list.isEmpty());
    assertEquals(1, MutableString.totalAllocations);

    MutableString m2 = list.add();
    m2.set("bar");
    assertEquals(2, list.size());
    assertEquals(2, MutableString.totalAllocations);

    assertEquals("foo", list.get(0).get());
    assertEquals("bar", list.get(1).get());

    // Clear the list, which should retain the allocated objects in the pool.
    list.clear();
    assertEquals(0, list.size());

    // Verify that objects are reused.
    MutableString m3 = list.add();
    // Reused objects are not cleared.
    assertEquals("foo", m3.get());
    m3.set("foo2");
    // Objects from the pool don't require allocation.
    assertEquals(2, MutableString.totalAllocations);

    MutableString m4 = list.add();
    assertEquals("bar", m4.get());
    m4.set("bar2");
    assertEquals(2, MutableString.totalAllocations);

    // Get another object, which will require allocation.
    MutableString m5 = list.add();
    assertEquals("", m5.get());
    m5.set("new");
    assertEquals(3, list.size());
    assertEquals(3, MutableString.totalAllocations);

    // Test resize to a larger size.
    list.resize(5);
    assertEquals(5, list.size());
    assertEquals(5, MutableString.totalAllocations);

    // Check list contents.
    assertEquals("foo2", list.get(0).get());
    assertEquals("bar2", list.get(1).get());
    assertEquals("new", list.get(2).get());
    assertEquals("", list.get(3).get());
    assertEquals("", list.get(4).get());

    // Verify that access beyond the end of the list throws an exception.
    try {
      list.get(5);
      fail("Expected IndexOutOfBoundsException");
    } catch (IndexOutOfBoundsException expected) {
    }

    // Resize to a smaller size.
    list.resize(2);
    assertEquals(2, list.size());

    // Verify swap().
    list.swap(0, 1);
    assertEquals("bar2", list.get(0).get());
    assertEquals("foo2", list.get(1).get());

    // Verify that access beyond the end of the list throws an exception, even when the underlying
    // array is large enough.
    try {
      list.get(3);
      fail("Expected IndexOutOfBoundsException");
    } catch (IndexOutOfBoundsException expected) {
    }

    // Resize back to a larger size, without requiring allocation.
    list.resize(5);
    assertEquals(5, list.size());
    assertEquals(5, MutableString.totalAllocations);
    list.get(4).set("last");

    // Verify iteration.
    Iterator<MutableString> it = list.iterator();
    // We know there are 5 elements.
    for (int i = 0; i < 5; i++) {
      assertTrue(it.hasNext());
      MutableString m = it.next();
      switch (i) {
        case 0:
          assertEquals("bar2", m.get());
          break;
        case 1:
          assertEquals("foo2", m.get());
          break;
        case 2:
          assertEquals("new", m.get());
          break;
        case 3:
          assertEquals("", m.get());
          break;
        case 4:
          assertEquals("last", m.get());
          break;
        default:
          fail("Unexpected iteration index: " + i);
      }
    }
    // The iterator should have been exhausted.
    assertFalse(it.hasNext());
    // Calling next() should throw an exception.
    try {
      it.next();
      fail("Expected NoSuchElementException");
    } catch (NoSuchElementException expected) {
    }

    // Sort the list alphabetically by comparing the Strings returned by the get() method.
    list.sort(Comparator.comparing(MutableString::get));
    assertEquals("", list.get(0).get());
    assertEquals("bar2", list.get(1).get());
    assertEquals("foo2", list.get(2).get());
    assertEquals("last", list.get(3).get());
    assertEquals("new", list.get(4).get());

    assertEquals(5, MutableString.totalAllocations);
  }

  @Test
  public void testRemove() {
    PooledList<MutableString> list = new PooledList<>(MutableString::new);
    list.add().set("zero");
    list.add().set("one");
    list.add().set("two");
    list.add().set("three");
    list.add().set("four");
    list.add().set("five");
    assertEquals(6, list.size());
    assertEquals("five", list.get(5).get());

    list.remove(2);
    assertEquals(5, list.size());
    assertEquals("five", list.get(4).get());

    Assert.assertThrows(IndexOutOfBoundsException.class, () -> list.remove(5));

    // Removing an element puts it back in the pool, so now adding another element won't require
    // allocation.
    int allocations = MutableString.totalAllocations;
    list.add().set("reused");
    assertEquals(allocations, MutableString.totalAllocations);
  }
}
