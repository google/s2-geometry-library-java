/*
 * Copyright 2011 Google Inc.
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

package com.google.testing.util;

import static java.util.Arrays.asList;

import com.google.common.base.Objects;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import java.util.Comparator;

import junit.framework.Assert;

import java.util.Iterator;
import java.util.List;

public final class MoreAsserts {

  private MoreAsserts() { }

  /**
   * Asserts that {@code actual} contains precisely the elements
   * {@code expected}, in any order.  Both collections may contain
   * duplicates, and this method will only pass if the quantities are
   * exactly the same.
   */
  public static void assertContentsAnyOrder(
      String message, Iterable<?> actual, Object... expected) {
    assertEqualsImpl(message,
        HashMultiset.create(asList(expected)), HashMultiset.create(actual));
  }

  /**
   * Variant of {@link #assertContentsAnyOrder(String,Iterable,Object...)}
   * using a generic message.
   */
  public static void assertContentsAnyOrder(
      Iterable<?> actual, Object... expected) {
    assertContentsAnyOrder((String) null, actual, expected);
  }

  /**
   * Asserts that {@code actual} contains precisely the elements
   * {@code expected}, in any order.  Both collections may contain
   * duplicates, and this method will only pass if the quantities are
   * exactly the same.  This method uses the user-provided Comparator
   * object for doing the object comparison, instead of relying on the
   * contents' implementation of {@link Object#equals(Object)}. It also takes
   * in the expected set of objects as an Iterable.
   * <p>
   * Note the different order of expected and actual from the other
   * {@link #assertContentsAnyOrder(String,Iterable,Object...)}
   */
  public static <T> void assertContentsAnyOrder(String message,
      Iterable<? extends T> expected, Iterable<? extends T> actual,
      Comparator<? super T> comparator) {
    // We should not iterate over an Iterable more than once. There's
    // no guarentees that Iterable.iterator() returns an iterator over the
    // entire collection every time.
    //
    // Why don't we use TreeMultiset? Unfortunately, TreeMultiset.toString()
    // produces really odd output for duplicates. In addition, our contract
    // states that we use the comparator to compare equality, not to order
    // items.
    ImmutableList<T> actualList = ImmutableList.copyOf(actual);
    ImmutableList<T> expectedList = ImmutableList.copyOf(expected);

    // First compare sizes to save ourselves on N X M operation.
    // This also handles the case where "expected" is a subset of "actual".
    if (actualList.size() != expectedList.size()) {
      failNotEqual(message, expectedList, actualList);
    }

    // Now for each expected value, iterate through actuals and delete entry
    // if found. We need to make another copy of the "actual" items because
    // we will be removing items from this list, and we need to keep the original
    // for the failure message.
    List<T> unfoundItems = Lists.newLinkedList(actualList);
    for (T ex : expectedList) {
      boolean found = false;
      Iterator<T> iter = unfoundItems.iterator();
      while (iter.hasNext()) {
        T ac = iter.next();
        if (comparator.compare(ex, ac) == 0) {
          iter.remove();
          found = true;
          break;
        }
      }
      if (!found) {
        failNotEqual(message, expectedList, actualList);
      }
    }
  }

  /**
   * Variant of {@link #assertContentsAnyOrder(String,Iterable,Object...)}
   * using a generic message.
   */
  public static <T> void assertContentsAnyOrder(
      Iterable<? extends T> expected, Iterable<? extends T> actual,
      Comparator<? super T> comparator) {
    assertContentsAnyOrder((String) null, expected, actual, comparator);
  }

  private static void failNotEqual(String message, Object expected,
      Object actual) {
    if ((expected != null) && (actual != null)
        && expected.toString().equals(actual.toString())) {
      failWithMessage(message, "expected:<("
          + expected.getClass().getName() + ") " + expected + "> but was:<("
          + actual.getClass().getName() + ") " + actual + ">");
    } else {
      failWithMessage(message, "expected:<" + expected + "> but was:<" + actual
          + ">");
    }
  }

  /**
   * Replacement of {@link Assert#assertEquals} which provides the same error
   * message in GWT and java.
   */
  private static void assertEqualsImpl(
      String message, Object expected, Object actual) {
    if (!Objects.equal(expected, actual)) {
      failWithMessage(
          message, "expected:<" + expected + "> but was:<" + actual + ">");
    }
  }

  private static void failWithMessage(String userMessage, String ourMessage) {
    Assert.fail((userMessage == null)
        ? ourMessage
        : userMessage + ' ' + ourMessage);
  }
}
