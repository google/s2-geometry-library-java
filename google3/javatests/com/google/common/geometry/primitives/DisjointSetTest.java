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
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.GeometryTestCase;
import com.google.common.geometry.S2Point;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link DisjointSet}. */
@RunWith(JUnit4.class)
public class DisjointSetTest extends GeometryTestCase {

  @Test
  public void testS2PointSetCompiles() {
    DisjointSet<S2Point> set = new DisjointSet<>();
    set.add(S2Point.X_POS);
    S2Point root = set.findRoot(S2Point.X_POS);
    assertNotNull(root);
    assertEquals(S2Point.X_POS, root);
    set.clear();
    assertNull(set.findRoot(S2Point.X_POS));
  }

  @Test
  public void testInsertMoreThanOnceFails() {
    DisjointSet<Integer> set = new DisjointSet<>();
    assertTrue(set.add(1));
    assertFalse(set.add(1));
    assertFalse(set.add(1));
  }

  @Test
  public void testFindRootWorks() {
    DisjointSet<Integer> set = new DisjointSet<>();
    set.add(1);
    Integer root = set.findRoot(1);
    assertNotNull(root);
    assertEquals(1, (int) root);
    assertNull(set.findRoot(2));
  }

  @Test
  public void testUnionWorks() {
    DisjointSet<Integer> set = new DisjointSet<>();
    for (int i = 0; i < 10; ++i) {
      assertTrue(set.add(i));
    }

    // Union into two disjoint sets
    for (int i = 0; i < 4; ++i) {
      assertTrue(set.union(i + 0, i + 1));
      assertTrue(set.union(i + 5, i + 6));
    }

    // The values in the two sets should all have the same root now.
    for (int i = 0; i < 5; ++i) {
      Integer root = set.findRoot(i);
      assertNotNull(root);
      assertEquals(0, (int) root);

      root = set.findRoot(i + 5);
      assertNotNull(root);
      assertEquals(5, (int) root);
    }

    // Check that attempting to union with something not in the set fails.
    assertFalse(set.union(0, 13));
    assertFalse(set.union(13, 0));
    assertFalse(set.union(12, 13));

    // Union of two elements from two sets should unify the sets
    assertTrue(set.union(3, 7));
    for (int i = 0; i < 10; ++i) {
      Integer root = set.findRoot(i);
      assertNotNull(root);
      assertEquals(0, (int) root);
    }

    // Union of two elements already in the same set should return true.
    assertTrue(set.union(3, 7));
  }

  @Test
  public void testSizeAndClearWorks() {
    DisjointSet<Integer> set = new DisjointSet<>(10);
    for (int i = 0; i < 10; ++i) {
      assertTrue(set.add(i));
    }

    assertEquals(10, set.size());
    for (int i = 0; i < set.size() - 1; ++i) {
      assertTrue(set.union(i, i + 1));
    }

    // Unioning doesn't change size.
    assertEquals(10, set.size());
    set.clear();
    assertEquals(0, set.size());
  }
}
