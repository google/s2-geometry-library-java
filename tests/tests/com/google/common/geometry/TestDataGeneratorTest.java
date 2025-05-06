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
package com.google.common.geometry;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class TestDataGeneratorTest {

  private final TestDataGenerator tdg = new TestDataGenerator();

  @Test
  public void testGetRandomCellIdAtLevel() {
    S2CellId l20 = tdg.getRandomCellId(20);
    // child level must be greater or equal to parent level.
    assertThrows(IllegalArgumentException.class, () -> tdg.getRandomCellId(l20, 19));

    checkRandomChildCellIds(0, 0);
    checkRandomChildCellIds(0, 1);
    checkRandomChildCellIds(0, S2CellId.MAX_LEVEL - 1);
    checkRandomChildCellIds(0, S2CellId.MAX_LEVEL);

    checkRandomChildCellIds(20, 20);
    checkRandomChildCellIds(20, 21);
    checkRandomChildCellIds(20, S2CellId.MAX_LEVEL - 1);
    checkRandomChildCellIds(20, S2CellId.MAX_LEVEL);

    checkRandomChildCellIds(S2CellId.MAX_LEVEL, S2CellId.MAX_LEVEL);
  }

  /** Check correctness of a few random cell ids for the given parent and child levels. */
  private void checkRandomChildCellIds(int parentLevel, int childLevel) {
    S2CellId parent = tdg.getRandomCellId(parentLevel);
    for (int i = 0; i < 10; i++) {
      S2CellId child = tdg.getRandomCellId(parent, childLevel);
      assertEquals(child.level(), childLevel);
      assertTrue(parent.contains(child));
    }
  }
}
