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

import static com.google.common.geometry.S2CellRangeIterator.makeS2CellRangeIterator;
import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.CellRelation;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class S2CellRangeIteratorTest {

  @Test
  public void testRelation() {
    // Create an index with one point each on S2CellId faces 0, 1, and 2.
    S2ShapeIndex index = makeIndexOrDie("0:0 | 0:90 | 90:0 # #");
    S2CellRangeIterator<Cell> it0 = makeS2CellRangeIterator(index);
    S2CellRangeIterator<Cell> it1 = makeS2CellRangeIterator(index);
    it1.next();
    assertEquals(-1, it0.relation(it1));
    assertEquals(+1, it1.relation(it0));
    it1.prev();
    assertEquals(0, it1.relation(it0));
    it1.finish();
    assertEquals(+1, it1.relation(it0));
  }

  @Test
  public void testNext() {
    // Create an index with one point each on S2CellId faces 0, 1, and 2.
    S2ShapeIndex index = makeIndexOrDie("0:0 | 0:90 | 90:0 # #");
    S2CellRangeIterator<Cell> it = makeS2CellRangeIterator(index);
    assertEquals(0, it.id().face());
    it.next();
    assertEquals(1, it.id().face());
    it.next();
    assertEquals(2, it.id().face());
    it.next();
    // S2Iterator.id() throws an assertion if it is done. If assertions are disabled, it returns a
    // sentinel value.
    Assert.assertThrows(AssertionError.class, it::id);
    assertTrue(it.done());
  }

  @Test
  public void testLocate() {
    // Create two indices with one point each on S2CellId faces 0, 1, and 2.
    S2ShapeIndex index0 = makeIndexOrDie("0:0 | 0:90 | 90:0 # #");
    S2ShapeIndex index1 = makeIndexOrDie("0:0 | 0:90 | 90:0 # #");
    S2CellRangeIterator<Cell> it0 = makeS2CellRangeIterator(index0);
    S2CellRangeIterator<Cell> it1 = makeS2CellRangeIterator(index1);
    it0.next();
    CellRelation unused = it1.locate(it0);
    assertEquals(it1.id(), it0.id());
  }

  @Test
  public void testEmptyIndex() {
    S2ShapeIndex empty = makeIndexOrDie("# #");
    S2ShapeIndex nonEmpty = makeIndexOrDie("0:0 # #");
    S2CellRangeIterator<Cell> emptyIt = makeS2CellRangeIterator(empty);
    S2CellRangeIterator<Cell> nonEmptyIt = makeS2CellRangeIterator(nonEmpty);
    assertFalse(nonEmptyIt.done());
    assertTrue(emptyIt.done());

    emptyIt.seekTo(nonEmptyIt);
    assertTrue(emptyIt.done());

    emptyIt.seekBeyond(nonEmptyIt);
    assertTrue(emptyIt.done());

    emptyIt.seekTo(emptyIt);
    assertTrue(emptyIt.done());

    emptyIt.seekBeyond(emptyIt);
    assertTrue(emptyIt.done());
  }
}
