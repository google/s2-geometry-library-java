/*
 * Copyright 2013 Google Inc.
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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Comparator;
import org.jspecify.annotations.Nullable;

/**
 * Contains test constants for the geometry package that require different values for the Java
 * client and Javascript client versions.
 *
 * <p>This contains the Java client version.
 */
class TestPlatform {

  /** Number of random cells to check neighbors for in S2CellIdTest#testNeighbors. */
  static final int S2_CELL_ID_TEST_RANDOM_NEIGHBORS_CHECKS = 1000;

  /**
   * Maximum amount to increase angle by in S2LoopTest.testAreaCentroid. This needs to be different
   * for the Java and Javascript test versions since this test does a bunch of math, and the
   * Javascript version is so much slower than the native Java version.
   */
  static final double S2_LOOP_TEST_MAX_ANGLE_DELTA = 1e-6;

  /** Number of iterations to run in S2LoopTest.testLoopRelations2. */
  static final int S2_LOOP_TEST_LOOP_RELATIONS2_ITERATIONS = 1000;

  /** Number of S2ShapeIndex edges to validate quadratically. */
  static final int S2_SHAPEINDEX_NUM_EDGES = 100;

  /**
   * Tests Java serialization of the object. If comparator is null, uses the object's equals method.
   */
  @SuppressWarnings("unchecked")
  static <T> void testSerialization(T obj, @Nullable Comparator<T> comparator)
      throws IOException, ClassNotFoundException {
    ByteArrayOutputStream byteStream = new ByteArrayOutputStream();
    ObjectOutputStream objectOutputStream = new ObjectOutputStream(byteStream);

    objectOutputStream.writeObject(obj);

    ObjectInputStream objectInputStream =
        new ObjectInputStream(new ByteArrayInputStream(byteStream.toByteArray()));

    T readObj = (T) objectInputStream.readObject();

    if (comparator != null) {
      assertEquals("Expected " + obj + ", got " + readObj, 0, comparator.compare(obj, readObj));
    } else {
      assertEquals(obj, readObj);
    }
  }

  /**
   * Asserts that two doubles are identical, treating different NaN representations as not
   * identical. Web cannot identify different NaN representations so super-sourced there with a
   * fallback implementation.
   */
  static void assertIdentical(double expected, double actual) {
    assertEquals(Double.doubleToRawLongBits(expected), Double.doubleToRawLongBits(expected));
  }
}
