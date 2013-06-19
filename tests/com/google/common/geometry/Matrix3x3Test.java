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

import com.google.common.annotations.GwtCompatible;

import junit.framework.TestCase;

/**
 * Test case for a simple 3x3 matrix class.
 */
@GwtCompatible
public class Matrix3x3Test extends TestCase {
  public void testCtor() {
    Matrix3x3 m = new Matrix3x3(2, 1, 2, 3, 4, 5, 6);
    assertEquals(2, m.cols());
    assertEquals(3, m.rows());
    assertEquals(1D, m.get(0, 0));
    assertEquals(6D, m.get(2, 1));
    assertEquals(4D, m.get(1, 1));
    m.set(1, 1, 1);
    assertEquals(1D, m.get(1, 1));
  }

  public void testMatrixTranspose() {
    Matrix3x3 m = new Matrix3x3(2, 1, 2, 3, 4);
    assertEquals(new Matrix3x3(2, 1, 3, 2, 4), m.transpose());
  }

  public void testMatrixMult() {
    Matrix3x3 a = new Matrix3x3(3, 1, 2, 3, 4, 5, 6);
    Matrix3x3 b = new Matrix3x3(1, 3, 2, 1);
    Matrix3x3 result = new Matrix3x3(1, 10, 28);
    assertEquals(result, a.mult(b));
  }
}
