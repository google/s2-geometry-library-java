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

import static com.google.common.geometry.GeometryTestCase.assertExactly;

import junit.framework.TestCase;

/** Unit tests for the {@link Matrix} class. */
public class MatrixTest extends TestCase {
  public void testCtor() {
    Matrix m = new Matrix(2, 1, 2, 3, 4, 5, 6);
    assertEquals(2, m.cols());
    assertEquals(3, m.rows());
    assertExactly(1D, m.get(0, 0));
    assertExactly(6D, m.get(2, 1));
    assertExactly(4D, m.get(1, 1));
    m.set(1, 1, 1);
    assertExactly(1D, m.get(1, 1));
  }

  public void testMatrixTranspose() {
    Matrix m = new Matrix(2, 1, 2, 3, 4);
    assertEquals(new Matrix(2, 1, 3, 2, 4), m.transpose());
  }

  public void testMatrixMult() {
    Matrix a = new Matrix(3, 1, 2, 3, 4, 5, 6);
    Matrix b = new Matrix(1, 3, 2, 1);
    Matrix result = new Matrix(1, 10, 28);
    assertEquals(result, a.mult(b));
  }
}
