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

import junit.framework.TestCase;

/**
 * Verifies R2Vector.
 */
public class R2VectorTest extends TestCase {
  public void testOrtho() {
    assertEquals(new R2Vector(1, 1), new R2Vector(1, -1).ortho());
    assertEquals(new R2Vector(1, -1), new R2Vector(-1, -1).ortho());
    assertEquals(new R2Vector(-1, -1), new R2Vector(-1, 1).ortho());
    assertEquals(new R2Vector(1, 1), new R2Vector(1, -1).ortho());
  }

  public void testSub() {
    assertEquals(new R2Vector(3, 1), R2Vector.sub(new R2Vector(4, 3), new R2Vector(1, 2)));
  }
}
