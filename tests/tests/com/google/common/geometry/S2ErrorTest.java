/*
 * Copyright 2014 Google Inc.
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

/** Verifies S2Error. */
public class S2ErrorTest extends TestCase {
  public void testBasic() {
    S2Error error = new S2Error();
    error.init(S2Error.Code.DUPLICATE_VERTICES, "Vertex %d is the same as vertex %d", 23, 47);
    // Prepend additional context to the message.
    error.init(error.code(), "Loop %d: %s", 5, error.text());
    assertEquals(S2Error.Code.DUPLICATE_VERTICES, error.code());
    assertEquals("Loop 5: Vertex 23 is the same as vertex 47", error.text());
  }
}
