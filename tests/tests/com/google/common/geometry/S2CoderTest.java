/*
 * Copyright 2023 Google Inc.
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

import com.google.common.geometry.PrimitiveArrays.Bytes;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link S2Coder}. */
@RunWith(JUnit4.class)
public class S2CoderTest extends GeometryTestCase {
  @Test
  public void testVarint() throws IOException {
    Long[] tests = {Long.MIN_VALUE, -1L, 0L, 1L, 4L, Long.MAX_VALUE};
    for (Long test : tests) {
      assertEquals(test, encodeDecode(S2Coder.UNBOXED_VARINT, test));
    }
    assertThrows(NullPointerException.class, () -> encodeDecode(S2Coder.UNBOXED_VARINT, null));
  }

  @Test
  public void testString() throws IOException {
    String[] tests = {"", "foo", "ăѣծềſģȟ"};
    for (String test : tests) {
      assertEquals(test, encodeDecode(S2Coder.STRING, test));
    }
  }

  private static <T> T encodeDecode(S2Coder<T> coder, T value) throws IOException {
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    coder.encode(value, bytes);
    return coder.decode(Bytes.fromByteArray(bytes.toByteArray()));
  }
}
