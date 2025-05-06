/*
 * Copyright 2018 Google Inc.
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

import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.io.BaseEncoding;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link VectorCoder}. */
@RunWith(JUnit4.class)
public class VectorCoderTest {

  @Test
  public void testDecodeFromByteString() throws IOException {
    List<String> expected = Lists.newArrayList("fuji", "mutsu");

    byte[] b = BaseEncoding.base16().decode("00000010040966756A696D75747375");
    PrimitiveArrays.Bytes data = PrimitiveArrays.Bytes.fromByteArray(b);
    int offset = 3;
    List<String> v = VectorCoder.STRING.decode(data, data.cursor(offset));
    assertEquals(expected, v);
  }

  @Test
  public void testEmpty() throws IOException {
    checkEncodedStringVector(Lists.newArrayList(), 1);
  }

  @Test
  public void testEmptyString() throws IOException {
    checkEncodedStringVector(Lists.newArrayList(""), 2);
  }

  @Test
  public void testRepeatedEmptyStrings() throws IOException {
    checkEncodedStringVector(Lists.newArrayList("", "", ""), 4);
  }

  @Test
  public void testOneString() throws IOException {
    checkEncodedStringVector(Lists.newArrayList("apples"), 8);
  }

  @Test
  public void testTwoStrings() throws IOException {
    checkEncodedStringVector(Lists.newArrayList("fuji", "mustu"), 12);
  }

  @Test
  public void testTwoBigStrings() throws IOException {
    checkEncodedStringVector(Lists.newArrayList("x".repeat(10000), "y".repeat(100000)), 110007);
  }

  private void checkEncodedStringVector(List<String> input, int expectedBytes) throws IOException {
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    VectorCoder.STRING.encode(input, output);
    assertEquals(expectedBytes, output.size());

    PrimitiveArrays.Bytes data = PrimitiveArrays.Bytes.fromByteArray(output.toByteArray());
    Cursor cursor = data.cursor();
    List<String> actual = VectorCoder.STRING.decode(data, cursor);
    assertEquals(input, actual);
    assertEquals(expectedBytes, cursor.position);
  }
}
