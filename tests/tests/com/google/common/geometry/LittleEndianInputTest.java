/*
 * Copyright 2016 Google Inc.
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

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;

/** Unit tests for {@link LittleEndianInput}. */
public class LittleEndianInputTest extends TestCase {
  public void testByte() throws IOException {
    LittleEndianInput input = input(0, -1);
    assertEquals(0, input.readByte());
    assertEquals(-1, input.readByte());
    expectEOF(input);
  }

  public void testBytes() throws IOException {
    LittleEndianInput input = input(0, 1, 2, 3, 4, 5);
    byte[][] expected = {Bytes.toArray(Ints.asList(0, 1, 2)), Bytes.toArray(Ints.asList(3, 4, 5))};
    assertTrue(Arrays.toString(expected[0]), Arrays.equals(expected[0], input.readBytes(3)));
    assertTrue(Arrays.toString(expected[1]), Arrays.equals(expected[1], input.readBytes(3)));
    expectEOF(input);
  }

  public void testInt() throws IOException {
    LittleEndianInput input = input(0, 0, 0, 0, -1, -1, -1, -1);
    assertEquals(0, input.readInt());
    assertEquals(-1, input.readInt());
    expectEOF(input);
  }

  public void testLong() throws IOException {
    LittleEndianInput input = input(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);
    assertEquals(0L, input.readLong());
    assertEquals(-1L, input.readLong());
    expectEOF(input);
  }

  public void testFloat() throws IOException {
    LittleEndianInput input = input(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -64, 127, 0, 0, 0x80, 0x7F);
    assertEquals(0F, input.readFloat(), 0F);
    assertEquals(Float.MIN_VALUE, input.readFloat(), 0F);
    assertEquals(0, Float.compare(Float.NaN, input.readFloat()), 0F);
    assertEquals(Float.POSITIVE_INFINITY, input.readFloat(), 0F);
    expectEOF(input);
  }

  public void testDouble() throws IOException {
    LittleEndianInput input =
        input(
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8, 127, 0, 0, 0, 0,
            0, 0, -16, 127);
    assertEquals(0D, input.readDouble(), 0D);
    assertEquals(Double.MIN_VALUE, input.readDouble(), 0D);
    assertEquals(0, Double.compare(Double.NaN, input.readDouble()));
    assertEquals(Double.POSITIVE_INFINITY, input.readDouble(), 0D);
    expectEOF(input);
  }

  public void testVarint32() throws IOException {
    LittleEndianInput input = input(0, 1, 127, 0xa2, 0x74);
    assertEquals(0, input.readVarint32());
    assertEquals(1, input.readVarint32());
    assertEquals(127, input.readVarint32());
    // 14882
    assertEquals((0x22 << 0) | (0x74 << 7), input.readVarint32());
    expectEOF(input);
  }

  public void testVarint64() throws IOException {
    LittleEndianInput input =
        input(
            0xbe, 0xf7, 0x92, 0x84, 0x0b, 0xbe, 0xf7, 0x92, 0x84, 0x1b, 0x80, 0xe6, 0xeb, 0x9c,
            0xc3, 0xc9, 0xa4, 0x49, 0x9b, 0xa8, 0xf9, 0xc2, 0xbb, 0xd6, 0x80, 0x85, 0xa6, 0x01);
    // 2961488830
    assertEquals(
        (0x3e << 0) | (0x77 << 7) | (0x12 << 14) | (0x04 << 21) | (0x0bL << 28),
        input.readVarint64());
    // 7256456126
    assertEquals(
        (0x3e << 0) | (0x77 << 7) | (0x12 << 14) | (0x04 << 21) | (0x1bL << 28),
        input.readVarint64());
    // 41256202580718336
    assertEquals(
        (0x00 << 0)
            | (0x66 << 7)
            | (0x6b << 14)
            | (0x1c << 21)
            | (0x43L << 28)
            | (0x49L << 35)
            | (0x24L << 42)
            | (0x49L << 49),
        input.readVarint64());
    // 11964378330978735131
    assertEquals(
        (0x1b << 0)
            | (0x28 << 7)
            | (0x79 << 14)
            | (0x42 << 21)
            | (0x3bL << 28)
            | (0x56L << 35)
            | (0x00L << 42)
            | (0x05L << 49)
            | (0x26L << 56)
            | (0x01L << 63),
        input.readVarint64());
    expectEOF(input);
  }

  private LittleEndianInput input(int... bytes) {
    return new LittleEndianInput(new ByteArrayInputStream(Bytes.toArray(Ints.asList(bytes))));
  }

  private void expectEOF(LittleEndianInput input) {
    try {
      byte unused = input.readByte();
      fail();
    } catch (IOException eof) {
      // Expected.
    }
  }
}
