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
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;

/** Unit tests for {@link LittleEndianOutput}. */
public class LittleEndianOutputTest extends TestCase {
  private ByteArrayOutputStream bytes;
  private LittleEndianOutput output;

  @Override
  protected void setUp() {
    bytes = new ByteArrayOutputStream();
    output = new LittleEndianOutput(bytes);
  }

  public void testByte() throws IOException {
    output.writeByte((byte) 0x00);
    output.writeByte((byte) 0xFF);
    expect(0, -1);
  }

  public void testBytes() throws IOException {
    byte[][] input = {Bytes.toArray(Ints.asList(0, 1, 2)), Bytes.toArray(Ints.asList(3, 4, 5))};
    output.writeBytes(input[0]);
    output.writeBytes(input[1]);
    expect(0, 1, 2, 3, 4, 5);
  }

  public void testInt() throws IOException {
    output.writeInt(0);
    output.writeInt(-1);
    expect(0, 0, 0, 0, -1, -1, -1, -1);
  }

  public void testLong() throws IOException {
    output.writeLong(0);
    output.writeLong(-1);
    expect(0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1);
  }

  public void testFloat() throws IOException {
    output.writeFloat(0F);
    output.writeFloat(Float.MIN_VALUE);
    output.writeFloat(Float.NaN);
    output.writeFloat(Float.POSITIVE_INFINITY);
    expect(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -64, 127, 0, 0, -128, 127);
  }

  public void testDouble() throws IOException {
    output.writeDouble(0D);
    output.writeDouble(Double.MIN_VALUE);
    output.writeDouble(Double.NaN);
    output.writeDouble(Double.POSITIVE_INFINITY);
    expect(
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8, 127, 0, 0, 0, 0, 0, 0,
        -16, 127);
  }

  public void testVarint32() throws IOException {
    output.writeVarint32(0);
    output.writeVarint32(1);
    output.writeVarint32(127);
    // 14882
    output.writeVarint32((0x22 << 0) | (0x74 << 7));
    expect(0, 1, 127, 0xa2, 0x74);
  }

  public void testVarint64() throws IOException {
    // 2961488830
    output.writeVarint64((0x3e << 0) | (0x77 << 7) | (0x12 << 14) | (0x04 << 21) | (0x0bL << 28));
    // 7256456126
    output.writeVarint64((0x3e << 0) | (0x77 << 7) | (0x12 << 14) | (0x04 << 21) | (0x1bL << 28));
    // 41256202580718336
    output.writeVarint64(
        (0x00 << 0)
            | (0x66 << 7)
            | (0x6b << 14)
            | (0x1c << 21)
            | (0x43L << 28)
            | (0x49L << 35)
            | (0x24L << 42)
            | (0x49L << 49));
    // 11964378330978735131
    output.writeVarint64(
        (0x1b << 0)
            | (0x28 << 7)
            | (0x79 << 14)
            | (0x42 << 21)
            | (0x3bL << 28)
            | (0x56L << 35)
            | (0x00L << 42)
            | (0x05L << 49)
            | (0x26L << 56)
            | (0x01L << 63));
    expect(
        0xbe, 0xf7, 0x92, 0x84, 0x0b, 0xbe, 0xf7, 0x92, 0x84, 0x1b, 0x80, 0xe6, 0xeb, 0x9c, 0xc3,
        0xc9, 0xa4, 0x49, 0x9b, 0xa8, 0xf9, 0xc2, 0xbb, 0xd6, 0x80, 0x85, 0xa6, 0x01);
  }

  private void expect(int... expected) {
    byte[] actual = bytes.toByteArray();
    assertTrue(
        Arrays.toString(actual), Arrays.equals(Bytes.toArray(Ints.asList(expected)), actual));
  }
}
