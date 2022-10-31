/*
 * Copyright 2019 Google Inc.
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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.primitives.ImmutableLongArray;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;

/** Unit tests for {@link PrimitiveArrays}. */
public class PrimitiveArraysTest extends GeometryTestCase {
  @GwtIncompatible("Uses ByteBuffer")
  public void testBytesFromByteBuffer() {
    byte[] b = new byte[] {0, 1, 2};
    Bytes data = Bytes.fromByteBuffer(ByteBuffer.wrap(b));

    assertEquals(3, b.length);
    for (int i = 0; i < b.length; i++) {
      assertEquals(b[i], data.get(i));
    }
  }

  public void testBytesFromByteArray() {
    byte[] b = new byte[] {0, 1, 2};
    Bytes data = Bytes.fromByteArray(b);

    assertEquals(3, b.length);
    for (int i = 0; i < b.length; i++) {
      assertEquals(b[i], data.get(i));
    }
  }

  public void testBytesToInputStream() throws IOException {
    byte[] b = new byte[] {0, 1, 2, 3, 4, 5};
    Bytes data = Bytes.fromByteArray(b);

    byte[] actual = new byte[b.length];
    assertEquals(b.length, data.toInputStream().read(actual));
    assertEquals(
        com.google.common.primitives.Bytes.asList(b),
        com.google.common.primitives.Bytes.asList(actual));
  }

  public void testBytesReadVarint64() {
    byte[] b = new byte[] {(byte) 0b10101100, (byte) 0b00000010};
    Bytes data = Bytes.fromByteArray(b);

    Cursor cursor = data.cursor();
    assertEquals(300, data.readVarint64(cursor));
    assertEquals(2, cursor.position);
  }

  public void testBytesReadVarint64_malformed() {
    byte[] b = new byte[11];
    Arrays.fill(b, (byte) 0b11111111);
    Bytes data = Bytes.fromByteArray(b);

    Cursor cursor = data.cursor();
    // No support for assertThrows in GWT.
    try {
      long unused = data.readVarint64(cursor);
      throw new RuntimeException();
    } catch (IllegalArgumentException e) {
      assertEquals("Malformed varint.", e.getMessage());
    }
  }

  public void testBytesReadUintWithLength() {
    byte[] b = new byte[] {(byte) 0b11111111, (byte) 0b11111111};
    Bytes data = Bytes.fromByteArray(b);

    Cursor cursor = data.cursor();
    assertEquals(65535, data.readUintWithLength(cursor, b.length));
    assertEquals(b.length, cursor.position);
  }

  public void testBytesReadLittleEndianDouble() {
    byte[] b = new byte[] {0, 0, 0, 0, 0, -64, 94, 64};
    Bytes data = Bytes.fromByteArray(b);

    int position = 0;
    double tolerance = 0;
    assertEquals(123.0, data.readLittleEndianDouble(position), tolerance);
  }

  public void testLongsFromImmutableLongArray() {
    ImmutableLongArray expected = ImmutableLongArray.of(1, 2, 3);
    Longs actual = Longs.fromImmutableLongArray(expected);

    assertEquals(expected.length(), actual.length());
    for (int i = 0; i < actual.length(); i++) {
      assertEquals(expected.get(i), actual.get(i));
    }
  }

  public void testLongsToIntArray() {
    Longs actual = Longs.fromImmutableLongArray(ImmutableLongArray.of(1, 2, 3));
    int[] expected = actual.toIntArray();

    assertEquals(expected.length, actual.length());
    for (int i = 0; i < actual.length(); i++) {
      assertEquals(expected[i], actual.get(i));
    }

    // No support for assertThrows in GWT.
    try {
      int[] unused = Longs.fromImmutableLongArray(
          ImmutableLongArray.of(Long.MAX_VALUE)).toIntArray();
      throw new RuntimeException();
    } catch (IllegalArgumentException e) {
      // Do nothing.
    }
  }
}
