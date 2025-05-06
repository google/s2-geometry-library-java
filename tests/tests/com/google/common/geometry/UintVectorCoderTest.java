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
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.io.BaseEncoding;
import com.google.common.primitives.ImmutableLongArray;
import com.google.common.primitives.UnsignedInteger;
import com.google.common.primitives.UnsignedInts;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for UintVectorCoder. */
@RunWith(JUnit4.class)
public class UintVectorCoderTest {

  List<Long> decodeLongsFromBytes(Bytes data, int offset) throws IOException {
    Longs longs = UintVectorCoder.UINT64.decode(data, data.cursor(offset));
    List<Long> actual = new ArrayList<>();
    for (int i = 0; i < longs.length(); i++) {
      actual.add(longs.get(i));
    }
    return actual;
  }

  @Test
  public void testDecodeLongsFromByteString() throws IOException {
    List<Long> expected = Lists.newArrayList(0L, 0L, 0L);

    // The first 3 bytes of 'bytes' are padding placed to ensure we are handling offsets correctly.
    List<Long> actual = decodeLongsFromBytes(
        Bytes.fromByteArray(BaseEncoding.base16().decode("00000018000000")), 3);

    assertEquals(expected, actual);
  }

  @Test
  public void testEmpty_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList();

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(1, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testEmpty_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList();

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(1, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testZero_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(0L);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(2, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testZero_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(0);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(2, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testRepeatedZeros_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(0L, 0L, 0L);

    byte[] data = encodeUint64Vector(expected);
    assertEquals(4, data.length);

    List<Long> actual = decodeUint64Vector(data);
    assertEquals(expected, actual);
  }

  @Test
  public void testRepeatedZeros_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(0, 0, 0);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(4, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testMaxInt_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(~0L);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(9, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testMaxInt_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(~0);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(5, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testOneByte_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(0L, 255L, 1L, 254L);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(5, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testOneByte_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(0, 255, 1, 254);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(5, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testTwoBytes_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(0L, 255L, 256L, 254L);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(9, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testTwoBytes_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(0, 255, 256, 254);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(9, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testThreeBytes_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(0xffffffL, 0x0102L, 0L, 0x050403L);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(13, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testThreeBytes_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(0xffffff, 0x0102, 0, 0x050403);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(13, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testFourBytes_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(0xffffffffL, (long) Integer.MAX_VALUE);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(9, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testFourBytes_uint32() throws IOException {
    List<Integer> expected = Lists.newArrayList(0xffffffff, Integer.MAX_VALUE);

    byte[] data = encodeUint32Vector(expected);
    List<Integer> actual = decodeUint32Vector(data);

    assertEquals(9, data.length);
    assertEquals(expected, actual);
  }

  @Test
  public void testEightBytes_uint64() throws IOException {
    List<Long> expected = Lists.newArrayList(~0L, 0L, 0x0102030405060708L, Long.MAX_VALUE);

    byte[] data = encodeUint64Vector(expected);
    List<Long> actual = decodeUint64Vector(data);

    assertEquals(33, data.length);
    assertEquals(expected, actual);
  }

  private static byte[] encodeUint64Vector(List<Long> values) throws IOException {
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    UintVectorCoder.UINT64.encode(
        Longs.fromImmutableLongArray(ImmutableLongArray.copyOf(values)), output);
    return output.toByteArray();
  }

  private static List<Long> decodeUint64Vector(byte[] bytes) throws IOException {
    Bytes data = Bytes.fromByteArray(bytes);
    Longs longs = UintVectorCoder.UINT64.decode(data, data.cursor());
    List<Long> result = new ArrayList<>();
    for (int i = 0; i < longs.length(); i++) {
      result.add(longs.get(i));
    }
    return result;
  }

  private static byte[] encodeUint32Vector(List<Integer> values) throws IOException {
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    ImmutableLongArray.Builder builder = ImmutableLongArray.builder();
    for (int value : values) {
      builder.add(UnsignedInteger.fromIntBits(value).longValue());
    }
    UintVectorCoder.UINT32.encode(Longs.fromImmutableLongArray(builder.build()), output);
    return output.toByteArray();
  }

  private static List<Integer> decodeUint32Vector(byte[] bytes) throws IOException {
    Bytes data = Bytes.fromByteArray(bytes);
    Longs longs = UintVectorCoder.UINT32.decode(data, data.cursor());
    List<Integer> result = new ArrayList<>();
    for (int i = 0; i < longs.length(); i++) {
      result.add(UnsignedInts.checkedCast(longs.get(i)));
    }
    return result;
  }
}
