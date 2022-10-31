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

import static com.google.common.geometry.EncodedInts.decodeZigZag32;
import static com.google.common.geometry.EncodedInts.decodeZigZag64;
import static com.google.common.geometry.EncodedInts.encodeZigZag32;
import static com.google.common.geometry.EncodedInts.encodeZigZag64;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import junit.framework.TestCase;

/** Unit tests for the {@link EncodedInts} class. */
public class EncodedIntsTest extends TestCase {
  public void testVarInts() throws IOException {
    ByteArrayOutputStream out = new ByteArrayOutputStream();
    long[] values = {Long.MIN_VALUE, -217, -1, 0, 1, 255, Long.MAX_VALUE};
    int expectedSize = 0;
    for (long v : values) {
      expectedSize += EncodedInts.varIntSize(v);
      EncodedInts.writeVarint64(out, v);
    }
    assertEquals("Expected different total size", expectedSize, out.size());
    InputStream in = new ByteArrayInputStream(out.toByteArray());
    for (long v : values) {
      assertEquals("Expected different value", v, EncodedInts.readVarint64(in));
    }
    assertTrue("Expected to reach EOF", in.read() < 0);
  }

  public void testEncodeZigZag32() {
    assertEquals(0, encodeZigZag32(0));
    assertEquals(1, encodeZigZag32(-1));
    assertEquals(2, encodeZigZag32(1));
    assertEquals(3, encodeZigZag32(-2));
    assertEquals(0x7FFFFFFE, encodeZigZag32(0x3FFFFFFF));
    assertEquals(0x7FFFFFFF, encodeZigZag32(0xC0000000));
    assertEquals(0xFFFFFFFE, encodeZigZag32(0x7FFFFFFF));
    assertEquals(0xFFFFFFFF, encodeZigZag32(0x80000000));

    // Some easier-to-verify round-trip tests. The inputs (other than 0, 1, -1) were chosen semi-
    // randomly via keyboard bashing.
    assertEquals(0, encodeZigZag32(decodeZigZag32(0)));
    assertEquals(1, encodeZigZag32(decodeZigZag32(1)));
    assertEquals(-1, encodeZigZag32(decodeZigZag32(-1)));
    assertEquals(14927, encodeZigZag32(decodeZigZag32(14927)));
    assertEquals(-3612, encodeZigZag32(decodeZigZag32(-3612)));
  }

  public void testEncodeZigZag64() {
    assertEquals(0, encodeZigZag64(0));
    assertEquals(1, encodeZigZag64(-1));
    assertEquals(2, encodeZigZag64(1));
    assertEquals(3, encodeZigZag64(-2));
    assertEquals(0x000000007FFFFFFEL, encodeZigZag64(0x000000003FFFFFFFL));
    assertEquals(0x000000007FFFFFFFL, encodeZigZag64(0xFFFFFFFFC0000000L));
    assertEquals(0x00000000FFFFFFFEL, encodeZigZag64(0x000000007FFFFFFFL));
    assertEquals(0x00000000FFFFFFFFL, encodeZigZag64(0xFFFFFFFF80000000L));
    assertEquals(0xFFFFFFFFFFFFFFFEL, encodeZigZag64(0x7FFFFFFFFFFFFFFFL));
    assertEquals(0xFFFFFFFFFFFFFFFFL, encodeZigZag64(0x8000000000000000L));

    // Some easier-to-verify round-trip tests. The inputs (other than 0, 1, -1) were chosen semi-
    // randomly via keyboard bashing.
    assertEquals(0, encodeZigZag64(decodeZigZag64(0)));
    assertEquals(1, encodeZigZag64(decodeZigZag64(1)));
    assertEquals(-1, encodeZigZag64(decodeZigZag64(-1)));
    assertEquals(14927, encodeZigZag64(decodeZigZag64(14927)));
    assertEquals(-3612, encodeZigZag64(decodeZigZag64(-3612)));

    assertEquals(856912304801416L, encodeZigZag64(decodeZigZag64(856912304801416L)));
    assertEquals(-75123905439571256L, encodeZigZag64(decodeZigZag64(-75123905439571256L)));
  }

  public void testDecodeZigZag32() {
    assertEquals(0, decodeZigZag32(0));
    assertEquals(-1, decodeZigZag32(1));
    assertEquals(1, decodeZigZag32(2));
    assertEquals(-2, decodeZigZag32(3));
    assertEquals(0x3FFFFFFF, decodeZigZag32(0x7FFFFFFE));
    assertEquals(0xC0000000, decodeZigZag32(0x7FFFFFFF));
    assertEquals(0x7FFFFFFF, decodeZigZag32(0xFFFFFFFE));
    assertEquals(0x80000000, decodeZigZag32(0xFFFFFFFF));
  }

  public void testDecodeZigZag64() {
    assertEquals(0, decodeZigZag64(0));
    assertEquals(-1, decodeZigZag64(1));
    assertEquals(1, decodeZigZag64(2));
    assertEquals(-2, decodeZigZag64(3));
    assertEquals(0x000000003FFFFFFFL, decodeZigZag64(0x000000007FFFFFFEL));
    assertEquals(0xFFFFFFFFC0000000L, decodeZigZag64(0x000000007FFFFFFFL));
    assertEquals(0x000000007FFFFFFFL, decodeZigZag64(0x00000000FFFFFFFEL));
    assertEquals(0xFFFFFFFF80000000L, decodeZigZag64(0x00000000FFFFFFFFL));
    assertEquals(0x7FFFFFFFFFFFFFFFL, decodeZigZag64(0xFFFFFFFFFFFFFFFEL));
    assertEquals(0x8000000000000000L, decodeZigZag64(0xFFFFFFFFFFFFFFFFL));
  }

  public void testInterleaveBits() {
    checkBits(0xC000000000000000L, 0x80000000, 0x80000000);
    checkBits(0xEAAAAAAA2AAAAA88L, 0x80000000, 0xFFFF7FFA);
    checkBits(0xEAAAAAAAAAAAA8AAL, 0x80000000, 0xFFFFFFEF);
    checkBits(0xEAAAAAAAAAAAAAAAL, 0x80000000, 0xFFFFFFFF);
    checkBits(0x4000000000000000L, 0x80000000, 0x00000000);
    checkBits(0x40000000000000A8L, 0x80000000, 0x0000000E);
    checkBits(0x6AAAAAAAAAAAAA22L, 0x80000000, 0x7FFFFFF5);
    checkBits(0xD555555515555544L, 0xFFFF7FFA, 0x80000000);
    checkBits(0xFFFFFFFF3FFFFFCCL, 0xFFFF7FFA, 0xFFFF7FFA);
    checkBits(0xFFFFFFFFBFFFFDEEL, 0xFFFF7FFA, 0xFFFFFFEF);
    checkBits(0xFFFFFFFFBFFFFFEEL, 0xFFFF7FFA, 0xFFFFFFFF);
    checkBits(0x5555555515555544L, 0xFFFF7FFA, 0x00000000);
    checkBits(0x55555555155555ECL, 0xFFFF7FFA, 0x0000000E);
    checkBits(0x7FFFFFFFBFFFFF66L, 0xFFFF7FFA, 0x7FFFFFF5);
    checkBits(0xD555555555555455L, 0xFFFFFFEF, 0x80000000);
    checkBits(0xFFFFFFFF7FFFFEDDL, 0xFFFFFFEF, 0xFFFF7FFA);
    checkBits(0xFFFFFFFFFFFFFCFFL, 0xFFFFFFEF, 0xFFFFFFEF);
    checkBits(0xFFFFFFFFFFFFFEFFL, 0xFFFFFFEF, 0xFFFFFFFF);
    checkBits(0x5555555555555455L, 0xFFFFFFEF, 0x00000000);
    checkBits(0x55555555555554FDL, 0xFFFFFFEF, 0x0000000E);
    checkBits(0x7FFFFFFFFFFFFE77L, 0xFFFFFFEF, 0x7FFFFFF5);
    checkBits(0xD555555555555555L, 0xFFFFFFFF, 0x80000000);
    checkBits(0xFFFFFFFF7FFFFFDDL, 0xFFFFFFFF, 0xFFFF7FFA);
    checkBits(0xFFFFFFFFFFFFFDFFL, 0xFFFFFFFF, 0xFFFFFFEF);
    checkBits(0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFF, 0xFFFFFFFF);
    checkBits(0x5555555555555555L, 0xFFFFFFFF, 0x00000000);
    checkBits(0x55555555555555FDL, 0xFFFFFFFF, 0x0000000E);
    checkBits(0x7FFFFFFFFFFFFF77L, 0xFFFFFFFF, 0x7FFFFFF5);
    checkBits(0x8000000000000000L, 0x00000000, 0x80000000);
    checkBits(0xAAAAAAAA2AAAAA88L, 0x00000000, 0xFFFF7FFA);
    checkBits(0xAAAAAAAAAAAAA8AAL, 0x00000000, 0xFFFFFFEF);
    checkBits(0xAAAAAAAAAAAAAAAAL, 0x00000000, 0xFFFFFFFF);
    checkBits(0x0000000000000000L, 0x00000000, 0x00000000);
    checkBits(0x00000000000000A8L, 0x00000000, 0x0000000E);
    checkBits(0x2AAAAAAAAAAAAA22L, 0x00000000, 0x7FFFFFF5);
    checkBits(0x8000000000000054L, 0x0000000E, 0x80000000);
    checkBits(0xAAAAAAAA2AAAAADCL, 0x0000000E, 0xFFFF7FFA);
    checkBits(0xAAAAAAAAAAAAA8FEL, 0x0000000E, 0xFFFFFFEF);
    checkBits(0xAAAAAAAAAAAAAAFEL, 0x0000000E, 0xFFFFFFFF);
    checkBits(0x0000000000000054L, 0x0000000E, 0x00000000);
    checkBits(0x00000000000000FCL, 0x0000000E, 0x0000000E);
    checkBits(0x2AAAAAAAAAAAAA76L, 0x0000000E, 0x7FFFFFF5);
    checkBits(0x9555555555555511L, 0x7FFFFFF5, 0x80000000);
    checkBits(0xBFFFFFFF7FFFFF99L, 0x7FFFFFF5, 0xFFFF7FFA);
    checkBits(0xBFFFFFFFFFFFFDBBL, 0x7FFFFFF5, 0xFFFFFFEF);
    checkBits(0xBFFFFFFFFFFFFFBBL, 0x7FFFFFF5, 0xFFFFFFFF);
    checkBits(0x1555555555555511L, 0x7FFFFFF5, 0x00000000);
    checkBits(0x15555555555555B9L, 0x7FFFFFF5, 0x0000000E);
    checkBits(0x3FFFFFFFFFFFFF33L, 0x7FFFFFF5, 0x7FFFFFF5);
  }

  private static void checkBits(long expected, int v1, int v2) {
    long bits = EncodedInts.interleaveBits(v1, v2);
    assertEquals(expected, bits);
    assertEquals(v1, EncodedInts.deinterleaveBits1(bits));
    assertEquals(v2, EncodedInts.deinterleaveBits2(bits));
  }

  public void testInterleaveBitPairs() {
    checkBitPairs(0xA000000000000000L, 0x80000000, 0x80000000);
    checkBitPairs(0xECCCCCCC4CCCCC88L, 0x80000000, 0xFFFF7FFA);
    checkBitPairs(0xECCCCCCCCCCCC8CCL, 0x80000000, 0xFFFFFFEF);
    checkBitPairs(0xECCCCCCCCCCCCCCCL, 0x80000000, 0xFFFFFFFF);
    checkBitPairs(0x2000000000000000L, 0x80000000, 0x00000000);
    checkBitPairs(0x20000000000000C8L, 0x80000000, 0x0000000E);
    checkBitPairs(0x6CCCCCCCCCCCCC44L, 0x80000000, 0x7FFFFFF5);
    checkBitPairs(0xB333333313333322L, 0xFFFF7FFA, 0x80000000);
    checkBitPairs(0xFFFFFFFF5FFFFFAAL, 0xFFFF7FFA, 0xFFFF7FFA);
    checkBitPairs(0xFFFFFFFFDFFFFBEEL, 0xFFFF7FFA, 0xFFFFFFEF);
    checkBitPairs(0xFFFFFFFFDFFFFFEEL, 0xFFFF7FFA, 0xFFFFFFFF);
    checkBitPairs(0x3333333313333322L, 0xFFFF7FFA, 0x00000000);
    checkBitPairs(0x33333333133333EAL, 0xFFFF7FFA, 0x0000000E);
    checkBitPairs(0x7FFFFFFFDFFFFF66L, 0xFFFF7FFA, 0x7FFFFFF5);
    checkBitPairs(0xB333333333333233L, 0xFFFFFFEF, 0x80000000);
    checkBitPairs(0xFFFFFFFF7FFFFEBBL, 0xFFFFFFEF, 0xFFFF7FFA);
    checkBitPairs(0xFFFFFFFFFFFFFAFFL, 0xFFFFFFEF, 0xFFFFFFEF);
    checkBitPairs(0xFFFFFFFFFFFFFEFFL, 0xFFFFFFEF, 0xFFFFFFFF);
    checkBitPairs(0x3333333333333233L, 0xFFFFFFEF, 0x00000000);
    checkBitPairs(0x33333333333332FBL, 0xFFFFFFEF, 0x0000000E);
    checkBitPairs(0x7FFFFFFFFFFFFE77L, 0xFFFFFFEF, 0x7FFFFFF5);
    checkBitPairs(0xB333333333333333L, 0xFFFFFFFF, 0x80000000);
    checkBitPairs(0xFFFFFFFF7FFFFFBBL, 0xFFFFFFFF, 0xFFFF7FFA);
    checkBitPairs(0xFFFFFFFFFFFFFBFFL, 0xFFFFFFFF, 0xFFFFFFEF);
    checkBitPairs(0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFF, 0xFFFFFFFF);
    checkBitPairs(0x3333333333333333L, 0xFFFFFFFF, 0x00000000);
    checkBitPairs(0x33333333333333FBL, 0xFFFFFFFF, 0x0000000E);
    checkBitPairs(0x7FFFFFFFFFFFFF77L, 0xFFFFFFFF, 0x7FFFFFF5);
    checkBitPairs(0x8000000000000000L, 0x00000000, 0x80000000);
    checkBitPairs(0xCCCCCCCC4CCCCC88L, 0x00000000, 0xFFFF7FFA);
    checkBitPairs(0xCCCCCCCCCCCCC8CCL, 0x00000000, 0xFFFFFFEF);
    checkBitPairs(0xCCCCCCCCCCCCCCCCL, 0x00000000, 0xFFFFFFFF);
    checkBitPairs(0x0000000000000000L, 0x00000000, 0x00000000);
    checkBitPairs(0x00000000000000C8L, 0x00000000, 0x0000000E);
    checkBitPairs(0x4CCCCCCCCCCCCC44L, 0x00000000, 0x7FFFFFF5);
    checkBitPairs(0x8000000000000032L, 0x0000000E, 0x80000000);
    checkBitPairs(0xCCCCCCCC4CCCCCBAL, 0x0000000E, 0xFFFF7FFA);
    checkBitPairs(0xCCCCCCCCCCCCC8FEL, 0x0000000E, 0xFFFFFFEF);
    checkBitPairs(0xCCCCCCCCCCCCCCFEL, 0x0000000E, 0xFFFFFFFF);
    checkBitPairs(0x0000000000000032L, 0x0000000E, 0x00000000);
    checkBitPairs(0x00000000000000FAL, 0x0000000E, 0x0000000E);
    checkBitPairs(0x4CCCCCCCCCCCCC76L, 0x0000000E, 0x7FFFFFF5);
    checkBitPairs(0x9333333333333311L, 0x7FFFFFF5, 0x80000000);
    checkBitPairs(0xDFFFFFFF7FFFFF99L, 0x7FFFFFF5, 0xFFFF7FFA);
    checkBitPairs(0xDFFFFFFFFFFFFBDDL, 0x7FFFFFF5, 0xFFFFFFEF);
    checkBitPairs(0xDFFFFFFFFFFFFFDDL, 0x7FFFFFF5, 0xFFFFFFFF);
    checkBitPairs(0x1333333333333311L, 0x7FFFFFF5, 0x00000000);
    checkBitPairs(0x13333333333333D9L, 0x7FFFFFF5, 0x0000000E);
    checkBitPairs(0x5FFFFFFFFFFFFF55L, 0x7FFFFFF5, 0x7FFFFFF5);
  }

  private static void checkBitPairs(long expected, int v1, int v2) {
    long bits = EncodedInts.interleaveBitPairs(v1, v2);
    assertEquals(expected, bits);
    assertEquals(v1, EncodedInts.deinterleaveBitPairs1(bits));
    assertEquals(v2, EncodedInts.deinterleaveBitPairs2(bits));
  }
}
