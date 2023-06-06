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

import static com.google.common.geometry.TestDataGenerator.makePoint;

import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.io.BaseEncoding;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import junit.framework.TestCase;

/** Tests for S2CellIdVectorCoder. */
public class S2CellIdVectorCoderTest extends TestCase {

  private static final S2CellId NONE = S2CellId.none();
  private static final S2CellId SENTINEL = S2CellId.sentinel();

  public void testDecodeFromByteString() throws IOException {
    List<S2CellId> expected = Lists.newArrayList(SENTINEL, SENTINEL);

    Bytes bytes = Bytes.fromByteArray(BaseEncoding.base16().decode("00000007FFFFFFFFFFFFFF10FFFF"));
    List<S2CellId> actual = S2CellIdVectorCoder.INSTANCE.decode(bytes, bytes.cursor(3));

    assertEquals(expected, actual);
  }

  public void testEmpty() throws IOException {
    checkEncodeDecode(Lists.newArrayList(), 2);
  }

  public void testNone() throws IOException {
    checkEncodeDecode(Lists.newArrayList(NONE.id()), 3);
  }

  public void testNoneNone() throws IOException {
    checkEncodeDecode(Lists.newArrayList(NONE.id(), NONE.id()), 4);
  }

  public void testSentinel() throws IOException {
    checkEncodeDecode(Lists.newArrayList(SENTINEL.id()), 10);
  }

  public void testMaximumShiftCell() throws IOException {
    // Tests the encoding of a single cell at level 2, which corresponds to the maximum encodable
    // shift value (56).
    checkEncodeDecode(Lists.newArrayList(S2CellId.fromDebugString("0/00").id()), 3);
  }

  public void testSentinelSentinel() throws IOException {
    checkEncodeDecode(Lists.newArrayList(SENTINEL.id(), SENTINEL.id()), 11);
  }

  public void testNoneSentinelNone() throws IOException {
    checkEncodeDecode(Lists.newArrayList(NONE.id(), SENTINEL.id(), NONE.id()), 26);
  }

  public void testInvalidCells() throws IOException {
    // Tests that cells with an invalid LSB can be encoded.
    checkEncodeDecode(Lists.newArrayList(0x6L, 0xeL, 0x7eL), 5);
  }

  public void testOneByteLeafCells() throws IOException {
    // Tests that (1) if all cells are leaf cells, the low bit is not encoded, and (2) this can be
    // indicated using the standard 1-byte header.
    checkEncodeDecode(Lists.newArrayList(0x3L, 0x7L, 0x177L), 5);
  }

  public void testOneByteLevel29Cells() throws IOException {
    // Tests that (1) if all cells are at level 29, the low bit is not encoded, and (2) this can be
    // indicated using the standard 1-byte header.
    checkEncodeDecode(Lists.newArrayList(0xcL, 0x1cL, 0x47cL), 5);
  }

  public void testOneByteLevel28Cells() throws IOException {
    // Tests that (1) if all cells are at level 28, the low bit is not encoded, and (2) this can be
    // indicated using the extended 2-byte header.
    checkEncodeDecode(Lists.newArrayList(0x30L, 0x70L, 0x1770L), 6);
  }

  public void testOneByteMixedCellLevels() throws IOException {
    // Tests that cells at mixed levels can be encoded in one byte.
    checkEncodeDecode(Lists.newArrayList(0x300L, 0x1c00L, 0x7000L, 0xff00L), 6);
  }

  public void testOneByteMixedCellLevelsWithPrefix() throws IOException {
    // Tests that cells at mixed levels can be encoded in one byte even when they share a multi-byte
    // prefix.
    checkEncodeDecode(
        Lists.newArrayList(
            0x1234567800000300L, 0x1234567800001c00L,
            0x1234567800007000L, 0x123456780000ff00L),
        10);
  }

  public void testOneByteRangeWithBaseValue() throws IOException {
    // Tests that cells can be encoded in one byte by choosing a base value whose bit range overlaps
    // the delta values.
    // 1 byte header, 3 bytes base, 1 byte size, 4 bytes deltas.
    checkEncodeDecode(
        Lists.newArrayList(
            0x00ffff0000000000L, 0x0100fc0000000000L,
            0x0100500000000000L, 0x0100330000000000L),
        9);
  }

  public void testSixFaceCells() throws IOException {
    List<Long> expected = new ArrayList<>();
    for (int face = 0; face < 6; face++) {
      expected.add(S2CellId.fromFace(face).id());
    }
    checkEncodeDecode(expected, 8);
  }

  public void testFourLevel10Children() throws IOException {
    List<Long> expected = new ArrayList<>();
    S2CellId parent = S2CellId.fromDebugString("3/012301230");
    for (S2CellId id = parent.childBegin(); !id.equals(parent.childEnd()); id = id.next()) {
      expected.add(id.id());
    }
    checkEncodeDecode(expected, 8);
  }

  public void testFractalS2ShapeIndexCells() throws IOException {
    S2FractalBuilder fractalBuilder = new S2FractalBuilder(new Random(1));
    fractalBuilder.setLevelForApproxMaxEdges(3 * 1024);
    S2Point center = makePoint("47.677:-122.206");
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(fractalBuilder.makeLoop(S2.getFrame(center), S1Angle.degrees(1)));
    List<Long> expected = new ArrayList<>();
    S2Iterator<S2ShapeIndex.Cell> it = index.iterator();
    do {
      expected.add(it.id().id());
      it.next();
    } while (!it.done());
    assertEquals(969, expected.size());
    checkEncodeDecode(expected, 2911);
  }

  public void testCoveringCells() throws IOException {
    List<Long> expected =
        Lists.newArrayList(
            0x414a617f00000000L,
            0x414a61c000000000L,
            0x414a624000000000L,
            0x414a63c000000000L,
            0x414a647000000000L,
            0x414a64c000000000L,
            0x414a653000000000L,
            0x414a704000000000L,
            0x414a70c000000000L,
            0x414a714000000000L,
            0x414a71b000000000L,
            0x414a7a7c00000000L,
            0x414a7ac000000000L,
            0x414a8a4000000000L,
            0x414a8bc000000000L,
            0x414a8c4000000000L,
            0x414a8d7000000000L,
            0x414a8dc000000000L,
            0x414a914000000000L,
            0x414a91c000000000L,
            0x414a924000000000L,
            0x414a942c00000000L,
            0x414a95c000000000L,
            0x414a96c000000000L,
            0x414ab0c000000000L,
            0x414ab14000000000L,
            0x414ab34000000000L,
            0x414ab3c000000000L,
            0x414ab44000000000L,
            0x414ab4c000000000L,
            0x414ab6c000000000L,
            0x414ab74000000000L,
            0x414ab8c000000000L,
            0x414ab94000000000L,
            0x414aba1000000000L,
            0x414aba3000000000L,
            0x414abbc000000000L,
            0x414abe4000000000L,
            0x414abec000000000L,
            0x414abf4000000000L,
            0x46b5454000000000L,
            0x46b545c000000000L,
            0x46b5464000000000L,
            0x46b547c000000000L,
            0x46b5487000000000L,
            0x46b548c000000000L,
            0x46b5494000000000L,
            0x46b54a5400000000L,
            0x46b54ac000000000L,
            0x46b54b4000000000L,
            0x46b54bc000000000L,
            0x46b54c7000000000L,
            0x46b54c8004000000L,
            0x46b54ec000000000L,
            0x46b55ad400000000L,
            0x46b55b4000000000L,
            0x46b55bc000000000L,
            0x46b55c4000000000L,
            0x46b55c8100000000L,
            0x46b55dc000000000L,
            0x46b55e4000000000L,
            0x46b5604000000000L,
            0x46b560c000000000L,
            0x46b561c000000000L,
            0x46ca424000000000L,
            0x46ca42c000000000L,
            0x46ca43c000000000L,
            0x46ca444000000000L,
            0x46ca45c000000000L,
            0x46ca467000000000L,
            0x46ca469000000000L,
            0x46ca5fc000000000L,
            0x46ca604000000000L,
            0x46ca60c000000000L,
            0x46ca674000000000L,
            0x46ca679000000000L,
            0x46ca67f000000000L,
            0x46ca684000000000L,
            0x46ca855000000000L,
            0x46ca8c4000000000L,
            0x46ca8cc000000000L,
            0x46ca8e5400000000L,
            0x46ca8ec000000000L,
            0x46ca8f0100000000L,
            0x46ca8fc000000000L,
            0x46ca900400000000L,
            0x46ca98c000000000L,
            0x46ca994000000000L,
            0x46ca99c000000000L,
            0x46ca9a4000000000L,
            0x46ca9ac000000000L,
            0x46ca9bd500000000L,
            0x46ca9e4000000000L,
            0x46ca9ec000000000L,
            0x46caf34000000000L,
            0x46caf4c000000000L,
            0x46caf54000000000L);
    checkEncodeDecode(expected, 488);
  }

  public void testShift() throws IOException {
    // Test shift values which are larger than int32 (to test code paths like 1 << shift, where
    // shift >= 32).
    checkEncodeDecode(Lists.newArrayList(0x1689100000000000L), 5);
  }

  public void testLowerBoundLimits() throws IOException {
    // Test seeking before the beginning and past the end of the vector.
    S2CellId first = S2CellId.begin(S2CellId.MAX_LEVEL);
    S2CellId last = S2CellId.end(S2CellId.MAX_LEVEL).prev();

    ByteArrayOutputStream output = new ByteArrayOutputStream();
    S2CellIdVectorCoder.INSTANCE.encode(Lists.newArrayList(first, last), output);
    PrimitiveArrays.Bytes bytes = PrimitiveArrays.Bytes.fromByteArray(output.toByteArray());
    S2CellIdVector cellIds = S2CellIdVectorCoder.INSTANCE.decode(bytes, bytes.cursor());

    assertEquals(0, cellIds.lowerBound(NONE));
    assertEquals(0, cellIds.lowerBound(first));
    assertEquals(1, cellIds.lowerBound(first.next()));
    assertEquals(1, cellIds.lowerBound(last.prev()));
    assertEquals(1, cellIds.lowerBound(last));
    assertEquals(2, cellIds.lowerBound(last.next()));
    assertEquals(2, cellIds.lowerBound(SENTINEL));
  }

  // See EncodedS2CellIdVectorInitNeverCrashesRegression in encoded_s2cell_id_vector_test.cc.
  public void testEncodedS2CellIdVectorInitNeverCrashesRegression() throws IOException {
    // Previously, this would overflow the 32-bit multiplies of (size * bytesPerWord) and/or
    // (position * bytesPerWord) in UintVectorCoder.decode(), causing an
    // ArrayIndexOutOfBoundsException: "Index -2147483642 out of bounds for length 10".
    byte[] overflowBytes = {
      (byte) 32,
      (byte) 135,
      (byte) 128,
      (byte) 128,
      (byte) 128,
      (byte) 48,
      (byte) 39,
      (byte) 132,
      (byte) 143,
      (byte) 84
    };

    PrimitiveArrays.Bytes bytes = PrimitiveArrays.Bytes.fromByteArray(overflowBytes);
    try {
      S2CellIdVector unused = S2CellIdVectorCoder.INSTANCE.decode(bytes, bytes.cursor());
      fail("Expected IOException, didn't get it.");
    } catch (IOException expected) {
      // Expect "IOException: Decoding from 'data' with length 10 bytes, but 12884901894 bytes are
      // required."
    }
  }

  // Encodes the given vector and checks that it has the expected size and contents.
  private static void checkEncodeDecode(List<Long> expected, int expectedBytes) throws IOException {
    List<S2CellId> cellIds = new ArrayList<>();
    for (long id : expected) {
      cellIds.add(new S2CellId(id));
    }
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    S2CellIdVectorCoder.INSTANCE.encode(cellIds, output);
    Bytes data = Bytes.fromByteArray(output.toByteArray());
    Cursor cursor = data.cursor();
    List<S2CellId> actual = S2CellIdVectorCoder.INSTANCE.decode(data, cursor);

    assertEquals(expectedBytes, output.size());
    assertEquals(expectedBytes, cursor.position);
    assertEquals(cellIds, actual);
  }
}
