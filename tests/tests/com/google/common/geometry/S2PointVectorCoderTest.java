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
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.io.BaseEncoding;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for S2PointVectorCoder. */
@RunWith(JUnit4.class)
public class S2PointVectorCoderTest extends GeometryTestCase {

  /** Number of deltas per block in the implementation. */
  private static final int BLOCK_SIZE = 16;

  @Test
  public void testDecodeFromByteString() throws IOException {
    List<S2Point> expected = Lists.newArrayList(new S2Point(1, 2, 3), new S2Point(4, 5, 6));

    Bytes data =
        Bytes.fromByteArray(
            BaseEncoding.base16()
                .decode(
                    "00000010000000000000F03F00000000000000400000000000000840000000000000104000"
                        + "000000000014400000000000001840"));

    List<S2Point> actual = S2PointVectorCoder.FAST.decode(data, data.cursor(3));
    assertEquals(expected, actual);
  }

  @Test
  public void testEmpty() throws IOException {
    checkEncodeDecode(Lists.newArrayList(), S2PointVectorCoder.FAST, 1);

    // Test that an empty vector uses the FAST encoding.
    checkEncodeDecode(Lists.newArrayList(), S2PointVectorCoder.COMPACT, 1);
  }

  @Test
  public void testOnePoint() throws IOException {
    checkEncodeDecode(Lists.newArrayList(new S2Point(1, 0, 0)), S2PointVectorCoder.FAST, 25);

    // Encoding: header (2 bytes), block count (1 byte), block lengths (1 byte), block header (1
    // byte), delta (1 byte).
    checkEncodeDecode(Lists.newArrayList(new S2Point(1, 0, 0)), S2PointVectorCoder.COMPACT, 6);
  }

  @Test
  public void testOnePointWithExceptionsNoOverlap() throws IOException {
    // Test encoding a block with one point when other blocks have exceptions (which changes the
    // encoding for all blocks). The case below yields deltaBits == 8 and overlapBits == 0.
    //
    // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
    // Block 0: block header (1 byte), 16 deltas (16 bytes), exception (24 bytes)
    // Block 1: block header (1 byte), delta (1 byte)
    S2Point a = new S2Point(1, 0, 0);
    List<S2Point> points = new ArrayList<>();
    points.add(new S2Point(1, 2, 3).normalize());
    for (int i = 0; i < 15; ++i) {
      points.add(a);
    }
    points.add(a); // Second block
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 48);
  }

  @Test
  public void testOnePointWithExceptionsWithOverlap() throws IOException {
    // Test encoding a block with one point when other blocks have exceptions (which changes the
    // encoding for all blocks). The case below yields deltaBits == 8 and overlapBits == 4.
    //
    // Encoding: header (2 bytes), base (2 bytes), block count (1 byte), block offsets (2 bytes)
    // Block 0: header (1 byte), offset (2 bytes), 16 deltas (16 bytes), exception (24 bytes)
    // Block 1: header (1 byte), offset (2 bytes), delta (1 byte)
    S2Point a = new S2CellId(0x946df618d0000000L).toPoint();
    S2Point b = new S2CellId(0x947209e070000000L).toPoint();
    List<S2Point> points = new ArrayList<>();
    points.add(new S2Point(1, 2, 3).normalize());
    for (int i = 0; i < 15; ++i) {
      points.add(a);
    }
    points.add(b);  // Second block
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 54);
  }

  @Test
  public void testTenPoints() throws IOException {
    List<S2Point> points = new ArrayList<>();
    for (int i = 0; i < 10; ++i) {
      points.add(new S2Point(1, i, 0).normalize());
    }
    checkEncodeDecode(points, S2PointVectorCoder.FAST, 241);
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 231);
  }

  @Test
  public void testCellIdWithException() throws IOException {
    // Test one point encoded as an S2CellId with one point encoded as an exception.
    //
    // Encoding: header (2 bytes), block count (1 byte), block lengths (1 byte), block header (1
    // byte), two deltas (2 bytes), exception (24 bytes).
    checkEncodeDecode(
        Lists.newArrayList(
            S2CellId.fromDebugString("1/23").toPoint(), new S2Point(0.1, 0.2, 0.3).normalize()),
        S2PointVectorCoder.COMPACT,
        31);
  }

  @Test
  public void testPointsAtMultipleLevels() throws IOException {
    // Test that when points at multiple levels are present, the level with the most points is
    // chosen (preferring the smallest level in case of ties). (All other points are encoded as
    // exceptions.)

    // In this example, the two points at level 5 (on face 1) should be encoded. It is possible to
    // tell which points are encoded by the length of the encoding (since different numbers of
    // "base" bytes are encoded).
    //
    // Encoding: header (2 bytes), base (1 byte), block count (1 byte), block lengths (1 byte),
    // block header (1 byte), 5 deltas (5 bytes), S2Point exceptions (72 bytes).
    checkEncodeDecode(
        Lists.newArrayList(
            S2CellId.fromDebugString("2/11001310230102").toPoint(),
            S2CellId.fromDebugString("1/23322").toPoint(),
            S2CellId.fromDebugString("3/3").toPoint(),
            S2CellId.fromDebugString("1/23323").toPoint(),
            S2CellId.fromDebugString("2/12101023022012").toPoint()),
        S2PointVectorCoder.COMPACT,
        83);
  }

  @Test
  public void testPointsOnDifferentFaces() throws IOException {
    checkEncodeDecode(
        Lists.newArrayList(
            new S2CellId(0x87f627880299f7dbL).toPoint(),
            new S2CellId(0x52b332cd52cb8949L).toPoint()),
        S2PointVectorCoder.COMPACT,
        21);
  }

  @Test
  public void testNoOverlapOrExtraDeltaBitsNeeded() throws IOException {
    // This function tests the case in GetBlockCodes() where values can be encoded using the minimum
    // number of delta bits and no overlap. From the comments there:
    //
    //   Example 1: bMin = 0x72, bMax = 0x7e. The range is 0x0c. This can be encoded using
    //   deltaBits = 4 and overlapBits = 0, which allows us to represent an offset of 0x70 and a
    //   maximum delta of 0x0f, so that we can encode values up to 0x7f.
    //
    // To set up this test, we need at least two blocks: one to set the global minimum value, and
    // the other to encode a specific range of deltas. To make things easier, the first block has a
    // minimum value of zero.
    //
    // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
    // Block 0: header (1 byte), 8 deltas (8 bytes)
    // Block 1: header (1 byte), offset (1 byte), 4 deltas (2 bytes)
    int level = 3;
    List<S2Point> points =
        new ArrayList<>(Collections.nCopies(BLOCK_SIZE, encodedValueToPoint(0, level)));
    points.add(encodedValueToPoint(0x72, level));
    points.add(encodedValueToPoint(0x74, level));
    points.add(encodedValueToPoint(0x75, level));
    points.add(encodedValueToPoint(0x7e, level));
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 10 + BLOCK_SIZE / 2);
  }

  @Test
  public void testOverlapNeeded() throws IOException {
    // Like the above, but tests the following case:
    //
    //   Example 2: bMin = 0x78, bMax = 0x84. The range is 0x0c, but in this case it is not
    //   sufficient to use deltaBits = 4 and overlapBits = 0 because we can again only represent
    //   an offset of 0x70, so the maximum delta of 0x0f only lets us encode values up to 0x7f.
    //   However if we increase the overlap to 4 bits then we can represent an offset of 0x78, which
    //   lets us encode values up to 0x78 + 0x0f = 0x87.
    //
    // Encoding: header (2 bytes), block count (1 byte), block lengths (2 bytes)
    // Block 0: header (1 byte), 8 deltas (8 bytes)
    // Block 1: header (1 byte), offset (1 byte), 4 deltas (2 bytes)
    int level = 3;
    List<S2Point> points =
        new ArrayList<>(Collections.nCopies(BLOCK_SIZE, encodedValueToPoint(0, level)));
    points.add(encodedValueToPoint(0x78, level));
    points.add(encodedValueToPoint(0x7a, level));
    points.add(encodedValueToPoint(0x7c, level));
    points.add(encodedValueToPoint(0x84, level));
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 10 + BLOCK_SIZE / 2);
  }

  @Test
  public void testExtraDeltaBitsNeeded() throws IOException {
    // Like the above, but tests the following case:
    //
    //   Example 3: bMin = 0x08, bMax = 0x104. The range is 0xfc, so we should be able to use
    //   8-bit deltas. But even with a 4-bit overlap, we can still only encode offset = 0 and a
    //   maximum value of 0xff. (We don't allow bigger overlaps because statistically they are not
    //   worthwhile). Instead we increase the delta size to 12 bits, which handles this case easily.
    //
    // Encoding: header (2 bytes), block count (1 byte), block lengths (2 bytes)
    // Block 0: header (1 byte), 8 deltas (8 bytes)
    // Block 1: header (1 byte), 4 deltas (6 bytes)
    int level = 3;
    List<S2Point> points =
        new ArrayList<>(Collections.nCopies(BLOCK_SIZE, encodedValueToPoint(0, level)));
    points.add(encodedValueToPoint(0x08, level));
    points.add(encodedValueToPoint(0x4e, level));
    points.add(encodedValueToPoint(0x82, level));
    points.add(encodedValueToPoint(0x104, level));
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 13 + BLOCK_SIZE / 2);
  }

  @Test
  public void testExtraDeltaBitsAndOverlapNeeded() throws IOException {
    // Like the above, but tests the following case:
    //
    //   Example 4: bMin = 0xf08, bMax = 0x1004. The range is 0xfc, so we should be able to use
    //   8-bit deltas. With 8-bit deltas and no overlap, we have offset = 0xf00 and a maximum
    //   encodable value of 0xfff. With 8-bit deltas and a 4-bit overlap, we still have offset =
    //   0xf00 and a maximum encodable value of 0xfff. Even with 12-bit deltas, we have offset = 0
    //   and we can still only represent 0xfff. However with deltaBits = 12 and overlapBits = 4, we
    //   can represent offset = 0xf00 and a maximum encodable value of 0xf00 + 0xfff = 0x1eff.
    //
    // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
    // Block 0: header (1 byte), 8 deltas (8 bytes)
    // Block 1: header (1 byte), offset (1 byte), 4 deltas (6 bytes)
    int level = 5;
    List<S2Point> points =
        new ArrayList<>(Collections.nCopies(BLOCK_SIZE, encodedValueToPoint(0, level)));
    points.add(encodedValueToPoint(0xf08, level));
    points.add(encodedValueToPoint(0xf4e, level));
    points.add(encodedValueToPoint(0xf82, level));
    points.add(encodedValueToPoint(0x1004, level));
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 14 + BLOCK_SIZE / 2);
  }

  @Test
  public void testSixtyFourBitOffset() throws IOException {
    // Tests a case where a 64-bit block offset is needed.
    //
    // Encoding: header (2 bytes), block count (1 byte), block offsets (2 bytes)
    // Block 0: header (1 byte), 8 deltas (8 bytes)
    // Block 1: header (1 byte), offset (8 bytes), 2 deltas (1 byte)
    int level = S2CellId.MAX_LEVEL;
    List<S2Point> points =
        new ArrayList<>(Collections.nCopies(BLOCK_SIZE, S2CellId.begin(level).toPoint()));
    points.add(S2CellId.end(level).prev().toPoint());
    points.add(S2CellId.end(level).prev().prev().toPoint());
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 16 + BLOCK_SIZE / 2);
  }

  @Test
  public void testAllExceptionsBlock() throws IOException {
    // The encoding consists of two blocks; the first contains 16 encodable values, while the second
    // contains two exceptions.
    List<S2Point> points =
        new ArrayList<>(
            Collections.nCopies(BLOCK_SIZE, encodedValueToPoint(0, S2CellId.MAX_LEVEL)));
    points.add(new S2Point(0.1, 0.2, 0.3).normalize());
    points.add(new S2Point(0.3, 0.2, 0.1).normalize());
    // Encoding: header (2 bytes), block count (1 byte), block lengths (2 bytes).
    // 1st block header (1 byte), 16 deltas (16 bytes).
    // 2nd block header (1 byte), 2 deltas (1 byte), 2 exceptions (48 bytes).
    checkEncodeDecode(points, S2PointVectorCoder.COMPACT, 72);

    // Encoding: header (2 bytes), 18 S2Points (432 bytes).
    checkEncodeDecode(points, S2PointVectorCoder.FAST, 434);
  }

  @Test
  public void testFirstAtAllLevels() throws IOException {
    // Test encoding the first S2CellId at each level (which also happens to have the maximum face,
    // si, and ti values). All such S2CellIds can be encoded in 6 bytes because most of the bits
    // are zero.
    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      checkEncodeDecode(
          Lists.newArrayList(S2CellId.begin(level).toPoint()), S2PointVectorCoder.COMPACT, 6);
    }
  }

  @Test
  public void testLastAtAllLevels() throws IOException {
    // Test encoding the last S2CellId at each level. It turns out that such S2CellIds have the
    // largest possible face and ti values, and the minimum possible si value at that level. Such
    // S2CellIds can be encoded in 6 to 13 bytes depending on the level.
    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      // Note that 8 bit deltas are used to encode blocks of size 1, which reduces the size of
      // "base" from ((level + 2) / 4) to (level / 4) bytes.
      int expectedSize = 6 + level / 4;
      checkEncodeDecode(
          Lists.newArrayList(S2CellId.end(level).prev().toPoint()),
          S2PointVectorCoder.COMPACT,
          expectedSize);
    }
  }

  @Test
  public void testMaxFaceSiTiAtAllLevels() throws IOException {
    // Similar to the test above, but tests encoding the S2CellId at each level whose face, si, and
    // ti values are all maximal. This turns out to be the S2CellId whose human-readable form is
    // 5/222...22 (0xb555555555555555), however for clarity we construct it using
    // S2CellId.fromFaceIJ.
    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      S2CellId id =
          S2CellId.fromFaceIJ(5, S2CellId.MAX_SIZE - 1, S2CellId.MAX_SIZE - 1).parent(level);

      // This encoding is one byte bigger than the previous test at levels 7, 11, 15, 19, 23, and
      // 27. This is because in the previous test, the odd-numbered value bits are all zero (except
      // for the face number), which reduces the number of base bits needed by exactly 1. The
      // encoding size at level==3 is unaffected because for singleton blocks, the lowest 8 value
      // bits are encoded in the delta.
      int expectedSize = (level < 4) ? 6 : 6 + (level + 1) / 4;
      checkEncodeDecode(Lists.newArrayList(id.toPoint()), S2PointVectorCoder.COMPACT, expectedSize);
    }
  }

  @Test
  public void testLastTwoPointsAtAllLevels() throws IOException {
    // Test encoding the last two S2CellIds at each level.
    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      S2CellId id = S2CellId.end(level).prev();
      // Notice that this costs only 4 bits more than encoding the last S2CellId
      // by itself (see LastAtAllLevels). This is because encoding a block of
      // size 1 uses 8-bit deltas (which reduces the size of "base" by 4 bits),
      // while this test uses two 4-bit deltas.
      int expectedSize = 6 + (level + 2) / 4;
      checkEncodeDecode(
          Lists.newArrayList(id.toPoint(), id.prev().toPoint()),
          S2PointVectorCoder.COMPACT,
          expectedSize);
    }
  }

  @Test
  public void testManyDuplicatePointsAtAllLevels() throws IOException {
    // Test encoding 32 copies of the last S2CellId at each level. This uses between 27 and 38 bytes
    // depending on the level. (Note that the encoding can use less than 1 byte per point in this
    // situation.)
    for (int level = 0; level <= S2CellId.MAX_LEVEL; level++) {
      S2CellId id = S2CellId.end(level).prev();
      // Encoding: header (2 bytes), base ((level + 2) / 4 bytes), block count (1 byte), block
      // lengths (2 bytes), block headers (2 bytes), 32 deltas (16 bytes). At level 30 the encoding
      // size goes up by 1 byte because we can't encode an 8 byte "base" value, so instead this case
      // uses a base of 7 bytes plus a one-byte offset in each of the 2 blocks.
      int expectedSize = 23 + (level + 2) / 4;
      if (level == 30) {
        expectedSize += 1;
      }
      List<S2Point> points = new ArrayList<>(Collections.nCopies(32, id.toPoint()));
      checkEncodeDecode(points, S2PointVectorCoder.COMPACT, expectedSize);
    }
  }

  @Test
  public void testSnappedFractalLoops() throws IOException {
    final int maxPoints = 3 << 14;
    data.resetSeed();
    for (int numPoints = 3; numPoints <= maxPoints; numPoints *= 4) {
      int s2polygonSize = 0;
      int laxPolygonSize = 0;
      for (int i = 0; i < 10; ++i) {
        S2FractalBuilder builder = new S2FractalBuilder(data.rand);
        builder.setLevelForApproxMaxEdges(numPoints);
        Matrix frame = data.getRandomFrame();
        S2Loop loop = builder.makeLoop(frame, kmToAngle(10));
        List<S2Point> points = new ArrayList<>();
        for (int j = 0; j < loop.numVertices(); ++j) {
          points.add(S2CellId.fromPoint(loop.vertex(j)).toPoint());
        }
        S2Polygon s2polygon = new S2Polygon(new S2Loop(points));
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        s2polygon.encode(output);
        s2polygonSize += output.size();
        // S2LaxPolygonShape has 2 extra bytes of overhead to encode one loop.
        laxPolygonSize += checkEncodeDecode(points, S2PointVectorCoder.COMPACT, -1) + 2;
      }
      System.err.println(Platform.formatString(
          "numPoint=%5d  s2Polygon=%9d  laxPolygon=%9d",
          numPoints, s2polygonSize, laxPolygonSize));
    }
  }

  @Test
  public void testSnappedFractalPolylines() throws IOException {
    int maxPoints = 3 << 14;
    for (int numPoints = 3; numPoints <= maxPoints; numPoints *= 4) {
      data.setSeed(numPoints);

      S2FractalBuilder builder = new S2FractalBuilder(data.rand);
      builder.setLevelForApproxMaxEdges(numPoints);
      Matrix frame = data.getRandomFrame();
      S2Loop loop = builder.makeLoop(frame, kmToAngle(10));
      List<S2Point> points = new ArrayList<>();
      for (int j = 0; j < loop.numVertices(); j++) {
        points.add(S2CellId.fromPoint(loop.vertex(j)).toPoint());
      }
      checkEncodeDecode(points, S2PointVectorCoder.COMPACT, -1);
    }
  }

  /** Returns the number of bytes used to encode the input. */
  @CanIgnoreReturnValue
  private static int checkEncodeDecode(
      List<S2Point> input, S2Coder<List<S2Point>> coder, int expectedBytes) throws IOException {
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    coder.encode(input, output);

    PrimitiveArrays.Bytes data = PrimitiveArrays.Bytes.fromByteArray(output.toByteArray());
    Cursor cursor = data.cursor();
    List<S2Point> actual = coder.decode(data, cursor);

    if (expectedBytes >= 0) {
      assertEquals(expectedBytes, output.size());
      assertEquals(expectedBytes, cursor.position);
    }
    assertEquals(input, actual);
    return output.size();
  }

  // In order to make it easier to construct tests that encode particular values, this function
  // duplicates the part of EncodedS2PointVector that converts an encoded 64-bit value back to an
  // S2Point.
  private static S2Point encodedValueToPoint(long value, int level) {
    int sj = EncodedInts.deinterleaveBits1(value);
    int tj = EncodedInts.deinterleaveBits2(value);
    int shift = S2CellId.MAX_LEVEL - level;
    int si = (((sj << 1) | 1) << shift) & 0x7fffffff;
    int ti = (((tj << 1) | 1) << shift) & 0x7fffffff;
    int face = ((sj << shift) >>> 30) | (((tj << (shift + 1)) >>> 29) & 4);
    return S2Projections.faceUvToXyz(
            face,
            S2Projections.stToUV(S2Projections.siTiToSt(si)),
            S2Projections.stToUV(S2Projections.siTiToSt(ti)))
        .normalize();
  }
}
