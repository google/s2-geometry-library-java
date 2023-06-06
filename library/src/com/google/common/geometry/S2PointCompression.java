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

import static java.lang.Math.min;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.geometry.S2Projections.FaceSiTi;
import com.google.common.primitives.Longs;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Given a sequence of S2Points assumed to be the center of level-k cells, compresses it into a
 * stream using the following method:
 *
 * <ul>
 *   <li>decompose the points into (face, si, ti) tuples (see {@link S2Projections.FaceSiTi}).
 *   <li>run-length encode the faces, combining face number and count into a varint32. See the
 *       {@link FaceRunCoder} class.
 *   <li>right shift the (si, ti) to remove the part that's constant for all cells of level-k. The
 *       result is called the (pi, qi) space.
 *   <li>2nd derivative encode the pi and qi sequences (linear prediction).
 *   <li>zig-zag encode all derivative values but the first, which cannot be negative.
 *   <li>interleave the zig-zag encoded values.
 *   <li>encode the first interleaved value in a fixed length encoding (varint would make this value
 *       larger).
 *   <li>encode the remaining interleaved values as varint64s, as the derivative encoding should
 *       make the values small.
 * </ul>
 *
 * <p>In addition, provides a lossless method to compress a sequence of points even if some points
 * are not the center of level-k cells. These points are stored exactly, using 3 double precision
 * values, after the above encoded string, together with their index in the sequence (this leads to
 * some redundancy - it is expected that only a small fraction of the points are not cell centers).
 *
 * <p>To encode leaf cells, this requires 8 bytes for the first vertex plus an average of 3.8 bytes
 * for each additional vertex.
 */
public final strictfp class S2PointCompression {

  private S2PointCompression() {}

  private static final int DERIVATIVE_ENCODING_ORDER = 2;

  /**
   * Encode a list of points into an efficient, lossless binary representation, which can be decoded
   * by calling {@link S2PointCompression#decodePointsCompressed}. The encoding is byte-compatible
   * with the C++ version of the S2 library.
   *
   * <p>Points that are snapped to the specified level will require approximately 4 bytes per point,
   * while other points will require 24 bytes per point.
   *
   * @param points The list of points to encode.
   * @param level The {@link S2Cell} level at which points should be encoded.
   * @param output The output stream into which the encoding should be written.
   * @throws IOException if there was a problem writing into the output stream.
   */
  public static void encodePointsCompressed(List<S2Point> points, int level, OutputStream output)
      throws IOException {
    encodePointsCompressed(points, level, new LittleEndianOutput(output));
  }

  /**
   * If the list of points is empty, no faces will be encoded, and no (pi, qi) coordinates, but
   * a single Varint32 of value zero will be written for the number of off-center points.
   */
  static void encodePointsCompressed(List<S2Point> points, int level, LittleEndianOutput encoder)
      throws IOException {
    // Convert the points to (face, pi, qi) coordinates.
    FaceRunCoder faces = new FaceRunCoder();
    int[] verticesPi = new int[points.size()];
    int[] verticesQi = new int[points.size()];
    List<Integer> offCenter = new ArrayList<>();
    for (int i = 0; i < points.size(); i++) {
      FaceSiTi faceSiTi = S2Projections.PROJ.xyzToFaceSiTi(points.get(i));
      faces.addFace(faceSiTi.face);
      verticesPi[i] = siTiToPiQi(faceSiTi.si, level);
      verticesQi[i] = siTiToPiQi(faceSiTi.ti, level);
      if (S2Projections.PROJ.levelIfCenter(faceSiTi, points.get(i)) != level) {
        offCenter.add(i);
      }
    }

    // Encode the runs of the faces.
    faces.encode(encoder);

    // Encode the (pi, qi) coordinates of all the points, in order.
    NthDerivativeCoder piCoder = new NthDerivativeCoder(DERIVATIVE_ENCODING_ORDER);
    NthDerivativeCoder qiCoder = new NthDerivativeCoder(DERIVATIVE_ENCODING_ORDER);
    for (int i = 0; i < verticesPi.length; i++) {
      int pi = piCoder.encode(verticesPi[i]);
      int qi = qiCoder.encode(verticesQi[i]);
      if (i == 0) {
        // The first point will be just the (pi, qi) coordinates of the S2Point. NthDerivativeCoder
        // will not save anything in that case, so we encode in fixed format rather than varint to
        // avoid the varint overhead.

        // Interleave to reduce overhead from two partial bytes to one.
        long interleavedPiQi = EncodedInts.interleaveBits(pi, qi);

        // Java uses big-endian representation, but the wire format requires little-endian.
        // Simultaneously, we only write as many bytes as are actually required for the given level,
        // i.e. we wish to truncate the byte representation of the long. We do this by reversing the
        // bytes of the value, then truncating the byte array representing the result to the
        // required length. Note that if level is zero, bytesRequired is zero, and nothing is
        // written.
        int bytesRequired = (level + 7) / 8 * 2;
        byte[] littleEndianInterleavedPiQiBytes =
            Arrays.copyOf(Longs.toByteArray(Long.reverseBytes(interleavedPiQi)), bytesRequired);
        encoder.writeBytes(littleEndianInterleavedPiQiBytes);
      } else {
        // ZigZagEncode, as varint requires the maximum number of bytes for negative numbers.
        int zigZagEncodedPi = EncodedInts.encodeZigZag32(pi);
        int zigZagEncodedQi = EncodedInts.encodeZigZag32(qi);

        // Interleave to reduce overhead from two partial bytes to one.
        long interleavedPiQi = EncodedInts.interleaveBits(zigZagEncodedPi, zigZagEncodedQi);
        encoder.writeVarint64(interleavedPiQi);
      }
    }

    // Encode the number of off-center points.
    encoder.writeVarint32(offCenter.size());

    // Encode the actual off-center points.
    for (int index : offCenter) {
      encoder.writeVarint32(index);
      points.get(index).encode(encoder);
    }
  }

  /**
   * Decode a list of points that were encoded using {@link
   * S2PointCompression#encodePointsCompressed}.
   *
   * <p>Points that are snapped to the specified level will require approximately 4 bytes per point,
   * while other points will require 24 bytes per point.
   *
   * @param numVertices The number of points to decode.
   * @param level The {@link S2Cell} level at which points are encoded.
   * @param input The input stream containing the encoded point data.
   * @return the list of decoded points.
   * @throws IOException if there was a problem reading from the input stream.
   */
  public static List<S2Point> decodePointsCompressed(int numVertices, int level, InputStream input)
      throws IOException {
    return decodePointsCompressed(numVertices, level, new LittleEndianInput(input));
  }

  static List<S2Point> decodePointsCompressed(int numVertices, int level, LittleEndianInput decoder)
      throws IOException {
    List<S2Point> vertices = new ArrayList<>(numVertices);
    if (level > S2CellId.MAX_LEVEL || level < 0) {
      throw new IllegalArgumentException("Invalid S2Cell level provided: " + level);
    }

    FaceRunCoder faces = new FaceRunCoder();
    faces.decode(numVertices, decoder);
    Iterator<Integer> faceIterator = faces.getFaceIterator();

    NthDerivativeCoder piCoder = new NthDerivativeCoder(DERIVATIVE_ENCODING_ORDER);
    NthDerivativeCoder qiCoder = new NthDerivativeCoder(DERIVATIVE_ENCODING_ORDER);
    for (int i = 0; i < numVertices; i++) {
      int pi;
      int qi;
      if (i == 0) {
        // The interleaved coordinates are stored in a truncated (depending on level) little-endian
        // representation, but we need big-endian for Java, so reconstruct the necessary bytes here.
        // Do this by reading in the bytes as-is, padding the end with zeros to get the full 8-byte
        // array, converting to a long (still little-endian), and finally reversing the byte
        // representation to get big-endian. Note that if level is zero, bytesRequired is zero, and
        // nothing is read.
        int bytesRequired = (level + 7) / 8 * 2;
        byte[] littleEndianBytes = decoder.readBytes(bytesRequired);
        long interleavedPiQi =
            Long.reverseBytes(Longs.fromByteArray(Arrays.copyOf(littleEndianBytes, Longs.BYTES)));
        pi = piCoder.decode(EncodedInts.deinterleaveBits1(interleavedPiQi));
        qi = qiCoder.decode(EncodedInts.deinterleaveBits2(interleavedPiQi));
      } else {
        long piqi = decoder.readVarint64();
        pi = piCoder.decode(EncodedInts.decodeZigZag32(EncodedInts.deinterleaveBits1(piqi)));
        qi = qiCoder.decode(EncodedInts.decodeZigZag32(EncodedInts.deinterleaveBits2(piqi)));
      }
      int face = faceIterator.next();
      vertices.add(facePiQiToXyz(face, pi, qi, level));
    }

    // Now decode the off-center points.
    int numOffCenter = decoder.readVarint32();
    if (numOffCenter > numVertices) {
      throw new IOException("Number of off-center points is greater than total amount of points.");
    }
    for (int i = 0; i < numOffCenter; i++) {
      int index = decoder.readVarint32();
      double x = decoder.readDouble();
      double y = decoder.readDouble();
      double z = decoder.readDouble();
      try {
        vertices.set(index, new S2Point(x, y, z));
      } catch (IndexOutOfBoundsException e) {
        throw new IOException("Insufficient or invalid data: ", e);
      }
    }

    return vertices;
  }

  private static int siTiToPiQi(long si, int level) {
    si = min(si, S2Projections.MAX_SITI - 1);
    return (int) (si >>> (S2CellId.MAX_LEVEL + 1 - level));
  }

  private static double piQiToST(int pi, int level) {
    // We want to recover the position at the center of the cell. If the point was snapped to the
    // center of the cell, then modf(s * 2^level) == 0.5. Inverting STtoPiQi gives:
    // s = (pi + 0.5) / 2^level.
    return (pi + 0.5) / (1 << level);
  }

  private static S2Point facePiQiToXyz(int face, int pi, int qi, int level) {
    return S2Projections.faceUvToXyz(
            face,
            S2Projections.PROJ.stToUV(piQiToST(pi, level)),
            S2Projections.PROJ.stToUV(piQiToST(qi, level))).normalize();
  }

  private static class FaceRunCoder {
    private static class FaceRun {
      public FaceRun(int face, int count) {
        this.face = face;
        this.count = count;
      }

      public int face;
      public int count;
    }

    private final List<FaceRun> faces = new ArrayList<>();

    public void addFace(int face) {
      FaceRun lastRun = !faces.isEmpty() ? Iterables.getLast(faces) : null;
      if (lastRun != null && lastRun.face == face) {
        lastRun.count += 1;
      } else {
        faces.add(new FaceRun(face, 1));
      }
    }

    /** Writes the list of FaceRuns to the encoder. If 'faces' is an empty list, writes nothing. */
    public void encode(LittleEndianOutput encoder) throws IOException {
      for (FaceRun run : faces) {
        // It isn't necessary to encode the number of faces left for the last run, but since this
        // would only help if there were more than 21 faces, it will be a small overall savings,
        // much smaller than the bound encoding.
        encoder.writeVarint64(S2CellId.NUM_FACES * (long) run.count + run.face);
      }
    }

    /** Reads 'vertices' FaceRuns from the decoder. If 'vertices' is zero, reads nothing. */
    public void decode(int vertices, LittleEndianInput decoder) throws IOException {
      int facesParsed = 0;
      while (facesParsed < vertices) {
        long faceAndCount = decoder.readVarint64();
        int face = (int) (faceAndCount % S2CellId.NUM_FACES);
        int count = (int) (faceAndCount / S2CellId.NUM_FACES);
        if (faceAndCount < 0) {
          throw new IOException("Invalid face: " + face + ", from faceAndCount: " + faceAndCount);
        }
        if (count < 0) {
          throw new IOException("Invalid count: " + count + ", from faceAndCount: " + faceAndCount);
        }
        FaceRun run = new FaceRun(face, count);
        faces.add(run);
        facesParsed += run.count;
      }
    }

    public Iterator<Integer> getFaceIterator() {
      final Iterator<FaceRun> faceRunIterator = faces.iterator();
      // Special case if there are not faces at all.
      if (!faceRunIterator.hasNext()) {
        return Collections.emptyIterator();
      }
      return new Iterator<Integer>() {
        private FaceRun currentFaceRun = faceRunIterator.next();
        private int usedCountForCurrentFaceRun = 0;

        @Override
        public boolean hasNext() {
          return usedCountForCurrentFaceRun < currentFaceRun.count || faceRunIterator.hasNext();
        }

        @Override
        public void remove() {
          throw new UnsupportedOperationException();
        }

        @Override
        public Integer next() {
          if (usedCountForCurrentFaceRun < currentFaceRun.count) {
            usedCountForCurrentFaceRun++;
          } else {
            usedCountForCurrentFaceRun = 1;
            currentFaceRun = faceRunIterator.next();
          }
          return currentFaceRun.face;
        }
      };
    }
  }

  @VisibleForTesting
  static final class NthDerivativeCoder {

    // The range of supported Ns is [N_MIN, N_MAX].
    public static final int N_MIN = 0;
    public static final int N_MAX = 10;

    // The derivative order of the coder (the N in NthDerivative).
    private final int n;

    // The derivative order in which to code the next value (ramps up to n).
    private int m;

    // Value memory, from oldest to newest.
    private final int[] memory;

    public NthDerivativeCoder(int n) {
      Preconditions.checkArgument(N_MIN <= n && n <= N_MAX, "Unsupported N: %s", n);
      this.n = n;
      memory = new int[N_MAX];
      reset();
    }

    public int getN() {
      return n;
    }

    public int encode(int k) {
      for (int i = 0; i < m; i++) {
        int delta = k - memory[i];
        memory[i] = k;
        k = delta;
      }
      if (m < n) {
        memory[m] = k;
        m++;
      }
      return k;
    }

    public int decode(int k) {
      if (m < n) {
        m++;
      }
      for (int i = m - 1; i >= 0; i--) {
        memory[i] += k;
        k = memory[i];
      }
      return k;
    }

    public void reset() {
      Arrays.fill(memory, 0, n, 0);
      m = 0;
    }
  }
}
