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

import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static com.google.common.geometry.TestDataGenerator.snapPointToLevel;
import static com.google.common.geometry.TestDataGenerator.snapPointsToLevel;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public strictfp class S2PointCompressionTest extends GeometryTestCase {

  // Four vertex loop near the corner of faces 0, 1, and 2.
  private List<S2Point> fourVertexLoop;

  // Four vertex loop near the corner of faces 0, 1, and 2; unsnapped.
  private List<S2Point> fourVertexUnsnappedLoop;

  // Four vertex loop near the corner of faces 0, 1, and 2; snapped to level 14.
  private List<S2Point> fourVertexLevel14Loop;

  // 100 vertex loop near the corner of faces 0, 1, and 2.
  private List<S2Point> oneHundredVertexLoop;

  // 100 vertex loop near the corner of faces 0, 1, and 2; unsnapped.
  private List<S2Point> oneHundredVertexUnsnappedLoop;

  // 100 vertex loop near the corner of faces 0, 1, and 2; 15 points snapped to S2CellId.MAX_LEVEL,
  // the others not snapped.
  private List<S2Point> oneHundredVertexMixed15Loop;

  // 100 vertex loop near the corner of faces 0, 1, and 2; 25 points snapped to S2CellId.MAX_LEVEL,
  // the others not snapped.
  private List<S2Point> oneHundredVertexMixed25Loop;

  // 100 vertex loop near the corner of faces 0, 1, and 2; snapped to level 22.
  private List<S2Point> oneHundredVertexLevel22Loop;

  // A loop with two vertices on each of three faces.
  private List<S2Point> multiFaceLoop;

  // A straight line of 100 vertices on face 0 that should compress well.
  private List<S2Point> straightLine;

  @Override
  public void setUp() {
    super.setUp();

    fourVertexLoop = makeRegularPoints(4, 0.1, S2CellId.MAX_LEVEL);

    fourVertexUnsnappedLoop = makeRegularPoints(4, 0.1, -1);

    // Radius is 100m, so points are about 141 meters apart. Snapping to level 14 will move them by
    // < 47m.
    fourVertexLevel14Loop = makeRegularPoints(4, 0.1, 14);

    oneHundredVertexLoop = makeRegularPoints(100, 0.1, S2CellId.MAX_LEVEL);

    oneHundredVertexUnsnappedLoop = makeRegularPoints(100, 0.1, -1);

    oneHundredVertexMixed15Loop = makeRegularPoints(100, 0.1, -1);
    for (int i = 0; i < 15; i++) {
      oneHundredVertexMixed15Loop.set(
          3 * i, snapPointToLevel(oneHundredVertexMixed15Loop.get(3 * i), S2CellId.MAX_LEVEL));
    }

    oneHundredVertexMixed25Loop = makeRegularPoints(100, 0.1, -1);
    for (int i = 0; i < 25; i++) {
      oneHundredVertexMixed25Loop.set(
          4 * i, snapPointToLevel(oneHundredVertexMixed15Loop.get(4 * i), S2CellId.MAX_LEVEL));
    }

    // Circumference is 628m, so points are about 6 meters apart. Snapping to level 22 will move
    // them by < 2m.
    oneHundredVertexLevel22Loop = makeRegularPoints(100, 0.1, 22);

    List<S2Point> multiFacePoints = new ArrayList<>(6);
    multiFacePoints.add(S2Projections.faceUvToXyz(0, -0.5, 0.5).normalize());
    multiFacePoints.add(S2Projections.faceUvToXyz(1, -0.5, 0.5).normalize());
    multiFacePoints.add(S2Projections.faceUvToXyz(0, 0.5, -0.5).normalize());
    multiFacePoints.add(S2Projections.faceUvToXyz(2, -0.5, 0.5).normalize());
    multiFacePoints.add(S2Projections.faceUvToXyz(2, 0.5, -0.5).normalize());
    multiFacePoints.add(S2Projections.faceUvToXyz(2, 0.5, 0.5).normalize());
    multiFaceLoop = snapPointsToLevel(multiFacePoints, S2CellId.MAX_LEVEL);

    List<S2Point> linePoints = new ArrayList<>(100);
    for (int i = 0; i < 100; i++) {
      double s = 0.01 + 0.005 * i;
      double t = 0.01 + 0.009 * i;
      double u = S2Projections.PROJ.stToUV(s);
      double v = S2Projections.PROJ.stToUV(t);
      linePoints.add(S2Projections.faceUvToXyz(0, u, v).normalize());
    }
    straightLine = snapPointsToLevel(linePoints, S2CellId.MAX_LEVEL);
  }

  // Make a regular loop around the corner of faces 0, 1, and 2 with the specified radius in meters
  // (on the earth) and number of vertices.
  private List<S2Point> makeRegularPoints(int numVertices, double radiusKm, int level) {
    S2Point center = new S2Point(1.0, 1.0, 1.0).normalize();
    S1Angle radiusAngle = kmToAngle(radiusKm);

    List<S2Point> points = S2Loop.makeRegularVertices(center, radiusAngle, numVertices);

    if (level < 0) {
      return points;
    } else {
      return snapPointsToLevel(points, level);
    }
  }

  private byte[] encode(List<S2Point> points, int level) throws IOException {
    ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
    S2PointCompression.encodePointsCompressed(points, level, outputStream);
    return outputStream.toByteArray();
  }

  private List<S2Point> decode(int numPoints, int level, byte[] encoded) throws IOException {
    List<S2Point> decodedPoints =
        S2PointCompression.decodePointsCompressed(
            numPoints, level, new ByteArrayInputStream(encoded));
    assertEquals(numPoints, decodedPoints.size());
    return decodedPoints;
  }

  private void checkRoundtrip(List<S2Point> points, int level) throws IOException {
    byte[] encoded = encode(points, level);
    List<S2Point> decodedPoints = decode(points.size(), level, encoded);
    assertEquals(points, decodedPoints);
  }

  public void testRoundtrips_empty() throws Exception {
    checkRoundtrip(new ArrayList<S2Point>(), S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_fourVertexLoop() throws Exception {
    checkRoundtrip(fourVertexLoop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_fourVertexUnsnappedLoop() throws Exception {
    checkRoundtrip(fourVertexUnsnappedLoop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_fourVertexLevel14Loop() throws Exception {
    checkRoundtrip(fourVertexLevel14Loop, 14);
  }

  public void testRoundtrip_oneHundredVertexLoop() throws Exception {
    checkRoundtrip(oneHundredVertexLoop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_oneHundredVertexUnsnappedLoop() throws Exception {
    checkRoundtrip(oneHundredVertexUnsnappedLoop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_oneHundredVertexMixed15Loop() throws Exception {
    checkRoundtrip(oneHundredVertexMixed15Loop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_oneHundredVertexMixed25Loop() throws Exception {
    checkRoundtrip(oneHundredVertexMixed25Loop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_oneHundredVertexLevel22Loop() throws Exception {
    checkRoundtrip(oneHundredVertexLevel22Loop, 22);
  }

  public void testRoundtrip_multiFaceLoop() throws Exception {
    checkRoundtrip(multiFaceLoop, S2CellId.MAX_LEVEL);
  }

  public void testRoundtrip_straightLine() throws Exception {
    checkRoundtrip(straightLine, S2CellId.MAX_LEVEL);
  }

  public void testSize_fourVertexLoop() throws Exception {
    byte[] encoded = encode(fourVertexLoop, S2CellId.MAX_LEVEL);
    // It would take 32 bytes uncompressed.
    assertEquals(39, encoded.length);
  }

  public void testSize_fourVertexLevel14Loop() throws Exception {
    byte[] encoded = encode(fourVertexLevel14Loop, 14);
    // It would take 4 bytes per vertex without compression.
    assertEquals(23, encoded.length);
  }

  public void testSize_oneHundredVertexLoop() throws Exception {
    byte[] encoded = encode(oneHundredVertexLoop, S2CellId.MAX_LEVEL);
    assertEquals(257, encoded.length);
  }

  public void testSize_oneHundredVertexUnsnappedLoop() throws Exception {
    byte[] encoded = encode(oneHundredVertexUnsnappedLoop, S2CellId.MAX_LEVEL);
    assertEquals(2756, encoded.length);
  }

  public void testSize_oneHundredVertexLevel22Loop() throws Exception {
    byte[] encoded = encode(oneHundredVertexLevel22Loop, 22);
    assertEquals(148, encoded.length);
  }

  public void testSize_straightLine() throws Exception {
    byte[] encoded = encode(straightLine, S2CellId.MAX_LEVEL);
    // About 1 byte / vertex.
    assertEquals(straightLine.size() + 17, encoded.length);
  }

  private static final int NUM_REGRESSION_CASES = 1000000;

  public void testNthDerivativeCoder_fixed() throws Exception {
    int[] input = {1, 5, 10, 15, 20, 23};
    int[] order0 = {1, 5, 10, 15, 20, 23};
    int[] order1 = {1, 4, 5, 5, 5, 3};
    int[] order2 = {1, 4, 1, 0, 0, -2};

    S2PointCompression.NthDerivativeCoder encoder0 = new S2PointCompression.NthDerivativeCoder(0);
    S2PointCompression.NthDerivativeCoder decoder0 = new S2PointCompression.NthDerivativeCoder(0);
    S2PointCompression.NthDerivativeCoder encoder1 = new S2PointCompression.NthDerivativeCoder(1);
    S2PointCompression.NthDerivativeCoder decoder1 = new S2PointCompression.NthDerivativeCoder(1);
    S2PointCompression.NthDerivativeCoder encoder2 = new S2PointCompression.NthDerivativeCoder(2);
    S2PointCompression.NthDerivativeCoder decoder2 = new S2PointCompression.NthDerivativeCoder(2);
    for (int i = 0; i < input.length; i++) {
      assertEquals(order0[i], encoder0.encode(input[i]));
      assertEquals(input[i], encoder0.encode(decoder0.decode(order0[i])));
      assertEquals(order1[i], encoder1.encode(input[i]));
      assertEquals(input[i], decoder1.decode(order1[i]));
      assertEquals(order2[i], encoder2.encode(input[i]));
      assertEquals(input[i], decoder2.decode(order2[i]));
    }
  }

  private void checkRegression(int order) {
    Random random = new Random();
    int[] raw = new int[NUM_REGRESSION_CASES];
    int[] encoded = new int[NUM_REGRESSION_CASES];
    int[] decoded = new int[NUM_REGRESSION_CASES];
    for (int i = 0; i < NUM_REGRESSION_CASES; i++) {
      raw[i] = random.nextInt();
    }

    S2PointCompression.NthDerivativeCoder encoder =
        new S2PointCompression.NthDerivativeCoder(order);
    for (int i = 0; i < NUM_REGRESSION_CASES; i++) {
      encoded[i] = encoder.encode(raw[i]);
    }

    S2PointCompression.NthDerivativeCoder decoder =
        new S2PointCompression.NthDerivativeCoder(order);
    for (int i = 0; i < NUM_REGRESSION_CASES; i++) {
      decoded[i] = decoder.decode(encoded[i]);
    }

    for (int i = 0; i < NUM_REGRESSION_CASES; i++) {
      assertEquals(raw[i], decoded[i]);
    }
  }

  public void testNthDerivativeCoder_regression() throws Exception {
    // Always test the lowest orders.
    for (int i = 0; i <= 2; i++) {
      checkRegression(i);
    }

    Random random = new Random();
    int skip = random.nextInt(1000);
    for (int i = 0; i < skip; i++) {
      random.nextInt();
    }

    for (int i = 0; i < 10; i++) {
      checkRegression(random.nextInt(S2PointCompression.NthDerivativeCoder.N_MAX + 1));
    }
  }
}
