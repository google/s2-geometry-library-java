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

import static com.google.common.geometry.TestDataGenerator.kmToAngle;
import static com.google.common.geometry.TestDataGenerator.makePoint;
import static com.google.common.geometry.TestDataGenerator.makePolyline;

import com.google.common.collect.ImmutableMap;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.io.BaseEncoding;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

/** Tests for S2ShapeIndexCoder. */
public class S2ShapeIndexCoderTest extends GeometryTestCase {
  // TODO(user): Add benchmarks.

  private static final S2Shape oneEdge = makePolyline("1:1, 2:2");

  public void testEmpty() throws IOException {
    S2ShapeIndex expected = new S2ShapeIndex();
    checkEncodeDecode(expected, 4);
  }

  public void testNull() throws IOException {
    byte[] b = BaseEncoding.base16().decode("080028000000");
    S2ShapeIndex index = decode(b);
    assertNull(index.shapes.get(0));
  }

  public void testOneEdge() throws IOException {
    S2ShapeIndex expected = new S2ShapeIndex();
    expected.add(oneEdge);
    S2ShapeIndexTest.checkIteratorMethods(expected);
    checkEncodeDecode(expected, 8);
  }

  public void testOneEdgeNull() throws IOException {
    byte[] b =
        BaseEncoding.base16()
            .decode(
                "100036020102000000B4825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF913FFC7FB8B805F6EF3F"
                    + "28516A6D8FDBA13F27DCF7C958DEA13F28C809010408020010");
    S2ShapeIndex index = decode(b);
    assertNull(index.shapes.get(0));
    assertEquals(oneEdge, index.shapes.get(1));
  }

  public void testSameShapeTwice() throws IOException {
    S2ShapeIndex expected = new S2ShapeIndex();
    S2Polyline polyline = makePolyline("1:1, 2:2");
    expected.add(polyline);
    expected.add(polyline);
    checkEncodeDecode(expected, 12);
  }

  public void testRegularLoops() throws IOException {
    Map<Integer, Integer> testCases =
        ImmutableMap.<Integer, Integer>builder()
            .put(4, 8)
            .put(8, 8)
            .put(16, 16)
            .put(64, 77)
            .put(256, 327)
            .put(4096, 8813)
            .put(65536, 168291)
            .buildOrThrow();
    for (int numEdges : testCases.keySet()) {
      S2ShapeIndex expected = new S2ShapeIndex();
      S2Polygon poly =
          new S2Polygon(
              S2Loop.makeRegularLoop(
                  new S2Point(3, 2, 1).normalize(), S1Angle.degrees(0.1), numEdges));
      expected.add(poly.shape());
      checkEncodeDecode(expected, testCases.get(numEdges));
    }
  }

  public void testOverlappingPointClouds() throws IOException {
    int[][] testCases = new int[][] {{1, 50, 73}, {2, 100, 591}, {4, 100, 1461}};

    S2Cap cap = S2Cap.fromAxisAngle(new S2Point(0.1, -0.4, 0.3).normalize(), S1Angle.degrees(1));
    for (int[] testCase : testCases) {
      int numShapes = testCase[0];
      int numPointsPerShape = testCase[1];
      int expectedBytes = testCase[2];

      data.setSeed(numShapes);

      S2ShapeIndex expected = new S2ShapeIndex();
      for (int i = 0; i < numShapes; i++) {
        List<S2Point> points = new ArrayList<>();
        for (int j = 0; j < numPointsPerShape; j++) {
          points.add(data.samplePoint(cap));
        }
        expected.add(S2Point.Shape.fromList(points));
      }
      checkEncodeDecode(expected, expectedBytes);
    }
  }

  public void testOverlappingPolylines() throws IOException {
    int[][] testCases = new int[][] {{2, 50, 136}, {10, 50, 1015}, {20, 50, 2191}};

    S2Cap cap = S2Cap.fromAxisAngle(new S2Point(-0.2, -0.3, 0.4).normalize(), S1Angle.degrees(0.1));
    for (int[] testCase : testCases) {
      int numShapes = testCase[0];
      int numShapeEdges = testCase[1];
      int expectedBytes = testCase[2];

      data.setSeed(numShapes);

      S1Angle edgeLen = S1Angle.radians(2 * cap.angle().radians() / numShapeEdges);
      S2ShapeIndex expected = new S2ShapeIndex();
      for (int i = 0; i < numShapes; i++) {
        S2Point a = data.samplePoint(cap);
        S2Point b = data.getRandomPoint();
        List<S2Point> vertices = new ArrayList<>();
        for (int j = 0; j < numShapeEdges; j++) {
          vertices.add(
              S2EdgeUtil.interpolateAtDistance(S1Angle.radians(edgeLen.radians() * j), a, b));
        }
        expected.add(new S2Polyline(vertices));
      }
      checkEncodeDecode(expected, expectedBytes);
    }
  }

  public void testOverlappingLoops() throws IOException {
    int[][] testCases = new int[][] {{2, 250, 737}, {5, 250, 784}, {25, 50, 3848}};

    S2Cap cap = S2Cap.fromAxisAngle(new S2Point(-0.1, 0.25, 0.2).normalize(), S1Angle.degrees(3));
    for (int[] testCase : testCases) {
      int numShapes = testCase[0];
      int maxEdgesPerLoop = testCase[1];
      int expectedBytes = testCase[2];

      data.setSeed(numShapes);

      S2ShapeIndex expected = new S2ShapeIndex();
      for (int i = 0; i < numShapes; i++) {
        S2Point center = data.samplePoint(cap);
        double radiusFraction = data.nextDouble();
        // Scale the number of edges so that they are all about the same length (similar to modeling
        // all geometry at a similar resolution).
        int numEdges = Integer.max(3, (int) (maxEdgesPerLoop * radiusFraction));
        S2Polygon polygon =
            new S2Polygon(
                S2Loop.makeRegularLoop(
                    center, S1Angle.radians(cap.angle().radians() * radiusFraction), numEdges));
        expected.add(polygon.shape());
      }
      checkEncodeDecode(expected, expectedBytes);
    }
  }

  public void testSnappedFractalPolylines() throws IOException {
    S2ShapeIndex expected = new S2ShapeIndex();

    for (int i = 0; i < 5; i++) {
      data.setSeed(i);

      S2FractalBuilder builder = new S2FractalBuilder(data.rand);
      builder.setLevelForApproxMaxEdges(3 * 256);
      Matrix frame = S2.getFrame(makePoint("10:" + i));
      S2Loop loop = builder.makeLoop(frame, S1Angle.degrees(0.1));
      S2Polygon polygon = new S2Polygon();
      polygon.initToSnapped(new S2Polygon(loop), S2CellId.MAX_LEVEL);
      expected.add(polygon.shape());
    }
    checkEncodeDecode(expected, 10132);
  }

  public void testShapesThreadSafe() throws Exception {
    // Try reading each shape of an EncodedS2ShapeIndex in multiple threads at the same time. Given
    // that ByteBuffer is not thread-safe, if the decoding of each S2Shape is not properly
    // synchronized, this test should fail or timeout.
    int numIterations = 100;
    int numShapes = 20;

    for (int i = 0; i < numIterations; i++) {
      S2ShapeIndex index = new S2ShapeIndex();
      for (int j = 0; j < numShapes; j++) {
        int numVertices = 4 * data.skewed(10);
        index.add(
            new S2Polygon(S2Loop.makeRegularLoop(data.getRandomPoint(), kmToAngle(5), numVertices))
                .shape());
      }
      List<S2Shape> shapes = index.getShapes();
      ByteArrayOutputStream output = new ByteArrayOutputStream();
      VectorCoder.FAST_SHAPE.encode(shapes, output);
      S2ShapeIndexCoder.INSTANCE.encode(index, output);
      S2ShapeIndex encodedIndex = decode(output.toByteArray());

      CyclicBarrier start = new CyclicBarrier(numShapes + 1);
      CyclicBarrier end = new CyclicBarrier(numShapes + 1);

      final AtomicBoolean failed = new AtomicBoolean(false);
      IntStream.range(0, numShapes)
          .forEach(
              j ->
                  new Thread(
                          () -> {
                            try {
                              start.await();
                              if (!S2ShapeUtil.equals(shapes.get(j), encodedIndex.shapes.get(j))) {
                                failed.set(true);
                              }
                              end.await();
                            } catch (Exception e) {
                              throw new RuntimeException(e);
                            }
                          },
                          "reader-" + j)
                      .start());

      start.await();
      end.await();
      assertFalse(failed.get());
    }
  }

  public void testIteratorThreadSafe() throws Exception {
    // Try iterating over each cell of an EncodedS2ShapeIndex in multiple readers at the same time.
    int numIterations = 100;
    int numReaders = 8;

    for (int i = 0; i < numIterations; i++) {
      S2ShapeIndex index = new S2ShapeIndex();
      int numVertices = 4 * data.skewed(10);
      index.add(
          new S2Polygon(S2Loop.makeRegularLoop(data.getRandomPoint(), kmToAngle(5), numVertices))
              .shape());
      ByteArrayOutputStream output = new ByteArrayOutputStream();
      VectorCoder.FAST_SHAPE.encode(index.getShapes(), output);
      S2ShapeIndexCoder.INSTANCE.encode(index, output);

      S2ShapeIndex encodedIndex = decode(output.toByteArray());

      CyclicBarrier start = new CyclicBarrier(numReaders + 1);
      CyclicBarrier end = new CyclicBarrier(numReaders + 1);

      final AtomicBoolean failed = new AtomicBoolean(false);
      for (int j = 0; j < numReaders; j++) {
        new Thread(
                () -> {
                  try {
                    start.await();
                    S2Iterator<S2ShapeIndex.Cell> indexIter;
                    S2Iterator<S2ShapeIndex.Cell> encodedIndexIter;
                    for (indexIter = index.iterator(), encodedIndexIter = encodedIndex.iterator();
                        !indexIter.done() && !encodedIndexIter.done(); ) {
                      if (!S2ShapeUtil.equals(indexIter.entry(), encodedIndexIter.entry())) {
                        failed.set(true);
                        break;
                      }
                      indexIter.next();
                      encodedIndexIter.next();
                    }
                    if (!indexIter.done() || !encodedIndexIter.done()) {
                      failed.set(true);
                    }
                    end.await();
                  } catch (Exception e) {
                    throw new RuntimeException(e);
                  }
                },
                "reader-" + j)
            .start();
      }

      start.await();
      end.await();
      assertFalse(failed.get());
    }
  }

  private static S2ShapeIndex decode(byte[] bytes) throws IOException {
    Bytes data = Bytes.fromByteArray(bytes);
    Cursor cursor = data.cursor();
    List<S2Shape> shapes = VectorCoder.FAST_SHAPE.decode(data, cursor);
    return new S2ShapeIndexCoder(shapes).decode(data, cursor);
  }

  private static void checkEncodeDecode(S2ShapeIndex expected, int expectedBytes)
      throws IOException {
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    VectorCoder.FAST_SHAPE.encode(expected.getShapes(), output);
    S2ShapeIndexCoder.INSTANCE.encode(expected, output);

    Bytes data = Bytes.fromByteArray(output.toByteArray());
    Cursor cursor = data.cursor();
    List<S2Shape> shapes = VectorCoder.FAST_SHAPE.decode(data, cursor);
    long shapeIndexOffset = cursor.position;
    S2ShapeIndex actual = new S2ShapeIndexCoder(shapes).decode(data, cursor);

    assertEquals(expectedBytes, cursor.position - shapeIndexOffset);
    assertEquals(expected.options().getMaxEdgesPerCell(), actual.options().getMaxEdgesPerCell());
    assertEquals(expected, actual);
    S2ShapeIndexTest.checkIteratorMethods(actual);
  }

  private static void assertEquals(S2ShapeIndex expected, S2ShapeIndex actual) {
    assertTrue(S2ShapeUtil.equals(expected, actual));
  }

  private static void assertEquals(S2Shape expected, S2Shape actual) {
    assertTrue(S2ShapeUtil.equals(expected, actual));
  }
}
