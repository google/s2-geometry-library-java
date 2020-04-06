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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2Point.Shape;
import com.google.common.io.BaseEncoding;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.List;

@GwtIncompatible("S2LaxPolylineShape and S2LaxPolygonShape")
public class S2TaggedShapeCoderTest extends GeometryTestCase {

  public void testMixedShapes() throws IOException {
    S2ShapeIndex index = makeIndex("0:0 | 0:1 # 1:1, 1:2, 1:3 # 2:2; 2:3, 2:4, 3:3");
    index.add(S2LaxPolylineShape.create(makePolyline("1:1, 1:2, 1:3").vertices()));

    index.add(S2LaxPolygonShape.EMPTY);
    index.add(S2LaxPolygonShape.FULL);
    index.add(S2LaxPolygonShape.create(makePolygon("0:0, 0:4, 4:4, 4:0")));
    ByteArrayOutputStream output = new ByteArrayOutputStream();

    for (S2Shape shape : index.getShapes()) {
      S2TaggedShapeCoder.FAST.encode(shape, output);
      S2TaggedShapeCoder.COMPACT.encode(shape, output);
    }

    Bytes data = Bytes.fromByteArray(output.toByteArray());
    Cursor cursor = data.cursor();

    for (S2Shape expected : index.getShapes()) {
      S2Shape shape1 = S2TaggedShapeCoder.FAST.decode(data, cursor);
      S2Shape shape2 = S2TaggedShapeCoder.COMPACT.decode(data, cursor);

      assertTrue(S2ShapeUtil.equals(expected, shape1));
      assertTrue(S2ShapeUtil.equals(expected, shape2));
    }

    assertEquals(output.size(), cursor.position);
  }

  public void testDecodeFromByteString() {
    byte[] bytes =
        BaseEncoding.base16()
            .decode(
                "2932007C00E4002E0192010310000000000000F03F000000000000000000000000000000008AAFF597"
                    + "C0FEEF3F1EDD892B0BDF913F00000000000000000418B4825F3C81FDEF3F27DCF7C958DE913F"
                    + "1EDD892B0BDF913FD44A8442C3F9EF3FCE5B5A6FA6DDA13F1EDD892B0BDF913FAE0218F586F3"
                    + "EF3F3C3F66D2BBCAAA3F1EDD892B0BDF913F05010220FC7FB8B805F6EF3F28516A6D8FDBA13F"
                    + "27DCF7C958DEA13F96E20626CAEFEF3F4BF8A48399C7AA3F27DCF7C958DEA13F96B6DB0611E7"
                    + "EF3FC0221C80C6D8B13F27DCF7C958DEA13FE2337CCA8FE9EF3F6C573C9B60C2AA3F0EC9EF48"
                    + "C7CBAA3F0C0001040418B4825F3C81FDEF3F27DCF7C958DE913F1EDD892B0BDF913FD44A8442"
                    + "C3F9EF3FCE5B5A6FA6DDA13F1EDD892B0BDF913FAE0218F586F3EF3F3C3F66D2BBCAAA3F1EDD"
                    + "892B0BDF913F05010120000000000000F03F00000000000000000000000000000000F6FF7071"
                    + "0BECEF3F28516A6D8FDBB13F00000000000000003C4A985423D8EF3F199E8D966CD0B13F2851"
                    + "6A6D8FDBB13FF6FF70710BECEF3F000000000000000028516A6D8FDBB13F28C8390001000301"
                    + "04030405040704000738070E1B24292B3213000009030002130000110300092B000100010000"
                    + "00010D00000223041004020400020113082106110A4113000111030101");

    S2ShapeIndex expected = new S2ShapeIndex();
    expected.add(S2Point.Shape.fromList(Lists.newArrayList(makePoint("0:0"), makePoint("0:1"))));
    expected.add(S2LaxPolylineShape.create(makePolyline("1:1, 1:2, 1:3").vertices()));
    expected.add(
        S2LaxPolygonShape.create(
            Lists.newArrayList(
                makePolyline("2:2").vertices(), makePolyline("2:3, 2:4, 3:3").vertices())));
    expected.add(S2LaxPolylineShape.create(makePolyline("1:1, 1:2, 1:3").vertices()));
    expected.add(S2LaxPolygonShape.create(makePolygon("0:0, 0:4, 4:4, 4:0")));

    PrimitiveArrays.Bytes data = PrimitiveArrays.Bytes.fromByteArray(bytes);
    Cursor cursor = data.cursor();
    List<S2Shape> shapes = VectorCoder.FAST_SHAPE.decode(data, cursor);
    S2ShapeIndex actual = new S2ShapeIndexCoder(shapes).decode(data, cursor);
    assertTrue(S2ShapeUtil.equals(expected, actual));
  }

  public void testToBuilder() throws IOException {
    S2Point.Shape shape =
        new Shape() {
          @Override
          public S2Point get(int index) {
            return null;
          }

          @Override
          public int size() {
            return 0;
          }
        };

    ByteArrayOutputStream output = new ByteArrayOutputStream();
    try {
      S2TaggedShapeCoder.FAST.encode(shape, output);
      fail();
    } catch (IOException | IllegalArgumentException e) {
      assertTrue(e.getMessage().contains("No S2Coder matched S2Shape with type"));
    }

    output = new ByteArrayOutputStream();
    int typeTag = 8192;
    S2TaggedShapeCoder coder =
        S2TaggedShapeCoder.FAST.toBuilder()
            .add(shape.getClass(), S2Point.Shape.Coder.FAST, typeTag)
            .build();
    coder.encode(shape, output);
    assertTrue(output.size() > 0);
  }
}
