/*
 * Copyright 2023 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.parsePointsOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2BuilderUtil.IndexedLayer;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2PointVectorLayer}. */
@RunWith(JUnit4.class)
public final class S2PointVectorLayerTest extends GeometryTestCase {

  private void verifyS2PointVectorLayerResults(
      IntVector labelSetIds,
      IdSetLexicon labelSetLexicon,
      List<S2Point> output,
      String strExpectedPoints,
      List<IntVector> expectedLabels) {
    List<S2Point> expectedPoints = parsePointsOrDie(strExpectedPoints);

    assertEquals(expectedLabels.size(), labelSetIds.size());
    for (int i = 0; i < output.size(); ++i) {
      // Verify the expected point #i is equal to the output point #i.
      assertTrue(expectedPoints.get(i).equalsPoint(output.get(i)));
      // Verify that the expectedLabels for point #i has the same labels in the same order as
      // the output labels for point #i.
      int expectedNumLabels = expectedLabels.get(i).size();
      assertEquals(expectedNumLabels, labelSetLexicon.idSet(labelSetIds.get(i)).size());
      int[] idSet = labelSetLexicon.idSet(labelSetIds.get(i)).toArray();
      for (int k = 0; k < expectedNumLabels; k++) {
        int actualLabel = idSet[k];
        int expectedLabel = expectedLabels.get(i).get(k);
        assertEquals(expectedLabel, actualLabel);
      }
    }
  }

  private void addPoint(S2Point p, S2Builder builder) {
    builder.addEdge(p, p);
  }

  private IntVector v(int... labels) {
    return IntVector.copyOf(labels);
  }

  @Test
  public void testMergeDuplicates() {
    S2Builder builder = new S2Builder.Builder().build();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    IntVector labelSetIds = new IntVector();
    S2PointVectorLayer layer =
        new S2PointVectorLayer(
            labelSetIds, labelSetLexicon, new S2PointVectorLayer.Options(DuplicateEdges.MERGE));
    builder.startLayer(layer);

    builder.setLabel(1);
    addPoint(makePointOrDie("0:1"), builder);
    addPoint(makePointOrDie("0:2"), builder);
    builder.setLabel(2);
    addPoint(makePointOrDie("0:1"), builder);
    addPoint(makePointOrDie("0:4"), builder);
    addPoint(makePointOrDie("0:5"), builder);
    builder.clearLabels();
    addPoint(makePointOrDie("0:5"), builder);
    addPoint(makePointOrDie("0:6"), builder);
    S2Error error = new S2Error();
    assertTrue(builder.build(error));

    ImmutableList<IntVector> expectedLabels = ImmutableList.of(v(1, 2), v(1), v(2), v(2), v());
    String expectedPoints = "0:1, 0:2, 0:4, 0:5, 0:6";

    verifyS2PointVectorLayerResults(
        labelSetIds, labelSetLexicon, layer.getPointVector(), expectedPoints, expectedLabels);
  }

  @Test
  public void testKeepDuplicates() {
    S2Builder builder = new S2Builder.Builder().build();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    IntVector labelSetIds = new IntVector();
    S2PointVectorLayer layer =
        new S2PointVectorLayer(
            labelSetIds, labelSetLexicon, new S2PointVectorLayer.Options(DuplicateEdges.KEEP));
    builder.startLayer(layer);

    builder.setLabel(1);
    addPoint(makePointOrDie("0:1"), builder);
    addPoint(makePointOrDie("0:2"), builder);
    builder.setLabel(2);
    addPoint(makePointOrDie("0:1"), builder);
    addPoint(makePointOrDie("0:4"), builder);
    addPoint(makePointOrDie("0:5"), builder);
    builder.clearLabels();
    addPoint(makePointOrDie("0:5"), builder);
    addPoint(makePointOrDie("0:6"), builder);
    S2Error error = new S2Error();
    assertTrue(builder.build(error));

    ImmutableList<IntVector> expectedLabels =
        ImmutableList.of(v(1), v(2), v(1), v(2), v(2), v(), v());
    String expectedPoints = "0:1, 0:1, 0:2, 0:4, 0:5, 0:5, 0:6";

    verifyS2PointVectorLayerResults(
        labelSetIds, labelSetLexicon, layer.getPointVector(), expectedPoints, expectedLabels);
  }

  @Test
  public void testError() {
    S2Builder builder = new S2Builder.Builder().build();
    S2PointVectorLayer layer =
        new S2PointVectorLayer(new S2PointVectorLayer.Options(DuplicateEdges.KEEP));
    builder.startLayer(layer);

    addPoint(makePointOrDie("0:1"), builder);
    builder.addEdge(makePointOrDie("0:3"), makePointOrDie("0:4"));
    addPoint(makePointOrDie("0:5"), builder);
    S2Error error = new S2Error();
    assertFalse(builder.build(error));
    assertEquals(S2Error.Code.INVALID_ARGUMENT, error.code());
    assertEquals("Found non-degenerate edges", error.text());

    List<S2Point> output = layer.getPointVector();
    assertEquals(2, output.size());
    assertEquals(makePointOrDie("0:1"), output.get(0));
    assertEquals(makePointOrDie("0:5"), output.get(1));
  }

  @Test
  public void testIndexedS2PointVectorLayerAddsShapes() {
    S2Builder builder = new S2Builder.Builder().build();
    S2ShapeIndex index = new S2ShapeIndex();
    builder.startLayer(new IndexedLayer<S2PointVectorLayer>(index, new S2PointVectorLayer()));
    String point0Str = "0:0";
    String point1Str = "2:2";
    builder.addPoint(makePointOrDie(point0Str));
    builder.addPoint(makePointOrDie(point1Str));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    assertEquals(1, index.getShapes().size());
    S2Point.Shape shape = (S2Point.Shape) index.getShapes().get(0);
    assertEquals(2, shape.numEdges());
    assertEquals(point0Str, S2TextFormat.toString(shape.get(0)));
    assertEquals(point1Str, S2TextFormat.toString(shape.get(1)));
  }

  @Test
  public void testIndexedS2PointVectorLayerDoesNotAddEmptyShape() {
    S2Builder builder = new S2Builder.Builder().build();
    S2ShapeIndex index = new S2ShapeIndex();
    builder.startLayer(new IndexedLayer<S2PointVectorLayer>(index, new S2PointVectorLayer()));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    assertEquals(0, index.getShapes().size());
  }
}
