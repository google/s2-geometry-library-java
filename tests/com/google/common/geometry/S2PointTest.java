/*
 * Copyright 2013 Google Inc.
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

import com.google.common.annotations.GwtCompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.List;

/** Verifies S2Point. */
@GwtCompatible
public class S2PointTest extends GeometryTestCase {
  public void testVectorMethods() {
    S2Point a = new S2Point(1, 2, 3);
    S2Point b = new S2Point(0, 1, 0);
    assertEquals(new S2Point(1, 3, 3), a.add(b));
    assertEquals(new S2Point(1, 1, 3), a.sub(b));
    assertEquals(new S2Point(2, 4, 6), a.mul(2));
    assertEquals(new S2Point(0.5, 1, 1.5), a.div(2));
    assertEquals(2.0, a.dotProd(b), 0);
    assertEquals(new S2Point(-3, 0, 1), a.crossProd(b));
    assertEquals(b, b.normalize());
    assertEquals(a, a.neg().fabs());
    assertEquals(S2Point.Y_NEG, S2Point.X_POS.ortho());
  }

  public void testRotate() {
    // Check simple axial cases.
    S2Point p = makePoint("0:0");
    assertEquals(makePoint("0:90"), p.rotate(S2Point.Z_POS, Math.toRadians(90)), 1e-15);
    assertEquals(makePoint("-90:0"), p.rotate(S2Point.Y_POS, Math.toRadians(90)), 1e-15);
    assertEquals(makePoint("0:0"), p.rotate(S2Point.X_POS, Math.toRadians(90)), 1e-15);
    // Verify rotations in the plane containing the origin and two random points.
    for (int i = 0; i < 1000; i++) {
      S2Point a = randomPoint();
      S2Point b = randomPoint();
      S2Point axis = S2Point.crossProd(a, b);
      double angle = a.angle(b);
      // Rotate 'a' onto 'b'.
      assertEquals(b, a.rotate(axis, angle), 1e-15);
      // Rotate 'b' onto 'a'.
      assertEquals(a, b.rotate(axis, -angle), 1e-15);
      // Rotate 'a' to the midpoint of 'a' and 'b'.
      assertEquals(S2Point.normalize(S2Point.add(a, b)), a.rotate(axis, angle / 2), 1e-14);
      // Rotate 'b' to the antipodal point of 'a'.
      assertEquals(S2Point.neg(a), b.rotate(axis, Math.PI - angle), 1e-14);
    }
  }

  public void testS2Region() {
    S2Point point = new S2Point(1, 0, 0);

    S2Cap expectedCapBound = S2Cap.fromAxisHeight(point, 0);
    assertEquals(expectedCapBound, point.getCapBound());

    S2LatLng ll = new S2LatLng(point);
    S2LatLngRect expectedRect = new S2LatLngRect(ll, ll);
    assertEquals(expectedRect, point.getRectBound());

    // The leaf cell containing a point is still much larger than the point.
    S2Cell cell = new S2Cell(point);
    assertFalse(point.contains(cell));
    assertTrue(point.mayIntersect(cell));
  }

  public void testShape() {
    List<S2Point> points = ImmutableList.of(makePoint("0:0"), makePoint("1:1"), makePoint("2:2"));
    for (int j = 0; j < points.size(); j++) {
      List<S2Point> subset = points.subList(0, j);
      S2Shape shape = S2Point.Shape.fromList(subset);
      assertFalse(shape.containsOrigin());
      assertFalse(shape.hasInterior());
      assertEquals(subset.size(), shape.numEdges());
      MutableEdge edge = new MutableEdge();
      for (int i = 0; i < subset.size(); i++) {
        shape.getEdge(i, edge);
        assertEquals(edge.a, edge.b);
        assertEquals(subset.get(i), edge.a);
      }
    }
  }

  public void testSerialization() throws IOException {
    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    S2Point.X_NEG.encode(bos);
    assertEquals(S2Point.X_NEG, S2Point.decode(new ByteArrayInputStream(bos.toByteArray())));
  }

  public void testBuilder_origin() {
    assertEquals(S2Point.ORIGIN, new S2Point.Builder().build());
  }

  public void testBuilder_single() {
    S2Point base = new S2Point(0.1, 0.2, 0.3);
    S2Point built = new S2Point.Builder().add(base).build();
    assertEquals(base, built);
  }

  public void testBuilder_multiple() {
    assertEquals(
        new S2Point(0.51, -2.23, 2.73),
        new S2Point.Builder()
            .add(new S2Point(1, -2, 3))
            .add(new S2Point(-0.5, -0.25, -0.3))
            .add(new S2Point(0.01, 0.02, 0.03))
            .build());
  }

  public void testToBuilder() {
    S2Point base = new S2Point(0.25, 0.3, 0.5);
    S2Point.Builder builder = base.toBuilder();
    builder.add(new S2Point(0.1, 0.25, -0.5));
    assertEquals(new S2Point(0.35, 0.55, 0.0), builder.build());
    assertEquals(new S2Point(0.25, 0.3, 0.5), base); // Builder has a copy, not original
  }
}
