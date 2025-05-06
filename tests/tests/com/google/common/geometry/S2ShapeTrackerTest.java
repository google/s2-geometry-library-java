/*
 * Copyright 2024 Google Inc.
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

import static com.google.common.geometry.S2Cell.Boundary.BOTTOM_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.LEFT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.RIGHT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.TOP_EDGE;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Collections2;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Cell.Boundary;
import com.google.common.geometry.S2RobustCellClipper.Crossing;
import com.google.common.geometry.S2RobustCellClipper.CrossingType;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class S2ShapeTrackerTest extends GeometryTestCase {

  @Test
  public void testPointWorks() {
    S2ShapeTracker tracker = new S2ShapeTracker(0, 3);
    assertFalse(tracker.finished());

    // Mark chains as tracked.
    tracker.markChain(0);
    tracker.markChain(1);
    tracker.markChain(2);
    assertTrue(tracker.finished());

    // Chain tracking is idempotent, it shouldn't matter if we mark them again.
    tracker.markChain(0);
    tracker.markChain(1);
    tracker.markChain(2);
    assertTrue(tracker.finished());
  }

  @Test
  public void testPolylineWorks() {
    S2ShapeTracker tracker = new S2ShapeTracker(1, 3);
    assertFalse(tracker.finished());

    for (int i = 0; i < 3; ++i) {
      tracker.markChain(i);
    }

    // Turn on a couple crossings
    tracker.addPoint(0, 1, 1, 11);
    tracker.delPoint(3, 0, 1337, 13);
    tracker.addPoint(4, 1, 3141, 17);
    assertFalse(tracker.finished());

    // Turn them off one by one.
    tracker.delPoint(0, 1, 1, 11);
    assertFalse(tracker.finished());

    tracker.addPoint(3, 0, 1337, 13);
    assertFalse(tracker.finished());

    tracker.delPoint(4, 1, 3141, 17);
    assertTrue(tracker.finished());
  }

  @Test
  public void testPolygonWorks() {
    S2ShapeTracker tracker = new S2ShapeTracker(2, 5);
    assertFalse(tracker.finished());

    for (int i = 0; i < 5; ++i) {
      tracker.markChain(i);
    }

    // Add a couple of intervals.
    tracker.addInterval(0, 1, 1000, 11, 22);
    tracker.addInterval(3, 0, 1337, 23, 36);
    tracker.addInterval(4, 1, 3141, 55, 72);
    assertFalse(tracker.finished());

    // Subtract them off one by one.
    tracker.addInterval(0, 1, 1000, 22, 18);
    tracker.addInterval(0, 1, 1000, 18, 11);
    assertFalse(tracker.finished());

    tracker.addInterval(3, 0, 1337, 25, 24);
    tracker.addInterval(3, 0, 1337, 24, 23);
    tracker.addInterval(3, 0, 1337, 30, 25);
    tracker.addInterval(3, 0, 1337, 36, 30);
    assertFalse(tracker.finished());

    tracker.addInterval(4, 1, 3141, 72, 55);
    assertTrue(tracker.finished());
  }

  @Test
  public void testFaceAxesCancel() {
    S2ShapeTracker tracker = new S2ShapeTracker(2, 1);
    tracker.markChain(0);
    assertTrue(tracker.finished());

    // Test that adjacent face axes cancel properly.
    S2Cell cell0 = S2Cell.fromFace(0);
    S2Cell cell1 = S2Cell.fromFace(1);

    tracker.addInterval(
        0, 0, S2Projections.LIMIT_IJ, //
        cell0.getIJCoordOfBoundary(BOTTOM_EDGE),
        cell0.getIJCoordOfBoundary(TOP_EDGE));
    assertFalse(tracker.finished());

    tracker.addInterval(
        0, 1, S2Projections.LIMIT_IJ, //
        cell0.getIJCoordOfBoundary(RIGHT_EDGE),
        cell0.getIJCoordOfBoundary(BOTTOM_EDGE));
    assertFalse(tracker.finished());

    tracker.addInterval(
        1, 1, 0, //
        cell1.getIJCoordOfBoundary(TOP_EDGE),
        cell1.getIJCoordOfBoundary(BOTTOM_EDGE));
    assertFalse(tracker.finished());

    tracker.addInterval(
        2, 0, 0, //
        cell1.getIJCoordOfBoundary(TOP_EDGE),
        cell1.getIJCoordOfBoundary(BOTTOM_EDGE));
    assertTrue(tracker.finished());
  }

  @Test
  public void testFaceCellsClose() {
    ImmutableList<Integer> faces = ImmutableList.of(0, 1, 2, 3, 4, 5);

    // Adding all 6 face cells should naturally sum to zero, and the order shouldn't matter.
    for (List<Integer> facePermutation : Collections2.permutations(faces)) {
      S2ShapeTracker tracker = new S2ShapeTracker(2, 1);
      assertFalse(tracker.finished());
      tracker.markChain(0);
      assertTrue(tracker.finished());

      int cnt = 0;
      for (int face : facePermutation) {
        tracker.addCellBoundary(S2Cell.fromFace(face));
        assertEquals(++cnt == 6, tracker.finished());
      }
    }
  }

  @Test
  public void testFaceChildrenClose() {
    S2ShapeTracker tracker = new S2ShapeTracker(2, 1);
    assertFalse(tracker.finished());
    tracker.markChain(0);
    assertTrue(tracker.finished());

    // Adding all the children of all the faces at a given level should close.
    for (int face = 0; face < 6; ++face) {
      S2CellId faceId = S2CellId.fromFace(face);
      for (int level = 1; level <= 2; ++level) {
        S2CellId beg = faceId.childBegin(level);
        S2CellId end = faceId.childEnd(level);

        for (S2CellId id = beg; !id.equals(end); id = id.next()) {
          tracker.addCellBoundary(new S2Cell(id));

          boolean done = face == 5 && level == 2 && id.next().equals(end);
          assertEquals(done, tracker.finished());
        }
      }
    }
  }

  @Test
  public void testFaceCornersClose() {
    // We'll add intervals that contain the corners of each face cell, they should sum to zero.
    S2ShapeTracker tracker = new S2ShapeTracker(2, 1);
    assertFalse(tracker.finished());
    tracker.markChain(0);
    assertTrue(tracker.finished());

    // We have to be careful to order the crossings correctly on each boundary.
    ImmutableList<Crossing> crossings = ImmutableList.of(
        crossing(BOTTOM_EDGE, CrossingType.INCOMING, -1, -0.75, 0),
        crossing(BOTTOM_EDGE, CrossingType.OUTGOING, -1, +0.75, 0),
        crossing(RIGHT_EDGE, CrossingType.INCOMING, +1, -0.75, 0),
        crossing(RIGHT_EDGE, CrossingType.OUTGOING, +1, +0.75, 0),
        crossing(TOP_EDGE, CrossingType.INCOMING, +1, +0.75, 0),
        crossing(TOP_EDGE, CrossingType.OUTGOING, +1, -0.75, 0),
        crossing(LEFT_EDGE, CrossingType.INCOMING, -1, +0.75, 0),
        crossing(LEFT_EDGE, CrossingType.OUTGOING, -1, -0.75, 0));

    for (int face = 0; face < 6; ++face) {
      tracker.processCrossings(S2Cell.fromFace(face), crossings);
      assertEquals(face == 5, tracker.finished());
    }
  }

  @Test
  public void testSmallIntervalsWork() {
    // We'll add very small intervals near the face corners.  The rounding should expand them to
    // unit length and they should cancel properly.
    //
    // This should work for both polylines and polygons.
    for (int dim = 1; dim <= 2; ++dim) {
      S2ShapeTracker tracker = new S2ShapeTracker(dim, 1);
      assertFalse(tracker.finished());
      tracker.markChain(0);
      assertTrue(tracker.finished());

      // We have to be careful to order the crossings correctly on each boundary.
      ImmutableList<Crossing> crossings =
          ImmutableList.of(
              crossing(BOTTOM_EDGE, CrossingType.INCOMING, -1, -(1 - 1e-15), 0),
              crossing(BOTTOM_EDGE, CrossingType.OUTGOING, -1, +(1 - 1e-15), 0),
              crossing(RIGHT_EDGE, CrossingType.INCOMING, +1, -(1 - 1e-15), 0),
              crossing(RIGHT_EDGE, CrossingType.OUTGOING, +1, +(1 - 1e-15), 0),
              crossing(TOP_EDGE, CrossingType.INCOMING, +1, +(1 - 1e-15), 0),
              crossing(TOP_EDGE, CrossingType.OUTGOING, +1, -(1 - 1e-15), 0),
              crossing(LEFT_EDGE, CrossingType.INCOMING, -1, +(1 - 1e-15), 0),
              crossing(LEFT_EDGE, CrossingType.OUTGOING, -1, -(1 - 1e-15), 0));

      for (int face = 0; face < 6; ++face) {
        tracker.processCrossings(S2Cell.fromFace(face), crossings);
        assertEquals(face == 5, tracker.finished());
      }
    }
  }

  @Test
  public void testSmallTeeIntervalsWork() {
    // Add a very small interval near the cell center. We'll add it in children 0 and 3 of face 1,
    // and to face 0 directly, this forms a T intersection, but the tiny interval should be rounded
    // outward and still cancel at that point.
    S2ShapeTracker tracker = new S2ShapeTracker(2, 1);
    assertFalse(tracker.finished());
    tracker.markChain(0);
    assertTrue(tracker.finished());

    S2Cell cell;
    cell = new S2Cell(S2CellId.fromFace(1).child(3));
    tracker.processCrossings(
        cell,
        ImmutableList.of(
            crossing(BOTTOM_EDGE, CrossingType.INCOMING, -1, -0.75, 0),
            crossing(LEFT_EDGE, CrossingType.OUTGOING, -1, +1e-15, 0)));
    assertFalse(tracker.finished());

    cell = new S2Cell(S2CellId.fromFace(1).child(0));
    tracker.processCrossings(
        cell,
        ImmutableList.of(
            crossing(TOP_EDGE, CrossingType.OUTGOING, +1, -0.75, 0),
            crossing(LEFT_EDGE, CrossingType.INCOMING, -1, -1e-15, 0)));
    assertFalse(tracker.finished());

    cell = S2Cell.fromFace(0);
    tracker.processCrossings(
        cell,
        ImmutableList.of(
            crossing(RIGHT_EDGE, CrossingType.OUTGOING, +1, -1e-15, 0),
            crossing(RIGHT_EDGE, CrossingType.INCOMING, +1, +1e-15, 0)));
    assertTrue(tracker.finished());
  }

  @Test
  public void testNoTrivialCollisions() {
    // Try to find some trivial collisions by just adding two intervals with different endpoints but
    // the same length. S2ShapeTracker shouldn't be fooled.

    for (int i = 0; i < 50000; ++i) {
      // Generate random face, coord and axis.
      int face = data.uniformInt(0, 6);
      int coord = data.uniformInt(0, S2Projections.LIMIT_IJ);
      int axis = data.uniformInt(0, 2);

      // Generate a random interval along coord.  Intervals can't be empty.
      int ij0;
      int ij1;
      int len;
      do {
        ij0 = data.uniformInt(0, S2Projections.LIMIT_IJ);
        ij1 = data.uniformInt(ij0, S2Projections.LIMIT_IJ);
        len = ij1 - ij0;
      } while (len == 0);

      // We'll pretend we're tracking a polygon with one chain which we'll mark seen in advance.
      S2ShapeTracker tracker = new S2ShapeTracker(2, 1);
      tracker.markChain(0);
      assertTrue(tracker.finished());

      // Add the whole interval in a positive sense, shape should be unfinished.
      tracker.addInterval(face, axis, coord, ij0, ij1);
      assertFalse(tracker.finished());

      // Pick another interval of the same length but not the same endpoints.
      int ij00;
      int ij10;
      do {
        ij00 = data.uniformInt(0, S2Projections.LIMIT_IJ - len);
        ij10 = ij00 + len;
        if (ij10 >= S2Projections.LIMIT_IJ) {
          continue;
        }
      } while (ij00 == ij0 || ij10 == ij1);

      // Pick a random split point in the interval.
      int mid;
      do {
        mid = data.uniformInt(0, len) + ij00;
      } while (mid <= ij00 || mid == ij10);

      // Subtract it out in pieces.
      tracker.addInterval(face, axis, coord, ij10, mid);
      assertFalse(tracker.finished());

      tracker.addInterval(face, axis, coord, mid, ij00);
      assertFalse(tracker.finished());
    }
  }

  private static Crossing crossing(
      Boundary boundary, CrossingType type, double coord, double intercept, int edgeIndex) {
    Crossing c = new Crossing();
    c.set(boundary, type, coord, intercept, edgeIndex);
    return c;
  }
}
