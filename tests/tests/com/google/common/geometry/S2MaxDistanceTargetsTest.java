/*
 * Copyright 2022 Google Inc.
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

import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2FurthestEdgeQuery} targets. */
@RunWith(JUnit4.class)
public final class S2MaxDistanceTargetsTest extends GeometryTestCase {
  /**
   * Verifies that the overloads of updateBestDistance on {@link S2FurthestEdgeQuery.PointTarget}
   * work as expected with S1ChordAngle.maxCollector for distances from points, edges, and cells.
   */
  @Test
  public void testUpdateMaxDistanceToPointTarget() {
    S2Point targetPoint = makePointOrDie("0:0");
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(targetPoint);
    DistanceCollector<S1ChordAngle> dist0 = S1ChordAngle.maxCollector();
    dist0.set(S1ChordAngle.fromDegrees(0));
    DistanceCollector<S1ChordAngle> dist10 = S1ChordAngle.maxCollector();
    dist10.set(S1ChordAngle.fromDegrees(10));

    // Update max distance to the point target at 0:0 from another point.
    S2Point p = makePointOrDie("1:0");
    assertTrue(target.updateBestDistance(p, dist0));
    assertEquals(1.0, dist0.distance().degrees(), 1e-15);
    assertFalse(target.updateBestDistance(p, dist10));

    // Reset dist0 which was updated above.
    dist0.set(S1ChordAngle.fromDegrees(0));
    // Update max distance to the point target from an edge.
    List<S2Point> edge = parsePointsOrDie("0:-1, 0:1");
    assertTrue(target.updateBestDistance(edge.get(0), edge.get(1), dist0));
    assertEquals(1.0, dist0.distance().degrees(), 1e-15);
    assertFalse(target.updateBestDistance(edge.get(0), edge.get(1), dist10));

    // Reset dist0 again.
    dist0.set(S1ChordAngle.fromDegrees(0));
    // Update max distance to the point target from a leaf cell at 0:0. The point target is inside
    // the leaf cell, almost exactly at the center, but the max distance is to one of the corners of
    // the cell.
    S2Cell cell = new S2Cell(makePointOrDie("0:0"));
    assertTrue(target.updateBestDistance(cell, dist0));
    // The distance from the center to the corner of a leaf cell is tiny compared to 10 degrees, so
    // this will produce no update.
    assertFalse(target.updateBestDistance(cell, dist10));
  }

  /**
   * Verifies that the overloads of updateBestDistance on {@link S2FurthestEdgeQuery.EdgeTarget}
   * work as expected with S1ChordAngle.maxCollector for distances from points, edges, and cells.
   */
  @Test
  public void testUpdateMaxDistanceToEdgeTarget() {
    List<S2Point> targetEdge = parsePointsOrDie("0:-1, 0:1");
    S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.EdgeTarget<>(targetEdge.get(0), targetEdge.get(1));
    DistanceCollector<S1ChordAngle> dist0 = S1ChordAngle.maxCollector();
    dist0.set(S1ChordAngle.fromDegrees(0));
    DistanceCollector<S1ChordAngle> dist10 = S1ChordAngle.maxCollector();
    dist10.set(S1ChordAngle.fromDegrees(10));

    // Update max distance to the edge target from a point.
    S2Point p = makePointOrDie("0:2");
    assertTrue(target.updateBestDistance(p, dist0));
    assertDoubleNear(3.0, dist0.distance().degrees(), 1e-15);
    assertFalse(target.updateBestDistance(p, dist10));

    // Reset dist0 which was updated.
    dist0.set(S1ChordAngle.fromDegrees(0));
    // Update max distance to the edge target from another edge.
    List<S2Point> testEdge = parsePointsOrDie("0:2, 0:3");
    assertTrue(target.updateBestDistance(testEdge.get(0), testEdge.get(1), dist0));
    assertDoubleNear(4.0, dist0.distance().degrees(), 1e-15);
    assertFalse(target.updateBestDistance(testEdge.get(0), testEdge.get(1), dist10));

    // Reset dist0 again.
    dist0.set(S1ChordAngle.fromDegrees(0));
    // Update max distance to the edge target from a leaf cell at 0:0. The edge target intersects
    // the leaf cell, going through the center. The max distance is from a corner of the cell to a
    // target edge endpoint on the other side, just over one degree.
    S2Cell cell = new S2Cell(makePointOrDie("0:0"));
    assertTrue(target.updateBestDistance(cell, dist0));
    // But that one degree distance is small compared to dist10, so that won't update.
    assertFalse(target.updateBestDistance(cell, dist10));
  }

  /**
   * Verifies that the overloads of updateBestDistance on {@link S2FurthestEdgeQuery.CellTarget}
   * work as expected with S1ChordAngle.maxCollector for distances from points, edges, and cells.
   */
  @Test
  public void testUpdateMaxDistanceToCellTarget() {
    S2Cell targetCell = new S2Cell(makePointOrDie("0:1"));
    S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.CellTarget<>(targetCell);
    DistanceCollector<S1ChordAngle> dist0 = S1ChordAngle.maxCollector();
    dist0.set(S1ChordAngle.fromDegrees(0));
    DistanceCollector<S1ChordAngle> dist10 = S1ChordAngle.maxCollector();
    dist10.set(S1ChordAngle.fromDegrees(10));

    // Update max distance to the cell target.
    S2Point p = makePointOrDie("0:0");
    assertTrue(target.updateBestDistance(p, dist0));
    assertFalse(target.updateBestDistance(p, dist10));

    // Reset dist0 which was updated above.
    dist0.set(S1ChordAngle.fromDegrees(0));
    // Update max distance to the cell target from an edge.
    List<S2Point> edge = parsePointsOrDie("0:2, 0:3");
    assertTrue(target.updateBestDistance(edge.get(0), edge.get(1), dist0));
    assertFalse(target.updateBestDistance(edge.get(0), edge.get(1), dist10));

    // Reset dist0 again.
    dist0.set(S1ChordAngle.fromDegrees(0));
    // Update max distance to the cell target from another leaf cell at 0:0. This will be a bit more
    // than one degree.
    S2Cell cell = new S2Cell(makePointOrDie("0:0"));
    assertTrue(target.updateBestDistance(cell, dist0));
    // But that distance is small compared to 10 degrees, so this will produce no update.
    assertFalse(target.updateBestDistance(cell, dist10));
  }

  /**
   * Verifies that updateBestDistance from an edge to a {@link S2FurthestEdgeQuery.PointTarget} with
   * a {@link S1ChordAngle#maxCollector()} only returns true when the new distance from the edge to
   * the point target is greater than the old distance, not greater than or equal to.
   */
  @Test
  public void testUpdateMaxDistanceFromEdgeToPointWhenEqual() {
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(makePointOrDie("1:0"));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    List<S2Point> edge = parsePointsOrDie("0:-1, 0:1");
    assertTrue(target.updateBestDistance(edge.get(0), edge.get(1), dist));
    assertFalse(target.updateBestDistance(edge.get(0), edge.get(1), dist));
  }

  /**
   * Verifies that updateBestDistance from an edge to a {@link S2FurthestEdgeQuery.PointTarget} with
   * a {@link S1ChordAngle#maxCollector()} only returns true when the new distance from the edge to
   * the point target is greater than the old distance, not greater than or equal to.
   */
  @Test
  public void testUpdateMaxDistanceFromCellToPointWhenEqual() {
    S2Point targetPoint = makePointOrDie("1:0");
    S2FurthestEdgeQuery.PointTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.PointTarget<>(targetPoint);

    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    S2Cell cell = new S2Cell(makePointOrDie("0:0"));
    assertTrue(target.updateBestDistance(cell, dist));
    assertFalse(target.updateBestDistance(cell, dist));
  }

  /**
   * Verifies that updateBestDistance from an edge to a {@link S2FurthestEdgeQuery.EdgeTarget} with
   * a {@link S1ChordAngle#maxCollector()} only returns true when the new distance from the edge to
   * the edge target is greater than the old distance, not greater than or equal to.
   */
  @Test
  public void testUpdateMaxDistanceFromEdgeToEdgeWhenEqual() {
    S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.EdgeTarget<>(makePointOrDie("1:0"), makePointOrDie("1:1"));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    List<S2Point> edge = parsePointsOrDie("0:-1, 0:1");
    assertTrue(target.updateBestDistance(edge.get(0), edge.get(1), dist));
    assertFalse(target.updateBestDistance(edge.get(0), edge.get(1), dist));
  }

  /**
   * Verifies that updateBestDistance from a degenerate edge (a point) to an antipodal {@link
   * S2FurthestEdgeQuery.EdgeTarget} with a {@link S1ChordAngle#maxCollector()} returns true, and
   * that the distance is the maximum.
   */
  @Test
  public void testUpdateMaxDistanceFromEdgeToEdgeAntipodal() {
    S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.EdgeTarget<>(makePointOrDie("0:89"), makePointOrDie("0:91"));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    List<S2Point> edge = parsePointsOrDie("1:-90, -1:-90");
    assertTrue(target.updateBestDistance(edge.get(0), edge.get(1), dist));
    assertEquals(dist.distance(), S1ChordAngle.STRAIGHT);
  }

  /**
   * Verifies that updateBestDistance from a cell to a {@link S2FurthestEdgeQuery.EdgeTarget} with a
   * {@link S1ChordAngle#maxCollector()} only returns true when the new distance from the edge to
   * the edge target is greater than the old distance, not greater than or equal to.
   */
  @Test
  public void testUpdateMaxDistanceFromCellToEdgeWhenEqual() {
    S2FurthestEdgeQuery.EdgeTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.EdgeTarget<>(makePointOrDie("1:0"), makePointOrDie("1:1"));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    S2Cell cell = new S2Cell(makePointOrDie("0:0"));
    assertTrue(target.updateBestDistance(cell, dist));
    assertFalse(target.updateBestDistance(cell, dist));
  }

  /**
   * Verifies the consistency of {@link S2FurthestEdgeQuery.CellTarget#getCapBound()} and {@link
   * S2Cell#getMaxDistance()}.
   */
  @Test
  public void testMaxDistanceGetCapBound() {
    for (int i = 0; i < 100; i++) {
      // Random cells at uniformly distributed random levels.
      S2Cell cell = new S2Cell(data.getRandomCellId());
      // We'll get the maximum distance to these randomly located and sized cells.
      S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target =
          new S2FurthestEdgeQuery.CellTarget<>(cell);

      // Get a cap bounding the area of best possible (maximum) distance to the target cell.
      S2Cap cap = target.getCapBound();

      // Test random points from the entire surface of the sphere.
      for (int j = 0; j < 100; j++) {
        S2Point pTest = data.getRandomPoint();
        // Check that points that are outside of the cap always have a max distance to the cell
        // that is less than the maximum max distance.
        if (!cap.contains(pTest)) {
          S1ChordAngle dist = cell.getMaxDistance(pTest);
          assertTrue(dist.lessThan(S1ChordAngle.STRAIGHT));
        }
      }
    }
  }

  /**
   * Verifies that updateBestDistance from an edge to a {@link S2FurthestEdgeQuery.CellTarget} with
   * a {@link S1ChordAngle#maxCollector()} only returns true when the new distance from the edge to
   * the cell target is greater than the old distance, not greater than or equal to.
   */
  @Test
  public void testUpdateMaxDistanceFromEdgeToCellWhenEqual() {
    S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.CellTarget<>(new S2Cell(makePointOrDie("0:1")));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    List<S2Point> edge = parsePointsOrDie("0:-1, 0:1");
    assertTrue(target.updateBestDistance(edge.get(0), edge.get(1), dist));
    assertFalse(target.updateBestDistance(edge.get(0), edge.get(1), dist));
  }

  /**
   * Verifies that updateBestDistance from a cell to a {@link S2FurthestEdgeQuery.CellTarget} with a
   * {@link S1ChordAngle#maxCollector()} only returns true when the new distance from the cell to
   * the cell target is greater than the old distance, not greater than or equal to.
   */
  @Test
  public void testUpdateMaxDistanceFromCellToCellWhenEqual() {
    S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.CellTarget<>(new S2Cell(makePointOrDie("0:1")));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    S2Cell cell = new S2Cell(makePointOrDie("0:0"));
    assertTrue(target.updateBestDistance(cell, dist));
    assertFalse(target.updateBestDistance(cell, dist));
  }

  /**
   * Verifies that updateBestDistance for the distance from a cell target to an antipodal cell is
   * the maximum distance.
   */
  @Test
  public void testUpdateMaxDistanceToCellAntipodal() {
    S2Point p = makePointOrDie("0:0");
    S2FurthestEdgeQuery.CellTarget<S1ChordAngle> target =
        new S2FurthestEdgeQuery.CellTarget<>(new S2Cell(p));
    DistanceCollector<S1ChordAngle> dist = S1ChordAngle.maxCollector();
    S2Cell cell = new S2Cell(p.neg());
    assertTrue(target.updateBestDistance(cell, dist));
    assertEquals(dist.distance(), S1ChordAngle.STRAIGHT);
    // Expect a second update to do nothing.
    assertFalse(target.updateBestDistance(cell, dist));
  }
}
