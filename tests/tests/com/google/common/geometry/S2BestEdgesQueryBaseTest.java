/*
 * Copyright 2021 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static java.util.Comparator.reverseOrder;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2BestEdgesQueryBase.Result;
import com.google.common.geometry.S2ShapeUtil.PointVisitor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.Comparator;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit test for S2BestEdgesQueryBase. */
@RunWith(JUnit4.class)
public final class S2BestEdgesQueryBaseTest extends GeometryTestCase {

  static class PointTarget extends S2BestEdgesQueryBase.PointTarget<S1ChordAngle>
      implements S2BestDistanceTarget<S1ChordAngle> {
    public PointTarget(S2Point targetPoint) {
      super(targetPoint);
    }

    @Override
    public S2Cap getCapBound() {
      return S2Cap.fromAxisChord(point.neg(), S1ChordAngle.ZERO);
    }

    @Override
    public int maxBruteForceIndexSize() {
      return 2; // Arbitrary.
    }
  }

  // A target of two points, for testing result deduplication.
  static class TwoPointTarget implements S2BestDistanceTarget<S1ChordAngle> {
    protected final S2Point a;
    protected final S2Point b;

    public TwoPointTarget(S2Point a, S2Point b) {
      this.a = a;
      this.b = b;
    }

    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Point p, DistanceCollector<S1ChordAngle> collector) {
      return (collector.update(a, p) || collector.update(b, p));
    }

    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(
        S2Point v0, S2Point v1, DistanceCollector<S1ChordAngle> collector) {
      return collector.update(a, v0, v1) || collector.update(b, v0, v1);
    }

    @CanIgnoreReturnValue
    @Override
    public boolean updateBestDistance(S2Cell cell, DistanceCollector<S1ChordAngle> collector) {
      return (collector.update(a, cell) || collector.update(b, cell));
    }

    @CanIgnoreReturnValue
    @Override
    public boolean visitConnectedComponentPoints(PointVisitor visitor) {
      if (!visitor.apply(a)) {
        return false;
      }
      return visitor.apply(b);
    }

    @Override
    public S2Cap getCapBound() {
      S2Cap bound = S2Cap.fromAxisChord(a.neg(), S1ChordAngle.ZERO);
      return bound.addPoint(b.neg());
    }

    @Override
    public int maxBruteForceIndexSize() {
      return 0; // Don't use brute force unless it's requested in the options.
    }
  }

  /** A minimal implementation of a FurthestEdgeQuery for testing. */
  static class FurthestEdgeTestQuery extends S2BestEdgesQueryBase<S1ChordAngle> {
    public FurthestEdgeTestQuery(Options<S1ChordAngle> options) {
      super(options);
    }

    static class Builder extends S2BestEdgesQueryBase.Builder<S1ChordAngle> {
      public Builder(int maxResults) {
        super(S1ChordAngle.ZERO, S1ChordAngle.ZERO);
        this.maxResults = maxResults;
      }

      @Override
      public FurthestEdgeTestQuery build() {
        return new FurthestEdgeTestQuery(new Options<>(this));
      }

      @Override
      public FurthestEdgeTestQuery build(S2ShapeIndex index) {
        FurthestEdgeTestQuery query = build();
        query.init(index);
        return query;
      }
    }

    @Override
    protected DistanceCollector<S1ChordAngle> newDistanceCollector() {
      return S1ChordAngle.maxCollector();
    }

    @Override
    protected boolean atBestLimit(DistanceCollector<S1ChordAngle> distanceCollector) {
      return distanceCollector.distance().getLength2() >= S1ChordAngle.MAX_LENGTH2;
    }

    @Override
    public Comparator<S1ChordAngle> distanceComparator() {
      return reverseOrder();
    }

    @Override
    public S1ChordAngle zeroDistance() {
      return S1ChordAngle.ZERO;
    }

    @Override
    public S1ChordAngle bestDistance() {
      return S1ChordAngle.STRAIGHT;
    }

    @Override
    public S1ChordAngle worstDistance() {
      return S1ChordAngle.ZERO;
    }

    @Override
    public S1ChordAngle beyondWorstDistance() {
      return S1ChordAngle.NEGATIVE;
    }

    @Override
    public S1ChordAngle errorBoundedDistance(S1ChordAngle value) {
      return S1ChordAngle.add(value, options().maxError);
    }

    @Override
    public S1ChordAngle searchCapRadius(S1ChordAngle antipodalCapRadius, S1ChordAngle minDistance) {
      return S1ChordAngle.add(
          antipodalCapRadius, S1ChordAngle.sub(S1ChordAngle.STRAIGHT, minDistance));
    }

    @Override
    protected boolean visitBestDistanceContainingShapes(
        S2BestDistanceTarget<S1ChordAngle> target, S2ContainsPointQuery.ShapeVisitor visitor) {
      S2ContainsPointQuery containsPointQuery = new S2ContainsPointQuery(index);
      return target.visitConnectedComponentPoints(
          targetPoint -> containsPointQuery.visitContainingShapes(targetPoint.neg(), visitor));
    }

    FurthestEdgeTestQuery(S2ShapeIndex index, Options<S1ChordAngle> options) {
      super(options);
      init(index);
    }

    public List<Result<S1ChordAngle>> findFurthestEdges(S2BestDistanceTarget<S1ChordAngle> target) {
      return findBestEdges(target);
    }

    public void findFurthestEdges(
        S2BestDistanceTarget<S1ChordAngle> target,
        S2BestEdgesQueryBase.ResultVisitor<S1ChordAngle> visitor) {
      findBestEdges(target, visitor);
    }
  }

  // Track the number of visits the ResultVisitor gets.
  private int visitCount;

  @Test
  public void testMaxDistance() {
    // An index with four points. This will be a single shape with four edges.
    S2ShapeIndex index = makeIndexOrDie("0:0 | 1:0 | 2:0 | 3:0 # #");
    assertEquals(1, index.getShapes().size());

    // A target of a single point. The furthest point in the index is at 0:0, distance 4.
    PointTarget target = new PointTarget(S2TextFormat.makePointOrDie("4:0"));

    // A test query to get the single furthest result.
    FurthestEdgeTestQuery.Builder builder = new FurthestEdgeTestQuery.Builder(1);
    FurthestEdgeTestQuery query = builder.build(index);

    // Find the furthest single result.
    List<Result<S1ChordAngle>> results = query.findFurthestEdges(target);

    // The single result should be the furthest index point: edge 0 of the single shape, at 0:0.
    assertEquals(1, results.size());
    Result<S1ChordAngle> result = results.get(0);
    assertEquals(0, result.shapeId());
    assertEquals(0, result.edgeId());
    assertDoubleNear(4, result.distance().degrees(), 1e-13);

    // Running a similar query but with maxResults = 10 and obtaining results with a visitor that
    // returns false on the third visit should get three of the four possible results.
    builder.setMaxResults(10);
    FurthestEdgeTestQuery query2 = builder.build(index);

    visitCount = 0;
    query2.findFurthestEdges(
        target,
        (distance, s, edgeId) -> {
          visitCount++;

          // Just check the distance against the edge id to ensure they're valid.
          if (result.edgeId() == 0) {
            assertDoubleNear(4, result.distance().degrees(), 1e-13);
          } else if (result.edgeId() == 1) {
            assertDoubleNear(3, result.distance().degrees(), 1e-13);
          } else if (result.edgeId() == 2) {
            assertDoubleNear(2, result.distance().degrees(), 1e-13);
          } else if (result.edgeId() == 3) {
            assertDoubleNear(1, result.distance().degrees(), 1e-13);
          }

          return visitCount < 3;
        });

    assertEquals(3, visitCount);
  }

  @Test
  public void testShapeInteriorsAreDeduplicated() {
    // An index with one polygon.
    S2ShapeIndex index = makeIndexOrDie("# # 0:0, 0:5, 5:5, 5:0");

    // Get up to 8 results with no minimum distance.
    FurthestEdgeTestQuery.Builder builder = new FurthestEdgeTestQuery.Builder(8);
    FurthestEdgeTestQuery query = builder.build(index);

    // A single target of two points that are antipodal to the indexed polygon
    TwoPointTarget target =
        new TwoPointTarget(
            S2TextFormat.makePointOrDie("1:1").neg(), S2TextFormat.makePointOrDie("2:2").neg());

    // Get up to 8 results. If shape interiors are not correctly deduplicated, the two points of the
    // target will produce two interior results. There should only be one interior result, but there
    // will also be results for the four polygon edges.
    List<Result<S1ChordAngle>> results = query.findFurthestEdges(target);
    int interiorResults = 0;
    int edgeResults = 0;

    for (Result<S1ChordAngle> result : results) {
      assertEquals(0, result.shapeId());
      if (result.isInterior()) {
        interiorResults++;
        assertEquals(S1ChordAngle.STRAIGHT, result.distance());
        assertEquals(-1, result.edgeId());
      } else {
        edgeResults++;
        assertTrue(S1ChordAngle.STRAIGHT.greaterThan(result.distance()));
        assertTrue(result.edgeId() >= 0);
      }
    }

    assertEquals(1, interiorResults);
    assertEquals(4, edgeResults);
  }
}
