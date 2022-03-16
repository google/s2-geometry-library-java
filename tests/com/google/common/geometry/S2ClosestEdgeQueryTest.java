package com.google.common.geometry;


import com.google.common.collect.Lists;

import java.util.List;
import java.util.Random;

import static com.google.common.geometry.S2.M_PI;
import static java.lang.Math.min;
import static java.lang.Math.pow;

/** Verifies {@link S2ClosestEdgeQuery}. */
public class S2ClosestEdgeQueryTest extends GeometryTestCase {

    private static final int kNumIndexes = 50;
    private static final int kNumEdges = 100;
    private static final int kNumQueries = 200;


    interface S2ShapeIndexFactory {
        void addEdges(S2Cap index_cap, int num_edges, S2ShapeIndex index);
    }

    // Generates a fractal loop that approximately fills the given S2Cap.
    class FractalLoopShapeIndexFactory implements S2ShapeIndexFactory {

        @Override
        public void addEdges(S2Cap index_cap, int num_edges, S2ShapeIndex index) {
            S2FractalBuilder builder = new S2FractalBuilder(rand);
            builder.setLevelForApproxMaxEdges(num_edges);
            index.add(builder.makeLoop(getRandomFrameAt(index_cap.axis()), index_cap.angle()));
        }

    }

    // Generates a regular loop that approximately fills the given S2Cap.
    //
    // Regular loops are nearly the worst case for distance calculations, since
    // many edges are nearly equidistant from any query point that is not
    // immediately adjacent to the loop.
    static class RegularLoopShapeIndexFactory implements S2ShapeIndexFactory {

        @Override
        public void addEdges(S2Cap index_cap, int num_edges, S2ShapeIndex index) {
            S2Loop loop = S2Loop.makeRegularLoop(index_cap.axis(), index_cap.angle(), num_edges);
            index.add(loop);
        }

    }

    // Generates a cloud of points that approximately fills the given S2Cap.
    class PointCloudShapeIndexFactory implements S2ShapeIndexFactory {

        @Override
        public void addEdges(S2Cap index_cap, int num_edges, S2ShapeIndex index) {
            List<S2Point> points = Lists.newArrayList();
            for(int i = 0; i < num_edges; ++i) {
                points.add(samplePoint(index_cap));
            }
            index.add(new S2PointVectorShape(points));
        }

    }

    public void testNoEdges() {
        S2ShapeIndex index = new S2ShapeIndex();
        S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index);
        S2Point target =new S2Point(1, 0, 0);
        S2ClosestEdgeQuery.Result edge = query.findClosestEdge(target);
        assertEquals(S1ChordAngle.INFINITY, edge.distance);
        assertNull(edge.shape);
        assertEquals(-1, edge.edgeId);
        assertFalse(edge.isInterior());
        assertTrue(edge.isEmpty());
        assertEquals(S1ChordAngle.INFINITY, query.getDistance(target));
    }

    public void testOptionsNotModified() {
        // Tests that FindClosestEdge(), GetDistance(), and IsDistanceLess() do not
        // modify query.options, even though all of these methods have their own
        // specific options requirements.
        S2ClosestEdgeQuery.Options options = new S2ClosestEdgeQuery.Options()
                .setMaxResults(3)
                .setMaxDistance(S1ChordAngle.fromS1Angle(S1Angle.degrees(3)))
                .setMaxError(S1ChordAngle.fromS1Angle(S1Angle.degrees(0.001)));
        S2ShapeIndex index = makeIndex("1:1 | 1:2 | 1:3 # #");
        S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index, options);
        S2Point target = makePoint("2:2");

        assertEquals(1, query.findClosestEdge(target).edgeId);
        assertDoubleNear(1.0, query.getDistance(target).toAngle().degrees(), 1e-15);
        assertTrue(query.isDistanceLess(target, S1ChordAngle.fromS1Angle(S1Angle.degrees(1.5))));

        // Verify that none of the options above were modified.
        assertEquals(options.getMaxResults(), query.getOptions().getMaxResults());
        assertEquals(options.getMaxDistance(), query.getOptions().getMaxDistance());
        assertEquals(options.getMaxError(), query.getOptions().getMaxError());
    }

    public void testDistanceEqualToLimit() {
        // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
        // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
        // the distance to the target exactly equals the chosen limit.
        S2Point p0 = makePoint("23:12");
        S2Point p1 = makePoint("47:11");
        List<S2Point> index_points = Lists.newArrayList(p0);
        S2ShapeIndex index = new S2ShapeIndex();
        index.add(new S2PointVectorShape(index_points));
        S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index);

        // Start with two identical points and a zero distance.
        S1ChordAngle dist0 = S1ChordAngle.ZERO;
        assertFalse(query.isDistanceLess(p0, dist0));
        assertTrue(query.isDistanceLessOrEqual(p0, dist0));
        assertTrue(query.isConservativeDistanceLessOrEqual(p0, dist0));

        // Now try two points separated by a non-zero distance.
        S1ChordAngle dist1 = new S1ChordAngle(p0, p1);
        assertFalse(query.isDistanceLess(p1, dist1));
        assertTrue(query.isDistanceLessOrEqual(p1, dist1));
        assertTrue(query.isConservativeDistanceLessOrEqual(p1, dist1));
    }

    public void testTrueDistanceLessThanS1ChordAngleDistance() {
        // Tests that IsConservativeDistanceLessOrEqual returns points where the
        // true distance is slightly less than the one computed by S1ChordAngle.
        //
        // The points below had the worst error from among 100,000 random pairs.
        S2Point p0 = new S2Point(0.78516762584829192, -0.50200400690845970, -0.36263449417782678);
        S2Point p1 = new S2Point(0.78563011732429433, -0.50187655940493503, -0.36180828883938054);

        // The S1ChordAngle distance is ~4 ulps greater than the true distance.
        S1ChordAngle dist1 = new S1ChordAngle(p0, p1);
        S1ChordAngle limit = dist1.predecessor().predecessor().predecessor().predecessor();
        assertTrue(S2Predicates.compareDistance(p0, p1, limit.getLength2()) < 0);

        // Verify that IsConservativeDistanceLessOrEqual() still returns "p1".
        List<S2Point> index_points = Lists.newArrayList(p0);
        S2ShapeIndex index = new S2ShapeIndex();
        index.add(new S2PointVectorShape(index_points));
        S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index);
        assertFalse(query.isDistanceLess(p1, limit));
        assertFalse(query.isDistanceLessOrEqual(p1, limit));
        assertTrue(query.isConservativeDistanceLessOrEqual(p1, limit));
    }

    public void testTargetPointInsideIndexedPolygon() {
        // Tests a target point in the interior of an indexed polygon.
        // (The index also includes a polyline loop with no interior.)
        S2ShapeIndex index = makeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
        S2ClosestEdgeQuery.Options options = new S2ClosestEdgeQuery.Options()
                .setIncludeInteriors(true)
                .setMaxDistance(S1Angle.degrees(1));
        S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index, options);
        S2Point target = makePoint("2:12");
        List<S2ClosestEdgeQuery.Result> results = query.findClosestEdges(target);
        assertEquals(1, results.size());
        assertEquals(S1ChordAngle.ZERO, results.get(0).distance);
        assertEquals(index.shapes.get(1), results.get(0).shape);
        assertEquals(-1, results.get(0).edgeId);
        assertTrue(results.get(0).isInterior());
        assertFalse(results.get(0).isEmpty());
    }

    public void testTargetPointOutsideIndexedPolygon() {
        // Tests a target point in the interior of a polyline loop with no
        // interior.  (The index also includes a nearby polygon.)
        S2ShapeIndex index = makeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
        S2ClosestEdgeQuery.Options options = new S2ClosestEdgeQuery.Options()
                .setIncludeInteriors(true)
                .setMaxDistance(S1Angle.degrees(1));
        S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index, options);
        S2Point target = makePoint("2:2");
        List<S2ClosestEdgeQuery.Result> results = query.findClosestEdges(target);
        assertEquals(0, results.size());
    }

    public void testIsConservativeDistanceLessOrEqual() {
        // Test
        int num_tested = 0;
        int num_conservative_needed = 0;
        for(int iter = 0; iter < 1000; ++iter) {
            rand.setSeed(iter + 1);  // Easier to reproduce a specific case.
            S2Point x = randomPoint();
            S2Point dir = randomPoint();
            S1Angle r = S1Angle.radians(M_PI * pow(1e-30, rand.nextDouble()));
            S2Point y = S2EdgeUtil.interpolateAtDistance(r, x, dir);
            S1ChordAngle limit = S1ChordAngle.fromS1Angle(r);
            if (S2Predicates.compareDistance(x, y, limit.getLength2()) <= 0) {
                S2ShapeIndex index = new S2ShapeIndex();
                index.add(new S2PointVectorShape(Lists.newArrayList(x)));
                S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(index);
                assertTrue(query.isConservativeDistanceLessOrEqual(y, limit));
                ++num_tested;
                if (!query.isDistanceLess(y, limit)) ++num_conservative_needed;
            }
        }

        System.out.println("Num tested = " + num_tested + ", Num conservative = " + num_conservative_needed);

        // Verify that in most test cases, the distance between the target points
        // was close to the desired value.  Also verify that at least in some test
        // cases, the conservative distance test was actually necessary.
        assertTrue(num_tested >= 300);
        assertTrue(num_tested <= 700);
        assertTrue(num_conservative_needed >= 20);
    }

    public void testCircleEdges() {
        testWithIndexFactory(
                new RegularLoopShapeIndexFactory(),
            kNumIndexes, kNumEdges, kNumQueries
        );
    }

    public void testFractalEdges() {
        testWithIndexFactory(
            new FractalLoopShapeIndexFactory(),
            kNumIndexes, kNumEdges, kNumQueries
        );
    }

    public void testPointCloudEdges() {
        testWithIndexFactory(
            new PointCloudShapeIndexFactory(),
            kNumIndexes, kNumEdges, kNumQueries
        );
    }

    public void testConservativeCellDistanceIsUsed() {
        // These specific test cases happen to fail if max_error() is not properly
        // taken into account when measuring distances to S2ShapeIndex cells.
        for (int seed : Lists.newArrayList(42, 681, 894, 1018, 1750, 1759, 2401)) {
            testWithIndexFactory(
                new FractalLoopShapeIndexFactory(),
                5, 100, 10,
                    seed
            );
        }
    }

    // The approximate radius of S2Cap from which query edges are chosen.
    private static final S1Angle kTestCapRadius = metersToAngle(10000);

    // An approximate bound on the distance measurement error for "reasonable"
    // distances (say, less than Pi/2) due to using S1ChordAngle.
    private static final double kTestChordAngleError = 1e-15;

    // Use "query" to find the closest edge(s) to the given target.  Verify that
    // the results satisfy the search criteria.
    public void getClosestEdges(S2MinDistanceTarget target, S2ClosestEdgeQuery query, List<S2ClosestEdgeQuery.Result> edges) {
        query.findClosestEdges(edges, target);
        assertTrue(edges.size() <= query.getOptions().getMaxResults());
        if (query.getOptions().getMaxDistance() == S1ChordAngle.INFINITY) {
            int min_expected = min(query.getOptions().getMaxResults(), S2ClosestEdgeQuery.countEdges(query.index()));
            if (!query.getOptions().isIncludeInteriors()) {
                // We can predict exactly how many edges should be returned.
                assertEquals(min_expected, edges.size());
            } else {
                // All edges should be returned, and possibly some shape interiors.
                assertTrue(edges.size() >= min_expected);
            }
        }
        for (S2ClosestEdgeQuery.Result edge : edges) {
            // Check that the edge satisfies the max_distance() condition.
            assertTrue(edge.distance.getLength2() < query.getOptions().getMaxDistance().getLength2());
        }
    }

    private S2ClosestEdgeQuery.Result testFindClosestEdges(S2MinDistanceTarget target, S2ClosestEdgeQuery query) {
        List<S2ClosestEdgeQuery.Result> expected = Lists.newArrayList();
        List<S2ClosestEdgeQuery.Result> actual = Lists.newArrayList();
        query.useBruteForce(true);
        getClosestEdges(target, query, expected);
        query.useBruteForce(false);
        getClosestEdges(target, query, actual);
        assertTrue(
            checkDistanceResults(
                expected,
                actual,
                query.getOptions().getMaxResults(),
                query.getOptions().getMaxDistance(),
                query.getOptions().getMaxError()
            )
        );

        if (expected.isEmpty()) return new S2ClosestEdgeQuery.Result(S1ChordAngle.INFINITY);

        // Note that when options.maxError > 0, expected[0].distance() may not
        // be the minimum distance.  It is never larger by more than max_error(),
        // but the actual value also depends on max_results().
        //
        // Here we verify that GetDistance() and IsDistanceLess() return results
        // that are consistent with the max_error() setting.
        S1ChordAngle max_error = query.getOptions().getMaxError();
        S1ChordAngle min_distance = expected.get(0).distance;
        assertTrue(query.getDistance(target).getLength2() <= S1ChordAngle.add(min_distance, max_error).getLength2());

        // Test IsDistanceLess().
        assertFalse(query.isDistanceLess(target, S1ChordAngle.sub(min_distance, max_error)));
        assertTrue(query.isConservativeDistanceLessOrEqual(target, min_distance));

        // Return the closest edge result so that we can also test Project.
        return expected.get(0);
    }

    // The running time of this test is proportional to
    //    (num_indexes + num_queries) * num_edges.
    // (Note that every query is checked using the brute force algorithm.)
    private void testWithIndexFactory(
            S2ShapeIndexFactory factory,
            int num_indexes,
            int num_edges,
            int num_queries
    ) {
        testWithIndexFactory(factory, num_indexes, num_edges, num_queries, 1);
    }
    private void testWithIndexFactory(
            S2ShapeIndexFactory factory,
            int num_indexes,
            int num_edges,
            int num_queries,
            int seed
        ) {
            // Build a set of MutableS2ShapeIndexes containing the desired geometry.
            List<S2Cap> index_caps = Lists.newArrayList();
            List<S2ShapeIndex> indexes = Lists.newArrayList();
            for(int i = 0 ; i < num_indexes; ++i) {
                rand.setSeed(seed + i);
                index_caps.add(S2Cap.fromAxisAngle(randomPoint(), kTestCapRadius));
                indexes.add(new S2ShapeIndex());
                factory.addEdges(index_caps.get(i), num_edges, indexes.get(i));
            }
            for (int i = 0; i <num_queries; ++i) {
                rand.setSeed(seed + i);
                int i_index = rand.nextInt(num_indexes);
                S2Cap index_cap = index_caps.get(i_index);

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                S1Angle query_radius = index_cap.angle().mul(2);
                S2Cap query_cap = S2Cap.fromAxisAngle(index_cap.axis(), query_radius);
                S2ClosestEdgeQuery query = new S2ClosestEdgeQuery(indexes.get(i_index));

                // Occasionally we don't set any limit on the number of result edges.
                // (This may return all edges if we also don't set a distance limit.)
                if (!oneIn(5)) {
                    query.getOptions().setMaxResults(1 + rand.nextInt(10));
                }
                // We set a distance limit 2/3 of the time.
                if (!oneIn(3)) {
                    query.getOptions().setMaxDistance(query_radius.mul(rand.nextDouble()));
                }
                if (oneIn(2)) {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.getOptions().setMaxError(
                        S1Angle.radians(
                            pow(
                                1e-4,
                                rand.nextDouble()
                            ) * query_radius.radians()
                        )
                    );
                }
                query.getOptions().setIncludeInteriors(oneIn(2));
                int target_type = rand.nextInt(2);

                if (target_type == 0) {
                    // Find the edges closest to a given point.
                    S2Point point = samplePoint(query_cap);
                    S2ClosestEdgeQuery.PointTarget target = new S2ClosestEdgeQuery.PointTarget(point);
                    S2ClosestEdgeQuery.Result closest = testFindClosestEdges(target, query);
                    if (!closest.distance.isInfinity()) {
                        // Also test the Project method.
                        assertDoubleNear(
                                new S1Angle(point, query.project(point, closest)).radians(),
                            closest.distance.toAngle().radians(),
                                kTestChordAngleError
                        );
                    }
                } else {
                    // Find the edges closest to a given edge.
                    S2Point a = samplePoint(query_cap);
                    S2Point b = samplePoint(
                        S2Cap.fromAxisAngle(
                            a,
                            query_radius.mul(pow(1e-4, rand.nextDouble()))
                        )
                    );
                    S2ClosestEdgeQuery.EdgeTarget target = new S2ClosestEdgeQuery.EdgeTarget(a, b);
                    testFindClosestEdges(target, query);
                }
            }
        }



    // Compare two sets of "closest" items, where "expected" is computed via brute
    // force (i.e., considering every possible candidate) and "actual" is computed
    // using a spatial data structure.  Here "max_size" is a bound on the maximum
    // number of items, "max_distance" is a limit on the distance to any item, and
    // "max_error" is the maximum error allowed when selecting which items are
    // closest (see S2ClosestEdgeQuery::Options::max_error).
    private static boolean checkDistanceResults(
            List<S2ClosestEdgeQuery.Result> expected,
            List<S2ClosestEdgeQuery.Result> actual,
            int max_size,
            S1ChordAngle max_distance,
            S1ChordAngle max_error
    ) {
        // This is a conservative bound on the error in computing the distance from
        // the target geometry to an S2Cell.  Such errors can cause candidates to be
        // pruned from the result set even though they may be slightly closer.
        S1ChordAngle kMaxPruningError = S1ChordAngle.fromS1Angle(S1Angle.radians(1e-15));
        return (checkResultSet(
                actual,
                expected,
                max_size,
                max_distance,
                max_error,
                kMaxPruningError,
                "Missing"
        ) & /*not &&*/
        checkResultSet(
                expected,
                actual,
                max_size,
                max_distance,
                max_error,
                S1ChordAngle.ZERO,
                "Extra"
        ));
    }


    // Check that result set "x" contains all the expected results from "y", and
    // does not include any duplicate results.
    private static boolean checkResultSet(
            List<S2ClosestEdgeQuery.Result> x,
            List<S2ClosestEdgeQuery.Result> y,
            int max_size,
            S1ChordAngle max_distance,
            S1ChordAngle max_error,
            S1ChordAngle max_pruning_error,
            String label
    ) {
        // Results should be sorted by distance, but not necessarily then by Id.
        for (int i = 1; i < x.size() ;++i) {
            assertTrue(x.get(i).compareTo(x.get(i - 1)) > 0);
        }

        // Result set X should contain all the items from Y whose distance is less
        // than "limit" computed below.
        S1ChordAngle limit = S1ChordAngle.ZERO;
        if (x.size() < max_size) {
            // Result set X was not limited by "max_size", so it should contain all
            // the items up to "max_distance", except that a few items right near the
            // distance limit may be missed because the distance measurements used for
            // pruning S2Cells are not conservative.
            if (max_distance == S1ChordAngle.INFINITY) {
                limit = max_distance;
            } else {
                limit = S1ChordAngle.sub(max_distance, max_pruning_error);
            }
        } else if (!x.isEmpty()) {
            // Result set X contains only the closest "max_size" items, to within a
            // tolerance of "max_error + max_pruning_error".
            limit = S1ChordAngle.sub(S1ChordAngle.sub(x.get(x.size() - 1).distance, max_error), max_pruning_error);
        }

        boolean result = true;
        for (S2ClosestEdgeQuery.Result yp : y) {
            // Note that this test also catches duplicate values.
            int count = (int) x.stream().filter(xp -> xp.equals(yp)).count();
            if (yp.distance.getLength2() < limit.getLength2() && count != 1) {
                result = false;
                System.out.println((count > 1 ? "Duplicate" : label) + " distance = ${yp.first}, id = ${yp.second}");
            }
        }

        return result;
    }

}
