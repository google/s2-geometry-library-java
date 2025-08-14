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

import static com.google.common.base.Preconditions.checkState;
import static com.google.common.geometry.S2CellId.MAX_LEVEL;
import static com.google.common.geometry.S2CellId.begin;
import static com.google.common.geometry.S2CellId.end;
import static com.google.common.geometry.S2CellId.fromFacePosLevel;
import static com.google.common.geometry.S2DensityClusterQueryTest.EditableData.WIDTH;
import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;
import static java.lang.Math.toIntExact;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.geometry.S2DensityClusterQuery.CellInterpolator;
import com.google.common.geometry.S2DensityClusterQuery.Cluster;
import com.google.common.geometry.S2DensityClusterQueryTest.EditableData.Feature;
import com.google.common.geometry.S2DensityClusterQueryTest.EditableData.ST;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder.CellWeightFunction;
import com.google.common.geometry.S2DensityTree.DecodedPath;
import com.google.common.geometry.S2DensityTree.TreeEncoder;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.function.Function;
import java.util.function.Supplier;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link S2DensityClusterQuery}. */
@RunWith(JUnit4.class)
public final class S2DensityClusterQueryTest extends GeometryTestCase {
  // Set clustering in the scenario where clusters one or more whole face cells.
  @Test
  public void testClustersAsFaces() {
    ImmutableMap<S2CellId, Long> leaves =
        ImmutableMap.of(
            S2CellId.fromFace(0), 100L,
            S2CellId.fromFace(1), 100L,
            S2CellId.fromFace(2), 100L,
            S2CellId.fromFace(3), 100L,
            S2CellId.fromFace(4), 100L,
            S2CellId.fromFace(5), 100L);

    // Build the tree
    S2DensityTree tree = density(leaves);

    // Cluster the tree
    S2DensityClusterQuery query100 = new S2DensityClusterQuery(100);
    List<S2CellUnion> wt100clusters = query100.tightCoverings(tree);

    // Verify clusters are each a single face.
    assertEquals(6, wt100clusters.size());
    for (S2CellUnion cluster : wt100clusters) {
      assertEquals(1, cluster.cellIds().size());
      S2CellId cell = cluster.cellIds().get(0);
      assertTrue(cell.isFace());
    }

    // If the clusters have double the weight, there should be only three, each covering two faces.
    S2DensityClusterQuery query200 = new S2DensityClusterQuery(200);
    List<S2CellUnion> wt200clusters = query200.tightCoverings(tree);

    // Verify clusters are each two faces.
    assertEquals(3, wt200clusters.size());
    for (S2CellUnion cluster : wt200clusters) {
      assertEquals(2, cluster.size());
      S2CellId cell0 = cluster.cellIds().get(0);
      assertTrue(cell0.isFace());
      S2CellId cell1 = cluster.cellIds().get(0);
      assertTrue(cell1.isFace());
    }
  }

  // Check clustering that requires descending into children but not subdividing tree cells.
  @Test
  public void testClustersDivideFaces() {
    ImmutableMap<S2CellId, Long> leaves =
        ImmutableMap.of(
            // Face 0 weight divided between two children
            S2CellId.fromFace(0).child(1), 50L,
            S2CellId.fromFace(0).child(2), 50L,
            // Face 1 divided between four grandchildren
            S2CellId.fromFace(1).child(0).child(0), 25L,
            S2CellId.fromFace(1).child(0).child(1), 25L,
            S2CellId.fromFace(1).child(0).child(2), 25L,
            S2CellId.fromFace(1).child(2).child(0), 25L);

    S2DensityTree tree = density(leaves);

    // Test that tree.normalize() and creating a new DecodedPath doesn't crash.
    DecodedPath unused = new DecodedPath(tree.normalize());

    // If divided into clusters of weight 50, faces 0 and 1 will both be split into two clusters.
    S2DensityClusterQuery query50 = new S2DensityClusterQuery(50);
    List<S2CellUnion> clusters50 = query50.tightCoverings(tree);

    // There should be a total of four clusters.
    assertEquals(4, clusters50.size());

    S2CellUnion cluster = clusters50.get(0);
    assertEquals(1, cluster.size());
    assertEquals(S2CellId.fromFace(0).child(1), cluster.cellId(0));

    cluster = clusters50.get(1);
    assertEquals(1, cluster.size());
    assertEquals(S2CellId.fromFace(0).child(2), cluster.cellId(0));

    cluster = clusters50.get(2);
    assertEquals(2, cluster.size());
    assertEquals(S2CellId.fromFace(1).child(0).child(0), cluster.cellId(0));
    assertEquals(S2CellId.fromFace(1).child(0).child(1), cluster.cellId(1));

    cluster = clusters50.get(3);
    assertEquals(2, cluster.size());
    // The cluster contains the last two weight 25 tree nodes.
    assertEquals(S2CellId.fromFace(1).child(0).child(2), cluster.cellId(0));
    assertEquals(S2CellId.fromFace(1).child(2).child(0), cluster.cellId(1));
  }

  @Test
  public void testUninterpolate() {
    S2CellId face = S2CellId.fromFace(0);
    S2CellId child1 = face.child(1);
    assertAlmostEquals(0.25, new CellInterpolator(face).uninterpolate(child1));
  }

  @Test
  public void testNormalizingTree() {
    ImmutableMap<S2CellId, Long> leaves =
        ImmutableMap.ofEntries(
            // In this scenario, there are a total of 8 features each of weight 100 intersecting
            // face 0.
            Map.entry(S2CellId.fromFace(0), 800L),

            // In child 0, four features, each intersect all four children of (face 0, child 0).
            Map.entry(S2CellId.fromFace(0).child(0), 400L),
            Map.entry(S2CellId.fromFace(0).child(0).child(0), 400L),
            Map.entry(S2CellId.fromFace(0).child(0).child(1), 400L),
            Map.entry(S2CellId.fromFace(0).child(0).child(2), 400L),
            Map.entry(S2CellId.fromFace(0).child(0).child(3), 400L),

            // In child 2, four features, each intersect one (different) child of (face 0, child 2).
            Map.entry(S2CellId.fromFace(0).child(2), 400L),
            Map.entry(S2CellId.fromFace(0).child(2).child(0), 100L),
            Map.entry(S2CellId.fromFace(0).child(2).child(1), 100L),
            Map.entry(S2CellId.fromFace(0).child(2).child(2), 100L),
            Map.entry(S2CellId.fromFace(0).child(2).child(3), 100L));

    // Note that density() is not used here. We are deliberately creating a (valid) S2DensityTree
    // where deeper nodes in the tree have redundant weight.
    TreeEncoder encoder = new TreeEncoder();
    leaves.forEach(encoder::put);
    S2DensityTree tree = encoder.build();

    // Create clusters of weight 400 using the query.
    S2DensityClusterQuery query400 = new S2DensityClusterQuery(400);
    List<S2CellUnion> clusters400 = query400.tightCoverings(tree);

    // If S2DensityClusterQuery does not call normalizedTree(), there would be five clusters of
    // weight 400: one for each of the four children of (face 0, child 0), and one for (face 0,
    // child 2). But, because normalizedTree removes the redundant weight, there will be two
    // clusters.
    assertEquals(2, clusters400.size());

    // Note that although clusters are formed from leaves of the density tree, they are normalized
    // S2CellUnions, so each set of (child0, child1, child2, child3) is replaced by a single cell.
    S2CellUnion cluster0 = clusters400.get(0);
    assertEquals(1, cluster0.size());
    assertEquals(S2CellId.fromFace(0).child(0), cluster0.cellId(0));

    S2CellUnion cluster1 = clusters400.get(1);
    assertEquals(1, cluster1.size());
    assertEquals(S2CellId.fromFace(0).child(2), cluster1.cellId(0));
  }

  // Test the subdivision of cells using fractions of rangeMin and rangeMax.
  public void checkInterpolationOnFace(int face) {
    // A tree with a single face cell of weight 100.
    ImmutableMap<S2CellId, Long> leaves =
        ImmutableMap.of(
            S2CellId.fromFace(face), 100L);

    // Build the tree
    S2DensityTree tree = density(leaves);

    // Cluster with size 50, which will split the face into two by interpolation.
    S2DensityClusterQuery query50 = new S2DensityClusterQuery(50);
    List<S2CellUnion> clusters50 = query50.tightCoverings(tree);

    // There should be two clusters, each with two level 1 cells.
    assertEquals(2, clusters50.size());

    S2CellUnion cluster0 = clusters50.get(0);
    assertEquals(2, cluster0.size());
    assertEquals(S2CellId.fromFace(face).child(0), cluster0.cellId(0));
    assertEquals(S2CellId.fromFace(face).child(1), cluster0.cellId(1));

    S2CellUnion cluster1 = clusters50.get(1);
    assertEquals(2, cluster1.size());
    assertEquals(S2CellId.fromFace(face).child(2), cluster1.cellId(0));
    assertEquals(S2CellId.fromFace(face).child(3), cluster1.cellId(1));

    // Cluster with size 32. This should produce three clusters, where the last cluster is slightly
    // larger but within the acceptable 20% tolerance.
    S2DensityClusterQuery query32 = new S2DensityClusterQuery(32);
    List<S2CellUnion> clusters32 = query32.tightCoverings(tree);
    assertEquals(3, clusters32.size());
  }

  @Test
  public void testInterpolationOnFaces() {
    checkInterpolationOnFace(0);
    checkInterpolationOnFace(1);
    checkInterpolationOnFace(4);
  }

  // Tests accessing tree nodes with DecodedPath.
  @Test
  public void testDecodedPath() {
    ImmutableMap<S2CellId, Long> leaves =
        ImmutableMap.of(
            // Face 2 weight divided between two children
            S2CellId.fromFace(2).child(1), 50L,
            S2CellId.fromFace(2).child(2), 50L,
            // Face 4 divided between four grandchildren
            S2CellId.fromFace(4).child(0).child(0), 25L,
            S2CellId.fromFace(4).child(2).child(1), 25L,
            S2CellId.fromFace(4).child(2).child(2), 25L,
            S2CellId.fromFace(4).child(3).child(0), 25L);

    S2DensityTree tree = density(leaves);

    // Set up a DecodedPath for the tree.
    S2DensityTree.DecodedPath decoder = new S2DensityTree.DecodedPath(tree);

    // Check that the leaves can be obtained by DecodedPath.
    for (Map.Entry<S2CellId, Long> leaf : leaves.entrySet()) {
      DecodedPath freshDecoder = new DecodedPath(tree);
      // Lookup in a newly created decoder
      S2DensityTree.Cell cell = freshDecoder.cell(leaf.getKey());
      assertFalse(cell.isEmpty());
      assertEquals(cell.weight, (long) leaf.getValue());
      // Lookup in an existing decoder
      cell = decoder.cell(leaf.getKey());
      assertFalse(cell.isEmpty());
      assertEquals(cell.weight, (long) leaf.getValue());
    }

    // Face nodes not in the tree should return an empty cell.
    assertTrue(decoder.cell(S2CellId.fromFace(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(1)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(3)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(5)).isEmpty());

    // Face 2 should be present with the sum of the leaf weights.
    assertEquals(100, decoder.cell(S2CellId.fromFace(2)).weight);

    // Level 1 nodes of a face that isn't in the tree aren't in the tree either.
    assertTrue(decoder.cell(S2CellId.fromFace(1).child(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(1).child(1)).isEmpty());

    // Level 1 nodes of face 2 that is in the tree
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(0)).isEmpty());
    assertEquals(0, decoder.cell(S2CellId.fromFace(2).child(0)).weight);
    assertFalse(decoder.cell(S2CellId.fromFace(2).child(1)).isEmpty());
    assertEquals(50, decoder.cell(S2CellId.fromFace(2).child(1)).weight);
    assertFalse(decoder.cell(S2CellId.fromFace(2).child(2)).isEmpty());
    assertEquals(50, decoder.cell(S2CellId.fromFace(2).child(2)).weight);
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(3)).isEmpty());
    assertEquals(0, decoder.cell(S2CellId.fromFace(2).child(3)).weight);
  }

  private S2DensityTree buildDensityTree(S2Polygon polygon, int treeSizeBytes, int maxLevel) {
    S2ShapeIndex index = new S2ShapeIndex();
    IdentityHashMap<S2Shape, S2Polygon> features = new IdentityHashMap<>();
    S2Shape shape = polygon.shape();
    index.add(shape);
    features.put(shape, polygon);
    S2DensityTree result = S2DensityTree.featureDensity(
        treeSizeBytes, maxLevel, index, features::get, weight -> 1);
    return result;
  }

  private final S2RegionSharder sharder(String stageName, S2DensityTree tree, int clusterWeight) {
    System.err.println(
        Strings.lenientFormat("%s/InputTreeNumLeaves=%s", stageName, tree.getLeaves().size()));
    System.err.println(
        Strings.lenientFormat("%s/ClusterTargetWeight=%s", stageName, clusterWeight));

    List<S2CellUnion> partitions = new S2DensityClusterQuery(clusterWeight).tightCoverings(tree);
    System.err.println(Strings.lenientFormat("%s/ShardCount=%s", stageName, partitions.size()));

    // Shards should never be empty.
    for (S2CellUnion partition : partitions) {
      assertFalse(partition.isEmpty());
    }

    return new S2RegionSharder(partitions);
  }

  /**
   * Reproduces a bug that used to occur when feature weights are low, like the default 1, resulting
   * normalizedTree setting normalized weights to zero.
   */
  @Test
  public void testClustering() {
    // Values from GeoDensityOp & GeoClusterOp
    final int defaultTreeSizeBytes = 1 << 20;
    final int defaultMaxLevel = 20;
    final int defaultClusterWeight = 64 * (1 << 20);

    Pair<S2Polygon, S2Polygon> polygons = intersectingPolygons();
    S2DensityTree leftTree =
        buildDensityTree(polygons.getFirst(), defaultTreeSizeBytes, defaultMaxLevel);
    S2DensityTree rightTree =
        buildDensityTree(polygons.getSecond(), defaultTreeSizeBytes, defaultMaxLevel);
    S2DensityTree sumDensity =
        S2DensityTree.sumDensity(
            ImmutableList.of(leftTree, rightTree), defaultTreeSizeBytes, defaultMaxLevel);

    S2RegionSharder sharder = sharder("testClustering", sumDensity, defaultClusterWeight);

    // Each polygon should intersect the single cluster.
    ImmutableList<Integer> intersectingShardsLeft = ImmutableList.copyOf(
        sharder.getIntersectingShards(polygons.getFirst()));
    assertEquals(1, intersectingShardsLeft.size());
    ImmutableList<Integer> intersectingShardsRight= ImmutableList.copyOf(
        sharder.getIntersectingShards(polygons.getSecond()));
    assertEquals(1, intersectingShardsRight.size());
  }

  /** Returns a pair of polygons created by buffering intersecting lines. */
  private Pair<S2Polygon, S2Polygon> intersectingPolygons() {
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxLevel(22).setMaxCells(1000).build();

    // These lines intersect.
    S2Polyline leftLine =
        new S2Polyline(
            ImmutableList.of(
                S2LatLng.fromDegrees(-7.285604594541237E-4, 1.4571098021440003E-4).toPoint(),
                S2LatLng.fromDegrees(-7.285604593975788E-4, 7.285604594564796E-4).toPoint()));
    S2Polyline rightLine =
        new S2Polyline(
            ImmutableList.of(
                S2LatLng.fromDegrees(-4.371346081495913E-4, 4.371346081623137E-4).toPoint(),
                S2LatLng.fromDegrees(-0.0010199885340521325, 4.371346081623137E-4).toPoint()));

    S2Polygon leftPolygon =
        S2Polygon.fromCellUnionBorder(
            coverer.getCovering(
                new S2ShapeIndexBufferedRegion(
                    S2ShapeIndex.fromShapes(leftLine), S2Earth.metersToAngle(1.0))));
    S2Polygon rightPolygon =
        S2Polygon.fromCellUnionBorder(
            coverer.getCovering(
                new S2ShapeIndexBufferedRegion(
                    S2ShapeIndex.fromShapes(rightLine), S2Earth.metersToAngle(1.0))));
    return Pair.of(leftPolygon, rightPolygon);
  }

  /** A weight of one shouldn't degenerate to weights of zero internally. */
  @Test
  public void testMinClusterBoundsAreNonZero() {
    S2DensityClusterQuery query = new S2DensityClusterQuery(1);
    assertEquals(1, query.minClusterWeight());
    assertEquals(1, query.maxClusterWeight());
  }

  @Test public void testClamp() {
    S2CellId face = S2CellId.fromFace(1);
    S2CellId begin = begin(MAX_LEVEL);
    S2CellId end = end(MAX_LEVEL);
    S2DensityTree density = density(ImmutableMap.of(
        face.rangeMin(), 10L,
        face.child(1).child(2), 11L,
        face.child(2).child(1), 12L,
        face.rangeMax(), 13L,
        end.prev(), 14L));
    assertRange(density, face.rangeMin(), end, begin, end);
    assertRange(density, end.prev(), end, end.prev(), end);
    assertRange(density, begin, begin, begin, face.rangeMin());
    assertRange(density, end, end, end, end);
    S2CellId face2 = face.next();
    assertRange(density, face2.rangeMin(), face2.rangeMin(), face2, face2.next());
    assertRange(
        density(ImmutableMap.of(begin, 1L)),
        begin.next(), begin.next(),
        begin.next(), end); // Doesn't miss the density until level 30
  }

  @Test public void testCluster() {
    S2CellId parent = S2CellId.fromFace(1).child(2);
    long weight = 10;
    S2DensityTree density = density(ImmutableMap.of(parent, weight));
    for (int level = 1; level < 4; level++) {
      double fraction = 1.0 / (1 << (2 * level));
      long wc = round(ceil(fraction * weight));
      for (S2CellId c = parent.childBegin(parent.level() + level);
          c.rangeMin().lessThan(parent.rangeMax());
          c = c.next()) {
        S2CellId d = c.next().rangeMin();
        assertEquals(wc, S2DensityClusterQuery.cluster(density, c.rangeMin(), d).weight());
      }
    }
    assertEquals(0L, S2DensityClusterQuery.cluster(density, begin(30), parent.rangeMin()).weight());
    assertEquals(10L, S2DensityClusterQuery.cluster(density, begin(30), end(30)).weight());
  }

  @Test public void testRegionKeys() {
    S2CellId a = begin(30).advance(10);
    S2CellId b = begin(30).advance(20);
    S2CellId c = begin(30).advance(30);
    S2CellId d = begin(30).advance(40);
    Function<S2Region, String> sharder = S2DensityClusterQuery.regionKeys(
        ImmutableMap.of(
            "Foo", new Cluster(a, b, 10),
            "Bar", new Cluster(c, d, 20)),
        "Baz");
    assertEquals("Baz", sharder.apply(new S2Cell(a.prev())));
    assertEquals("Foo", sharder.apply(new S2Cell(a)));
    assertEquals("Foo", sharder.apply(new S2Cell(b.prev()))); // prev gets inside half-open range
    assertEquals("Baz", sharder.apply(new S2Cell(c.prev())));
    assertEquals("Bar", sharder.apply(new S2Cell(c)));
    assertEquals("Bar", sharder.apply(new S2Cell(d.prev()))); // prev gets inside half-open range
    assertEquals("Baz", sharder.apply(new S2Cell(begin(30))));
    assertEquals("Baz", sharder.apply(new S2Cell(end(30).prev()))); // prev gets inside valid range
  }

  @Test public void testEdit() {
    EditableData data = new EditableData(4, 6);
    // No records results in no clusters.
    data.assertEdits("");
    data.reset();
    // One record results in one cluster.
    data.assertEdits("0+0:1 -> 0:1,0:3=1", 0, 0, 1);
    data.reset();
    // Insert a cluster of 1 record, remove the record, and edit should remove the cluster.
    data.assertEdits("0+0:1 -> 0:1,0:3=1", 0, 0, 1);
    data.records.clear();
    data.assertEdits("0+0:1 -> null");
    data.reset();
    // Add a cluster, then grow it to the max size. It almost splits.
    data.assertEdits("0+0:1 -> 0:1,0:7=2", 0, 1, 1, 0, 5, 1);
    data.assertEdits("0+0:1 -> 0:1,0:7=6", 0, 3, 4);
    data.reset();
    // Add a cluster of 2 records, then add a 3rd so it does split.
    data.assertEdits("0+0:1 -> 0:1,0:7=2", 0, 1, 1, 0, 5, 1);
    data.assertEdits("0+0:1 -> 0:1,0:5=6 | 1+0:5 -> 0:5,0:7=1", 0, 3, 5);
    data.reset();
    // Add a cluster of 1 record, and then 2 clusters of 1 record before it.
    data.assertEdits("0+0:5 -> 0:5,0:7=4", 0, 5, 4);
    data.assertEdits("1+0:1 -> 0:1,0:3=5 | 1+0:3 -> 0:3,0:5=5", 0, 1, 5, 0, 3, 5);
    data.reset();
    // Add a cluster of 1 record, and then 2 clusters of 1 record after it.
    data.assertEdits("0+0:1 -> 0:1,0:3=5", 0, 1, 5);
    data.assertEdits("1+0:3 -> 0:3,0:5=5 | 1+0:5 -> 0:5,0:7=5", 0, 3, 5, 0, 5, 5);
    data.reset();
  }

  /**
   * A more extended test that verifies the process ramps up and ramps down, and exercises the most
   * common edge cases with a variety of distributions. Each distribution adds and removes a random
   * number of features until the collection has 10k features, and then adds and removes a random
   * number of features until the collection has 0 features. When ramping up, the added features are
   * on average greater than the removed features, and vice-versa when ramping down. The boundary
   * cases these tests hit include:
   * <ul>
   * <li>Cluster edits to empty tables, and edits that create empty tables.
   * <li>Repeated splits and merges over time.
   * <li>Density tree leaf cells that are larger than the desired cluster weight.
   * <li>Density tree leaves that are far smaller than the desired cluster weight
   * <li>Cluster edits that "steal" from the head of the next cluster or gap.
   * <li>Gaps that are large enough to get assigned to both neighboring clusters.
   * </ul>
   *
   * <p>This test also prints a CSV file on stdout that can be pasted into a spreadsheet to graph
   * the clusters as a function of time. It provides the following columns:
   * <ul>
   * <li>distribution: the sampling distribution used to place new features
   * <li>time: the integer time in the current distribution
   * <li>count: the number of clusters at this time
   * <li>sum: the number of features at this time
   * <li>min: the minimum number of features in any cluster
   * <li>max: the maximum number of features in any cluster
   * <li>mean: the average number of features across the clusters
   * <li>stdev: the standard deviation of the number of features across the clusters
   * <li>hilbert: a graph of the hilbert curve, where empty Hilbert regions print '.' and clustered
   * Hilbert regions print a single character to indicate the order in which the cluster in that
   * region was introduced, counting from 0-9 and then A-Z. It's only able to reliably show clusters
   * that span more than 1/64th of the Hilbert range.
   * </ul>
   */
  @Test
  public void testClusterEditDistributions() {
    EditableData data = new EditableData(800, 1200);
    S2CellId begin = S2CellId.begin(30); // start of hilbert curve
    data.printHeader("distribution", System.out);
    runDistribution(data, "Uniform", () -> begin.advance(nextLong(data.r, WIDTH)));
    runDistribution(data, "GaussianCentered",
        () -> begin.advance(WIDTH / 2 + round(data.r.nextGaussian() / 20 * WIDTH)));
    runDistribution(data, "GaussianAtBegin",
        () -> begin.advance(abs(round(data.r.nextGaussian() / 10 * WIDTH))));
    // Note the skewed distribution requires a much larger density tree, so the clusters are trash.
    runDistribution(data, "Skewed",
        () -> begin.advance(data.r.nextLong() & ((1L << nextLong(data.r, 60)) - 1)));
  }

  /** Returns a random long with at most the given magnitude. J2CL does not have this method. */
  private static long nextLong(Random r, long max) {
    return abs(r.nextLong()) % max;
  }

  private void runDistribution(EditableData data, String label, Supplier<S2CellId> distribution) {
    data.reset();
    Map<ST, Character> ids = new HashMap<>(); // unique IDs for each cluster, easy to print
    boolean rampUp = true;
    for (data.time = 0; data.time == 0 || !data.records.isEmpty(); data.time++) {
      rampUp &= data.records.size() < 10_000;
      int numRemoved = rampUp ? data.r.nextInt(10) : data.r.nextInt(100);
      int numAdded = rampUp ? max(data.r.nextInt(100), numRemoved) : data.r.nextInt(20);
      for (int feature = 0; feature < numAdded; feature++) {
        data.records.add(new Feature(distribution.get(), 1));
      }
      Collections.shuffle(data.records, data.r);
      data.records.subList(max(0, data.records.size() - numRemoved), data.records.size()).clear();
      data.editClusters();
      if (!data.clusters.isEmpty()) {
        data.printStats(label, ids, System.out);
      }
    }
    assertEquals(ImmutableMap.of(), data.clusters);
  }

  static class EditableData {
    /** The width of the hilbert curve in units we can pass to {@link S2CellId#advance}. */
    public static final long WIDTH = S2CellId.end(30).distanceFromBegin();

    private final S2DensityClusterQuery query;
    private final Random r = new Random(0); // a shared source of randomness if needed
    private final Map<ST, Cluster> clusters = new HashMap<>();
    private final List<Feature> records = new ArrayList<>();
    private int time;

    EditableData(int min, int max) {
      this.query = new S2DensityClusterQuery(min, max);
    }

    private void reset() {
      clusters.clear();
      records.clear();
      time = 0;
    }

    /** Edits 'clusters' with the new density of 'records'. */
    private void editClusters() {
      // Use just enough bytes for the tree that we can create the expected clusters in each test.
      // Making this value larger would not change the result, and making it smaller would only make
      // for incorrect clusters.
      BreadthFirstTreeBuilder builder = new BreadthFirstTreeBuilder(400, MAX_LEVEL);
      S2CellIndex index = new S2CellIndex();
      for (int i = 0; i < records.size(); i++) {
        index.add(records.get(i).cell, i);
      }
      index.build();
      S2CellUnion union = new S2CellUnion();
      CellWeightFunction weigher = id -> {
        union.cellIds().clear();
        union.cellIds().add(id);
        long sum = 0;
        for (int label : index.getIntersectingLabels(union)) {
          sum += records.get(label).weight;
        }
        return sum;
      };
      S2DensityTree density = builder.build(weigher);
      query.edit(clusters, density, c -> new ST(time, c.begin()));
    }

    /** Checks that cluster edits are expected after adding tuples of face,pos,weight. */
    private void assertEdits(String expectedStrings, long ... tuples) {
      for (int i = 0; i < tuples.length; i += 3) {
        S2CellId cell = fromFacePosLevel(toIntExact(tuples[i]), tuples[i + 1], MAX_LEVEL);
        records.add(new Feature(cell, tuples[i + 2]));
      }
      Map<ST, Cluster> copy = new HashMap<>(clusters);
      editClusters();
      List<String> strings = new ArrayList<>();
      for (Entry<ST, Cluster> cluster : clusters.entrySet()) {
        if (!cluster.getValue().equals(copy.get(cluster.getKey()))) {
          strings.add(cluster.getKey() + " -> " + cluster.getValue());
        }
      }
      for (Entry<ST, Cluster> cluster : copy.entrySet()) {
        if (!clusters.containsKey(cluster.getKey())) {
          strings.add(cluster.getKey() + " -> " + null);
        }
      }
      strings.sort(String::compareTo);
      assertEquals(expectedStrings, Joiner.on(" | ").join(strings));
      time++;
    }

    /** Prints the CSV stats header to 'out'. */
    public void printHeader(String label, PrintStream out) {
      out.println(label + ",time,count,sum,min,max,mean,stdev,hilbert");
    }

    /** Prints the current cluster stats at 'time' to the given output. */
    public void printStats(String runLabel, Map<ST, Character> clusterLabels, PrintStream out) {
      // Compute the histogram of cluster positions along the hilbert curve.
      List<Entry<ST, Cluster>> sorted = new ArrayList<>(clusters.entrySet());
      sorted.sort(Comparator.comparing(e -> e.getValue().begin()));
      int n = 0;
      char[] histogram = new char[64];
      double scale = 1.0 * histogram.length / WIDTH;
      for (var e : sorted) {
        Cluster cluster = e.getValue();
        int a = (int) round(floor(scale * cluster.begin().distanceFromBegin()));
        int b = (int) round(floor(scale * cluster.end().distanceFromBegin()));
        Arrays.fill(histogram, n, a, '.');
        Arrays.fill(histogram, a, b, clusterLabels.computeIfAbsent(e.getKey(), k -> {
          String base32 = Integer.toString(clusterLabels.size() + 1, 32);
          checkState(base32.length() == 1);
          return base32.charAt(0);
        }));
        n = b;
      }
      Arrays.fill(histogram, n, histogram.length, '.');

      // Compute stats on the clusters.
      long[] counts = new long[clusters.size()];
      long sum = 0;
      long min = Long.MAX_VALUE;
      long max = Long.MIN_VALUE;
      for (int i = 0; i < sorted.size(); i++) {
        Cluster cluster = sorted.get(i).getValue();
        counts[i] = records.stream().filter(f -> f.intersects(cluster)).count();
        sum += counts[i];
        min = min(min, counts[i]);
        max = max(max, counts[i]);
      }
      double mean = 1.0 * sum / sorted.size();
      double diff2 = 0;
      for (long count : counts) {
        double diff = count - mean;
        diff2 += diff * diff;
      }
      double stdev = sqrt(diff2 / sorted.size());

      // Print the csv for this increment.
      out.println(Joiner.on(',').join(ImmutableList.of(runLabel, time, sorted.size(),
          sum, min, max, round(mean), round(stdev), new String(histogram))));
    }

    /** A simple space+time key. */
    static final class ST {
      private final long time;
      private final S2CellId cell;

      /** Creates a new ST key. */
      public ST(long time, S2CellId cell) {
        this.time = time;
        this.cell = cell;
      }

      @Override public int hashCode() {
        return cell.hashCode() ^ (int) time;
      }

      @Override public boolean equals(Object other) {
        if (other instanceof ST) {
          ST st = (ST) other;
          return cell.equals(st.cell) && time == st.time;
        } else {
          return false;
        }
      }

      @Override public String toString() {
        // A real implementation would probably allow any cell level, but it's convenient for testing
        // if we restrict to level 30.
        assertEquals(30, cell.level());
        return time + "+" + cell.face() + ":" + Long.toHexString(cell.pos());
      }
    }

    /** A simple spatially positioned feature for testing. */
    static class Feature {
      private final S2CellId cell;
      private final long weight;

      /** Creates a new feature. */
      public Feature(S2CellId cell, long weight) {
        this.cell = cell;
        this.weight = weight;
      }

      /** Returns the S2CellId of this feature. */
      public S2CellId getCell() {
        return cell;
      }

      /** Returns the weight of this feature. */
      public long getWeight() {
        return weight;
      }

      /** Returns true if this feature intersects the given cluster. */
      public boolean intersects(Cluster cluster) {
        return cell.greaterOrEquals(cluster.begin()) && cell.lessThan(cluster.end());
      }

      @Override public int hashCode() {
        return cell.hashCode() ^ (int) weight;
      }

      @Override public boolean equals(Object other) {
        if (other instanceof Feature) {
          Feature f = (Feature) other;
          return cell.equals(f.cell) && weight == f.weight;
        } else {
          return false;
        }
      }
    }
  }

  private void assertRange(S2DensityTree density,
      S2CellId outBegin, S2CellId outEnd,
      S2CellId inBegin, S2CellId inEnd) {
    Cluster cluster = new Cluster(inBegin.rangeMin(), inEnd.rangeMin(), 0).clamp(density);
    assertEquals(0, cluster.weight());
    assertEquals(outBegin, cluster.begin());
    assertEquals(outEnd, cluster.end());
  }

  /** Returns a tree that sums the weights of all entries into all parent cells of those entries. */
  private static S2DensityTree density(Map<S2CellId, Long> bases) {
    TreeEncoder encoder = new TreeEncoder();
    Map<S2CellId, Long> map = new HashMap<>();
    for (Map.Entry<S2CellId, Long> entry : bases.entrySet()) {
      S2CellId cell = entry.getKey();
      long weight = entry.getValue();
      for (int level = 0; level <= cell.level(); level++) {
        S2CellId parent = cell.parent(level);
        map.put(parent, map.getOrDefault(parent, 0L) + weight);
      }
    }
    map.forEach(encoder::put);
    return encoder.build();
  }
}
