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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Strings;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.geometry.S2DensityClusterQuery.CellInterpolator;
import com.google.common.geometry.S2DensityTree.DecodedPath;
import com.google.common.geometry.S2DensityTree.TreeEncoder;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link S2DensityClusterQuery}. */
@RunWith(JUnit4.class)
public final class S2DensityClusterQueryTest extends GeometryTestCase {

  /**
   * Given a Map of S2CellId to weights for the leaves of a density tree, this function returns an
   * expanded map, adding all the required ancestors of those leaves, assigning them the weight of
   * the sum of the weights of all the leaves below them.
   */
  private static Map<S2CellId, Long> sumToRoot(Map<S2CellId, Long> leafWeights) {
    return leafWeights.keySet().stream()
        .flatMap(
            cell ->
                IntStream.range(0, cell.level() + 1)
                    .mapToObj(level -> Map.entry(cell.parent(level), leafWeights.get(cell))))
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, Long::sum));
  }

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
    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(leaves).forEach(encoder::put);
    S2DensityTree tree = encoder.build();

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

    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(leaves).forEach(encoder::put);
    S2DensityTree tree = encoder.build();

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

    // Note that sumToRoot() is not used here. We are deliberately creating a (valid) S2DensityTree
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
    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(leaves).forEach(encoder::put);
    S2DensityTree tree = encoder.build();

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

    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(leaves).forEach(encoder::put);
    S2DensityTree tree = encoder.build();

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
}
