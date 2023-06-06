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

import com.google.common.collect.ImmutableMap;
import com.google.common.geometry.S2DensityTree.DecodedPath;
import com.google.common.geometry.S2DensityTree.TreeEncoder;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class S2DensityClusterQueryTest extends GeometryTestCase {

  /**
   * Returns the map that contains the sum of cells in 'bases' and the parents of bases with the
   * corresponding weights.
   */
  private static Map<S2CellId, Long> sumToRoot(Map<S2CellId, Long> bases) {
    return bases.keySet().stream()
        .flatMap(
            cell ->
                IntStream.range(0, cell.level() + 1)
                    .mapToObj(level -> Map.entry(cell.parent(level), bases.get(cell))))
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, Long::sum));
  }

  // Set clustering in the scenario where clusters one or more whole face cells.
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

    // Test low level operations
    DecodedPath path = new DecodedPath(S2DensityClusterQuery.normalizedTree(tree));
    assertEquals(S2CellId.fromFace(0), S2DensityClusterQuery.firstLeafInTree(path));
    assertEquals(S2CellId.fromFace(5), S2DensityClusterQuery.lastLeafInTree(path));
    assertEquals(
        S2CellId.fromFace(1), S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(0)));
    assertEquals(
        S2CellId.fromFace(5), S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(4)));
    assertEquals(
        S2CellId.sentinel(), S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(5)));

    // Cluster the tree
    S2DensityClusterQuery query100 = new S2DensityClusterQuery(100);
    List<S2CellUnion> wt100clusters = query100.clusters(tree);

    // Verify clusters are each a single face.
    assertEquals(6, wt100clusters.size());
    for (S2CellUnion cluster : wt100clusters) {
      assertEquals(1, cluster.cellIds().size());
      S2CellId cell = cluster.cellIds().get(0);
      assertTrue(cell.isFace());
    }

    // If the clusters have double the weight, there should be only three, each covering two faces.
    S2DensityClusterQuery query200 = new S2DensityClusterQuery(200);
    List<S2CellUnion> wt200clusters = query200.clusters(tree);

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

    // Test low level operations
    DecodedPath path = new DecodedPath(S2DensityClusterQuery.normalizedTree(tree));
    assertEquals(S2CellId.fromFace(0).child(1),
        S2DensityClusterQuery.firstLeafInTree(path));
    assertEquals(S2CellId.fromFace(1).child(2).child(0),
        S2DensityClusterQuery.lastLeafInTree(path));

    // There are two leaves under face 0.
    S2CellUnion leavesUnder = new S2CellUnion();
    S2DensityClusterQuery.collectLeavesInTree(path, leavesUnder, S2CellId.fromFace(0));
    assertEquals(2, leavesUnder.size());
    assertEquals(S2CellId.fromFace(0).child(1), leavesUnder.cellIds().get(0));
    assertEquals(S2CellId.fromFace(0).child(2), leavesUnder.cellIds().get(1));

    // And four leaves under face 1.
    leavesUnder.cellIds().clear();
    S2DensityClusterQuery.collectLeavesInTree(path, leavesUnder, S2CellId.fromFace(1));
    assertEquals(4, leavesUnder.size());
    assertEquals(S2CellId.fromFace(1).child(0).child(0), leavesUnder.cellIds().get(0));
    assertEquals(S2CellId.fromFace(1).child(0).child(1), leavesUnder.cellIds().get(1));
    assertEquals(S2CellId.fromFace(1).child(0).child(2), leavesUnder.cellIds().get(2));
    assertEquals(S2CellId.fromFace(1).child(2).child(0), leavesUnder.cellIds().get(3));

    // Test nextLeafInTree for internal nodes.

    // Test nextLeafInTree() for each leaf in the tree.
    assertEquals(
        S2CellId.fromFace(0).child(2),
        S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(0).child(1)));
    assertEquals(
        S2CellId.fromFace(1).child(0).child(0),
        S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(0).child(2)));
    assertEquals(
        S2CellId.fromFace(1).child(0).child(1),
        S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(1).child(0).child(0)));
    assertEquals(
        S2CellId.fromFace(1).child(2).child(0),
        S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(1).child(0).child(2)));
    assertEquals(
        S2CellId.sentinel(),
        S2DensityClusterQuery.nextLeafInTree(path, S2CellId.fromFace(1).child(2).child(0)));

    // If divided into clusters of weight 50, faces 0 and 1 will both be split into two clusters.
    S2DensityClusterQuery query50 = new S2DensityClusterQuery(50);
    List<S2CellUnion> clusters50 = query50.clusters(tree);

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
    List<S2CellUnion> clusters400 = query400.clusters(tree);

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
    List<S2CellUnion> clusters50 = query50.clusters(tree);

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
    List<S2CellUnion> clusters32 = query32.clusters(tree);
    assertEquals(3, clusters32.size());
  }

  public void testInterpolationOnFaces() {
    checkInterpolationOnFace(0);
    checkInterpolationOnFace(1);
    checkInterpolationOnFace(4);
  }

  // Tests accessing tree nodes with DecodedPath.
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
}
