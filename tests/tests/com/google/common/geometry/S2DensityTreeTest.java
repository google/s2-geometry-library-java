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

import static com.google.common.geometry.S2CellId.FACE_CELLS;
import static com.google.common.geometry.S2CellId.MAX_LEVEL;
import static com.google.common.geometry.S2DensityTree.VERSION;
import static org.junit.Assert.assertThrows;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder.CellWeightFunction;
import com.google.common.geometry.S2DensityTree.Cell;
import com.google.common.geometry.S2DensityTree.CellVisitor.Action;
import com.google.common.geometry.S2DensityTree.DecodedPath;
import com.google.common.geometry.S2DensityTree.IndexCellWeightFunction;
import com.google.common.geometry.S2DensityTree.IntersectionDensity;
import com.google.common.geometry.S2DensityTree.SumDensity;
import com.google.common.geometry.S2DensityTree.TreeEncoder;
import com.google.common.geometry.S2DensityTree.TreeEncoder.ReversedLengthsWriter;
import com.google.common.geometry.S2DensityTree.TreeEncoder.ReversibleBytes;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class S2DensityTreeTest extends GeometryTestCase {
  private final S2Shape point = S2Point.Shape.singleton(
      S2CellId.fromFacePosLevel(3, S2CellId.POS_BITS / 2, 30).toPoint());
  private final S2Shape line = S2LaxPolylineShape.create(
      ImmutableList.of(S2Point.X_POS, S2Point.Y_POS, S2Point.X_NEG));
  private final S2Shape polygon = S2LaxPolygonShape.create(
      ImmutableList.of(ImmutableList.of(S2Point.X_POS, S2Point.Y_POS, S2Point.Z_POS)));
  private final ImmutableMap<S2CellId, Long> weights =
      ImmutableMap.<S2CellId, Long>builder()
          .put(FACE_CELLS[1], 3L)
          .put(FACE_CELLS[1].child(1), 1L)
          .put(FACE_CELLS[1].child(2), 2L)
          .put(FACE_CELLS[2], 4L)
          .put(FACE_CELLS[3], 2L)
          .put(FACE_CELLS[3].child(0), 2L)
          .buildOrThrow();

  private final S2CellId cell22 = S2CellId.fromFace(2).child(2);

  /** A simple S2DensityTree containing two leaves that are children of cell22. */
  private final S2DensityTree cell22Tree = encode(sumToRoot(ImmutableMap.of(
      cell22.child(2), 100L,
      cell22.child(3), 120L)));

  public void testReversibleBytes() {
    ReversibleBytes bytes = new ReversibleBytes();
    bytes.write(VERSION.getBytes(StandardCharsets.ISO_8859_1));
    assertEquals(VERSION, new String(bytes.toByteArray(), StandardCharsets.ISO_8859_1));
    // Reversing twice still gives us the right value.
    bytes.reverseFrom(0);
    assertEquals(VERSION, new String(bytes.reversedCopy(), StandardCharsets.ISO_8859_1));
  }

  public void testReversibleBytesClear() {
    ReversibleBytes bytes = new ReversibleBytes();
    bytes.write(1);
    assertEquals(1, bytes.size());
    bytes.clear();
    assertEquals(0, bytes.size());
  }

  public void testReversedLengthsWriter() throws IOException {
    // Create an output with some data in it before we get started.
    ReversibleBytes out = new ReversibleBytes();
    out.write((byte) 1);
    // Write blocks of length 0, 2, 0, 1.
    ReversedLengthsWriter lengths = new ReversedLengthsWriter(out);
    lengths.next();
    out.write((byte) 2);
    out.write((byte) 3);
    lengths.next();
    lengths.next();
    out.write((byte) 4);
    lengths.next();
    // Write the offsets we have collected, and a mask. The children have no meaning here, and we
    // already know which children are present, so the mask value doesn't matter.
    long mask = 123;
    lengths.write(mask);

    // Read from a reversed copy. We leave off the last block, so we expect lengths 1, 0, 2.
    InputStream in = new ByteArrayInputStream(out.reversedCopy());
    assertEquals("Unexpected mask", mask, EncodedInts.readVarint64(in));
    assertEquals("Unexpected non-empty block length", 1, EncodedInts.readVarint64(in));
    assertEquals("Unexpected non-empty block length", 0, EncodedInts.readVarint64(in));
    assertEquals("Unexpected non-empty block length", 2, EncodedInts.readVarint64(in));
    // Verify bytes of the first block (the last one written).
    assertEquals(4, in.read());
    // Verify bytes of the second block (the first one written).
    assertEquals(3, in.read());
    assertEquals(2, in.read());
    // At this point we're outside the lengths. We wrote a 1 before we started, so make sure that's
    // all we have left.
    assertEquals(1, in.read());
    assertEquals(-1, in.read());
  }

  public void testEncoderWeightTooLarge() {
    assertThrows("Weight out of range", IllegalArgumentException.class,
        () -> new BreadthFirstTreeBuilder(100, 1).build(cell -> Long.MAX_VALUE));
  }

  public void testEncodeEmpty() {
    checkEncoding(ImmutableMap.of());
  }

  public void testEncodeOneFace() {
    checkEncoding(ImmutableMap.of(FACE_CELLS[3], 17L));
  }

  public void testEncodeOneLeaf() {
    checkEncoding(sumToRoot(ImmutableMap.of(S2CellId.fromPoint(S2Point.Y_POS), 123L)));
  }

  public void testEncodeOneBranch() {
    // A pair of level 20 cells with a shared ancestor at level 10.
    S2CellId split = S2CellId.fromFaceIJ(1, 1 << 10, 2 << 10).parent(10);
    checkEncoding(sumToRoot(ImmutableMap.of(split.childBegin(20), 1L, split.childEnd(20), 17L)));
  }

  public void testEncodeEachFace() {
    checkEncoding(IntStream.range(0, 6).boxed()
        .collect(Collectors.toMap(i -> FACE_CELLS[i], i -> 10L + i)));
  }

  public void testVisitWhile() {
    S2DensityTree tree = encode(sumToRoot(ImmutableMap.of(
        S2CellId.fromPoint(new S2Point(1, 2, 3).normalize()), 1L)));
    List<S2CellId> cells = new ArrayList<>();
    assertFalse(tree.visitWhile((cell, weight) -> cell.level() < 4 && cells.add(cell)));
    assertEquals(4, cells.size());
  }

  public void testSelect() {
    // Create a base map of 6*4^4 entries of weight 1.
    Map<S2CellId, Long> base = new TreeMap<>();
    S2CellId end = S2CellId.end(4);
    for (S2CellId cell = S2CellId.begin(4); cell.compareTo(end) < 0; cell = cell.next()) {
      base.put(cell, 1L);
    }
    // Encode the sum of the base nodes into a density tree.
    S2DensityTree tree = encode(sumToRoot(base));
    // Select weight 100 cells. With root weight of 1536, weight 100 will be between the 6 face
    // cells and the 24 level 1 cells, so we'll get the latter.
    List<S2CellId> partitions = tree.getPartitioning(100);
    // Cells are all level 1.
    partitions.forEach(cell -> assertEquals(1, cell.level()));
    // No duplicates and as many level 1 cells as there are.
    assertEquals(24, partitions.stream().distinct().count());
  }

  /** Verifies that random more heavily populated maps work, and that a builder can be reused. */
  public void testEncodeRandomBranches() {
    TreeEncoder encoder = new TreeEncoder();
    for (long weight = 1; weight < 1000; weight++) {
      Map<S2CellId, Long> expected = new HashMap<>();
      for (int j = 0; j < 50; j++) {
        expected.put(data.randomCell(), weight);
      }
      expected = sumToRoot(expected);
      expected.forEach(encoder::put);
      assertEquals(expected, encoder.build().decode());
    }
  }

  /** Verifies an empty shape. */
  public void testShapeIndexEmpty() {
    checkCoverings(ImmutableMap.of(S2Loop.empty(), 1L));
  }

  /** Verifies a point. */
  public void testShapeIndexPoint() {
    checkCoverings(ImmutableMap.of(point, 1L));
  }

  /** Verifies a face center point. */
  public void testShapeIndexFacePoint() {
    checkCoverings(ImmutableMap.of(S2Point.Shape.singleton(FACE_CELLS[2].toPoint()), 1L));
  }

  /** Verifies a line. */
  public void testShapeIndexLine() {
    checkCoverings(ImmutableMap.of(line, 1L));
  }

  /** Verifies a polygon. */
  public void testShapeIndexPolygon() {
    checkCoverings(ImmutableMap.of(polygon, 1L));
  }

  /** Verifies a collection of multiple types. */
  public void testShapeIndexMultiple() {
    checkCoverings(ImmutableMap.of(line, 1L, polygon, 2L));
  }

  /**
   * Verifies the given shape/weight pairs' density measures against the intersects/contains results
   * from S2ShapeIndexRegion, which has the same semantics.
   */
  private void checkCoverings(Map<S2Shape, Long> shapeWeights) {
    S2ShapeIndex index = new S2ShapeIndex();
    shapeWeights.keySet().forEach(index::add);
    S2Region region = new S2ShapeIndexRegion(index);
    S2RegionCoverer coverer = S2RegionCoverer.builder().setMaxCells(64).build();
    S2CellUnion cover = coverer.getCovering(region);
    // Verify cover cells always intersect.
    IndexCellWeightFunction measure = new IndexCellWeightFunction(index, shapeWeights::get);
    for (S2CellId cell : cover) {
      long weight = weight(shapeWeights, cell);
      if (region.contains(new S2Cell(cell))) {
        weight = -weight;
      }
      assertEquals(weight, measure.applyAsLong(cell));
    }
    // Verify the complement does not intersect.
    S2CellUnion complement = new S2CellUnion();
    complement.getDifference(coverer.getCovering(S2Loop.full()), cover);
    for (S2CellId cell : complement) {
      // Unless the cell is on the border between 'cover' and its complement, since both
      // S2ShapeIndexRegion and S2DensityTree use a small amount of outward padding.
      long expected = cover.intersects(cell) ? weight(shapeWeights, cell) : 0;
      assertEquals(expected, measure.applyAsLong(cell));
    }
  }

  /**
   * Returns the sum of weights that intersect 'cell' using separate S2ShapeIndexRegion instances,
   * which must be done to test against that particular mayIntersect class's robustness properties.
   */
  private long weight(Map<S2Shape, Long> shapeWeights, S2CellId cell) {
    long sum = 0;
    for (S2Shape shape : shapeWeights.keySet()) {
      S2ShapeIndex index = new S2ShapeIndex();
      index.add(shape);
      if (new S2ShapeIndexRegion(index).mayIntersect(new S2Cell(cell))) {
        sum += shapeWeights.get(shape);
      }
    }
    return sum;
  }

  public void testFeatureDensity() {
    S2ShapeIndex shapes = new S2ShapeIndex();
    IdentityHashMap<S2Shape, String> features = new IdentityHashMap<>();
    Map<String, Integer> weights = new HashMap<>();
    S2Point p = S2LatLng.fromDegrees(5, 5).toPoint();
    S2Point q = S2LatLng.fromDegrees(-5, 5).toPoint();
    // Define alice.
    weights.put("Alice", 1);
    S2Shape alice1 = S2Point.Shape.singleton(p);
    shapes.add(alice1);
    features.put(alice1, "Alice");
    S2Shape alice2 = S2Point.Shape.singleton(q);
    shapes.add(alice2);
    features.put(alice2, "Alice");
    // Define bob.
    weights.put("Bob", 5);
    S2Shape bob1 = S2Point.Shape.singleton(p);
    shapes.add(bob1);
    features.put(bob1, "Bob");
    int maxBytes = 100;
    int maxLevel = 1;
    S2DensityTree tree = S2DensityTree.featureDensity(
        maxBytes, maxLevel, shapes, features::get, weights::get);
    assertEquals(
        ImmutableMap.builder()
            .put(S2CellId.fromPoint(p).parent(0), 6L)
            .put(S2CellId.fromPoint(p).parent(1), 6L)
            .put(S2CellId.fromPoint(q).parent(1), 1L)
            .buildOrThrow(),
        tree.decode());
  }

  /** Tests accessing tree leaves with DecodedPath. */
  public void testDecodedPathAccessLeaves() {
    S2DensityTree.DecodedPath decoder = new S2DensityTree.DecodedPath(cell22Tree);

    // Check that all nodes can be obtained by DecodedPath.
    for (Map.Entry<S2CellId, Long> leaf : cell22Tree.decode().entrySet()) {
      S2DensityTree.Cell cell = decoder.cell(leaf.getKey());
      assertEquals(cell.weight, (long) leaf.getValue());
    }
  }

  /** Verifies that DecodedPath returns empty nodes when cells are not present in the tree. */
  public void testDecodedPathNodesNotPresentAreEmpty() {
    S2DensityTree.DecodedPath decoder = new S2DensityTree.DecodedPath(cell22Tree);

    // Face nodes not in the tree should return an empty cell.
    assertTrue(decoder.cell(S2CellId.fromFace(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(1)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(5)).isEmpty());

    // Level 1 nodes of a face that isn't in the tree
    assertTrue(decoder.cell(S2CellId.fromFace(1).child(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(1).child(1)).isEmpty());

    // Level 1 nodes of the face that is in the tree
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(1)).isEmpty());

    // A child of a level 1 node that is not present, is also not present.
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(3).child(0)).isEmpty());

    // cell22 and its children in positions 2 and 3 are in the tree, all others are not.
    assertTrue(decoder.cell(cell22.child(0)).isEmpty());
    assertTrue(decoder.cell(cell22.child(1)).isEmpty());

    assertTrue(decoder.cell(cell22.child(1).child(0)).isEmpty());
    assertTrue(decoder.cell(cell22.child(1).child(1)).isEmpty());
  }

  /** Verifies that DecodedPath returns the expected weights for nodes that are in the tree. */
  public void testDecodedPathNodesHaveExpectedWeights() {
    S2DensityTree.DecodedPath decoder = new S2DensityTree.DecodedPath(cell22Tree);

    // Face 2 should be present with the sum of the leaf weights.
    assertEquals(220, decoder.cell(S2CellId.fromFace(2)).weight);
    assertEquals(220, decoder.cell(S2CellId.fromFace(2).child(2)).weight);

    assertEquals(100, decoder.cell(cell22.child(2)).weight);
    assertEquals(120, decoder.cell(cell22.child(3)).weight);
  }

  /** Verifies that we get no density from no inputs. */
  public void testEmptyInputs() {
    checkSum(ImmutableMap.of(), weights);
    checkIntersection(ImmutableMap.of(), weights);
  }

  /**
   * Verifies that we get the same tree back for the sum of one tree, or intersection of a tree with
   * itself (but with twice the weight).
   */
  public void testCombineOne() {
    checkSum(
        ImmutableMap.of(FACE_CELLS[1], 3L,
            FACE_CELLS[1].child(1), 1L,
            FACE_CELLS[1].child(2), 2L),
        weights,
        FACE_CELLS[1]);

    checkIntersection(
        ImmutableMap.of(
            FACE_CELLS[1], 6L,
            FACE_CELLS[1].child(1), 2L,
            FACE_CELLS[1].child(2), 4L),
        weights,
        FACE_CELLS[1],
        FACE_CELLS[1]);
  }

  /** Verifies that we get proper results when one tree contains the other. */
  public void testSumNested() {
    checkSum(
        ImmutableMap.of(FACE_CELLS[1], 4L, // 3 + 1
            FACE_CELLS[1].child(1), 2L,    // 1 + 1
            FACE_CELLS[1].child(2), 2L),   // 2
        weights,
        FACE_CELLS[1], FACE_CELLS[1].child(1));

    checkIntersection(
        ImmutableMap.of(
            FACE_CELLS[1], 4L,            // 3 + 1
            FACE_CELLS[1].child(1), 2L),  // 1 + 1
            // face 1, child 2 isn't in the intersection.
        weights,
        FACE_CELLS[1], FACE_CELLS[1].child(1));
  }

  /** Tests intersection of leaves of two trees. */
  public void testIntersectionOfLeavesTwoWay() {
    ImmutableMap<S2CellId, Long> weightsLeft =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 10L)
            .put(FACE_CELLS[1].child(1), 6L)
            .put(FACE_CELLS[1].child(2), 4L)
            .put(FACE_CELLS[1].child(3), 4L)
        .buildOrThrow();

    ImmutableMap<S2CellId, Long> weightsRight =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 10L)
            .put(FACE_CELLS[1].child(2), 6L)
            .put(FACE_CELLS[1].child(3), 4L)
            .put(FACE_CELLS[1].child(4), 4L)
        .buildOrThrow();

    List<S2DensityTree> left = buildTrees(weightsLeft, FACE_CELLS[1]);
    List<S2DensityTree> right = buildTrees(weightsRight, FACE_CELLS[1]);
    ImmutableList<S2DensityTree> trees = ImmutableList.of(left.get(0), right.get(0));
    S2DensityTree.IntersectionDensity intersectionDensity
        = S2DensityTree.IntersectionDensity.create(trees);
    S2CellUnion intersection = intersectionDensity.getIntersectingLeaves();

    assertTrue(intersection.cellIds().contains(FACE_CELLS[1].child(2)));
    assertTrue(intersection.cellIds().contains(FACE_CELLS[1].child(3)));
  }

  /** Verifies that combining trees caps cell weights at MAX_WEIGHT. */
  public void testCombinationExceedingMaxWeight() {
    // A map of cells to weights that are each less than the per-cell MAX_WEIGHT, but large enough
    // that adding together the weights will exceed the max.
    long max = Long.MAX_VALUE >>> 4;
    ImmutableMap<S2CellId, Long> heavyWeights =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], max / 2)
            .put(FACE_CELLS[1].child(1), max / 2)
            .put(FACE_CELLS[1].child(2), max / 2)
            .buildOrThrow();

    checkSum(
        ImmutableMap.of(FACE_CELLS[1], max,
            FACE_CELLS[1].child(1), max,
            FACE_CELLS[1].child(2), max / 2), // Child 2 is only in one of the trees.
        heavyWeights,
        FACE_CELLS[1], FACE_CELLS[1].child(1), FACE_CELLS[1].child(1));

    checkIntersection(
        ImmutableMap.of(
            FACE_CELLS[1], max,
            FACE_CELLS[1].child(1), max),
            // face 1, child 2 isn't in the intersection.
        heavyWeights,
        FACE_CELLS[1], FACE_CELLS[1].child(1), FACE_CELLS[1].child(1));
  }

  /** Verifies two disjoint trees don't affect each other. */
  public void testCombineDisjoint() {
    checkSum(
        ImmutableMap.of(
            FACE_CELLS[2], 4L,
            FACE_CELLS[3], 2L, FACE_CELLS[3].child(0), 2L),
        weights,
        FACE_CELLS[2], FACE_CELLS[3]);

    checkIntersection(
        ImmutableMap.of(),
        weights,
        FACE_CELLS[2], FACE_CELLS[3]);
  }

  /** Verifies the intersection of two trees is empty if they only have internal nodes in common. */
  public void testIntersection() {
    // Despite both having weight in face cell 1, the intersection of the trees is empty.
    checkIntersection(
        ImmutableMap.of(),
        weights,
        FACE_CELLS[1].child(1), FACE_CELLS[1].child(2));
  }

  public void testSumMaxLevel() {
    S2CellId cell = S2CellId.FACE_CELLS[5].child(2).child(1).child(0);
    for (int maxLevel = 0; maxLevel <= cell.level(); maxLevel++) {
      assertEquals(
          sumToRoot(ImmutableMap.of(cell.parent(maxLevel), 1L)),
          new BreadthFirstTreeBuilder(1_000, maxLevel, new TreeEncoder())
              .build(x -> x.intersects(cell) ? 1L : 0L).decode());
    }
  }

  public void testDecoderWithMissingUncle() {
    TreeEncoder encoder = new TreeEncoder();
    // Create two trees with 4000 weight total in the same level 2 cell. The last child at level 1
    // (the uncle) doesn't intersect, so level 2 processing must reload their actual parent.
    S2CellId leaf = S2CellId.fromToken("55");
    sumToRoot(ImmutableMap.of(leaf, 2052L)).forEach(encoder::put);
    S2DensityTree a = encoder.build();
    sumToRoot(ImmutableMap.of(leaf, 1948L)).forEach(encoder::put);
    S2DensityTree b = encoder.build();
    S2DensityTree sum = S2DensityTree.sumDensity(ImmutableList.of(a, b), 1 << 20, 18);
    assertTrue(a.decode().containsKey(leaf));
    assertTrue(b.decode().containsKey(leaf));
    assertTrue(sum.decode().containsKey(leaf));
  }

  /**
   * Creates a density tree for each root using the given map of weights, and returns those trees.
   */
  private static List<S2DensityTree> buildTrees(Map<S2CellId, Long> weights, S2CellId ... roots) {
    TreeEncoder encoder = new TreeEncoder();
    List<S2DensityTree> trees = new ArrayList<>();
    for (S2CellId root : roots) {
      // Insert weights in 'weights' under the given roots.
      insert(weights, encoder, root);
      // And then if this isn't a leaf cell, insert the root weight into all cells above it.
      long weight = weights.get(root);
      while (root.level() > 0) {
        root = root.parent();
        encoder.put(root, weight);
      }
      trees.add(encoder.build());
    }
    return trees;
  }

  /**
   * Creates a density tree for each root in the given 'roots' using the given 'weights', sums the
   * trees using SumDensity, and verifies that the sum matches the given 'expected' weights.
   */
  private static void checkSum(
      Map<S2CellId, Long> expected, Map<S2CellId, Long> weights, S2CellId ... roots) {
    List<S2DensityTree> trees = buildTrees(weights, roots);
    SumDensity sum = new SumDensity();
    sum.set(trees);
    TreeEncoder encoder = new TreeEncoder();
    S2DensityTree sumTree = new BreadthFirstTreeBuilder(1_000, 20, encoder).build(sum);
    assertEquals(expected, sumTree.decode());
  }

  /**
   * Creates a density tree for each root in the given 'roots' using the given 'weights', aggregates
   * the trees using IntersectionDensity, and verifies that the combined tree matches the given
   * 'expected' weights.
   */
  private static void checkIntersection(
      Map<S2CellId, Long> expected, Map<S2CellId, Long> weights, S2CellId... roots) {
    List<S2DensityTree> trees = buildTrees(weights, roots);
    IntersectionDensity intersection = IntersectionDensity.create(trees);
    intersection.set(trees);
    TreeEncoder encoder = new TreeEncoder();
    Map<S2CellId, Long> actual =
        new BreadthFirstTreeBuilder(1_000, 20, encoder).build(intersection).decode();

    StringBuilder message = new StringBuilder();
    message.append("Input Trees: \n");
    for (int i = 0; i < roots.length; i++) {
      message
          .append("  tree #")
          .append(i)
          .append(" to string = ")
          .append(toString(trees.get(i)));
      message
          .append("  tree #")
          .append(i)
          .append(" leaves = ")
          .append(toTokenSet(trees.get(i).getLeaves()))
          .append("\n\n");
    }
    message
        .append("\nIntersection of Input Tree Leaves: \n")
        .append(toTokenSet(intersection.getIntersectingLeaves()))
        .append("\n");
    message.append("\nExpected Weights: \n").append(decodedTreeToString(expected)).append("\n");
    message.append("\nActual Weights: \n").append(decodedTreeToString(actual)).append("\n");
    assertEquals(message.toString(), expected, actual);
  }

  private static void insert(Map<S2CellId, Long> weights, TreeEncoder encoder, S2CellId cell) {
    if (weights.containsKey(cell)) {
      encoder.put(cell, weights.get(cell));
      for (int i = 0; i < 4; i++) {
        insert(weights, encoder, cell.child(i));
      }
    }
  }

  /**
   * Dilate a tree with level 2 leaves, with dilation constrained to level 2, which are much
   * larger cells than necessary for the dilation radius.
   */
  public void testSmallDilationConstrainedToLeafLevel() {
    // A density tree with two leaves at level 2.
    final ImmutableMap<S2CellId, Long> twoLevelTwo =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 4L)
            .put(FACE_CELLS[1].child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(3), 2L)
            .put(FACE_CELLS[1].child(3).child(3), 2L)
            .buildOrThrow();

    S2DensityTree tree = encode(twoLevelTwo);

    // Dilate the tree by 1 km, which is very small relative to the cell size. Use maxLevelDiff 0,
    // so the dilated tree will have weight in the level 2 cells neighboring each of the two level 2
    // leaves. They are in face corners and have 7 neighbors each, so the resulting dilated tree has
    // 14 added level 2 leaf cells as well as the original two.
    S2DensityTree dilatedTree = S2DensityTree.dilate(tree, metersToAngle(1000), 0);
    ImmutableSet<String> expectedNodes = ImmutableSet.of(
        "1",
        /**/ "14",
        /**/ /**/ "15",
        /**/ /**/ "17",
        "3",
        /**/ "2c",
        /**/ /**/ "29",
        /**/ /**/ "2b", // Input node
        /**/ /**/ "2d",
        /**/ /**/ "2f",
        /**/ "3c",
        /**/ /**/ "39",
        /**/ /**/ "3b",
        /**/ /**/ "3d",
        /**/ /**/ "3f", // Input node
        "5",
        /**/ "44",
        /**/ /**/ "41",
        /**/ /**/ "43",
        "7",
        /**/ "6c",
        /**/ /**/ "69",
        /**/ /**/ "6b",
        "b",
        /**/ "ac",
        /**/ /**/ "ab",
        /**/ /**/ "ad"
    );
    Set<String> actualNodes = toTokenSet(dilatedTree);
    assertEquals(expectedNodes, actualNodes);
  }

  /**
   * Dilate a tree with level 2 leaves with dilation constrained to level 3, which are still much
   * larger cells than necessary for the dilation radius.
   */
  public void testSmallDilationRelativeToLeafSize() {
    // A density tree with two leaves at level 2.
    final ImmutableMap<S2CellId, Long> twoLevelTwo =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 4L)
            .put(FACE_CELLS[1].child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(3), 2L)
            .put(FACE_CELLS[1].child(3).child(3), 2L)
            .buildOrThrow();

    S2DensityTree tree = encode(twoLevelTwo);

    // Dilate the tree by 1 km, which is very small relative to the cell size. Use maxLevelDiff 1,
    // so the dilated tree will have weight in the two level 3 cells neighboring each of the four
    // sides of the two level 2 leaves and three (these are at the cube corners) level 3 diagonal
    // corners. The dilated tree ends up with 11 additional level 3 cells around each of the two
    // original level 2 leaf cells.
    S2DensityTree dilatedTree = S2DensityTree.dilate(tree, metersToAngle(1000), 1);
    assertEquals(24, dilatedTree.getLeaves().size());
  }

  /**
   * Dilate two versions of a tree with level 2 leaves that share a neighbor, and make sure that
   * the maximum weight is used in dilated nodes, regardless of what order they are added during
   * dilation.
   */
  public void testDilationUsesMaximum() {
    // Two density trees with two leaves at level 2 with a common neighbor "3b" not in the tree.
    // The only difference is the distribution of weight.
    final ImmutableMap<S2CellId, Long> tree1 =
        ImmutableMap.<S2CellId, Long>builder()
            .put(S2CellId.fromToken("3"), 10L)     // Face
            .put(S2CellId.fromToken("3c"), 2L)     // Level 1
            .put(S2CellId.fromToken("3d"), 2L)     // Level 2 child of 3c
            .put(S2CellId.fromToken("34"), 8L)     // Level 1
            .put(S2CellId.fromToken("31"), 8L)     // Level 2 child of 34
            .buildOrThrow();

    final ImmutableMap<S2CellId, Long> tree2 =
        ImmutableMap.<S2CellId, Long>builder()
            .put(S2CellId.fromToken("3"), 10L)     // Face
            .put(S2CellId.fromToken("3c"), 8L)     // Level 1
            .put(S2CellId.fromToken("3d"), 8L)     // Level 2 child of 3c
            .put(S2CellId.fromToken("34"), 2L)     // Level 1
            .put(S2CellId.fromToken("31"), 2L)     // Level 2 child of 34
            .buildOrThrow();

    // Dilate both with maxlevelDiff 0, so cell "3b" is added.
    S2DensityTree dilatedTree1 = S2DensityTree.dilate(encode(tree1), metersToAngle(1000), 0);
    S2DensityTree dilatedTree2 = S2DensityTree.dilate(encode(tree2), metersToAngle(1000), 0);

    // Check that the common neighbor "3b" gets the maximum weight 8 in both trees.
    dilatedTree1.visitCells((cellId, node) -> {
      if (cellId.equals(S2CellId.fromToken("3b"))) {
        assertEquals(8, node.weight());
      }
      return S2DensityTree.CellVisitor.Action.ENTER_CELL;
    });
    dilatedTree2.visitCells((cellId, node) -> {
      if (cellId.equals(S2CellId.fromToken("3b"))) {
        assertEquals(8, node.weight());
      }
      return S2DensityTree.CellVisitor.Action.ENTER_CELL;
    });
  }

  /**
   * Dilate a tree with level 5 leaves, which average about 288 km across, by 1000 km, which
   * requires using dilation with cells at level 2.
   */
  public void testDilationLargerThanLeafSize() {
    // A density tree with two leaves at level 5.
    final ImmutableMap<S2CellId, Long> twoLevelFive =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 4L)
            .put(FACE_CELLS[1].child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1).child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1).child(1).child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(3), 2L)
            .put(FACE_CELLS[1].child(1).child(3).child(3), 2L)
            .put(FACE_CELLS[1].child(1).child(3).child(3).child(3), 2L)
            .put(FACE_CELLS[1].child(1).child(3).child(3).child(3).child(3), 2L)
            .buildOrThrow();

    S2DensityTree tree = encode(twoLevelFive);

    // maxLevelDiff is set to 4, but the dilation level is limited by the dilation radius and will
    // actually be 2. The result is that the two level 2 nodes in the tree have all their level 2
    // neighbors added, and the nodes at higher levels are dropped. It is helpful to visualize the
    // results in supermario/s2viewer.
    S2DensityTree dilatedTree =
        S2DensityTree.dilate(tree, metersToAngle(1000 * 1000), 4);

    ImmutableSet<String> expectedNodes = ImmutableSet.of(
        "3",
        /**/ "24",
        /**/ /**/ "25",
        /**/ /**/ "27",
        /**/ "2c",
        /**/ /**/ "29",
        /**/ /**/ "2b",
        /**/ /**/ "2d",
        /**/ /**/ "2f",
        /**/ "34",
        /**/ /**/ "31",
        /**/ /**/ "33",
        /**/ "3c",
        /**/ /**/ "3b",
        "7",
        /**/ "6c",
        /**/ /**/ "69",
        /**/ /**/ "6b",
        "b",
        /**/ "ac",
        /**/ /**/ "ab",
        /**/ /**/ "ad"
    );
    Set<String> actualNodes = toTokenSet(dilatedTree);
    assertEquals(expectedNodes, actualNodes);
  }

  /**
   * Checks sumDensity and sumIntersectionDensity with dilated trees, with dilation constrained to
   * at most the leaf level of the input trees (maxLevelDiff == 0).
   */
  public void testDilateAtLeafLevelAndCombineDensity() {
    // A density tree with one leaf at level 2.
    final ImmutableMap<S2CellId, Long> oneLevelTwo =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 2L) // 3
            .put(FACE_CELLS[1].child(2), 2L) // 34
            .put(FACE_CELLS[1].child(2).child(1), 2L) // 33
            .buildOrThrow();

    // A density tree with two leaves at level 2.
    final ImmutableMap<S2CellId, Long> twoLevelTwo =
        ImmutableMap.<S2CellId, Long>builder()
            .put(FACE_CELLS[1], 4L)
            .put(FACE_CELLS[1].child(1), 2L)
            .put(FACE_CELLS[1].child(1).child(1), 2L)
            .put(FACE_CELLS[1].child(3), 2L)
            .put(FACE_CELLS[1].child(3).child(3), 2L)
            .buildOrThrow();

    S2DensityTree tree1 = encode(oneLevelTwo);
    S2DensityTree tree2 = encode(twoLevelTwo);

    assertEquals(toTokenSet(tree1.getLeaves()), ImmutableSet.of("33"));
    assertEquals(toTokenSet(tree2.getLeaves()), ImmutableSet.of("3f", "2b"));

    // Dilate both trees.
    ImmutableList<S2DensityTree> trees = ImmutableList.of(
        S2DensityTree.dilate(tree1, metersToAngle(1000), 0),
        S2DensityTree.dilate(tree2, metersToAngle(1000), 0));

    String message = String.format(
        "Pre-Dilation Tree1 Leaves: \n%s\n"
        + "\nDilated Tree1 Leaves: \n%s\n"
        + "\nPre-Dilation Tree2 Leaves: \n%s\n"
        + "\nDilated Tree2 Leaves: \n%s\n",
        toTokenSet(tree1.getLeaves()),
        toTokenSet(trees.get(0).getLeaves()),
        toTokenSet(tree2.getLeaves()),
        toTokenSet(trees.get(1).getLeaves()));

    // Dilated tree 1 has leaves that are the nine level-2 cells of "33" (roughly the Philippenes)
    // and its 8 neighbors. Every node has weight 2.
    assertEquals(
        message,
        toTokenSet(trees.get(0).getLeaves()),
        ImmutableSet.of("61", "2d", "35", "69", "33", "2f", "37", "67", "31"));
    trees.get(0).visitCells((cellId, node) -> {
      assertEquals("Node " + cellId.toToken(), 2, node.weight());
      return S2DensityTree.CellVisitor.Action.ENTER_CELL;
    });

    // Dilated tree 2 has eight level-2 cells for "3f" and "2b" (roughly Iran and western Australia)
    // and their *seven* neighbors each (they are at the corners of their face). They don't overlap.
    // All the nodes have weight 2, except for the face with weight 4,
    assertEquals(
        message,
        toTokenSet(trees.get(1).getLeaves()),
        ImmutableSet.of(
            "3b", "69", "ad", "2b", "41", "15", "3d", "2d",
            "43", "6b", "17", "39", "3f", "29", "2f", "ab"));
    trees.get(1).visitCells((cellId, node) -> {
      int expectedWeight = cellId.toToken().equals("3") ? 4 : 2;
      assertEquals("Node " + cellId.toToken(), expectedWeight, node.weight());
      return S2DensityTree.CellVisitor.Action.ENTER_CELL;
    });

    // Combine the dilated trees with sumDensity. The limit on the size of the combined tree is
    // large enough to not be a factor.
    S2DensityTree dilatedSum = S2DensityTree.sumDensity(trees, 16384, 2);

    // The three dilated level-2 cells south of "33" (the Philippenes) are 2f, 2d, and 69. These
    // are also the three dilated cells north of "2b" (Western Australia). So the summed tree has a
    // total of 16 + 9 - 3 = 22 leaves, where those three have weight 4 and all the others have
    // weight 2.
    // Face 0
    Set<String> actualLeaves = toTokenSet(dilatedSum.getLeaves());
    ImmutableSet<String> expectedUnionLeaves = ImmutableSet.of(
            "2d", "2f", "69", // The three neighbors shared between "33" and "2b"
            "33", "31", "61", "35", "37", "67", // "33" and its non-overlapping 5 neighbors
            "29", "6b", "2b", "ad", "ab", // "2b" and its non-overlapping 4 neighbors
            "41", "43",  "15", "3f", "39",  "17", "3d", "3b" // "3f" and its 7 neighbors.
            );

    assertEquals(
        message + "\nDilated Sum Leaves: " + actualLeaves, expectedUnionLeaves, actualLeaves);

    // Cells not listed here have weight 2.
    ImmutableMap<S2CellId, Long> expectedWeights =
      ImmutableMap.<S2CellId, Long>builder()
          // The three overlapping leaves each have summed weight 4.
          .put(S2CellId.fromToken("2d"), 4L)
          .put(S2CellId.fromToken("2f"), 4L)
          .put(S2CellId.fromToken("69"), 4L)
          // Face 3 has weight 6, as the sum of the weights of the input tree faces.
          .put(S2CellId.fromToken("3"), 6L)
          // Level 1 cells "2c" and "6c" have weight 4, propagated up from their child leaves.
          .put(S2CellId.fromToken("2c"), 4L)
          .put(S2CellId.fromToken("6c"), 4L)
          // Face 7 has weight 4, propagated up from child "69".
          .put(S2CellId.fromToken("7"), 4L)
          .buildOrThrow();

    dilatedSum.visitCells(
        (cellId, node) -> {
          if (expectedWeights.containsKey(cellId)) {
            assertEquals(cellId.toToken(), (long) expectedWeights.get(cellId), node.weight());
          } else {
            assertEquals("Other Leaf " + cellId.toToken(), 2, node.weight());
          }
          return S2DensityTree.CellVisitor.Action.ENTER_CELL;
        });

    // Intersect the dilated trees. Only the three overlapping cells are expected as leaves, with
    // weight 4.
    S2DensityTree dilatedIntersection = S2DensityTree.intersectionDensity(trees, 16384, 2);
    actualLeaves = toTokenSet(dilatedIntersection.getLeaves());
    ImmutableSet<String> expectedIntersectionLeaves = ImmutableSet.of("2d", "2f", "69");

    assertEquals(message + "\nDilated Intersection Leaves: " + actualLeaves,
        actualLeaves, expectedIntersectionLeaves);
    dilatedIntersection.visitCells((cellId, node) -> {
      assertEquals(
            decodedTreeToString(dilatedSum.decode()) + "\nNon-leaf " + cellId.toToken(),
            (long) expectedWeights.get(cellId),
            node.weight());
      return S2DensityTree.CellVisitor.Action.ENTER_CELL;
    });
  }

  /** A test case with geometry from HD roads. */
  public void testIntersectionFromIndexes() {
    final int maxLevel = 20;

    // The left and right shapes are about 300 meters apart.
    S2ShapeIndex leftShapes =
        S2TextFormat.makeIndexOrDie(
            "# 0:0.0089989, 0:0.0179977 | 0:0.0179977, 0:0.0269966 "
                + "# 3.15184885545761e-05:0.0179977305305156, "
                + "3.15184897208095e-05:0.00899890901688273, "
                + "-3.15184897208095e-05:0.00899890901688273, "
                + "-3.15184885545761e-05:0.0179977305305156 "
                + "| 3.15184866108635e-05:0.0269965337367333, "
                + "3.15184885545761e-05:0.0179977305305156, "
                + "-3.15184885545761e-05:0.0179977305305156, "
                + "-3.15184866108635e-05:0.0269965337367333");
    S2DensityTree leftTree = S2DensityTree.vertexDensity(leftShapes, 2048, maxLevel);

    S2ShapeIndex rightShapes =
        S2TextFormat.makeIndexOrDie(
            "# 0.00224971648787036:0.00449943297574072, 0.00224971648787036:0.0134982989272222 #");
    S2DensityTree rightTree = S2DensityTree.vertexDensity(rightShapes, 2048, maxLevel);

    // Dilate both trees by at least half the distance between them.
    S1Angle dilationAngle = metersToAngle(150);
    ImmutableList<S2DensityTree> dilatedTrees =
        ImmutableList.of(
            S2DensityTree.dilate(leftTree, dilationAngle, 0),
            S2DensityTree.dilate(rightTree, dilationAngle, 0));

    String message =
        String.format(
            "Pre-Dilation LeftTree Leaves: \n%s\n"
                + "\nDilated LeftTree Leaves: \n%s\n"
                + "\nPre-Dilation RightTree Leaves: \n%s\n"
                + "\nDilated RightTree Leaves: \n%s\n",
            toTokenSet(leftTree.getLeaves()),
            toTokenSet(dilatedTrees.get(0).getLeaves()),
            toTokenSet(rightTree.getLeaves()),
            toTokenSet(dilatedTrees.get(1).getLeaves()));

    // The intersection of the dilated trees should not be empty.
    S2DensityTree dilatedIntersection = S2DensityTree.intersectionDensity(
        dilatedTrees, 16384, maxLevel);
    assertFalse(message, dilatedIntersection.getLeaves().isEmpty());
  }

  public void testVisitCellsAtFaceCenter() {
    // Two Level 16 cells near the center of face 0, in different level 1 cells.
    String[] tokens = { "0ffffffd5", "10000002b" };

    // Turn the tokens into an S2DensityTree, using them as leaves with weight 1.
    IdentityHashMap<S2CellId, Long> cellWeights = new IdentityHashMap<>();
    for (String token : tokens) {
      cellWeights.put(S2CellId.fromToken(token), 1L);
    }
    S2DensityTree tree = encode(sumToRoot(cellWeights));

    // Get the set of cells in the tree.
    Set<S2CellId> treeCellIds = tree.decode().keySet();

    // Validate that all the ancestors of the leaves are in the tree, and gather the set of cells
    // we expect to visit.
    Set<S2CellId> expected = new HashSet<>();
    for (String token : tokens) {
      S2CellId cellId = S2CellId.fromToken(token);
      while (cellId.level() >= 0) {
        if (cellId.level() <= 14) {
          expected.add(cellId);
        }
        assertTrue(treeCellIds.contains(cellId));
        cellId = cellId.parent();
      }
    }

    // Visit cells with children down to level 14, keeping track of what was actually visited.
    Set<S2CellId> actual = new HashSet<>();
    tree.visitCells(
        (cellId, node) -> {
          actual.add(cellId);
          if (node.hasChildren() && cellId.level() < 14) {
            return S2DensityTree.CellVisitor.Action.ENTER_CELL;
          }
          return S2DensityTree.CellVisitor.Action.SKIP_CELL;
        });

    assertEquals(toTokenSet(expected), toTokenSet(actual));
  }

  public void testDilationAtFaceCenter() {
    // Two Level 16 cells near the center of face 0, in different level 1 cells.
    String[] tokens = { "0ffffffd5", "10000002b" };

    // Turn the tokens into an S2DensityTree, using them as leaves with weight 1.
    IdentityHashMap<S2CellId, Long> cellWeights = new IdentityHashMap<>();
    for (String token : tokens) {
      cellWeights.put(S2CellId.fromToken(token), 1L);
    }
    S2DensityTree tree = encode(sumToRoot(cellWeights));

    // Dilate the tree by 300 meters, which requires using cells at level 14.
    S2DensityTree dilated = S2DensityTree.dilate(tree, metersToAngle(300), 0);

    // The level 14 parents of the two leaves are 0ffffffd and 10000003. They are adjacent, and
    // surrounded by 10 more level 14 cells, forming a 4 by 3 grid of cells which should be the
    // leaves of the dilated tree.
    String[] expectedDilatedLeaves = {
      "0fffffe5", "0fffffe3", "1000001d", "1000001b",
      "0ffffffb", "0ffffffd", "10000003", "10000005",
      "0ffffff9", "0fffffff", "10000001", "10000007"
    };
    Set<String> expected = new HashSet<>(Arrays.asList(expectedDilatedLeaves));

    String message =
        String.format(
            "Pre-dilated tree = %s \nwith leaves: %s\nDilated tree = %s \nwith leaves: %s\n",
            toTokenSet(tree),
            toTokenSet(tree.getLeaves()),
            toTokenSet(dilated),
            toTokenSet(dilated.getLeaves()));

    Set<String> actual = toTokenSet(dilated.getLeaves());
    assertEquals(message, expected, actual);
  }

  /**
   * Verifies cells with weight intersect, and unsplit nodes with weight are contained by the
   * S2Region for a tree.
   */
  public void testDensityRegion() {
    S2DensityTree tree = encode(weights);
    S2Region region = tree.asRegion();
    tree.visitCells((S2CellId cell, Cell node) -> {
      assertEquals(node.weight() > 0, region.mayIntersect(new S2Cell(cell)));
      assertEquals(node.weight() > 0 && !node.hasChildren(), region.contains(new S2Cell(cell)));
      return Action.ENTER_CELL;
    });
  }

  /** Verifies a large random set of cells. */
  public void testDensitySize() {
    // Build a large random weights map and a function that measures intersecting weights from it.
    // We may get slightly fewer than the desired number of cells since duplicates overwrite.
    TreeMap<S2CellId, Long> weights = new TreeMap<>();
    for (int i = 0; i < 1_000; i++) {
      S2CellId parent = data.randomCell();
      for (S2CellId child : parent.childrenAtLevel(Math.min(MAX_LEVEL, parent.level() + 5))) {
        weights.put(child, i + 500L);
      }
    }
    CellWeightFunction weigher = target -> sum(weights, target);

    // Sum up the requested vs. actual sizes for a variety of requested sizes, not limiting by cell
    // level at all.
    long requestedSum = 0;
    long actualSum = 0;
    for (int requestedSize = 1; requestedSize < 10_000_000; requestedSize *= 2) {
      BreadthFirstTreeBuilder builder = new BreadthFirstTreeBuilder(requestedSize, MAX_LEVEL);
      S2DensityTree result = builder.build(weigher);
      long actualSize = result.size();
      requestedSum += requestedSize;
      actualSum += actualSize;
    }

    // The ratio can be as little as 0.5 or as great as 4, in general, although very small maps may
    // have much worse ratios (due to the header and face encodings) and this sample we've measured
    // has been averaged quite a bit, and isn't as spread out as real trees are, so the actual value
    // computed here is about 0.7.
    double ratio = 1.0 * actualSum / requestedSum;
    assertTrue("Ratio was %s" + ratio, ratio >= 0.5 && ratio <= 4);
  }

  /** Verifies that the sum operator treats leaves and only leaves as contained. */
  public void testSumAtLeaf() {
    S2DensityTree tree = encode(weights);
    SumDensity sum = new SumDensity();
    sum.set(ImmutableList.of(tree));
    DecodedPath path = new DecodedPath(tree);
    weights.forEach((cell, weight) -> {
      Cell node = path.cell(cell);
      assertEquals(node.hasChildren(), sum.applyAsLong(cell) > 0);
    });
  }

  /** Returns the result of encoding the given weights into a tree. */
  private static S2DensityTree encode(Map<S2CellId, Long> weights) {
    TreeEncoder encoder = new TreeEncoder();
    weights.forEach(encoder::put);
    return encoder.build();
  }

  /** Returns the sum of weights that intersect 'query'. */
  private static long sum(SortedMap<S2CellId, Long> weights, S2CellId query) {
    long[] sum = {0};
    weights.subMap(query.rangeMin(), query.rangeMax().next())
        .forEach((cell, weight) -> sum[0] += weight);
    return sum[0];
  }

  /**
   * Returns the map that contains the sum of cells in 'bases' and the parents of bases with the
   * corresponding weights.
   */
  private static Map<S2CellId, Long> sumToRoot(Map<S2CellId, Long> bases) {
    return bases.keySet().stream()
        .flatMap(cell -> IntStream.range(0, cell.level() + 1)
            .mapToObj(level -> Map.entry(cell.parent(level), bases.get(cell))))
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, Long::sum));
  }

  /** Encodes entries in the map with a continuous parent/child connection to a face cell. */
  private void checkEncoding(Map<S2CellId, Long> entries) {
    assertEquals(entries(entries), entries(encode(entries).decode()));
  }

  private static List<Entry<S2CellId, Long>> entries(Map<S2CellId, Long> cells) {
    Comparator<S2CellId> byCell = Comparator.comparingInt(S2CellId::level).thenComparing(x -> x);
    return ImmutableList.sortedCopyOf(Entry.comparingByKey(byCell), cells.entrySet());
  }
}
