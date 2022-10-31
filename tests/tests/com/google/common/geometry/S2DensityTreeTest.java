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
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder.CellWeightFunction;
import com.google.common.geometry.S2DensityTree.Cell;
import com.google.common.geometry.S2DensityTree.CellVisitor.Action;
import com.google.common.geometry.S2DensityTree.IndexCellWeightFunction;
import com.google.common.geometry.S2DensityTree.SumDensity;
import com.google.common.geometry.S2DensityTree.TreeEncoder;
import com.google.common.geometry.S2DensityTree.TreeEncoder.ReversedLengthsWriter;
import com.google.common.geometry.S2DensityTree.TreeEncoder.ReversibleBytes;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
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
  private final Map<S2CellId, Long> weights = ImmutableMap.<S2CellId, Long>builder()
      .put(FACE_CELLS[1], 3L)
          .put(FACE_CELLS[1].child(1), 1L)
          .put(FACE_CELLS[1].child(2), 2L)
      .put(FACE_CELLS[2], 4L)
      .put(FACE_CELLS[3], 2L)
          .put(FACE_CELLS[3].child(0), 2L)
      .build();

  public void testReversibleBytes() {
    ReversibleBytes bytes = new ReversibleBytes();
    bytes.write(VERSION.getBytes(StandardCharsets.ISO_8859_1));
    assertEquals(VERSION, new String(bytes.toByteArray(), StandardCharsets.ISO_8859_1));
    // Reversing twice still gives us the right value.
    bytes.reverseFrom(0);
    assertEquals(VERSION, new String(bytes.reversedCopy(), StandardCharsets.ISO_8859_1));
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
    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(ImmutableMap.of(S2CellId.fromPoint(new S2Point(1, 2, 3).normalize()), 1L))
        .forEach(encoder::put);
    S2DensityTree tree = encoder.build();
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
    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(base).forEach(encoder::put);
    S2DensityTree tree = encoder.build();
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
    Random r = new Random(0);
    TreeEncoder encoder = new TreeEncoder();
    for (long weight = 1; weight < 1000; weight++) {
      Map<S2CellId, Long> expected = new HashMap<>();
      for (int j = 0; j < 50; j++) {
        expected.put(randomCell(r), weight);
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
            .build(),
        tree.decode());
  }

  // Tests accessing tree nodes with DecodedPath.
  public void testDecodedPath() {
    // Create an S2DensityTree with the given leaf cells and weights.
    S2CellId cell22 = S2CellId.fromFace(2).child(2);
    ImmutableMap<S2CellId, Long> leaves = ImmutableMap.of(
        cell22.child(2), 100L,
        cell22.child(3), 120L);

    TreeEncoder encoder = new TreeEncoder();
    sumToRoot(leaves).forEach(encoder::put);
    S2DensityTree tree = encoder.build();

    S2DensityTree.DecodedPath decoder = new S2DensityTree.DecodedPath();
    decoder.set(tree);

    // Check that the leaves can be obtained by DecodedPath.
    for (Map.Entry<S2CellId, Long> leaf : leaves.entrySet()) {
      S2DensityTree.Cell cell = decoder.cell(leaf.getKey());
      assertNotNull(cell);
      assertEquals(cell.weight, (long) leaf.getValue());
    }

    // Face nodes not in the tree should return an empty cell.
    assertTrue(decoder.cell(S2CellId.fromFace(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(1)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(5)).isEmpty());

    // Face 2 should be present with the sum of the leaf weights.
    assertEquals(220, decoder.cell(S2CellId.fromFace(2)).weight);

    // Level 1 nodes of a face that isn't in the tree
    assertTrue(decoder.cell(S2CellId.fromFace(1).child(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(1).child(1)).isEmpty());

    // Level 1 nodes of the face that is in the tree
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(0)).isEmpty());
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(1)).isEmpty());
    assertEquals(220, decoder.cell(S2CellId.fromFace(2).child(2)).weight);
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(3)).isEmpty());

    // Level 2 nodes
    assertTrue(decoder.cell(S2CellId.fromFace(2).child(3).child(0)).isEmpty());

    assertTrue(decoder.cell(cell22.child(0)).isEmpty());
    assertTrue(decoder.cell(cell22.child(1)).isEmpty());
    assertEquals(100, decoder.cell(cell22.child(2)).weight);
    assertEquals(120, decoder.cell(cell22.child(3)).weight);

    // Level 3 nodes are not present.
    assertTrue(decoder.cell(cell22.child(1).child(0)).isEmpty());
    assertTrue(decoder.cell(cell22.child(1).child(1)).isEmpty());
    assertTrue(decoder.cell(cell22.child(1).child(2)).isEmpty());
    assertTrue(decoder.cell(cell22.child(1).child(3)).isEmpty());

    assertTrue(decoder.cell(cell22.child(3).child(0)).isEmpty());
    assertTrue(decoder.cell(cell22.child(3).child(1)).isEmpty());
    assertTrue(decoder.cell(cell22.child(3).child(2)).isEmpty());
    assertTrue(decoder.cell(cell22.child(3).child(3)).isEmpty());
  }

  /** Verifies that we get no density from no inputs. */
  public void testSumEmpty() {
    checkSum(ImmutableMap.of());
  }

  /** Verifies that we get no density from one input with no cells. */
  public void testSumOne() {
    checkSum(
        ImmutableMap.of(FACE_CELLS[1], 3L,
            FACE_CELLS[1].child(1), 1L,
            FACE_CELLS[1].child(2), 2L),
        FACE_CELLS[1]);
  }

  /** Verifies that we get proper results when one tree contains the other. */
  public void testSumNested() {
    checkSum(
        // 3, 1, 2 + 1, 1, 0 = 4, 2, 2
        ImmutableMap.of(FACE_CELLS[1], 4L,
            FACE_CELLS[1].child(1), 2L,
            FACE_CELLS[1].child(2), 2L),
        FACE_CELLS[1], FACE_CELLS[1].child(1));
  }

  /** Verifies two disjoint trees don't affect each other. */
  public void testSumDisjoint() {
    checkSum(
        ImmutableMap.of(
            FACE_CELLS[2], 4L,
            FACE_CELLS[3], 2L, FACE_CELLS[3].child(0), 2L),
        FACE_CELLS[2], FACE_CELLS[3]);
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

  /** Creates a density tree for each root using the given weights, and verifies the sum. */
  private void checkSum(Map<S2CellId, Long> expected, S2CellId ... roots) {
    TreeEncoder encoder = new TreeEncoder();
    List<S2DensityTree> trees = new ArrayList<>();
    for (S2CellId root : roots) {
      // Insert weights in 'weights' under the given roots.
      insert(encoder, root);
      // And then if this isn't a leaf cell, insert the root weight into all cells above it.
      long weight = weights.get(root);
      while (root.level() > 0) {
        root = root.parent();
        encoder.put(root, weight);
      }
      trees.add(encoder.build());
    }
    SumDensity sum = new SumDensity();
    sum.set(trees);
    assertEquals(expected, new BreadthFirstTreeBuilder(1_000, 20, encoder).build(sum).decode());
  }

  private void insert(TreeEncoder encoder, S2CellId cell) {
    if (weights.containsKey(cell)) {
      encoder.put(cell, weights.get(cell));
      for (int i = 0; i < 4; i++) {
        insert(encoder, cell.child(i));
      }
    }
  }

  /** Verifies cells with weight intersect, and unsplit nodes with weight are contained. */
  public void testDensityRegion() {
    TreeEncoder encoder = new TreeEncoder();
    weights.forEach(encoder::put);
    S2DensityTree tree = encoder.build();
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
    Random r = new Random(0);
    TreeMap<S2CellId, Long> weights = new TreeMap<>();
    for (int i = 0; i < 1_000; i++) {
      S2CellId parent = randomCell(r);
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

  /** Returns the sum of weights that intersect 'query'. */
  private static long sum(SortedMap<S2CellId, Long> weights, S2CellId query) {
    long[] sum = {0};
    weights.subMap(query.rangeMin(), query.rangeMax().next())
        .forEach((cell, weight) -> sum[0] += weight);
    return sum[0];
  }

  private static S2CellId randomCell(Random r) {
    int face = r.nextInt(6);
    long pos = r.nextLong() & S2CellId.POS_BITS;
    int level = r.nextInt(30);
    return S2CellId.fromFacePosLevel(face, pos, level);
  }

  /**
   *  Returns the map that contains the sum of cells in 'bases' and the parents of bases with the
   *  corresponding weights.
   */
  private static Map<S2CellId, Long> sumToRoot(Map<S2CellId, Long> bases) {
    return bases.keySet().stream()
        .flatMap(cell -> IntStream.range(0, cell.level() + 1)
            .mapToObj(level -> Map.entry(cell.parent(level), bases.get(cell))))
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, Long::sum));
  }

  /** Encodes entries in the map with a continuous parent/child connection to a face cell. */
  private void checkEncoding(Map<S2CellId, Long> entries) {
    TreeEncoder encoder = new TreeEncoder();
    entries.forEach(encoder::put);
    assertEquals(entries(entries), entries(encoder.build().decode()));
  }

  private static List<Entry<S2CellId, Long>> entries(Map<S2CellId, Long> cells) {
    Comparator<S2CellId> byCell = Comparator.comparingInt(S2CellId::level).thenComparing(x -> x);
    return ImmutableList.sortedCopyOf(Entry.comparingByKey(byCell), cells.entrySet());
  }
}
