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

import static com.google.common.geometry.EncodedInts.varIntSize;
import static com.google.common.geometry.EncodedInts.writeVarint64;
import static com.google.common.geometry.S2CellId.FACE_CELLS;
import static com.google.common.geometry.S2CellId.MAX_LEVEL;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.toIntExact;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.collect.Maps;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2DensityTree.BreadthFirstTreeBuilder.CellWeightFunction;
import com.google.common.geometry.S2ShapeIndexRegion.ShapeVisitor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.IOException;
import java.io.OutputStream;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.ToLongFunction;
import jsinterop.annotations.JsEnum;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * A density tree is a kind of spatial histogram, computed over a collection of shapes in an
 * S2ShapeIndex. It's most often used to do spatial clustering into equal-size shards, where spatial
 * datasets are well-known for their skew, which would defeat a regular tiling of the surface. But
 * with detailed spatial density, we can easily divide the surface into equal-weight shards.
 *
 * <p>A density tree is more formally a map from S2CellIds to weights of objects that intersect or
 * contain that cell, with these properties:
 *
 * <ul>
 *   <li>If any cell is present in the tree, then all of its ancestors are present in the tree.
 *   <li>If a cell has any descendants in the tree, then the descendants together intersect and
 *       include the weight of all objects intersecting the parent.
 * </ul>
 *
 * <p>Some consequences of these properties are:
 *
 * <ul>
 *   <li>If a shape intersects a density tree, it intersects at least one leaf cell in the tree.
 *   <li>All the shapes are covered by leaf cells of the tree.
 * </ul>
 *
 * <p>Note that the weight of a parent cell is generally not the sum of the weights of its children.
 * It may be less, as an object may intersect multiple children.
 *
 * <p>Weights have no units. Their meaning is defined by the operation that computes them. For
 * example, the weight of a shape could be the number of vertices in each shape. Generally a client
 * will compute {@link #shapeDensity} for a whole set of shapes at once by putting them in an {@link
 * S2ShapeIndex} and providing a {@link ShapeWeightFunction} to map a shape to the weight of that
 * shape. If there are multiple indices to weigh, intermediate trees may be combined with {@link
 * #sumDensity} or {@link #intersectionDensity}. The convenience function {@link #vertexDensity} is
 * useful when the measure of a shape's weight is simply the number of vertices in the shape.
 *
 * <p>As a limitation, if a shape overlaps several cells at a given level, S2DensityTree does not
 * keep track of the fact that the weights in those cells are correlated. So while it is capable of
 * accurately measuring the weight in individual cells, it can very significantly overcount the
 * weight in a collection of cells. This can cause the S2CellIds returned by {@link
 * #getPartitioning(int)} to have significantly less than the targeted weight.
 *
 * <p>It's strongly encouraged to use the bulk versions of the above methods if they will be called
 * more than once, as they involve a lot of temporary allocations that may be pooled by {@link
 * #createShapeDensityOp}, {@link #createVertexDensityOp}, {@link #createSumDensityOp}, and {@link
 * #createIntersectionDensityOp}.
 *
 * <p>Density trees are fast and cheap to initialize from a {@link Bytes buffer}, but the entries
 * are parsed lazily. If lookups will visited more often than there are entries, it may be faster to
 * {@link #decode} into a cell,weight map that has faster lookups, but note the in-memory size is
 * several times larger.
 *
 * <p>S2DensityTrees can be created with a restricted size threshold:
 *
 * {@snippet :
 * S2ShapeIndex index = new S2ShapeIndex();
 * index.add(...);
 * S2DensityTree tree = S2DensityTree.vertexDensity(index, 10_000, 15);
 * }
 *
 * <p>An example of summing trees:
 *
 * {@snippet :
 * List<S2DensityTree> trees = ...;
 * S2DensityTree sumTree = S2DensityTree.sumDensity(trees, 100_000, 15);
 * }
 *
 * <p>Example of getting approximately 100MiB disjoint shards from a given tree (assuming weights
 * are bytes):
 *
 * {@snippet :
 * S2DensityTree tree = ...;
 * List<S2CellId> partitions = tree.getPartitioning(100 * 2 << 10);
 * }
 */
@JsType
@SuppressWarnings("Assertion")
public class S2DensityTree {
  /** The current version string that prefixes all density trees, in the latin1 charset. */
  public static final String VERSION = "S2DensityTree0";

  /** The {@link #VERSION} as latin1 bytes. */
  private static final byte[] VERSION_BYTES = VERSION.getBytes(StandardCharsets.ISO_8859_1);

  /** The {@link #VERSION_BYTES} backwards since we write them this way before reversing. */
  private static final byte[] REVERSED_VERSION_BYTES = new byte[VERSION_BYTES.length];

  static {
    for (int i = 0, j = VERSION_BYTES.length - 1; j >= 0; i++, j--) {
      REVERSED_VERSION_BYTES[i] = VERSION_BYTES[j];
    }
  }

  /** An {@link S2Coder} of density trees. */
  public static final S2Coder<S2DensityTree> CODER =
      new S2Coder<S2DensityTree>() {
        @JsIgnore // OutputStream is not available in J2CL.
        @Override
        public void encode(S2DensityTree value, OutputStream output) throws IOException {
          value.encoded.writeTo(output);
        }

        @Override
        public S2DensityTree decode(Bytes data, Cursor cursor) throws IOException {
          return S2DensityTree.decode(data, cursor);
        }

        @Override
        public boolean isLazy() {
          return true;
        }
      };

  private static final int CHILD_MASK_BITS = 4;
  private static final int CHILD_MASK = (1 << CHILD_MASK_BITS) - 1;
  private static final long MAX_WEIGHT = Long.MAX_VALUE >>> CHILD_MASK_BITS;

  private final Bytes encoded;
  private final long[] facePositions;

  /** Creates a new instance, where 'encoded' must have no other references. */
  private S2DensityTree(Bytes encoded, Cursor cursor) throws IOException {
    this.encoded = encoded;
    this.facePositions = decodeHeader(encoded, cursor);
  }

  /** Returns the encoded size of this tree in bytes. */
  public long size() {
    return encoded.length();
  }

  /** Returns the bytes underlying this tree. This is a cheap operation. */
  public Bytes encoded() {
    return encoded;
  }

  /** As {@link #decode(Bytes, Cursor)} with a default cursor over the whole set of bytes. */
  @JsIgnore // No method overloading in J2CL. Use decode(bytes, bytes.cursor()).
  public static S2DensityTree decode(Bytes bytes) throws IOException {
    return new S2DensityTree(bytes, bytes.cursor());
  }

  /** Reads a tree from the given source of bytes. This does not fully decode the tree. */
  public static S2DensityTree decode(Bytes bytes, Cursor cursor) throws IOException {
    return new S2DensityTree(bytes, cursor);
  }

  /**
   * Returns a fully-decoded map of this tree. This is only useful if far more lookups will be done
   * than there are entries in the map, and the lookups are sufficiently random that a {@link
   * DecodedPath} is not sufficient. In that case fully decoding all the cells may be faster than
   * lazily-decoding cells on each lookup.
   *
   * <p>Since the result is space-inefficient, considering using visitors to collect decoded results
   * instead of this method.
   */
  public Map<S2CellId, Long> decode() {
    TreeMap<S2CellId, Long> results = new TreeMap<>();
    visitAll(results::put);
    return results;
  }

  /** Visits all nodes in this tree in depth-first order. */
  public void visitAll(CellVisitor.All visitor) {
    visitCells(visitor);
  }

  /**
   * Visits nodes in depth-first order while the visitor returns true, and returns false iff visitor
   * returns false.
   */
  public boolean visitWhile(CellVisitor.While visitor) {
    return visitCells(
        (cell, node) ->
            visitor.visit(cell, node.weight)
                ? CellVisitor.Action.ENTER_CELL
                : CellVisitor.Action.STOP);
  }

  /** Returns an S2Region view of the density tree. */
  public S2Region asRegion() {
    DecodedPath path = new DecodedPath(this);
    return new S2Region() {
      @Override
      public S2Cap getCapBound() {
        return S2Cap.full();
      }

      @Override
      public S2LatLngRect getRectBound() {
        return S2LatLngRect.full();
      }

      @Override
      public boolean mayIntersect(S2Cell cell) {
        // If the cell is in the tree, or its last ancestor in the tree is a leaf, then the objects
        // weighed by the tree may intersect the given cell. Otherwise, the spatial extent of the
        // objects weighed by the tree cannot intersect the cell.
        Cell node = path.cell(cell.id());
        return !node.isEmpty();
      }

      @Override
      public boolean contains(S2Cell cell) {
        return contains(cell.id());
      }

      @Override
      public boolean contains(S2Point p) {
        return contains(S2CellId.fromPoint(p));
      }

      private boolean contains(S2CellId id) {
        // If the cell for the given id (or an ancestor) is in the tree and is a leaf, it is
        // contained by the objects weighed by the tree (or the tree has insufficient detail to
        // differentiate between containment and intersection, which we treat as containment).
        // Otherwise the node is either empty, indicating the cell id is disjoint from the tree, or
        // it has children, indicating the cell id is not contained by all the objects weighed by
        // the tree.
        Cell node = path.cell(id);
        return !node.isEmpty() && !node.hasChildren();
      }
    };
  }

  /**
   * Density trees map cells to weights that intersect that cell, but larger features may intersect
   * many density cells, so the deeper into the tree we look, the more a feature will contribute
   * undetected duplicates of its weight many times.
   *
   * <p>So we offer building a normalized version of the tree where every node's weight is scaled by
   * (its parent's weight / the sum of weights of the node and its siblings). This makes the weight
   * of a parent equal to the sum of its children.
   */
  public S2DensityTree normalize() {
    // Visit cells in depth-first order.
    Map<S2CellId, Long> normalizedWeights = new HashMap<>();
    DecodedPath path = new DecodedPath(this);
    visitAll(
        (cell, weight) -> {
          if (!cell.isFace()) {
            S2CellId parent = cell.parent();
            long parentWeight = normalizedWeights.get(parent);
            long siblingWeight = 0;
            for (int i = 0; i < 4; i++) {
              siblingWeight += path.weight(parent.child(i));
            }
            // Rescale weight by the ratio of parentWeight/siblingWeight, which will be less than 1.
            if (parentWeight < Long.MAX_VALUE / weight) {
              // The product weight*parentWeight can't overflow, so use the fast formula.
              weight = (weight * parentWeight - 1) / siblingWeight + 1;
            } else {
              // The product very likely overflows, so use the more expensive formula.
              // This is likely to be used in the first few levels of a density tree for large
              // datasets, but the fast formula should be used for most cells.
              weight =
                  BigInteger.valueOf(weight)
                      .multiply(BigInteger.valueOf(parentWeight))
                      .subtract(BigInteger.ONE)
                      .divide(BigInteger.valueOf(siblingWeight))
                      .add(BigInteger.ONE)
                      .longValue();
            }
          }
          normalizedWeights.put(cell, weight);
        });
    TreeEncoder encoder = new TreeEncoder();
    normalizedWeights.forEach(encoder::put);
    return encoder.build();
  }

  /**
   * Returns a selection of disjoint cells with generally 25% to 100% of 'maxWeight' weight per
   * cell. The cells as a set cover the spatial extent of the objects whose weight is captured in
   * the tree. If the density tree's leaves are ever heavier than 'maxWeight' then those cells will
   * be heavier.
   */
  public List<S2CellId> getPartitioning(int maxWeight) {
    List<S2CellId> clusters = new ArrayList<>();
    visitCells(
        (cell, node) -> {
          if (node.weight > maxWeight && node.hasChildren()) {
            return CellVisitor.Action.ENTER_CELL;
          } else {
            clusters.add(cell);
            return CellVisitor.Action.SKIP_CELL;
          }
        });
    return clusters;
  }

  /** Returns the S2CellIds of the leaves of this tree as an S2CellUnion. */
  public S2CellUnion getLeaves() {
    S2CellUnion result = new S2CellUnion();
    result.initRawCellIds(getLeavesAsList());
    return result;
  }

  /** Returns the S2CellIds of the leaves of this tree as an ArrayList. */
  public ArrayList<S2CellId> getLeavesAsList() {
    ArrayList<S2CellId> leaves = new ArrayList<>();
    visitCells(
        (cell, node) -> {
          if (node.hasChildren()) {
            return CellVisitor.Action.ENTER_CELL;
          }
          leaves.add(cell);
          return CellVisitor.Action.SKIP_CELL;
        });
    return leaves;
  }

  /** Returns the maximum level of any node in this tree. */
  public int getMaxLevel() {
    final int[] maxLevel = new int[1]; // Array is workaround for Java lambdas requiring final
    visitCells(
        (cell, node) -> {
          maxLevel[0] = max(maxLevel[0], cell.level());
          return CellVisitor.Action.ENTER_CELL;
        });
    return maxLevel[0];
  }

  /**
   * Visits each cell,weight pair in depth-first order, halting early and returning false if {@link
   * CellVisitor#visit} ever returns {@link CellVisitor.Action#STOP}.
   */
  @CanIgnoreReturnValue
  public boolean visitCells(CellVisitor visitor) {
    Cursor cursor = encoded.cursor();
    Cell cell = new Cell();
    for (int i = 0; i < FACE_CELLS.length; i++) {
      cursor.position = facePositions[i];
      if (cursor.position >= 0 && !visitCells(cell, cursor, FACE_CELLS[i], visitor)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Decodes and visits cells beginning from the given 'cell' in depth-first order. Returns true if
   * we should continue to visit cells.
   */
  private boolean visitCells(Cell decoded, Cursor cursor, S2CellId cell, CellVisitor visitor) {
    decoded.decode(encoded, cursor);
    switch (visitor.visit(cell, decoded)) {
      case SKIP_CELL:
        return true;
      case ENTER_CELL:
        // Copy the positions before we recurse so we can reuse one Cell.
        int p0 = decoded.positions[0];
        int p1 = decoded.positions[1];
        int p2 = decoded.positions[2];
        int p3 = decoded.positions[3];
        // Visit the children. Skip unset children. Stop early if requested from child processing.
        // Rely on short-circuited logic to halt early if any visit call returns false.
        return (p0 < 0 || visitCells(decoded, cursor.seek(p0), cell.child(0), visitor))
            && (p1 < 0 || visitCells(decoded, cursor.seek(p1), cell.child(1), visitor))
            && (p2 < 0 || visitCells(decoded, cursor.seek(p2), cell.child(2), visitor))
            && (p3 < 0 || visitCells(decoded, cursor.seek(p3), cell.child(3), visitor));
      case STOP:
        return false;
      default:
        throw new IllegalStateException("Unexpected next state");
    }
  }

  /** A visitor to stream through the cell,weight pairs of the density tree in depth-first order. */
  @JsType
  public interface CellVisitor {
    /**
     * @param cell the cell
     * @param node the weight and set of available children in 'cell'
     * @return whether to skip children, visit children, or stop altogether
     */
    Action visit(S2CellId cell, Cell node);

    /** The action requested by the visitor. */
    @JsEnum
    enum Action {
      /** Continue visitation but skip the children of this node. */
      SKIP_CELL,
      /** Continue visitation with the children of this node. */
      ENTER_CELL,
      /** Abort visitation. {@link S2DensityTree#visitCells(CellVisitor)} will return false. */
      STOP
    }

    /** A visitor of all nodes. */
    @JsType
    interface All extends CellVisitor {
      @Override
      default Action visit(S2CellId cell, Cell node) {
        visit(cell, node.weight);
        return Action.ENTER_CELL;
      }

      /** Visits the next cell. */
      @JsIgnore // J2CL has no method overloads.
      void visit(S2CellId cell, long weight);
    }

    /** A visitor of nodes while {@link #visit} returns true. */
    @JsType
    interface While extends CellVisitor {
      @Override
      default Action visit(S2CellId cell, Cell node) {
        return visit(cell, node.weight) ? Action.ENTER_CELL : Action.STOP;
      }

      /** Visits the next cell. */
      @JsIgnore // J2CL has no method overloads.
      boolean visit(S2CellId cell, long weight);
    }
  }

  /**
   * A function that returns the weight of S2Shapes.
   *
   * <p>See {@link S2DensityTree class comments} for the meaning of weights).
   */
  @JsType
  public interface ShapeWeightFunction extends ToLongFunction<S2Shape> {
    /**
     * Returns the weight of the given shape, which must be in the closed range [0, {@link
     * #MAX_WEIGHT}].
     */
    @Override
    long applyAsLong(S2Shape shape);

    /**
     * Resets this weigher for the next cell. This is useful to override if multiple shapes map to
     * the same underlying row and we want to measure the row once per cell.
     */
    default void reset() {}
  }

  /**
   * Returns the density tree of a given index, with the given weigher of shapes in the index. If
   * multiple indices must have their density computed, use {@link #createShapeDensityOp} to get a
   * reusable {@link ShapeDensityOp} that is often much faster for repeated calls.
   *
   * @param approximateSizeBytes an approximate limit in the tree size, which is grown an entire
   *     level at a time until this limit is reached, so in practice the actual size may be up to 4x
   *     this limit in the worst-case
   * @param maxLevel the max level to create in the tree, where most paths in the tree will reach
   *     this limit except cells contained by polygons, which can save a lot of space in polygonal
   *     datasets
   */
  public static S2DensityTree shapeDensity(
      S2ShapeIndex index, ShapeWeightFunction weigher, int approximateSizeBytes, int maxLevel) {
    return createShapeDensityOp(approximateSizeBytes, maxLevel).apply(index, weigher);
  }

  /**
   * As {@link #shapeDensity}, where the weigher is {@link S2ShapeUtil#numVertices}. If multiple
   * indices must have their density computed, use {@link #createVertexDensityOp} to get a reusable
   * {@link VertexDensityOp} that is often much faster for repeated calls.
   *
   * @param approximateSizeBytes an approximate limit in the tree size, which is grown an entire
   *     level at a time until this limit is reached, so in practice the actual size may be up to 4x
   *     this limit in the worst-case
   * @param maxLevel the max level to create in the tree, where most paths in the tree will reach
   *     this limit except cells contained by polygons, which can save a lot of space in polygonal
   *     datasets
   */
  public static S2DensityTree vertexDensity(
      S2ShapeIndex index, int approximateSizeBytes, int maxLevel) {
    return createVertexDensityOp(approximateSizeBytes, maxLevel).apply(index);
  }

  /**
   * Returns the feature density, where each feature that intersects a cell is weighed once whether
   * it has one shape that intersects the feature or many.
   *
   * @param shapes an index that gets the shapes that intersect a cell to weigh
   * @param features a lookup function that gets the feature corresponding to each shape
   * @param weights a lookup function that gets the weight of each feature
   * @param <F> the type of feature
   */
  public static <F> S2DensityTree featureDensity(
      int approximateSizeBytes,
      int maxLevel,
      S2ShapeIndex shapes,
      FeatureLookup<F> features,
      FeatureWeigher<F> weights) {
    FeatureDensityOp<F> op = createFeatureDensityOp(approximateSizeBytes, maxLevel);
    return op.apply(shapes, features, weights);
  }

  /**
   * Given an input S2DensityTree, returns a "dilated" tree, i.e. a new density tree that has weight
   * in all areas within the given distance "radius" of the spatial boundary of the leaves of the
   * tree.
   *
   * <p>Dilation is required for creating density trees for clustering data for cases like a
   * D-Within filter join, so that clusters created from the density tree will contain all pairs of
   * features within distance "radius" of each other.
   *
   * <p>Dilation is done by adding weights to neighbors of leaves of the tree. If the radius is
   * small relative to the area covered by the density tree, the required size of the neighboring
   * cells could also be small, resulting in an exponential blowup in the number of cells in the
   * tree. The 'maxLevelDiff' parameter prevents this by limiting the level of neighboring added
   * cells, relative to the size of the leaf being dilated. Therefore there is a tradeoff between
   * dilation accuracy and tree size: a higher maxLevelDiff will add more, smaller cells.
   *
   * <p>If the dilation radius is large compared to cells in the tree, then dilation discards nodes
   * with a cell level higher than the cell level required for dilation, and instead adds neighbors
   * to nodes in the tree that are at the dilation cell level.
   */
  public static S2DensityTree dilate(S2DensityTree tree, S1Angle radius, int maxLevelDiff) {
    // Get the leaves of the density tree as an S2CellUnion, convert them to an S2Polygon, then to
    // an S2LaxPolylineShape of the boundary of the tree. Then put the boundary in an S2ShapeIndex
    // and finally get a S2ShapeIndexRegion. Checking intersection with the boundary region is used
    // to determine which nodes in the density tree are on the spatial boundary, and thus need to be
    // dilated.
    S2CellUnion leaves = tree.getLeaves();

    // Create an empty Map<S2CellId, Long> weights in which we will compute the new density tree.
    final Map<S2CellId, Long> weights = new HashMap<>();

    // 'radiusLevel' is the highest cell level (smallest possible cells) that could be used for
    // dilation that will cover all areas within 'radius' of the tree. Existing cells in the tree
    // at higher levels will be dropped.
    int radiusLevel = S2Projections.MIN_WIDTH.getMaxLevel(radius.radians());

    // Compute an S2CellUnion 'dilationCells' of the 'halo' around the existing leaves that will
    // contain the cells added by dilation, and does not contain the existing leaves. Note that
    // expansion may use a level lower than radiusLevel.
    S2CellUnion expandedLeaves = S2CellUnion.copyFrom(leaves);
    expandedLeaves.expand(radius, maxLevelDiff);
    S2CellUnion dilationCells = new S2CellUnion();
    dilationCells.getDifference(expandedLeaves, leaves);

    // Visit all the cells in the current tree, copying them to the dilated tree, and where required
    // also adding neighboring cells for dilation.
    tree.visitCells(
        (S2CellId cellId, Cell node) -> {
          // Add this existing node to the output tree. If it is already present in the tree (as a
          // neighbor of a previously visited and dilated cell), update its dilated weight to the
          // maximum of the current dilated weight and the pre-dilation weight.
          long dilatedWeight = max(weights.getOrDefault(cellId, 0L), node.weight());
          weights.put(cellId, dilatedWeight);

          // If this node has children, and has a level less than radiusLevel, it will become an
          // internal node of the dilated tree, so dilation will occur at its children.
          if (node.hasChildren() && cellId.level() < radiusLevel) {
            return CellVisitor.Action.ENTER_CELL;
          }

          // Otherwise, this cell is either a leaf of the input tree, or has level "radiusLevel".
          // Either way, it will be a leaf of the dilated tree. It may need to be dilated.
          // (Determining precisely which leaves must be dilated could be done better at the cost of
          // some additional complexity: the current implementation does more dilation than needed.)

          // To dilate this output leaf, we want to use the highest level possible for accuracy,
          // which is radiusLevel, but the cell level to use for dilating may not be more than the
          // current cell level plus maxLevelDiff.
          int dilateLevel = min(radiusLevel, maxLevelDiff + cellId.level());

          // For each neighbor that intersects the dilation area, and is not already present with
          // sufficient weight, insert it or increase the weight to 'dilatedWeight'. Also ensure
          // that all its ancestors also have weight at least 'dilatedWeight' so the tree is valid.
          cellId.visitNeighbors(
              dilateLevel,
              nbr -> {
                if (!dilationCells.intersects(nbr)) {
                  // This neighbor isn't in the area of dilation, ignore it but keep going.
                  return true;
                }
                if (weights.getOrDefault(nbr, 0L) < dilatedWeight) {
                  weights.put(nbr, dilatedWeight);
                  for (;
                      nbr.level() > 0 && weights.getOrDefault(nbr.parent(), 0L) < dilatedWeight;
                      nbr = nbr.parent()) {
                    weights.put(nbr.parent(), dilatedWeight);
                  }
                }
                return true;
              });

          return CellVisitor.Action.SKIP_CELL;
        });

    // Convert the map into a density tree.
    TreeEncoder encoder = new TreeEncoder();
    weights.forEach(encoder::put);
    return encoder.build();
  }

  /**
   * Returns a new density tree that contains the summed weights across the cells in the given input
   * trees. This sum may be a lossy operation if any of the input trees are more detailed than the
   * other limit parameters allow. If this will be called multiple times, use {@link
   * #createSumDensityOp} to get a reusable {@link AggregateDensityOp} that is often much faster for
   * repeated calls.
   *
   * @param approximateSizeBytes an approximate limit in the tree size, which is grown an entire
   *     level at a time until this limit is reached, so in practice the actual size may be up to 4x
   *     this limit in the worst-case
   * @param maxLevel the max level to create in the tree, where most paths in the tree will reach
   *     this limit except cells contained by polygons, which can save a lot of space in polygonal
   *     datasets
   */
  @JsIgnore // Iterable<S2DensityTree> "is not usable by JavaScript" but not clear why.
  public static S2DensityTree sumDensity(
      Iterable<S2DensityTree> trees, int approximateSizeBytes, int maxLevel) {
    return createSumDensityOp(approximateSizeBytes, maxLevel).apply(trees);
  }

  /**
   * Returns a new density tree that contains summed weights of cells in the given input trees where
   * *all* the input trees contribute non-zero weight *at their leaves* to every output cell. Cells
   * where any input tree has no weight at the leaves are pruned from the output tree.
   *
   * <p>For example, given two input trees, one of which covers a collection of features in Portland
   * and one of which covers a collection of features in Seattle, the intersection of the two trees
   * is empty, even though both trees have the level 5 cell that covers both cities in common, along
   * with all the ancestors of that cell.
   *
   * <p>This sum may be a lossy operation if any of the input trees are more detailed than the other
   * limit parameters allow. If this will be called multiple times, use {@link
   * #createIntersectionDensityOp} to get a reusable {@link AggregateDensityOp} that is often much
   * faster for repeated calls.
   *
   * @param trees S2DensityTrees to combine
   * @param approximateSizeBytes an approximate limit in the tree size, which is grown an entire
   *     level at a time until this limit is reached, so in practice the actual size may be up to 4x
   *     this limit in the worst-case
   * @param maxLevel the max level to create in the tree, where most paths in the tree will reach
   *     this limit except cells contained by polygons, which can save a lot of space in polygonal
   *     datasets
   */
  @JsIgnore // Iterable<S2DensityTree> "is not usable by JavaScript" but not clear why.
  public static S2DensityTree intersectionDensity(
      Iterable<S2DensityTree> trees, int approximateSizeBytes, int maxLevel) {
    return createIntersectionDensityOp(trees, approximateSizeBytes, maxLevel).apply(trees);
  }

  /**
   * A density function captures state needed to compute density and is not thread-safe, but greatly
   * accelerates repeated density computations with the same settings.
   */
  @JsType
  public interface ShapeDensityOp
      extends BiFunction<S2ShapeIndex, ShapeWeightFunction, S2DensityTree> {
    /** Returns the density of {@code index} using {@code weigher} at each cell. */
    @Override
    S2DensityTree apply(S2ShapeIndex index, ShapeWeightFunction weigher);
  }

  /**
   * A vertex density function captures state needed to compute vertex density and is not thread-
   * safe, but greatly accelerates repeated density computations with the same settings.
   */
  @JsType
  public interface VertexDensityOp extends Function<S2ShapeIndex, S2DensityTree> {
    /** Returns the density of {@code index}, weighing the number of vertices at each cell. */
    @Override
    S2DensityTree apply(S2ShapeIndex index);
  }

  /** A density function that collects the sum of feature weights with any shapes in each cell. */
  @JsType
  public interface FeatureDensityOp<F> {
    /**
     * Returns the density of {@code shapes}, by converting shapes that intersect each cell to
     * features with {@code features}, and weighing each <strong>instance</strong> of the feature
     * with {@code weights} just once.
     */
    S2DensityTree apply(S2ShapeIndex shapes, FeatureLookup<F> features, FeatureWeigher<F> weights);
  }

  /** A density function that aggregates values in each cell in the input trees into a new tree. */
  @JsType
  public interface AggregateDensityOp extends Function<Iterable<S2DensityTree>, S2DensityTree> {
    @Override
    S2DensityTree apply(Iterable<S2DensityTree> trees);
  }

  /**
   * As {@link #shapeDensity} but returns a thread-unsafe function that pools its allocations, so
   * testing many indices is faster than repeatedly calling shapeDensity.
   */
  public static ShapeDensityOp createShapeDensityOp(int approximateSizeBytes, int maxLevel) {
    BreadthFirstTreeBuilder acc = new BreadthFirstTreeBuilder(approximateSizeBytes, maxLevel);
    return (index, weigher) -> acc.build(new IndexCellWeightFunction(index, weigher));
  }

  /**
   * As {@link #vertexDensity} but returns a thread-unsafe function that pools its allocations, so
   * computing the density of many indices is faster than repeatedly calling vertexDensity.
   */
  public static VertexDensityOp createVertexDensityOp(int approximateSizeBytes, int maxLevel) {
    BreadthFirstTreeBuilder acc = new BreadthFirstTreeBuilder(approximateSizeBytes, maxLevel);
    IdentityHashMap<S2Shape, Integer> vertices = new IdentityHashMap<>();
    return index -> {
      vertices.clear();
      for (S2Shape s : index.getShapes()) {
        vertices.put(s, S2ShapeUtil.numVertices(s));
      }
      return acc.build(new IndexCellWeightFunction(index, vertices::get));
    };
  }

  /**
   * As {@link #featureDensity} but returns a thread-unsafe function that pools its allocations, so
   * computing the density of many indices is faster than repeatedly calling featureDensity.
   */
  public static <F> FeatureDensityOp<F> createFeatureDensityOp(
      int approximateSizeBytes, int maxLevel) {
    BreadthFirstTreeBuilder acc = new BreadthFirstTreeBuilder(approximateSizeBytes, maxLevel);
    FeatureDensityWeigher weigher = new FeatureDensityWeigher();
    return (shapes, features, weights) -> {
      weigher.init(shapes, features, weights);
      return acc.build(weigher);
    };
  }

  /** A lookup API that finds the feature corresponding to a shape. */
  @JsType
  public interface FeatureLookup<T> {
    /** Returns the feature for {@code shape}. */
    T feature(S2Shape shape);
  }

  /** A weigher of a given feature. */
  @JsType
  public interface FeatureWeigher<T> {
    /** Returns the positive non-zero weight of {@code feature}. */
    long weight(T feature);
  }

  /** A weigher of feature weights, where each feature may have multiple shapes. */
  @JsType
  public static class FeatureDensityWeigher implements CellWeightFunction {
    /** The region to test each cell against. */
    private S2ShapeIndexRegion region;

    /** Maps shape IDs to the feature ID that owns the shape, i.e. {@code id = ids[shapeId]}. */
    private int[] ids;

    /** Maps feature IDs to weights, i.e. {@code weight = weights[ids[shapeId]]}. */
    private int[] weights;

    /** The last call each ID was seen in, i.e. {@code call = lastCall[ids[shapeId]]}. */
    private int[] lastCall;

    /** The next call number. Increment to separate previous results from new results. */
    private int nextCall;

    /** The sum of weight seen in the current cell. */
    private long sum;

    /** If true then all features contained the current cell. */
    private boolean allContained;

    /**
     * Initialize this weigher.
     *
     * @param shapes an index that gets the shapes that intersect a cell to weigh
     * @param featureLookup a lookup function that gets the features corresponding to each shape
     * @param weightLookup a lookup function that gets the weight of each feature
     */
    public <T> void init(
        S2ShapeIndex shapes, FeatureLookup<T> featureLookup, FeatureWeigher<T> weightLookup) {
      int[] ids = new int[shapes.shapes.size()];
      // Assign consecutive integer ids to features in shape ID order.
      IdentityHashMap<T, Integer> map = Maps.newIdentityHashMap();
      for (int i = 0; i < ids.length; i++) {
        T feature = featureLookup.feature(shapes.shapes.get(i));
        ids[i] = map.computeIfAbsent(feature, s -> map.size());
      }
      // Compute the weights of each feature.
      int[] weights = new int[map.size()];
      map.forEach(
          (feature, id) -> {
            weights[id] = toIntExact(weightLookup.weight(feature));
          });
      initFromArrays(shapes, ids, weights);
    }

    /**
     * Inits this weigher to the given literal values. This allows minimizing copying but is more
     * difficult to use, so consider {@link #init} if such inputs are already available.
     *
     * @param shapes the shapes to test against each weighted cell
     * @param ids {@code ids[shapeID]} contains the row ID of each shape in 'shapes'
     * @param weights {@code weights[rowID]} contains the weight of each row ID in 'ids'
     */
    public void initFromArrays(S2ShapeIndex shapes, int[] ids, int[] weights) {
      this.region = new S2ShapeIndexRegion(shapes);
      this.ids = ids;
      this.weights = weights;
      this.lastCall = new int[weights.length];
      this.nextCall = 1;
    }

    @Override
    public long applyAsLong(S2CellId cell) {
      sum = 0;
      allContained = true;
      region.visitIntersectingShapes(new S2Cell(cell), this::sumFeatures);
      // sumFeatures uses lastCall[id]==nextCall to detect whether it's the first time we have seen
      // a feature *for this cell*. Rather than clear a set of IDs, just increment nextCall so all
      // lastCall elements can be seen the first time. If the nextCall value reaches its max value,
      // only then do we reset things to the starting position, but unlike Set.clear() after each
      // call to this method, rolling over the nextCall should be incredibly rare.
      if (++nextCall == Integer.MAX_VALUE) {
        nextCall = 1;
        Arrays.fill(lastCall, 0);
      }
      return allContained ? -sum : sum;
    }

    /** Adds the weight of this shape if it belongs to a feature we haven't yet seen. */
    private boolean sumFeatures(int shapeId, boolean containsTarget) {
      int featureId = ids[shapeId];
      if (lastCall[featureId] != nextCall) {
        lastCall[featureId] = nextCall;
        sum += weights[featureId];
        allContained &= containsTarget;
      }
      return true;
    }
  }

  /**
   * As {@link #sumDensity} but returns a thread-unsafe function that pools its allocations, so
   * computing the density of many sets of trees is faster than repeatedly calling sumDensity.
   */
  public static AggregateDensityOp createSumDensityOp(int approximateSizeBytes, int maxLevel) {
    BreadthFirstTreeBuilder acc = new BreadthFirstTreeBuilder(approximateSizeBytes, maxLevel);
    SumDensity sum = new SumDensity();
    return trees -> {
      sum.set(trees);
      return acc.build(sum);
    };
  }

  /**
   * Returns a thread-unsafe function that combines density trees. The function returns a tree that
   * contains the sum of the density of the trees provided to the function, but only where the sum
   * intersects the intersection of the leaves of the "intersectionTrees" provided here. Other parts
   * of the summed tree are discarded.
   *
   * <p>The resulting tree may not be as detailed as the input trees, but will be as detailed as
   * possible given the limits passed to this method and the available detail of the inputs.
   */
  public static AggregateDensityOp createIntersectionDensityOp(
      Iterable<S2DensityTree> intersectionTrees, int approximateSizeBytes, int maxLevel) {
    BreadthFirstTreeBuilder acc = new BreadthFirstTreeBuilder(approximateSizeBytes, maxLevel);
    IntersectionDensity intersection = IntersectionDensity.create(intersectionTrees);
    return trees -> {
      intersection.set(trees);
      return acc.build(intersection);
    };
  }

  /**
   * A reusable measurer of shape index density that computes the weight at each cell the index
   * intersects, as the sum of the weigher for each shape in the index that intersects a cell.
   *
   * <p>This class is not thread safe, but reusing an instance saves a lot of allocation.
   */
  public static class IndexCellWeightFunction
      implements BreadthFirstTreeBuilder.CellWeightFunction {
    private final S2ShapeIndex index;
    private final S2ShapeIndexRegion region;
    private final ShapeWeightFunction weigher;

    /**
     * @param index The index to locate shapes within that intersect cells being weighed
     * @param weigher The weigher of each shape in the given index
     */
    public IndexCellWeightFunction(S2ShapeIndex index, ShapeWeightFunction weigher) {
      this.index = index;
      this.region = new S2ShapeIndexRegion(index);
      this.weigher = weigher;
    }

    /**
     * Returns the total weight of all shapes that intersect the given cell. Negates the returned
     * total if the cell is fully contained by all shapes that intersect it.
     */
    @Override
    public long applyAsLong(S2CellId cell) {
      class IntersectingVisitor implements ShapeVisitor {
        long sum = 0;
        long limit = MAX_WEIGHT;
        boolean allContained = true;

        @Override
        public boolean test(int shapeId, boolean containsTarget) {
          long weight = weigher.applyAsLong(index.getShapes().get(shapeId));
          Preconditions.checkState(weight >= 0, "Weight must not be negative.");
          limit -= weight;
          Preconditions.checkState(limit >= 0, "Weight must be at most " + MAX_WEIGHT);
          sum += weight;
          allContained &= containsTarget;
          return true;
        }

        /** Returns the sum of the shape weights, negated if all shapes contained the cell. */
        long finalWeight() {
          return allContained ? -sum : sum;
        }
      }
      IntersectingVisitor visitor = new IntersectingVisitor();
      region.visitIntersectingShapes(new S2Cell(cell), visitor);
      weigher.reset();
      return visitor.finalWeight();
    }
  }

  /** A weigher of features described by cells. */
  public static class CellDensityFunction implements CellWeightFunction {
    private final S2CellIndex index;
    private final S2CellUnion[] coverings;
    private final long[] weights;
    private final S2CellUnion query = new S2CellUnion();

    /**
     * Creates a new cell density function. The index must have been created from 'coverings', as
     * done by {@link #fromArrays}.
     *
     * @param index the index of cell->row, which must contain 'coverings'
     * @param coverings an array of coverings indexed by label in 'index'
     * @param weights an array of weights indexed by label in 'index'
     */
    public CellDensityFunction(S2CellIndex index, S2CellUnion[] coverings, long[] weights) {
      Preconditions.checkArgument(coverings.length == weights.length);
      this.index = index;
      this.coverings = coverings;
      this.weights = weights;
    }

    /**
     * Creates a new cell density function. A cell index is created from 'coverings' with the
     * corresponding weights. If the index of 'coverings' has already been created, consider calling
     * the {@link CellDensityFunction () constructor} instead.
     *
     * @param coverings an array of coverings, pairwise with 'weights'
     * @param weights an array of weights, pairwise with 'coverings'
     */
    public static CellDensityFunction fromArrays(S2CellUnion[] coverings, long[] weights) {
      S2CellIndex index = new S2CellIndex();
      for (int i = 0; i < coverings.length; i++) {
        index.add(coverings[i], i);
      }
      index.build();
      return new CellDensityFunction(index, coverings, weights);
    }

    @Override
    public long applyAsLong(S2CellId cell) {
      query.cellIds().add(cell);
      long sum = 0;
      boolean contained = true;
      for (int label : index.getIntersectingLabels(query)) {
        sum += weights[label];
        contained &= coverings[label].contains(query);
      }
      query.cellIds().clear();
      return contained ? -sum : sum;
    }

    /** Returns a weigher of features, given functions to extract the covering and weight. */
    public static <T> CellWeightFunction fromFeatures(
        List<T> features, Function<T, S2CellUnion> covering, ToLongFunction<T> weight) {
      S2CellUnion[] coverings = new S2CellUnion[features.size()];
      long[] weights = new long[features.size()];
      for (int i = 0; i < features.size(); i++) {
        T feature = features.get(i);
        coverings[i] = covering.apply(feature);
        weights[i] = weight.applyAsLong(feature);
      }
      return fromArrays(coverings, weights);
    }

    /** Returns density of 'features', using given functions to get covering and weight. */
    public static <T> S2DensityTree density(
        List<T> features,
        Function<T, S2CellUnion> covering,
        ToLongFunction<T> weight,
        BreadthFirstTreeBuilder builder) {
      return builder.build(fromFeatures(features, covering, weight));
    }
  }

  /** An abstract combiner of {@link S2DensityTree}s. */
  abstract static class AggregateDensity implements BreadthFirstTreeBuilder.CellWeightFunction {
    /** The decoded path for each input. */
    protected final List<DecodedPath> decodedPaths = new ArrayList<>();

    /** The number of trees being combined. */
    protected int size;

    /**
     * Combines the input trees. We stream through the cells of the input trees together in one pass
     * so we don't have the fully-decoded results in memory at the same time, but the given sequence
     * must have few enough trees that we can hold them all in memory at the same time.
     */
    public void set(Iterable<S2DensityTree> trees) {
      size = Iterables.size(trees);
      while (decodedPaths.size() < size) {
        decodedPaths.add(new DecodedPath());
      }
      int i = 0;
      for (S2DensityTree tree : trees) {
        decodedPaths.get(i++).set(tree);
      }
    }
  }

  /**
   * A combiner of {@link S2DensityTree} that respects a given max size in bytes. If the set of
   * trees remains under the max size, then the highest level cells will be the most detailed paths
   * in any of the input trees. If the set of trees exceeds the max size, then we will prefer to
   * discard higher level cells from the resulting density map, so space remains bounded while we
   * gradually lose precision.
   */
  public static class SumDensity extends AggregateDensity {
    /**
     * Sums the weights for the given cell id over all the input trees. A positive weight is
     * returned if at least one of the trees has children at the cell with the given id. A negative
     * weight is returned if the given cell id is a leaf or a child of a leaf in all the trees.
     */
    @Override
    public long applyAsLong(S2CellId id) {
      long sum = 0;
      boolean allContained = true;
      for (int i = 0; i < size; i++) {
        Cell node = decodedPaths.get(i).cell(id);
        sum += node.weight;
        allContained &= !node.hasChildren();
      }
      sum = min(MAX_WEIGHT, sum);
      return allContained ? -sum : sum;
    }
  }

  /**
   * A combiner of {@link S2DensityTree} that sums the weights of the input density trees in cells
   * that intersect the leaves of all the provided density trees. Cells that do not intersect every
   * tree at the leaves are pruned from the combined tree. Like {@link SumDensity}, respects a given
   * max size in bytes. If the set of trees remains under the max size, then the highest level cells
   * will be the most detailed paths in any of the input trees. If the set of trees exceeds the max
   * size, then we will prefer to discard higher level cells from the resulting density map, so
   * space remains bounded while we gradually lose precision.
   */
  @VisibleForTesting
  static class IntersectionDensity extends SumDensity {
    private final S2CellUnion intersectingLeaves;

    /** Constructs a combiner of input densities for the given trees. */
    public static IntersectionDensity create(Iterable<S2DensityTree> trees) {
      Iterator<S2DensityTree> iterator = trees.iterator();

      // Compute the intersection of the leaves of the given S2DensityTrees as an S2CellUnion.
      // *All* of the input trees must have leaves covering a cell for that cell to appear in
      // 'leaves'.
      S2CellUnion leaves = iterator.hasNext() ? iterator.next().getLeaves() : new S2CellUnion();
      while (iterator.hasNext()) {
        S2CellUnion temp = new S2CellUnion();
        temp.getIntersection(leaves, iterator.next().getLeaves());
        leaves = temp;
      }

      return new IntersectionDensity(leaves);
    }

    private IntersectionDensity(S2CellUnion intersectingLeaves) {
      this.intersectingLeaves = intersectingLeaves;
    }

    @VisibleForTesting
    public S2CellUnion getIntersectingLeaves() {
      return intersectingLeaves;
    }

    /**
     * Sums the weights for the given cell id over the intersection of the input trees. If any of
     * the input trees doesn't have a leaf intersecting the given cell, returns zero. Otherwise,
     * behaves the same as SumDensity.
     */
    @Override
    public long applyAsLong(S2CellId id) {
      return intersectingLeaves.intersects(id) ? super.applyAsLong(id) : 0;
    }
  }

  /**
   * A random access cell decoder that only decodes as much of the path from the root down to the
   * specified cell id passed to {@link #weight} or {@link #cell} as is different from the last
   * call. The path's current {@link #cell} is either from the tree, or has weight 0 and no children
   * if the current path is disjoint from the tree. Greatly accelerates calls to many nearby cells.
   */
  public static class DecodedPath implements CellWeightFunction {
    private S2DensityTree tree;
    private Cursor cursor;
    private final List<Cell> stack = new ArrayList<>();

    /** The id of the cell last actually decoded, or null. Not an empty cell. */
    private S2CellId lastDecoded;

    private static final Cell EMPTY_CELL;

    static {
      EMPTY_CELL = new Cell();
      EMPTY_CELL.setNoChildren();
    }

    /** Default constructor. {@link #set(S2DensityTree)} must be called before use. */
    public DecodedPath() {}

    /** Constructor that calls set(). */
    public DecodedPath(S2DensityTree tree) {
      set(tree);
    }

    /** Sets the tree to decode from. Must be called before {@link #weight} or {@link #cell}. */
    public void set(S2DensityTree tree) {
      this.tree = tree;
      this.cursor = tree.encoded.cursor();
      this.lastDecoded = null;
    }

    // Implements the CellWeightFunction API, so we can use trees as input to further weighting.
    @Override public long applyAsLong(S2CellId id) {
      Cell node = cell(id);
      return node.hasChildren() ? node.weight() : -node.weight();
    }

    /**
     * Returns the weight of the Cell at 'cell' or the weight of the leaf containing 'cell', if
     * present in the tree. Otherwise returns weight 0. The stack remembers the path to the last
     * decoded cell in order to skip decoding unchanged parents.
     */
    public long weight(S2CellId id) {
      return cell(id).weight;
    }

    /**
     * Returns the Cell at 'cell' or that is a leaf containing 'cell', if present in the tree.
     * Otherwise, 'cell' has weight 0, so an {@link Cell#isEmpty() empty cell} is returned. The
     * stack remembers the path to the last decoded cell in order to skip decoding unchanged
     * parents.
     */
    public Cell cell(S2CellId cell) {
      int decodedLevel = decodedLevel(cell);
      if (decodedLevel < 0) {
        Cell node = loadFace(cell.face());
        if (node.isEmpty()) {
          return node;
        }
        decodedLevel = 0;
      }
      return loadCell(decodedLevel, cell);
    }

    /**
     * Returns the level of the deepest already-decoded parent of the given cell id, or -1 if none
     * are decoded.
     */
    private int decodedLevel(S2CellId cell) {
      if (lastDecoded == null) {
        return -1;
      } else if (lastDecoded.contains(cell)) {
        return lastDecoded.level();
      } else {
        return lastDecoded.getCommonAncestorLevel(cell);
      }
    }

    /**
     * Returns a decoded face cell from the stack for the given face, in [0..5]. If the face cell is
     * not present in the tree, an {@link Cell#isEmpty() empty cell} is returned.
     */
    private Cell loadFace(int face) {
      Cell cell = get(0);
      if (tree.facePositions[face] < 0) {
        // The face cell isn't present in this tree, so return no weight or children.
        cell.weight = 0;
        cell.setNoChildren();
        lastDecoded = null;
      } else {
        // The face cell is present so load the weight and child positions.
        cursor.position = tree.facePositions[face];
        cell.decode(tree.encoded, cursor);
        lastDecoded = FACE_CELLS[face];
      }
      return cell;
    }

    /**
     * Returns a decoded Cell for 'cell', or if that Cell is not in the tree, the deepest ancestor
     * of 'id' in the tree, if that ancestor is a leaf. Otherwise, returns an empty Cell. Only
     * levels above 'decodedLevel' must be loaded.
     */
    private Cell loadCell(int decodedLevel, S2CellId cell) {
      Cell decodedCell = get(decodedLevel);
      for (int i = decodedLevel + 1; i <= cell.level(); i++) {
        cursor.position = decodedCell.positions[cell.childPosition(i)];
        if (cursor.position < 0) {
          // This tree doesn't contain the requested cell. Return the last available ancestor, if
          // that ancestor is a leaf, otherwise an empty cell.
          return decodedCell.hasChildren() ? EMPTY_CELL : decodedCell;
        }
        decodedCell = get(i);
        // We have a child to decode, so decode it.
        decodedCell.decode(tree.encoded, cursor);
        lastDecoded = cell.parent(i);
      }
      return decodedCell;
    }

    /** Returns the cell at the given level in the stack. */
    private Cell get(int level) {
      // Breadth-first processing should only ever increase the level 1 at a time.
      if (stack.size() == level) {
        stack.add(new Cell());
      }
      return stack.get(level);
    }
  }

  /**
   * A builder of density trees that visits a 'weigher' function in breadth-first order, starting
   * with the top level {@link S2CellId#FACE_CELLS face cells}, and visits every S2 cell where the
   * weigher result is greater than 0, the cell level is below the given max level, and the final
   * tree size is estimated to be within the max size.
   */
  @JsType
  public static class BreadthFirstTreeBuilder {
    private final int approximateSizeBytes;
    private final int maxLevel;
    private final TreeEncoder encoder;
    private long[] ranges = new long[8];
    private int rangesSize;
    private long[] nextLevelRanges = new long[8];
    private int nextLevelRangesSize;

    /**
     * @param approximateSizeBytes The max size in bytes to aim for, actual usage may be 4x higher,
     *     but averages around this size
     * @param maxLevel The max level to compute density information for
     */
    @JsIgnore
    public BreadthFirstTreeBuilder(int approximateSizeBytes, int maxLevel) {
      this(approximateSizeBytes, maxLevel, new TreeEncoder());
    }

    /** As {@link #BreadthFirstTreeBuilder(int, int)} but with a reusable encoder. */
    public BreadthFirstTreeBuilder(int approximateSizeBytes, int maxLevel, TreeEncoder encoder) {
      this.approximateSizeBytes = approximateSizeBytes;
      this.maxLevel = maxLevel;
      this.encoder = encoder;
      clear();
    }

    /** Resets this builder. */
    public void clear() {
      rangesSize = 0;
      nextLevelRangesSize = 0;
      encoder.clear();
    }

    /**
     * Builds the density tree that forms from a breadth-first visitation of the given weigher,
     * which must produce values greater than 0 for cells to have their children weighed.
     *
     * @param weigher A producer of the absolute weight of each cell, which must produce 0 if the
     *     cell is disjoint, a weight between 1 and {@link #MAX_WEIGHT} if the cell intersects, and
     *     a negation of the weight if the cell is contained.
     */
    public S2DensityTree build(CellWeightFunction weigher) {
      clear();

      // Weigh cells in breadth first order until the resulting cell/weight pairs exceed the size
      // limit or we hit the max cell level.
      int sizeEstimateBytes = 0;
      ranges[rangesSize++] = S2CellId.begin(MAX_LEVEL).id();
      ranges[rangesSize++] = S2CellId.end(MAX_LEVEL).id();
      for (int level = 0; level <= maxLevel && sizeEstimateBytes < approximateSizeBytes; level++) {
        long lastHilbertEnd = S2CellId.sentinel().id();
        for (int i = 0; i < rangesSize; i += 2) {
          S2CellId start = new S2CellId(ranges[i]);
          S2CellId end = new S2CellId(ranges[i + 1]);
          for (S2CellId cell = start.parent(level); cell.lessThan(end); cell = cell.next()) {
            // Get the weight and skip this cell unless it's larger than 0.
            long weight = weigher.applyAsLong(cell);
            if (weight == 0) {
              // Skip disjoint cells.
              continue;
            } else if (weight < 0) {
              // Get the absolute weight and skip searching children, as there are none.
              weight = -weight;
            } else {
              // Add this hilbert range to the ranges to scan at the next level.
              long hilbertStart = cell.childBegin(MAX_LEVEL).id();
              long hilbertEnd = cell.childEnd(MAX_LEVEL).id();
              if (hilbertStart == lastHilbertEnd) {
                // Extend the existing range.
                nextLevelRanges[nextLevelRangesSize - 1] = hilbertEnd;
              } else {
                // Add a new range. Note the full test is equality because the array length is even.
                if (nextLevelRangesSize == nextLevelRanges.length) {
                  nextLevelRanges = Arrays.copyOf(nextLevelRanges, nextLevelRanges.length * 2);
                }
                nextLevelRanges[nextLevelRangesSize++] = hilbertStart;
                nextLevelRanges[nextLevelRangesSize++] = hilbertEnd;
              }
              lastHilbertEnd = hilbertEnd;
            }

            // Save the weight for repacking later and estimate the size it will consume.
            Preconditions.checkArgument(weight <= MAX_WEIGHT, "Weigher exceeded max weight");
            encoder.put(cell, weight);
            sizeEstimateBytes += estimateSize(weight);
          }
        }

        // Swap the next level ranges into the ranges, and clear the next level ranges.
        long[] tempArray = ranges;
        ranges = nextLevelRanges;
        nextLevelRanges = tempArray;
        rangesSize = nextLevelRangesSize;
        nextLevelRangesSize = 0;
      }

      S2DensityTree tree = encoder.build();
      clear();
      return tree;
    }

    /**
     * A function that produces the weight that intersects the given cell in an abstract collection
     * of cells with weights, such as a shape index with a weighing function, or a density tree. The
     * returned weight is 0 if the cell is disjoint from the collection. The weight is normally the
     * sum of the weights of features that intersect the cell. However, the returned weight is
     * negated if the weight is constant throughout the cell, such as when the cell is entirely
     * contained by a polygon, or when details of the distribution of weight in the cell are not
     * available.
     */
    @JsType
    public interface CellWeightFunction extends ToLongFunction<S2CellId> {
      @Override
      long applyAsLong(S2CellId cell);
    }
  }

  /** Loads face positions. The cursor must be positioned at the start of the version string. */
  static long[] decodeHeader(Bytes encoded, Cursor cursor) throws IOException {
    // Verify the version string.
    try {
      for (byte expected : VERSION_BYTES) {
        byte actual = encoded.get(cursor.position++);
        if (expected != actual) {
          throw new IOException("Invalid version byte, expected " + expected + ", got " + actual);
        }
      }
    } catch (ArrayIndexOutOfBoundsException e) {
      throw new IOException("Insufficient or invalid input bytes: ", e);
    }

    // Read the lengths of each cube face except the last.
    long[] offsets = new long[6];
    long faceMask = encoded.readVarint64(cursor);
    // There are faceBits-1 lengths to read.
    int lengthsRemaining = Long.bitCount(faceMask) - 1;
    int offset = 0;
    for (int i = 0; i < FACE_CELLS.length; i++) {
      if ((faceMask & (1L << i)) == 0) {
        offsets[i] = -1;
      } else {
        offsets[i] = offset;
        if (lengthsRemaining-- > 0) {
          offset += checkLength(encoded.readVarint64(cursor));
        }
      }
    }

    for (int i = 0; i < FACE_CELLS.length; i++) {
      if (offsets[i] >= 0) {
        offsets[i] += cursor.position;
      }
    }

    return offsets;
  }

  /** The decoded weight and positions of encoded child cells. */
  @JsType
  public static class Cell {
    private static final int[] NO_CHILDREN = {-1, -1, -1, -1};

    /** The weight of this cell. */
    long weight;

    /** The position of each child in {@link S2CellId#childPosition} order, or -1 if child unset. */
    final int[] positions = new int[4];

    /**
     * Returns true iff this cell is empty. Empty cells are used to indicate nodes that are not
     * present in a decoded tree. Empty cells have zero weight (and no children).
     */
    public boolean isEmpty() {
      return this.weight == 0L;
    }

    /** Unsets all children. */
    public void setNoChildren() {
      System.arraycopy(NO_CHILDREN, 0, positions, 0, positions.length);
    }

    /** Returns true iff there is a child available under this node. */
    public boolean hasChildren() {
      return !Arrays.equals(positions, NO_CHILDREN);
    }

    /** Loads this cell with the weight and child positions at the given cursor position. */
    public void decode(Bytes encoded, Cursor cursor) {
      long bits = encoded.readVarint64(cursor);
      weight = bits >>> CHILD_MASK_BITS;
      int childMask = (int) bits & CHILD_MASK;

      // Read the cumulative offsets for each child in the child mask.
      int offset = 0;
      for (int i = 0; i < 4; i++) {
        if ((childMask & 1) != 0) {
          // Child is set. Set the current offset.
          positions[i] = offset;
          if (childMask > 1) {
            // There's another child, so increment the cumulative offset by its offset.
            offset += checkLength(encoded.readVarint64(cursor));
          }
        } else {
          // Child is unset.
          positions[i] = -1;
        }
        // Slide the child mask down to test the next child.
        childMask >>>= 1;
      }

      // Offsets are relative to the position right after them, so increment them by that position.
      for (int i = 0; i < 4; i++) {
        if (positions[i] >= 0) {
          positions[i] = checkLength(positions[i] + cursor.position);
        }
      }
    }

    /** Returns the weight in this node. */
    public long weight() {
      return weight;
    }

    /** For debugging or logging. */
    @Override
    public String toString() {
      return Platform.formatString(
          "Cell[weight=%d, positions=%s]", weight, Arrays.toString(positions));
    }
  }

  /** Reads a length and checks that it is sensible. This is lightweight corruption detection. */
  private static int checkLength(long length) {
    Preconditions.checkState(length > 0, "Invalid length, corrupt bytes");
    Preconditions.checkState(length <= Integer.MAX_VALUE, "Invalid length, corrupt bytes");
    return (int) length;
  }

  /**
   * Returns the estimated encoded size of 'weight', where we can precompute the weight size, but
   * must guess that the offset to the weight has to skip an average of two siblings with that
   * weight.
   */
  public static int estimateSize(long weight) {
    int weightSize = varIntSize(weight << CHILD_MASK_BITS | 0XF);
    return weightSize + 2 * varIntSize(weightSize);
  }

  /**
   * A collector of cell/weight pairs and an encoder of them. Since the face and cell encodings both
   * require varint offsets to the children in the header, we write everything backwards, and then
   * reverse the array at the end.
   */
  @JsType
  static class TreeEncoder {
    private final ReversibleBytes output = new ReversibleBytes();
    // The cell/weight pairs that have been added to the encoder, in the order they were added. Only
    // the entries from 0 to size are used.
    private WeightedCell[] weights = new WeightedCell[8];
    private int size;

    TreeEncoder() {
      for (int i = 0; i < weights.length; i++) {
        weights[i] = new WeightedCell();
      }
      clear();
    }

    /** Clears any partially encoded results. */
    public void clear() {
      size = 0;
      output.clear();
    }

    /** Inserts the given cell/weight pair into the current encoder. */
    public void put(S2CellId cell, long weight) {
      if (size == weights.length) {
        weights = Arrays.copyOf(weights, weights.length * 2);
        for (int i = size; i < weights.length; i++) {
          weights[i] = new WeightedCell();
        }
      }
      weights[size++].set(cell, weight);
    }

    /**
     * Orders by rangeMin and breaks ties by reversed cell level, so reading from the end is a
     * preorder S2Cell traversal where faces and children are reversed.
     */
    private static final Comparator<WeightedCell> ENCODER_ORDER =
        (a, b) -> {
          if (b.cell.contains(a.cell)) {
            return -1;
          } else if (a.cell.contains(b.cell)) {
            return 1;
          } else {
            return a.cell.compareTo(b.cell);
          }
        };

    /** Encodes the current set of cell/weight pairs and returns the new density tree. */
    public S2DensityTree build() {
      try {
        // Encode weighted cells in order.
        Arrays.sort(weights, 0, size, ENCODER_ORDER);
        output.clear();
        ReversedLengthsWriter lengths = new ReversedLengthsWriter(output);
        int faceMask = 0;
        while (size > 0) {
          WeightedCell weighted = weights[--size];
          assert weighted.cell.level() == 0 : "First cell of a face must be a face cell";
          int face = weighted.cell.face();
          encodeCellReverse(0, weighted);
          lengths.next();
          faceMask |= 1 << face;
        }
        lengths.write(faceMask);
        output.write(REVERSED_VERSION_BYTES);
        return S2DensityTree.decode(Bytes.fromByteArray(output.reversedCopy()));
      } catch (IOException e) {
        throw new IllegalStateException("ByteStack EOF", e);
      } finally {
        clear();
      }
    }

    /**
     * Writes the encoded cell format in reverse and then clears the cell ID so 'weights' only
     * retains the heap that is useful to the next call.
     */
    private void encodeCellReverse(int level, WeightedCell parent) {
      int childMask = 0;
      ReversedLengthsWriter lengths = new ReversedLengthsWriter(output);
      while (size > 0 && parent.cell.contains(weights[size - 1].cell)) {
        WeightedCell child = weights[--size];
        assert parent.cell.equals(child.cell.parent()) : "Must be a direct descendant";
        int childLevel = level + 1;
        childMask |= 1 << child.cell.childPosition(childLevel);
        encodeCellReverse(childLevel, child);
        lengths.next();
      }
      // TODO(eengle): use 0s to refer to the parent weight
      lengths.write((parent.weight << CHILD_MASK_BITS) | childMask);
      parent.cell = null;
    }

    /** A cell and the weight measured at that cell. */
    private static class WeightedCell {
      S2CellId cell;
      long weight;

      void set(S2CellId cell, long weight) {
        this.cell = cell;
        this.weight = weight;
      }
    }

    /**
     * Records the lengths of each in a sequence of chunks written by the caller in reverse order,
     * and then writes a mask followed by the lengths of each chunk in reverse order. This allows
     * unpredictable varints to be written at the start of each parent without having a more
     * complicated scheme than simply reversing all the bytes after encoding, or a wasteful scheme
     * like double-buffering.
     */
    @VisibleForTesting
    static class ReversedLengthsWriter {
      final ReversibleBytes output;
      final int[] lengths = new int[6];
      int size;
      int start;

      ReversedLengthsWriter(ReversibleBytes output) {
        this.output = output;
        this.start = output.size;
      }

      /** Records the size of another encoded chunk. The caller should encode chunks in reverse. */
      void next() {
        lengths[size++] = output.size - start;
        start = output.size;
      }

      /**
       * Writes the mask and lengths in normal order and then reverses them, where the last length
       * is not included because the position of whatever comes after can be inferred from a higher-
       * level data structure (most often the encoded parent cell).
       */
      void write(long mask) {
        try {
          writeVarint64(output, mask);
          for (int i = size - 1; i > 0; i--) {
            writeVarint64(output, lengths[i]);
          }
        } catch (IOException e) {
          throw new IllegalStateException("Can't be thrown by ReversibleBytes: ", e);
        }
        output.reverseFrom(start);
        start = output.size;
      }
    }

    /** An output stream we can do reversals on. */
    @VisibleForTesting
    static final class ReversibleBytes extends OutputStream {
      private byte[] bytes = new byte[256];
      private int size;

      /** Resets the output counter to 0. Retains previously-allocated memory. */
      public void clear() {
        size = 0;
      }

      /** Returns the current number of bytes in this output. */
      public int size() {
        return size;
      }

      @Override
      public void write(int b) {
        maybeResize(1);
        bytes[size++] = (byte) b;
      }

      @Override
      public void write(byte[] b) {
        write(b, 0, b.length);
      }

      @Override
      public void write(byte[] b, int off, int len) {
        maybeResize(len);
        for (int i = 0; i < len; i++) {
          bytes[size++] = b[off++];
        }
      }

      private void maybeResize(int increase) {
        while (size + increase > bytes.length) {
          bytes = Arrays.copyOf(bytes, bytes.length * 2);
        }
      }

      /** Reverses the range of bytes from start (inclusive) to the end of the output. */
      public void reverseFrom(int start) {
        for (int i = start, j = size - 1; i < j; i++, j--) {
          byte b = bytes[i];
          bytes[i] = bytes[j];
          bytes[j] = b;
        }
      }

      /** Returns a copy of the bytes. */
      public byte[] toByteArray() {
        return Arrays.copyOf(bytes, size);
      }

      /** Returns a copy of the bytes in reverse order. */
      public byte[] reversedCopy() {
        byte[] reversed = new byte[size];
        for (int i = 0, j = size - 1; i < size; i++, j--) {
          reversed[i] = bytes[j];
        }
        return reversed;
      }
    }
  }
}
