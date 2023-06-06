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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.geometry.S2DensityTree.Cell;
import com.google.common.geometry.S2DensityTree.DecodedPath;
import com.google.common.geometry.S2DensityTree.TreeEncoder;
import java.util.ArrayList;
import java.util.List;

/**
 * A query that computes clusters from density that are within 40% of the target cluster size, and
 * within 20% of each of the upper and lower target cluster boundaries. The clusters computed are
 * as detailed as the provided S2DensityTree, and even more detailed (via interpolation) if needed
 * to achieve the target cluster sizes. As such, the clusters may be large, generally containing
 * (in aggregate) as many S2CellIds as the provided S2DensityTree has leaf nodes. However, the
 * clusters (S2CellUnions) are normalized, so groups of 4 child cells will be replaced by their
 * parent cell.
 *
 * <p>The clusters are built from leaf nodes of the density tree, so only the parts of the Hilbert
 * curve range that actually have weight in the S2DensityTree are covered by clusters.
 *
 * <p>If and when the density tree does not have enough detail to build clusters sufficiently close
 * to the target weight, interpolation within S2Cells in [rangeMin() ... rangeMax()] of density tree
 * leaf nodes is used to divide weight between clusters.
 */
public class S2DensityClusterQuery {
  // The desired weight of each cluster.
  private final long clusterSize;
  // While searching for a cluster boundary, the sum of nodes included so far.
  private long sum;
  // How much tolerance we have for cluster boundaries vs. the target.
  private long tolerance;

  /** The face cell ids in decreasing order. */
  private static final S2CellId[] REVERSE_FACE_CELLS = new S2CellId[6];
  static {
    for (int face = 0; face < 6; face++) {
      REVERSE_FACE_CELLS[5 - face] = S2CellId.fromFace(face);
    }
  }

  /**
   * Creates a clusterer for a given clusterSize in the same unit of weight used to compute density.
   */
  public S2DensityClusterQuery(long clusterSize) {
    this.clusterSize = clusterSize;
    this.tolerance = Math.round(0.2 * clusterSize);
  }

  /**
   * Density trees map cells to weights that intersect that cell, but larger features may intersect
   * many density cells, so the deeper into the tree we look, the more a feature will contribute
   * undetected duplicates of its weight many times.
   *
   * <p>So, first we build a normalized version of the tree where every node's weight is scaled by
   * (its parent's weight / the sum of weights of the node and its siblings). This makes the weight
   * of a parent equal to the sum of its children.
   */
  @VisibleForTesting
  static S2DensityTree normalizedTree(S2DensityTree tree) {
    DecodedPath path = new DecodedPath(tree);
    TreeEncoder encoder = new TreeEncoder();
    // Visit cells in depth-first order.
    tree.visitAll((cell, weight) -> {
      if (!cell.isFace()) {
        S2CellId parent = cell.parent();
        long parentWeight = path.weight(parent);
        long siblingWeight = 0;
        for (int i = 0; i < 4; i++) {
          siblingWeight += path.weight(parent.child(i));
        }
        weight = Math.round(1.0 * parentWeight * weight / siblingWeight);
      }
      encoder.put(cell, weight);
    });
    return encoder.build();
  }

  /**
   * Returns the leftmost leaf node's id from the part of the density tree that intersects 'ids'.
   * The given 'ids' must be in increasing order. Returns sentinel() if the intersection of the
   * tree with 'ids' is empty.
   */
  private static S2CellId firstLeafInTree(DecodedPath path, S2CellId... ids) {
    for (S2CellId id : ids) {
      Cell node = path.cell(id);
      if (node.isEmpty()) {
        continue;
      }
      if (!node.hasChildren()) {
        return id;
      }
      return firstLeafInTree(path, id.child(0), id.child(1), id.child(2), id.child(3));
    }
    return S2CellId.sentinel();
  }

  /** Returns the leftmost leaf node's cell id, or sentinel() if the tree is empty. */
  @VisibleForTesting
  static S2CellId firstLeafInTree(DecodedPath path) {
    return firstLeafInTree(path, S2CellId.FACE_CELLS);
  }

  /**
   * Returns the rightmost leaf node's id from the part of the density tree that intersects 'ids'.
   * The given 'ids' must be in decreasing order. Returns sentinel() if the intersection of the tree
   * with 'ids' is empty.
   */
  private static S2CellId lastLeafInTree(DecodedPath path, S2CellId... ids) {
    for (S2CellId id : ids) {
      Cell node = path.cell(id);
      if (node.isEmpty()) {
        continue;
      }
      if (!node.hasChildren()) {
        return id;
      }
      return lastLeafInTree(path, id.child(3), id.child(2), id.child(1), id.child(0));
    }
    return S2CellId.sentinel();
  }

  /** Returns the rightmost leaf node's cell id, or sentinel() if the tree is empty. */
  @VisibleForTesting
  static S2CellId lastLeafInTree(DecodedPath path) {
    return lastLeafInTree(path, REVERSE_FACE_CELLS);
  }

  /**
   * Returns the leftmost leaf cell that is entirely right of 'prev', in the part of the density
   * tree that intersects 'ids'. Returns sentinel() if no part of the density tree under 'ids' is
   * right of 'prev'.
   *
   * <p>The given 'ids' must be in increasing order.
   */
  private static S2CellId nextLeafInTree(DecodedPath path, S2CellId prev, S2CellId... ids) {
    for (S2CellId id : ids) {
      if (id.rangeMax().lessThan(prev.rangeMin())) {
        continue;  // id and its descendants are all completely left of 'prev'.
      }
      // id is at least partly right of 'prev', so id or some descendant of id might be the first
      // leaf right of 'prev'.
      Cell node = path.cell(id);
      if (node.isEmpty()) {
        continue;
      }
      if (node.hasChildren()) {
        // TODO(torrey): Search the tree directly, looking at node.children instead of intersecting
        // the S2CellId hierarchy with the density tree as this implementation does.
        // One of id's descendants may be right of 'prev', otherwise we need to try the next id.
        S2CellId n = nextLeafInTree(path, prev, id.child(0), id.child(1), id.child(2), id.child(3));
        if (n.isValid()) {
          return n;
        }
      } else {
        // This is a leaf node.
        if (id.rangeMin().greaterThan(prev.rangeMax())) {
          // This leaf node is entirely right of 'prev'. It is the leftmost.
          return id;
        }
      }
    }
    // All ids were tested and none had a leaf node descendent right of 'prev'.
    return S2CellId.sentinel();
  }

  /**
   * Returns the leftmost leaf cell in the density tree that is right of 'prev', or sentinel() if
   * there is none.
   */
  @VisibleForTesting
  static S2CellId nextLeafInTree(DecodedPath path, S2CellId prev) {
    return nextLeafInTree(path, prev, S2CellId.FACE_CELLS);
  }

  /** Recursively add all the leaves in the density tree under 'ids' to 'cellUnion'. */
  @VisibleForTesting
  static void collectLeavesInTree(
      DecodedPath path, S2CellUnion cellUnion, S2CellId... ids) {
    for (S2CellId id : ids) {
      Cell node = path.cell(id);
      if (node.isEmpty()) {
        continue;
      }
      if (node.hasChildren()) {
        collectLeavesInTree(path, cellUnion, id.child(0), id.child(1), id.child(2), id.child(3));
      } else {
        cellUnion.cellIds().add(id);
      }
    }
  }

  /**
   * Returns the clusters for the given tree.
   *
   * <p>Clusters are a collection of disjoint S2CellUnions that together cover the density tree's
   * spatial extent, i.e. every node of the density tree is contained by one of the S2CellUnions.
   * This method attempts to build clusters that each have weight approximately equal to
   * 'clusterSize'.
   */
  public List<S2CellUnion> clusters(S2DensityTree tree) {
    DecodedPath path = new DecodedPath(normalizedTree(tree));

    // Normalizing weights removes much of the redundancy in the tree, but some remains concentrated
    // in the higher cell levels. So knowing that cluster N should roughly cover the byte range
    // [N*clusterSize, (N+1)*clusterSize), we find cluster boundaries by searching for a density
    // node at a multiple of the clusterSize, stopping as early in the tree as possible so that we
    // don't encounter too many duplicates of a feature's weight.
    List<S2CellUnion> clusters = new ArrayList<>();

    // As clusters are formed, 'first' is updated to always be the leftmost cell id in the current
    // cluster. It is possible that 'first' is not a node in the tree: it may be interpolated within
    // the range of a leaf node of a tree. 'first' begins as the leftmost leaf of the density tree.
    S2CellId first = firstLeafInTree(path);

    // The running total of weight of the tree included in all clusters so far.
    sum = 0;
    // The target weight for the cluster currently being formed, as a sum from the beginning of the
    // tree including all previous clusters.
    long target = clusterSize;

    do {
      // TODO(eengle): Simplify this class with use of S2DensityTree.asRegion.
      // Create and grow another cluster, obtaining the start of the next cluster.
      S2CellUnion cluster = new S2CellUnion();
      first = growCluster(cluster, target, path, first, S2CellId.FACE_CELLS);
      cluster.normalize();
      clusters.add(cluster);
      target += clusterSize;
    } while (first.lessThan(S2CellId.sentinel()));

    return clusters;
  }

  /**
   * Adds cells to the given 'cluster' from the part of the tree intersecting 'ids' and beginning
   * with 'start' or a descendent of 'start'.
   *
   * Returns a S2CellId 'next', which will be the first cell of the next cluster, or
   * {@link S2CellId#sentinel()} if the cluster was grown to the end of the density tree under
   * 'ids'. In a recursive call, sentinel() indicates that the search should continue at the next
   * level up, but at the top level, sentinel() indicates that this was the last cluster.
   *
   * <p>If the current cluster ends with a node in the tree and it is not the last leaf in the tree,
   * 'next' will be the leftmost node in the tree that is right of the end of the cluster. If the
   * current cluster ends at a cell interpolated within the rangeMin() ... rangeMax() of a cell in
   * the tree, 'next' will be the next() cell in that range.
   *
   * <p>The given ids must be in increasing order.
   */
  public S2CellId growCluster(
      S2CellUnion cluster, long target, DecodedPath path, S2CellId start, S2CellId... ids) {
    // Add the weight of each 'id' from 'start' to 'sum' until (target +/- tolerance) is reached.
    for (S2CellId id : ids) {
      // If 'id' ends before start, it was included in previous clusters or doesn't intersect.
      if (id.rangeMax().lessThan(start.rangeMin())) {
        continue;
      }

      // If 'id' isn't in the tree, it has no weight, skip it.
      Cell node = path.cell(id);
      if (node.isEmpty()) {
        continue;
      }

      // At this point, part or all of the S2CellId range covered by 'id' will be part of the
      // cluster. If 'id' is not a parent of 'start', the entire node for 'id' is a candidate to be
      // included. Otherwise, we must descend into the tree under 'id', as the previous cluster did.
      if (id.equals(start) || !id.contains(start)) {
        // Consider adding the entire 'node'.
        if (node.weight() + sum <= target - tolerance) {
          // Including 'node' and its children (if any) will still be well under the desired target.
          // Add it and keep going to the next id.
          // Rather than adding 'id' itself, which could include chunks of Hilbert curve space that
          // have no weight in the tree, we collect the leaves in the tree under node.
          collectLeavesInTree(path, cluster, id);
          sum += node.weight();
          continue;
        } else if (node.weight() + sum <= target + tolerance) {
          // Including node, or the descendants of node (if any) gets us to within tolerance of the
          // target, which is close enough.
          sum += node.weight();
          collectLeavesInTree(path, cluster, id);
          // The next cluster, if any, will start at the next leaf after id, if there is one.
          return nextLeafInTree(path, id);
        }
        // Including 'node' would produce a total weight more than 120% of the target.
      }

      // The node at 'id' must be subdivided, either because it is too heavy, or because the last
      // cluster ended part way through the descendants of 'id'.
      if (node.hasChildren()) {
        // Recursively explore the children of 'node' to grow the cluster.
        S2CellId next = growCluster(
            cluster, target, path, start, id.child(0), id.child(1), id.child(2), id.child(3));
        if (!next.equals(S2CellId.sentinel())) {
          // The target weight was reached somewhere before the last descendant of 'id', we're done.
          return next;
        }
        // sentinel() indicates that the cluster has grown to the end of the tree under 'id'.
        // There are two cases:
        if (sum < target - tolerance) {
          // Case 1: We haven't reached the target weight yet. Keep going with the next id.
          continue;
        } else {
          // Case 2: The target weight was reached simultaneously with adding the last leaf under
          // 'id'. So the next cluster, if any, starts at the next leaf after 'id'.
          return nextLeafInTree(path, id);
        }

      } else {
        // The node must be subdivided but has no children. We will interpolate in the Hilbert curve
        // range of 'id', assuming the node weight is uniformly distributed across that range.
        long rangeMin = id.rangeMin().id();
        long rangeMax = id.rangeMax().id();
        long rangeSize = rangeMax - rangeMin;

        // TODO(torrey): document how this interpolation works with signed longs, as it's confusing.

        // TODO(torrey): Rather than going all the way down to level 30 leaf cells, consider how
        // much weight is needed to add to 'sum' to get to 'target', and the total weight available
        // in 'node'. Divide just finely enough to get to 'target' within 'tolerance'. Going down
        // four levels gives gives increments of (weight/128) which would likely be enough.

        // Case 1. This node has not been subdivided yet. 'start' equals id, or is not contained by
        // it. Begin subdividing the range, taking a chunk of the target size.
        if (id.equals(start) || !id.contains(start)) {
          // The next cluster will continue in this interpolated range.
          return addSubrangeToCluster(cluster, target, rangeMin, node.weight(), rangeSize);
        }

        // Case 2. This node was subdivided already and part of it was in the last cluster. We're
        // continuing from 'start', which is a leaf cell in (rangeMin .. rangeMax). We will take
        // another chunk of the range, beginning with 'start'.
        Preconditions.checkState(id.contains(start) && start.isLeaf());

        // How much of this nodes' weight is left?
        long rangeRemains = rangeMax - start.prev().id();
        double remainingWeight = 1.0 * node.weight() * rangeRemains / rangeSize;

        if (sum + remainingWeight > target + tolerance) {
          // There's too much weight remaining for this cluster. Take another chunk of the range,
          // and continue the range with the next cluster.
          return addSubrangeToCluster(cluster, target, start.id(), node.weight(), rangeSize);
        } else {
          // Add all the remaining part of the range to the current cluster.
          sum += (long) remainingWeight;

          S2CellUnion range = new S2CellUnion();
          range.initFromMinMax(start, new S2CellId(rangeMax));
          cluster.cellIds().addAll(range.cellIds());

          if (sum >= target - tolerance) {
            // That was enough weight for this cluster. The next cluster starts at the next leaf.
            return nextLeafInTree(path, new S2CellId(rangeMax));
          }
          // Otherwise, there was not enough weight left in the node of this id. Keep going with
          // the next id.
        }
      }
    }

    // We ran out of density tree intersecting 'ids' before reaching the target weight. Keep going
    // at the next level up, if there is a next level up.
    return S2CellId.sentinel();
  }

  // Adds cells from the leaf cell range from baseId to targetId to the cluster, adjusting targetId
  // if needed. Returns the next leaf cell in the range.
  private S2CellId addSubrangeToCluster(
      S2CellUnion cluster, long target, long baseId, long nodeWeight, long rangeSize) {
    double targetFraction = (1.0 * target - sum) / nodeWeight;
    // TODO(torrey): If I don't use leaf cells, I should be more accurate here.
    sum += (long) (targetFraction * nodeWeight);
    Preconditions.checkState(target - tolerance <= sum && sum <= target + tolerance);

    long targetId = (long) Math.floor(baseId + targetFraction * rangeSize);
    S2CellId interpolated = new S2CellId(targetId);

    // Leaf cell ids are two units apart, so 'interpolated' has only a 50% chance of being a
    // leaf. If it isn't a valid leaf, the preceding long value will be.
    if (!(interpolated.isValid() && interpolated.isLeaf())) {
      interpolated = new S2CellId(targetId - 1);
      Preconditions.checkState(interpolated.isValid());
      Preconditions.checkState(interpolated.isLeaf());
    }

    // Add a minimal set of S2CellIds covering [rangeMin .. interpolated] to the cluster.
    S2CellUnion range = new S2CellUnion();
    range.initFromMinMax(new S2CellId(baseId), interpolated);
    cluster.cellIds().addAll(range.cellIds());

    // The next cluster starts on the leaf cell id after 'interpolated'.
    return interpolated.next();
  }
}
