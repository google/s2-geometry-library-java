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

import static com.google.common.collect.ImmutableList.of;
import static com.google.common.collect.ImmutableList.toImmutableList;
import static com.google.common.collect.Streams.stream;
import static com.google.common.geometry.S2CellId.MAX_LEVEL;
import static com.google.common.geometry.S2CellId.begin;
import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2DensityTree.CellVisitor.Action;
import java.util.Iterator;
import java.util.stream.Stream;

/**
 * A query that computes clusters from density that are each {@link #S2DensityClusterQuery(long)
 * near} a target cluster weight, or {@link #S2DensityClusterQuery(long, long) within a range} of
 * acceptable cluster weights. The clusters may be computed in several ways:
 *
 * <ul>
 *   <li>{@link #clusters}: Returns a {@link Cluster} for each cluster, which provides the cluster
 *       {@link Cluster#weight weight} and S2CellId range, from inclusive {@link Cluster#begin} to
 *       exclusive {@link Cluster#end}). A {@link Cluster#covering} method provides easy access to
 *       the range as a covering suitable for e.g. {@link S2CellIndex}.
 *   <li>{@link #coverings}: A convenience method that returns the {@link Cluster#covering
 *       coverings} of each cluster when the weight and range aren't directly needed. The union of
 *       these coverings spans the world, so consider:
 *   <li>{@link #tightCoverings}: A convenience method that returns the {@link Cluster#covering
 *       coverings} of each cluster intersected with the leaves of the density tree. This is the
 *       tightest covering of each cluster that can be produced from the density tree, but it can be
 *       far more complex than desirable, requiring as many cells across all clusters as the density
 *       tree has leaf cells. To set a tradeoff between tightness and complexity, consider:
 *   <li>{@link #approximateCoverings}: A convenience method that returns the tight coverings of
 *       each cluster re-covered with a coverer, which defaults to {@link #DEFAULT_COVERER} to
 *       produce 256-cell coverings of each cluster. This essentially re-covers the tight coverings,
 *       but faster and with much less RAM than would be needed to generate the tight coverings only
 *       to approximate them. However, these coverings may overlap each other slightly, and since
 *       it's difficult to configure a coverer to generate non-overlapping coverings, consider:
 *   <li>{@link #clusterRegions}: A convenience method that returns the clusters as a region, which
 *       can be used to test whether any cell intersects the cluster without having to generate a
 *       covering for it. Regions are tight like {@link #tightCoverings}, lightweight like {@link
 *       #approximateCoverings}, but do not intersect each other (no S2Cell is contained by more
 *       than one of the resulting regions).
 * </ul>
 *
 * <p>In each case above, if the provided density tree has leaves that are larger than the maximum
 * cluster weight, this class will linearly interpolate clusters within the density tree leaf as
 * though the density is linearly distributed throughout the Hilbert range of the density leaf cell.
 * Thus it isn't necessary to ensure all leaf cells are smaller than the desired cluster weight, but
 * rather to ensure the density tree is detailed enough that the distribution has become roughly
 * linear in the leaf cells (this is not always possible since geospatial datasets are so skewed,
 * but each additional level of density tree cells significantly resolves the skew, so the
 * distribution becomes linear more quickly than might be expected).
 */
public class S2DensityClusterQuery {
  /** A good default coverer for complex clusters. */
  public static final S2RegionCoverer DEFAULT_COVERER =
      S2RegionCoverer.builder().setMaxCells(256).build();

  private final long minClusterWeight;
  private final long maxClusterWeight;

  /** Creates a new clusterer with clusters +/- 20% of the given desired cluster weight. */
  public S2DensityClusterQuery(long desiredClusterWeight) {
    this(max(1, desiredClusterWeight * 8 / 10), max(1, desiredClusterWeight * 12 / 10));
  }

  /**
   * Creates a new clusterer with finer control than {@link #S2DensityClusterQuery(long)}.
   *
   * @param minClusterWeight the smallest desired cluster weight; aim for 80% of maxClusterWeight
   * @param maxClusterWeight the largest desired cluster weight; aim for 125% of minClusterWeight
   */
  public S2DensityClusterQuery(long minClusterWeight, long maxClusterWeight) {
    Preconditions.checkArgument(minClusterWeight > 0);
    Preconditions.checkArgument(minClusterWeight <= maxClusterWeight);
    this.minClusterWeight = minClusterWeight;
    this.maxClusterWeight = maxClusterWeight;
  }

  /** Returns the min weight of resulting clusters. */
  public long minClusterWeight() {
    return minClusterWeight;
  }

  /** Returns the max weight of resulting clusters. */
  public long maxClusterWeight() {
    return maxClusterWeight;
  }

  /**
   * Returns coverings of {@link #clusters} built from 'density'. This simply collects the {@link
   * Cluster#covering coverings} without any further processing.
   */
  public ImmutableList<S2CellUnion> coverings(S2DensityTree density) {
    return stream(clusters(density)).map(Cluster::covering).collect(toImmutableList());
  }

  /**
   * Returns the {@link Cluster#intersection} of each cluster covering with the {@link
   * S2DensityTree#getLeaves leaves} of 'density', which provides the tightest covering that can be
   * produced for each cluster. This results in the fewest false positives in an {@link S2CellIndex}
   * or similar lookup system, but these coverings are usually far larger than desirable and can (in
   * aggregate) have as many cells as the density tree has leaves.
   */
  public ImmutableList<S2CellUnion> tightCoverings(S2DensityTree density) {
    S2CellUnion leaves = density.getLeaves();
    leaves.normalize();
    return stream(clusters(density)).map(c -> c.intersection(leaves)).collect(toImmutableList());
  }

  /**
   * As {@link #approximateCoverings(S2RegionCoverer, S2DensityTree)} with {@link #DEFAULT_COVERER}
   * as the coverer.
   */
  public ImmutableList<S2CellUnion> approximateCoverings(S2DensityTree density) {
    return approximateCoverings(DEFAULT_COVERER, density);
  }

  /**
   * Returns the coverings from 'coverer' for the {@link S2RegionIntersection} between each {@link
   * Cluster#covering cluster covering} and the {@link S2DensityTree#getLeaves leaves} of 'density'.
   * This results in coverings that are equivalent to re-covering {@link #tightCoverings} to make
   * them simpler, but is far more efficient. The resulting coverings are approximate and so may
   * overlap each other, but avoid being either too complex or disjoint from the density tree up to
   * the level the coverer is able to produce. Clusters with around a million features are usually
   * pretty well described with a covering of 256 cells.
   */
  public ImmutableList<S2CellUnion> approximateCoverings(
      S2RegionCoverer coverer, S2DensityTree density) {
    return regions(density).map(coverer::getCovering).collect(toImmutableList());
  }

  /**
   * Returns the clusters from 'density' as regions that intersect the range of each cluster and the
   * leaves of 'density'. The resulting regions can be {@link S2RegionCoverer#getCovering covered}
   * to produce coverings like {@link #approximateCoverings}, but it's usually more useful to work
   * with cluster boundaries as regions because they can be tested at any cell level.
   */
  public ImmutableList<S2Region> clusterRegions(S2DensityTree density) {
    return regions(density).collect(toImmutableList());
  }

  /** Returns a stream of regions that intersect the density tree leaves and the cluster. */
  private Stream<S2Region> regions(S2DensityTree density) {
    S2Region region = density.asRegion();
    return stream(clusters(density)).map(c -> new S2RegionIntersection(of(c.covering(), region)));
  }

  /**
   * Returns {@link Cluster clusters} computed from 'density'. The resulting clusters collectively
   * cover the density tree's spatial extent, i.e. every node of the density tree is contained by
   * one of the S2CellUnions. The resulting clusters will not overlap each other. This method builds
   * clusters that each have weight approximately equal to the desired cluster weight.
   */
  public Iterable<Cluster> clusters(S2DensityTree density) {
    S2DensityTree tree = density.normalize();
    return () ->
        new Iterator<>() {
          Cluster c = S2DensityClusterQuery.this.next(tree, begin(MAX_LEVEL));

          @Override
          public boolean hasNext() {
            return c.weight > 0;
          }

          @Override
          public Cluster next() {
            Cluster result = c;
            c = S2DensityClusterQuery.this.next(tree, c.end);
            return result;
          }
        };
  }

  /** Returns the next cluster beginning from 'start'. */
  private Cluster next(S2DensityTree tree, S2CellId start) {
    Cluster cluster = new Cluster(start);
    tree.visitCells(
        (cell, node) -> {
          if (cell.rangeMax().lessThan(cluster.begin)) {
            // Current node is entirely before this cluster starts, so move forward.
            return Action.SKIP_CELL;
          }
          CellInterpolator range = new CellInterpolator(cell);
          double ratio =
              node.hasChildren() ? 0 : max(0, min(1, range.uninterpolate(cluster.begin)));
          long sum = cluster.weight + (long) ceil((1 - ratio) * node.weight());
          if (sum < minClusterWeight) {
            // Current node is too small to finish the current cluster, so move forward.
            cluster.weight = sum;
            return Action.SKIP_CELL;
          }
          if (sum <= maxClusterWeight) {
            // Current node is large enough to finish the current cluster, so stop.
            cluster.weight = sum;
            cluster.end = cell.rangeMax().next();
            return Action.STOP;
          }
          if (node.hasChildren()) {
            // Current node is too large but has children, so explore the children.
            return Action.ENTER_CELL;
          }
          // Current node is too large and has no children, so prorate the Hilbert range.
          long missingWeight = (minClusterWeight + maxClusterWeight + 1) / 2 - cluster.weight;
          cluster.weight += missingWeight;
          ratio += 1.0 * missingWeight / node.weight();
          cluster.end = range.interpolate(ratio).rangeMin();
          return Action.STOP;
        });
    return cluster;
  }

  /** An interpolator of cells. */
  @VisibleForTesting
  static class CellInterpolator {
    private final S2CellId parent;
    private final int level;
    private final long begin;
    private final long length;

    CellInterpolator(S2CellId parent) {
      // A cell's Hilbert range can be larger than 2^52 at level 30, which floating point math can't
      // handle without underflow. So choose a cell level up to 26 levels deeper, so the range can
      // be no larger than 2^(26*2) bits.
      this.parent = parent;
      this.level = min(MAX_LEVEL, parent.level() + 26);
      this.begin = parent.childBegin(level).distanceFromBegin();
      this.length = parent.childEnd(level).distanceFromBegin() - begin;
    }

    /** Returns a cell that is 'ratio' percent of the way along the Hilbert range of 'parent'. */
    public S2CellId interpolate(double ratio) {
      long steps = (long) ceil(ratio * length);
      return parent.childBegin(level).advance(steps);
    }

    /** Returns fraction along 'parent' that 'child' begins. */
    public double uninterpolate(S2CellId child) {
      Preconditions.checkArgument(child.level() >= parent.level());
      S2CellId adjusted = level <= child.level() ? child.parent(level) : child.childBegin(level);
      return 1.0 * (adjusted.distanceFromBegin() - begin) / length;
    }
  }

  /** A cluster is a cell range and an expected weight based on the density it was computed from. */
  public static final class Cluster {
    /** The cluster boundary as an S2CellId Hilbert range. */
    private S2CellId begin;

    private S2CellId end = S2CellId.end(MAX_LEVEL);

    /** The cluster weight. */
    private long weight;

    /** Initialize a cluster with weight 0, from 'begin' to end of the Hilbert range. */
    private Cluster(S2CellId begin) {
      Preconditions.checkArgument(begin.level() == MAX_LEVEL);
      this.begin = begin;
    }

    /** Returns the level 30 {@link S2CellId#begin beginning} of this cluster range. */
    public S2CellId begin() {
      return begin;
    }

    /** Returns the level 30 {@link S2CellId#end end} of this cluster range. */
    public S2CellId end() {
      return end;
    }

    /** Returns the covering of this cluster. */
    public S2CellUnion covering() {
      S2CellUnion covering = new S2CellUnion();
      covering.initFromBeginEnd(begin, end);
      return covering;
    }

    /** Returns the intersection of this cluster's range with the given covering. */
    public S2CellUnion intersection(S2CellUnion covering) {
      S2CellUnion intersection = new S2CellUnion();
      intersection.getIntersection(covering, covering());
      return intersection;
    }

    /** Returns the weight of this cluster. */
    public long weight() {
      return weight;
    }
  }
}
