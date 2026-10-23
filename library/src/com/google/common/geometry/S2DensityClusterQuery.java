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

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkState;
import static com.google.common.collect.Comparators.max;
import static com.google.common.collect.Comparators.min;
import static com.google.common.collect.ImmutableList.of;
import static com.google.common.collect.ImmutableList.toImmutableList;
import static com.google.common.collect.Maps.immutableEntry;
import static com.google.common.collect.Streams.stream;
import static com.google.common.geometry.S2CellId.MAX_LEVEL;
import static com.google.common.geometry.S2CellId.begin;
import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.geometry.S2DensityTree.CellVisitor.Action;
import com.google.common.geometry.S2DensityTree.DecodedPath;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.function.Function;
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
  private final boolean interpolate;

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
    this(minClusterWeight, maxClusterWeight, true);
  }

  private S2DensityClusterQuery(long minClusterWeight, long maxClusterWeight, boolean interpolate) {
    Preconditions.checkArgument(minClusterWeight > 0);
    Preconditions.checkArgument(minClusterWeight <= maxClusterWeight);
    this.minClusterWeight = minClusterWeight;
    this.maxClusterWeight = maxClusterWeight;
    this.interpolate = interpolate;
  }

  /** Returns a copy of this query with interpolation set to the given value. */
  public S2DensityClusterQuery withInterpolation(boolean interpolate) {
    return new S2DensityClusterQuery(minClusterWeight, maxClusterWeight, interpolate);
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

  /**
   * {@link #edits Edits} 'clusters' in the range {@code [begin, end)} due to the new 'density',
   * creating new keys with 'newKey'. Clusters outside the range are removed. Clusters that cross
   * the range are clipped to the range.
   */
  public <K> void edit(
      Collection<Entry<K, Cluster>> clusters,
      S2DensityTree density,
      Function<Cluster, K> newKey,
      S2CellId begin,
      S2CellId end) {
    List<Entry<K, Cluster>> edits = new ArrayList<>();
    edits(edits, clusters, density, begin, end);
    clusters.clear();
    for (Entry<K, Cluster> edit : edits) {
      K key = edit.getKey();
      Cluster cluster = edit.getValue();
      if (cluster.weight() > 0) {
        if (key == null) {
          edit = immutableEntry(newKey.apply(cluster), cluster);
        }
        clusters.add(edit);
      }
    }
  }

  /**
   * Adds 'edits' that will update old 'clusters' in the range {@code [begin, end)} due to the new
   * 'density'. New clusters have a null key, removed clusters have a 0 weight, and other clusters
   * have been updated. Clusters outside the range are removed, and clusters that cross the range
   * are clipped to the range. The caller should apply these edits to the source of 'clusters' prior
   * to calling this method again. The caller may provide these edits to {@link #regionKeys} to set
   * keys of features in each 'edits' key.
   *
   * <p>The edits preserve the relative order of keys (but may delete or insert around that order).
   *
   * @param edits the collection to populate with edits
   * @param clusters the existing clusters to update, which must be disjoint and in spatial order
   * @param density the density tree to cluster
   * @param begin the start of the cell range (inclusive)
   * @param end the end of the cell range (exclusive)
   */
  public <K> void edits(
      Collection<Entry<K, Cluster>> edits,
      Collection<? extends Entry<K, Cluster>> clusters,
      S2DensityTree density,
      S2CellId begin,
      S2CellId end) {
    Entry<K, Cluster> tail = immutableEntry(null, new Cluster(begin, begin, 0));
    for (Entry<K, Cluster> entry : clusters) {
      // Delete clusters outside the range (set weight to 0).
      Cluster c = entry.getValue();
      if (c.end().lessOrEquals(begin) || c.begin().greaterOrEquals(end)) {
        edits.add(immutableEntry(entry.getKey(), c.withWeight(0)));
        continue;
      }

      // Compute the range of this interval.
      S2CellId cBegin = max(c.begin(), begin);
      S2CellId cEnd = min(c.end(), end);
      checkArgument(tail.getValue().end().lessOrEquals(cBegin), "Not disjoint and sorted.");

      // Process gap before cluster if present.
      if (tail.getValue().end().lessThan(cBegin)) {
        tail = processInterval(edits, density, null, tail.getValue().end(), cBegin, tail);
      }
      // Process cluster.
      tail = processInterval(edits, density, entry.getKey(), cBegin, cEnd, tail);
    }

    // Process trailing gap if present.
    if (tail.getValue().end().lessThan(end)) {
      tail = processInterval(edits, density, null, tail.getValue().end(), end, tail);
    }

    // Flush tail if active.
    if (tail.getKey() != null || tail.getValue().weight() > 0) {
      edits.add(immutableEntry(tail.getKey(), tail.getValue().clamp(density)));
    }
  }

  /**
   * Processes a single Hilbert interval {@code [begin, end)}, merging with any active {@code tail}
   * from preceding intervals, emitting completed clusters directly to {@code edits}, splitting
   * oversized ranges, and returning the new {@code tail}.
   */
  private <K> Entry<K, Cluster> processInterval(
      Collection<Entry<K, Cluster>> edits,
      S2DensityTree density,
      K key,
      S2CellId begin,
      S2CellId end,
      Entry<K, Cluster> tail) {
    K tailKey = tail.getKey();
    Cluster tailCluster = tail.getValue();

    // Merge with active tail if present, deleting the smaller key if both are present.
    if (tailKey != null || tailCluster.weight() > 0) {
      Cluster union = cluster(density, tailCluster.begin(), end);
      if (union.weight() <= maxClusterWeight) {
        begin = tailCluster.begin();
        if (tailKey != null && key != null) {
          if (tailCluster.weight() > cluster(density, tailCluster.end(), end).weight()) {
            edits.add(immutableEntry(key, new Cluster(begin, begin, 0)));
            key = tailKey;
          } else {
            edits.add(immutableEntry(tailKey, new Cluster(begin, begin, 0)));
          }
        } else if (tailKey != null) {
          key = tailKey;
        }
        tail = immutableEntry(null, new Cluster(begin, begin, 0));
      } else {
        Cluster head = next(density, tailCluster.begin(), end, 0);
        if (head.end().lessThan(end) && head.weight() <= maxClusterWeight) {
          Cluster rem = cluster(density, head.end(), end);
          if (rem.weight() >= minClusterWeight) {
            edits.add(immutableEntry(tailKey, head.clamp(density)));
            begin = head.end();
            tail = immutableEntry(null, new Cluster(begin, begin, 0));
          }
        }
        if (tail.getKey() != null || tail.getValue().weight() > 0) {
          edits.add(immutableEntry(tail.getKey(), tail.getValue().clamp(density)));
          tail = immutableEntry(null, new Cluster(begin, begin, 0));
        }
      }
    }

    // Process and split current interval.
    Cluster cluster = cluster(density, begin, end);
    if (cluster.weight() == 0 && key == null) {
      return immutableEntry(null, new Cluster(end, end, 0));
    }
    if (cluster.weight() > maxClusterWeight) {
      S2CellId cur = begin;
      while (cur.lessThan(end)) {
        Cluster split = next(density, cur, end);
        if (split.end().equals(end)) {
          return immutableEntry(key, split);
        }
        edits.add(immutableEntry(key, split.clamp(density)));
        key = null;
        cur = split.end();
      }
      return immutableEntry(null, new Cluster(end, end, 0));
    }
    if (cluster.weight() < minClusterWeight) {
      return immutableEntry(key, cluster);
    }
    edits.add(immutableEntry(key, cluster.clamp(density)));
    return immutableEntry(null, new Cluster(end, end, 0));
  }

  /**
   * Returns the weight of tree nodes contained by the range {@code [begin, end)}, or interpolates
   * the contained portion of tree leaves.
   */
  @VisibleForTesting
  static Cluster cluster(S2DensityTree density, S2CellId begin, S2CellId end) {
    return new Cluster(begin, end, weight(density, begin, end));
  }

  /** Returns the weight of tree nodes contained by the range {@code [begin, end)}. */
  public static long weight(S2DensityTree density, S2CellId begin, S2CellId end) {
    long[] sum = {0};
    density.visitCells((cell, node) -> {
      // Skip disjoint cells.
      if (begin.greaterThan(cell.rangeMax()) || end.lessOrEquals(cell.rangeMin())) {
        return Action.SKIP_CELL;
      }
      // Use the entire weight of contained cells.
      if (begin.lessOrEquals(cell.rangeMin()) && end.greaterThan(cell.rangeMax())) {
        sum[0] += node.weight();
        return Action.SKIP_CELL;
      }
      // Recurse if we can since children may be less ambiguous.
      if (node.hasChildren()) {
        return Action.ENTER_CELL;
      }
      // Otherwise determine the portion of 'cell' that intersects the begin,end range.
      CellInterpolator range = new CellInterpolator(cell);
      double t1 = max(0, range.uninterpolate(begin));
      double t2 = end.isValid() ? min(1, range.uninterpolate(end)) : 1;
      sum[0] += round(ceil(node.weight() * (t2 - t1)));
      return Action.SKIP_CELL;
    });
    return sum[0];
  }

  /**
   * Returns a function that maps a region to the key it should be assigned to as part of cluster
   * updates. This function must be called on all inserted and updated records, as well as unedited
   * records keyed by any key in {@code clusters}. After such records have been given keys by this
   * function, and the cluster edits in {@code clusters} have been persisted as well, the dataset
   * will be in a properly clustered state except for any exceptional results added to the default
   * key, which allows further clustering.
   *
   * @param clusters the edited clusters, for updated/removed/inserted features
   * @param defaultKey a default key to use if a record has no good spatial match
   */
  public static <K> Function<S2Region, K> regionKeys(
      Iterable<Entry<K, Cluster>> clusters, K defaultKey) {
    // Make parallel lists of keys and coverings, to select the key for the best covering by index.
    List<K> keys = new ArrayList<>();
    List<S2CellUnion> coverings = new ArrayList<>();
    for (Entry<K, Cluster> edit : clusters) {
      K key = edit.getKey();
      Cluster cluster = edit.getValue();
      if (cluster.weight() > 0) {
        keys.add(key);
        coverings.add(cluster.covering());
      }
    }
    // Use the first non-covering index as the default and put the default key in that index.
    int defaultKeyIndex = keys.size();
    keys.add(defaultKey);

    // Return a sharder that maps a region to the proper key.
    // TODO(eengle): Have S2RegionSharder return the default if a region intersects *many* clusters.
    S2RegionSharder sharder = new S2RegionSharder(coverings);
    return region -> keys.get(sharder.getMostIntersectingShard(region, defaultKeyIndex));
  }

  /** Returns the next cluster in 'density' by weighing cells from 'begin'. */
  Cluster next(S2DensityTree density, S2CellId begin) {
    return next(density, begin, S2CellId.end(MAX_LEVEL));
  }

  /** Returns the next cluster in 'density' by weighing cells from 'begin' up to 'end'. */
  Cluster next(S2DensityTree density, S2CellId begin, S2CellId end) {
    return next(density, begin, end, 0);
  }

  /** Returns the next cluster in 'density' by incrementing 'weight' from 'begin' up to 'end'. */
  Cluster next(S2DensityTree density, S2CellId begin, S2CellId end, long weight) {
    WeightVisitor visitor = new WeightVisitor(begin, end, weight);
    density.visitCells(visitor);
    return visitor.result();
  }

  /** A weigher of density tree cells. */
  private class WeightVisitor implements S2DensityTree.CellVisitor {
    private final S2CellId begin;
    private S2CellId end;
    private long weight;

    /**
     * Creates a new weight visitor.
     * @param begin the first S2CellId in the range to weigh
     * @param end the first S2CellId after the range to weigh
     * @param weight the starting weight before 'begin'
     **/
    WeightVisitor(S2CellId begin, S2CellId end, long weight) {
      this.begin = begin;
      this.end = end;
      this.weight = weight;
    }

    @Override public Action visit(S2CellId cell, S2DensityTree.Cell node) {
      if (cell.rangeMax().lessThan(begin)) {
        return Action.SKIP_CELL;
      }
      if (cell.rangeMin().greaterOrEquals(end)) {
        return Action.STOP;
      }
      CellInterpolator range = new CellInterpolator(cell);
      double ratio = node.hasChildren() ? 0 : max(0, min(1, range.uninterpolate(begin)));
      long sum = weight + (long) ceil((1 - ratio) * node.weight());
      if (sum < minClusterWeight && cell.rangeMax().lessThan(end)) {
        weight = sum;
        return Action.SKIP_CELL;
      }
      if (sum <= maxClusterWeight && cell.rangeMax().lessThan(end)) {
        weight = sum;
        end = cell.rangeMax().next();
        return Action.STOP;
      }
      if (node.hasChildren()) {
        return Action.ENTER_CELL;
      }
      if (interpolate) {
        // Current node is too large and has no children, so prorate the Hilbert range.
        long missing = (minClusterWeight + maxClusterWeight + 1) / 2 - weight;
        double targetRatio = ratio + 1.0 * missing / node.weight();
        if (targetRatio < 1.0) {
          S2CellId nextEnd = range.interpolate(targetRatio).rangeMin();
          if (nextEnd.lessThan(end)) {
            weight += missing;
            end = nextEnd;
            return Action.STOP;
          }
        }
        double t = range.uninterpolate(end);
        weight += (long) ceil((max(0.0, min(1.0, t)) - ratio) * node.weight());
      } else if (cell.lessThan(end)) {
        weight = sum;
        end = cell.rangeMax().next();
      } else {
        end = max(begin, cell.rangeMin());
      }
      return Action.STOP;
    }

    /** Returns the cluster result that we just weighed. */
    Cluster result() {
      return new Cluster(begin, end, weight);
    }
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
    /** A coder for one {@link Cluster}. */
    public static final S2Coder<Cluster> CODER = new S2Coder<>() {
      @Override public void encode(Cluster c, OutputStream out) throws IOException {
        LittleEndianOutput.writeLong(out, c.begin().id());
        EncodedInts.writeVarint64(out, c.end().distanceFromBegin() - c.begin().distanceFromBegin());
        EncodedInts.writeVarint64(out, c.weight());
      }
      @Override public Cluster decode(Bytes data, Cursor cursor) throws IOException {
        try (InputStream is = data.toInputStream(cursor)) {
          S2CellId begin = new S2CellId(LittleEndianInput.readLong(is));
          S2CellId end = begin.advance(EncodedInts.readVarint64(is));
          long weight = EncodedInts.readVarint64(is);
          return new Cluster(begin, end, weight);
        }
      }
      @Override public boolean isLazy() {
        return false;
      }
    };

    /** A coder for {@link Cluster clusters}. Compresses much better than {@link #CODER}. */
    public static final S2Coder<List<Cluster>> VECTOR_CODER = new S2Coder<>() {
      @Override public void encode(List<Cluster> clusters, OutputStream out) throws IOException {
        // Encode the begins and ends using the existing S2CellIdVectorCoder.
        S2CellIdVectorCoder.INSTANCE.encode(Lists.transform(clusters, Cluster::begin), out);
        S2CellIdVectorCoder.INSTANCE.encode(Lists.transform(clusters, Cluster::end), out);
        // Clusters were created to be nearly the same size, so encode the deviation from the mean.
        long mean = round(clusters.stream().mapToLong(Cluster::weight).average().orElse(0));
        Longs deltas = new Longs() {
          @Override public int length() {
            return clusters.size();
          }
          @Override public long get(int index) {
            return EncodedInts.encodeZigZag64(clusters.get(index).weight() - mean);
          }
        };
        EncodedInts.writeVarint64(out, mean);
        UintVectorCoder.UINT64.encode(deltas, out);
      }
      @Override public List<Cluster> decode(Bytes data, Cursor cursor) throws IOException {
        List<S2CellId> begins = S2CellIdVectorCoder.INSTANCE.decode(data, cursor);
        List<S2CellId> ends = S2CellIdVectorCoder.INSTANCE.decode(data, cursor);
        long mean = data.readVarint64(cursor);
        Longs deltas = UintVectorCoder.UINT64.decode(data, cursor);
        return new AbstractList<>() {
          @Override public int size() {
            return begins.size();
          }
          @Override public Cluster get(int index) {
            long weight = mean + EncodedInts.decodeZigZag64(deltas.get(index));
            return new Cluster(begins.get(index), ends.get(index), weight);
          }
        };
      }
    };

    private final S2CellId begin;
    private final S2CellId end;
    private final long weight;

    /** Creates a cluster with the given range of {@link S2CellId#isLeaf leaf cells} and weight. */
    public Cluster(S2CellId begin, S2CellId end, long weight) {
      Preconditions.checkArgument(begin.isLeaf() && end.isLeaf());
      Preconditions.checkArgument(begin.lessOrEquals(end));
      this.begin = begin;
      this.end = end;
      this.weight = weight;
    }

    /** Returns the level 30 {@link S2CellId#begin beginning} of this cluster range. */
    public S2CellId begin() {
      return begin;
    }

    /** Returns the level 30 {@link S2CellId#end end} of this cluster range. */
    public S2CellId end() {
      return end;
    }

    /** Returns the weight of this cluster. */
    public long weight() {
      return weight;
    }

    /** Returns the updated weight of this cluster. */
    public long weight(S2DensityTree updated) {
      return S2DensityClusterQuery.weight(updated, begin, end);
    }

    /** Returns a new cluster with the given weight. */
    public Cluster withWeight(long weight) {
      return new Cluster(begin, end, weight);
    }

    /** Returns whether the given cell intersects the range of this cluster. */
    public boolean intersects(S2CellId id) {
      return id.rangeMax().greaterOrEquals(begin) && id.rangeMin().lessThan(end);
    }

    /** Returns whether the given cluster intersects the range of this cluster. */
    public boolean intersects(Cluster other) {
      return this.end.greaterThan(other.begin) && this.begin.lessThan(other.end);
    }

    /** Returns true if this cluster range contains the other cluster range. */
    public boolean contains(Cluster other) {
      return this.begin.lessOrEquals(other.begin) && this.end.greaterOrEquals(other.end);
    }

    /** Returns true if this cluster range is contained by the other cluster range. */
    public boolean within(Cluster other) {
      return other.contains(this);
    }

    /** Returns true if this cluster range overlaps the other range but is neither within nor contains it. */
    public boolean crosses(Cluster other) {
      return intersects(other) && !within(other) && !contains(other);
    }

    @Override public int hashCode() {
      return begin.hashCode() ^ end.hashCode() ^ (int) (weight ^ 31);
    }

    @Override
    public boolean equals(Object o) {
      if (o instanceof Cluster) {
        Cluster c = (Cluster) o;
        return begin.equals(c.begin) && end.equals(c.end) && weight == c.weight;
      } else {
        return false;
      }
    }

    @Override public String toString() {
      // The level is always 30, so just show the face:pos instead of the always-long token.
      return begin.face() + ":" + Long.toHexString(begin.pos()) + ","
          + end.face() + ":" + Long.toHexString(end.pos()) + "="
          + weight;
    }

    /**
     * Returns a new cluster that is the intersection between this cluster and the non-zero leaves
     * of 'density'. The weight is the same, but the range may be far smaller.
     */
    public Cluster clamp(S2DensityTree density) {
      DecodedPath path = new DecodedPath();
      path.set(density);
      S2CellId clampedBegin = clamp(path, S2CellId.fromFace(0), 6);
      S2CellId clampedEnd = clamp(path, S2CellId.fromFace(5), -6);
      if (!clampedBegin.isValid()) {
        return new Cluster(begin, begin, 0); // No density at all
      }
      assert clampedEnd.isValid(); // clampedBegin.isValid()==clampedEnd.isValid()
      clampedBegin = max(begin, clampedBegin.rangeMin());
      clampedEnd = min(end, clampedEnd.next().rangeMin());
      checkState(clampedBegin.lessOrEquals(clampedEnd));
      return new Cluster(clampedBegin, clampedEnd, weight);
    }

    /** Returns the min S2CellId with weight in 'path' that comes at or after 'begin', or 'end'. */
    private S2CellId clamp(DecodedPath path, S2CellId id, int steps) {
      for (int dir = Integer.signum(steps); steps != 0; steps -= dir, id = id.advance(dir)) {
        // Skip IDs disjoint from the cluster range.
        if (id.rangeMax().lessThan(begin) || id.rangeMin().greaterOrEquals(end)) {
          continue;
        }
        // Skip IDs disjoint from the density tree's leaves.
        S2DensityTree.Cell cell = path.cell(id);
        if (cell.weight() == 0) {
          continue;
        }
        // Return leaves of the density tree.
        if (!cell.hasChildren()) {
          return id;
        }
        // Otherwise search the children, returning the first valid one.
        S2CellId child = clamp(path, id.child(dir < 0 ? 3 : 0), 4 * dir);
        if (child.isValid()) {
          return child;
        }
      }
      return S2CellId.sentinel();
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
  }
}
