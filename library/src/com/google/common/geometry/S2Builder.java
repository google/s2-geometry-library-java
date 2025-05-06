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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2RobustCrossProd.robustCrossProd;
import static java.lang.Math.acos;
import static java.lang.Math.max;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import com.google.common.base.Preconditions;
import com.google.common.base.Strings;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderUtil.GraphShape;
import com.google.common.geometry.S2CrossingEdgesQuery.CrossingType;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.IntComparator;
import com.google.common.geometry.primitives.Ints.IntList;
import com.google.common.geometry.primitives.Ints.OfInt;
import com.google.common.geometry.primitives.Sorter;
import com.google.common.geometry.primitives.Sorter.SortableCollection;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

/**
 * S2Builder is a tool for assembling geometry from edges. Here are some of the things it is
 * designed for:
 *
 * <ol>
 *   <li>Building polygons, polylines, and polygon meshes from unsorted collections of edges.
 *   <li>Snapping geometry to discrete representations (such as S2CellId centers or E7 lat/lng
 *       coordinates) while preserving the input topology and with guaranteed error bounds.
 *   <li>Simplifying geometry (e.g. for indexing, display, or storage).
 *   <li>Importing geometry from other formats, including repairing geometry that has errors.
 *   <li>As a tool for implementing more complex operations such as polygon intersections and
 *       unions.
 * </ol>
 *
 * <p>The implementation is based on the framework of "snap rounding". Unlike most snap rounding
 * implementations, S2Builder defines edges as geodesics on the sphere (straight lines) and uses the
 * topology of the sphere (i.e., there are no "seams" at the poles or 180th meridian). The algorithm
 * is designed to be 100% robust for arbitrary input geometry. It offers the following properties:
 *
 * <ul>
 *   <li>Guaranteed bounds on how far input vertices and edges can move during the snapping process
 *       (i.e., at most the given "snapRadius").
 *   <li>Guaranteed minimum separation between edges and vertices other than their endpoints
 *       (similar to the goals of Iterated Snap Rounding). In other words, edges that do not
 *       intersect in the output are guaranteed to have a minimum separation between them.
 *   <li>Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the input already meets
 *       the output criteria then it will not be modified.
 *   <li>Preservation of the input topology (up to the creation of degeneracies). This means that
 *       there exists a continuous deformation from the input to the output such that no vertex
 *       crosses an edge. In other words, self-intersections won't be created, loops won't change
 *       orientation, etc.
 *   <li>The ability to snap to arbitrary discrete point sets (such as S2CellId centers, E7 lat/lng
 *       points on the sphere, or simply a subset of the input vertices), rather than being limited
 *       to an integer grid.
 * </ul>
 *
 * <p>Here are some of its other features:
 *
 * <ul>
 *   <li>It can handle both directed and undirected edges. Undirected edges can be useful for
 *       importing data from other formats, e.g. where loops have unspecified orientations.
 *   <li>It can eliminate self-intersections by finding all edge pairs that cross and adding a new
 *       vertex at each intersection point.
 *   <li>It can simplify polygons to within a specified tolerance. For example, if two vertices are
 *       close enough they will be merged, and if an edge passes nearby a vertex then it will be
 *       rerouted through that vertex. Optionally, it can also detect nearly straight chains of
 *       short edges and replace them with a single long edge, while maintaining the same accuracy,
 *       separation, and topology guarantees ("simplifyEdgeChains").
 *   <li>It supports many different output types through the concept of "layers" (polylines,
 *       polygons, polygon meshes, etc.) You can build multiple layers at once in order to ensure
 *       that snapping does not create intersections between different objects (for example, you can
 *       simplify a set of contour lines without the risk of having them cross each other).
 *   <li>It supports edge labels, which allow you to attach arbitrary information to edges and have
 *       it preserved during the snapping process. (This can also be achieved using layers, at a
 *       coarser level of granularity.)
 * </ul>
 *
 * <p>Every edge can have a set of non-negative integer labels attached to it, assigned when the
 * edge (or shape) is added to the builder. When used with an appropriate layer type, you can then
 * retrieve the labels associated with each output edge. This can be useful when merging or
 * combining data from several sources. (Note that in many cases it is easier to use separate output
 * layers rather than labels.)
 *
 * <p>Labels are 32-bit non-negative integers. To support other label types, you can use
 * ValueLexicon to store the set of unique labels seen so far:
 *
 * {@snippet :
 *   ValueLexicon<MyLabelType> myLabelLexicon = new ValueLexicon<>();
 *   builder.setLabel(myLabelLexicon.add(label));
 * }
 *
 * <p>The current set of labels is represented as a stack. This makes it easy to add and remove
 * labels hierarchically (e.g., polygon 5, loop 2). Use setLabel() and clearLabels() if you need at
 * most one label per edge.
 *
 * <p>Caveat: Because S2Builder only works with edges, it cannot distinguish between the empty and
 * full polygons. If your application can generate both the empty and full polygons, you must
 * implement logic outside of this class.
 *
 * <p>Example showing how to snap a polygon to E7 coordinates:
 *
 * {@snippet :
 * S2Polygon inputPolygon = ...;
 * S2Error error = new S2Error();
 *
 * // Create a builder with a SnapFunction for E7, and a polygon layer.
 * S2Builder builder = new S2Builder.Builder(new IntLatLngSnapFunction(7)).build();
 * S2PolygonLayer layer = new S2BuilderUtil.S2PolygonLayer();
 * builder.startLayer(layer);
 *
 * // Add the input polygon to the builder and build. Check that the build was successful.
 * builder.addPolygon(inputPolygon);
 * if (!builder.build(error)) {
 *   logger.atError().log(error.text());
 *    ...
 *  }
 * }</pre>
 */
@SuppressWarnings("Assertion")
public class S2Builder {

  /** The options for this S2Builder. */
  private final Builder options;

  /**
   * The maximum distance (inclusive) that a vertex can move when snapped, equal to
   * S1ChordAngle(options.snapFunction().snapRadius()).
   */
  private final S1ChordAngle siteSnapRadiusChordAngle;

  /**
   * The maximum distance (inclusive) that an edge can move when snapping to a snap site. It can be
   * slightly larger than the site snap radius when edges are being split at crossings.
   */
  private final S1ChordAngle edgeSnapRadiusChordAngle;

  /**
   * True if we need to check that snapping has not changed the input topology around any vertex
   * (i.e. Voronoi site). Normally this is only necessary for forced vertices, but if the snap
   * radius is very small (e.g., zero) and splitCrossingEdges() is true then we need to do this for
   * all vertices. In all other situations, any snapped edge that crosses a vertex will also be
   * closer than minEdgeVertexSeparation() to that vertex, which will cause us to add a separation
   * site anyway.
   */
  private final boolean checkAllSiteCrossings;

  /**
   * The maximum distance that a vertex can be separated from an edge while still affecting how
   * that edge is snapped.
   */
  private final S1Angle maxEdgeDeviation;

  /**
   * Sites within a distance of edgeSiteQueryRadiusChordAngle from an edge are candidates for
   * snapping and/or avoidance.
   */
  private final S1ChordAngle edgeSiteQueryRadiusChordAngle;

  /**
   * The maximum edge length such that even if both endpoints move by the maximum distance allowed
   * (i.e. edgeSnapRadius), the center of the edge will still move by less than maxEdgeDeviation.
   */
  private final S1ChordAngle minEdgeLengthToSplitChordAngle;

  /**
   * minSiteSeparation is set from the snap function's {@link SnapFunction#minVertexSeparation()}.
   */
  private final S1Angle minSiteSeparation;

  /**
   * To implement idempotency, we check whether the input geometry could possibly be the output of a
   * previous S2Builder invocation. This involves testing whether any site/site or edge/site pairs
   * are too close together. This checking is done using exact predicates, which requires converting
   * separation values to S1ChordAngles.
   *
   * <p>This is the minimum separation between sites (minSiteSeparation) as a chord angle.
   */
  private final S1ChordAngle minSiteSeparationChordAngle;
  /** The minimum separation between edges and sites as a chord angle. */
  private final S1ChordAngle minEdgeSiteSeparationChordAngle;

  /**
   * An upper bound on the distance computed by S2ClosestPointQuery where the true distance might
   * be less than minEdgeSiteSeparationChordAngle.
   */
  private final S1ChordAngle minEdgeSiteSeparationChordAngleLimit;

  /**
   * The maximum possible distance between two sites whose Voronoi regions touch, increased to
   * account for errors. (The maximum radius of each Voronoi region is edgeSnapRadius.)
   */
  private final S1ChordAngle maxAdjacentSiteSeparationChordAngle;

  /**
   * The squared sine of the edge snap radius. This is equivalent to the snap radius (squared) for
   * distances measured through the interior of the sphere to the plane containing an edge. This
   * value is used only when interpolating new points along edges (see
   * {@link #getSeparationSite(S2Point, S2Point, S2Point, int)}).
   */
  private final double edgeSnapRadiusSin2;

  /** A reference to the argument to build(). */
  private S2Error error;

  /**
   * True if snapping was requested. This is true if either snapRadius() is positive, or
   * splitCrossingEdges() is true (which implicitly requests snapping to ensure that both crossing
   * edges are snapped to the intersection point).
   */
  private final boolean snappingRequested;

  /**
   * Initially false, and set to true when it is discovered that at least one input vertex or edge
   * does not meet the output guarantees (e.g., that vertices are separated by at least
   * snapFunction.minVertexSeparation).
   */
  private boolean snappingNeeded;

  /**
   * The list of layers for this S2Builder. The last layer in the list is the current layer. All
   * edges in the S2Builder are assigned to an S2BuilderLayer. Specifically, edges are assigned to
   * the current layer when the edge is added. The specific implementation of S2BuilderLayer and the
   * corresponding options in layerOptions control how the edges in the layer are assembled into
   * output geometry.
   */
  private final List<S2BuilderLayer> layers = new ArrayList<>();

  /**
   * Each S2BuilderLayer in the 'layers' list has a corresponding GraphOptions in the 'layerOptions'
   * list that specifies how the geometry assigned to that layer is to be processed. Note that the
   * layer options may be modified by S2Builder in some cases, see
   * {@link S2BuilderGraph#processEdges()} for details.
   */
  private final List<GraphOptions> layerOptions = new ArrayList<>();

  /**
   * Each S2BuilderLayer in the 'layers' list has a corresponding entry in 'layerBegins' that
   * records the index in 'inputEdges' of the first edge assigned to the layer.
   */
  private final IntVector layerBegins = new IntVector();

  /**
   * Each S2BuilderLayer in the 'layers' list has a corresponding entry in
   * 'layerIsFullPolygonPredicate' that defines whether a polygon in that layer with no edges is
   * full or empty. The default predicate assigned by {@link #startLayer(S2BuilderLayer)} returns
   * false (empty) for all polygons. That default can be changed with
   * {@link #addIsFullPolygonPredicate(IsFullPolygonPredicate)} if desired.
   */
  private final List<IsFullPolygonPredicate> layerIsFullPolygonPredicates = new ArrayList<>();

  /**
   * A flag indicating whether labelSet has been modified since the last time labelSetId was
   * computed.
   */
  private boolean labelSetModified;

  /**
   * The input geometry to the S2Builder consists of vertices and edges. All the vertices in all
   * layers are stored in the 'inputVertices' list. Note that this list is modified as the
   * algorithm runs to include output sites.
   */
  private final List<S2Point> inputVertices = new ArrayList<>();

  /**
   * All the edges in all layers are stored in 'inputEdges', as pairs of input vertex ids, which
   * are simply indices into the 'inputVertices' list.
   */
  private final EdgeList inputEdges = new EdgeList();

  /**
   * If labels are used, then each input edge has a corresponding "label set id" (an int)
   * identifying the set of labels attached to that edge. The 'labelSetIds' defines that
   * correspondence, mapping input edge ids to label set ids. Label set ids are keys for the
   * 'labelSetLexicon'. If labels are not used, this vector will not be defined.
   */
  private final IntVector labelSetIds = new IntVector();

  /** Stores sets of labels that are assigned to one or more input edges. */
  private final IdSetLexicon labelSetLexicon = new IdSetLexicon();

  /** The current set of labels (represented as a stack). */
  private final IntVector labelSet = new IntVector();

  /**
   * The labelSetId is the integer id of the current label set, computed on demand by adding it to
   * the labelSetLexicon.
   */
  private int labelSetId;

  /**
   * The number of sites specified using forceVertex(). These sites are always at the beginning of
   * the sites list.
   */
  private int numForcedSites;

  /** The set of snapped vertex locations ("sites"). */
  private final List<S2Point> sites = new ArrayList<>();

  /** A reusable crossing edges query. Only constructed if needed, in addEdgeCrossings. */
  private S2CrossingEdgesQuery crossingEdgesQuery;

  /**
   * A map from each input edge id to the set of sites "nearby" that edge, defined as the set of
   * sites that are candidates for snapping and/or avoidance. Sites are kept sorted by
   * increasing distance from the origin of the input edge.
   *
   * <p>Once snapping is finished, this field is discarded unless edge chain simplification was
   * requested, in which case instead the sites are filtered by removing the ones that each edge
   * was snapped to, leaving only the "sites to avoid" (needed for simplification).
   */
  private final EdgeSites edgeSites = new EdgeSites();

  /** The S2Builder constructor takes a Builder. Clients should use Builder.build(). */
  private S2Builder(Builder builder) {
    options = builder;

    SnapFunction snapFunction = builder.snapFunction();
    S1Angle snapRadius = snapFunction.snapRadius();
    Preconditions.checkArgument(snapRadius.lessOrEquals(S2BuilderSnapFunctions.maxSnapRadius()));

    // Convert the snap radius to an S1ChordAngle. This is the "true snap radius" used when
    // evaluating exact predicates (S2Predicates).
    siteSnapRadiusChordAngle = S1ChordAngle.fromS1Angle(snapRadius);

    // When intersectionTolerance() is non-zero we need to use a larger snap radius for edges than
    // for vertices to ensure that both edges are snapped to the edge intersection location. This is
    // because the computed intersection point is not exact; it may be up to intersectionTolerance()
    // away from its true position. The computed intersection point might then be snapped to some
    // other vertex up to snapRadius away. So to ensure that both edges are snapped to a common
    // vertex, we need to increase the snap radius for edges to at least the sum of these two values
    // (calculated conservatively).
    S1Angle edgeSnapRadius = options.edgeSnapRadius();
    edgeSnapRadiusChordAngle = roundUp(edgeSnapRadius);
    snappingRequested = edgeSnapRadius.greaterThan(S1Angle.ZERO);

    // Compute the maximum distance that a vertex can be separated from an edge while still
    // affecting how that edge is snapped.
    maxEdgeDeviation = options.maxEdgeDeviation();
    edgeSiteQueryRadiusChordAngle = S1ChordAngle.fromS1Angle(
        maxEdgeDeviation.add(snapFunction.minEdgeVertexSeparation()));

    // Compute the maximum edge length such that even if both endpoints move by the maximum distance
    // allowed (i.e. edgeSnapRadius), the center of the edge will still move by less than
    // maxEdgeDeviation(). This saves us a lot of work since then we don't need to check the actual
    // deviation.
    if (!snappingRequested) {
      minEdgeLengthToSplitChordAngle = S1ChordAngle.INFINITY;
    } else {
      // This value varies between 30 and 50 degrees depending on the snap radius.
      minEdgeLengthToSplitChordAngle = S1ChordAngle.fromRadians(
          2 * acos(sin(edgeSnapRadius.radians()) / sin(maxEdgeDeviation.radians())));
    }

    // In rare cases we may need to explicitly check that the input topology is preserved, i.e. that
    // edges do not cross vertices when snapped. This is only necessary (1) for vertices added using
    // forceVertex(), and (2) when the snap radius is smaller than intersectionTolerance() (which is
    // typically either zero or S2.kIntersectionError, about 9e-16 radians). This condition arises
    // because when a geodesic edge is snapped, the edge center can move further than its endpoints.
    // This can cause an edge to pass on the wrong side of an input vertex. (Note that this could
    // not happen in a planar version of this algorithm.) Usually we don't need to consider this
    // possibility explicitly, because if the snapped edge passes on the wrong side of a vertex then
    // it is also closer than minEdgeVertexSeparation() to that vertex, which will cause a
    // separation site to be added.
    //
    // If the condition below is true then we need to check all sites (i.e. snapped input vertices)
    // for topology changes. However, this is almost never the case because:
    // <pre>
    //            maxEdgeDeviation() == 1.1 * edgeSnapRadius()
    //      and   minEdgeVertexSeparation() >= 0.219 * snapRadius()
    // </pre>
    // for all currently implemented snap functions. The condition below is only true when
    // intersectionTolerance() is non-zero (which causes edgeSnapRadius() to exceed snapRadius()
    // by S2.kIntersectionError) and snapRadius() is very small (at most
    // S2.kIntersectionError / 1.19).
    checkAllSiteCrossings =
        options
            .maxEdgeDeviation()
            .greaterThan(options.edgeSnapRadius().add(snapFunction.minEdgeVertexSeparation()));

    if (options.intersectionTolerance().lessOrEquals(S1Angle.ZERO)) {
      Preconditions.checkState(!checkAllSiteCrossings);
    }

    // To implement idempotency, we check whether the input geometry could possibly be the output
    // of a previous S2Builder invocation. This involves testing whether any site/site or edge/
    // site pairs are too close together. This is done using exact predicates, which require
    // converting the minimum separation values to an S1ChordAngle.
    minSiteSeparation = snapFunction.minVertexSeparation();
    minSiteSeparationChordAngle = S1ChordAngle.fromS1Angle(minSiteSeparation);
    minEdgeSiteSeparationChordAngle = S1ChordAngle.fromS1Angle(
        snapFunction.minEdgeVertexSeparation());

    // This is an upper bound on the distance computed by S2ClosestPointQuery where the true
    // distance might be less than minEdgeSiteSeparation.
    minEdgeSiteSeparationChordAngleLimit = addPointToEdgeError(minEdgeSiteSeparationChordAngle);

    // Compute the maximum possible distance between two sites whose Voronoi regions touch. (The
    // maximum radius of each Voronoi region is edgeSnapRadius). Then increase this bound to
    // account for errors.
    maxAdjacentSiteSeparationChordAngle = addPointToPointError(roundUp(edgeSnapRadius.mul(2)));

    // Finally, we also precompute sin^2(edgeSnapRadius), which is simply the squared distance
    // between a vertex and an edge measured perpendicular to the plane containing the edge, and
    // increase this value by the maximum error in the calculation to compare this distance
    // against the bound.
    double d = sin(edgeSnapRadius.radians());
    edgeSnapRadiusSin2 = (d * d)
      + ((9.5 * d + 2.5 + 2 * sqrt(3)) * d + 9 * DBL_EPSILON) * DBL_EPSILON;

    // Initialize the current label set.
    labelSetId = IdSetLexicon.EMPTY_SET_ID;
    labelSetModified = false;

    // If snapping was requested, we try to determine whether the input geometry already meets the
    // output requirements. This is necessary for idempotency, and can also save work. If we
    // discover any reason that the input geometry needs to be modified, snappingNeeded is set to
    // true.
    snappingNeeded = false;
  }

  // TODO(torrey): Add getters to S2Builder so specific options can be obtained without needing to
  // call toBuilder.

  /** Returns a new S2Builder.Builder initialized with the options used by this S2Builder. */
  public Builder toBuilder() {
    return new Builder(options);
  }

  /** Returns a predicate that returns a constant value (true or false); */
  public static IsFullPolygonPredicate isFullPolygon(boolean isFull) {
    return graph -> isFull;
  }

  private static S1ChordAngle roundUp(S1Angle a) {
    S1ChordAngle ca = S1ChordAngle.fromS1Angle(a);
    return ca.plusError(ca.getS1AngleConstructorMaxError());
  }

  private static S1ChordAngle addPointToPointError(S1ChordAngle ca) {
    return ca.plusError(ca.getS2PointConstructorMaxError());
  }

  private static S1ChordAngle addPointToEdgeError(S1ChordAngle ca) {
    return ca.plusError(S2EdgeUtil.getMinDistanceMaxError(ca));
  }

  /** Convenience method that returns the current layer, which is the last one. */
  public S2BuilderLayer currentLayer() {
    return layers.get(layers.size() - 1);
  }

  /** Convenience method that returns the options for the current layer. */
  public GraphOptions currentLayerOptions() {
    return layerOptions.get(layers.size() - 1);
  }

  /**
   * Starts a new output layer. This method must be called before adding any edges to the S2Builder.
   * You may call this method multiple times to build multiple geometric objects that are snapped to
   * the same set of sites.
   *
   * <p>For example, if you have a set of contour lines, then you could put each contour line in a
   * separate layer. This keeps the contour lines separate from each other, while also ensuring that
   * no crossing edges are created when they are snapped and/or simplified. (This is not true if the
   * contour lines are snapped or simplified independently.)
   *
   * <p>Similarly, if you have a set of polygons that share common boundaries (e.g., countries), you
   * can snap and/or simplify them at the same time by putting them in different layers, while
   * ensuring that their boundaries remain consistent (i.e., no crossing edges or T-vertices are
   * introduced).
   *
   * <p>Ownership of the layer is transferred to the S2Builder. Example usage:
   *
   * {@snippet :
   * S2Polyline line1;
   * S2Polyline line2;
   * builder.startLayer(new S2PolylineLayer(line1));
   * // ... Add edges using builder.addEdge(), etc.
   * builder.startLayer(new S2PolylineLayer(line2));
   * // ... Add edges using builder.addEdge(), etc.
   * S2Error error = new S2Error();
   * builder.build(error);  // Builds "line1" & "line2"
   * Preconditions.check(error.ok());
   * }
   */
  public void startLayer(S2BuilderLayer layer) {
    layerOptions.add(layer.graphOptions());
    layerBegins.add(inputEdges.size());
    layerIsFullPolygonPredicates.add(isFullPolygon(false));
    layers.add(layer);
  }

  /**
   * Input vertices are stored in a vector, with some removal of duplicates. Edges are represented
   * as (VertexId, VertexId) pairs. All edges are stored in a single vector; each layer corresponds
   * to a contiguous range. Returns the input vertex id.
   */
  @CanIgnoreReturnValue
  private int addVertex(S2Point v) {
    // Remove duplicate vertices that follow the pattern AB, BC, CD. If we want to do anything
    // more sophisticated, either use a ValueLexicon, or sort the vertices once they have all
    // been added, remove duplicates, and update the edges.
    if (inputVertices.isEmpty() || !v.equalsPoint(inputVertices.get(inputVertices.size() - 1))) {
      inputVertices.add(v);
    }
    return inputVertices.size() - 1;
  }

  /** The first 'numForcedSites' in the sites list are forced sites. */
  boolean isForced(int siteId) {
    return siteId < numForcedSites;
  }

  /**
   * Adds a point to the current layer, as a degenerate edge. The sementics depend on the layer:
   *
   * <p>If the current layer will build points, this adds a point.
   *
   * <p>If the current layer will build lines, then this point represents a degenerate 1-dimensional
   * edge. Note that line-building layers may discard degenerate edges.
   *
   * <p>If the current layer will build polygons, then this point represents a degenerate
   * 2-dimensional loop. (Whether it is a degenerate shell or hole depends on the layer: for the
   * S2LaxPolygonLayer, it depends on if it is contained by another loop). Again, note that polygon
   * building layers may discard degenerate edges.
   */
  public void addPoint(S2Point v) {
    addEdge(v, v);
  }

  /** Adds the given edge to the current layer. */
  public void addEdge(S2Point v0, S2Point v1) {
    Preconditions.checkState(!layers.isEmpty(), "Call startLayer before adding any edges.");
    if (v0.equalsPoint(v1)
        && (currentLayerOptions().degenerateEdges() == GraphOptions.DegenerateEdges.DISCARD)) {
      return;
    }

    int j0 = addVertex(v0);
    int j1 = addVertex(v1);
    inputEdges.add(j0, j1);

    // If there are any labels, then attach them to this input edge.
    if (labelSetModified) {
      if (labelSetIds.isEmpty()) {
        // Populate the missing entries of previously added edges with empty label sets.
        // (labelSetId is still EMPTY_SET_ID in this case).
        labelSetIds.enlarge(inputEdges.size() - 1);
        labelSetIds.fill(labelSetId);
      }
      // Add the new labelSet to the lexicon and get its new labelSetId.
      labelSetId = labelSetLexicon.add(labelSet);
      // Map the newly added edge to the new labelSetId.
      labelSetIds.push(labelSetId);
      labelSetModified = false;
    } else if (!labelSetIds.isEmpty()) {
      // Map the newly added edge id to the current labelSetId.
      labelSetIds.push(labelSetId);
    }
  }

  /**
   * Adds the edges in the given polyline to the current layer. Note that an S2Polyline with one
   * vertex defines a degenerate edge.
   */
  public void addPolyline(S2Polyline polyline) {
    int n = polyline.numVertices();
    if (n == 1) {
      addPoint(polyline.vertex(0));
      return;
    }
    for (int i = 1; i < n; ++i) {
      addEdge(polyline.vertex(i - 1), polyline.vertex(i));
    }
  }

  /**
   * Adds the edges in the given loop. Note that a loop consisting of one vertex adds a single
   * degenerate edge.
   *
   * <p>If the sign() of an S2Loop is negative (i.e. the loop represents a hole within a polygon),
   * the edge directions are automatically reversed to ensure that the polygon interior is always
   * to the left of every edge.
   */
  public void addLoop(S2Loop loop) {
    // Ignore loops that do not have a boundary.
    if (loop.isEmptyOrFull()) {
      return;
    }

    // For loops that represent holes, we add the edge from vertex n-1 to vertex n-2 first. This is
    // because these edges will be assembled into a clockwise loop, which will eventually be
    // normalized in S2Polygon by calling S2Loop.invert(). S2Loop.invert() reverses the order of the
    // vertices, so to end up with the original vertex order (0, 1, ..., n-1) we need to build a
    // clockwise loop with vertex order (n-1, n-2, ..., 0). This is done by adding the edge
    // (n-1, n-2) first, and then ensuring that nuild() assembles loops starting from edges in the
    // order they were added.
    for (int i = 0; i < loop.numVertices(); ++i) {
      addEdge(loop.orientedVertex(i), loop.orientedVertex(i + 1));
    }
  }

  /**
   * Adds the loops in the given polygon to the current layer. Loops representing holes have their
   * edge directions automatically reversed as described for addLoop(). Note that this method does
   * not distinguish between the empty and full polygons, i.e. adding a full polygon has the same
   * effect as adding an empty one.
   */
  public void addPolygon(S2Polygon polygon) {
    for (int i = 0; i < polygon.numLoops(); ++i) {
      addLoop(polygon.loop(i));
    }
  }

  /** Adds the edges of the given shape to the current layer. */
  public void addShape(S2Shape shape) {
    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
    for (int e = 0, n = shape.numEdges(); e < n; ++e) {
      shape.getEdge(e, edge);
      addEdge(edge.a, edge.b);
    }
  }

  /**
   * If the given "vertex" is the intersection point of two edges AB and CD (as computed by {@link
   * S2EdgeUtil#getIntersection(S2Point, S2Point, S2Point, S2Point)}), this method ensures that AB
   * and CD snap to a common vertex. (Note that the common vertex may be different than "vertex" in
   * order to ensure that no pair of vertices is closer than the given snap radius). Unlike {@link
   * Builder#splitCrossingEdges()}, this method may be used to split crossing edge pairs
   * selectively.
   *
   * <p>This method can also be used to tessellate edges using {@link
   * S2EdgeUtil#interpolateAtDistance(S1Angle, S2Point, S2Point, S1Angle)} or {@link
   * S2EdgeUtil#project(S2Point, S2Point, S2Point)} provided that a suitable intersection tolerance
   * is specified (see {@link Builder#intersectionTolerance()} for details).
   *
   * <p>This method implicitly overrides the {@link Builder#idempotent()} option, since adding an
   * intersection point implies a desire to have nearby edges snapped to it even if these edges
   * already satisfy the S2Builder output guarantees. (Otherwise for example edges would never be
   * snapped to nearby intersection points when the snap radius is zero).
   *
   * <p>Note that unlike {@link #forceVertex(S2Point)}, this method maintains all S2Builder
   * guarantees regarding minimum vertex-vertex separation, minimum edge-vertex separation, and edge
   * chain simplification.
   *
   * <p>REQUIRES: options().intersectionTolerance() > S1Angle.Zero()
   *
   * <p>REQUIRES: "vertex" was computed by {@link S2EdgeUtil#getIntersection(S2Point, S2Point,
   * S2Point, S2Point)} (in order to guarantee that both edges snap to a common vertex)
   */
  public void addIntersection(S2Point vertex) {
    // It is an error to call this method without first setting intersectionTolerance() to a
    // non-zero value.
    Preconditions.checkState(options.intersectionTolerance().greaterOrEquals(S1Angle.ZERO));

    // Calling this method also overrides the idempotent() option.
    snappingNeeded = true;
    addVertex(vertex);
  }

  /**
   * For layers that are assembled into polygons, this method specifies a predicate that is called
   * when the output consists entirely of degenerate edges and/or sibling pairs. The predicate is
   * given an S2BuilderGraph containing the output edges (if any) and is responsible for deciding
   * whether this graph represents the empty polygon (possibly with degenerate shells) or the full
   * polygon (possibly with degenerate holes). Note that this cannot be determined from the output
   * edges alone; it also requires knowledge of the input geometry. (Also see
   * {@link IsFullPolygonPredicate} above.)
   *
   * <p>This method should be called at most once per layer; additional calls simply overwrite the
   * previous value for the current layer.
   *
   * <p>The default predicate simply returns false (i.e., degenerate polygons are assumed to be
   * empty). Arguably it would better to return an error in this case, but the fact is that
   * relatively few clients need to be able to construct full polygons, and it is unreasonable to
   * expect all such clients to supply an appropriate predicate.
   *
   * <p>The reason for having a predicate rather than a boolean value is that the predicate is
   * responsible for determining whether the output polygon is empty or full. In general the input
   * geometry is not degenerate, but rather collapses into a degenerate configuration due to
   * snapping and/or simplification.
   *
   * <p>TODO(user): Provide standard predicates to handle common cases, e.g. valid input
   * geometry that becomes degenerate due to snapping.
   */
  public void addIsFullPolygonPredicate(IsFullPolygonPredicate predicate) {
    layerIsFullPolygonPredicates.set(layerIsFullPolygonPredicates.size() - 1, predicate);
  }

  /**
   * Forces a vertex to be located at the given position. This can be used to prevent certain input
   * vertices from moving. However if you are trying to preserve input edges, be aware that this
   * option does not prevent edges from being split by new vertices.
   *
   * <p>Forced vertices are subject to the following limitations:
   * <ul>
   * <li> Forced vertices are never snapped. This is true even when the given position is not
   *      allowed by the given snap function (e.g. you can force a vertex at a non-S2CellId center
   *      when using S2CellIdSnapFunction). If you want to ensure that forced vertices obey the
   *      snap function restrictions, you must call snapFunction().snapPoint() explicitly.
   * <li> There is no guaranteed minimum separation between pairs of forced vertices, i.e.
   *      snapFunction().minVertexSeparation() does not apply. (This must be true because forced
   *      vertices can be placed arbitrarily.)
   * <li> There is no guaranteed minimum separation between forced vertices and non-incident edges,
   *      i.e. snapFunction().minEdgeVertexSeparation() does not apply.
   * <li> Forced vertices are never simplified away (i.e. when simplification is requested using
   *      options().simplifyEdgeChains()).
   * </li>
   * </ul>
   * <p>All other guarantees continue to hold, e.g. the input topology will always be preserved.
   */
  public void forceVertex(S2Point vertex) {
    sites.add(vertex);
  }

  /** Convenience function that clears the stack and adds a single label, which must be >= 0. */
  public void setLabel(int label) {
    Preconditions.checkArgument(label >= 0);
    labelSet.clear();
    labelSet.push(label);
    labelSetModified = true;
  }

  /** Clear the stack of labels. */
  public void clearLabels() {
    labelSet.clear();
    labelSetModified = true;
  }

  /** Add a label to the stack. The label must be >= 0. */
  public void pushLabel(int label) {
    Preconditions.checkArgument(label >= 0);
    labelSet.push(label);
    labelSetModified = true;
  }

  /** Remove a label from the stack. */
  public void popLabel() {
    labelSet.pop();
    labelSetModified = true;
  }

  /**
   * Performs the requested edge splitting, snapping, simplification, etc, and then assembles the
   * resulting edges into the requested output layers.
   *
   * <p>Returns true if all edges were assembled; otherwise sets "error" appropriately and returns
   * false. Depending on the error, some or all output layers may have been created. Automatically
   * clears the S2Builder state so that it can be reused.
   */
  public boolean build(S2Error error) {
    // Check that an S2Error object was provided here, because this is friendlier than crashing on
    // the "error.clear()" call below. It would be easy to allow (error == null) by declaring a
    // local "tmpError", but it seems better to make clients think about error handling.
    Preconditions.checkNotNull(error);
    this.error = error;
    this.error.clear();

    // Mark the end of the last layer.
    layerBegins.push(inputEdges.size());

    // See the algorithm overview at the top of this file.
    if (snappingRequested && !options.idempotent()) {
      snappingNeeded = true;
    }
    try {
      chooseSites();
    } catch (IllegalArgumentException e) {
      // If the input geometry has a NaN value in the S2Points, the S1ChordAngle constructor will
      // throw an IllegalArgumentException with the message "Invalid length2: NaN". In this case
      // we want to set an S2Error rather than throwing an exception. For consistency with the C++
      // S2 API, we use S2Error.Code.BUILDER_SNAP_RADIUS_TOO_SMALL.
      if (e.getMessage() != null && e.getMessage().contains("NaN")) {
        error.init(S2Error.Code.BUILDER_SNAP_RADIUS_TOO_SMALL, "NaN in input geometry?");
        clear();
        return false;
      }
      throw e;
    }

    buildLayers();
    clear();

    return error.ok();
  }

  /** Clears all input data so the builder may be reused. Any options specified are preserved. */
  public void clear() {
    // Note that these calls do not change vector capacities.
    inputVertices.clear();
    inputEdges.clear();
    layers.clear();
    layerOptions.clear();
    layerBegins.clear();
    layerIsFullPolygonPredicates.clear();
    labelSetIds.clear();
    labelSetLexicon.clear();
    labelSet.clear();
    labelSetModified = false;
    sites.clear();
    edgeSites.clear();
    snappingNeeded = false;
  }

  private void chooseSites() {
    if (inputVertices.isEmpty()) {
      return;
    }

    // Note that although we always create an S2ShapeIndex, often it is not actually built (because
    // this happens lazily).
    S2ShapeIndex inputEdgeIndex = new S2ShapeIndex();

    // Add all the input vertices and edges to the new index as a single new shape.
    inputEdgeIndex.add(new GraphShape(inputEdges, inputVertices));
    if (options.splitCrossingEdges()) {
      addEdgeCrossings(inputEdgeIndex);
    }

    if (snappingRequested) {
      S2PointIndex<Integer> siteIndex = new S2PointIndex<>();
      addForcedSites(siteIndex);
      chooseInitialSites(siteIndex);
      collectSiteEdges(siteIndex);
    }

    if (snappingNeeded) {
      addExtraSites(inputEdgeIndex);
    } else {
      chooseAllVerticesAsSites();
    }
  }

  /**
   * Sorts the input vertices, discards duplicates, and uses the result as the list of sites. (We
   * sort in the same order used by chooseInitialSites() to avoid inconsistencies in tests.) We also
   * assign the result back to inputVertices and update the input edges to use the new vertex
   * numbering so that InputVertexId == SiteId. This simplifies the implementation of snapEdge() for
   * this case.
   */
  private void chooseAllVerticesAsSites() {
    sites.clear();

    // Sorts the vertex ids in the order we want as described below.
    final IntVector sortedInputVertexIds = sortInputVertices();
    // A map from input vertex id to site id.
    final IntVector vmap = IntVector.ofSize(inputVertices.size());

    // Consider the sorted vertices and add them as sites if they are not duplicates.
    S2Point previousSite = S2Point.ORIGIN; // won't match any vertex.
    for (int in = 0; in < sortedInputVertexIds.size(); in++) {
      // Get the next input vertex to consider as a candidate site.
      int vertexId = sortedInputVertexIds.get(in);
      S2Point candidateSite = inputVertices.get(vertexId);

      // Check the previously added site to see if it is exactly equal. Since sortedInputVertexIds
      // is sorted first by leaf S2CellId and then by actual S2Point, exactly equal points must be
      // adjacent in the sorted list.
      if (!previousSite.equalsPoint(candidateSite)) {
        // Not a duplicate. Add this unique vertex as a site.
        sites.add(candidateSite);
        previousSite = candidateSite;
      }

      // Map the current input vertex id to the last added site id.
      vmap.set(vertexId, sites.size() - 1);
    }

    // Update the input vertices list to be the list of sites, and then renumber the input edges to
    // match.
    inputVertices.clear();
    inputVertices.addAll(sites);
    for (int edgeId = 0; edgeId < inputEdges.size(); ++edgeId) {
      int newSrcId = vmap.get(inputEdges.getSrcId(edgeId));
      int newDstId = vmap.get(inputEdges.getDstId(edgeId));
      inputEdges.set(edgeId, newSrcId, newDstId);
    }
  }

  /**
   * Sort all the input vertices in the order that we wish to consider them as candidate Voronoi
   * sites in {@link #chooseInitialSites}. Also used by {@link #chooseAllVerticesAsSites}, which
   * deduplicates the input vertices by checking consecutive sorted vertices for equality, and so
   * currently requires that input vertices be sorted so that equal vertices are adjacent.
   *
   * <p>For {@link #chooseInitialSites}, any sort order would produce correct output, so we have
   * complete flexibility in choosing the sort key. We could even leave them unsorted, although this
   * would have the disadvantage that changing the order of the input edges could cause S2Builder to
   * snap to a different set of Voronoi sites.
   *
   * <p>We have chosen to sort them primarily by S2CellId since this improves the performance of
   * many S2Builder phases (due to better spatial locality). It also allows the possibility of
   * replacing the current S2PointIndex approach with a more efficient recursive divide-and-conquer
   * algorithm.
   *
   * <p>However, sorting by leaf S2CellId alone has two small disadvantages in the case where the
   * candidate sites are densely spaced relative to the snap radius (e.g., when using the
   * IdentitySnapFunction, or when snapping to E6/E7 near the poles, or snapping to S2CellId/E6/E7
   * using a snap radius larger than the minimum value required):
   *
   * <p>First, it tends to bias the Voronoi site locations towards points that are earlier on the
   * S2CellId Hilbert curve. For example, suppose that there are two parallel rows of input vertices
   * on opposite sides of the edge between two large S2Cells, and the rows are separated by less
   * than the snap radius. Then only vertices from the cell with the smaller S2CellId are selected,
   * because they are considered first and prevent us from selecting the sites from the other cell
   * (because they are closer than "snapRadius" to an existing site).
   *
   * <p>Second, it tends to choose more Voronoi sites than necessary, because at each step we choose
   * the first site along the Hilbert curve that is at least "snapRadius" away from all previously
   * selected sites. This tends to yield sites whose "coverage discs" overlap quite a bit, whereas
   * it would be better to cover all the input vertices with a smaller set of coverage discs that
   * don't overlap as much. (This is the "geometric set cover problem", which is NP-hard.)
   *
   * <p>It is not worth going to much trouble to fix these problems, because they really aren't that
   * important (and don't affect the guarantees made by the algorithm), but here are a couple of
   * heuristics that might help:
   *
   * <p>1. Sort the input vertices by S2CellId at a coarse level (down to cells that are
   * O(snapRadius) in size), and then sort by a fingerprint of the S2Point coordinates (i.e., quasi-
   * randomly). This would retain most of the advantages of S2CellId sorting, but makes it more
   * likely that we will select sites that are further apart.
   *
   * <p>2. Rather than choosing the first uncovered input vertex and snapping it to obtain the next
   * Voronoi site, instead look ahead through following candidates in S2CellId order and choose the
   * furthest candidate whose snapped location covers all previous uncovered input vertices.
   */
  private IntVector sortInputVertices() {
    // The vertex ids, which will be sorted by cell id and returned.
    final IntVector sortedVertexIds = new IntVector();
    // The leaf cell ids corresponding to the vertex point. Only needed for sorting.
    final long[] cellIds = new long[inputVertices.size()];

    // Initialize.
    sortedVertexIds.ensureCapacity(inputVertices.size());
    for (int inputVertexId = 0; inputVertexId < inputVertices.size(); ++inputVertexId) {
      cellIds[inputVertexId] = S2CellId.fromPoint(inputVertices.get(inputVertexId)).id();
      sortedVertexIds.add(inputVertexId);
    }

    // Wrap the IntVector and ArrayList in a SortableCollection so they can be sorted together.
    // Orders by cellId, breaking ties by comparing vertex points as described above.
    Sorter.sort(new SortableCollection() {
      @Override
      public int size() {
        return sortedVertexIds.size();
      }

      @Override
      public void truncate(int start) {
        sortedVertexIds.truncate(start);
      }

      @Override
      public void swap(int a, int b) {
        sortedVertexIds.swap(a, b);
        long t = cellIds[a];
        cellIds[a] = cellIds[b];
        cellIds[b] = t;
      }

      @Override
      public boolean less(int a, int b) {
        int c = UnsignedLongs.compare(cellIds[a], cellIds[b]);
        if (c != 0) {
          return c < 0;
        }
        // Break ties by comparing vertices, so equal vertices are adjacent in the sorted list.
        // The chooseAllVerticesAsSites() method relies on this.
        S2Point vertexA = inputVertices.get(sortedVertexIds.get(a));
        S2Point vertexB = inputVertices.get(sortedVertexIds.get(b));
        return vertexA.compareTo(vertexB) < 0;
      }
    });
    return sortedVertexIds;
  }

  /**
   * Check all edge pairs for crossings, and add the corresponding intersection points to
   * inputVertices. The intersection points will be snapped and merged with the other vertices
   * during site selection.
   */
  private void addEdgeCrossings(S2ShapeIndex inputEdgeIndex) {
    inputEdgeIndex.applyUpdates();

    ArrayList<S2Point> newVertices = new ArrayList<>();
    if (crossingEdgesQuery == null) {
      crossingEdgesQuery = new S2CrossingEdgesQuery(CrossingType.INTERIOR);
    }

    crossingEdgesQuery.visitCrossingEdgePairs(
        inputEdgeIndex,
        (shapeA, edgeIdA, aEdgeSrc, aEdgeDst, shapeB, edgeIdB, bEdgeSrc, bEdgeDst, isInterior) -> {
          newVertices.add(S2EdgeUtil.getIntersection(aEdgeSrc, aEdgeDst, bEdgeSrc, bEdgeDst));
          return true;
        });

    if (newVertices.isEmpty()) {
      return;
    }
    snappingNeeded = true;
    inputVertices.addAll(newVertices);
  }

  private void addForcedSites(S2PointIndex<Integer> siteIndex) {
    // TODO(torrey): Change 'sites' to be a wrapped array of S2Points, sort() followed by
    // deduplicate is not very efficient.
    Collections.sort(sites);
    S2BuilderUtil.deduplicateSortedList(sites);

    // Add the forced sites to the site index.
    for (int siteId = 0; siteId < sites.size(); ++siteId) {
      siteIndex.add(sites.get(siteId), siteId);
    }
    numForcedSites = sites.size();
  }

  private void chooseInitialSites(S2PointIndex<Integer> siteIndex) {
    // Prepare to find all points whose distance is <= minSiteSeparationChordAngle.
    S2ClosestPointQuery<Integer> siteQuery = new S2ClosestPointQuery<>(siteIndex);
    siteQuery.setConservativeMaxDistance(minSiteSeparationChordAngle);
    List<S2ClosestPointQuery.Result<Integer>> nearbyExistingSites = new ArrayList<>();

    // Apply the snapFunction() to each input vertex, then check whether any existing site is closer
    // than minVertexSeparation(). If not, then add a new site.
    //
    // NOTE(ericv): There are actually two reasonable algorithms, which we call "snap first" (the
    // one above) and "snap last". The latter checks for each input vertex whether any existing
    // site is closer than snapRadius(), and only then applies the snapFunction() and adds a new
    // site. "Snap last" can yield slightly fewer sites in some cases, but it is also more expensive
    // and can produce surprising results. For example, if you snap the polyline "0:0, 0:0.7" using
    // IntLatLngSnapFunction(0), the result is "0:0, 0:0" rather than the expected "0:0, 0:1",
    // because the snap radius is approximately sqrt(2) degrees and therefore it is legal to snap
    // both input points to "0:0". "Snap first" produces "0:0, 0:1" as expected.

    IntVector inputVertexIds = sortInputVertices();
    for (int i = 0; i < inputVertexIds.size(); ++i) {
      int inputVertexId = inputVertexIds.get(i);
      S2Point vertex = inputVertices.get(inputVertexId);
      S2Point site = snapSite(vertex);
      // If any vertex moves when snapped, the output cannot be idempotent.
      snappingNeeded = snappingNeeded || !site.equalsPoint(vertex);

      boolean addSite = true;
      if (siteSnapRadiusChordAngle.isZero()) {
        // If the snap radius is zero, always add, except avoid obvious duplicates.
        addSite = sites.isEmpty() || !site.equalsPoint(sites.get(sites.size() - 1));
      } else {
        // findClosestPoints() measures distances conservatively, so we need to recheck the
        // distances using exact predicates.
        //
        // NOTE(ericv): When the snap radius is large compared to the average vertex spacing, we
        // could possibly avoid the call the findClosestPoints by checking whether sites.back() is
        // close enough.
        // TODO(torrey): Update ClosestPointQuery to return nearbyExistingSites via a visitor.
        nearbyExistingSites.clear();
        // Check if site is close to previously added sites
        siteQuery.findClosestPoints(nearbyExistingSites, site);
        for (S2ClosestPointQuery.Result<Integer> result : nearbyExistingSites) {
          if (S2Predicates.compareDistance(site, result.entry().point(),
                  minSiteSeparationChordAngle.getLength2()) <= 0) {
            addSite = false;
            // This pair of sites is too close. If the sites are distinct, then the output cannot
            // be idempotent.
            snappingNeeded = snappingNeeded || !site.equalsPoint(result.entry().point());
          }
        }
      }
      if (addSite) {
        siteIndex.add(site, sites.size());
        sites.add(site);
        siteQuery.reset();
      }
    }
  }

  private S2Point snapSite(S2Point point) {
    if (!snappingRequested) {
      return point;
    }
    S2Point site = options.snapFunction().snapPoint(point);
    S1ChordAngle distMoved = new S1ChordAngle(site, point);
    if (distMoved.greaterThan(siteSnapRadiusChordAngle)) {
      error.init(
          S2Error.Code.BUILDER_SNAP_RADIUS_TOO_SMALL,
          "Snap function moved vertex (%s, %s, %s) by %s, which is more than the "
              + "specified snap radius of %s",
          point.x,
          point.y,
          point.z,
          distMoved.radians(),
          siteSnapRadiusChordAngle.radians());
    }
    return site;
  }

  /**
   * For each edge, find all sites within edgeSiteQueryRadiusChordAngle and store them in edgeSites.
   * Also, to implement idempotency this method also checks whether the input vertices and edges may
   * already satisfy the output criteria. If any problems are found then snappingNeeded is set to
   * true.
   */
  private void collectSiteEdges(S2PointIndex<Integer> siteIndex) {
    // Find all points whose distance is <= edgeSiteQueryRadiusChordAngle.
    S2ClosestPointQuery<Integer> siteQuery = new S2ClosestPointQuery<>(siteIndex);
    siteQuery.setConservativeMaxDistance(edgeSiteQueryRadiusChordAngle);

    // Nearby site list is reused for each edge.
    ArrayList<S2ClosestPointQuery.Result<Integer>> nearbySites = new ArrayList<>();

    for (int inputEdgeId = 0; inputEdgeId < inputEdges.size(); ++inputEdgeId) {
      // Find points near this input edge.
      S2Point v0 = inputVertices.get(inputEdges.getSrcId(inputEdgeId));
      S2Point v1 = inputVertices.get(inputEdges.getDstId(inputEdgeId));

      // TODO(torrey): Change ClosestPointQuery to return results with a visitor.
      nearbySites.clear(); // S2ClosestPointQuery doesn't clear before appending!
      siteQuery.findClosestPointsToEdge(nearbySites, v0, v1);

      // Gather the site ids of those points, and determine if the points need to be snapped.
      IntVector sitesNearEdge = new IntVector();
      sitesNearEdge.ensureCapacity(nearbySites.size());
      for (S2ClosestPointQuery.Result<Integer> result : nearbySites) {
        sitesNearEdge.add(result.entry().data());
        if (!snappingNeeded
            && result.distance().lessThan(minEdgeSiteSeparationChordAngleLimit)
            && !result.entry().point().equalsPoint(v0)
            && !result.entry().point().equalsPoint(v1)) {

          // Here, we may need to change snappingNeeded to true if the current graph could not be
          // the output of a previous S2Builder operation. The conservative result from the
          // S2ClosestPointQuery indicates the nearby site, which is not an edge endpoint, _may_ be
          // too close to the edge. Use an exact predicate to determine if the nearby site is
          // actually too close to the edge.
          snappingNeeded = S2Predicates.compareEdgeDistance(
              result.entry().point(), v0, v1, minEdgeSiteSeparationChordAngle.getLength2()) < 0;
        }
      }

      // Sort the sites near the edge by increasing distance from the edge source vertex.
      sitesNearEdge.sort(new SiteIdDistanceComparator(v0));
      // Store the sorted sites in edgeSites, which takes ownership of sitesNearEdge.
      edgeSites.add(inputEdgeId, sitesNearEdge);
    }
  }

  /** Sort a list of site IDs by increasing distance from a specified point. */
  private class SiteIdDistanceComparator implements IntComparator {
    private final S2Point xPoint;

    public SiteIdDistanceComparator(S2Point p) {
      xPoint = p;
    }

    @Override
    public int compare(int a, int b) {
      S2Point aPoint = sites.get(a);
      S2Point bPoint = sites.get(b);
      // Returns -1, 0, or +1 according to whether AX < BX, A == B, or AX > BX respectively.
      return S2Predicates.compareDistances(xPoint, aPoint, bPoint);
    }
  }

  /**
   * There are two situations where we need to add extra Voronoi sites in order to ensure that the
   * snapped edges meet the output requirements:
   * <pre>
   *
   *  (1) If a snapped edge deviates from its input edge by more than maxEdgeDeviation(), we add a
   *      new site on the input edge near the middle of the snapped edge. This causes the snapped
   *      edge to split into two pieces, so that it follows the input edge more closely.
   *
   *  (2) If a snapped edge is closer than minEdgeVertexSeparation() to any nearby site (the "site
   *      to avoid") or passes on the wrong side of it relative to the input edge, then we add a
   *      new site (the "separation site") along the input edge near the site to avoid. This causes
   *      the snapped edge to follow the input edge more closely, so that it is guaranteed to pass
   *      on the correct side of the site to avoid with a separation of at least the required
   *      distance.
   *
   * </pre>
   *
   * We check these conditions by snapping all the input edges to a chain of Voronoi sites and then
   * testing each edge in the chain. If a site needs to be added, we mark all nearby edges for
   * re-snapping.
   */
  private void addExtraSites(S2ShapeIndex inputEdgeIndex) {
    // Note that we could save some work in addSnappedEdges() by saving the snapped edge chains in
    // a vector, but currently this is not worthwhile since snapEdge() accounts for less than 5% of
    // the runtime (in C++). TODO(torrey): Runtime analysis for Java and update this comment.

    // A set of input edge ids that need to be re-snapped.
    IntSet edgesToResnap = new IntOpenHashSet();
    IntVector siteIdChain = new IntVector(); // Temporary storage.

    // The first pass is different because we snap every edge. In the following passes we only snap
    // edges that are near the extra sites that were added.
    for (int inputEdgeId = 0; inputEdgeId < inputEdges.size(); ++inputEdgeId) {
      snapEdge(inputEdgeId, siteIdChain);
      edgesToResnap.remove(inputEdgeId);
      maybeAddExtraSites(inputEdgeId, siteIdChain, inputEdgeIndex, edgesToResnap);
    }

    while (!edgesToResnap.isEmpty()) {
      HashSet<Integer> edgesToSnap = new HashSet<>(edgesToResnap);
      edgesToResnap.clear();
      for (int inputEdgeId : edgesToSnap) {
        snapEdge(inputEdgeId, siteIdChain);
        edgesToResnap.remove(inputEdgeId);
        maybeAddExtraSites(inputEdgeId, siteIdChain, inputEdgeIndex, edgesToResnap);
      }
    }
  }

  /**
   * The input edge identified by 'inputEdgeId' has been snapped to a series of one or more sites,
   * identified by site ids in 'siteIdsChain'. Here, we check for various problems that can arise
   * due to snapping, and fix them by adding additional new sites to the snapped edge. Those
   * problems are: 1. The snapped edge may come too close to some other site and/or pass on the
   * wrong side of some other site. 2. The snapped edge may deviate too far from the original edge.
   */
  private void maybeAddExtraSites(
      int inputEdgeId,
      IntVector siteIdsChain,
      S2ShapeIndex inputEdgeIndex,
      IntSet edgesToResnap) {
    // If the input includes NaN vertices, snapping can produce an empty chain.
    if (siteIdsChain.isEmpty()) {
      return;
    }

    // The snapped edge 'siteIdsChain' is always a subsequence of the nearby sites (edgeSites), so
    // we walk through the two arrays in parallel looking for sites that weren't snapped. These are
    // the "sites to avoid". We also keep track of the current snapped edge, since it is the only
    // edge that can be too close or pass on the wrong side of a site to avoid. Vertices beyond the
    // chain endpoints in either direction can be ignored because only the interiors of chain edges
    // can be too close to a site to avoid.

    // The input vertices.
    S2Point a0 = inputVertices.get(inputEdges.getSrcId(inputEdgeId));
    S2Point a1 = inputVertices.get(inputEdges.getDstId(inputEdgeId));

    IntList nearbySites = edgeSites.get(inputEdgeId);
    for (int i = 0, j = 0; j < nearbySites.size(); ++j) {
      int siteId = nearbySites.get(j);
      if (siteId == siteIdsChain.get(i)) {
        // This site is a vertex of the snapped edge chain. Advance to the next nearby site.
        if (++i == siteIdsChain.size()) {
          break; // Sites beyond the end of the snapped chain can be ignored.
        }

        // Check whether this snapped edge deviates too far from its original position. If so, we
        // split the edge by adding an extra site.
        S2Point v0 = sites.get(siteIdsChain.get(i - 1));
        S2Point v1 = sites.get(siteIdsChain.get(i));
        if (new S1ChordAngle(v0, v1).lessThan(minEdgeLengthToSplitChordAngle)) {
          // This snapped edge is too short to deviate too far.
          continue;
        }

        // TODO(torrey): Possible optimization. With the identity snap function, it may often(?) be
        // that a0==v0 and a1==v1 and we could skip the distance calculation?

        if (!S2EdgeUtil.isEdgeBNearEdgeA(a0, a1, v0, v1, maxEdgeDeviation)) {
          // Some point on snapped edge (v0, v1) is too far from the input edge (a0, a1).
          // Add a new site on the input edge, positioned so that it splits the snapped edge into
          // two approximately equal pieces. Then we find all the edges near the new site
          // (including this one) and add them to the edgesToResnap queue.
          //
          // Note that with large snap radii, it is possible that the snapped edge wraps around the
          // sphere the "wrong way". To handle this we find the preferred split location by
          // projecting both endpoints of the snapped edge onto the input edge and taking their
          // midpoint.
          S2Point mid =
              S2Point.add(S2EdgeUtil.project(v0, a0, a1), S2EdgeUtil.project(v1, a0, a1))
                  .normalize();
          S2Point newSite = getSeparationSite(mid, v0, v1, inputEdgeId);
          addExtraSite(newSite, inputEdgeIndex, edgesToResnap);

          // In the case where the edge wrapped around the sphere the "wrong way", it is not safe
          // to continue checking this edge. It will be marked for resnapping, and we will come back
          // to it in the next pass.
          return;
        }
      } else {
        // This site is near the input edge but is not part of the snapped chain.
        if (i == 0) {
          continue; // Sites before the start of the chain can be ignored.
        }
        // We need to ensure that non-forced sites are separated by at least
        // minEdgeVertexSeparation() from the snapped chain. This happens automatically as part of
        // the algorithm except where there are portions of the input edge that are not within
        // edgeSnapRadius() of any site. These portions of the original edge are called "coverage
        // gaps". Therefore, if we find that a site to avoid that is too close to the snapped edge
        // chain, we can fix the problem by adding a new site (the "separation site") on the input
        // edge in the corresponding coverage gap located as closely as possible to the site to
        // avoid. This technique is guaranteed to produce the required minimum separation, and
        // the entire process of adding separation sites is guaranteed to terminate.
        S2Point siteToAvoid = sites.get(siteId);
        S2Point v0 = sites.get(siteIdsChain.get(i - 1));
        S2Point v1 = sites.get(siteIdsChain.get(i));
        boolean addSeparationSite = false;
        if (!isForced(siteId)
            && minEdgeSiteSeparationChordAngle.greaterThan(S1ChordAngle.ZERO)
            && S2Predicates.compareEdgeDistance(
                    siteToAvoid, v0, v1, minEdgeSiteSeparationChordAngle.getLength2())
                < 0) {
          addSeparationSite = true;
        }

        // Similarly, we also add a separation site whenever a snapped edge passes on the wrong
        // side of a site to avoid. Normally we don't need to worry about this, since if an edge
        // passes on the wrong side of a nearby site then it is also too close to it. However if the
        // snap radius is very small and intersectionTolerance() is non-zero then we need to check
        // this condition explicitly (see the "checkAllSiteCrossings" flag for details). We also
        // need to check this condition explicitly for forced vertices. Again, we can solve this
        // problem by adding a "separation site" in the corresponding coverage gap located as
        // closely as possible to the site to avoid.
        //
        // It is possible to show that when all points are projected onto the great circle through
        // (a0, a1), no improper crossing occurs unless the site to avoid is located between a0 and
        // a1, and also between v0 and v1.
        // TODO(torrey): Verify whether all these checks are necessary.
        if (!addSeparationSite
            && (isForced(siteId) || checkAllSiteCrossings)
            && (S2Predicates.sign(a0, a1, siteToAvoid) != S2Predicates.sign(v0, v1, siteToAvoid))
            && S2Predicates.compareEdgeDirections(a0, a1, a0, siteToAvoid) > 0
            && S2Predicates.compareEdgeDirections(a0, a1, siteToAvoid, a1) > 0
            && S2Predicates.compareEdgeDirections(a0, a1, v0, siteToAvoid) > 0
            && S2Predicates.compareEdgeDirections(a0, a1, siteToAvoid, v1) > 0) {
          addSeparationSite = true;
        }

        if (addSeparationSite) {
          // We add a new site (the separation site) in the coverage gap, on the input edge,
          // located as closely as possible to the site to avoid. Then we find all the edges near
          // the new site (including this one) and add them to the edgesToResnap queue.
          S2Point newSite = getSeparationSite(siteToAvoid, v0, v1, inputEdgeId);
          assert !siteToAvoid.equalsPoint(newSite);
          addExtraSite(newSite, inputEdgeIndex, edgesToResnap);

          // Skip the remaining sites near this chain edge, and then continue scanning this chain.
          // Note that this is safe even though the call to addExtraSite() above added a new site
          // to "nearbySites".
          for (; nearbySites.get(j + 1) != siteIdsChain.get(i); ++j) {}
        }
      }
    }
  }

  /**
   * Adds a new site, then updates "edgeSites" for all edges near the new site and adds them to
   * "edgesToResnap" for resnapping.
   */
  private void addExtraSite(S2Point newSite, S2ShapeIndex inputEdgeIndex, IntSet edgesToResnap) {
    assert sites.isEmpty() || !newSite.equalsPoint(sites.get(sites.size() - 1));
    int newSiteId = sites.size();
    sites.add(newSite);

    // Find all edges whose distance is <= edgeSiteQueryRadiusChordAngle.
    S2ClosestEdgeQuery.Builder builder = new S2ClosestEdgeQuery.Builder();
    builder.setConservativeMaxDistance(edgeSiteQueryRadiusChordAngle);
    builder.setIncludeInteriors(false);

    // Note that inputEdgeIndex contains a single shape with all the edges.
    inputEdgeIndex.applyUpdates();
    S2ClosestEdgeQuery<S1ChordAngle> query = builder.build(inputEdgeIndex);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> newSiteTarget =
        new S2ClosestEdgeQuery.PointTarget<>(newSite);
    query.findClosestEdges(newSiteTarget, (distance, shapeId, edgeId) -> {
      S2Point v0 = inputVertices.get(inputEdges.getSrcId(edgeId));
      edgeSites.insertSiteIdByDistance(edgeId, newSiteId, v0);
      edgesToResnap.add(edgeId);
      return true;
    });
  }

  /**
   * Returns a new snap site which is as close as possible to the projection of 'siteToAvoid' onto
   * the input edge, in the "coverage gap" between v0 and v1. The new site should be included in the
   * snapped site id chain for inputEdgeId between the snapped sites v0 and v1.
   */
  private S2Point getSeparationSite(S2Point siteToAvoid, S2Point v0, S2Point v1, int inputEdgeId) {
    // Define the "coverage disc" of a site S to be the disc centered at S with radius "snapRadius".
    // Similarly, define the "coverage interval" of S for an edge XY to be the intersection of XY
    // with the coverage disc of S. The SnapFunction implementations guarantee that the only way
    // that a snapped edge can be closer than minEdgeVertexSeparation() to a non-snapped site (i.e.
    // siteToAvoid) is if there is a gap in the coverage of XY near this site. We can fix this
    // problem simply by adding a new site to fill this gap, located as closely as possible to the
    // site to avoid.
    //
    // To calculate the coverage gap, we look at the two snapped sites on either side of
    // siteToAvoid, and find the endpoints of their coverage intervals. Then we place a new site in
    // the gap, located as closely as possible to the site to avoid. Note that the new site may
    // move when it is snapped by the snapFunction, but it is guaranteed not to move by more than
    // snapRadius and therefore its coverage interval will still intersect the gap.
    S2Point x = inputVertices.get(inputEdges.getSrcId(inputEdgeId));
    S2Point y = inputVertices.get(inputEdges.getDstId(inputEdgeId));
    S2Point xyDir = y.sub(x);
    S2Point n = robustCrossProd(x, y);
    // Find the point as close as possible to 'siteToAvoid' on XY.
    S2Point newSite = S2EdgeUtil.project(siteToAvoid, x, y, n);
    // Ensure the new site is in the coverage gap between v0 and v1.
    S2Point gapMin = getCoverageEndpoint(v0, n);
    S2Point gapMax = getCoverageEndpoint(v1, n.neg());
    if (newSite.sub(gapMin).dotProd(xyDir) < 0) {
      newSite = gapMin;
    } else if (gapMax.sub(newSite).dotProd(xyDir) < 0) {
      newSite = gapMax;
    }
    newSite = snapSite(newSite);
    assert !v0.equalsPoint(newSite);
    assert !v1.equalsPoint(newSite);
    return newSite;
  }

  /**
   * Given a site P and a normal N for an edge XY, this method computes the intersection of XY with
   * the disc of radius snapRadius() around P, and returns the intersection point that is further
   * along the edge XY toward Y. Note that the endpoints X and Y are not required.
   */
  private S2Point getCoverageEndpoint(S2Point p, S2Point n) {
    // Consider the plane perpendicular to P that cuts off a spherical cap of radius snapRadius().
    // This plane intersects the plane through the edge XY (perpendicular to N) along a line, and
    // that line intersects the unit sphere at two points Q and R, and we want to return the point
    // R that is further along the edge XY toward Y.
    //
    // Let M be the midpoint of QR. This is the point along QR that is closest to P. We can now
    // express R as the sum of two perpendicular vectors OM and MR in the plane XY. Vector MR is in
    // the direction N x P, while vector OM is in the direction (N x P) x N, where N = X x Y.
    //
    // The length of OM can be found using the Pythagorean theorem on triangle OPM, and the length
    // of MR can be found using the Pythagorean theorem on triangle OMR.
    //
    // In the calculations below, we save some work by scaling all the vectors by
    // n.crossProd(p).norm2(), and normalizing at the end.
    double n2 = n.norm2();
    double nDp = n.dotProd(p);
    S2Point nXp = n.crossProd(p);
    S2Point nXpXn = p.mul(n2).sub(n.mul(nDp));
    S2Point om = nXpXn.mul(sqrt(1 - edgeSnapRadiusSin2));
    double mr2 = edgeSnapRadiusSin2 * n2 - nDp * nDp;

    // MR is constructed so that it points toward Y (rather than X).
    S2Point mr = nXp.mul(sqrt(max(0.0, mr2)));
    return om.add(mr).normalize();
  }

  /**
   * Given an inputEdgeId, appends to the given 'siteIdChain' the site ids of the vertices of the
   * chain of output edges that represent that inputEdgeId in the S2Builder output.
   */
  private void snapEdge(int inputEdgeId, IntVector siteIdChain) {
    siteIdChain.clear();

    int srcSiteId = inputEdges.getSrcId(inputEdgeId);
    int dstSiteId = inputEdges.getDstId(inputEdgeId);
    // If no snapping is needed, the output chain is simply the input edge source and destination
    // vertices.
    if (!snappingNeeded) {
      // Note that the input vertices have been renumbered such that InputVertexId and SiteId are
      // the same (see chooseAllVerticesAsSites).
      siteIdChain.push(srcSiteId);
      siteIdChain.push(dstSiteId);
      return;
    }
    // XY is the edge to snap.
    S2Point x = inputVertices.get(srcSiteId);
    S2Point y = inputVertices.get(dstSiteId);

    // TODO(torrey): Eric V. noted the following possibilities for optimization, applicable here
    // and in C++.
    // Optimization: if there is only one nearby site, return.
    // Optimization: if there are exactly two nearby sites, and one is close enough to each vertex,
    // then return.

    // Now iterate through the sites, which are sorted in order of distance from 'x', the origin of
    // the inputEdgeId. We keep track of the sequence of sites that are visited.
    IntList candidates = edgeSites.get(inputEdgeId);
    for (OfInt siteIter = candidates.intIterator(); siteIter.hasNext(); ) {
      int siteId = siteIter.nextInt();
      // 'c' is the next candidate site to consider.
      S2Point c = sites.get(siteId);
      // Skip any sites that are too far away from the edge to snap. (There will be some of these,
      // because we also keep track of "sites to avoid".)  Note that some sites may be close enough
      // to the line containing the edge, but not to the edge itself, so we can just use the dot
      // product with the edge normal.
      if (S2Predicates.compareEdgeDistance(c, x, y, edgeSnapRadiusChordAngle.getLength2()) > 0) {
        continue;
      }

      // Site C will be added to the chain, unless the new site C is excluded by the previous site
      // B. It is also possible for B to be excluded by C, in which case B is removed, and
      // the loop repeats with the next-previous site, and so on.
      boolean addSiteC = true;
      // Consider the sites added so far, if any, starting with the last and working backward.
      for (; !siteIdChain.isEmpty(); siteIdChain.pop()) {
        S2Point b = sites.get(siteIdChain.peek());

        // First, check whether B and C are so far apart that their clipped Voronoi regions can't
        // intersect.
        S1ChordAngle bc = new S1ChordAngle(b, c);
        if (bc.greaterOrEquals(maxAdjacentSiteSeparationChordAngle)) {
          // Neither site could exclude the other. Add site C.
          break;
        }

        // Otherwise, we want to check whether site C prevents the Voronoi region of B from
        // intersecting XY, or vice versa. This can be determined by computing the "coverage
        // interval" (the segment of XY intersected by the coverage disc of radius snapRadius) for
        // each site. If the coverage interval  of one site contains the coverage interval of the
        // other, then the contained site can be excluded.
        S2Predicates.Excluded result =
            S2Predicates.getVoronoiSiteExclusion(b, c, x, y, edgeSnapRadiusChordAngle.getLength2());
        if (result == S2Predicates.Excluded.FIRST) {
          // Site B is excluded by C. Remove it from siteIdChain and consider C vs. the previous
          // site B in the chain.
          continue;
        }
        if (result == S2Predicates.Excluded.SECOND) {
          // Site C is excluded by B. Don't add site C, continue on to the next candidate.
          addSiteC = false;
          break;
        }

        // Neither B nor C excludes the other.
        assert result == S2Predicates.Excluded.NEITHER;
        // Otherwise, if there is a site A preceding B and C, check whether A is close enough to B
        // and C that it might further clip the Voronoi region of B, and result in excluding B.
        if (siteIdChain.size() < 2) {
          // Nothing precedes B and C. Add C and continue.
          break;
        }
        S2Point a = sites.get(siteIdChain.get(siteIdChain.size() - 2));
        S1ChordAngle ac = new S1ChordAngle(a, c);
        if (ac.greaterOrEquals(maxAdjacentSiteSeparationChordAngle)) {
          // A is too far from C to be relevant. Add C and continue.
          break;
        }

        // If triangles ABC and XYB have the same orientation, the circumcenter Z of ABC is
        // guaranteed to be on the same side of XY as B.
        int xyb = S2Predicates.sign(x, y, b);
        if (S2Predicates.sign(a, b, c) == xyb) {
          break; // The circumcenter is on the same side as B but further away. Add C and continue.
        }
        // TODO(torrey): Other possible optimizations:
        //  - if AB > maxAdjacentSiteSeparationChordAngle then keep B.
        //  - if d(B, XY) < 0.5 * min(AB, BC) then keep B.

        // If the circumcenter of ABC is on the same side of XY as B, then B is excluded by A and C
        // combined. Otherwise B is needed and we can exit.
        if (S2Predicates.edgeCircumcenterSign(x, y, a, b, c) != xyb) {
          break; // B is needed. Add C and continue.
        }
        // Site B is excluded by A and C combined. Remove it and consider the previous site B.
      }
      if (addSiteC) {
        siteIdChain.push(siteId);
      }
    }
    // The endpoints of the edge must snap to at least one site.
    assert !siteIdChain.isEmpty();

    // Verify that the last site on the chain is the closest candidate site to the edge endpoint Y.
    S2Point lastSite = sites.get(siteIdChain.peek());
    assert checkSnappingInvariant(y, lastSite, candidates);
  }

  private boolean checkSnappingInvariant(S2Point endpoint, S2Point lastSite, IntList candidates) {
    // If the distance from the input edge endpoint to the last site on the chain is greater
    // than the distance from the endpoint to any other candidate site, the invariant is broken.
    for (OfInt i = candidates.intIterator(); i.hasNext(); ) {
      if (S2Predicates.compareDistances(endpoint, lastSite, sites.get(i.nextInt())) > 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Snaps and simplifies the edges for all layers. Then for each layer, constructs an
   * S2BuilderGraph with the output edges and vertices for that layer, and passes the graph to
   * layer.build to construct the layer output.
   */
  private void buildLayers() {
    // For each layer, the list of output edges in that layer.
    ArrayList<EdgeList> layerEdges = new ArrayList<>();

    // For each layer, maps output edge id (in the layer) to an "input edge id set id" (an int)
    // representing the set of input edge ids that were snapped to this edge.
    ArrayList<IntVector> layerInputEdgeIdSetIds = new ArrayList<>();

    // For all layers, the actual InputEdgeIds assigned to an edge can be retrieved from
    // "inputEdgeIdSetLexicon" with the "input edge id set id" for that edge.
    IdSetLexicon inputEdgeIdSetLexicon = new IdSetLexicon();

    // If vertex filtering is used, then for each layer, the set of vertices in that layer.
    ArrayList<List<S2Point>> layerVertices = null;

    // Snaps and possibly simplifies the edges for each layer, filling in the provided lists.
    buildLayerEdges(layerEdges, layerInputEdgeIdSetIds, inputEdgeIdSetLexicon);

    // If there are a large number of layers, then we build a minimal subset of vertices for each
    // layer. This ensures that layer types that iterate over vertices will run in time proportional
    // to the size of that layer rather than the size of all layers combined.
    // TODO(torrey): Determine an appropriate value for this constant using a benchmark.
    final int kMinLayersForVertexFiltering = 10;
    if (layers.size() >= kMinLayersForVertexFiltering) {
      // Disable vertex filtering if it is disallowed by any layer. (This could be optimized, but
      // in current applications either all layers allow filtering or none of them do.)
      boolean allowVertexFiltering = true;
      for (final GraphOptions options : layerOptions) {
        allowVertexFiltering &= options.allowVertexFiltering();
      }
      if (allowVertexFiltering) {
        layerVertices = new ArrayList<>();
        layerVertices.ensureCapacity(layers.size());
        // Temporary storage reused by filterVertices.
        IntVector tmp1 = new IntVector();
        IntVector tmp2 = new IntVector();
        for (int i = 0; i < layers.size(); ++i) {
          layerVertices.add(S2BuilderGraph.filterVertices(sites, layerEdges.get(i), tmp1, tmp2));
        }
      }
    }

    // TODO(torrey): Reuse the same S2BuilderGraph for building all layers. Or, the last thing
    // buildLayerEdges() does is create an S2BuilderGraph and process the layer edges with it
    // for each layer. Could that step be factored out to reuse the graphs here?
    for (int i = 0; i < layers.size(); ++i) {
      List<S2Point> vertices = (layerVertices == null ? sites : layerVertices.get(i));
      S2BuilderGraph graph =
          new S2BuilderGraph(
              layerOptions.get(i),
              vertices,
              layerEdges.get(i),
              layerInputEdgeIdSetIds.get(i),
              inputEdgeIdSetLexicon,
              labelSetIds,
              labelSetLexicon,
              layerIsFullPolygonPredicates.get(i));
      // Building a layer may fail, in which case 'error' is filled in.
      boolean unused = layers.get(i).build(graph, error);
    }
  }

  /**
   * Snaps and possibly simplifies the edges for each layer, populating the given output arguments.
   * The resulting edges can be used to construct an S2BuilderGraph directly (no further processing
   * is necessary). The output arguments are:
   * <ul>
   * <li>For each layer, an EdgeList containing the output edges in that layer.</li>
   * <li>For each layer, an IntVector mapping edge id to an id in the lexicon, containing sets of
   * input edge ids, thus providing the input edge ids mapped to each output edge.</li>
   * </ul>
   * <p>Note S2BuilderGraph.processEdges can modify layerOptions in some cases, changing undirected
   * edges to directed ones.
   */
  private void buildLayerEdges(
      ArrayList<EdgeList> layerEdges,
      ArrayList<IntVector> layerInputEdgeIdSetIds,
      IdSetLexicon inputEdgeIdSetLexicon) {
    // Edge chains are simplified only when a non-zero snap radius is specified.
    boolean simplify = snappingNeeded && options.simplifyEdgeChains();

    // If simplifying, we build a map from each site to the set of input vertices that snapped to
    // that site. The common case is exactly one input vertex per site. TODO(torrey): An IntVector
    // for every site, commonly holding a single input vertex id is inefficient, fix that.
    ArrayList<IntVector> siteVertices = new ArrayList<>();
    if (simplify) {
      // siteVertices is initialized here and filled in by addSnappedEdges.
      for (int i = 0; i < sites.size(); ++i) {
        siteVertices.add(new IntVector());
      }
    }

    for (int i = 0; i < layers.size(); ++i) {
      layerEdges.add(new EdgeList());
      layerInputEdgeIdSetIds.add(new IntVector());
    }

    for (int layer = 0; layer < layers.size(); ++layer) {
      // Snap the edges for layer 'i', filling in layerEdges and layerInputEdgeIdSetIds, while
      // updating inputEdgeIdSetLexicon and filling in siteVertices if required.
      addSnappedEdges(
          layerBegins.get(layer),
          layerBegins.get(layer + 1),
          layerOptions.get(layer),
          layerEdges.get(layer),
          layerInputEdgeIdSetIds.get(layer),
          inputEdgeIdSetLexicon,
          siteVertices);
    }
    // dumpLayerEdges("Snapped", layerEdges, layerInputEdgeIdSetIds, inputEdgeIdSetLexicon);

    // We simplify edge chains before processing the per-layer GraphOptions because simplification
    // can create duplicate edges and/or sibling edge pairs which may need to be removed.
    if (simplify) {
      simplifyEdgeChains(siteVertices, layerEdges, layerInputEdgeIdSetIds, inputEdgeIdSetLexicon);
      siteVertices.clear();
    }
    // dumpLayerEdges("Simplified", layerEdges, layerInputEdgeIdSetIds, inputEdgeIdSetLexicon);

    // At this point we have no further need for nearby site data, so we clear it to save space.
    // We keep inputVertices and inputEdges so that S2Builder.Layer implementations can access them
    // if desired. (This is useful for determining how snapping has changed the input geometry.)
    edgeSites.clear();

    for (int i = 0; i < layers.size(); ++i) {
      // The errors generated by processEdges are really warnings, so we simply record them and
      // continue.
      S2BuilderGraph.processEdges(
          layerOptions.get(i),
          layerEdges.get(i),
          layerInputEdgeIdSetIds.get(i),
          inputEdgeIdSetLexicon,
          error);
    }
    // dumpLayerEdges("Processed", layerEdges, layerInputEdgeIdSetIds, inputEdgeIdSetLexicon);
  }

  /** Simplifies edge chains, updating its input/output arguments as necessary. */
  private void simplifyEdgeChains(
      ArrayList<IntVector> siteVertices,
      ArrayList<EdgeList> layerEdges,
      ArrayList<IntVector> layerInputEdgeIdSetIds,
      IdSetLexicon inputEdgeIdSetLexicon) {
    if (layers.isEmpty()) {
      return;
    }

    // An EdgeList with edges from all layers, ordered lexicographically.
    EdgeList mergedEdges = new EdgeList();
    // Maps edges in mergedEdges to the corresponding IdSet id in the lexicon.
    IntVector mergedInputEdgeIdSetIds = new IntVector();
    // Maps edges in mergedEdges to the layer the edge is from.
    IntVector mergedEdgeLayers = new IntVector();

    // Merge the edges from all layers (in order to build a single graph).
    mergeLayerEdges(
        layerEdges, layerInputEdgeIdSetIds, mergedEdges, mergedInputEdgeIdSetIds, mergedEdgeLayers);

    // The following fields will be reconstructed by EdgeChainSimplifier.
    for (EdgeList edges : layerEdges) {
      edges.clear();
    }
    for (IntVector inputEdgeIds : layerInputEdgeIdSetIds) {
      inputEdgeIds.clear();
    }

    // The graph options are irrelevant for edge chain simplification, but we try to set them
    // appropriately anyway.
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.DIRECTED, DegenerateEdges.KEEP, DuplicateEdges.KEEP, SiblingPairs.KEEP);

    // Construct a new S2BuilderGraph for the EdgeChainSimplifier, containing the merged edges from
    // all the layers in lexicographic order.
    S2BuilderGraph mergedGraph =
        new S2BuilderGraph(
            graphOptions,
            sites,
            mergedEdges,
            mergedInputEdgeIdSetIds,
            inputEdgeIdSetLexicon,
            null,
            null,
            isFullPolygon(false));

    EdgeChainSimplifier simplifier =
        new EdgeChainSimplifier(
            this,
            mergedGraph,
            mergedEdgeLayers,
            siteVertices,
            layerEdges,
            layerInputEdgeIdSetIds,
            inputEdgeIdSetLexicon);

    simplifier.run();
  }

  /**
   * Merges the edges from all layers and sorts them in lexicographic order so that we can construct
   * a single graph. The sort is stable because all ties are broken using edgeIds if needed, which
   * means that any duplicate edges within each layer will still be sorted by InputEdgeId.
   *
   * <p>The input parameters are layerEdges, which is an EdgeList for each layer, and
   * layerInputEdgeIdSetIds, with a map from EdgeList index to IdSet id for each layer.
   *
   * <p>Fills in the provided output parameters mergedEdges, mergedInputEdgeIdSetIds, and
   * mergedEdgelayers with the data from the edges in all layers, in lexicographic order.
   */
  private void mergeLayerEdges(
      // Input parameters.
      ArrayList<EdgeList> layerEdges,
      ArrayList<IntVector> layerInputEdgeIdSetIds,
      // Output parameters.
      EdgeList mergedEdges,
      IntVector mergedInputEdgeIdSetIds,
      IntVector mergedEdgeLayers) {
    // A LayerEdgeList containing (layer, edgeIndex) for every edge in every layer's edgeList.
    LayerEdgeIdList order = new LayerEdgeIdList(layerEdges);
    for (int layer = 0; layer < layerEdges.size(); ++layer) {
      for (int edgeIndex = 0; edgeIndex < layerEdges.get(layer).size(); ++edgeIndex) {
        order.add(layer, edgeIndex);
      }
    }

    // Sort by source vertex id, then destination vertex id, then layer id, then edge index.
    order.sort();

    mergedEdges.ensureCapacity(order.size());
    mergedInputEdgeIdSetIds.ensureCapacity(order.size());
    mergedEdgeLayers.ensureCapacity(order.size());

    // Iterate over the LayerEdgeId pairs in sorted order and fill in the output parameters.
    for (int i = 0; i < order.size(); i++) {
      mergedEdges.add(order.getSrcId(i), order.getDstId(i));
      mergedInputEdgeIdSetIds.add(
          layerInputEdgeIdSetIds.get(order.getLayerId(i)).get(order.getEdgeId(i)));
      mergedEdgeLayers.add(order.getLayerId(i));
    }
  }

  /**
   * Snaps all the input edges for a given layer, populating the given output arguments. If
   * siteVertices is non-empty then it is updated so that siteVertices.get(site) contains a list of
   * all input vertices that were snapped to that site.
   */
  private void addSnappedEdges(
      // Input parameters
      int beginInputEdgeId,
      int endInputEdgeId,
      GraphOptions options,
      // Output parameters
      EdgeList edges,
      IntVector inputEdgeIdSetIds,
      IdSetLexicon inputEdgeIdSetLexicon,
      ArrayList<IntVector> siteVertices) {
    boolean discardDegenerateEdges =
        options.degenerateEdges() == GraphOptions.DegenerateEdges.DISCARD;
    IntVector siteIdChain = new IntVector();
    for (int inputEdgeId = beginInputEdgeId; inputEdgeId < endInputEdgeId; ++inputEdgeId) {
      int inputEdgeIdSetId = inputEdgeIdSetLexicon.addSingleton(inputEdgeId);
      // Snap the input edge to a chain of sites.
      snapEdge(inputEdgeId, siteIdChain);
      if (siteIdChain.isEmpty()) {
        continue;
      }
      // If siteVertices is non-empty, edges will be simplified. Then ensure siteVertices for the
      // first site in the chain has the input edge source vertex.
      maybeAddInputVertex(inputEdges.getSrcId(inputEdgeId), siteIdChain.get(0), siteVertices);
      if (siteIdChain.size() == 1) {
        if (discardDegenerateEdges) {
          continue;
        }
        // Add the (degenerate) single edge in the chain to 'edges' and 'inputEdgeIdSetIds'.
        addSnappedEdge(
            siteIdChain.get(0),
            siteIdChain.get(0),
            inputEdgeIdSetId,
            options.edgeType(),
            edges,
            inputEdgeIdSetIds);
      } else {
        // If siteVertices is non-empty, ensure siteVertices for the last site in the chain has the
        // input edge destination vertex.
        maybeAddInputVertex(inputEdges.getDstId(inputEdgeId), siteIdChain.peek(), siteVertices);
        // Add each edge in the chain to 'edges' and 'inputEdgeIdSetIds'.
        for (int i = 1; i < siteIdChain.size(); ++i) {
          addSnappedEdge(
              siteIdChain.get(i - 1),
              siteIdChain.get(i),
              inputEdgeIdSetId,
              options.edgeType(),
              edges,
              inputEdgeIdSetIds);
        }
      }
    }
  }

  /**
   * If "siteVertices" is non-empty, ensures that siteVertices.get(siteId) contains "inputVertexId".
   * Duplicate entries are allowed. The purpose of this function is to build a map from site id to
   * input vertices, so that simplifyEdgeChains() can quickly find all the input vertices that
   * snapped to a particular site.
   */
  private void maybeAddInputVertex(
      int inputVertexId, int siteId, ArrayList<IntVector> siteVertices) {
    if (siteVertices.isEmpty()) {
      return;
    }

    // Optimization: check if we just added this vertex. This is worthwhile because the input edges
    // usually form a continuous chain, i.e. the destination of one edge is the same as the source
    // of the next edge.
    IntVector vertices = siteVertices.get(siteId);
    if (vertices.isEmpty() || vertices.peek() != inputVertexId) {
      vertices.push(inputVertexId);
    }
  }

  /**
   * Adds the given edge to "edges" and "inputEdgeIdSetIds". If undirected edges are being used,
   * also adds an edge in the opposite direction.
   */
  private void addSnappedEdge(
      int srcSiteId,
      int dstSiteId,
      int inputEdgeIdSetId,
      EdgeType edgeType,
      EdgeList edges,
      IntVector inputEdgeIdSetIds) {
    edges.add(srcSiteId, dstSiteId);
    inputEdgeIdSetIds.add(inputEdgeIdSetId);
    if (edgeType == EdgeType.UNDIRECTED) {
      edges.add(dstSiteId, srcSiteId);
      // Automatically created edges do not have input edge ids or labels. This can be used to
      // distinguish the original direction of the undirected edge.
      inputEdgeIdSetIds.add(IdSetLexicon.EMPTY_SET_ID);
    }
  }

  /**
   * EdgeSites is conceptually a map from edge ids to sets of site ids "nearby" that edge, defined
   * as the set of sites that are candidates for snapping and/or avoidance. Sites are kept sorted by
   * increasing distance from the origin of the input edge. Most edges have just two nearby sites,
   * i.e. the edge end points, so EdgeSites is optimized for that common case.
   */
  private class EdgeSites {
    public EdgeSites() {}
    // TODO(torrey): Implementation. Eric's comments: "we probably do want an array of some kind of
    // value, where the array index is the input edge ID, and the array values are describing the
    // endpoints of the edge and the nearby vertices.
    //
    // E.g. we can have this be a single IntVector, where the first 2*numInputEdges elements are
    // the IDs of the start/end vertices, and if the closest edge query produces more than 2 nearby
    // vertices, then the first entry can contain the first result and the second entry can contain
    // the remaining results at the end of the single IntVector (as size followed by that many
    // remaining entries). That eliminates any objects per input edge."
    // If the data structure and add() method are changed to copy the provided sitesNearEdge, then
    // collectSiteEdges can be changed to reuse an IntVector instead of allocating each one.
    private final HashMap<Integer, IntVector> map = new HashMap<>();

    /** Clears the contents of this EdgeSites. */
    public void clear() {
      map.clear();
    }

    /**
     * Stores the provided list of site ids 'sitesNearEdge' for the provided 'inputEdgeId'. The
     * 'sitesNearEdge' must already be sorted. Currently, takes ownership of the provided IntVector
     * which may not be further modified.
     */
    public void add(int inputEdgeId, IntVector sitesNearEdge) {
      map.put(inputEdgeId, sitesNearEdge);
    }

    /**
     * Inserts 'newSiteId' into the list of site ids near the edge 'edgeId'. The list is already
     * sorted by distance from point 'x' and the sorted order is maintained.
     */
    public void insertSiteIdByDistance(int edgeId, int newSiteId, S2Point x) {
      IntVector sitesNearEdge = map.get(edgeId);
      // Returns the lowest index containing a site id greater than or equal to newSiteId, or
      // sitesNearEdge.size() if all existing site ids are lower (the normal case when adding new
      // sites.)
      int pos = sitesNearEdge.lowerBound(newSiteId, new SiteIdDistanceComparator(x));
      // The newSiteId must not have already been added.
      Preconditions.checkState(pos == sitesNearEdge.size() || sitesNearEdge.get(pos) != newSiteId);
      sitesNearEdge.insert(pos, newSiteId);
    }

    /** Returns the sorted list of site ids for the given 'edgeId'. */
    public IntList get(int edgeId) {
      if (!map.containsKey(edgeId)) {
        map.put(edgeId, new IntVector());
      }
      return map.get(edgeId);
    }
  }

  /** A class that encapsulates the state needed for simplifying edge chains. */
  private static class EdgeChainSimplifier {
    /** The S2Builder that owns this EdgeChainSimplifier. */
    private final S2Builder builder;

    /**
     * The S2BuilderGraph that this EdgeChainSimplifier is working on, containing the merged edges
     * from all layers.
     */
    private final S2BuilderGraph graph;

    /** Constructed by this EdgeChainSimplifier, maps graph vertices to their incoming edge ids. */
    private final S2BuilderGraph.VertexInMap vertexInMap;

    /** Constructed by this EdgeChainSimplifier, maps graph vertices to their outgoing edge ids. */
    private final S2BuilderGraph.VertexOutMap vertexOutMap;

    /** Maps edgeIds to the layer they came from. */
    private final IntVector edgeLayers;

    /** Maps site ids to all the input vertex ids that snapped to the site. */
    private final List<IntVector> siteVertices;

    /** For each layer, the edges in the layer. Modified by simplification. */
    private final List<EdgeList> layerEdges;

    /** For each layer, maps edges in the EdgeList to the IdSet in the lexicon of input edge ids. */
    private final List<IntVector> layerInputEdgeIdSetIds;

    /** For all layers, sets of input edge ids mapped to output edges. */
    private final IdSetLexicon inputEdgeIdSetLexicon;

    /** Convenience member copied from `builder`. */
    private final IntVector layerBegins;

    /**
     * Constants for clarity with the 'isInterior' and 'used' IntVectors, which are storing
     * booleans as ints.
     */
    private static final int TRUE = 1;
    private static final int FALSE = 0;

    /**
     * isInterior[v]==TRUE indicates that VertexId "v" is eligible to be an interior vertex of a
     * simplified edge chain. You can think of it as vertex whose indegree and outdegree are both 1
     * (although the actual definition is a bit more complicated because of duplicate edges and
     * layers).
     */
    private final IntVector isInterior;

    /** used[e]==TRUE indicates that EdgeId "e" has already been processed. */
    private final IntVector used;

    // Temporary objects declared here to avoid repeated allocation.
    private final IntVector tmpVertices = new IntVector();
    private final IntVector tmpEdges = new IntVector();
    // Temporary storage used only in simplifyChain().
    private final IntSet usedVertices = new IntOpenHashSet();

    /** New output edges created by simplification. */
    private final EdgeList newEdges = new EdgeList();

    /** Maps new edges by index in 'newEdges' to ids of sets of input edge ids in the lexicon. */
    private final IntVector newInputEdgeIdSetIds = new IntVector();

    /** Maps new edges by index in 'newEdges' to the layer of the new edge. */
    private final IntVector newEdgeLayers = new IntVector();

    /**
     * Constructs an EdgeChainSimplifier. The graph "g" contains all edges from all layers.
     * "edgeLayers" indicates the original layer for each edge. "siteVertices" is a map from SiteId
     * to the set of InputVertexIds that were snapped to that site. "layerEdges" and
     * "layerInputEdgeIdSetIds" are output arguments where the simplified edge chains will be
     * placed. The input and output edges are not sorted.
     */
    public EdgeChainSimplifier(
        // Input parameters.
        S2Builder builder,
        S2BuilderGraph graph,
        IntVector edgeLayers,
        List<IntVector> siteVertices,
        // Output parameters.
        List<EdgeList> layerEdges, // Simplified edges for each layer.
        List<IntVector> layerInputEdgeIdSetIds, // InputEdgeIdSetIds for the new edges.
        IdSetLexicon inputEdgeIdSetLexicon) { // New lexicon entries are added, if needed.
      this.builder = builder;
      this.graph = graph;
      this.edgeLayers = edgeLayers;

      // Unmodified, only used to ensure the simplified edge is close enough to the input vertices.
      this.siteVertices = siteVertices;

      vertexInMap = new S2BuilderGraph.VertexInMap(graph);
      vertexOutMap = new S2BuilderGraph.VertexOutMap(graph);

      this.layerEdges = layerEdges;
      this.layerInputEdgeIdSetIds = layerInputEdgeIdSetIds;
      this.inputEdgeIdSetLexicon = inputEdgeIdSetLexicon;
      this.layerBegins = builder.layerBegins;

      isInterior = IntVector.ofSize(graph.numVertices());
      used = IntVector.ofSize(graph.numEdges());

      newEdges.ensureCapacity(graph.numEdges());
      newInputEdgeIdSetIds.ensureCapacity(graph.numEdges());
      newEdgeLayers.ensureCapacity(graph.numEdges());
    }

    public void run() {
      // Determine which vertices can be interior vertices of an edge chain.
      for (int vertexId = 0; vertexId < graph.numVertices(); ++vertexId) {
        isInterior.set(vertexId, isInterior(vertexId));
      }

      // Attempt to simplify all edge chains that start from a non-interior vertex. (This takes
      // care of all chains except loops.)
      for (int edgeId = 0; edgeId < graph.numEdges(); ++edgeId) {
        int srcVertex = graph.edgeSrcId(edgeId);
        int dstVertex = graph.edgeDstId(edgeId);
        if (used.get(edgeId) == TRUE) {
          continue;
        }
        if (isInterior.get(srcVertex) == TRUE) {
          continue;
        }
        if (isInterior.get(dstVertex) == FALSE) {
          outputEdge(edgeId); // An edge between two non-interior vertices.
        } else {
          simplifyChain(srcVertex, dstVertex);
        }
      }

      // If there are any edges left, they form one or more disjoint loops where all vertices are
      // interior vertices.
      //
      // TODO(ericv): It would be better to start from the edge with the smallest minInputEdgeId(),
      // since that would make the output more predictable for testing purposes. It also means that
      // we won't create an edge that spans the start and end of a polyline if the polyline is
      // snapped into a loop. (Unfortunately there are pathological examples that prevent us from
      // guaranteeing this in general, e.g. there could be two polylines in different layers that
      // snap to the same loop but start at different positions. In general we only consider input
      // edge ids to be a hint towards the preferred output ordering.)
      for (int edgeId = 0; edgeId < graph.numEdges(); ++edgeId) {
        if (used.get(edgeId) == TRUE) {
          continue;
        }

        int srcVertex = graph.edgeSrcId(edgeId);
        int dstVertex = graph.edgeDstId(edgeId);
        if (srcVertex == dstVertex) {
          // Note that it is safe to output degenerate edges as we go along, because this vertex has
          // at least one non-degenerate outgoing edge and therefore we will (or just did) start an
          // edge chain here.
          outputEdge(edgeId);
        } else {
          simplifyChain(srcVertex, dstVertex);
        }
      }

      // TODO(torrey): The graph is not needed past here, so we could save some memory by clearing
      // the underlying Edge and InputEdgeIdSetId vectors.

      // Finally, copy the output edges into the appropriate layers. They don't need to be sorted
      // because the input edges were also unsorted.
      for (int index = 0; index < newEdges.size(); ++index) {
        int layer = newEdgeLayers.get(index);
        int src = newEdges.getSrcId(index);
        int dst = newEdges.getDstId(index);
        int newInputEdgeIdSetId = newInputEdgeIdSetIds.get(index);
        layerEdges.get(layer).add(src, dst);
        layerInputEdgeIdSetIds.get(layer).add(newInputEdgeIdSetId);
      }
    }

    /** Copies the given edge to the output and marks it as used. */
    private void outputEdge(int edgeId) {
      newEdges.add(graph.edges().getSrcId(edgeId), graph.edges().getDstId(edgeId));
      newInputEdgeIdSetIds.push(graph.inputEdgeIdSetId(edgeId));
      newEdgeLayers.push(edgeLayers.get(edgeId));
      used.set(edgeId, TRUE);
    }

    /** Returns the layer that a given graph edge belongs to. */
    private int graphEdgeLayer(int edgeId) {
      return edgeLayers.get(edgeId);
    }

    /** Returns the layer than a given input edge belongs to. */
    private int inputEdgeLayer(int inputEdgeId) {
      // TODO(torrey): If this method shows up in profiling, the result could be stored with each
      // edge (i.e., edgeLayers and newEdgeLayers).
      assert inputEdgeId >= 0;
      return layerBegins.upperBound(inputEdgeId) - (layerBegins.get(0) + 1);
    }

    /**
     * Returns TRUE if the given vertexId can be an interior vertex of a simplified edge chain. See
     * the InteriorVertexMatcher class for what this implies.
     *
     * <p>(Uses integer constants TRUE and FALSE because the result is stored in an IntVector.)
     */
    private int isInterior(int vertexId) {
      // Check a few simple prerequisites.
      if (vertexOutMap.outDegree(vertexId) == 0) {
        return FALSE;
      }
      if (vertexOutMap.outDegree(vertexId) != vertexInMap.inDegree(vertexId)) {
        return FALSE;
      }
      if (builder.isForced(vertexId)) {
        return FALSE; // Keep forced vertices.
      }

      // Sort the edge ids so that they are grouped by layer.
      IntVector edgeIds = tmpEdges; // Avoid allocating each time.
      edgeIds.clear();
      vertexOutMap.edgeIds(vertexId).forEach(edgeIds::push);
      vertexInMap.edgeIds(vertexId).forEach(edgeIds::push);
      edgeIds.sort((x, y) -> Integer.compare(graphEdgeLayer(x), graphEdgeLayer(y)));

      // Now feed the edges in each layer to the InteriorVertexMatcher.
      InteriorVertexMatcher matcher = new InteriorVertexMatcher(vertexId);
      for (int e = 0; e < edgeIds.size(); ) {
        int edgeId = edgeIds.get(e);
        int layer = graphEdgeLayer(edgeId);
        matcher.startLayer();
        for (; e < edgeIds.size() && graphEdgeLayer(edgeIds.get(e)) == layer; e++) {
          edgeId = edgeIds.get(e);
          int srcVertex = graph.edgeSrcId(edgeId);
          int dstVertex = graph.edgeDstId(edgeId);
          if (srcVertex == vertexId) {
            matcher.tally(dstVertex, /* outgoing= */ true);
          }
          if (dstVertex == vertexId) {
            matcher.tally(srcVertex, /* outgoing= */ false);
          }
        }

        if (!matcher.matches()) {
          return FALSE;
        }
      }
      return TRUE;
    }

    /**
     * Follows the edge chain starting with (v0, v1) until either we find a non-interior vertex or
     * we return to the original vertex v0. At each vertex we simplify a subchain of edges that is
     * as long as possible.
     */
    private void simplifyChain(int v0, int v1) {
      // Avoid allocating "chain" each time by reusing it.
      IntVector chain = tmpVertices;

      // Contains the set of vertices that have either been avoided or added to the chain so far.
      // This is necessary so that avoidSites() doesn't try to avoid vertices that have already
      // been added to the chain.

      S2PolylineSimplifier simplifier = new S2PolylineSimplifier();
      int vStart = v0;
      boolean done = false;
      do {
        usedVertices.clear();
        // Simplify a subchain of edges starting with (v0, v1).
        chain.push(v0);
        usedVertices.add(v0);
        simplifier.init(graph.vertex(v0));
        // Note that if the first edge (v0, v1) is longer than the maximum length allowed for
        // simplification, then avoidSites() will return false, and we exit the loop below after the
        // first iteration.
        boolean simplify = avoidSites(v0, v0, v1, usedVertices, simplifier);
        do {
          chain.push(v1);
          usedVertices.add(v1);
          done = isInterior.get(v1) == FALSE || v1 == vStart;
          if (done) {
            break;
          }

          // Attempt to extend the chain to the next vertex.
          int vPrev = v0;
          v0 = v1;
          // Get the next vertex in the chain after edge (vPrev, v0).
          v1 = followChain(vPrev, v0);
        } while (simplify
            // Constrain the simplified edge to be close enough to the input vertices snapped to it.
            && targetInputVertices(v0, simplifier)
            // Constrain the simplified edge to avoid nearby sites.
            && avoidSites(chain.get(0), v0, v1, usedVertices, simplifier)
            // The edge from the *original* v0, can be extended to v1 within the constraints.
            && simplifier.extend(graph.vertex(v1)));


        // 'chain' is now a series of sites and the edges between them that can be replaced by a
        // single edge from the first to last vertices of the chain ... except that there may be
        // many copies of 'chain', possibly in different layers, so there will be multiple copies
        // of the simplified edge, and we need to track input edge ids snapped to the simplified
        // output edge(s). The 'mergeChain' method handles those details.
        if (chain.size() == 2) {
          outputAllEdges(chain.get(0), chain.get(1)); // Could not simplify.
        } else {
          mergeChain(chain);
        }
        // Note that any degenerate edges that were not merged into a chain are output by
        // EdgeChainSimplifier.run().
        chain.clear();
      } while (!done);
    }

    /**
     * Given an edge (v0, v1) where v1 is an interior vertex, returns the (unique) next vertex id
     * in the edge chain.
     */
    private int followChain(int v0, int v1) {
      assert isInterior.get(v1) == TRUE;
      for (OfInt edgeIter = vertexOutMap.edgeIds(v1).intIterator(); edgeIter.hasNext(); ) {
        int edgeId = edgeIter.nextInt();
        int edgeDstVertex = graph.edgeDstId(edgeId);
        if (edgeDstVertex != v0 && edgeDstVertex != v1) {
          return edgeDstVertex;
        }
      }

      throw new IllegalStateException("Could not find next edge in edge chain");
    }

    /** Copies all input edges between v0 and v1 (in both directions) to the output. */
    private void outputAllEdges(int v0, int v1) {
      vertexOutMap.edgeIds(v0, v1).forEach(this::outputEdge);
      vertexOutMap.edgeIds(v1, v0).forEach(this::outputEdge);
    }

    /**
     * Ensures that the simplified edge passes within "edgeSnapRadius" of all the *input* vertices
     * that snapped to the given vertex "vertexId".
     */
    private boolean targetInputVertices(int vertexId, S2PolylineSimplifier simplifier) {
      for (OfInt iter = siteVertices.get(vertexId).intIterator(); iter.hasNext(); ) {
        int inputVertexId = iter.nextInt();
        if (!simplifier.targetDisc(
            builder.inputVertices.get(inputVertexId), builder.edgeSnapRadiusChordAngle)) {
          return false;
        }
      }
      return true;
    }

    /**
     * Given the starting vertex v0 and last edge (v1, v2) of an edge chain, restricts the allowable
     * range of angles in order to ensure that all sites near the edge (v1, v2) are avoided by at
     * least minEdgeVertexSeparation.
     *
     * <p>Returns true if sites near the edge can be avoided given previous constraints, and the
     * edge can be simplified, or false otherwise.
     */
    private boolean avoidSites(
        int v0,
        int v1,
        int v2,
        IntSet usedVertices,
        S2PolylineSimplifier simplifier) {
      S2Point p0 = graph.vertex(v0);
      S2Point p1 = graph.vertex(v1);
      S2Point p2 = graph.vertex(v2);
      S1ChordAngle r1 = new S1ChordAngle(p0, p1);
      S1ChordAngle r2 = new S1ChordAngle(p0, p2);

      // The distance from the start of the edge chain must increase monotonically for each vertex,
      // since we don't want to simplify chains that backtrack on themselves (we want a parametric
      // approximation, not a geometric one).
      if (r2.lessThan(r1)) {
        return false;
      }

      // We also limit the maximum edge length in order to guarantee that the simplified edge stays
      // with maxEdgeDeviation() of all the input edges that snap to it.
      if (r2.greaterOrEquals(builder.minEdgeLengthToSplitChordAngle)) {
        return false;
      }

      // Otherwise it is sufficient to consider the nearby sites (edgeSites) for a single input
      // edge that snapped to (v1, v2) or (v2, v1). This is because each edge has a list of all
      // sites within (maxEdgeDeviation + minEdgeVertexSeparation), and since the output edge
      // is within maxEdgeDeviation of all input edges, this list includes all sites within
      // minEdgeVertexSeparation of the output edge.
      //
      // Usually there is only one edge to choose from, but it's not much more effort to choose the
      // edge with the shortest list of edgeSites.
      int bestInputEdgeId = -1;
      EdgeSites edgeSites = builder.edgeSites;
      for (OfInt edgeIter = vertexOutMap.edgeIds(v1, v2).intIterator(); edgeIter.hasNext(); ) {
        int edgeId = edgeIter.nextInt();
        for (OfInt inputEdgeIter = graph.inputEdgeIds(edgeId).intIterator();
            inputEdgeIter.hasNext(); ) {
          int inputEdgeId = inputEdgeIter.nextInt();
          if (bestInputEdgeId < 0
              || edgeSites.get(inputEdgeId).size() < edgeSites.get(bestInputEdgeId).size()) {
            bestInputEdgeId = inputEdgeId;
          }
        }
      }
      for (OfInt edgeIter = vertexOutMap.edgeIds(v2, v1).intIterator(); edgeIter.hasNext(); ) {
        int edgeId = edgeIter.nextInt();
        for (OfInt inputEdgeIter = graph.inputEdgeIds(edgeId).intIterator();
            inputEdgeIter.hasNext(); ) {
          int inputEdgeId = inputEdgeIter.nextInt();
          if (bestInputEdgeId < 0
              || edgeSites.get(inputEdgeId).size() < edgeSites.get(bestInputEdgeId).size()) {
            bestInputEdgeId = inputEdgeId;
          }
        }
      }
      assert bestInputEdgeId >= 0; // Because there is at least one outgoing edge.

      // Consider the sites near the best input edge snapped to (v1,v2).
      for (OfInt iter = edgeSites.get(bestInputEdgeId).intIterator(); iter.hasNext(); ) {
        // The id and location of the nearby site.
        int vertexId = iter.nextInt();
        S2Point p = graph.vertex(vertexId);
        // Distance from the starting vertex to this nearby site.
        S1ChordAngle r = new S1ChordAngle(p0, p);

        // Sites whose distance from "p0" is at least "r2" are not relevant yet.
        if (r.greaterOrEquals(r2)) {
          continue;
        }

        // The following test prevents us from avoiding previous vertices of the edge chain that
        // also happen to be nearby the current edge. (It also happens to ensure that each vertex is
        // avoided at most once, but this is just an optimization.)
        // TODO(torrey): Would it be faster to do this check first, before considering distance?
        if (!usedVertices.add(vertexId)) {
          continue;
        }

        // We need to figure out whether this site is to the left or right of the edge chain. For
        // the first edge this is easy. Otherwise, since we are only considering sites in the radius
        // range (r1, r2), we can do this by checking whether the site is to the left of the wedge
        // (p0, p1, p2).
        boolean discOnLeft = (v1 == v0)
            ? (S2Predicates.sign(p1, p2, p) > 0)
            : S2Predicates.orderedCCW(p0, p2, p, p1);
        if (!simplifier.avoidDisc(p, builder.minEdgeSiteSeparationChordAngle, discOnLeft)) {
          // This site cannot be avoided given previous constraints.
          return false;
        }
      }
      return true;
    }

    /**
     * Clears the IntVectors in the provided 'mergedInputIds' list, makes sure there are at least
     * 'entries' of them, and each has at least 'capacity'.
     */
    private void clearMergedInputIds(List<IntVector> mergedInputIds, int entries, int capacity) {
      for (IntVector vec : mergedInputIds) {
        vec.ensureCapacity(capacity);
        vec.clear();
      }
      for (int j = mergedInputIds.size(); j < entries; j++) {
        mergedInputIds.add(IntVector.ofCapacity(capacity));
      }
    }

    /**
     * Given the vertices in a simplified edge chain, adds the corresponding simplified edge(s) to
     * the output. The simplified edges are between the first and last vertices in the chain. Note
     * that (1) the edge chain may exist in multiple layers, (2) the edge chain may exist in both
     * directions, (3) there may be more than one copy of an edge chain (in either direction)
     * within a single layer.
     */
    private void mergeChain(IntVector vertices) {
      // Suppose that all interior vertices have M outgoing edges and N incoming edges. Our goal is
      // to group the edges into M outgoing chains and N incoming chains, and then replace each
      // chain by a single edge.

      // For each edge in the chain, a list of the input edge ids snapped and simplified to it.
      final List<IntVector> mergedInputIds = new ArrayList<>();
      // The input edge ids snapped and simplified to all the degenerate edges along interior
      // vertices of the chain.
      final IntVector degenerateIds = new IntVector();

      // Consider each consecutive pair of vertices along the chain.
      for (int i = 1; i < vertices.size(); ++i) {
        int v0 = vertices.get(i - 1);
        int v1 = vertices.get(i);
        // Find all the graph edges in each direction between the two vertices in the chain.
        S2BuilderGraph.VertexOutEdgeIds v0v1Edges = vertexOutMap.edgeIds(v0, v1);
        S2BuilderGraph.VertexOutEdgeIds v1v0Edges = vertexOutMap.edgeIds(v1, v0);

        if (i == 1) {
          // One-time initialization at the start of the loop: Allocate space to store the input
          // edge ids associated with each edge.
          int numEdges = v0v1Edges.size() + v1v0Edges.size();
          // Clear the IntVectors in mergedInputIds and ensure there is one for each edge.
          clearMergedInputIds(mergedInputIds, numEdges, vertices.size() - 1);
        } else {
          // For each interior vertex, we build a list of input edge ids associated with degenerate
          // edges. Each input edge ids will be assigned to one of the output edges later. (Normally
          // there are no degenerate edges at all since most layer types don't want them.)
          assert isInterior.get(v0) == TRUE;
          // Find all loops at v0 and add the input edge ids associated with them to degenerateIds.
          vertexOutMap.edgeIds(v0, v0).forEach(edgeId -> {
            graph.inputEdgeIds(edgeId).forEach(degenerateIds::add);
            used.set(edgeId, TRUE);
          });
        }

        // Because the edges were created in layer order, and all sorts used are stable, the edges
        // are still in layer order. Therefore, we can simply merge together all the edges in the
        // same relative position.
        int j = 0; // For tracking an invariant: The number of edges seen between v0 and v1.

        // For all graph edges from v0 to v1 in order, gather the input edge ids and add them to the
        // mergedInputIds for that graph edge.
        for (OfInt edgeIter = v0v1Edges.intIterator(); edgeIter.hasNext(); ) {
          int edgeId = edgeIter.nextInt();
          IntVector merged = mergedInputIds.get(j);
          for (OfInt idIter = graph.inputEdgeIds(edgeId).intIterator(); idIter.hasNext(); ) {
            merged.add(idIter.nextInt());
          }
          used.set(edgeId, TRUE);
          ++j;
        }

        // Similarly for graph edges from v1 to v0.
        for (OfInt edgeIter = v1v0Edges.intIterator(); edgeIter.hasNext(); ) {
          int edgeId = edgeIter.nextInt();
          IntVector merged = mergedInputIds.get(j);
          for (OfInt idIter = graph.inputEdgeIds(edgeId).intIterator(); idIter.hasNext(); ) {
            merged.add(idIter.nextInt());
          }
          used.set(edgeId, TRUE);
          ++j;
        }

        // Invariant: there are merged ids for each graph edge seen.
        assert mergedInputIds.size() == j;
      }

      // All input edge ids for edges along the chain have been gathered into mergedInputIds.

      // If there were input edge ids assigned to degenerate edges, assign those input edge ids to
      // output edges.
      if (!degenerateIds.isEmpty()) {
        degenerateIds.sort();
        assignDegenerateEdges(degenerateIds, mergedInputIds);
      }

      // Output the merged edges.
      int v0 = vertices.get(0); // First vertex in the chain
      int v1 = vertices.get(1); // Second vertex in the chain
      int vb = vertices.peek(); // Last vertex in the chain.

      // For every non-simplified edge from v0 to v1, output a simplified edge from v0 to vb in the
      // corresponding layer.
      vertexOutMap.edgeIds(v0, v1).forEach(edgeId -> {
        newEdges.add(v0, vb);
        newEdgeLayers.add(graphEdgeLayer(edgeId));
      });
      // And the same from v1 to v0.
      vertexOutMap.edgeIds(v1, v0).forEach(edgeId -> {
        newEdges.add(vb, v0);
        newEdgeLayers.add(graphEdgeLayer(edgeId));
      });

      // Add the sets of merged input edge ids to the lexicon, and map the new edges to them.
      for (IntVector ids : mergedInputIds) {
        newInputEdgeIdSetIds.add(inputEdgeIdSetLexicon.add(ids));
      }
    }

    /**
     * Given "degenerateIds", a list of the input edge ids associated with degenerate edges in the
     * interior of an edge chain, assigns each such input edge id to one of the output edges.
     *
     * <p>"mergedIds" maps output edge id to the input edge ids mapped to that output edge.
     */
    private void assignDegenerateEdges(IntList degenerateIds, List<IntVector> mergedIds) {
      // Each degenerate edge is assigned to an output edge in the appropriate layer. If there is
      // more than one candidate, we use heuristics so that if the input consists of a chain of
      // edges provided in consecutive order (some of which became degenerate), then all those input
      // edges are assigned to the same output edge. For example, suppose that one output edge is
      // labeled with input edges 3,4,7,8, while another output edge is labeled with input edges
      // 50,51,54,57. Then if we encounter degenerate edges labeled with input edges 5 and 6, we
      // would prefer to assign them to the first edge (yielding the continuous range 3,4,5,6,7,8).
      //
      // The heuristic below is only smart enough to handle the case where the candidate output
      // edges have non-overlapping ranges of input edges. (Otherwise there is probably not a good
      // heuristic for assigning the degenerate edges in any case.)

      // Duplicate edge ids don't affect the heuristic below, so we don't bother removing them.
      // (They will be removed by IdSetLexicon.add().)
      // Sort each list of input edge ids.
      for (IntVector ids : mergedIds) {
        ids.sort();
      }

      // Sort the output edges by their minimum input edge id. This is sufficient for the purpose of
      // determining which layer they belong to. With EdgeType.UNDIRECTED, some edges might not have
      // any input edge ids (i.e., if they consist entirely of siblings of input edges). We simply
      // remove such edges from the lists of candidates.

      // Gather the output edge ids that have at least one input edge assigned to them in 'order'.
      IntVector order = IntVector.ofCapacity(mergedIds.size());
      for (int i = 0; i < mergedIds.size(); ++i) {
        if (!mergedIds.get(i).isEmpty()) {
          order.add(i);
        }
      }

      // Sort by the minimum input edge id mapped to that edge.
      order.sort(
          (int i, int j) ->  Integer.compare(mergedIds.get(i).get(0), mergedIds.get(j).get(0)));

      // Now determine where each degenerate edge should be assigned.
      for (OfInt degenerateEdgeIter = degenerateIds.intIterator(); degenerateEdgeIter.hasNext(); ) {
        int degenerateEdgeId = degenerateEdgeIter.nextInt();
        int layer = inputEdgeLayer(degenerateEdgeId);
        // Find the first output edge whose range of input edge ids starts after "degenerateEdgeId".
        // If the *previous* edge belongs to the correct layer, then we assign the degenerate edge
        // that edge.
        int index = order.upperBound(degenerateEdgeId,
            // x is degenerateEdgeId, an input edge id. y is an output edge id, and
            // mergedIds.get(y).get(0) is the minimum input edge id assigned to y.
            (int x, int y) -> Integer.compare(x, mergedIds.get(y).get(0)));

        if (index != 0) {
          if (mergedIds.get(order.get(index - 1)).get(0) >= layerBegins.get(layer)) {
            index--;
          }
        }
        int edgeId = order.get(index);
        assert layer == inputEdgeLayer(mergedIds.get(edgeId).get(0));
        mergedIds.get(edgeId).add(degenerateEdgeId);
      }
    }

    /**
     * A helper class for determining whether a vertex can be an interior vertex of a simplified
     * edge chain. Such a vertex must be adjacent to exactly two vertices (across all layers
     * combined), and in each layer the number of incoming edges from one vertex must equal the
     * number of outgoing edges to the other vertex (in both directions). Furthermore the vertex
     * cannot have any degenerate edges in a given layer unless it has at least one non-degenerate
     * edge in that layer as well. (Note that usually there will not be any degenerate edges at
     * all, since most layer types discard them.)
     *
     * <p>The last condition is necessary to prevent the following: suppose that one layer has a
     * chain ABC and another layer has a degenerate edge BB (with no other edges incident to B).
     * Then we can't simplify ABC to AC because there would be no suitable replacement for the edge
     * BB (since the input edge that mapped to BB can't be replaced by any of the edges AA, AC, or
     * CC without moving further than snapRadius).
     */
    private static class InteriorVertexMatcher {
      /** The id of the vertex that we are checking. */
      private final int v0;
      /** The vertex id of the first neighbor of v0 encountered. */
      private int v1;
      /** The vertex id of the second neighbor of v0 encountered. */
      private int v2;
      /** The number of degenerate edges at v0. */
      private int n0;
      /** The number of edges connecting v0 to v1 */
      private int n1;
      /** Number of edges connecting v0 to v2. */
      private int n2;
      /** outdegree(v0) - indegree(v0) */
      private int excessOut;
      /** Have we seen more than two vertices neighboring v0? */
      private boolean tooManyEndpoints;

      /**
       * Constructs an InteriorVertexMatcher to check if the vertex with vertex id 'v0' can be an
       * interior vertex of an edge chain.
       */
      public InteriorVertexMatcher(int v0) {
        this.v0 = v0;
        v1 = -1;
        v2 = -2;
        n0 = 0;
        n1 = 0;
        n2 = 0;
        excessOut = 0;
        tooManyEndpoints = false;
      }

      /** Starts analyzing the edges of a new layer. */
      public void startLayer() {
        excessOut = 0;
        n0 = 0;
        n1 = 0;
        n2 = 0;
      }

      /**
       * This method should be called for each edge incident to "v0" in a given layer. (For
       * degenerate edges, it should be called twice.)
       */
      public void tally(int v, boolean outgoing) {
        excessOut += outgoing ? 1 : -1;  // outdegree - indegree
        if (v == v0) {
          ++n0;  // Counts both endpoints of each degenerate edge.
        } else {
          // We keep track of the total number of edges (incoming or outgoing) connecting v0 to up
          // to two adjacent vertices.
          if (v1 < 0) {
            v1 = v;
          }
          if (v1 == v) {
            ++n1;
          } else {
            if (v2 < 0) {
              v2 = v;
            }
            if (v2 == v) {
              ++n2;
            } else {
              tooManyEndpoints = true;
            }
          }
        }
      }

      /**
       * This method should be called after processing the edges for each layer. It returns true if
       * "v0" is an interior vertex based on the edges so far.
       */
      public boolean matches() {
        // We check that there are the same number of incoming and outgoing edges in each direction
        // by verifying that (1) indegree(v0) == outdegree(v0) and (2) the total number of edges
        // (incoming and outgoing) to "v1" and "v2" are equal.  We also check the condition on
        // degenerate edges that is documented above.
        return !tooManyEndpoints && excessOut == 0 && n1 == n2 && (n0 == 0 || n1 > 0);
      }
    }
  }

  /**
   * A LayerEdgeIdList is conceptually a List of pairs of LayerId and EdgeId, for every edge in
   * every layer of the graph. The reason for it to exist is that it supports sorting the edges in
   * the order required by mergeLayerEdges(). The pairs are sorted lexicographically, i.e. first by
   * source vertex id, then by destination vertex id, then layer id and lastly edge id.
   */
  private static class LayerEdgeIdList implements Sorter.SortableCollection {
    /** layerEdges is provided to the constructor and maps layerId to an EdgeList of layer edges. */
    private final ArrayList<EdgeList> layerEdges;
    /** layerIds contains the "layerId" part of the (layerId, edgeId) pairs. */
    private final IntVector layerIds = new IntVector();
    /** edgeIds contains the "edgeId" part of the (layerId, edgeId) pairs. */
    private final IntVector edgeIds = new IntVector();

    /**
     * Constructs a new LayerEdgeIdList that uses the given layerEdges to look up source and
     * destination vertex ids.
     */
    public LayerEdgeIdList(ArrayList<EdgeList> layerEdges) {
      this.layerEdges = layerEdges;
    }

    /**
     * Adds the given layerId and edgeId to this list, assigning it the next available sequential
     * index.
     */
    public void add(int layerId, int edgeId) {
      layerIds.add(layerId);
      edgeIds.add(edgeId);
    }

    /** Returns the number of (layer id, edge id) pairs in this LayerEdgeIdList. */
    @Override
    public int size() {
      return layerIds.size();
    }

    @Override
    public void truncate(int start) {
      layerIds.truncate(start);
      edgeIds.truncate(start);
    }

    /** Returns the layer id of the (layer id, edge id) pair at the given index. */
    public int getLayerId(int index) {
      return layerIds.get(index);
    }

    /** Returns the edge id of the (layer id, edge id) pair at the given index. */
    public int getEdgeId(int index) {
      return edgeIds.get(index);
    }

    /** Returns the id of the source vertex of the edge at the given index. */
    public int getSrcId(int index) {
      return layerEdges.get(layerIds.get(index)).getSrcId(edgeIds.get(index));
    }

    /** Returns the id of the destination vertex of the edge at the given index. */
    public int getDstId(int index) {
      return layerEdges.get(layerIds.get(index)).getDstId(edgeIds.get(index));
    }

    /** Returns true if the (layerId, edgeId) pair at indexA precedes the pair at indexB. */
    @Override
    public boolean less(int indexA, int indexB) {
      int layerA = layerIds.get(indexA);
      int layerB = layerIds.get(indexB);
      int edgeIdA = edgeIds.get(indexA);
      int edgeIdB = edgeIds.get(indexB);
      EdgeList layerEdgesA = layerEdges.get(layerA);
      EdgeList layerEdgesB = layerEdges.get(layerB);
      int srcA = layerEdgesA.getSrcId(edgeIdA);
      int srcB = layerEdgesB.getSrcId(edgeIdB);
      // First compare source vertex ids.
      if (srcA < srcB) {
        return true;
      } else if (srcB < srcA) {
        return false;
      }

      // Next compare destination vertex ids.
      int dstA = layerEdgesA.getDstId(edgeIdA);
      int dstB = layerEdgesB.getDstId(edgeIdB);
      if (dstA < dstB) {
        return true;
      } else if (dstB < dstA) {
        return false;
      }

      // Tiebreaking between duplicate edges order must be stable. Compare the layer edge ids pairs
      // themselves: layer id and then edge id.
      if (layerA < layerB) {
        return true;
      } else if (layerB < layerA) {
        return false;
      }
      return edgeIdA < edgeIdB;
    }

    @Override
    public void swap(int indexA, int indexB) {
      layerIds.swap(indexA, indexB);
      edgeIds.swap(indexA, indexB);
    }
  }

  /**
   * EdgeType indicates whether the input edges are undirected. Typically this is specified for each
   * output layer (e.g., {@link S2PolygonLayer}).
   *
   * <p>Directed edges are preferred, since otherwise the output is ambiguous. For example, output
   * polygons may be the *inverse* of the intended result (e.g., a polygon intended to represent the
   * world's oceans may instead represent the world's land masses). Directed edges are also somewhat
   * more efficient.
   *
   * <p>However even with undirected edges, most S2Builder layer types try to preserve the input
   * edge direction whenever possible. Generally, edges are reversed only when it would yield a
   * simpler output. For example, S2PolygonLayer assumes that polygons created from undirected edges
   * should cover at most half of the sphere. Similarly, S2PolylineVectorLayer assembles edges into
   * as few polylines as possible, even if this means reversing some of the "undirected" input
   * edges.
   *
   * <p>For shapes with interiors, directed edges should be oriented so that the interior is to the
   * left of all edges. This means that for a polygon with holes, the outer loops ("shells") should
   * be directed counter-clockwise while the inner loops ("holes") should be directed clockwise.
   * Note that S2Builder.addPolygon() follows this convention automatically.
   */
  public enum EdgeType {
    DIRECTED,
    UNDIRECTED
  }

  /**
   * A SnapFunction restricts the locations of the output vertices. For example, there are
   * predefined snap functions that require vertices to be located at S2CellId centers or at
   * E5/E6/E7 coordinates. The SnapFunction can also specify a minimum spacing between vertices (the
   * "snap radius"). SnapFunction implementations should be immutable. S2Builder and other users
   * of SnapFunctions may assume they are immutable, and have undefined behavior if they are changed
   * while being used.
   *
   * <p>A SnapFunction defines the following methods:
   *
   * <ol>
   *   <li>The snapPoint() method, which snaps a point P to a nearby point (the "candidate snap
   *       site"). Any point may be returned, including P itself (this is the "identity snap
   *       function").
   *   <li>snapRadius(), the maximum distance that vertices can move when snapped. The snapRadius
   *       must be at least as large as the maximum distance between P and snapPoint(P) for any
   *       point P.
   *       <p>Note that the maximum distance that edge interiors can move when snapped is slightly
   *       larger than "snapRadius", and is returned by the function {@link
   *       S2Builder.Builder#maxEdgeDeviation()} (see there for details).
   *   <li>minVertexSeparation(), the guaranteed minimum distance between vertices in the output.
   *       This is generally a fraction of "snapRadius" where the fraction depends on the snap
   *       function.
   *   <li>minEdgeVertexSeparation(), the guaranteed minimum distance between edges and non-incident
   *       vertices in the output. This is generally a fraction of "snapRadius" where the fraction
   *       depends on the snap function.
   * </ol>
   *
   * <p>It is important to note that snapPoint() does not define the actual mapping from input
   * vertices to output vertices, since the points it returns (the candidate snap sites) are further
   * filtered to ensure that they are separated by at least the snap radius. For example, if you
   * specify E7 coordinates (2cm resolution) and a snap radius of 10m, then a subset of points
   * returned by SnapPoint will be chosen (the "snap sites"), and each input vertex will be mapped
   * to the closest site. Therefore you cannot assume that P is necessarily snapped to snapPoint(P).
   *
   * <p>S2Builder makes the following guarantees:
   *
   * <ol>
   *   <li>Every vertex is at a location returned by snapPoint().
   *   <li>Vertices are within "snapRadius" of the corresponding input vertex.
   *   <li>Edges are within "maxEdgeDeviation" of the corresponding input edge (a distance slightly
   *       larger than "snapRadius").
   *   <li>Vertices are separated by at least "minVertexSeparation" (a fraction of "snapRadius" that
   *       depends on the snap function).
   *   <li>Edges and non-incident vertices are separated by at least "minEdgeVertexSeparation" (a
   *       fraction of "snapRadius").
   *   <li>Vertex and edge locations do not change unless one of the conditions above is not already
   *       met (idempotency / stability).
   *   <li>The topology of the input geometry is preserved (up to the creation of degeneracies).
   *       This means that there exists a continuous deformation from the input to the output such
   *       that no vertex crosses an edge.
   * </ol>
   */
  public interface SnapFunction {
    /**
     * The maximum distance that vertices can move when snapped. The snap radius can be any value
     * between zero and {@link S2BuilderSnapFunctions#maxSnapRadius()}.
     *
     * <p>If the snap radius is zero, then vertices are snapped together only if they are identical.
     * Edges will not be snapped to any vertices other than their endpoints, even if there are
     * vertices whose distance to the edge is zero, unless splitCrossingEdges() is true (see below).
     *
     * <p>REQUIRES: {@code snapRadius() <= maxSnapRadius()}.
     */
    public abstract S1Angle snapRadius();

    /**
     * The guaranteed minimum distance between vertices in the output. This is generally some
     * fraction of "snapRadius".
     */
    public abstract S1Angle minVertexSeparation();

    /**
     * The guaranteed minimum spacing between edges and non-incident vertices in the output. This is
     * generally some fraction of "snapRadius".
     */
    public abstract S1Angle minEdgeVertexSeparation();

    /**
     * Returns a candidate snap site for the given point. The final vertex locations are a subset of
     * the snap sites returned by this function (spaced at least "minVertexSeparation" apart).
     *
     * <p>The only requirement is that snapPoint(x) must return a point whose distance from "x" is
     * no greater than "snapRadius".
     */
    public abstract S2Point snapPoint(S2Point point);
  }

  /**
   * The S2Builder.Builder is a mutable container for S2Builder's options, and allows for
   * constructing new S2Builders with those options. Once an S2Builder has been constructed, its
   * options are immutable. A new S2Builder.Builder can be obtained from an S2Builder with a copy of
   * its options.
   *
   * <p>The default options are documented with their getter and setter methods. TODO(torrey):
   * consider moving those here.
   */
  public static class Builder {
    private SnapFunction snapFunction =
        new S2BuilderSnapFunctions.IdentitySnapFunction(S1Angle.ZERO);
    private boolean splitCrossingEdges = false;
    private S1Angle intersectionTolerance = S1Angle.ZERO;
    private boolean simplifyEdgeChains = false;
    private boolean idempotent = true;

    /**
     * Constructs a new Builder with default options. In short, those are to use an
     * identitySnapFunction with snapRadius 0, do not split crossing edges or simplify edge chains,
     * and intersection tolerance of zero.
     */
    public Builder() {}

    /** Copy constructor. */
    public Builder(Builder options) {
      setSnapFunction(options.snapFunction);
      setIdempotent(options.idempotent);
      setSimplifyEdgeChains(options.simplifyEdgeChains);
      setSplitCrossingEdges(options.splitCrossingEdges);
      setIntersectionTolerance(options.intersectionTolerance);
    }

    /** Convenience constructor that calls setSnapFunction(). */
    public Builder(SnapFunction snapFunction) {
      setSnapFunction(snapFunction);
    }

    /** Constructs and returns a new S2Builder using the options in this Builder. */
    public S2Builder build() {
      return new S2Builder(this);
    }

    /**
     * Sets the desired snap function. Note that if your input data includes vertices that were
     * created using {@link S2EdgeUtil#getIntersection(S2Point, S2Point, S2Point, S2Point)}, then
     * you should use a "snapRadius" of at least {@link S2EdgeUtil#INTERSECTION_MERGE_RADIUS}, e.g.
     * by calling
     *
     * {@snippet :
     * options.setSnapFunction(new IdentitySnapFunction(S2.INTERSECTION_MERGE_RADIUS));
     * }
     *
     * <p>DEFAULT: {@code S2BuilderSnapFunctions.IdentitySnapFunction(S1Angle.ZERO)}. This does no
     * snapping and preserves all input vertices exactly.
     */
    @CanIgnoreReturnValue
    public Builder setSnapFunction(SnapFunction snapFunction) {
      this.snapFunction = snapFunction;
      return this;
    }

    /** Returns the current SnapFunction. See {@link setSnapFunction(SnapFunction)}. */
    public SnapFunction snapFunction() {
      return snapFunction;
    }

    /**
     * The maximum distance from snapped edge vertices to the original edge. This is
     * snapFunction().snapRadius() plus {@link #intersectionTolerance()}.
     *
     * <p>When the splitCrossingEdges() option is true (see below), the value of
     * intersectionTolerance() is the larger of {@link S2EdgeUtil#INTERSECTION_ERROR} and the value
     * from the Builder options, which defaults to zero.
     */
    public S1Angle edgeSnapRadius() {
      return snapFunction().snapRadius().add(intersectionTolerance());
    }

    /**
     * The maximum distance that any point along an edge can move when snapped. It is slightly
     * larger than edgeSnapRadius() because when a geodesic edge is snapped, the edge center moves
     * further than its endpoints. S2Builder ensures that this distance is at most 10% larger than
     * edgeSnapRadius().
     */
    public S1Angle maxEdgeDeviation() {
      // We want maxEdgeDeviation() to be large enough compared to snapRadius() that edge splitting
      // is rare.
      //
      // Using spherical trigonometry, if the endpoints of an edge of length L move by at most a
      // distance R, the center of the edge moves by at most asin(sin(R) / cos(L / 2)). Thus the
      // (maxEdgeDeviation / snapRadius) ratio increases with both the snap radius R and the edge
      // length L.
      //
      // We arbitrarily limit the edge deviation to be at most 10% more than the snap radius. With
      // the maximum allowed snap radius of 70 degrees, this means that edges up to 30.6 degrees
      // long are never split. For smaller snap radii, edges up to 49 degrees long are never split.
      // (Edges of any length are not split unless their endpoints move far enough so that the
      // actual edge deviation exceeds the limit; in practice, splitting is rare even with long
      // edges.) Note that it is always possible to split edges when maxEdgeDeviation() is exceeded;
      // see maybeAddExtraSites().
      Preconditions.checkState(
          snapFunction().snapRadius().lessOrEquals(S2BuilderSnapFunctions.maxSnapRadius()));
      final double kMaxEdgeDeviationRatio = 1.1;
      return edgeSnapRadius().mul(kMaxEdgeDeviationRatio);
    }

    /**
     * If true, then detect all pairs of crossing edges and eliminate them by adding a new vertex at
     * their intersection point. See also the {@link #addIntersection(S2Point)} method which allows
     * intersection points to be added selectively.
     *
     * <p>When this option if true, {@link #intersectionTolerance()} is automatically set to a
     * minimum of {@link S2EdgeUtil#INTERSECTION_ERROR} (see {@link #intersectionTolerance()} for
     * why this is necessary). Note that this means that edges can move by up to {@link
     * S2EdgeUtil#INTERSECTION_ERROR} even when the specified snap radius is zero. The exact
     * distance that edges can move is always given by {@link #maxEdgeDeviation()} defined above.
     *
     * <p>Undirected edges should always be used when the output is a polygon, since splitting a
     * directed loop at a self-intersection converts it into two loops that don't define a
     * consistent interior according to the "interior is on the left" rule. (On the other hand, it
     * is fine to use directed edges when defining a polygon *mesh* because in that case the input
     * consists of sibling edge pairs.)
     *
     * <p>Self-intersections can also arise when importing data from a 2D projection. You can
     * minimize this problem by subdividing the input edges so that the S2 edges (which are
     * geodesics) stay close to the original projected edges (which are curves on the sphere). This
     * can be done using S2EdgeTessellator, for example.
     *
     * <p>DEFAULT: false
     */
    public boolean splitCrossingEdges() {
      return splitCrossingEdges;
    }

    /** Sets the splitCrossingEdges option. See {@link #splitCrossingEdges()}. */
    @CanIgnoreReturnValue
    public Builder setSplitCrossingEdges(boolean splitCrossingEdges) {
      this.splitCrossingEdges = splitCrossingEdges;
      return this;
    }

    /**
     * Specifies the maximum allowable distance between a vertex added by addIntersection() and the
     * edge(s) that it is intended to snap to. The intersection tolerance must be set before
     * addIntersection() can be used. It has the effect of increasing the snap radius for edges (but
     * not vertices) by the given distance.
     *
     * <p>The intersection tolerance should be set to the maximum error in the intersection
     * calculation used. For example, if {@link S2EdgeUtil#getIntersection(S2Point, S2Point,
     * S2Point, S2Point)} is used then the error should be set to {@link
     * S2EdgeUtil#INTERSECTION_ERROR}. If {@link S2EdgeUtil#getPointOnLine(S2Point, S2Point,
     * S2Point)} is used then the error should be set to {@link S2EdgeUtil#GET_POINT_ON_LINE_ERROR}.
     * If {@link S2EdgeUtil#project(S2Point, S2Point, S2Point)} is used then the error should be set
     * to {@link S2EdgeUtil#GET_POINT_ON_RAY_PERPENDICULAR_ERROR}. If more than one method is used
     * then the intersection tolerance should be set to the maximum such error.
     *
     * <p>The reason this option is necessary is that computed intersection points are not exact.
     * For example, {@link S2EdgeUtil#getIntersection(a, b, c, d)} returns a point up to {@link
     * S2EdgeUtil#INTERSECTION_ERROR} away from the true mathematical intersection of the edges AB
     * and CD. Furthermore such intersection points are subject to further snapping in order to
     * ensure that no pair of vertices is closer than the specified snap radius. For example,
     * suppose the computed intersection point X of edges AB and CD is 1 nanonmeter away from both
     * edges, and the snap radius is 1 meter. In that case X might snap to another vertex Y exactly
     * 1 meter away, which would leave us with a vertex Y that could be up to 1.000000001 meters
     * from the edges AB and/or CD. This means that AB and/or CD might not snap to Y leaving us with
     * two edges that still cross each other.
     *
     * <p>However if the intersection tolerance is set to 1 nanometer then the snap radius for edges
     * is increased to 1.000000001 meters ensuring that both edges snap to a common vertex even in
     * this worst case. (Tthis technique does not work if the vertex snap radius is increased as
     * well; it requires edges and vertices to be handled differently.)
     *
     * <p>Note that this option allows edges to move by up to the given intersection tolerance even
     * when the snap radius is zero. The exact distance that edges can move is always given by
     * maxEdgeDeviation() defined above.
     *
     * <p>When splitCrossingEdges() is true, the intersection tolerance is automatically set to a
     * minimum of {@link S2EdgeUtil#INTERSECTION_ERROR}. A larger value can be specified by calling
     * this method explicitly.
     *
     * <p>DEFAULT: {@link S1Angle#ZERO}
     */
    public S1Angle intersectionTolerance() {
      if (!splitCrossingEdges()) {
        return intersectionTolerance;
      }
      return S1Angle.max(intersectionTolerance, S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR));
    }

    /** Sets the intersectionTolerance option. See {@link intersectionTolerance()}. */
    @CanIgnoreReturnValue
    public Builder setIntersectionTolerance(S1Angle intersectionTolerance) {
      Preconditions.checkArgument(intersectionTolerance.greaterThan(S1Angle.ZERO));
      this.intersectionTolerance = intersectionTolerance;
      return this;
    }

    /**
     * If true, then simplify the output geometry by replacing nearly straight chains of short edges
     * with a single long edge.
     *
     * <p>The combined effect of snapping and simplifying will not change the input by more than the
     * guaranteed tolerances (see the list documented with the SnapFunction class). For example,
     * simplified edges are guaranteed to pass within snapRadius() of the *original* positions of
     * all vertices that were removed from that edge. This is a much tighter guarantee than can be
     * achieved by snapping and simplifying separately.
     *
     * <p>However, note that this option does not guarantee idempotency. In other words, simplifying
     * geometry that has already been simplified once may simplify it further. (This is unavoidable,
     * since tolerances are measured with respect to the original geometry, which is no longer
     * available when the geometry is simplified a second time.)
     *
     * <p>When the output consists of multiple layers, simplification is guaranteed to be
     * consistent: for example, edge chains are simplified in the same way across layers, and
     * simplification preserves topological relationships between layers (e.g., no crossing edges
     * will be created). Note that edge chains in different layers do not need to be identical (or
     * even have the same number of vertices, etc) in order to be simplified together. All that is
     * required is that they are close enough together so that the same simplified edge can meet all
     * of their individual snapping guarantees.
     *
     * <p>Note that edge chains are approximated as parametric curves rather than point sets. This
     * means that if an edge chain backtracks on itself (for example, ABCDEFEDCDEFGH) then such
     * backtracking will be preserved to within snapRadius() (for example, if the preceding point
     * were all in a straight line then the edge chain would be simplified to ACFCFH, noting that C
     * and F have degree > 2 and therefore can't be simplified away).
     *
     * <p>Simplified edges are assigned all labels associated with the edges of the simplified
     * chain.
     *
     * <p>For this option to have any effect, a SnapFunction with a non-zero snapRadius() must be
     * specified. Also note that vertices specified using forceVertex are never simplified away.
     *
     * <p>DEFAULT: false
     */
    public boolean simplifyEdgeChains() {
      return simplifyEdgeChains;
    }

    /** Sets the simplifyEdgeChains option. See {@link #simplifyEdgeChains()}. */
    @CanIgnoreReturnValue
    public Builder setSimplifyEdgeChains(boolean simplifyEdgeChains) {
      this.simplifyEdgeChains = simplifyEdgeChains;

      // Simplification requires a non-zero snap radius, and while it might be possible to do some
      // simplifying without snapping, it is much simpler to always snap (even if the input
      // geometry already meets the other output requirements). We need to compute edgeSites in
      // order to avoid approaching non-incident vertices too closely, for example.
      setIdempotent(false);
      return this;
    }

    /**
     * If true, then snapping occurs only when the input geometry does not already meet the
     * S2Builder output guarantees (see the SnapFunction class description for details). This means
     * that if all input vertices are at snapped locations, all vertex pairs are separated by at
     * least minVertexSeparation(), and all edge-vertex pairs are separated by at least
     * minEdgeVertexSeparation(), then no snapping is done.
     *
     * <p>If false, then all vertex pairs and edge-vertex pairs closer than "snapRadius" will be
     * considered for snapping. This can be useful, for example, if you know that your geometry
     * contains errors and you want to make sure that features closer together than "snapRadius" are
     * merged.
     *
     * <p>This option is automatically turned off by {@link #simplifyEdgeChains()}, since
     * simplifying edge chains is never guaranteed to be idempotent.
     *
     * <p>DEFAULT: true
     */
    public boolean idempotent() {
      return idempotent;
    }

    /** Sets the idempotent option. See {@link #idempotent()}. */
    @CanIgnoreReturnValue
    public Builder setIdempotent(boolean idempotent) {
      this.idempotent = idempotent;
      return this;
    }
  }

  /**
   * The GraphOptions class is only needed by S2BuilderLayer implementations. A layer is responsible
   * for assembling an S2BuilderGraph of snapped edges into the desired output format (e.g., an
   * S2Polygon). The GraphOptions class allows each Layer type to specify requirements on its input
   * graph: for example, if DegenerateEdges.DISCARD is specified, then S2Builder will ensure that
   * all degenerate edges are removed before passing the graph to the S2Layer.build() method.
   */
  public static class GraphOptions {
    private EdgeType edgeType;
    private DegenerateEdges degenerateEdges;
    private DuplicateEdges duplicateEdges;
    private SiblingPairs siblingPairs;
    private boolean allowVertexFiltering;

    /**
     * All S2BuilderLayer subtypes should specify GraphOptions explicitly using this constructor,
     * rather than relying on default values.
     */
    public GraphOptions(
        EdgeType edgeType,
        DegenerateEdges degenerateEdges,
        DuplicateEdges duplicateEdges,
        SiblingPairs siblingPairs) {
      this.edgeType = edgeType;
      this.degenerateEdges = degenerateEdges;
      this.duplicateEdges = duplicateEdges;
      this.siblingPairs = siblingPairs;
      this.allowVertexFiltering = true;
    }

    /**
     * The default GraphOptions constructor specifies that all edges should be kept, since this
     * produces the least surprising output and makes it easier to diagnose the problem when an
     * option is left unspecified.
     */
    public GraphOptions() {
      this.edgeType = EdgeType.DIRECTED;
      this.degenerateEdges = DegenerateEdges.KEEP;
      this.duplicateEdges = DuplicateEdges.KEEP;
      this.siblingPairs = SiblingPairs.KEEP;
      this.allowVertexFiltering = true;
    }

    /** A copy constructor for GraphOptions. */
    public GraphOptions(GraphOptions other) {
      this.edgeType = other.edgeType();
      this.degenerateEdges = other.degenerateEdges();
      this.duplicateEdges = other.duplicateEdges();
      this.siblingPairs = other.siblingPairs();
      this.allowVertexFiltering = other.allowVertexFiltering();
    }

    @Override
    public String toString() {
      return Strings.lenientFormat(
          "EdgeType %s, DegenerateEdges %s, DuplicateEdges %s, SiblingPairs %s,"
              + " AllowVertexFiltering %s",
          edgeType, degenerateEdges, duplicateEdges, siblingPairs, allowVertexFiltering);
    }

    /**
     * Specifies whether the S2Builder input edges should be treated as undirected. If true, then
     * all input edges are duplicated into pairs consisting of an edge and a sibling (reverse) edge.
     * Note that the automatically created sibling edge has an empty set of labels and does not have
     * an associated inputEdgeId.
     *
     * <p>The layer implementation is responsible for ensuring that exactly one edge from each pair
     * is used in the output, i.e. *only half* of the graph edges will be used. (Note that some
     * values of the siblingPairs() option automatically take care of this issue by removing half of
     * the edges and changing edgeType() to DIRECTED.)
     *
     * <p>DEFAULT: EdgeType.DIRECTED
     */
    public EdgeType edgeType() {
      return edgeType;
    }

    /** Sets the edgeType for this GraphOptions. See {@link edgeType()}. */
    public void setEdgeType(EdgeType edgeType) {
      this.edgeType = edgeType;
    }

    /**
     * Controls how degenerate edges (i.e., an edge from a vertex to itself) are handled. Such edges
     * may be present in the input, or they may be created when both endpoints of an edge are
     * snapped to the same output vertex.
     *
     * <p>The options available are:
     *
     * <p>DISCARD: Discards all degenerate edges. This is useful for layers that do not support
     * degeneracies, such as S2PolygonLayer.
     *
     * <p>DISCARD_EXCESS: Discards all degenerate edges that are connected to non-degenerate edges
     * and merges any remaining duplicate degenerate edges. This is useful for simplifying polygons
     * while ensuring that loops that collapse to a single point do not disappear.
     *
     * <p>KEEP: Keeps all degenerate edges. Be aware that this may create many redundant edges when
     * simplifying geometry (e.g., a polyline of the form AABBBBBCCCCCCDDDD). DegenerateEdges.KEEP
     * is mainly useful for algorithms that require an output edge for every input edge.
     *
     * <p>DEFAULT: DegenerateEdges.KEEP
     */
    public enum DegenerateEdges {
      DISCARD,
      DISCARD_EXCESS,
      KEEP,
    };

    /** Returns the DegenerateEdges option for this GraphOptions. See {@link DegenerateEdges}. */
    public DegenerateEdges degenerateEdges() {
      return degenerateEdges;
    }

    /** Sets the DegenerateEdges option for this GraphOptions. See {@link DegenerateEdges}. */
    public void setDegenerateEdges(DegenerateEdges degenerateEdges) {
      this.degenerateEdges = degenerateEdges;
    }

    /**
     * Controls how duplicate edges (i.e., edges that are present multiple times) are handled. Such
     * edges may be present in the input, or they can be created when vertices are snapped together.
     * When several edges are merged, the result is a single edge labelled with all of the original
     * input edge ids.
     *
     * <p>DEFAULT: DuplicateEdges.KEEP
     */
    public enum DuplicateEdges {
      MERGE,
      KEEP,
    };

    /** Returns the DuplicateEdges option for this GraphOptions. See {@link DuplicateEdges}. */
    public DuplicateEdges duplicateEdges() {
      return duplicateEdges;
    }

    /** Sets the DuplicateEdges option for this GraphOptions. See {@link DuplicateEdges}. */
    public void setDuplicateEdges(DuplicateEdges duplicateEdges) {
      this.duplicateEdges = duplicateEdges;
    }

    /**
     * Controls how sibling edge pairs (i.e., pairs consisting of an edge and its reverse edge) are
     * handled. Layer types that define an interior (e.g., polygons) normally discard such edge
     * pairs since they do not affect the result (i.e., they define a "loop" with no interior). The
     * various options include:
     *
     * <ul>
     *   <li>DISCARD: Discards all sibling edge pairs.
     *   <li>DISCARD_EXCESS: Like DISCARD, except that a single sibling pair is kept if the result
     *       would otherwise be empty. This is useful for polygons with degeneracies
     *       (S2LaxPolygonShape), and for simplifying polylines while ensuring that they are not
     *       split into multiple disconnected pieces.
     *   <li>KEEP: Keeps sibling pairs. This can be used to create polylines that double back on
     *       themselves, or degenerate loops (with a layer type such as S2LaxPolygonShape).
     *   <li>REQUIRE: Requires that all edges have a sibling (and returns an error otherwise). This
     *       is useful with layer types that create a collection of adjacent polygons (a polygon
     *       mesh).
     *   <li>CREATE: Ensures that all edges have a sibling edge by creating them if necessary. This
     *       is useful with polygon meshes where the input polygons do not cover the entire sphere.
     *       Such edges always have an empty set of labels and do not have an associated
     *       InputEdgeId.
     * </ul>
     *
     * <p>If edgeType() is EdgeType.UNDIRECTED, a sibling edge pair is considered to consist of four
     * edges (two duplicate edges and their siblings), since only two of these four edges will be
     * used in the final output.
     *
     * <p>Furthermore, since the options REQUIRE and CREATE guarantee that all edges will have
     * siblings, S2Builder implements these options for undirected edges by discarding half of the
     * edges in each direction and changing the edgeType() to EdgeType::DIRECTED. For example, two
     * undirected input edges between vertices A and B would first be converted into two directed
     * edges in each direction, and then one edge of each pair would be discarded leaving only one
     * edge in each direction.
     *
     * <p>Degenerate edges are considered not to have siblings. If such edges are present, they are
     * passed through unchanged by SiblingPairs.DISCARD. For SiblingPairs.REQUIRE or
     * SiblingPairs.CREATE with undirected edges, the number of copies of each degenerate edge is
     * reduced by a factor of two.
     *
     * <p>Any of the options that discard edges (DISCARD, DISCARD_EXCESS, and REQUIRE/CREATE in the
     * case of undirected edges) have the side effect that when duplicate edges are present, all of
     * the corresponding edge labels are merged together and assigned to the remaining edges. (This
     * avoids the problem of having to decide which edges are discarded.) Note that this merging
     * takes place even when all copies of an edge are kept. For example, consider the graph {AB1,
     * AB2, AB3, BA4, CD5, CD6} (where XYn denotes an edge from X to Y with label "n"). With
     * SiblingPairs.DISCARD, we need to discard one of the copies of AB. But which one? Rather than
     * choosing arbitrarily, instead we merge the labels of all duplicate edges (even ones where no
     * sibling pairs were discarded), yielding {AB123, AB123, CD45, CD45} (assuming that duplicate
     * edges are being kept).
     *
     * <p>Notice that the labels of duplicate edges are merged even if no siblings were discarded
     * (such as CD5, CD6 in this example), and that this would happen even with duplicate degenerate
     * edges (e.g. the edges EE7, EE8).
     *
     * <p>DEFAULT: SiblingPairs.KEEP
     */
    public enum SiblingPairs {
      DISCARD,
      DISCARD_EXCESS,
      KEEP,
      REQUIRE,
      CREATE,
    };

    /** Returns the SiblingPairs option for this GraphOptions. See {@link SiblingPairs}. */
    public SiblingPairs siblingPairs() {
      return siblingPairs;
    }

    /** Sets the SiblingPairs option for this GraphOptions. See {@link SiblingPairs}. */
    public void setSiblingPairs(SiblingPairs siblingPairs) {
      this.siblingPairs = siblingPairs;
      // TODO(torrey): Implement and test the documented interactions between options. For instance,
      // setting certain values of siblingPairs changes edgeType.
    }

    /**
     * This is a specialized option that is only needed by clients that want to work with the graphs
     * for multiple layers at the same time (e.g., in order to check whether the same edge is
     * present in two different graphs). Note that if you need to do this, usually it is easier just
     * to build a single graph with suitable edge labels.
     *
     * <p>When there are a large number of layers, then by default S2Builder builds a minimal
     * subgraph for each layer containing only the vertices needed by the edges in that layer. This
     * ensures that layer types that iterate over the vertices run in time proportional to the size
     * of that layer rather than the size of all layers combined. (For example, if there are a
     * million layers with one edge each, then each layer would be passed a graph with 2 vertices
     * rather than 2 million vertices.)
     *
     * <p>If this option is set to false, this optimization is disabled. Instead the graph passed to
     * this layer will contain the full set of vertices. (This is not recommended when the number of
     * layers could be large.)
     *
     * <p>DEFAULT: true
     */
    public boolean allowVertexFiltering() {
      return allowVertexFiltering;
    }

    public void setAllowVertexFiltering(boolean allowVertexFiltering) {
      this.allowVertexFiltering = allowVertexFiltering;
    }

    @Override
    public boolean equals(Object other) {
      if (!(other instanceof GraphOptions)) {
        return false;
      }
      GraphOptions otherGraphOptions = (GraphOptions) other;
      return this.edgeType == otherGraphOptions.edgeType
          && this.degenerateEdges == otherGraphOptions.degenerateEdges
          && this.duplicateEdges == otherGraphOptions.duplicateEdges
          && this.siblingPairs == otherGraphOptions.siblingPairs
          && this.allowVertexFiltering == otherGraphOptions.allowVertexFiltering;
    }

    @Override
    public int hashCode() {
      return Objects.hash(
          edgeType, degenerateEdges, duplicateEdges, siblingPairs, allowVertexFiltering);
    }
  }

  /**
   * For output layers that represent polygons, there is an ambiguity inherent in spherical geometry
   * that does not exist in planar geometry. Namely, if a polygon has no edges, does it represent
   * the empty polygon (containing no points) or the full polygon (containing all points)? This
   * ambiguity also occurs for polygons that consist only of degeneracies, e.g. a degenerate loop
   * with only two edges could be either a degenerate shell in the empty polygon or a degenerate
   * hole in the full polygon.
   *
   * <p>To resolve this ambiguity, an IsFullPolygonPredicate may be specified for each output layer
   * (see addIsFullPolygonPredicate below). If the output after snapping consists only of degenerate
   * edges and/or sibling pairs (including the case where there are no edges at all), then the layer
   * implementation calls the given predicate to determine whether the polygon is empty or full
   * except for those degeneracies. The predicate is given an S2BuilderGraph containing the output
   * edges, but note that in general the predicate must also have knowledge of the input geometry in
   * order to determine the correct result.
   *
   * <p>This predicate is only needed by layers that are assembled into polygons. It is not used by
   * other layer types.
   */
  public interface IsFullPolygonPredicate {
    public boolean test(S2BuilderGraph graph);
  }
}
