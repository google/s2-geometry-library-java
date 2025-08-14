/*
 * Copyright 2023 Google Inc.
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

import static com.google.common.base.Preconditions.checkNotNull;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2Builder.SnapFunction;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderGraph.VertexInMap;
import com.google.common.geometry.S2BuilderGraph.VertexOutMap;
import com.google.common.geometry.S2CrossingEdgesQuery.CrossingType;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntPairVector;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * This class implements boolean operations (intersection, union, difference, and symmetric
 * difference) for regions whose boundaries are defined by geodesic edges.
 *
 * <p>S2BooleanOperation operates on exactly two input regions at a time. Each region is represented
 * as an S2ShapeIndex and may contain any number of points, polylines, and polygons. The region is
 * essentially the union of these objects, except that polygon interiors must be disjoint from all
 * other geometry (including other polygon interiors). If the input geometry for a region does not
 * meet this condition, it can be normalized by computing its union first. Duplicate polygon edges
 * are not allowed (even among different polygons), however polylines may have duplicate edges and
 * may even be self-intersecting. Note that points or polylines are allowed to coincide with the
 * boundaries of polygons.
 *
 * <p>Degeneracies are fully supported. Supported degeneracy types include the following:
 *
 * <ul>
 *   <li>Point polylines consisting of a single degenerate edge AA.
 *   <li>Point loops consisting of a single vertex A. Such loops may represent either shells or
 *       holes according to whether the loop adds to or subtracts from the surrounding region of the
 *       polygon.
 *   <li>Sibling edge pairs of the form {AB, BA}. Such sibling pairs may represent either shells or
 *       holes according to whether they add to or subtract from the surrounding region. The edges
 *       of a sibling pair may belong to the same polygon loop (e.g. a loop AB) or to different
 *       polygon loops or polygons (e.g. the polygons {ABC, CBD}).
 * </ul>
 *
 * <p>A big advantage of degeneracy support is that geometry may be simplified without completely
 * losing small details. For example, if a polygon representing a land area with many lakes and
 * rivers is simplified using a tolerance of 1 kilometer, every water feature in the input is
 * guaranteed to be within 1 kilometer of some water feature in the output (even if some lakes and
 * rivers are merged and/or reduced to degenerate point or sibling edge pair holes). Mathematically
 * speaking, degeneracy support allows geometry to be simplified while guaranteeing that the
 * Hausdorff distance between the boundaries of the original and simplified geometries is at most
 * the simplification tolerance. It also allows geometry to be simplified without changing its
 * dimension, thus preserving boundary semantics. (Note that the boundary of a polyline ABCD is
 * {A,D}, whereas the boundary of a degenerate shell ABCDCB is its entire set of vertices and
 * edges.)
 *
 * <p>Points and polyline edges are treated as multisets: if the same point or polyline edge appears
 * multiple times in the input, it will appear multiple times in the output. For example, the union
 * of a point with an identical point consists of two points. This feature is useful for modeling
 * large sets of points or polylines as a single region while maintaining their distinct identities,
 * even when the points or polylines intersect each other. It is also useful for reconstructing
 * polylines that loop back on themselves (e.g., time series such as GPS tracks). If duplicate
 * geometry is not desired, it can easily be removed by choosing the appropriate S2Builder output
 * layer options.
 *
 * <p>Self-intersecting polylines can be manipulated without materializing new vertices at the
 * self-intersection points. This feature is important when processing polylines with large numbers
 * of self-intersections such as GPS tracks (e.g., consider the path of a race car in the Indy 500).
 *
 * <p>Polylines are always considered to be directed. Polyline edges between the same pair of
 * vertices are defined to intersect even if the two edges are in opposite directions. (Undirected
 * polylines can be modeled by specifying GraphOptions.EdgeType.UNDIRECTED in the S2Builder output
 * layer.)
 *
 * <p>The output of each operation is sent to an S2BuilderLayer provided by the client. This allows
 * clients to build any representation of the geometry they choose. It also allows the client to do
 * additional postprocessing of the output before building data structures; for example, the client
 * can easily discard degeneracies or convert them to another data type.
 *
 * <p>The boundaries of polygons and polylines can be modeled as open, semi-open, or closed.
 * Polyline boundaries are controlled by the PolylineModel class, whose options are as follows:
 *
 * <ul>
 *   <li>In the OPEN model, polylines do not contain their first or last vertex except for one
 *       special case: namely, if the polyline forms a loop and the polylineLoopsHaveBoundaries()
 *       option is set to false, then the first/last vertex is contained.
 *   <li>In the SEMI_OPEN model, polylines contain all vertices except the last. Therefore if one
 *       polyline starts where another polyline stops, the two polylines do not intersect.
 *   <li>In the CLOSED model, polylines contain all of their vertices.
 * </ul>
 *
 * <p>When multiple polylines are present, they are processed independently and have no effect on
 * each other. For example, in the OPEN boundary model the polyline ABC contains the vertex B, while
 * set of polylines {AB, BC} does not. (If you want to treat the polylines as a union instead, with
 * boundaries merged according to the "mod 2" rule, this can be achieved by reassembling the edges
 * into maximal polylines using S2PolylineVectorLayer with EdgeType.UNDIRECTED,
 * DuplicateEdges.MERGE, and PolylineType.WALK.)
 *
 * <p>Polygon boundaries are controlled by the PolygonModel class, which has the following options:
 *
 * <ul>
 *   <li>In the OPEN model, polygons do not contain their vertices or edges. This implies that a
 *       polyline that follows the boundary of a polygon will not intersect it.
 *   <li>In the SEMI_OPEN model, polygon point containment is defined such that if several polygons
 *       tile the region around a vertex, then exactly one of those polygons contains that vertex.
 *       Similarly polygons contain all of their edges, but none of their reversed edges. This
 *       implies that a polyline and polygon edge with the same endpoints intersect if and only if
 *       they are in the same direction. (This rule ensures that if a polyline is intersected with a
 *       polygon and its complement, the two resulting polylines do not have any edges in common.)
 *   <li>In the CLOSED model, polygons contain all their vertices, edges, and reversed edges. This
 *       implies that a polyline that shares an edge (in either direction) with a polygon is defined
 *       to intersect it. Similarly, this is the only model where polygons that touch at a vertex or
 *       along an edge intersect.
 * </ul>
 *
 * <p>Note that PolylineModel and PolygonModel are defined as separate classes in order to allow for
 * possible future extensions.
 *
 * <p>Operations between geometry of different dimensions are defined as follows:
 *
 * <ul>
 *   <li>For UNION, the higher-dimensional shape always wins. For example the union of a closed
 *       polygon A with a polyline B that coincides with the boundary of A consists only of the
 *       polygon A.
 *   <li>For INTERSECTION, the lower-dimensional shape always wins. For example, the intersection of
 *       a closed polygon A with a point B that coincides with a vertex of A consists only of the
 *       point B.
 *   <li>For DIFFERENCE, higher-dimensional shapes are not affected by subtracting lower-dimensional
 *       shapes. For example, subtracting a point or polyline from a polygon A yields the original
 *       polygon A. This rule exists because in general, it is impossible to represent the output
 *       using the specified boundary model(s). (Consider subtracting one vertex from a
 *       PolylineModel.CLOSED polyline, or subtracting one edge from a PolygonModel.CLOSED polygon.)
 *       If you want to perform operations like this, consider representing all boundaries
 *       explicitly (topological boundaries) using OPEN boundary models. Another option for polygons
 *       is to subtract a degenerate loop, which yields a polygon with a degenerate hole (see
 *       S2LaxPolygonShape).
 * </ul>
 *
 * <p>Note that in the case of Precision.EXACT operations, the above remarks only apply to the
 * output before snapping. Snapping may cause nearby distinct edges to become coincident, e.g. a
 * polyline may become coincident with a polygon boundary. However also note that S2BooleanOperation
 * is perfectly happy to accept such geometry as input.
 *
 * <p>Example usage:
 *
 * {@snippet :
 * S2ShapeIndex a, b;   // Input geometry, e.g. containing polygons.
 * S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder();
 * builder.setSnapFunction(snapFunction);
 * S2PolygonLayer layer = new S2PolygonLayer();
 * S2BooleanOperation op = builder.build(OpType.INTERSECTION, layer, options);
 * S2Error error = new S2Error;
 * if (!op.build(a, b, error)) {
 *   logger.atError(error);
 *   return ...;
 * }
 * S2Polygon polygon = layer.getPolygon();   // Output geometry.
 * ...
 * }
 *
 * <p>If the output includes objects of different dimensions, they can be assembled into different
 * layers with code like this:
 *
 * {@snippet :
 * S2BooleanOperation op = builder.build(
 *     OpType.UNION,
 *     new S2PointVectorLayer(),
 *     new S2PolylineVectorLayer(),
 *     new S2PolygonLayer());
 * }
 *
 * <p>Boolean operations are implemented by constructing the boundary of the result and then using
 * S2Builder to assemble the edges. The boundary is obtained by clipping each of the two input
 * regions to the interior or exterior of the other region. For example, to compute the union of A
 * and B, we clip the boundary of A to the exterior of B and the boundary of B to the exterior of A;
 * the resulting set of edges defines the union of the two regions.
 *
 * <p>We use exact predicates, but inexact constructions (e.g. computing the intersection point of
 * two edges). Nevertheless, the following algorithm is guaranteed to be 100% robust, in that the
 * computed boundary stays within a small tolerance (snapRadius + S2EdgeUtil.INTERSECTION_ERROR) of
 * the exact result, and also preserves the correct topology (i.e., no crossing edges).
 *
 * <p>Unfortunately this robustness cannot quite be achieved using the strategy outlined above
 * (clipping the two input regions and assembling the resulting edges). Since computed intersection
 * points are not exact, the input geometry passed to S2Builder might contain self-intersections,
 * and these self-intersections cannot be eliminated reliably by snap rounding.
 *
 * <p>So instead, we pass S2Builder the entire set of input edges where at least some portion of
 * each edge belongs to the output boundary. We allow S2Builder to compute the intersection points
 * and snap round the edges (which it does in a way that is guaranteed to preserve the input
 * topology). Then once this is finished, we remove the portions of each edge that would have been
 * clipped if we had done the clipping first. This step only involves deciding whether to keep or
 * discard each edge in the output, since all intersection points have already been resolved, and
 * therefore there is no risk of creating new self-intersections.
 */
@SuppressWarnings("Assertion")
public class S2BooleanOperation {

  /** A bit mask representing all six faces of the S2 cube. */
  private static final byte ALL_FACES_MASK = 0x3f;

  /** The Options in use for this boolean operation. */
  private final Options options;

  /** The type of boolean operation being performed. */
  private final OpType opType;

  /** The two input regions. */
  private final S2ShapeIndex[] regions = new S2ShapeIndex[2];

  /**
   * The output of S2BooleanOperation consists either of zero layers (in which case resultEmpty is
   * used), one layer for all the output geometry, or three layers for the output geometry by
   * dimension.
   */
  private final List<S2BuilderLayer> layers;

  /**
   * If there are zero output layers, resultEmpty will be set to indicate if the output is empty.
   * Otherwise it is unused.
   */
  private boolean resultEmpty;

  // A collection of special InputEdgeIds that allow the GraphEdgeClipper state modifications to be
  // inserted into the list of edge crossings.

  /** If true, specifies that edges are retained. Otherwise, "outside", edges are discarded. */
  private static final int SET_INSIDE = -1;

  /** If true, specifies that edges should be clipped to the exterior of the other region. */
  private static final int SET_INVERT_B = -2;

  /** If true, specifies that edges should be reversed before emitting them. */
  private static final int SET_REVERSE_A = -3;

  /** The supported operation types. */
  public enum OpType {
    UNION, // Contained by either region.
    INTERSECTION, // Contained by both regions.
    DIFFERENCE, // Contained by the first region but not the second.
    SYMMETRIC_DIFFERENCE // Contained by one region but not the other.
  }

  /**
   * Defines whether polygons are considered to contain their vertices and/or edges (see definitions
   * above).
   */
  public enum PolygonModel {
    OPEN,
    SEMI_OPEN,
    CLOSED
  }

  /**
   * Defines whether polylines are considered to contain their endpoints (see definitions above).
   */
  public enum PolylineModel {
    OPEN,
    SEMI_OPEN,
    CLOSED
  }

  /**
   * With Precision.EXACT, the operation is evaluated using the exact input geometry. Predicates
   * that use this option will produce exact results; for example, they can distinguish between a
   * polyline that barely intersects a polygon from one that barely misses it. Constructive
   * operations (ones that yield new geometry, as opposed to predicates) are implemented by
   * computing the exact result and then snap rounding it according to the given snapFunction() (see
   * below). This is as close as it is possible to get to the exact result while requiring that
   * vertex coordinates have type "double".
   *
   * <p>Note that Precision.SNAPPED is not yet implemented. When it is implemented, the behavior
   * will be that the input regions are snapped together *before* the operation is evaluated. So for
   * example, two polygons that overlap slightly will be treated as though they share a common
   * boundary, and similarly two polygons that are slightly separated from each other will be
   * treated as though they share a common boundary. Snapped results are useful for dealing with
   * points, since in S2 the only points that lie exactly on a polyline or polygon edge are the
   * endpoints of that edge.
   *
   * <p>Conceptually, the difference between these two options is that with Precision.SNAPPED, the
   * inputs are snap rounded (together), whereas with Precision.EXACT only the result is snap
   * rounded.
   */
  public enum Precision {
    EXACT,
    SNAPPED // Not yet implemented. See above.
  }

  /**
   * A SourceId with a non-negative edgeId identifies an edge from one of the two input
   * S2ShapeIndexes. A SourceId consists of a region id (0 or 1), a shape id within that region's
   * S2ShapeIndex, and an edge id within that shape. SourceIds don't actually exist as objects, but
   * are encoded and decoded as long values.
   *
   * <p>SourceIds with negative edgeIds are "special", and used to encode state changes for the
   * GraphEdgeClipper.
   */
  private static class SourceId {
    private static final long REGION_MASK = 1L << 63;
    private static final long SHAPE_MASK = (REGION_MASK - 1) << 32;
    private static final long EDGE_MASK = (1L << 32) - 1;

    /** Encodes a SourceId consisting of the given regionId, shapeId, and edgeId into a long. */
    private static long encode(int regionId, int shapeId, int edgeId) {
      Preconditions.checkArgument(edgeId >= 0);
      return (((long) regionId) << 63) | (((long) shapeId) << 32) | edgeId;
    }

    /** Returns an encoded SourceId with a special edgeIds that control the GraphEdgeClipper. */
    public static long special(int specialEdgeId) {
      Preconditions.checkArgument(specialEdgeId < 0);
      return specialEdgeId; // encode(0, 0, specialEdgeId);
    }
  }

  /**
   * Options for S2BooleanOperation.
   */
  public static final class Options {
    private SnapFunction snapFunction;
    private PolygonModel polygonModel;
    private PolylineModel polylineModel;
    private boolean polylineLoopsHaveBoundaries;
    private boolean splitAllCrossingPolylineEdges;
    // Currently, Precision cannot be changed from the default.
    private final Precision precision;
    // Currently, conservativeOutput cannot be changed from the default.
    private final boolean conservativeOutput;

    /** The default SnapFunction. See {@link #snapFunction()}. */
    public static final SnapFunction DEFAULT_SNAP_FUNCTION =
        new S2BuilderSnapFunctions.IdentitySnapFunction();

    /** The default PolygonModel. See {@link #polygonModel()}. */
    public static final PolygonModel DEFAULT_POLYGON_MODEL = PolygonModel.SEMI_OPEN;

    /** The default PolylineModel. See {@link #polylineModel()}. */
    public static final PolylineModel DEFAULT_POLYLINE_MODEL = PolylineModel.CLOSED;

    /** The default value of {@link #polylineLoopsHaveBoundaries()}. */
    public static final boolean DEFAULT_POLYLINE_LOOPS_HAVE_BOUNDARIES = true;

    /** The default value of {@link #splitAllCrossingPolylineEdges()}. */
    public static final boolean DEFAULT_SPLIT_ALL_CROSSING_POLYLINE_EDGES = false;

    /** The default value of {@link #precision()}. */
    public static final Precision DEFAULT_PRECISION = Precision.EXACT;

    /** The default value of {@link #conservativeOutput()}. */
    public static final boolean DEFAULT_CONSERVATIVE_OUTPUT = false;

    /** An immutable instance of Options with the default values. */
    public static final Options DEFAULT = new Options();

    /** The Options constructor initializes the default options. */
    public Options() {
      this.snapFunction = DEFAULT_SNAP_FUNCTION;
      this.polygonModel = DEFAULT_POLYGON_MODEL;
      this.polylineModel = DEFAULT_POLYLINE_MODEL;
      this.polylineLoopsHaveBoundaries = DEFAULT_POLYLINE_LOOPS_HAVE_BOUNDARIES;
      this.splitAllCrossingPolylineEdges = DEFAULT_SPLIT_ALL_CROSSING_POLYLINE_EDGES;
      this.precision = DEFAULT_PRECISION;
      this.conservativeOutput = DEFAULT_CONSERVATIVE_OUTPUT;
    }

    /** A copy constructor for internal use. */
    Options(Options other) {
      this.snapFunction = other.snapFunction;
      this.polygonModel = other.polygonModel;
      this.polylineModel = other.polylineModel;
      this.polylineLoopsHaveBoundaries = other.polylineLoopsHaveBoundaries;
      this.splitAllCrossingPolylineEdges = other.splitAllCrossingPolylineEdges;
      this.precision = other.precision;
      this.conservativeOutput = other.conservativeOutput;
    }

    /**
     * Returns the function to be used for snap rounding. The default is {@link
     * S2BuilderSnapFunctions.IdentitySnapFunction} with a snapRadius of zero. This does no snapping
     * and preserves all input vertices exactly unless there are crossing edges, in which case the
     * snap radius is increased to the maximum intersection point error,
     * {@link S2EdgeUtil#INTERSECTION_ERROR}.
     */
    public SnapFunction snapFunction() {
      return snapFunction;
    }

    /**
     * Returns whether polygons are considered to contain their vertices and/or edges (see comments
     * above). The default is PolygonModel.SEMI_OPEN.
     */
    public PolygonModel polygonModel() {
      return polygonModel;
    }

    /**
     * Returns whether polylines are considered to contain their vertices (see comments above). The
     * default is PolylineModel.CLOSED.
     */
    public PolylineModel polylineModel() {
      return polylineModel;
    }

    /**
     * Returns whether a polyline loop is considered to have a non-empty boundary. By default this
     * option is true, meaning that even if the first and last vertices of a polyline are the same,
     * the polyline is considered to have a well-defined "start" and "end". For example, if the
     * polyline boundary model is OPEN then the polyline loop would not include the start/end
     * vertices. These are the best semantics for most applications, such as GPS tracks or road
     * network segments.
     *
     * <p>If the polyline forms a loop and this option is set to false, then instead the first and
     * last vertices are considered to represent a single vertex in the interior of the polyline. In
     * this case the boundary of the polyline is empty, meaning that the first/last vertex will be
     * contained by the polyline even if the boundary model is OPEN. (Note that this option also has
     * a small effect on the CLOSED boundary model, because the first/last vertices of a polyline
     * loop are considered to represent one vertex rather than two.)
     *
     * <p>The main reason for this option is to implement the "mod 2 union" boundary semantics of
     * the OpenGIS Simple Features spec. This can be achieved by making sure that all polylines are
     * constructed using S2Builder.Graph.PolylineType.WALK (which ensures that all polylines are as
     * long as possible), and then setting this option to false.
     */
    public boolean polylineLoopsHaveBoundaries() {
      return polylineLoopsHaveBoundaries;
    }

    /**
     * Specifies that a new vertex should be added whenever a polyline edge crosses another polyline
     * edge. Note that this can cause the size of polylines with many self-intersections to increase
     * quadratically.
     *
     * <p>If false, new vertices are added only when a polyline from one input region cross a
     * polyline from the other input region. This allows self-intersecting input polylines to be
     * modified as little as possible.
     */
    public boolean splitAllCrossingPolylineEdges() {
      return splitAllCrossingPolylineEdges;
    }

    /**
     * Specifies whether the operation should use the exact input geometry (Precision.EXACT), or
     * whether the two input regions should be snapped together first (Precision.SNAPPED).
     */
    public Precision precision() {
      return precision;
    }

    /**
     * If true, the input geometry is interpreted as representing nearby geometry that has been
     * snapped or simplified. It then outputs a conservative result based on the value of
     * polygonModel() and polylineModel(). For the most part, this only affects the handling of
     * degeneracies.
     *
     * <p>- If the model is OPEN, the result is as open as possible. For example, the intersection
     * of two identical degenerate shells is empty under PolygonModel.OPEN because they could have
     * been disjoint before snapping. Similarly, two identical degenerate polylines have an empty
     * intersection under PolylineModel.OPEN.
     *
     * <p>- If the model is CLOSED, the result is as closed as possible. In the case of the
     * DIFFERENCE operation, this is equivalent to evaluating A - B as Closure(A) - Interior(B). For
     * other operations, it affects only the handling of degeneracies. For example, the union of two
     * identical degenerate holes is empty under PolygonModel.CLOSED (i.e., the hole disappears)
     * because the holes could have been disjoint before snapping.
     *
     * <p>- If the model is SEMI_OPEN, the result is as degenerate as possible. New degeneracies
     * will not be created, but all degeneracies that coincide with the opposite region's boundary
     * are retained unless this would cause a duplicate polygon edge to be created. This model is is
     * very useful for working with input data that has both positive and negative degeneracies
     * (i.e., degenerate shells and holes).
     */
    public boolean conservativeOutput() {
      return conservativeOutput;
    }

    /**
     * If specified, then each output edge will be labelled with one or more SourceIds indicating
     * which input edge(s) it corresponds to. This can be useful if your input geometry has
     * additional data that needs to be propagated from the input to the output (e.g., elevations).
     *
     * <p>You can access the labels by using an S2BuilderLayer type that supports labels, such as
     * S2PolygonLayer. The layer outputs a "labelSetLexicon" and an "labelSetId" for each edge. You
     * can then look up the source information for each edge like this:
     *
     * {@snippet :
     * for (int label : labelSetLexicon.idSet(labelSetId)) {
     * SourceId src = sourceIdLexicon.value(label);
     * // regionId() specifies which S2ShapeIndex the edge is from (0 or 1).
     * doSomething(src.regionId(), src.shapeId(), src.edgeId());
     * }
     * }
     */
    // ValueLexicon<SourceId> sourceIdLexicon();  TODO(user): Implement

    /** For debugging. */
    @Override
    public String toString() {
      return "S2BooleanOperation.Options{snapFunction="
          + snapFunction()
          + ", polygonModel="
          + polygonModel()
          + ", polylineModel="
          + polylineModel()
          + ", polylineLoopsHaveBoundaries="
          + polylineLoopsHaveBoundaries()
          + ", splitAllCrossingPolylineEdges="
          + splitAllCrossingPolylineEdges()
          + "}";
    }
  }

  /**
   * The S2BooleanOperation.Builder provides a mutable interface to S2BooleanOperation's Options,
   * and allows for constructing new S2BooleanOperations with those options.
   */
  public static class Builder  {
    /** These options are mutated by the Builder. */
    private Options options = new Options();

    /** Constructs a Builder with default {@link Options} values. */
    public Builder() {
      options.snapFunction = Options.DEFAULT_SNAP_FUNCTION;
    }

    /**
     * Convenience constructor that uses the given snapFunction and otherwise uses default Options.
     */
    public Builder(SnapFunction snapFunction) {
      options.snapFunction = snapFunction;
    }

    /** Constructor that accepts provided Options. */
    public Builder(Options options) {
      this.options = options;
    }

    /** Returns a snapshot of the Options currently set on this Builder. */
    public Options options() {
      return new Options(options);
    }

    /**
     * Using a snapshot of the Options currently set for this Builder, constructs and returns an
     * S2BooleanOperation for which the output boundary edges will be sent to the three specified
     * layers according to their dimension. Points (represented by degenerate edges) are sent to
     * layer0, polyline edges are sent to layer1, and polygon edges are sent to layer2.
     *
     * <p>The dimension of an edge is defined as the minimum dimension of the two input edges that
     * produced it. For example, the intersection of two crossing polyline edges is a considered to
     * be a degenerate polyline rather than a point, so it is sent to layer 1. Clients can easily
     * reclassify such polylines as points if desired, but this rule makes it easier for clients
     * that want to process point, polyline, and polygon inputs differently.
     *
     * <p>The layers are always built in the order 0, 1, 2, and all arguments to the layer.build()
     * calls are guaranteed to be valid until the last call returns. All Graph objects have the same
     * set of vertices and the same lexicon objects, in order to make it easier to write classes
     * that process all the edges in parallel.
     */
    public S2BooleanOperation build(
        OpType opType,
        S2BuilderLayer layer0,
        S2BuilderLayer layer1,
        S2BuilderLayer layer2) {
      return new S2BooleanOperation(opType, ImmutableList.of(layer0, layer1, layer2), options());
    }

    /**
     * Using a snapshot of the Options currently defined set for Builder, constructs and returns an
     * S2BooleanOperation for which the output boundary edges should be sent to the single provided
     * S2Builder layer. This version can be used when the dimension of the output geometry is known
     * (e.g., intersecting two polygons to yield a third polygon).
     */
    public S2BooleanOperation build(OpType opType, S2BuilderLayer layer) {
      return new S2BooleanOperation(opType, ImmutableList.of(layer), options());
    }

    /** For debugging. */
    @Override
    public String toString() {
      return "S2BooleanOperation.Builder{snapFunction="
          + options.snapFunction()
          + ", polygonModel="
          + options.polygonModel()
          + ", polylineModel="
          + options.polylineModel()
          + ", polylineLoopsHaveBoundaries="
          + options.polylineLoopsHaveBoundaries()
          + ", splitAllCrossingPolylineEdges="
          + options.splitAllCrossingPolylineEdges()
          + "}";
    }

    /**
     * Specifies the function to be used for snap rounding. See the comments above for
     * {@link Options#snapFunction()}.
     */
    @CanIgnoreReturnValue
    public Builder setSnapFunction(SnapFunction snapFunction) {
      options.snapFunction = snapFunction;
      return this;
    }

    /**
     * Specifies whether polygons are considered to contain their vertices and/or edges. See the
     * comments above for {@link Options#polygonModel()}.
     */
    @CanIgnoreReturnValue
    public Builder setPolygonModel(PolygonModel model) {
      options.polygonModel = model;
      return this;
    }

    /**
     * Specifies whether polylines are considered to contain their vertices. See the
     * comments above for {@link Options#polylineModel()}.
     */
    @CanIgnoreReturnValue
    public Builder setPolylineModel(PolylineModel model) {
      options.polylineModel = model;
      return this;
    }

    /**
     * Specifies whether a polyline loop is considered to have a non-empty boundary. See the
     * comments above for {@link Options#polylineLoopsHaveBoundaries()}.
     */
    public void setPolylineLoopsHaveBoundaries(boolean value) {
      options.polylineLoopsHaveBoundaries = value;
    }

    /**
     * Specifies if crossing polyline edges should be split by adding a new vertex at each crossing.
     * See the comments above for {@link Options#splitAllCrossingPolylineEdges()}.
     */
    @CanIgnoreReturnValue
    public Builder setSplitAllCrossingPolylineEdges(boolean value) {
      options.splitAllCrossingPolylineEdges = value;
      return this;
    }

    // TODO(user): Precision.SNAPPED is not yet implemented in C++ or Java.
    // public Builder setPrecision(Precision precision) {
    //   options.precision = precision;
    //   return this;
    // }

    // TODO(user): Implement support for conservative output = true.
    // public Builder setConservativeOutput(boolean conservative) {
    //  options.conservativeOutput = conservative;
    //  return this;
    // }

    // TODO(user): Support sourceIdLexicon.
    // public Builder setSourceIdLexicon(ValueLexicon<SourceId> sourceIdLexicon) {
    //  options.sourceIdLexicon = sourceIdLexicon;
    //  return this;
    // }
  }

  /**
   * Internal constructor for an S2BooleanOperation. Uses the provided Options. Either zero, one, or
   * three layers must be provided. If three layers are provided, the S2BooleanOperation will send
   * output boundary edges to the three specified layers according to their dimension. Points
   * (represented by degenerate edges) are sent to layer 0, polyline edges are sent to layer 1, and
   * polygon edges are sent to layer 2.
   *
   * <p>Otherwise, if one layer is provided, the S2BooleanOperation will send all output boundary
   * edges to the single specified layer. If zero layers are provided, the S2BooleanOperation can be
   * used for efficiently testing boolean relationships without producing output boundary edges at
   * all. (See {@link #isEmpty(OpType, S2ShapeIndex, S2ShapeIndex, Options)} below).
   *
   * <p>For the case with three output layers:
   *
   * <p>The dimension of an edge is defined as the minimum dimension of the two input edges that
   * produced it. For example, the intersection of two crossing polyline edges is a considered to be
   * a degenerate polyline rather than a point, so it is sent to layer 1. Clients can easily
   * reclassify such polylines as points if desired, but this rule makes it easier for clients that
   * want to process point, polyline, and polygon inputs differently.
   *
   * <p>The layers are always built in the order 0, 1, 2, and all arguments to the build() calls are
   * guaranteed to be valid until the last call returns. All Graph objects have the same set of
   * vertices and the same lexicon objects, in order to make it easier to write classes that process
   * all the edges in parallel.
   */
  private S2BooleanOperation(OpType opType, List<S2BuilderLayer> layers, Options options) {
    this.options = options;
    this.opType = opType;
    this.layers = layers;
  }

  /** Returns the OpType of this S2BooleanOperation. */
  public OpType opType() {
    return opType;
  }

  /** Returns this {@code S2BooleanOperation}'s immutable Options. */
  public Options options() {
    return options;
  }

  /**
   * Executes the operation specified to the constructor on the provided inputs. Returns true on
   * success, and otherwise sets "error" appropriately.
   *
   * <p>S2BooleanOperation does not set any errors itself, but the S2BuilderLayer or layers may.
   * Generally, errors should not occur for valid inputs, and as of mid-2025, S2BuilderLayer
   * implementations provided in the S2 library, (outside of unit test implementations), only set
   * errors in the following cases:
   *
   * <ul>
   *   <li>S2LaxPolygonLayer doesn't support undirected edges yet, and sets an error if they are
   *       requested.
   *   <li>S2LaxPolylineLayer and S2PolylineLayer set an error if a single polyline was requested,
   *       but the given edges cannot be assembled into a single polyline.
   *   <li>S2PointVectorLayer was given edges that are not points.
   * </ul>
   */
  public boolean build(S2ShapeIndex a, S2ShapeIndex b, S2Error error) {
    regions[0] = checkNotNull(a);
    regions[1] = checkNotNull(b);

    Impl impl = new Impl(this);
    return impl.build(error);
  }

  /**
   * Executes the operation specified to the constructor on the provided inputs. This method wraps
   * {@link #build(S2ShapeIndex, S2ShapeIndex, S2Error)} as a convenience for clients who either do
   * not wish to handle errors, or prefer using exception handling.
   *
   * @throws S2Exception wrapping the underlying {@link S2Error}, if any error occurs.
   */
  public void buildUnsafe(S2ShapeIndex a, S2ShapeIndex b) {
    S2Error error = new S2Error();
    if (!build(a, b, error)) {
      throw new S2Exception(error);
    }
  }

  /**
   * Returns true if the result of the given operation is empty, using the provided indexes and
   * Options. Much faster than actually computing the output.
   */
  public static boolean isEmpty(OpType opType, S2ShapeIndex a, S2ShapeIndex b, Options options) {
    S2BooleanOperation op = new S2BooleanOperation(opType, ImmutableList.of(), options);
    S2Error error = new S2Error();
    boolean ok = op.build(a, b, error);
    assert ok;
    return op.resultEmpty;
  }

  /**
   * Returns true if the result of the given operation is empty, using the provided indexes and
   * default Options. Much faster than actually computing the output.
   */
  public static boolean isEmpty(OpType opType, S2ShapeIndex a, S2ShapeIndex b) {
    return isEmpty(opType, a, b, Options.DEFAULT);
  }

  /** Convenience method that returns true if A intersects B, using the given Options. */
  public static boolean intersects(S2ShapeIndex a, S2ShapeIndex b, Options options) {
    return !isEmpty(OpType.INTERSECTION, b, a, options);
  }

  /** Convenience method that returns true if A intersects B, using default Options. */
  public static boolean intersects(S2ShapeIndex a, S2ShapeIndex b) {
    return !isEmpty(OpType.INTERSECTION, b, a);
  }

  /**
   * Convenience method that returns true if A contains B, i.e., if the difference (B - A) is empty.
   */
  public static boolean contains(S2ShapeIndex a, S2ShapeIndex b, Options options) {
    return isEmpty(OpType.DIFFERENCE, b, a, options);
  }

  /** Convenience method that returns true if A contains B, using default Options. */
  public static boolean contains(S2ShapeIndex a, S2ShapeIndex b) {
    return isEmpty(OpType.DIFFERENCE, b, a);
  }

  /**
   * Convenience method that returns true if the symmetric difference of A and B is empty. (Note
   * that A and B may still not be identical, e.g. A may contain two copies of a polyline while B
   * contains one.)
   */
  public static boolean equals(S2ShapeIndex a, S2ShapeIndex b, Options options) {
    return isEmpty(OpType.SYMMETRIC_DIFFERENCE, b, a, options);
  }

  /**
   * Convenience method that returns true if the symmetric difference of A and B is empty, using
   * default Options. (Note that A and B may still not be identical, e.g. A may contain two copies
   * of a polyline while B contains one.)
   */
  public static boolean equals(S2ShapeIndex a, S2ShapeIndex b) {
    return isEmpty(OpType.SYMMETRIC_DIFFERENCE, b, a);
  }

  /**
   * CrossingInputEdge represents an input edge B that crosses some other input edge A. It stores
   * the input edge id of edge B and also whether it crosses edge A from left to right (or vice
   * versa).
   *
   * <p>TODO: Consider removing this class. CrossingInputEdge objects are used in two
   * places: (1) CrossingInputEdgePair is an int inputEdgeId and a CrossingInputEdge. That could
   * just be two ints (the input edge ids) and a boolean leftToRight, and perhaps encoded with
   * leftToRight as a sign bit. (2) CrossingInputEdgeVector, in GraphEdgeClipper, is just a List of
   * CrossingInputEdge with a lowerBound() operation. That could be an IntVector, encoding the
   * leftToRight as the sign bit.
   */
  private static class CrossingInputEdge {
    private final int inputEdgeId;
    private final boolean leftToRight;

    /**
     * Constructs a CrossingInputEdge for input edge B, crossing the other input edge A leftToRight.
     */
    public CrossingInputEdge(int inputEdgeId, boolean leftToRight) {
      this.inputEdgeId = inputEdgeId;
      this.leftToRight = leftToRight;
    }

    /** The input edge id of the edge B that crosses the other input edge A. */
    public int inputEdgeId() {
      return inputEdgeId;
    }

    /** Whether the input edge B crosses the other input edge A from left to right. */
    public boolean leftToRight() {
      return leftToRight;
    }

    @Override
    public String toString() {
      if (inputEdgeId >= 0) {
        return Platform.formatString("CrossingInputEdge(%d, %s)", inputEdgeId, leftToRight);
      }
      switch (inputEdgeId) {
        case SET_REVERSE_A:
          return Platform.formatString("CrossingInputEdge(SET_REVERSE_A = %s", leftToRight);
        case SET_INVERT_B:
          return Platform.formatString("CrossingInputEdge(SET_INVERT_B = %s", leftToRight);
        case SET_INSIDE:
          return Platform.formatString("CrossingInputEdge(SET_INSIDE = %s", leftToRight);
        default:
          throw new IllegalArgumentException("Unknown edgeId: " + inputEdgeId);
      }
    }
  }

  /**
   * Like a {@code Pair<InputEdgeId, CrossingInputEdge>}.
   *
   * <p>TODO: Consider removing this class. CrossingInputEdgePairs are only used in
   * InputEdgeCrossings, which is just an ArrayList of CrossingInputEdgePairs. This could be
   * replaced with parallel IntVectors (simple) or a single IntVector (a bit more complicated, but
   * better locality).
   */
  private static class CrossingInputEdgePair {
    private final int inputEdgeId;
    private final CrossingInputEdge crossingInputEdge;

    /**
     * Constructs a CrossingInputEdgePair for an input edge id, and some other input edge that
     * crosses it.
     */
    CrossingInputEdgePair(int inputEdgeId, CrossingInputEdge crossingInputEdge) {
      this.inputEdgeId = inputEdgeId;
      this.crossingInputEdge = crossingInputEdge;
    }

    /** The input edge id of the first of the crossing edge pair. */
    public int first() {
      return inputEdgeId;
    }

    /** The input edge id of the second of the crossing edge pair, and the crossing direction. */
    public CrossingInputEdge second() {
      return crossingInputEdge;
    }

    @Override
    public String toString() {
      return Platform.formatString("CrossingInputEdgePair(%d, %s)", inputEdgeId, crossingInputEdge);
    }
  }

  /**
   * InputEdgeCrossings represents all pairs of intersecting input edges and also certain
   * GraphEdgeClipper state modifications (SET_INSIDE, etc). It is sorted lexicographically except
   * for entries representing state modifications, which are sorted by the first InputEdgeId only.
   */
  private static class InputEdgeCrossings extends ArrayList<CrossingInputEdgePair> {}

  /**
   * Given two input edges A and B that intersect, suppose that A maps to a chain of snapped edges
   * A_0, A_1, ..., A_m and B maps to a chain of snapped edges B_0, B_1, ..., B_n. CrossingGraphEdge
   * represents an edge from chain B that shares a vertex with chain A. It is used as a temporary
   * data representation while processing chain A. The contents are:
   *
   * <ol>
   *   <li>"edgeId" - the Graph.EdgeId (int) of an edge from chain B.
   *   <li>"aIndex" - the index of the vertex (A_i) that is shared with chain A.
   *   <li>"outgoing" - true if the shared vertex is the first vertex of the B edge.
   *   <li>"otherVertexId" - the Graph.VertexId (int) of the vertex that is not shared with chain A.
   * </ol>
   *
   * Note that if an edge from the B chain shares both vertices with the A chain, there will be two
   * entries: an outgoing edge that treats its first vertex as being shared, and an incoming edge
   * that treats its second vertex as being shared.
   */
  private static class CrossingGraphEdge {
    final int edgeId;
    final int aIndex;
    final boolean outgoing;
    final int otherVertexId;

    /** Constructs a CrossingGraphEdge with contents as described above. */
    public CrossingGraphEdge(int edgeId, int aIndex, boolean outgoing, int otherVertexId) {
      this.edgeId = edgeId;
      this.aIndex = aIndex;
      this.outgoing = outgoing;
      this.otherVertexId = otherVertexId;
    }

    @Override
    public String toString() {
      return Platform.formatString(
          "CrossingGraphEdge(edgeId=%d, aIndex=%d, outgoing=%b, otherVertexId=%d)",
          edgeId, aIndex, outgoing, otherVertexId);
    }
  }

  // TODO: Reconsider implementation of CrossingGraphEdge and CrossingGraphEdgeVector.
  // C++ uses an InlinedVector, so this is likely to be a small list.
  //
  // CrossingGraphEdgeVector can probably be implemented as a wrapper around a single IntVector,
  // and CrossingGraphEdge can be removed entirely, replaced with an index into the vector.
  //
  // CrossingGraphEdges are created in gatherIncidentEdges and immediately added to one of the
  //   CrossingGraphVectors in the bEdges array of CrossingGraphVectors.
  // In getCrossedVertexIndex, the elements of a CrossingGraphEdgeVector are iterated over.
  // In getVertexRank, the .edgeId and .outgoing fields are used to compute rank of a
  //   CrossingGraphEdge.
  //
  // CrossingGraphEdgeVectors exist only in an array, 'bEdges', created in GraphEdgeClipper.run(),
  //   and reinitialized for each edge as run() loops over the graph edges.
  // The initialization of bEdges[] is done in gatherIncidentEdges().
  // Later in run(), individual CrossingGraphEdgeVectors in bEdges are passed to
  //   getCrossedVertexIndex().
  //

  /**
   * A List of CrossingGraphEdge, which represents all the crossing input edges of an input edge A.
   * Each CrossingGraphEdge represents one edge from a chain B that shares a vertex with chain A
   * (see above).
   */
  private static class CrossingGraphEdgeVector extends ArrayList<CrossingGraphEdge> {
    /** Returns the first CrossingGraphEdge in the vector. */
    public CrossingGraphEdge first() {
      return get(0);
    }

    /** Returns the last CrossingGraphEdge in the vector. */
    public CrossingGraphEdge last() {
      return get(size() - 1);
    }

    /** Returns the 'aIndex' value of the first CrossingGraphEdge in the vector. */
    public int firstAIndex() {
      return get(0).aIndex;
    }

    /** Returns the 'aIndex' value of the last CrossingGraphEdge in the vector. */
    public int lastAIndex() {
      return get(size() - 1).aIndex;
    }

    @Override
    public String toString() {
      StringBuilder sb = new StringBuilder();
      for (int i = 0; i < size(); ++i) {
        sb.append(get(i));
        if (i < size() - 1) {
          sb.append(", ");
        }
      }
      return sb.toString();
    }
  }

  /** A vector of CrossingGraphEdgeVector. */
  private static class CrossingGraphEdgeVectorVector {
    private CrossingGraphEdgeVector[] vectors = new CrossingGraphEdgeVector[32];
    private int size = 0;

    public int size() {
      return size;
    }

    public CrossingGraphEdgeVector get(int index) {
      return vectors[index];
    }

    /**
     * Adds a CrossingGraphEdge constructed from the provided values to the end of the
     * CrossingGraphEdgeVector that is at position 'vectorIndex' in this vector of vectors.
     */
    public void addCrossingGraphEdge(
        int vectorIndex, int edgeId, int aIndex, boolean outgoing, int otherVertexId) {
      vectors[vectorIndex].add(
          new CrossingGraphEdge(edgeId, aIndex, outgoing, otherVertexId));
    }

    /** Clears all the contents and reinitializes to contain 'size' empty CrossingGraphVectors. */
    public void reinitialize(int size) {
      this.size = size;

      // Only reallocate a larger array if required.
      if (vectors.length < size) {
        vectors = Arrays.copyOf(vectors, size);
      }
      // Clear the contents of each vector.
      for (int v = 0; v < size; ++v) {
        if (vectors[v] == null) {
          vectors[v] = new CrossingGraphEdgeVector();
        } else {
          vectors[v].clear();
        }
      }
    }
  }

  /**
   * Returns a vector of EdgeIds sorted by input edge id. When more than one output edge has the
   * same input edge id (i.e., the input edge snapped to a chain of edges), the edges are sorted so
   * that they form a directed edge chain.
   *
   * <p>This function could possibly be moved to S2BuilderGraph, but note that it has special
   * requirements. Namely, duplicate edges and sibling pairs must be kept in order to ensure that
   * every output edge corresponds to exactly one input edge. See also {@link
   * S2BuilderGraph#getInputEdgeOrder()}.
   */
  private static IntArrayList getInputEdgeChainOrder(S2BuilderGraph g, IntArrayList inputEdgeIds) {
    Preconditions.checkState(g.options().edgeType() == EdgeType.DIRECTED);
    Preconditions.checkState(g.options().duplicateEdges() == DuplicateEdges.KEEP);
    Preconditions.checkState(g.options().siblingPairs() == SiblingPairs.KEEP);

    // First, sort the edges so that the edges corresponding to each input edge are consecutive.
    // (Each input edge was snapped to a chain of output edges, or two chains in the case of
    // undirected input edges.)
    IntArrayList order = g.getInputEdgeOrder(inputEdgeIds);

    // Now sort the group of edges corresponding to each input edge in edge chain order (e.g. AB,
    // BC, CD).
    IntPairVector vmap = new IntPairVector(); // Map from source vertex id to edge id.
    // TODO(torrey): Consider using a BitSet for indegree, as the only thing that matters is if it
    // is nonzero or not.
    int[] indegree = new int[g.numVertices()]; // Restricted to current input edge.
    for (int end, begin = 0; begin < order.size(); begin = end) {
      // Gather the edges that came from a single input edge.
      int inputEdgeId = inputEdgeIds.getInt(order.getInt(begin));
      for (end = begin; end < order.size(); ++end) {
        if (inputEdgeIds.getInt(order.getInt(end)) != inputEdgeId) {
          break;
        }
      }
      if (end - begin == 1) {
        continue;
      }

      // Build a map from the source vertex of each edge to its edge id, and also compute the
      // indegree at each vertex considering only the edges that came from the current input edge.
      for (int i = begin; i < end; ++i) {
        int edgeId = order.getInt(i);
        vmap.add(g.edgeSrcId(edgeId), edgeId);
        indegree[g.edgeDstId(edgeId)] += 1;
      }
      vmap.sort();

      // Find the starting edge for building the edge chain.
      int nextEdgeId = g.numEdges();
      for (int i = begin; i < end; ++i) {
        int edgeId = order.getInt(i);
        if (indegree[g.edgeSrcId(edgeId)] == 0) {
          nextEdgeId = edgeId;
        }
      }
      // Build the edge chain.
      for (int i = begin; ; ) {
        order.set(i, nextEdgeId);
        int vertexId = g.edgeDstId(nextEdgeId);
        indegree[vertexId] = 0; // Clear as we go along.
        if (++i == end) {
          break;
        }
        int outIndex = vmap.lowerBound(vertexId, 0);
        assert vertexId == vmap.getFirst(outIndex);
        nextEdgeId = vmap.getSecond(outIndex);
      }
      vmap.clear();
    }
    return order;
  }

  /**
   * Given a set of clipping instructions encoded as an InputEdgeCrossings, GraphEdgeClipper
   * determines which graph edges correspond to clipped portions of input edges and removes them.
   *
   * <p>The clipping model is as follows. The input consists of edge chains. The clipper maintains
   * an "inside" boolean state as it clips each chain, and toggles this state whenever an input edge
   * is crossed. Any edges that are deemed to be "outside" after clipping are removed.
   *
   * <p>The "inside" state can be reset when necessary (e.g., when jumping to the start of a new
   * chain) by adding a special crossing marked SET_INSIDE. There are also two other special
   * "crossings" that modify the clipping parameters: SET_INVERT_B specifies that edges should be
   * clipped to the exterior of the other region, and SET_REVERSE_A specifies that edges should be
   * reversed before emitting them (which is needed to implement difference operations).
   */
  private static class GraphEdgeClipper {

    /**
     * Each CrossingInputEdge represents an input edge B that crosses some other input edge A. The
     * CrossingInputEdgeVector is a list of these, constructed in order by input edge id.
     */
    private static class CrossingInputEdgeVector extends ArrayList<CrossingInputEdge> {
      // TODO(torrey): Replace CrossingInputEdgeVector with an IntVector that contains the input
      // edge ids, and encodes the leftToRight in the sign bit.
      /**
       * Returns the index of the first CrossingInputEdge with an input edge id greater than or
       * equal to the specified input edge id. Requires that this CrossingInputEdgeVector is in
       * sorted order by input edge id.
       */
      public int lowerBound(int inputEdgeId) {
        return S2ShapeUtil.lowerBound(0, size(), target -> get(target).inputEdgeId() < inputEdgeId);
      }
    }

    /** The S2BuilderGraph being processed. */
    private final S2BuilderGraph g;

    /** VertexInMap of the S2BuilderGraph being processed. Maps vertex ids to incoming edge ids. */
    private final VertexInMap in;

    /** VertexOutMap of the S2BuilderGraph being processed. Maps vertex ids to outgoing edge ids. */
    private final VertexOutMap out;

    /** The dimension of each input edge (0, 1 or 2). */
    private final IntVector inputDimensions;

    /**
     * The set of all crossings to be used when clipping the edges of "g", sorted in lexicographic
     * order.
     */
    private final InputEdgeCrossings inputCrossings;

    /**
     * The clipped set of edges and their corresponding set of input edge ids. (These can be used to
     * construct a new S2BuilderGraph.)
     */
    private final EdgeList newEdges;

    /** The corresponding input edge ids of the clipped set of edges in 'newEdges'. */
    private final IntArrayList newInputEdgeIds;

    /**
     * Every graph edge is associated with exactly one input edge in our case, which means that
     * inputEdgeIds is a vector of InputEdgeIds rather than a vector of InputEdgeIdSetIds. (This
     * also relies on the fact that IdSetLexicon represents a singleton set as the value of its
     * single element.)
     */
    private final IntArrayList inputEdgeIds;

    /**
     * Graph edge ids sorted in input edge id order. When more than one output edge has the same
     * input edge id (i.e., the input edge snapped to a chain of edges), the edges are sorted so
     * that they form a directed edge chain.
     */
    private final IntArrayList order;

    /** The rank of each graph edge is its index in the sorted "order". */
    private final int[] rank;

    /** For using an IntVector to store a boolean 'true'. */
    private static final int TRUE = 1;

    /** For using an IntVector to store a boolean 'false'. */
    private static final int FALSE = 0;

    /**
     * Constructs a new GraphEdgeClipper.
     *
     * <p>"inputDimensions" is a vector specifying the dimension of each input edge (0, 1, or 2).
     * "inputCrossings" is a list of CrossingInputEdgePair, containing special instructions to
     * control state (SET_INSIDE, SET_REVERSE_A, SET_INVERT_B) as well as all crossings to be used
     * when clipping the edges of "g", sorted in lexicographic order.
     *
     * <p>The clipped set of edges and their corresponding set of input edge ids are returned in
     * "newEdges" and "newInputEdgeIds". (These can be used to construct a new S2BuilderGraph.)
     */
    public GraphEdgeClipper(
        S2BuilderGraph g,
        IntVector inputDimensions,
        InputEdgeCrossings inputCrossings,
        EdgeList newEdges,
        IntArrayList newInputEdgeIds) {
      this.g = g;
      this.in = new VertexInMap(g);
      this.out = new VertexOutMap(g);
      this.inputDimensions = inputDimensions;
      this.inputCrossings = inputCrossings;
      this.newEdges = newEdges;
      this.newInputEdgeIds = newInputEdgeIds;
      this.inputEdgeIds = g.inputEdgeIdSetIds();
      order = getInputEdgeChainOrder(g, inputEdgeIds);
      rank = new int[order.size()];
      // order.forEach((k, v) -> rank[v] = k);
      for (int k = 0; k < order.size(); k++) {
        rank[order.getInt(k)] = k;
      }

      // newEdges is obtained by filtering the graph edges and therefore the number of graph edges
      // is an upper bound on its size.
      newEdges.ensureCapacity(g.numEdges());
      newInputEdgeIds.ensureCapacity(g.numEdges());
    }

    /** Performs the graph clipping operation. */
    @SuppressWarnings("IncrementInForLoopAndHeader")
    public void run() {
      // Vectors are initialized once here and reused to avoid reallocation.

      // Vertex ids that the current input edge was snapped to.
      IntVector aVertices = new IntVector();
      // The sum of signed crossings at each A vertex.
      IntVector aNumCrossings = new IntVector();

      // Is the A vertex "isolated" (TRUE)? Otherwise FALSE. A vertex is isolated if it has at least
      // one crossing, but no adjacent emitted edge.
      IntVector aIsolated = new IntVector(); // Actually a vector of boolean

       // A list of all input edges that cross the single input edge corresponding to the graph
       // edge(s) being processed in the current iteration of the outer for loop.
      CrossingInputEdgeVector bInputEdges = new CrossingInputEdgeVector();

      // A list of CrossingGraphEdgeVector, each of which contains CrossingGraphEdges.
      // TODO(torrey): Document the semantics of bEdges, here and in C++.
      CrossingGraphEdgeVectorVector bEdges = new CrossingGraphEdgeVectorVector();

      boolean inside = false;
      boolean invertB = false;
      boolean reverseA = false;

      // Iterates through the inputCrossings in lexicographic order.
      int inputCrossingsIndex = 0;

      // Iterate over the graph edges in order of input edge id, as specified by 'order'. Each
      // iteration of this for loop processes all the snapped graph edges corresponding to a single
      // input edge 'A'.

      // 'i' is advanced to a new input edge at the top of the loop, and also advanced within the
      // loop over the other snapped graph edges for input edge 'A'.
      for (int i = 0; i < order.size(); ++i) {
        // edge0 is the first graph edge for input edge A.
        int edge0Id = order.getInt(i);
        int aInputEdgeId = inputEdgeIds.getInt(edge0Id);
        int edge0SrcId = g.edgeSrcId(edge0Id);
        int edge0DstId = g.edgeDstId(edge0Id);

        // For the current input edge (the "A" input edge), gather all the input edges that cross it
        // (the "B" input edges), as well as any the state change instructions that apply. Note that
        // bInputEdges is constructed in sorted order by input edge id, as inputCrossings is also
        // sorted in lexicographic order.
        bInputEdges.clear();
        for (; inputCrossingsIndex < inputCrossings.size(); ++inputCrossingsIndex) {
          CrossingInputEdgePair crossingPair = inputCrossings.get(inputCrossingsIndex);

          if (crossingPair.first() != aInputEdgeId) {
            break; // No more crossings or state changes for this input edge id.
          }
          if (crossingPair.second().inputEdgeId() >= 0) {
            // A crossing of aInputEdgeId by some other input edge.
            bInputEdges.add(crossingPair.second());
          } else if (crossingPair.second().inputEdgeId() == SET_INSIDE) {
            inside = crossingPair.second().leftToRight();
          } else if (crossingPair.second().inputEdgeId() == SET_INVERT_B) {
            invertB = crossingPair.second().leftToRight();
          } else {
            assert crossingPair.second().inputEdgeId() == SET_REVERSE_A;
            reverseA = crossingPair.second().leftToRight();
          }
        }

        // Optimization for degenerate edges.
        // TODO(ericv): If the output layer for this edge dimension specifies
        // DegenerateEdges.DISCARD, then remove the edge here.
        if (edge0SrcId == edge0DstId) {
          inside ^= ((bInputEdges.size() & 1) == 1);
          addEdge(edge0SrcId, edge0DstId, aInputEdgeId);
          continue;
        }

        // Optimization for the case where there are no crossings of this edge.
        if (bInputEdges.isEmpty()) {
          // In general the caller only passes edges that are part of the output (i.e., we could
          // assert(inside) here). The one exception is for polyline/polygon operations, where the
          // polygon edges are needed to compute the polyline output but are not emitted themselves.
          if (inside) {
            if (reverseA) {
              addEdge(edge0DstId, edge0SrcId, aInputEdgeId);
            } else {
              addEdge(edge0SrcId, edge0DstId, aInputEdgeId);
            }
          }
          continue;
        }

        // Walk along the chain of snapped edges for input edge A, and at each vertex collect all
        // the incident snapped edges that belong to one of the crossing edge chains (the "B" input
        // edges).

        // aVertices are the vertex ids of the snapped edges for input edge A. Add the first one.
        aVertices.clear();
        aVertices.add(edge0SrcId);

        // Get the snapped crossing edges incident to vertex 0 of the snapped edge of input edge A.
        bEdges.reinitialize(bInputEdges.size());
        gatherIncidentEdges(aVertices, 0, bInputEdges, bEdges);
        // Walk along the snapped edge chain for input edge A, adding their endpoints to aVertices
        // and getting the snapped crossing edges incident to each of those vertices.
        for (; i < order.size() && inputEdgeIds.getInt(order.getInt(i)) == aInputEdgeId; ++i) {
          aVertices.add(g.edgeDstId(order.getInt(i)));
          gatherIncidentEdges(aVertices, aVertices.size() - 1, bInputEdges, bEdges);
        }

        --i; // back up so we're still on input edge A

        // Now for each B edge chain, decide which vertex of the A chain it crosses, and keep track
        // of the sum of signed crossings at each A vertex. The sign of a crossing depends on
        // whether the other edge crosses from left to right or right to left.
        //
        // This would not be necessary if all calculations were done in exact arithmetic, because
        // crossings would have strictly alternating signs. But because we have already snapped the
        // result, some crossing locations are ambiguous, and getCrossedVertexIndex() handles this
        // by choosing a candidate vertex arbitrarily. The end result is that rarely, we may see two
        // crossings in a row with the same sign. We correct for this by adding extra output edges
        // that essentially link up the crossings in the correct (alternating sign) order. Compared
        // to the "correct" behavior, the only difference is that we have added some extra sibling
        // pairs (consisting of an edge and its corresponding reverse edge) which do not affect the
        // result.
        aNumCrossings.clear();
        aNumCrossings.resize(aVertices.size());
        aIsolated.clear();
        aIsolated.resize(aVertices.size());
        for (int bi = 0; bi < bInputEdges.size(); ++bi) {
          CrossingInputEdge bCrossingEdge = bInputEdges.get(bi);
          int bInputEdgeId = bCrossingEdge.inputEdgeId();
          boolean leftToRight = bCrossingEdge.leftToRight();
          int aIndex = getCrossedVertexIndex(aVertices, bEdges.get(bi), leftToRight);
          if (aIndex < 0) {
            // TODO(user): This fails in some cases. Need to get to the bottom of it.
            throw new IllegalStateException("Failed to get crossed vertex index.");
          }

          // Keep track of the number of signed crossings (see above).
          boolean isLine = inputDimensions.get(bInputEdgeId) == 1;
          int sign = isLine ? 0 : (leftToRight == invertB) ? -1 : 1;
          aNumCrossings.incrementAt(aIndex, sign);

          // Any polyline or polygon vertex that has at least one crossing but no adjacent
          // emitted edge may be emitted as an isolated vertex.
          aIsolated.set(aIndex, TRUE);
        }

        // Finally, we iterate through the A edge chain, keeping track of the number of signed
        // crossings as we go along. The "multiplicity" is defined as the cumulative number of
        // signed crossings, and indicates how many edges should be output (and in which direction)
        // in order to link up the edge crossings in the correct order. (The multiplicity is almost
        // always either 0 or 1 except in very rare cases.)
        int multiplicity = (inside ? 1 : 0) + aNumCrossings.get(0);

        for (int ai = 1; ai < aVertices.size(); ++ai) {
          if (multiplicity != 0) {
            aIsolated.set(ai - 1, FALSE);
            aIsolated.set(ai, FALSE);
          }
          int edgeCount = reverseA ? -multiplicity : multiplicity;
          // Output any forward edges required.
          for (int ei = 0; ei < edgeCount; ++ei) {
            addEdge(aVertices.get(ai - 1), aVertices.get(ai), aInputEdgeId);
          }
          // Output any reverse edges required.
          for (int ei = edgeCount; ei < 0; ++ei) {
            addEdge(aVertices.get(ai), aVertices.get(ai - 1), aInputEdgeId);
          }
          multiplicity += aNumCrossings.get(ai);
        }
        // Multiplicities other than 0 or 1 can only occur in the edge interior.
        // If this assertion fails, it is likely due to invalid inputs.
        assert multiplicity == 0 || multiplicity == 1;
        inside = (multiplicity != 0);

        // Output any isolated polyline vertices.
        // TODO(ericv): Only do this if an output layer wants degenerate edges.
        if (inputDimensions.get(aInputEdgeId) != 0) {
          for (int ai = 0; ai < aVertices.size(); ++ai) {
            if (aIsolated.get(ai) == TRUE) {
              addEdge(aVertices.get(ai), aVertices.get(ai), aInputEdgeId);
            }
          }
        }
      }
    }

    private void addEdge(int edgeSrcId, int edgeDstId, int inputEdgeId) {
      newEdges.add(edgeSrcId, edgeDstId);
      newInputEdgeIds.add(inputEdgeId);
    }

    /**
     * Given the vertices of the snapped edge chain for an input edge A, and the list of input edges
     * B that cross input edge A, ordered by input edge id, this method gathers all the snapped
     * edges of B that are incident to a given snapped vertex of A. The incident edges for each
     * input edge of B are appended to a separate output vector. (A and B can refer to either the
     * input edge or the corresponding snapped edge chain.)
     *
     * @param a Snapped edge chain for an input edge A, as vertex ids
     * @param ai The given snapped vertex of A, as an index into 'a'
     * @param bInputEdges All the B input edges that cross input edge A
     * @param bEdges Output: A list of vectors of snapped edges of B incident to vertex ai
     */
    private void gatherIncidentEdges(
        IntVector a,
        int ai,
        CrossingInputEdgeVector bInputEdges,
        CrossingGraphEdgeVectorVector bEdges) {
      // Examine all the snapped edges incident to the given vertex of A. If any edge comes from a B
      // input edge, append it to the appropriate CrossingGraphEdgeVector.
      assert bInputEdges.size() == bEdges.size();
      in.edgeIds(a.get(ai)) // Incoming snapped edge ids at the snapped vertex of A, from A and B
          .forEach(
              snappedEdgeId -> {
                // The input edge id of this incoming snapped edge to ai.
                int inputEdgeId = inputEdgeIds.getInt(snappedEdgeId);
                // Index of the first input edge of B that crosses this input edge, if any
                int index = bInputEdges.lowerBound(inputEdgeId);
                // A matching, crossing, input edge of B crosses here, at snapped vertex 'ai'
                if (index != bInputEdges.size()
                    && bInputEdges.get(index).inputEdgeId() == inputEdgeId) {
                  // Add a CrossingGraphEdge at the index of that input edge.
                  bEdges.addCrossingGraphEdge(
                      index, snappedEdgeId, ai, false, g.edgeSrcId(snappedEdgeId));
                }
              });
      out.edgeIds(a.get(ai))
          .forEach(
              snappedEdgeId -> {
                int inputEdgeId = inputEdgeIds.getInt(snappedEdgeId);
                int index = bInputEdges.lowerBound(inputEdgeId);
                if (index != bInputEdges.size()
                    && bInputEdges.get(index).inputEdgeId() == inputEdgeId) {
                  bEdges.addCrossingGraphEdge(
                      index, snappedEdgeId, ai, true, g.edgeDstId(snappedEdgeId));
                }
              });
    }

    /**
     * Given an edge chain A that is crossed by another edge chain B (where "leftToRight" indicates
     * whether B crosses A from left to right), this method decides which vertex of A the crossing
     * takes place at. The parameters are the vertices of the A chain ("a") and the set of edges in
     * the B chain ("b") that are incident to vertices of A. The B chain edges are sorted in
     * increasing order of (aIndex, outgoing) tuple.
     */
    private int getCrossedVertexIndex(IntVector a, CrossingGraphEdgeVector b, boolean leftToRight) {
      if (a.isEmpty() || b.isEmpty()) {
        throw new IllegalStateException(
            "GraphEdgeClipper.getCrossedVertexIndex called with "
                + a.size()
                + " vertex ids and "
                + b.size()
                + " crossing graph edges.");
      }

      // TODO(user): Rework this algorithm.

      // The reason this calculation is tricky is that after snapping, the A and B chains may meet
      // and separate several times. For example, if B crosses A from left to right, then B may
      // touch A, make an excursion to the left of A, come back to A, then make an excursion to the
      // right of A and come back to A again, like this:
      //
      //  *--B--*-\             /-*-\
      //           B-\       /-B     B-\      6     7     8     9
      //  *--A--*--A--*-A,B-*--A--*--A--*-A,B-*--A--*--A--*-A,B-*
      //  0     1     2     3     4     5      \-B     B-/
      //                                          \-*-/
      //
      // (where "*" is a vertex, and "A" and "B" are edge labels). Note that B may also follow A for
      // one or more edges whenever they touch (e.g. between vertices 2 and 3). In this case the
      // only vertices of A where the crossing could take place are 5 and 6, i.e. after all
      // excursions of B to the left of A, and before all excursions of B to the right of A.
      //
      // Other factors to consider are that the portion of B before and/or after the crossing may be
      // degenerate, and some or all of the B edges may be reversed relative to the A edges.

      // First, check whether edge A is degenerate.
      int n = a.size();
      if (n == 1) {
        return 0;
      }

      // If edge chain B is incident to only one vertex of A, we're done.
      if (b.firstAIndex() == b.lastAIndex()) {
        return b.firstAIndex();
      }

      // Determine whether the B chain visits the first and last vertices that it shares with the A
      // chain in the same order or the reverse order. This is only needed to implement one special
      // case (see below).
      boolean bReversed = getVertexRank(b.first()) > getVertexRank(b.last());

      // Examine each incident B edge and use it to narrow the range of positions where the crossing
      // could occur in the B chain. Vertex positions are represented as a range [lo, hi] of vertex
      // ranks in the B chain (see getVertexRank).
      //
      // Note that if an edge of B is incident to the first or last vertex of A, we can't test which
      // side of the A chain it is on. (An S2Predicates.Sign test doesn't work; e.g. if the B edge
      // is XY and the first edge of A is YZ, then snapping can change the sign of XYZ while
      // maintaining topological guarantees.) There can be up to 4 such edges (one incoming and one
      // outgoing edge at each endpoint of A). Two of these edges logically extend past the end of
      // the A chain and place no restrictions on the crossing vertex. The other two edges define
      // the ends of the subchain where B shares vertices with A. We save these edges in order to
      // handle a special case (see below).

      // Vertex ranks of acceptable crossings
      int lo = -1;
      int hi = order.size();

      // "b" subchain connecting "a" endpoints
      int bFirstEdgeId = -1;
      int bLastEdgeId = -1;

      for (CrossingGraphEdge e : b) {
        int ai = e.aIndex;
        if (ai == 0) {
          if (e.outgoing != bReversed && e.otherVertexId != a.get(1)) {
            bFirstEdgeId = e.edgeId;
          }
        } else if (ai == n - 1) {
          if (e.outgoing == bReversed && e.otherVertexId != a.get(n - 2)) {
            bLastEdgeId = e.edgeId;
          }
        } else {
          // This B edge is incident to an interior vertex of the A chain. First check whether this
          // edge is identical (or reversed) to an edge in the A chain, in which case it does not
          // create any restrictions.
          if (e.otherVertexId == a.get(ai - 1) || e.otherVertexId == a.get(ai + 1)) {
            continue;
          }

          // Otherwise we can test which side of the A chain the edge lies on.
          boolean onLeft =
              S2Predicates.orderedCCW(
                  g.vertex(a.get(ai + 1)),
                  g.vertex(e.otherVertexId),
                  g.vertex(a.get(ai - 1)),
                  g.vertex(a.get(ai)));

          // Every B edge that is incident to an interior vertex of the A chain places some
          // restriction on where the crossing vertex could be.
          if (leftToRight == onLeft) {
            // This is a pre-crossing edge, so the crossing cannot be before the destination vertex
            // of this edge. (For example, the input B edge crosses the input A edge from left to
            // right and this edge of the B chain is to the left of the A chain.)
            lo = max(lo, rank[e.edgeId] + 1);
          } else {
            // This is a post-crossing edge, so the crossing cannot be after the source vertex of
            // this edge.
            hi = min(hi, rank[e.edgeId]);
          }
        }
      }

      // There is one special case. If a subchain of B connects the first and last vertices of A,
      // then together with the edges of A this forms a loop whose orientation can be tested to
      // determine whether B is on the left or right side of A. This is only possible (and only
      // necessary) if the B subchain does not include any interior vertices of A, since otherwise
      // the B chain might cross from one side of A to the other.
      //
      // Note that it would be possible to avoid this test in some situations by checking whether
      // either endpoint of the A chain has two incident B edges, in which case we could check which
      // side of the B chain the A edge is on and use this to limit the possible crossing locations.
      if (bFirstEdgeId >= 0 && bLastEdgeId >= 0) {
        // Swap the edges if necessary so that they are in B chain order.
        if (bReversed) {
          int tmp = bFirstEdgeId;
          bLastEdgeId = bFirstEdgeId;
          bFirstEdgeId = tmp;
        }

        // The B subchain connects the first and last vertices of A. We test whether the chain
        // includes any interior vertices of A by iterating through the incident B edges again,
        // looking for ones that belong to the B subchain and are not incident to the first or last
        // vertex of A.
        boolean hasInteriorVertex = false;
        for (CrossingGraphEdge e : b) {
          if (e.aIndex > 0
              && e.aIndex < n - 1
              && rank[e.edgeId] >= rank[bFirstEdgeId]
              && rank[e.edgeId] <= rank[bLastEdgeId]) {
            hasInteriorVertex = true;
            break;
          }
        }
        if (!hasInteriorVertex) {
          // The B subchain is not incident to any interior vertex of A.
          boolean onLeft = edgeChainOnLeft(a, bFirstEdgeId, bLastEdgeId);
          if (leftToRight == onLeft) {
            lo = max(lo, rank[bLastEdgeId] + 1);
          } else {
            hi = min(hi, rank[bFirstEdgeId]);
          }
        }
      }

      // Otherwise we choose the smallest shared VertexId in the acceptable range, in order to
      // ensure that both chains choose the same crossing vertex.
      int best = -1;
      assert lo <= hi;
      for (CrossingGraphEdge e : b) {
        int ai = e.aIndex;
        int vrank = getVertexRank(e);
        if (vrank >= lo && vrank <= hi && (best < 0 || a.get(ai) < a.get(best))) {
          best = ai;
        }
      }
      return best;
    }

    /**
     * Returns the "vertex rank" of the shared vertex associated with the given CrossingGraphEdge.
     * Recall that graph edges are sorted in input edge order, and that the rank of an edge is its
     * position in this order (rank[e]). VertexRank(e) is defined such that VertexRank(e.src) ==
     * rank[e] and VertexRank(e.dst) == rank[e] + 1. Note that the concept of "vertex rank" is only
     * defined within a single edge chain (since different edge chains can have overlapping vertex
     * ranks).
     */
    private int getVertexRank(CrossingGraphEdge e) {
      return rank[e.edgeId] + (!e.outgoing ? 1 : 0);
    }

    /**
     * Given edge chains A and B that form a loop (after possibly reversing the direction of chain
     * B), returns true if chain B is to the left of chain A. Chain A is given as a sequence of
     * vertices, while chain B is specified as the first and last edges of the chain.
     */
    private boolean edgeChainOnLeft(IntVector a, int bFirstEdgeId, int bLastEdgeId) {
      // TODO(user): Rather than collecting all the vertices in a list, accumulate the angle.

      // Gather all the interior vertices of the B subchain.
      IntVector loop = new IntVector(); // vertex ids
      for (int i = rank[bFirstEdgeId]; i < rank[bLastEdgeId]; ++i) {
        loop.push(g.edgeDstId(order.getInt(i)));
      }
      // Possibly reverse the chain so that it forms a loop when "a" is appended.
      if (g.edgeDstId(bLastEdgeId) != a.get(0)) {
        loop.reverse();
      }
      // Append all the elements of 'a' to the loop.
      loop.insert(loop.size(), a);
      // Append duplicates of the first two vertices to simplify vertex indexing.
      for (int j = 0; j < 2; j++) {
        loop.insert(loop.size(), loop.get(j));
      }
      // Now B is to the left of A if and only if the loop is counterclockwise.
      double sum = 0;
      for (int i = 2; i < loop.size(); ++i) {
        sum +=
            S2.turnAngle(
                g.vertex(loop.get(i - 2)), g.vertex(loop.get(i - 1)), g.vertex(loop.get(i)));
      }
      return sum > 0;
    }
  }

  /**
   * Given a set of clipping instructions encoded as a set of intersections between input edges,
   * EdgeClippingLayer determines which graph edges correspond to clipped portions of input edges
   * and removes them. It assembles the remaining edges into a new S2BuilderGraph and passes the
   * result to the given output layer or layers for assembly.
   */
  private static class EdgeClippingLayer implements S2BuilderLayer {
    /** The output layer or layers. */
    private final List<S2BuilderLayer> layers;

    /** The dimensions of the input edges. */
    private final IntVector inputDimensions;

    /** The crossings between the input edges, along with GraphEdgeClipper state modifications. */
    private final InputEdgeCrossings inputCrossings;

    /**
     * Constructs a new EdgeClippingLayer that will output the remaining edges after clipping to the
     * given layers, according to their dimension.
     */
    public EdgeClippingLayer(
        List<S2BuilderLayer> layers, IntVector inputDimensions, InputEdgeCrossings inputCrossings) {
      this.layers = layers;
      this.inputDimensions = inputDimensions;
      this.inputCrossings = inputCrossings;
    }

    @Override
    public String toString() {
      return "S2BooleanOperation.EdgeClippingLayer with " + layers.size() + " output layers.";
    }

    /**
     * GraphOptions required by EdgeClippingLayer. We keep all edges, including degenerate ones, so
     * that we can figure out the correspondence between input edge crossings and output edge
     * crossings.
     */
    @Override
    public S2Builder.GraphOptions graphOptions() {
      return new S2Builder.GraphOptions(
          EdgeType.DIRECTED, DegenerateEdges.KEEP, DuplicateEdges.KEEP, SiblingPairs.KEEP);
    }

    /**
     * Builds this layer. Runs the GraphEdgeClipper, and then sends the remaining edges to the
     * output layer or layers, and then builds that layer or layers in turn.
     */
    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      // The bulk of the work is handled by GraphEdgeClipper.
      EdgeList newEdges = new EdgeList();
      IntArrayList newInputEdgeIds = new IntArrayList();

      GraphEdgeClipper gec =
          new GraphEdgeClipper(g, inputDimensions, inputCrossings, newEdges, newInputEdgeIds);
      gec.run();
      // TODO(torrey): Consider making the S2BooleanOperation.Impl, S2Builder, EdgeClippingLayer,
      // and GraphEdgeClipper reusable for repeated S2BooleanOperation operations. Otherwise, it
      // may be worth triggering garbage collection here, below, and/or when an operation finishes.
      gec = null;

      // Construct one or more subgraphs from the clipped edges and pass them to the given output
      // layer(s). We start with a copy of the input graph's IdSetLexicon because this is necessary
      // in general, even though in this case it is guaranteed to be empty because no edges have
      // been merged.
      IdSetLexicon newInputEdgeIdSetLexicon = g.inputEdgeIdSetLexicon();
      if (layers.size() == 1) {
        S2BuilderGraph newGraph =
            g.makeSubgraph(
                layers.get(0).graphOptions(),
                newEdges,
                newInputEdgeIds,
                newInputEdgeIdSetLexicon,
                g.isFullPolygonPredicate(),
                error);
        boolean unused = layers.get(0).build(newGraph, error);
      } else {
        // The Graph objects must be valid until the last build() call completes, so we store all
        // of the graph data in arrays with 3 elements.
        Preconditions.checkState(layers.size() == 3);
        EdgeList[] layerEdges = new EdgeList[3];
        IntArrayList[] layerInputEdgeIds = new IntArrayList[3];
        for (int d = 0; d < 3; ++d) {
          layerEdges[d] = new EdgeList();
          layerInputEdgeIds[d] = new IntArrayList();
        }

        // Separate the edges according to their dimension.
        for (int i = 0; i < newEdges.size(); ++i) {
          int newInputEdgeId = newInputEdgeIds.getInt(i);
          int d = inputDimensions.get(newInputEdgeId);
          layerEdges[d].copyEdge(newEdges, i);
          layerInputEdgeIds[d].add(newInputEdgeId);
        }

        // TODO(torrey): As above, consider options regarding clearing variables to save space.
        newEdges.clear();
        newInputEdgeIds.clear();

        S2BuilderGraph[] layerGraphs = new S2BuilderGraph[3];
        for (int d = 0; d < 3; ++d) {
          layerGraphs[d] =
              g.makeSubgraph(
                  layers.get(d).graphOptions(),
                  layerEdges[d],
                  layerInputEdgeIds[d],
                  newInputEdgeIdSetLexicon,
                  g.isFullPolygonPredicate(),
                  error);
          boolean unused = layers.get(d).build(layerGraphs[d], error);
        }
      }
      return error.ok();
    }
  }

  /**
   * ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex, consisting of a shapeId,
   * edgeId pair. The shapeId is the index of the shape in the S2ShapeIndex shape list.
   *
   * <p>TODO: Consider if this class is actually required, or if we really only need collections of
   * (shapeId, edgeId) that can be stored more efficiently. ShapeEdgeIds are used in:
   *
   * <ol>
   *   <li>class ShapeEdge, along with the two S2Points
   *   <li>a {@code List<ShapeEdgeId>} for chainStarts, which is iterated over, ends with SENTINEL
   *   <li>An IndexCrossing has two ShapeEdgeIds.
   *   <li>Part of the elements iterated over by a CrossingIterator
   *   <li>The minimum of two ShapeEdgeIds is selected during iteration and referenced
   *   <li>Constructed while processing edges of the chain
   *   <li>Stored as keys in a HashMap to Boolean, isDegenerateHole
   * </ol>
   */
  private static class ShapeEdgeId implements Comparable<ShapeEdgeId> {
    /** SENTINEL is used to mark the end of a collection of ShapeEdgeIds. */
    public static final ShapeEdgeId SENTINEL = new ShapeEdgeId(Integer.MAX_VALUE, 0);

    /** An instance value that indicates "no edge". */
    public static final ShapeEdgeId NONE = new ShapeEdgeId(-1, -1);

    /** The shape id of this ShapeEdgeId. */
    public final int shapeId;

    /** The edge id of this ShapeEdgeId. */
    public final int edgeId;

    public ShapeEdgeId(int shapeId, int edgeId) {
      this.shapeId = shapeId;
      this.edgeId = edgeId;
    }

    public boolean isEqualTo(ShapeEdgeId other) {
      return this.shapeId == other.shapeId && this.edgeId == other.edgeId;
    }

    @Override
    public int hashCode() {
      return Objects.hash(shapeId, edgeId);
    }

    @Override
    public boolean equals(Object obj) {
      return obj instanceof ShapeEdgeId && isEqualTo((ShapeEdgeId) obj);
    }

    /** Orders ShapeEdgeIds by shapeId, breaking ties with edgeId. */
    @Override
    public int compareTo(ShapeEdgeId other) {
      int result = Integer.compare(this.shapeId, other.shapeId);
      if (result != 0) {
        return result;
      }
      return Integer.compare(this.edgeId, other.edgeId);
    }

    public static ShapeEdgeId min(ShapeEdgeId a, ShapeEdgeId b) {
      return (a.compareTo(b) < 0) ? a : b;
    }

    @Override
    public String toString() {
      if (shapeId == Integer.MAX_VALUE && edgeId == 0) {
        return "ShapeEdgeId.SENTINEL";
      }
      return "ShapeEdgeId(" + shapeId + ", " + edgeId + ")";
    }
  }

  /**
   * A ShapeEdgeId (shape id and edge id) along with a MutableEdge containing the edge endpoints.
   */
  private static class ShapeEdge {
    public final ShapeEdgeId id;
    public final S2Point v0;
    public final S2Point v1;

    public ShapeEdge(int shapeId, int edgeId, S2Point v0, S2Point v1) {
      this.id = new ShapeEdgeId(shapeId, edgeId);
      this.v0 = v0;
      this.v1 = v1;
    }

    @Override
    public String toString() {
      return "ShapeEdge(" + id + ", " + v0 + ", " + v1 + ")";
    }
  }

  /**
   * Provides the actual implementation of S2BooleanOperation, clipping each of the two regions to
   * the other region.
   */
  private static class Impl {
    /** The S2BooleanOperation that this Impl is performing. */
    private final S2BooleanOperation op;

    /** The S2Builder options used to construct the output. */
    private S2Builder.Builder builderOptions;

    /**
     * If isBooleanOutput, no output geometry is created, and only the existence of some output
     * geometry is determined.
     */
    private final boolean isBooleanOutput;

    // The S2Builder used to construct the output. Note that the S2Builder object is created only
    // when isBooleanOutput is false.
    private S2Builder builder;

    // A vector specifying the dimension of each edge added to S2Builder.
    // TODO(torrey): Consider creating a ByteVector
    private final IntVector inputDimensions = new IntVector();

    // The set of all input edge crossings, which is used by EdgeClippingLayer to construct the
    // clipped output polygon.
    private final InputEdgeCrossings inputCrossings = new InputEdgeCrossings();

    // A vector containing all pairs of crossing edges from the two input regions (including edge
    // pairs that share a common vertex). The first element of each pair is an edge from
    // "indexCrossingsFirstRegionId", while the second element of each pair is an edge from the
    // other region.
    private final IndexCrossings indexCrossings = new IndexCrossings();

    // Indicates that the first element of each crossing edge pair in "indexCrossings" corresponds
    // to an edge from the given region. This field is negative if indexCrossings has not been
    // computed yet.
    private int indexCrossingsFirstRegionId;

    // Temporary storage used in getChainStarts(), declared here to avoid repeatedly allocating
    // memory.
    private final IndexCrossings tmpCrossings = new IndexCrossings();

    /** Constructs an S2BooleanOperation.Impl for the given S2BooleanOperation. */
    public Impl(S2BooleanOperation op) {
      this.op = op;
      isBooleanOutput = op.layers.isEmpty();
      indexCrossingsFirstRegionId = -1;
    }

    /** Builds the output of the S2BooleanOperation. Returns true if the operation succeeded. */
    public boolean build(S2Error error) {
      error.clear();
      doBuild(error);
      return error.ok();
    }

    /**
     * Clips the boundary of A to the interior of the opposite region B and adds the resulting edges
     * to the output. Optionally, any combination of region A, region B, and the result may be
     * inverted, which allows operations such as union and difference to be implemented.
     *
     * <p>Note that when an input region is inverted with respect to the output (e.g., invertA !=
     * invertResult), all polygon edges are reversed and all points and polylines are discarded,
     * since the complement of such objects cannot be represented. (If you want to compute the
     * complement of points or polylines, you can use S2LaxPolygonShape to represent your geometry
     * as degenerate polygons instead.)
     *
     * <p>This method must be called an even number of times (first to clip A to B and then to clip
     * B to A), calling doneBoundaryPair() after each pair.
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean addBoundary(
        int aRegionId,
        boolean invertA,
        boolean invertB,
        boolean invertResult,
        List<ShapeEdgeId> aChainStarts,
        CrossingProcessor crossingProcessor) {
      S2ShapeIndex aIndex = op.regions[aRegionId];
      S2ShapeIndex bIndex = op.regions[1 - aRegionId];
      // Initialize indexCrossings to the set of crossing edge pairs such that the first element of
      // each pair is an edge from "aRegionId".
      if (!initIndexCrossings(aRegionId)) {
        return false;
      }
      // Allocated once here, reused in the loop below.
      S2Shape.ChainPosition chainPosition = new S2Shape.ChainPosition();

      // Start building the boundary of A, clipped by B
      crossingProcessor.startBoundary(aRegionId, invertA, invertB, invertResult);

      // Walk the boundary of region A and build a list of all edge crossings. We also keep track of
      // whether the current vertex is "inside" region B.

      // Iterate over the shape edges of A that start "inside" B.
      Iterator<ShapeEdgeId> chainStartIterator = aChainStarts.iterator();
      // Get the first shape edge id "inside" B. If there is none, SENTINEL is returned.
      ShapeEdgeId chainStart = chainStartIterator.next();

      CrossingIterator nextCrossingIterator =
          new CrossingIterator(bIndex, indexCrossings, /* crossingsComplete == */ true);

      // The first ShapeEdgeId to process is the smaller of the current chainStart, or the current
      // (i.e. first) edge 'a' of the current (a, b) crossing.
      ShapeEdgeId nextId = ShapeEdgeId.min(chainStart, nextCrossingIterator.aId());

      // Keep advancing the iterators until both return SENTINEL.
      List<S2Shape> aShapes = aIndex.getShapes();
      while (!nextId.isEqualTo(ShapeEdgeId.SENTINEL)) {
        // Get the 'aShape' that the edge nextId belongs to, and start processing it.
        int aShapeId = nextId.shapeId;
        S2Shape aShape = aShapes.get(aShapeId);
        crossingProcessor.startShape(aShape);

        // Process edges for aShape.
        while (nextId.shapeId == aShapeId) {
          // TODO(torrey): Eric V. suggested special handling of an aShape of dimension 0. Could
          // omit most of this code, including the inner loop, since all chains are of length 1.

          // Get the chain information for the edge chain from aShape containing nextId.
          int edgeId = nextId.edgeId;
          aShape.getChainPosition(edgeId, chainPosition);
          int chainId = chainPosition.chainId;
          int chainLimit = aShape.getChainStart(chainId) + aShape.getChainLength(chainId);

          // If we're processing the current chainStart, we're "inside" B, and we advance
          // chainStartIterator.
          boolean startInside = nextId.isEqualTo(chainStart);
          if (startInside) {
            chainStart = chainStartIterator.next();
          }

          // Process edges of the aShape edge chain.
          crossingProcessor.startChain(chainId, aShape, startInside);

          while (edgeId < chainLimit) {
            // The next edge to process from the chain.
            ShapeEdgeId aId = new ShapeEdgeId(aShapeId, edgeId);

            // Loop invariant: We are "inside" B, or processing the set of crossings of the
            // particular edge of A that nextCrossingIterator is currently on.
            assert crossingProcessor.inside() || nextCrossingIterator.aId().isEqualTo(aId);

            // Process edge aId, advancing nextCrossingIterator over all crossings of aId (if any)
            // by edges from region B.
            if (!crossingProcessor.processEdge(aId, nextCrossingIterator)) {
              return false;
            }
            if (crossingProcessor.inside()) {
              // If we're still "inside" the output area B after processing crossings of edge aId,
              // go to the next edge of the current chain.
              ++edgeId;
            } else
            // Otherwise, after processing crossings of edge aId, we're no longer "inside" B. If
            // nextCrossingIterator is still on the current aShape, and hasn't reached the end of
            // the chain ...
            if (nextCrossingIterator.aId().shapeId == aShapeId
                && nextCrossingIterator.aId().edgeId < chainLimit) {
              // ... then the next edge of the current chain to process is the next one from
              // nextCrossingIterator.
              edgeId = nextCrossingIterator.aId().edgeId;
            } else {
              break; // Done with this chain.
            }
          }

          // We've finished a chain. Determine which edge to process next, which is the smaller
          // of the current chainStart, or the A edge of the current (a, b) crossing.
          nextId = ShapeEdgeId.min(chainStart, nextCrossingIterator.aId());
        }
      }
      return true;
    }

    /**
     * Returns the first edge of each edge chain from "aRegionId" whose first vertex is contained by
     * opposite region's polygons (using the semi-open boundary model). Each input region and the
     * result region are inverted as specified (invertA, invertB, and invertResult) before testing
     * for containment. The algorithm uses these "chain starts" in order to clip the boundary of A
     * to the interior of B in an output-sensitive way.
     *
     * <p>This method supports "early exit" in the case where a boolean predicate is being evaluated
     * and the algorithm discovers that the result region will be non-empty.
     */
    private boolean getChainStarts(
        int aRegionId,
        boolean invertA,
        boolean invertB,
        boolean invertResult,
        CrossingProcessor crossingProcessor,
        List<ShapeEdgeId> chainStarts) {
      S2ShapeIndex aIndex = op.regions[aRegionId];
      S2ShapeIndex bIndex = op.regions[1 - aRegionId];

      if (isBooleanOutput) {
        // If boolean output is requested, then we use the CrossingProcessor to determine whether
        // the first edge of each chain will be emitted to the output region. This lets us terminate
        // the operation early in many cases.
        crossingProcessor.startBoundary(aRegionId, invertA, invertB, invertResult);
      }

      // If region B has no two-dimensional shapes and is not inverted, then by definition no chain
      // starts are contained. However, if boolean output is requested then we check for containment
      // anyway, since as a side effect we may discover that the result region is non-empty and
      // terminate the entire operation early.
      boolean bHasInterior = S2ShapeUtil.hasInterior(bIndex);
      if (bHasInterior || invertB || isBooleanOutput) {
        S2ContainsPointQuery bQuery = new S2ContainsPointQuery(bIndex);
        MutableEdge reusableEdge = new MutableEdge();

        int numShapeIds = aIndex.getShapes().size();
        for (int shapeId = 0; shapeId < numShapeIds; ++shapeId) {
          S2Shape aShape = aIndex.getShapes().get(shapeId);

          // If region A is being subtracted from region B, points and polylines in region A can be
          // ignored since these shapes never contribute to the output (they can only remove edges
          // from region B).
          if (invertA != invertResult && aShape.dimension() < 2) {
            continue;
          }
          if (isBooleanOutput) {
            crossingProcessor.startShape(aShape);
          }

          for (int chainId = 0; chainId < aShape.numChains(); ++chainId) {
            List<S2Point> chain = aShape.chain(chainId);
            if (chain.isEmpty()) {
              continue;
            }

            boolean inside = (bHasInterior && bQuery.contains(chain.get(0))) != invertB;
            int startEdgeId = aShape.getChainStart(chainId);
            if (inside) {
              chainStarts.add(new ShapeEdgeId(shapeId, startEdgeId));
            }
            if (isBooleanOutput) {
              crossingProcessor.startChain(chainId, aShape, inside);
              aShape.getEdge(startEdgeId, reusableEdge);
              ShapeEdge a =
                  new ShapeEdge(
                      shapeId, startEdgeId, reusableEdge.getStart(), reusableEdge.getEnd());
              if (!processIncidentEdges(a, bIndex, bQuery, crossingProcessor)) {
                return false;
              }
            }
          }
        }
      }

      chainStarts.add(ShapeEdgeId.SENTINEL);
      return true;
    }

    /**
     * Get edges from the given bIndex that are incident to the ShapeEdge 'a', and process those
     * with the given CrossingProcessors, adding them to the output if required.
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean processIncidentEdges(
        ShapeEdge a,
        S2ShapeIndex bIndex,
        S2ContainsPointQuery bQuery,
        CrossingProcessor crossingProcessor) {
      tmpCrossings.clear();
      // Find shape edges in bIndex that are incident to ShapeEdge 'a'.
      bQuery.visitIncidentEdges(
          a.v0,
          (bShapeId, bEdgeId, bStart, bEnd) -> {
            ShapeEdge b = new ShapeEdge(bShapeId, bEdgeId, bStart, bEnd);
            return addIndexCrossing(a, b, /* isInterior= */ false, tmpCrossings);
          },
          new MutableEdge());

      // Fast path for the common case where there are no incident edges. We return false
      // (terminating early) if the first chain edge will be emitted.
      if (tmpCrossings.isEmpty()) {
        return !crossingProcessor.inside();
      }
      // Otherwise we invoke the full CrossingProcessor logic to determine whether the first chain
      // edge will be emitted.
      if (tmpCrossings.size() > 1) {
        tmpCrossings.sort();
      }
      tmpCrossings.add(new IndexCrossing(ShapeEdgeId.SENTINEL, ShapeEdgeId.SENTINEL));
      CrossingIterator nextCrossing =
          new CrossingIterator(bIndex, tmpCrossings, /* crossingsComplete= */ false);
      return crossingProcessor.processEdge(a.id, nextCrossing);
    }

    private boolean addIndexCrossing(
        ShapeEdge a, ShapeEdge b, boolean isInterior, IndexCrossings crossings) {
      IndexCrossing crossing = new IndexCrossing(a.id, b.id);
      crossings.add(crossing);

      if (isInterior) {
        crossing.isInteriorCrossing = true;
        if (S2Predicates.sign(a.v0, a.v1, b.v0) > 0) {
          crossing.leftToRight = true;
        }
        builder.addIntersection(S2EdgeUtil.getIntersection(a.v0, a.v1, b.v0, b.v1));
      } else {
        // TODO(ericv): This field isn't used unless one shape is a polygon and the other is a
        // polyline or polygon, but we don't have the shape dimension information readily available
        // here.
        if (S2EdgeUtil.vertexCrossing(a.v0, a.v1, b.v0, b.v1)) {
          crossing.isVertexCrossing = true;
        }
      }
      return true; // Continue visiting.
    }

    /**
     * Initialize indexCrossings to the set of crossing edge pairs such that the first element of
     * each pair is an edge from "regionId".
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean initIndexCrossings(int regionId) {
      if (regionId == indexCrossingsFirstRegionId) {
        // Already initialized.
        return true;
      }
      if (indexCrossingsFirstRegionId < 0) {
        Preconditions.checkState(regionId == 0); // For efficiency, not correctness.

        // TODO(torrey): ericv noted this would be more efficient (also in C++) if
        // visitCrossingEdgePairs() returned the sign (+1 or -1) of the interior crossing, i.e.
        // "int interiorCrossingSign" rather than "boolean isInterior".
        S2CrossingEdgesQuery query = new S2CrossingEdgesQuery(CrossingType.ALL);
        if (!query.visitCrossingEdgePairs(
            op.regions[0],
            op.regions[1],
            (aShapeId, aEdgeId, aEdgeSrc, aEdgeDst, // from regions[0]
             bShapeId, bEdgeId, bEdgeSrc, bEdgeDst, // from regions[1]
             isInterior) -> {
              // For all supported operations (union, intersection, and difference), if the
              // input edges have an interior crossing then the output is guaranteed to have at
              // least one edge.
              if (isInterior && isBooleanOutput) {
                return false;
              }
              // TODO: addIndexCrossing should take the shapeid and edge ids and edge endpoints
              // directly, avoiding building ShapeEdge objects.
              ShapeEdge a = new ShapeEdge(aShapeId, aEdgeId, aEdgeSrc, aEdgeDst);
              ShapeEdge b = new ShapeEdge(bShapeId, bEdgeId, bEdgeSrc, bEdgeDst);
              return addIndexCrossing(a, b, isInterior, indexCrossings);
            })) {
          return false;
        }
        if (indexCrossings.size() > 1) {
          indexCrossings.sort(); // Don't need to remove duplicates as the iterator skips them.
        }
        // Add a sentinel value to simplify the loop logic.
        indexCrossings.add(new IndexCrossing(ShapeEdgeId.SENTINEL, ShapeEdgeId.SENTINEL));
        indexCrossingsFirstRegionId = 0;
      }

      if (regionId != indexCrossingsFirstRegionId) {
        indexCrossings.reverseAll();
        indexCrossings.sort();
        indexCrossingsFirstRegionId = regionId;
      }
      return true;
    }

    /**
     * Adds the boundary of each region to the output, clipped to the "inside" of the other region.
     * Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty. Otherwise, adds edges to the output S2Builder.
     *
     * <p>addBoundaryPair is the core of the S2BooleanOperation implementation. Each of the specific
     * boolean operations (UNION, INTERSECTION, ...) is implemented as a single call to
     * addBoundaryPair with the appropriate invertA, invertB, invertResult flags, except for
     * SYMMETRIC_DIFFERENCE, which requires two calls.
     */
    private boolean addBoundaryPair(
        boolean invertA,
        boolean invertB,
        boolean invertResult,
        CrossingProcessor crossingProcessor) {
      // Optimization: if the operation is DIFFERENCE or SYMMETRIC_DIFFERENCE, it is worthwhile
      // checking whether the two regions are identical (in which case the output is empty).
      OpType type = op.opType;
      if (type == OpType.DIFFERENCE || type == OpType.SYMMETRIC_DIFFERENCE) {
        if (areRegionsIdentical()) {
          return true;
        }
      } else if (isBooleanOutput) {
        // TODO(ericv): When boolean output is requested there are other quick checks that could be
        // done here, such as checking whether a full cell from one S2ShapeIndex intersects a non-
        // empty cell of the other S2ShapeIndex.
      }

      // TODO(torrey): Better data structures? ShapeEdgeId is a pair of ints.
      List<ShapeEdgeId> aStarts = new ArrayList<>();
      List<ShapeEdgeId> bStarts = new ArrayList<>();

      // For INTERSECTION, where invertA, invertB, and invertResult are all false:
      // 1. Get the chain starts for region A: the first edge of each edge chain in region A whose
      //    first vertex is contained by region B's polygons. These edges will be part of the
      //    boundary of the intersection of the two regions, because they are inside one and on
      //    the boundary of the other.
      // 2. Similarly, get the chain starts for region B: the first edge of each edge chain in
      //    region B whose first vertex is contained by region A's polygons. These edges are also
      //    inside one region and on the boundary of the other, so will be part of the boundary of
      //    the intersection.
      // 3. If there are no chain starts, the intersection is empty. Otherwise, walk the boundaries,
      //    clipping them to the interior of the other region, and adding the resulting edges to the
      //    output.
      //
      // For UNION, everything is inverted. So the chain starts are the first edges of chains that
      // are NOT contained by the other region. (Edges that are contained by the other region are
      // not required in the output). When walking boundaries, edges are clipped to the exterior of
      // the other region.
      //
      // For DIFFERENCE, region B is inverted to subtract it from region A. So the chain starts are
      // the first edges of edge chains of A that are NOT contained by region B, so not subtracted,
      // and the first edges of edge chains of B that ARE contained by region A, and are subtracting
      // from A. When walking chains of A, they are clipped to the exterior of B, and when walking
      // chains of B, they are clipped to the interior of A.

      if (!getChainStarts(0, invertA, invertB, invertResult, crossingProcessor, aStarts)
          || !getChainStarts(1, invertB, invertA, invertResult, crossingProcessor, bStarts)
          || !addBoundary(0, invertA, invertB, invertResult, aStarts, crossingProcessor)
          || !addBoundary(1, invertB, invertA, invertResult, bStarts, crossingProcessor)) {
        Preconditions.checkState(isBooleanOutput);
        return false;
      }
      if (!isBooleanOutput) {
        crossingProcessor.doneBoundaryPair();
      }
      return true;
    }

    /**
     * When subtracting regions, we can save a lot of work by detecting the relatively common case
     * where the two regions are identical.
     */
    private boolean areRegionsIdentical() {
      return S2ShapeUtil.equals(op.regions[0].getShapes(), op.regions[1].getShapes());
    }

    /**
     * Builds the specific boolean operation (UNION, INTERSECTION, ...) on the provided inputs,
     * adding edges to the output S2Builder if output is needed (not boolean results).
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean buildOpType(OpType opType) {
      // CrossingProcessor does the real work of emitting the output edges.
      CrossingProcessor crossingProcessor =
          new CrossingProcessor(
              op.options.polygonModel(),
              op.options.polylineModel(),
              op.options.polylineLoopsHaveBoundaries(),
              builder,
              inputDimensions,
              inputCrossings);
      switch (opType) {
        case UNION:
          // A | B == ~(~A & ~B)
          return addBoundaryPair(true, true, true, crossingProcessor);

        case INTERSECTION:
          // A & B
          return addBoundaryPair(false, false, false, crossingProcessor);

        case DIFFERENCE:
          // A - B = A & ~B
          //
          // Note that degeneracies are implemented such that the symmetric operation (-B + A) also
          // produces correct results. This can be tested by swapping op.regions[0, 1] and calling
          // addBoundaryPair(true, false, false), which computes (~B & A).
          return addBoundaryPair(false, true, false, crossingProcessor);

        case SYMMETRIC_DIFFERENCE:
          // Compute the union of (A - B) and (B - A).
          return (addBoundaryPair(false, true, false, crossingProcessor)
              && addBoundaryPair(true, false, false, crossingProcessor));
      }
      throw new IllegalStateException("Invalid S2BooleanOperation.OpType" + opType);
    }

    /** Returns a bit mask indicating which of the 6 S2 cube faces intersect the index contents. */
    private static byte getFaceMask(S2ShapeIndex index) {
      byte mask = 0;
      S2Iterator<S2ShapeIndex.Cell> it = index.iterator();
      while (!it.done()) {
        int face = it.id().face();
        mask = (byte) (mask | 1 << face);
        it.seek(S2CellId.fromFace(face + 1).rangeMin());
      }
      return mask;
    }

    /**
     * Given a polygon edge graph containing only degenerate edges and sibling edge pairs, the
     * purpose of this function is to decide whether the polygon is empty or full except for the
     * degeneracies, i.e. whether the degeneracies represent shells or holes. The given graph may be
     * null, which indicates a graph with no edges.
     */
    private boolean isFullPolygonResult(S2BuilderGraph unused) {
      // TODO(user): Refactor isFullPolygonResult to be an S2Builder.IsFullPolygonPredicate.
      //
      // If there are no edges of dimension 2, the result could be either the empty polygon or the
      // full polygon. Note that this is harder to determine than you might think due to snapping.
      // For example, the union of two non-empty polygons can be empty, because both polygons
      // consist of tiny loops that are eliminated by snapping. Similarly, even if two polygons both
      // contain a common point their intersection can still be empty.
      //
      // We distinguish empty from full results using two heuristics:
      //
      //  1. We compute a bit mask representing the subset of the six S2 cube faces intersected by
      //     each input geometry, and use this to determine if only one of the two results is
      //     possible. (This test is very fast.) Note that snapping will never cause the result to
      //     cover an entire extra cube face because the maximum allowed snap radius is too small,
      //     (it is checked below to be at most 70 degrees).
      //
      //  2. We compute the area of each input geometry, and use this to bound the minimum and
      //     maximum area of the result. If only one of {0, 4*Pi} is possible then we are done. If
      //     neither is possible then we choose the one that is closest to being possible (since
      //     snapping can change the result area). Both results are possible only when computing the
      //     symmetric difference of two regions of area 2*Pi each, in which case we must resort to
      //     additional heuristics (see below).
      //
      // TODO(ericv): Implement a predicate that uses the results of edge snapping directly, rather
      // than computing areas. This would not only be much faster but would also allows all cases
      // to be handled 100% robustly.

      assert S2BuilderSnapFunctions.maxSnapRadius().degrees() <= 70;

      S2ShapeIndex a = op.regions[0];
      S2ShapeIndex b = op.regions[1];
      switch (op.opType) {
        case UNION:
          return isFullPolygonUnion(a, b);

        case INTERSECTION:
          return isFullPolygonIntersection(a, b);

        case DIFFERENCE:
          return isFullPolygonDifference(a, b);

        case SYMMETRIC_DIFFERENCE:
          return isFullPolygonSymmetricDifference(a, b);
      }
      throw new IllegalStateException("Invalid S2BooleanOperation.OpType");
    }

    private boolean isFullPolygonUnion(S2ShapeIndex a, S2ShapeIndex b) {
      // See comments in isFullPolygonResult(). The most common case is that neither input polygon
      // is empty but the result is empty due to snapping.

      // The result can be full only if the union of the two input geometries intersects all six
      // faces of the S2 cube. This test is fast.
      if ((getFaceMask(a) | getFaceMask(b)) != ALL_FACES_MASK) {
        return false;
      }

      // The union area satisfies:
      //
      //   max(A, B) <= Union(A, B) <= min(4*Pi, A + B)
      //
      // where A, B can refer to a polygon or its area. We then choose the result that assumes the
      // smallest amount of error.
      double aArea = S2ShapeIndexMeasures.area(a);
      double bArea = S2ShapeIndexMeasures.area(b);
      double minArea = max(aArea, bArea);
      double maxArea = min(4 * PI, aArea + bArea);
      return minArea > 4 * PI - maxArea;
    }

    private boolean isFullPolygonIntersection(S2ShapeIndex a, S2ShapeIndex b) {
      // See comments in isFullPolygonResult(). By far the most common case is that the result is
      // empty.

      // The result can be full only if each of the two input geometries intersects all six faces of
      // the S2 cube. This test is fast.
      if ((getFaceMask(a) & getFaceMask(b)) != ALL_FACES_MASK) {
        return false;
      }

      // The intersection area satisfies:
      //
      //   max(0, A + B - 4*Pi) <= Intersection(A, B) <= min(A, B)
      //
      // where A, B can refer to a polygon or its area. We then choose the result that assumes the
      // smallest amount of error.
      double aArea = S2ShapeIndexMeasures.area(a);
      double bArea = S2ShapeIndexMeasures.area(b);
      double minArea = max(0.0, aArea + bArea - 4 * PI);
      double maxArea = min(aArea, bArea);
      return minArea > 4 * PI - maxArea;
    }

    private boolean isFullPolygonDifference(S2ShapeIndex a, S2ShapeIndex b) {
      // See comments in isFullPolygonResult(). By far the most common case is that the result is
      // empty.

      // The result can be full only if each cube face is intersected by the first geometry. (The
      // second geometry is irrelevant, since for example it could consist of a tiny loop on each S2
      // cube face.)  This test is fast.
      if (getFaceMask(a) != ALL_FACES_MASK) {
        return false;
      }

      // The difference area satisfies:
      //
      //   max(0, A - B) <= Difference(A, B) <= min(A, 4*Pi - B)
      //
      // where A, B can refer to a polygon or its area. We then choose the result that assumes the
      // smallest amount of error.
      double aArea = S2ShapeIndexMeasures.area(a);
      double bArea = S2ShapeIndexMeasures.area(b);
      double minArea = max(0.0, aArea - bArea);
      double maxArea = min(aArea, 4 * PI - bArea);
      return minArea > 4 * PI - maxArea;
    }

    private boolean isFullPolygonSymmetricDifference(S2ShapeIndex a, S2ShapeIndex b) {
      // See comments in isFullPolygonResult(). By far the most common case is that the result is
      // empty.

      // The result can be full only if the union of the two input geometries intersects all six
      // faces of the S2 cube. This test is fast.
      byte aMask = getFaceMask(a);
      byte bMask = getFaceMask(b);
      if ((aMask | bMask) != ALL_FACES_MASK) {
        return false;
      }

      // The symmetric difference area satisfies:
      //
      //   |A - B| <= SymmetricDifference(A, B) <= 4*Pi - |4*Pi - (A + B)|
      //
      // where A, B can refer to a polygon or its area.
      double aArea = S2ShapeIndexMeasures.area(a);
      double bArea = S2ShapeIndexMeasures.area(b);
      double minArea = abs(aArea - bArea);
      double maxArea = 4 * PI - abs(4 * PI - (aArea + bArea));

      // Now we choose the result that assumes the smallest amount of error (minArea in the empty
      // case, and (4*Pi - maxArea) in the full case). However in the case of symmetric difference
      // these two errors may be equal, meaning that the result is ambiguous. This happens when both
      // polygons have area 2*Pi. Furthermore, this can happen even when the areas are not exactly
      // 2*Pi due to snapping and area calculation errors.
      //
      // To determine whether the result is ambiguous, we compute a rough estimate of the maximum
      // expected area error (including errors due to snapping), using the worst-case error bound
      // for a hemisphere defined by 4 vertices.
      S1Angle edgeSnapRadius = builderOptions.edgeSnapRadius();
      double hemisphereAreaError = 2 * PI * edgeSnapRadius.radians() + 40 * S2.DBL_EPSILON;

      // The following sign is the difference between the error needed for an empty result and the
      // error needed for a full result. It is negative if an empty result is possible, positive if
      // a full result is possible, and zero if both results are possible.
      double errorSign = minArea - (4 * PI - maxArea);
      if (abs(errorSign) <= hemisphereAreaError) {
        // Handling the ambiguous case correctly requires a more sophisticated algorithm (see
        // below), but we can at least handle the simple cases by testing whether both input
        // geometries intersect all 6 cube faces. If not, then the result is definitely full.
        if ((aMask & bMask) != ALL_FACES_MASK) {
          return true;
        }

        // Otherwise both regions have area 2*Pi and intersect all 6 cube faces. We choose "empty"
        // in this case under the assumption that it is more likely that the user is computing the
        // difference between two nearly identical polygons.
        //
        // TODO(ericv): Implement a robust algorithm based on examining the edge snapping results
        // directly, or alternatively add another heuristic (such as testing containment of random
        // points, or using a larger bit mask in the tests above, e.g. a 24-bit mask representing
        // all level 1 cells).
        return false;
      }
      return errorSign > 0;
    }

    /**
     * Creates S2Builder options, and an S2Builder with an EdgeClippingLayer if required (i.e.
     * isBooleanOutput is false). Then runs buildOpType. If an error occurs, the given error is set.
     */
    private void doBuild(S2Error error) {
      builderOptions = new S2Builder.Builder(op.options.snapFunction());
      builderOptions.setIntersectionTolerance(S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR));
      if (op.options.splitAllCrossingPolylineEdges()) {
        builderOptions.setSplitCrossingEdges(true);
      }
      // TODO(ericv): Ideally idempotent() should be true, but existing clients expect vertices
      // closer than the full "snapRadius" to be snapped.
      builderOptions.setIdempotent(false);

      if (isBooleanOutput) {
        // buildOpType() returns true if and only if the result has no edges.
        op.resultEmpty = buildOpType(op.opType) && !isFullPolygonResult(null);
        return;
      }

      builder = builderOptions.build();
      builder.startLayer(new EdgeClippingLayer(op.layers, inputDimensions, inputCrossings));

      // Add a predicate that decides whether a result with no polygon edges should be interpreted
      // as the empty polygon or the full polygon.
      builder.addIsFullPolygonPredicate(this::isFullPolygonResult);

      boolean unused = buildOpType(op.opType);
      unused = builder.build(error);
    }
  }

  /**
   * PointCrossingResult describes the relationship between a point from region A and a set of
   * crossing edges from region B. For example, "matchesPolygon" indicates whether a polygon vertex
   * from region B matches the given point.
   */
  private static class PointCrossingResult {
    // Note that "matchesPolyline" is true only if the point matches a polyline vertex of B *and*
    // the polyline contains that vertex, whereas "matchesPolygon" is true if the point matches any
    // polygon vertex.

    /** The point from region A matches a point from region B. */
    boolean matchesPoint = false;

    /** The point from revion A matches a _contained_ vertex of a polyline from region B. */
    boolean matchesPolyline = false;

    /** The point from region A matches a polygon vertex from region B. */
    boolean matchesPolygon = false;
  }

  /**
   * A helper class for iterating through the edges from region B that cross a particular edge from
   * region A. It caches information from the current shape, chain, and edge so that it doesn't need
   * to be looked up repeatedly. Typical usage:
   *
   * {@snippet :
   * void someFunction(ShapeEdgeId aId, CrossingIterator it) {
   *   // Iterate through the edges that cross edge "aId".
   *   for (; !it.done(aId); it.next()) {
   *     ... use it.bShape(), it.bEdge(), etc ...
   *   }
   * }
   * }
   */
  private static class CrossingIterator {
    private final S2ShapeIndex bIndex;
    private final IndexCrossings.Iterator it;
    private S2Shape bShape;
    private int bShapeId;
    private int bDimension;
    private final boolean crossingsComplete;

    // The chain Id to which the bEdge belongs. Computed on demand and cached.
    private int bEdgeChainId;

    /**
     * Creates an iterator over crossing edge pairs (a, b) where "b" is an edge from "bIndex".
     * "crossingsComplete" indicates that "crossings" contains all edge crossings between the two
     * regions (rather than a subset).
     */
    public CrossingIterator(
        S2ShapeIndex bIndex, IndexCrossings crossings, boolean crossingsComplete) {
      this.bIndex = bIndex;
      this.it = crossings.iterator(); // Starts at the first element, or SENTINEL if empty.
      bShapeId = -1;
      this.crossingsComplete = crossingsComplete;
      update();
    }

    /** Advances to the next crossing edge pair. */
    public void next() {
      it.next();
      update();
    }

    /**
     * Returns true if the current aId() is not equal to the given id, i.e. if we have been
     * processing edges from the given id, we are now finished.
     */
    public boolean done(ShapeEdgeId id) {
      return !aId().isEqualTo(id);
    }

    /** True if all edge crossings are available (see above). */
    public boolean crossingsComplete() {
      return crossingsComplete;
    }

    /** True if this crossing occurs at a point interior to both edges. */
    public boolean isInteriorCrossing() {
      return it.current.isInteriorCrossing;
    }

    /**
     * Equal to S2EdgeUtil.vertexCrossing(aEdge, bEdge), provided that aEdge and bEdge have exactly
     * one vertex in common and neither edge is degenerate.
     */
    public boolean isVertexCrossing() {
      return it.current.isVertexCrossing;
    }

    /** True if aEdge crosses bEdge from left to right (for interior crossings). */
    public boolean leftToRight() {
      return it.current.leftToRight;
    }

    /** Returns the S2ShapeIndex edge (shape id and edge id) from region A. */
    public ShapeEdgeId aId() {
      return it.current.aEdge;
    }

    /** Returns the S2ShapeIndex edge (shape id and edge id) from region B. */
    public ShapeEdgeId bId() {
      return it.current.bEdge;
    }

    /** Returns the S2ShapeIndex of region B. */
    public S2ShapeIndex bIndex() {
      return bIndex;
    }

    /** Returns the S2Shape from region B that contains the bEdge. */
    public S2Shape bShape() {
      return bShape;
    }

    /** Returns the dimension of the S2Shape from region B that contains the bEdge. */
    public int bDimension() {
      return bDimension;
    }

    /** Returns the id of the S2Shape from region B that contains the bEdge. */
    public int bShapeId() {
      return bShapeId;
    }

    /** Returns the edge id of the bEdge within its shape. */
    public int bEdgeId() {
      return bId().edgeId;
    }

    /** Fills the given MutableEdge with the endpoints of the bEdge. */
    public void getBEdge(MutableEdge edge) {
      bShape.getEdge(bEdgeId(), edge); // Opportunity to cache this.
    }

    /** Returns the chain id of the bEdge within its shape. */
    public int bEdgeChainId() {
      if (bEdgeChainId < 0) {
        // TODO(torrey): Look at how this code is used, consider improvements, avoid allocation
        S2Shape.ChainPosition position = new S2Shape.ChainPosition();
        bShape.getChainPosition(bEdgeId(), position);
        bEdgeChainId = position.chainId;
      }
      return bEdgeChainId;
    }

    /** Updates information about the B shape whenever it changes. */
    private void update() {
      if (!aId().isEqualTo(ShapeEdgeId.SENTINEL) && bId().shapeId != bShapeId) {
        bShapeId = bId().shapeId;
        bShape = bIndex.getShapes().get(bShapeId);
        bDimension = bShape.dimension();
        bEdgeChainId = -1; // Computed on demand.
      }
    }
  }

  /**
   * An IndexCrossing represents a pair of intersecting S2ShapeIndex edges ("aEdge" and "bEdge"). We
   * store all such intersections because the algorithm needs them twice, once when processing the
   * boundary of region A and once when processing the boundary of region B.
   */
  private static class IndexCrossing {
    private static final IndexCrossing SENTINEL =
        new IndexCrossing(ShapeEdgeId.SENTINEL, ShapeEdgeId.SENTINEL);

    /** The intersecting shape index edge from region A. */
    private ShapeEdgeId aEdge;
    /** The intersecting shape index edge from region B. */
    private ShapeEdgeId bEdge;

    // TODO(torrey): C++ packs these three booleans into one int.

    /** True if S2.crossingSign(aEdge, bEdge) > 0. */
    boolean isInteriorCrossing;

    /**
     * True if "aEdge" crosses "bEdge" from left to right. Undefined if isInteriorCrossing is
     * false.
     */
    boolean leftToRight;

    /**
     * Equal to S2.vertexCrossing(aEdge, bEdge). Undefined if "aEdge" and "bEdge" do not share
     * exactly one vertex or either edge is degenerate.
     */
    boolean isVertexCrossing;

    /**
     * Constructs a new IndexCrossing. All flags are "false" by default.
     */
    public IndexCrossing(ShapeEdgeId a, ShapeEdgeId b) {
      this.aEdge = a;
      this.bEdge = b;
      this.isInteriorCrossing = false;
      this.leftToRight = false;
      this.isVertexCrossing = false;
    }

    /** Reverses the crossing by swapping aEdge and bEdge. */
    public void reverse() {
      ShapeEdgeId tmp = aEdge;
      aEdge = bEdge;
      bEdge = tmp;
      // The following predicates get inverted when the edges are swapped.
      leftToRight ^= true;
      isVertexCrossing ^= true;
    }

    @Override
    public int hashCode() {
      return aEdge.hashCode() * 7 + bEdge.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
      if (!(obj instanceof IndexCrossing)) {
        return false;
      }
      IndexCrossing other = (IndexCrossing) obj;
      return this.isEqualTo(other);
    }

    public boolean isEqualTo(IndexCrossing other) {
      return this.aEdge.shapeId == other.aEdge.shapeId
          && this.aEdge.edgeId == other.aEdge.edgeId
          && this.bEdge.shapeId == other.bEdge.shapeId
          && this.bEdge.edgeId == other.bEdge.edgeId;
    }
  }

  /**
   * A list of IndexCrossings that may be sorted and iterated. Each IndexCrossing contains two
   * ShapeEdgeIds and details of type type of crossing.
   */
  private static class IndexCrossings {
    // TODO(torrey): This is likely worth optimizing.
    private final ArrayList<IndexCrossing> list = new ArrayList<>();

    public static final Comparator<IndexCrossing> COMPARATOR =
        (x, y) -> {
          int result = Integer.compare(x.aEdge.shapeId, y.aEdge.shapeId);
          if (result != 0) {
            return result;
          }
          result = Integer.compare(x.aEdge.edgeId, y.aEdge.edgeId);
          if (result != 0) {
            return result;
          }
          result = Integer.compare(x.bEdge.shapeId, y.bEdge.shapeId);
          if (result != 0) {
            return result;
          }
          return Integer.compare(x.bEdge.edgeId, y.bEdge.edgeId);
        };

    /**
     * An iterator over IndexCrossings that skips duplicate entries. Users access the 'current'
     * entry directly.
     */
    public class Iterator {
      // The position of the current element.
      private int position = 0;
      // The current IndexCrossing, or SENTINEL if no more are available.
      public IndexCrossing current = null;

      /** Constructs a new Iterator positioned on the first IndexCrossing. */
      public Iterator() {
        if (list.isEmpty()) {
          current = IndexCrossing.SENTINEL;
        } else {
          current = list.get(0);
        }
      }

      /** Returns true if there is another IndexCrossing. */
      public boolean hasNext() {
        return !current.isEqualTo(IndexCrossing.SENTINEL);
      }

      /** Advances to the next IndexCrossing. */
      public void next() {
        assert hasNext();
        IndexCrossing last = list.get(position);
        position++;
        while (last.isEqualTo(list.get(position)) && position < list.size()) {
          position++;
        }
        current = (position == list.size()) ? IndexCrossing.SENTINEL : list.get(position);
      }
    }

    /** Returns an IndexCrossings.Iterator over this list of IndexCrossings. */
    public IndexCrossings.Iterator iterator() {
      return new IndexCrossings.Iterator();
    }

    /** Clears this list of IndexCrossings. */
    public void clear() {
      list.clear();
    }

    /** Reverses all the IndexCrossings in this list. */
    public void reverseAll() {
      for (IndexCrossing crossing : list) {
        crossing.reverse();
      }
    }

    /** Sorts this list of IndexCrossings according to the IndexCrossings.COMPARATOR. */
    public void sort() {
      list.sort(COMPARATOR);
    }

    public void add(IndexCrossing crossing) {
      list.add(crossing);
    }

    public int size() {
      return list.size();
    }

    public boolean isEmpty() {
      return list.isEmpty();
    }
  }

  /**
   * CrossingProcessor processes edge crossings and maintains the necessary state in order to clip
   * the boundary of one region to the interior or exterior of the other region. It processes all
   * the edges from one region that cross a specific edge of the other region, outputs the
   * appropriate edges to an S2Builder, and outputs other information required by GraphEdgeClipper
   * to the given vectors.
   */
  private static class CrossingProcessor {
    // Constructor parameters:

    private final PolygonModel polygonModel;
    private final PolylineModel polylineModel;
    private final boolean polylineLoopsHaveBoundaries;

    /**
     * The output of the CrossingProcessor consists of a subset of the input edges that are emitted
     * to "builder", and some auxiliary information that allows GraphEdgeClipper to determine which
     * segments of those input edges belong to the output. The auxiliary information consists of the
     * dimension of each input edge, and set of input edges from the other region that cross each
     * input edge. These output fields are provided to the constructor. This builder is null if
     * boolean output was requested.
     */
    private final S2Builder builder;
    /** The dimension of each input edge. */
    private final IntVector inputDimensions;
    /** The set of input edges from the other region that cross each input edge. */
    private final InputEdgeCrossings inputCrossings;

    // Fields set by startBoundary:

    private int aRegionId;
    private int bRegionId;
    private boolean invertA;
    private boolean invertB;
    private boolean invertResult;
    /** True if this is a UNION operation. */
    private boolean isUnion;

    // Fields set by startShape:

    private S2Shape aShape;
    private int aDimension;

    // Fields set by startChain:

    private int chainId;
    private int chainStart;
    private int chainLimit;

    // Fields updated by processEdge:

    private final SourceEdgeCrossings sourceEdgeCrossings = new SourceEdgeCrossings();

    /**
     * A set of edges that cross the current edge being processed by processEdge() but that have not
     * yet been associated with a particular S2Builder edge. This is necessary because processEdge
     * can create up to three S2Builder edges per input edge: one to represent the edge interior,
     * and up to two more to represent an isolated start and/or end vertex. The crossing edges must
     * be associated with the S2Builder edge that represents the edge interior, and they are stored
     * here until that edge is created.
     */
    private final List<SourceEdgeCrossing> pendingSourceEdgeCrossings = new ArrayList<>();

    // TODO(torrey): Consider alternatives to HashMap. C++ uses
    // absl::btree_map<SourceId, InputEdgeId>;

    /**
     * A map that translates from SourceId (the (regionId, shapeId, edgeId) triple that identifies
     * an S2ShapeIndex edge, encoded into a Long) to InputEdgeId (the sequentially increasing
     * numbers assigned to input edges by S2Builder).
     */
    private static class SourceIdMap extends HashMap<Long, Integer> {}
    private final SourceIdMap sourceIdMap = new SourceIdMap();

    // TODO(torrey): Consider alternatives to HashMap. C++ uses a flat_hash_map<ShapeEdgeId, bool>.
    // A ShapeEdgeId is just two ints.

    /**
     * For each edge in region B that defines a degenerate loop (either a point loop or a sibling
     * pair), indicates whether that loop represents a shell or a hole. This information is used
     * during the second pass of addBoundaryPair() to determine the output for degenerate edges.
     */
    private final Map<ShapeEdgeId, Boolean> isDegenerateHole = new HashMap<>();

    /**
     * Indicates whether the point being processed along the current edge chain is in the polygonal
     * interior of the opposite region, using semi-open boundaries. If "invertB" is true then this
     * field is inverted.
     *
     * <p>Equal to: bIndex.contains(current point) ^ invertB
     */
    private boolean inside;

    /**
     * The value that "inside" would have just before the end of the previous edge added to
     * S2Builder. This value is used to determine whether the GraphEdgeClipper state needs to be
     * updated when jumping from one edge chain to another.
     */
    private boolean prevInside;

    /**
     * The maximum edge id of any edge in the current chain whose v0 vertex has already been
     * emitted. This is used to determine when an isolated vertex needs to be emitted, e.g. when two
     * closed polygons share only a vertex.
     */
    private int v0EmittedMaxEdgeId;

    /**
     * True if the first vertex of the current chain has been emitted. This is used when processing
     * loops in order to determine whether the first/last vertex of the loop should be emitted as an
     * isolated vertex.
     */
    private boolean chainV0Emitted;

    /**
     * Prepares to build output for the given polygon and polyline boundary models. Edges are
     * emitted to "builder", while other auxiliary data is appended to the given vectors.
     *
     * <p>If a predicate is being evaluated (i.e., we do not need to construct the actual result),
     * then "builder" and the various output vectors should all be null.
     */
    public CrossingProcessor(
        PolygonModel polygonModel,
        PolylineModel polylineModel,
        boolean polylineLoopsHaveBoundaries,
        S2Builder builder,
        IntVector inputDimensions,
        InputEdgeCrossings inputCrossings) {
      this.polygonModel = polygonModel;
      this.polylineModel = polylineModel;
      this.polylineLoopsHaveBoundaries = polylineLoopsHaveBoundaries;
      this.builder = builder;
      this.inputDimensions = inputDimensions;
      this.inputCrossings = inputCrossings;
      this.prevInside = false;
    }

    /**
     * Starts processing edges from the given region. "invertA", "invertB", and "invertResult"
     * indicate whether region A, region B, and/or the result should be inverted, which allows
     * operations such as union and difference to be implemented. For example, union is ~(~A & ~B).
     *
     * <p>This method should be called in pairs, once to process the edges from region A and once to
     * process the edges from region B.
     */
    public void startBoundary(
        int aRegionId, boolean invertA, boolean invertB, boolean invertResult) {
      this.aRegionId = aRegionId;
      this.bRegionId = 1 - aRegionId;
      this.invertA = invertA;
      this.invertB = invertB;
      this.invertResult = invertResult;
      this.isUnion = invertB && invertResult;

      // Specify to GraphEdgeClipper how these edges should be clipped.
      setClippingState(SET_REVERSE_A, invertA != invertResult);
      setClippingState(SET_INVERT_B, invertB);
    }

    /** Starts processing edges from the given shape. */
    public void startShape(S2Shape aShape) {
      this.aShape = aShape;
      aDimension = aShape.dimension();
    }

    /**
     * Starts processing edges from the given chain. If 'inside', the first vertex of the chain is
     * "inside" region b (after allowing for inversions).
     */
    public void startChain(int chainId, S2Shape shape, boolean inside) {
      this.chainId = chainId;
      // TODO(torrey): If this shows up in profiling, consider how startChain is called and perhaps
      // find ways to do less work here. In C++, chain.start and chain.length are members of the
      // Chain object, so the next two lines are very fast. In Java, implementations of
      // getChainStart() are constant time except for S2Polygon, which requires iterating over the
      // chains. Implementations of getChainLength are constant time.
      this.chainStart = shape.getChainStart(chainId);
      this.chainLimit = chainStart + shape.getChainLength(chainId);
      this.inside = inside;
      this.v0EmittedMaxEdgeId = chainStart - 1; // No edges emitted yet.
      this.chainV0Emitted = false;
    }

    // Allocated once here and used only in processEdge().
    private final MutableEdge edgeA = new MutableEdge();

    /**
     * Processes the given edge "aId" of the current chain. The iterator "it" should be positioned
     * to the set of edges from the other region that cross "aId" (if any). The iterator will be
     * advanced over all those edges.
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    public boolean processEdge(ShapeEdgeId aId, CrossingIterator it) {
      // Get the endpoints of edge aId.
      aShape.getChainEdge(chainId, aId.edgeId - chainStart, edgeA);
      if (aDimension == 0) {
        return processEdge0(aId, edgeA, it);
      } else if (aDimension == 1) {
        return processEdge1(aId, edgeA, it);
      } else {
        assert aDimension == 2;
        return processEdge2(aId, edgeA, it);
      }
    }

    /**
     * Processes an edge of dimension 0 (i.e., a point) from region A.
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean processEdge0(ShapeEdgeId aId, MutableEdge a, CrossingIterator it) {
      assert a.getStart().equalsPoint(a.getEnd());
      // When a region is inverted, all points and polylines are discarded.
      if (invertA != invertResult) {
        // Advance the crossing iterator over all crossings with aId.
        skipCrossings(aId, it);
        return true;
      }
      // Advance the crossing iterator over all crossings with edge aId, and determine the
      // relationship of the starting point of edge aId, and the shapes in region B.
      PointCrossingResult r = processPointCrossings(aId, a.getStart(), it);

      // "contained" indicates whether the current point is inside the polygonal interior of the
      // opposite region, using semi-open boundaries.
      boolean contained = inside ^ invertB;
      if (r.matchesPolygon && polygonModel != PolygonModel.SEMI_OPEN) {
        contained = (polygonModel == PolygonModel.CLOSED);
      }
      if (r.matchesPolyline) {
        contained = true;
      }

      // The output of UNION includes duplicate values, so ensure that points are not suppressed by
      // other points.
      if (r.matchesPoint && !isUnion) {
        contained = true;
      }

      // Test whether the point is contained after region B is inverted.
      if (contained == invertB) {
        return true; // Don't exit early.
      }
      return addPointEdge(a.getStart(), 0);
    }

    /**
     * Advances over any edge crossings with aId, as they are not needed to determine the result.
     * TODO(torrey): This could be a method on CrossingIterator, like: it.advancePast(aId)
     */
    private void skipCrossings(ShapeEdgeId aId, CrossingIterator it) {
      while (!it.done(aId)) {
        it.next();
      }
    }

    /**
     * Returns a summary of the relationship between a point from region A and a set of crossing
     * edges from region B (see PointCrossingResult). Advances the given CrossingIterator over all
     * crossings with edge aId. TODO(torrey): Only called from one place. Perhaps inline it to
     * avoid constructing and returning a PointCrossingResult?
     */
    private PointCrossingResult processPointCrossings(
        ShapeEdgeId aId, S2Point a0, CrossingIterator it) {
      PointCrossingResult r = new PointCrossingResult();
      for (; !it.done(aId); it.next()) {
        if (it.bDimension() == 0) {
          r.matchesPoint = true;
        } else if (it.bDimension() == 1) {
          if (polylineEdgeContainsVertex(a0, it, 0)) {
            r.matchesPolyline = true;
          }
        } else {
          r.matchesPolygon = true;
        }
      }
      return r;
    }

    /**
     * EdgeCrossingResult describes the relationship between an edge (a0, a1) from region A and a
     * set of crossing edges from region B. For example, "matchesPolygon" indicates whether (a0, a1)
     * matches a polygon edge from region B.
     */
    private static class EdgeCrossingResult {
      /** This field indicates that (a0, a1) exactly matches an edge of B, in either direction. */
      public boolean matchesPolyline = false;

      /**
       * This field indicates that a B polyline contains the degenerate polyline (a0, a0). This is
       * identical to whether the B polyline contains the point a0 except when the B polyline is
       * degenerate, since a degenerate polyline VV contains itself in all boundary models but
       * contains the point V only in the CLOSED polyline model.
       */
      public boolean a0MatchesPolyline = false;

      /**
       * This field indicates that a B polyline contains the degenerate polyline (a1, a1). As with
       * the previous field, this is identical to whether the B polyline contains the point a1
       * except when the B polyline is degenerate.
       */
      public boolean a1MatchesPolyline = false;

      /**
       * This field indicates that vertex a0 of (a0, a1) matches a polygon vertex of B. (Unlike with
       * polylines, the polygon may not contain that vertex.)
       */
      public boolean a0MatchesPolygon = false;

      /**
       * This field indicates that vertex a1 of (a0, a1) matches a polygon vertex of B. (Unlike with
       * polylines, the polygon may not contain that vertex.)
       */
      public boolean a1MatchesPolygon = false;

      /** When a0 != a1, identifies any B polygon edge that exactly matches (a0, a1). */
      public ShapeEdgeId polygonMatchId = ShapeEdgeId.NONE;

      /**
       * When a0 != a1, identifies any B polygon edge that exactly matches sibling edge (a1, a0).
       */
      public ShapeEdgeId siblingMatchId = ShapeEdgeId.NONE;

      /**
       * When a0 == a1, identifies any B polygon (degenerate) edge that exactly matches (a0, a0).
       */
      public ShapeEdgeId a0LoopMatchId = ShapeEdgeId.NONE;

      /** Convenience function to test whether a matching polygon edge was found. */
      public boolean matchesPolygon() {
        return polygonMatchId.edgeId >= 0;
      }

      /** Convenience function to test whether a matching sibling edge was found. */
      public boolean matchesSibling() {
        return siblingMatchId.edgeId >= 0;
      }

      /** Convenience function to test whether a matching degenerate polygon loop was found. */
      public boolean loopMatchesA0() {
        return a0LoopMatchId.edgeId >= 0;
      }

      // These fields count the number of edge crossings at a0, a1, and the interior of (a0, a1).

      /** Count of polygon crossings at a0. */
      public int a0Crossings = 0;

      /** Count of polygon crossings at a1. */
      public int a1Crossings = 0;

      /** Count of polygon crossings in edge interior. */
      public int interiorCrossings = 0;
    }

    /**
     * Processes an edge of dimension 1 (i.e., a polyline edge) from region A.
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean processEdge1(ShapeEdgeId aId, MutableEdge a, CrossingIterator it) {
      // When a region is inverted, all points and polylines are discarded.
      if (invertA != invertResult) {
        skipCrossings(aId, it);
        return true;
      }

      // Evaluate whether the start vertex should belong to the output, in case it needs to be
      // emitted as an isolated vertex.
      EdgeCrossingResult r = processEdgeCrossings(aId, a, it);
      boolean a0Inside = isPolylineVertexInside(r.a0MatchesPolyline, r.a0MatchesPolygon);

      // Test whether the entire polyline edge should be emitted (or not emitted) because it matches
      // a polyline or polygon edge.
      boolean isDegenerate = edgeIsPoint(a);
      inside ^= ((r.a0Crossings & 1) == 1);
      if (inside != isPolylineEdgeInside(r, isDegenerate)) {
        inside ^= true; // Invert the inside state.
        ++r.a1Crossings; // Restore the correct (semi-open) state later.
      }

      // If neither edge adjacent to v0 was emitted, and this polyline contains v0, and the other
      // region contains v0, then emit an isolated vertex.
      MutableEdge chainEdge = new MutableEdge();
      aShape.getChainEdge(chainId, chainLimit - chainStart - 1, chainEdge);
      if (!polylineLoopsHaveBoundaries
          && aId.edgeId == chainStart
          && a.getStart().equalsPoint(chainEdge.getEnd())) {
        // This is the first vertex of a polyline loop, so we can't decide if it is isolated until
        // we process the last polyline edge.
        chainV0Emitted = inside;
      } else if (isV0Isolated(aId)
          && !isDegenerate
          && polylineContainsV0(aId.edgeId, chainStart)
          && a0Inside) {
        if (!addPointEdge(a.getStart(), 1)) {
          return false;
        }
      }

      // Test whether the entire edge or any part of it belongs to the output.
      if (inside || r.interiorCrossings > 0) {
        // Note: updates "inside" to correspond to the state just before a1.
        if (!addEdge(aId, a.getStart(), a.getEnd(), /* dimension= */ 1, r.interiorCrossings)) {
          return false;
        }
      }
      // Remember whether the edge portion just before "a1" was emitted, so that
      // we can decide whether "a1" need to be emitted as an isolated vertex.
      if (inside) {
        v0EmittedMaxEdgeId = aId.edgeId + 1;
      }

      // Verify that edge crossings are being counted correctly.
      inside ^= (r.a1Crossings & 1) == 1;
      // This is a releatively expensive check and must only be done in debug mode. One of the
      // reasons this can fail is if one of the S2BooleanOperation inputs has overlapping polygon
      // interiors, which is not allowed.
      assert !it.crossingsComplete()
          || (new S2ContainsPointQuery(it.bIndex()).contains(a.getEnd()) == inside ^ invertB);

      // Special case to test whether the last vertex of a polyline should be emitted as an isolated
      // vertex.
      aShape.getChainEdge(chainId, chainStart, chainEdge);
      if (it.crossingsComplete()
          && !isDegenerate
          && isChainLastVertexIsolated(aId)
          && (polylineModel == PolylineModel.CLOSED
              || (!polylineLoopsHaveBoundaries && a.getEnd().equalsPoint(chainEdge.getStart())))
          && isPolylineVertexInside(r.a1MatchesPolyline, r.a1MatchesPolygon)) {
        if (!addPointEdge(a.getEnd(), 1)) {
          return false;
        }
      }
      return true;
    }

    /**
     * Returns true if the current point being processed (which must be a polyline vertex) is
     * contained by the opposite region (after inversion if "invertB" is true). "matchesPolyline"
     * and "matchesPolygon" indicate whether the vertex matches a polyline/ polygon vertex of the
     * opposite region.
     */
    private boolean isPolylineVertexInside(boolean matchesPolyline, boolean matchesPolygon) {
      // Initially "contained" indicates whether the current point is inside the polygonal interior
      // of region B using semi-open boundaries.
      boolean contained = inside ^ invertB;

      // For UNION the output includes duplicate polylines. The test below ensures that isolated
      // polyline vertices are not suppressed by other polyline vertices in the output.
      if (matchesPolyline && !isUnion) {
        contained = true;
      } else if (matchesPolygon && polygonModel != PolygonModel.SEMI_OPEN) {
        contained = (polygonModel == PolygonModel.CLOSED);
      }
      // Finally, invert the result if the opposite region should be inverted.
      return contained ^ invertB;
    }

    /**
     * Returns true if the current polyline edge is contained by the opposite region (after
     * inversion if "invertB" is true).
     */
    private boolean isPolylineEdgeInside(EdgeCrossingResult r, boolean isDegenerate) {
      // Initially "contained" indicates whether the current point (just past a0) is inside the
      // polygonal interior of region B using semi-open boundaries.
      boolean contained = inside ^ invertB;

      // Note that if r.matchesPolyline and isUnion is true, then "contained" will be false (unless
      // there is also a matching polygon edge) since polyline edges are not allowed in the interior
      // of B. In this case we leave "contained" as false since it causes both matching edges to be
      // emitted.
      if (r.matchesPolyline && !isUnion) {
        contained = true;
      } else if (isDegenerate) {
        // First allow the polygon boundary model to override the semi-open rules. Note that a
        // polygon vertex (dimension 2) is considered to completely contain degenerate OPEN and
        // SEMI_OPEN polylines (dimension 1) even though the latter do not contain any points. This
        // is because dimension 2 points are considered to be a strict superset of dimension 1
        // points.
        if (polygonModel != PolygonModel.SEMI_OPEN && r.a0MatchesPolygon) {
          contained = (polygonModel == PolygonModel.CLOSED);
        }
        // Note that r.a0MatchesPolyline is true if and only if some B polyline contains the
        // degenerate polyline (a0, a0).
        if (r.a0MatchesPolyline && !isUnion) {
          contained = true;
        }
      } else if (r.matchesPolygon()) {
        // In the SEMI_OPEN model, polygon sibling pairs cancel each other and have no effect on
        // point or edge containment.
        if (!(polygonModel == PolygonModel.SEMI_OPEN && r.matchesSibling())) {
          contained = (polygonModel != PolygonModel.OPEN);
        }
      } else if (r.matchesSibling()) {
        contained = (polygonModel == PolygonModel.CLOSED);
      }
      // Finally, invert the result if the opposite region should be inverted.
      return contained ^ invertB;
    }

    /**
     * Processes an edge of dimension 2 (i.e., a polygon edge) from region A.
     *
     * <p>Supports "early exit" in the case of boolean results by returning false as soon as the
     * result is known to be non-empty.
     */
    private boolean processEdge2(ShapeEdgeId aId, MutableEdge a, CrossingIterator it) {
      // Whenever the two regions contain the same edge, or opposite edges of a sibling pair, or one
      // region contains a point loop while the other contains a matching vertex, then in general
      // the result depends on whether one or both sides represent a degenerate shell or hole.
      //
      // In each pass it is easy to determine whether edges in region B represent degenerate
      // geometry, and if so whether they represent a shell or hole, since this can be determined
      // from the inside state and the matchesPolygon() / matchesSibling() methods of
      // EdgeCrossingResult. However, this information is not readily available for region A.
      //
      // We handle this by saving the shell/hole status of each degenerate loop in region B during
      // the first pass, and deferring the processing of any edges that meet the criteria above
      // until the second pass. (Note that regions A,B correspond to regions 0,1 respectively in the
      // first pass whereas they refer to regions 1,0 respectively in the second pass.)
      //
      // The first pass ignores:
      //  - degenerate edges of A that are incident to any edge of B
      //  - non-degenerate edges of A that match or are siblings to an edge of B
      //
      // The first pass also records the shell/hole status of:
      //  - degenerate edges of B that are incident to any edge of A
      //  - sibling pairs of B where either edge matches an edge of A
      //
      // The second pass processes and perhaps outputs:
      //  - degenerate edges of B that are incident to any edge of A
      //  - non-degenerate edges of B that match or are siblings to an edge of A
      //
      // The following flag indicates that we are in the second pass described above, i.e. that we
      // are emitting any necessary edges that were ignored by the first pass.
      boolean emitShared = (aRegionId == 1);

      // Degeneracies such as isolated vertices and sibling pairs can only be created by
      // intersecting CLOSED polygons or unioning OPEN polygons.
      boolean createDegen =
          (polygonModel == PolygonModel.CLOSED && !invertA && !invertB)
              || (polygonModel == PolygonModel.OPEN && invertA && invertB);

      // In addition, existing degeneracies are kept when an open boundary is subtracted. Note that
      // "keepDegenB" is only defined for completeness. It is needed to ensure that the "reverse
      // subtraction operator" (B - A) preserves degeneracies correctly, however in practice this
      // operator is only used internally to implement symmetric difference, and in that situation
      // the preserved degeneracy is always removed from the final result because it overlaps other
      // geometry.
      boolean keepDegenA = (polygonModel == PolygonModel.OPEN && invertB);
      boolean keepDegenB = (polygonModel == PolygonModel.OPEN && invertA);

      EdgeCrossingResult r = processEdgeCrossings(aId, a, it);
      assert !r.matchesPolyline;

      // If only one region is inverted, matching/sibling relations are reversed.
      if (invertA != invertB) {
        ShapeEdgeId tmp = r.polygonMatchId;
        r.polygonMatchId = r.siblingMatchId;
        r.siblingMatchId = tmp;
      }

      boolean isPoint = a.getStart().equalsPoint(a.getEnd());
      if (!emitShared) {
        // Remember the shell/hole status of degenerate B edges that are incident to any edge of A.
        // (We don't need to do this for vertex a1 since it is the same as vertex a0 of the
        // following A loop edge.)
        if (r.loopMatchesA0()) {
          isDegenerateHole.put(r.a0LoopMatchId, inside);
          if (isPoint) {
            return true;
          }
        }

        // Point loops are handled identically to points in the semi-open model, and are easier to
        // process in the first pass (since otherwise in the r.a0MatchesPolygon case we would need
        // to remember the containment status of the matching vertex). Otherwise we defer processing
        // such loops to the second pass so that we can distinguish whether the degenerate edge
        // represents a hole or shell.
        if (polygonModel != PolygonModel.SEMI_OPEN) {
          if (isPoint && r.a0MatchesPolygon) {
            return true;
          }
        }
      }

      inside ^= (r.a0Crossings & 1) == 1;
      if (!emitShared) {
        // Defer processing A edges that match or are siblings to an edge of B.
        if (r.matchesPolygon() || r.matchesSibling()) {
          // For sibling pairs, also remember their shell/hole status.
          if (r.matchesPolygon() && r.matchesSibling()) {
            isDegenerateHole.put(r.polygonMatchId, inside);
            isDegenerateHole.put(r.siblingMatchId, inside);
          }
          assert r.interiorCrossings == 0;
          inside ^= (r.a1Crossings & 1) == 1;
          return true;
        }
      }

      // Remember whether the B geometry represents a sibling pair hole.
      boolean isBHole = r.matchesPolygon() && r.matchesSibling() && inside;

      // WARNING: From here forward, "inside" changes meaning to "this edge should be emitted":
      //
      //   At this point, "inside" indicates whether the initial part of the A edge is contained by
      //   the B geometry using semi-open rules. The following code implements the various other
      //   polygon boundary rules by changing the value of "inside" when necessary to indicate
      //   whether the current A edge should be emitted to the output or not. "semiOpenInside"
      //   remembers the true value of "inside" so that it can be restored later.

      // TODO(torrey): Consider restructuring: instead of changing the meaning of inside and
      // introducing "semiOpenInside", introduce a "shouldEmit" boolean? Keep the C++ consistent.

      boolean semiOpenInside = inside;
      if (isPoint) {
        if (r.loopMatchesA0()) {
          // Both sides are point loops. The edge is kept only:
          //  - for closed intersection, open union, and open difference;
          //  - if A and B are both holes or both shells.
          inside = createDegen || keepDegenA || (inside == isDegenerateHole.get(r.a0LoopMatchId));
        } else if (r.a0MatchesPolygon) {
          // A point loop in A matches a polygon vertex in B. Note that this code can emit an extra
          // isolated vertex if A represents a point hole, but this doesn't matter (see comments on
          // the call to AddPointEdge below).
          if (polygonModel != PolygonModel.SEMI_OPEN) {
            inside = createDegen || keepDegenA;
          }
        }
      } else if (r.matchesPolygon()) {
        if (isDegenerate(aId)) {
          // The A edge has a sibling. The edge is kept only:
          //  - for closed intersection, open union, and open difference;
          //  - if the A sibling pair is a hole and the B edge has no sibling; or
          //  - if the B geometry is also a sibling pair and A and B are both holes or both shells.
          inside =
              createDegen
                  || keepDegenA
                  || (!r.matchesSibling() || inside) == isDegenerateHole.get(aId);
        } else {
          // Matching edges are kept unless the B geometry is a sibling pair, in which case it is
          // kept only for closed intersection, open union, and open difference.
          if (!r.matchesSibling() || createDegen || keepDegenB) {
            inside = true;
          }
        }
      } else if (r.matchesSibling()) {
        if (isDegenerate(aId)) {
          // The A edge has a sibling. The edge is kept only if A is a sibling pair shell and the
          // operation is closed intersection, open union, or open difference.
          inside = (createDegen || keepDegenA) && !isDegenerateHole.get(aId);
        } else {
          inside = createDegen;
        }
      }
      if (inside != semiOpenInside) {
        ++r.a1Crossings; // Restores the correct (semi-open) state later.
      }

      // Test whether the first vertex of this edge should be emitted as an isolated degenerate
      // vertex. This is only needed in the second pass when:
      //  - a0 matches a vertex of the B polygon;
      //  - the initial part of the A edge will not be emitted; and
      //  - the operation is closed intersection or open union, or open difference and the B
      //    geometry is a point loop.
      //
      // The logic does not attempt to avoid redundant extra vertices (e.g. the extra code in
      // processEdge1() that checks whether the vertex is the endpoint of the preceding emitted
      // edge) since these will be removed during S2BuilderGraph creation by
      // DegenerateEdges.DISCARD or DISCARD_EXCESS (which are necessary in any case due to
      // snapping).
      if (emitShared
          && r.a0MatchesPolygon
          && !inside
          && (createDegen || (keepDegenB && r.loopMatchesA0()))) {
        if (!addPointEdge(a.getStart(), 2)) {
          return false;
        }
      }

      // Since we skipped edges in the first pass that only had a sibling pair match in the B
      // geometry, we sometimes need to emit the sibling pair of an edge in the second pass. This
      // happens only if:
      //  - the operation is closed intersection, open union, or open difference;
      //  - the A geometry is not a sibling pair (since otherwise we will process that edge as
      //    well); and
      //  - the B geometry is not a sibling pair hole (since then only one edge should be emitted).
      if (r.matchesSibling() && (createDegen || keepDegenB) && !isDegenerate(aId) && !isBHole) {
        if (!addEdge(r.siblingMatchId, a.getEnd(), a.getStart(), /* dimension= */ 2, 0)) {
          return false;
        }
      }

      // Test whether the entire edge or any part of it belongs to the output.
      if (inside || r.interiorCrossings > 0) {
        // Note: updates "inside" to correspond to the state just before a1.
        if (!addEdge(aId, a.getStart(), a.getEnd(), /* dimension= */ 2, r.interiorCrossings)) {
          return false;
        }
      }
      inside ^= (r.a1Crossings & 1) == 1;

      // Verify that edge crossings are being counted correctly. This is a relatively expensive
      // check and should only be done in debug mode.
      assert !it.crossingsComplete()
          || (new S2ContainsPointQuery(it.bIndex()).contains(a.getEnd()) == inside ^ invertB);
      return true;
    }

    /**
     * Returns a summary of the relationship between a test edge from region A and a set of crossing
     * edges from region B (see EdgeCrossingResult). Advances the iterator over all those crossing
     * edges.
     *
     * <p>NOTE(ericv): We could save a bit of work when matching polygon vertices by passing in a
     * flag saying whether this information is needed. For example it is only needed in processEdge2
     * when (emitShared && createDegenerate).
     */
    private EdgeCrossingResult processEdgeCrossings(
        ShapeEdgeId aId, MutableEdge a, CrossingIterator it) {
      pendingSourceEdgeCrossings.clear();
      EdgeCrossingResult r = new EdgeCrossingResult();
      if (it.done(aId)) {
        return r;
      }

      MutableEdge b = new MutableEdge();
      for (; !it.done(aId); it.next()) {
        // Polyline and polygon "inside" states are not affected by point geometry.
        if (it.bDimension() == 0) {
          continue;
        }
        it.getBEdge(b);
        if (it.isInteriorCrossing()) {
          // The crossing occurs in the edge interior. The condition below says that (1) polyline
          // crossings don't affect the polygon "inside" state, and (2) subtracting a crossing
          // polyline from a polyline does not affect its "inside" state. (Note that vertices are
          // still created at the intersection points.)
          if (aDimension <= it.bDimension() && !(invertB != invertResult && it.bDimension() == 1)) {
            long srcId = SourceId.encode(bRegionId, it.bShapeId(), it.bEdgeId());
            addInteriorCrossing(new SourceEdgeCrossing(srcId, it.leftToRight()));
          }
          r.interiorCrossings += (it.bDimension() == 1) ? 2 : 1;
        } else if (it.bDimension() == 1) {
          // The polygon "inside" state is not affected by polyline geometry.
          if (aDimension == 2) {
            continue;
          }

          if (a.isEqualTo(b) || a.isSiblingOf(b)) {
            r.matchesPolyline = true;
          }
          if (b.hasEndpoint(a.getStart()) && polylineEdgeContainsVertex(a.getStart(), it, 1)) {
            r.a0MatchesPolyline = true;
          }
          if (b.hasEndpoint(a.getEnd()) && polylineEdgeContainsVertex(a.getEnd(), it, 1)) {
            r.a1MatchesPolyline = true;
          }
        } else {
          assert it.bDimension() == 2;
          if (edgeIsPoint(a) || edgeIsPoint(b)) {
            // There are no edge crossings since at least one edge is degenerate.
            if (edgeEqualsPoint(b, a.getStart())) {
              r.a0LoopMatchId = it.bId();
            }
          } else if (a.isEqualTo(b)) {
            ++r.a0Crossings;
            r.polygonMatchId = it.bId();
          } else if (a.isSiblingOf(b)) {
            ++r.a0Crossings;
            r.siblingMatchId = it.bId();
          } else if (it.isVertexCrossing()) {
            if (b.hasEndpoint(a.getStart())) {
              ++r.a0Crossings;
            } else {
              ++r.a1Crossings;
            }
          }
          if (b.hasEndpoint(a.getStart())) {
            r.a0MatchesPolygon = true;
          }
          if (b.hasEndpoint(a.getEnd())) {
            r.a1MatchesPolygon = true;
          }
        }
      }
      return r;
    }

    /**
     * Returns true if the vertex "v" is contained by the polyline edge referred to by the
     * CrossingIterator "it", taking into account the PolylineModel. "dimension" is 0 or 1 according
     * to whether "v" should be modeled as a point or as a degenerate polyline. (This only makes a
     * difference when the containing polyline is degenerate, since the polyline AA contains itself
     * in all boundary models but contains the point A only in the CLOSED model.)
     *
     * <p>REQUIRES: it.bDimension() == 1
     *
     * <p>REQUIRES: "v" is an endpoint of it.bEdge()
     */
    private boolean polylineEdgeContainsVertex(S2Point v, CrossingIterator it, int dimension) {
      assert it.bDimension() == 1;
      assert dimension == 0 || dimension == 1;

      // Closed polylines contain all their vertices.
      if (polylineModel == PolylineModel.CLOSED) {
        return true;
      }

      MutableEdge bEdge = new MutableEdge();
      it.getBEdge(bEdge);
      assert bEdge.hasEndpoint(v);

      int bEdgeId = it.bEdgeId();
      int bEdgeChainId = it.bEdgeChainId();
      int bChainStart = it.bShape().getChainStart(bEdgeChainId);
      int bChainLength = it.bShape().getChainLength(bEdgeChainId);
      int bChainLimit = bChainStart + bChainLength; // exclusive

      // A polyline contains its last vertex only when the polyline is degenerate (v0 == v1) and "v"
      // is modeled as a degenerate polyline (dimension == 1). This corresponds to the fact that the
      // polyline AA contains itself in all boundary models, but contains the point A only in the
      // CLOSED model.
      if (bEdgeId == bChainLimit - 1 // the last edge in the chain
          && v.equalsPoint(bEdge.getEnd()) // the last vertex in the chain
          && (dimension == 0 || bEdgeId > 0 || !v.equalsPoint(bEdge.getStart()))) {
        return false;
      }

      // Otherwise all interior vertices are contained. The first polyline vertex is contained if
      // either the polyline model is not OPEN, or the polyline forms a loop and
      // polylineLoopsHaveBoundaries is false.
      if (polylineContainsV0(bEdgeId, bChainStart)) {
        return true;
      }
      if (!v.equalsPoint(bEdge.getStart())) {
        return true;
      }
      if (polylineLoopsHaveBoundaries) {
        return false;
      }
      // Reuse bEdge to get the endpoint of the last chain edge to detect a loop.
      // TODO(torrey): Could this be more efficient with getChainVertex()?
      it.bShape().getChainEdge(bEdgeChainId, bChainLength - 1, bEdge);
      return v.equalsPoint(bEdge.getEnd());
    }

    /**
     * This method should be called after each pair of calls to startBoundary. (The only operation
     * that processes more than one pair of boundaries is SYMMETRIC_DIFFERENCE, which computes the
     * union of A-B and B-A.)
     *
     * <p>Resets the state of the CrossingProcessor. Translates the temporary representation of
     * crossing edges (SourceId) into the format expected by EdgeClippingLayer (InputEdgeId).
     */
    public void doneBoundaryPair() {
      // Add entries that translate the "special" crossings.
      sourceIdMap.put(SourceId.special(SET_INSIDE), SET_INSIDE);
      sourceIdMap.put(SourceId.special(SET_INVERT_B), SET_INVERT_B);
      sourceIdMap.put(SourceId.special(SET_REVERSE_A), SET_REVERSE_A);

      sourceEdgeCrossings.forEach(
          (firstInputEdgeId, secondSourceEdgeCrossing) -> {
            long sourceId = secondSourceEdgeCrossing.sourceId();

            assert sourceIdMap.containsKey(sourceId);
            int secondInputEdgeId = sourceIdMap.get(sourceId);

            inputCrossings.add(
                new CrossingInputEdgePair(
                    firstInputEdgeId,
                    new CrossingInputEdge(
                        secondInputEdgeId, secondSourceEdgeCrossing.leftToRight())));
          });

      sourceEdgeCrossings.clear();
      sourceIdMap.clear();
    }

    /**
     * Indicates whether the point being processed along the current edge chain is in the polygonal
     * interior of the opposite region, using semi-open boundaries. If "invertB" is true then this
     * field is inverted.
     *
     * <p>This value along with the set of incident edges can be used to compute whether the
     * opposite region contains this point under any of the supported boundary models
     * (PolylineModel.CLOSED, etc).
     */
    public boolean inside() {
      return inside;
    }

    /**
     * SourceEdgeCrossing represents edge crossings and changes in clipping state. Some
     * SourceEdgeCrossings use special SourceIds for SET_REVERSE_A, SET_INVERT_B, or SET_INSIDE,
     * while the common case is to represent an input that crosses some other edge, and the
     * direction of that crossing. If leftToRight is true, the input edge crosses the other edge
     * from left to right.
     */
    private static class SourceEdgeCrossing {

      /** The SourceId encoded into a long. */
      private final long srcId;

      /** True if the edge represented by srcId crosses the other edge from left to right. */
      private final boolean leftToRight;

      public SourceEdgeCrossing(long srcId, boolean leftToRight) {
        this.srcId = srcId;
        this.leftToRight = leftToRight;
      }

      public long sourceId() {
        return srcId;
      }

      public boolean leftToRight() {
        return leftToRight;
      }

      @Override
      public String toString() {
        return "SourceEdgeCrossing(" + srcId + ", " + leftToRight + ")";
      }
    }

    private static interface SourceEdgeCrossingVisitor {
      public void visit(int i, SourceEdgeCrossing crossing);
    }

    /**
     * A temporary representation of inputCrossings that is used internally until all necessary
     * edges from *both* polygons have been emitted to the S2Builder. This field is then converted
     * by doneBoundaryPair() into the InputEdgeCrossings format expected by GraphEdgeClipper.
     *
     * <p>The reason that we can't construct inputCrossings directly is that it uses InputEdgeIds to
     * identify the edges from both polygons, and when we are processing edges from the first
     * polygon, InputEdgeIds have not yet been assigned to the second polygon. So instead this field
     * identifies edges from the first polygon using an InputEdgeId (int), and edges from the second
     * polygon using a (regionId, shapeId, edgeId) tuple (i.e., a SourceId).
     *
     * <p>All crossings are represented twice, once to indicate that an edge from polygon 0 is
     * crossed by an edge from polygon 1, and once to indicate that an edge from polygon 1 is
     * crossed by an edge from polygon 0. The entries are sorted lexicographically by their eventual
     * InputEdgeIds except for GraphEdgeClipper state modifications, which are sorted by the first
     * InputEdgeId only.
     */
    private static class SourceEdgeCrossings {
      // TODO(torrey): This is a placeholder implementation. Consider a more memory efficient
      // way to store this data.

      /**
       * A crossing edge. The "int inputEdgeId" is the input edge id of the edge in the first
       * polygon, and the SourceEdgeCrossing is the SourceId of the edge in the second polygon,
       * along with the crossing direction.
       */
      private static class InputSourceEdgeCrossing {
        public final int inputEdgeId;
        public final SourceEdgeCrossing sourceEdgeCrossing;

        public InputSourceEdgeCrossing(int inputEdgeId, SourceEdgeCrossing crossing) {
          this.inputEdgeId = inputEdgeId;
          this.sourceEdgeCrossing = crossing;
        }
      }

      /** The list of crossing edges. */
      final ArrayList<InputSourceEdgeCrossing> list = new ArrayList<>();

      /**
       * Adds a pair of crossing edges where the first edge is the input edge id 'i' from the first
       * polygon, and the second edge and crossing direction are in the given SourceEdgeCrossing.
       */
      public void add(int inputEdgeId, SourceEdgeCrossing crossing) {
        list.add(new InputSourceEdgeCrossing(inputEdgeId, crossing));
      }

      /** Visits each of the input crossings. */
      public void forEach(SourceEdgeCrossingVisitor visitor) {
        for (InputSourceEdgeCrossing crossing : list) {
          visitor.visit(crossing.inputEdgeId, crossing.sourceEdgeCrossing);
        }
      }

      /** Clears the list of crossings. */
      public void clear() {
        list.clear();
      }
    }

    private int inputEdgeId() {
      // The addEdge method keeps essentially parallel arrays of values including the input
      // dimension. So when crossings of that edge are added, the edge ID that crossings are being
      // added to is just the most recent ID that addEdge was working with, which is the number of
      // addEdge calls so far, and therefore also the size of any of the parallel arrays.
      return inputDimensions.size();
    }

    /**
     * Returns true if the edges on either side of the first vertex of the current edge have not
     * been emitted.
     *
     * <p>REQUIRES: This method is called just after updating "inside" for "v0".
     */
    private boolean isV0Isolated(ShapeEdgeId aId) {
      return !inside && v0EmittedMaxEdgeId < aId.edgeId;
    }

    /**
     * Returns true if "aId" is the last edge of the current chain, and the edges on either side of
     * the last vertex have not been emitted (including the possibility that the chain forms a
     * loop).
     */
    private boolean isChainLastVertexIsolated(ShapeEdgeId aId) {
      return (aId.edgeId == chainLimit - 1 && !chainV0Emitted && v0EmittedMaxEdgeId <= aId.edgeId);
    }

    /**
     * Returns true if the given polyline edge contains "v0", taking into account the specified
     * PolylineModel.
     */
    private boolean polylineContainsV0(int edgeId, int chainStart) {
      return (polylineModel != PolylineModel.OPEN || edgeId > chainStart);
    }

    private boolean isDegenerate(ShapeEdgeId aId) {
      return isDegenerateHole.containsKey(aId);
    }

    private void addCrossing(SourceEdgeCrossing crossing) {
      sourceEdgeCrossings.add(inputEdgeId(), crossing);
    }

    private void addInteriorCrossing(SourceEdgeCrossing crossing) {
      // Crossing edges are queued until the S2Builder edge that they are supposed to be associated
      // with is created (see addEdge() and pendingsourceEdgeCrossings for details).
      pendingSourceEdgeCrossings.add(crossing);
    }

    private void setClippingState(int parameter, boolean state) {
      addCrossing(new SourceEdgeCrossing(SourceId.special(parameter), state));
    }

    /**
     * Supports "early exit" in the case of boolean results by returning false as soon as the result
     * is known to be non-empty.
     */
    private boolean addEdge(
        ShapeEdgeId aId, S2Point start, S2Point end, int dimension, int interiorCrossings) {
      if (builder == null) {
        return false; // Boolean output.
      }
      if (interiorCrossings > 0) {
        // Add the edges that cross this edge to the output so that GraphEdgeClipper can find them.
        for (SourceEdgeCrossing crossing : pendingSourceEdgeCrossings) {
          sourceEdgeCrossings.add(inputEdgeId(), crossing);
        }
        // Build a map that translates temporary edge ids (SourceId) to the representation used by
        // EdgeClippingLayer (InputEdgeId).
        long srcId = SourceId.encode(aRegionId, aId.shapeId, aId.edgeId);
        sourceIdMap.put(srcId, inputEdgeId());
      }

      // Set the GraphEdgeClipper's "inside" state to match ours.
      if (inside != prevInside) {
        setClippingState(SET_INSIDE, inside);
      }
      inputDimensions.add(dimension);
      builder.addEdge(start, end);
      inside ^= (interiorCrossings & 1) == 1;
      prevInside = inside;
      return true;
    }

    /**
     * Supports "early exit" in the case of boolean results by returning false as soon as the result
     * is known to be non-empty.
     */
    private boolean addPointEdge(S2Point p, int dimension) {
      if (builder == null) {
        return false; // Boolean output.
      }
      if (!prevInside) {
        setClippingState(SET_INSIDE, true);
      }

      inputDimensions.add(dimension);
      builder.addEdge(p, p);
      prevInside = true;
      return true;
    }
  }

  /** Returns true iff the endpoints of the given edge are exactly equal. */
  public static boolean edgeIsPoint(MutableEdge e) {
    return e.getStart().equalsPoint(e.getEnd());
  }

  /** Returns true iff both endpoints of the given edge are equal to the given point. */
  public static boolean edgeEqualsPoint(MutableEdge e, S2Point point) {
    return point.equalsPoint(e.getStart()) && point.equalsPoint(e.getEnd());
  }
}
