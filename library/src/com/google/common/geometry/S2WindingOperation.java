/*
 * Copyright 2025 Google Inc.
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

import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2Builder.IsFullPolygonPredicate;
import com.google.common.geometry.S2Builder.SnapFunction;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.primitives.IdSetLexicon;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.util.List;

/**
 * Given a set of possibly self-intersecting closed loops, this class computes a partitioning of the
 * sphere into regions of constant winding number and returns the subset of those regions selected
 * by a given winding rule. This functionality can be used to implement N-way boolean polygon
 * operations including union, intersection, and symmetric difference, as well as more exotic
 * operations such as buffering and Minkowski sums.
 *
 * <p>Recall that in the plane, the winding number of a closed curve around a point is an integer
 * representing the number of times the curve travels counterclockwise around that point. For
 * example, points in the interior of a planar simple polygon with counterclockwise boundary have a
 * winding number of +1, while points outside the polygon have a winding number of 0. If the
 * boundary is clockwise then points in the interior have a winding number of -1.
 *
 * <p>The interior of a complex self-intersecting closed boundary curve may be defined by choosing
 * an appropriate winding rule. For example, the interior could be defined as all regions whose
 * winding number is positive, or all regions whose winding number is non-zero, or all regions whose
 * winding number is odd. Different winding rules are useful for implementing the various operations
 * mentioned above (union, symmetric difference, etc).
 *
 * <p>Unfortunately the concept of winding number on the sphere is not well-defined, due to the fact
 * that a given closed curve does not partition the sphere into regions of constant winding number.
 * This is because the winding number changes when a point crosses either the given curve or that
 * curve's reflection through the origin.
 *
 * <p>Instead we engage in a slight abuse of terminology by modifying the concept of winding number
 * on the sphere to be a relative one. Given a set of closed curves on the sphere and the winding
 * number of a reference point R, we define the winding number of every other point P by counting
 * signed edge crossings. When the edge RP crosses from the right side of a curve to its left the
 * winding number increases by one, and when the edge RP crosses from the left side of a curve to
 * its right the winding number decreases by one.
 *
 * <p>This definition agrees with the classical one in the plane when R is taken to be the point at
 * infinity with a winding number of zero. It also agrees with the rule used by the S2 library to
 * define polygon interiors, namely that the interior of a loop is the region to its left. And most
 * importantly, it satisfies the property that a closed curve partitions the sphere into regions of
 * constant winding number.
 */
@SuppressWarnings("Assertion")
public class S2WindingOperation {

  /** Options for an S2WindingOperation. */
  public static class Options {
    private SnapFunction snapFunction;
    private boolean includeDegeneracies = false;

    public Options() {
      this.snapFunction = new IdentitySnapFunction(S1Angle.ZERO);
    }

    // Convenience constructor that calls setSnapFunction().
    public Options(SnapFunction snapFunction) {
      this.snapFunction = snapFunction;
    }

    // Copy constructor
    public Options(Options options) {
      this.snapFunction = options.snapFunction;
      this.includeDegeneracies = options.includeDegeneracies;
    }

    /**
     * Sets the SnapFunction to be used for snap rounding the output during the call to build(). The
     * default is an IdentitySnapFunction(S1Angle.ZERO).
     */
    public void setSnapFunction(SnapFunction snapFunction) {
      this.snapFunction = snapFunction;
    }

    /** Returns the SnapFunction used for snap rounding the output during the call to build(). */
    public SnapFunction snapFunction() {
      return snapFunction;
    }

    /**
     * This option can be enabled to provide limited support for degeneracies (i.e., sibling edge
     * pairs and isolated vertices). By default the output does not include such edges because they
     * do not bound any interior. If this option is true, then the output includes additional
     * degeneracies as follows:
     *
     * <ul>
     *   <li>{@code WindingRule.ODD} outputs degeneracies whose multiplicity is odd.
     *   <li>All other winding rules output degeneracies contained by regions whose winding number
     *       is zero.
     * </ul>
     *
     * These rules are sufficient to implement the following useful operations:
     *
     * <ul>
     *   <li>{@code WindingRule.ODD} can be used to compute the N-way symmetric difference of any
     *       combination of points, polylines, and polygons.
     *   <li>{@code WindingRule.POSITIVE} can be used to compute the N-way union of any combination
     *       of points, polylines, and polygons.
     * </ul>
     *
     * These statements only apply when closed boundaries are being used (see {@code
     * S2BooleanOperation.{Polygon,Polyline}Model.CLOSED}) and require the client to convert points
     * and polylines to degenerate loops and then back again (e.g. using {@code
     * s2builderutil.ClosedSetNormalizer}). Note that the handling of degeneracies is NOT sufficient
     * for other polygon operations (e.g. intersection) or other boundary models (e.g, semi-open).
     *
     * <p>DEFAULT: false
     */
    public void setIncludeDegeneracies(boolean includeDegeneracies) {
      this.includeDegeneracies = includeDegeneracies;
    }

    /**
     * Returns whether degeneracies should be included in the output. See {@link
     * #setIncludeDegeneracies}.
     */
    public boolean includeDegeneracies() {
      return includeDegeneracies;
    }
  }

  private Options options;
  private S2Builder builder;
  private int refInputEdgeId;
  private int refWindingIn;
  private WindingRule rule;

  /** Constructor that does not set the result layer or options, requires init() to be called. */
  public S2WindingOperation() {}

  /** Convenience constructor that calls init(). */
  public S2WindingOperation(S2BuilderLayer resultLayer, Options options) {
    init(resultLayer, options);
  }

  // Starts an operation that sends its output to the given S2Builder layer.
  // This method may be called more than once.
  public void init(S2BuilderLayer resultLayer, Options options) {
    this.options = options;
    S2Builder.Builder builderOptions = new S2Builder.Builder(options.snapFunction());
    builderOptions.setSplitCrossingEdges(true);

    builder = builderOptions.build();
    WindingLayer windingLayer = new WindingLayer(this, resultLayer);
    builder.startLayer(windingLayer);
  }

  public void init(S2BuilderLayer resultLayer) {
    init(resultLayer, new Options());
  }

  public Options options() {
    return options;
  }

  /**
   * Adds a loop to the set of loops used to partition the sphere. The given loop may be degenerate
   * or self-intersecting.
   */
  public void addLoop(List<S2Point> loop) {
    builder.addLoop(loop);
  }

  /**
   * Specifies the winding rule used to determine which regions belong to the result.
   *
   * <p>Note that additional winding rules may be created by adjusting the winding number of the
   * reference point. For example, to select regions whose winding number is at least 10, simply
   * subtract 9 from refWinding and then use WindingRule.POSITIVE. (This can be used to implement
   * N-way polygon intersection). Similarly, additional behaviors can be obtained by reversing some
   * of the input loops (e.g., this can be used to compute one polygon minus the union of several
   * other polygons).
   */
  public enum WindingRule {
    POSITIVE, // Winding number > 0
    NEGATIVE, // Winding number < 0
    NON_ZERO, // Winding number != 0
    ODD // Winding number % 2 == 1
  }

  /**
   * Given a reference point {@code refP} whose winding number is given to be {@code refWinding},
   * snaps all the given input loops together using the given {@code snapFunction()} and then
   * partitions the sphere into regions of constant winding number. As discussed above, the winding
   * number increases by one when any loop edge is crossed from right to left, and decreases by one
   * when any loop edge is crossed from left to right. Each partition region is classified as
   * belonging to the interior of the result if and only if its winding number matches the given
   * rule (e.g. {@code WindingRule.POSITIVE}).
   *
   * <p>The output consists of the boundary of this interior region plus possibly certain
   * degeneracies (as controlled by the {@code includeDegeneracies()} option). The boundary edges
   * are sent to the S2Builder result layer specified in the constructor, along with an appropriate
   * IsFullPolygonPredicate that can be used to distinguish whether the result is empty or full
   * (even when degeneracies are present). Note that distinguishing empty from full results is a
   * problem unique to spherical geometry.
   *
   * @param refP The reference point.
   * @param refWinding The winding number of the reference point.
   * @param rule The winding rule to use.
   * @param error If an error occurs, error will be initialized with the details.
   * @return true if the build was successful, otherwise initializes the error and returns false.
   */
  public boolean build(S2Point refP, int refWinding, WindingRule rule, S2Error error) {
    // The reference point must be an S2Builder input vertex in order to determine how its winding
    // number is affected by snapping. We record the input edge id of this edge so that we can find
    // it later.
    this.refInputEdgeId = builder.numInputEdges();
    builder.addPoint(refP);
    this.refWindingIn = refWinding;
    this.rule = rule;
    return builder.build(error);
  }

  /** The actual winding number operation is implemented as an S2BuilderLayer. */
  private static class WindingLayer implements S2BuilderLayer {
    private final S2WindingOperation op;
    private final S2BuilderLayer resultLayer;

    // The graph data that will be sent to resultLayer.
    private final EdgeList resultEdges = new EdgeList();
    private final IntArrayList resultInputEdgeIds = new IntArrayList();

    public WindingLayer(S2WindingOperation op, S2BuilderLayer resultLayer) {
      this.op = op;
      this.resultLayer = resultLayer;
    }

    @Override
    public GraphOptions graphOptions() {
      // The algorithm below has several steps with different graph requirements, however the first
      // step is to determine how snapping has affected the winding number of the reference point.
      // This requires keeping all degenerate edges and sibling pairs. We also keep all duplicate
      // edges since this makes it easier to remove the reference edge.
      return new GraphOptions(
          EdgeType.DIRECTED, DegenerateEdges.KEEP, DuplicateEdges.KEEP, SiblingPairs.KEEP);
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      if (!error.ok()) {
        return false; // Abort if an error already exists.
      }

      // The WindingOracle computes the winding number of any point on the sphere. It requires
      // knowing the winding number at a reference point which must be an input vertex to S2Builder.
      // This is achieved by adding a degenerate edge from the reference point to itself
      // (refInputEdgeId). Once we have computed the change in winding number, we create a new graph
      // with this edge removed (since it should not be emitted to the result layer).
      WindingOracle oracle = new WindingOracle(op.refInputEdgeId, op.refWindingIn, op.builder, g);
      assert error.ok() : error.text(); // No errors are possible.

      // Now that we have computed the change in winding number, we create a new graph with the
      // reference edge removed.
      EdgeList newEdges = new EdgeList();
      IntArrayList newInputEdgeIds = new IntArrayList();
      IdSetLexicon newInputEdgeIdSetLexicon = g.inputEdgeIdSetLexicon();

      // Copy all of the edges except the reference edge.
      for (int e = 0; e < g.numEdges(); ++e) {
        if (g.inputEdgeIds(e).first() == op.refInputEdgeId) {
          continue;
        }
        newEdges.copyEdge(g.edges(), e);
        newInputEdgeIds.add(g.inputEdgeIdSetId(e));
      }

      // Our goal is to assemble the given edges into loops that partition the sphere. In order to
      // do this we merge duplicate edges and create sibling edges so that every region can have its
      // own directed boundary loop.
      //
      // Isolated degenerate edges and sibling pairs are preserved in order to provide limited
      // support for working with geometry of dimensions 0 and 1 (i.e., points and polylines).
      // Clients can simply convert these objects to degenerate loops and then convert these
      // degenerate loops back to points/polylines in the output using S2ClosedSetNormalizer.
      GraphOptions newOptions =
          new GraphOptions(
              EdgeType.DIRECTED,
              DegenerateEdges.DISCARD_EXCESS,
              DuplicateEdges.MERGE,
              SiblingPairs.CREATE);
      S2BuilderGraph newGraph =
          g.makeSubgraph(
              newOptions, newEdges, newInputEdgeIds, newInputEdgeIdSetLexicon, null, error);
      if (!error.ok()) {
        return false;
      }

      // Now visit each connected component of the graph and assemble its edges into loops. For each
      // loop we determine the winding number of the region to its left, and if the winding number
      // matches the given rule then we add its edges to resultEdges.
      computeBoundary(newGraph, oracle, error);

      // Now we construct yet another S2BuilderGraph by starting with the edges that bound the
      // desired output region and then processing them according to the client's requested
      // GraphOptions. (Note that processEdges() can change these options in certain cases; see
      // S2Builder.GraphOptions for details.)
      //
      // The IsFullPolygonPredicate below allows clients to distinguish full from empty results
      // (including cases where there are degeneracies). Note that we can use the winding number at
      // any point on the sphere for this purpose.
      IsFullPolygonPredicate isFullPolygonPredicate = gr -> matchesRule(oracle.currentRefWinding());

      IdSetLexicon resultIdSetLexicon = newGraph.inputEdgeIdSetLexicon();
      S2BuilderGraph resultGraph =
          newGraph.makeSubgraph(
              resultLayer.graphOptions(),
              resultEdges,
              resultInputEdgeIds,
              resultIdSetLexicon,
              isFullPolygonPredicate,
              error);
      boolean unused = resultLayer.build(resultGraph, error);
      return error.ok();
    }

    private void computeBoundary(S2BuilderGraph g, WindingOracle oracle, S2Error error) {
      // We assemble the edges into loops using an algorithm similar to
      // S2BuilderGraph.getDirectedComponents(), except that we also keep track of winding numbers
      // and output the relevant edges as we go along.
      //
      // The following accounts for siblingMap, leftTurnMap, and edgeWinding, which have
      // g.numEdges() elements each.
      IntArrayList siblingMap = g.getSiblingMap();
      int[] leftTurnMap = new int[g.numEdges()];
      g.fillLeftTurnMap(siblingMap, leftTurnMap, error);
      if (!error.ok()) {
        return;
      }

      // A map from EdgeId to the winding number of the region it bounds.
      int[] edgeWinding = new int[g.numEdges()];
      IntArrayList frontier = new IntArrayList(); // Unexplored sibling edge ids.
      for (int edgeIdMin = 0; edgeIdMin < g.numEdges(); ++edgeIdMin) {
        if (leftTurnMap[edgeIdMin] < 0) {
          continue; // Already visited.
        }

        // We have found a new connected component of the graph. Each component consists of a set of
        // loops that partition the sphere. We start by choosing an arbitrary vertex "v0" and
        // computing the winding number at this vertex. Recall that point containment is defined
        // such that when a set of loops partition the sphere, every point is contained by exactly
        // one loop. Therefore the winding number at "v0" is the same as the winding number of the
        // adjacent loop that contains it. We choose "e0" to be an arbitrary edge of this loop (it
        // is the incoming edge to "v0").
        int vertexId0 = g.edgeDstId(edgeIdMin);

        int edge0Id = getContainingLoopEdgeId(vertexId0, edgeIdMin, g, leftTurnMap, siblingMap);
        edgeWinding[edge0Id] = oracle.getWindingNumber(g.vertex(vertexId0));

        // Now visit all loop edges in this connected component of the graph.
        // "frontier" is a stack of unexplored siblings of the edges visited far.
        frontier.add(edge0Id);
        while (!frontier.isEmpty()) {
          int edgeId = frontier.popInt();
          if (leftTurnMap[edgeId] < 0) {
            continue; // Already visited.
          }

          // Visit all edges of the loop starting at "edgeId".
          int winding = edgeWinding[edgeId];
          for (int nextEdgeId; leftTurnMap[edgeId] >= 0; edgeId = nextEdgeId) {
            // Count signed edge crossings to determine the winding number of the sibling region.
            // Input edges that snapped to "e" decrease the winding number by one (since we cross
            // them from left to right), while input edges that snapped to its sibling edge increase
            // the winding number by one (since we cross them from right to left).
            int siblingEdgeId = siblingMap.getInt(edgeId);
            int windingMinus = g.inputEdgeIds(edgeId).size();
            int windingPlus = g.inputEdgeIds(siblingEdgeId).size();
            int siblingWinding = winding - windingMinus + windingPlus;

            // Output all edges that bound the region selected by the winding rule, plus certain
            // degenerate edges.
            if ((matchesRule(winding) && !matchesRule(siblingWinding))
                || matchesDegeneracy(winding, windingMinus, windingPlus)) {
              outputEdge(g, edgeId);
            }
            nextEdgeId = leftTurnMap[edgeId];
            leftTurnMap[edgeId] = -1;
            // If the sibling hasn't been visited yet, add it to the frontier.
            if (leftTurnMap[siblingEdgeId] >= 0) {
              edgeWinding[siblingEdgeId] = siblingWinding;
              frontier.add(siblingEdgeId);
            }
          }
        }
      }
    }

    /**
     * Given an incoming "startEdgeId" to "vertexId", returns an edge of the loop that contains
     * "vertexId" with respect to the usual semi-open boundary rules.
     */
    private int getContainingLoopEdgeId(
        int vertexId,
        int startEdgeId,
        S2BuilderGraph g,
        int[] leftTurnMap,
        IntArrayList siblingMap) {
      // If the given edge is degenerate, this is an isolated vertex.
      assert g.edgeDstId(startEdgeId) == vertexId;
      if (g.edgeSrcId(startEdgeId) == vertexId) {
        return startEdgeId;
      }

      // Otherwise search for the loop that contains "v".
      int edgeId0 = startEdgeId;
      for (; ; ) {
        int edgeId1 = leftTurnMap[edgeId0];
        assert g.edgeDstId(edgeId0) == vertexId;
        assert g.edgeSrcId(edgeId1) == vertexId;

        // The first test below handles the case where there are only two edges incident to this
        // vertex (i.e., the vertex angle is 360 degrees).
        if (g.edgeSrcId(edgeId0) == g.edgeDstId(edgeId1)
            || S2Predicates.angleContainsVertex(
                g.vertex(g.edgeSrcId(edgeId0)),
                g.vertex(vertexId),
                g.vertex(g.edgeDstId(edgeId1)))) {
          return edgeId0;
        }
        edgeId0 = siblingMap.getInt(edgeId1);
        assert edgeId0 != startEdgeId;
      }
    }

    private boolean matchesRule(int winding) {
      switch (op.rule) {
        case POSITIVE:
          return winding > 0;
        case NEGATIVE:
          return winding < 0;
        case NON_ZERO:
          return winding != 0;
        case ODD:
          return (winding & 1) != 0;
      }
      throw new IllegalStateException("Unknown WindingRule " + op.rule);
    }

    private boolean matchesDegeneracy(int winding, int windingMinus, int windingPlus) {
      if (!op.options.includeDegeneracies()) {
        return false;
      }

      // A degeneracy is either a self-loop or a sibling pair where equal numbers of input edges
      // snapped to both edges of the pair. The test below covers both cases (because self-loops are
      // their own sibling).
      if (windingMinus != windingPlus) {
        return false;
      }

      if (op.rule == WindingRule.ODD) {
        // Any degeneracy whose multiplicity is odd should be part of the result independent of the
        // winding number of the region that contains it. This rule allows computing symmetric
        // differences of any combination of points, polylines, and polygons (where the first two
        // are represented as degenerate loops).
        return (windingPlus & 1) != 0;
      } else {
        // For all other winding rules we output degeneracies only if they are contained by a region
        // of winding number zero. Even though the interface to this class does not provide enough
        // information to allow consistent handling of degeneracies in general, this rule is
        // sufficient for several important cases. Specifically it allows computing N-way unions of
        // any mixture of points, polylines, and polygons by converting the points and polylines to
        // degenerate loops. In this case all input loops are degenerate or CCW, and the boundary of
        // the result can be computed using WindingRule.POSITIVE. Since there are no clockwise
        // loops, all degeneracies contained by a region of winding number zero represent degenerate
        // shells and should be emitted. (They can be converted back to points/polylines
        // using S2ClosedSetNormalizer.)
        //
        // Similarly, this heuristic is sufficient to compute unions of points, polylines, and
        // polygons where all boundaries are clockwise (by using WindingRule.NEGATIVE) or where all
        // boundaries are of an unknown but consistent orientation (by using WindingRule.NON_ZERO).
        return winding == 0;
      }
    }

    private void outputEdge(S2BuilderGraph g, int edgeId) {
      resultEdges.copyEdge(g.edges(), edgeId);
      resultInputEdgeIds.add(g.inputEdgeIdSetId(edgeId));
    }
  }

  /**
   * The purpose of WindingOracle is to compute winding numbers with respect to a set of input loops
   * after snapping. It is given the input edges (via S2Builder), the output eges (an
   * S2Builder.Graph), and the winding number at a reference point with respect to the input edges
   * ("refInputEdgeId" and "refWindingIn"). getWindingNumber() can then be called to compute the
   * winding number at any point with respect to the snapped loops. The main algorithm uses this to
   * compute the winding number at one arbitrary vertex of each connected component of the S2Builder
   * output graph.
   */
  private static class WindingOracle {
    private final S2BuilderGraph g;

    // The current reference point. Initially this is simply the snapped position of the input
    // reference point, but it is updated on each call to getWindingNumber (see below).
    private S2Point refP;

    // The winding number at the current reference point.
    private int refWinding;

    // For each connected component of the S2Builder output graph, we determine the winding number
    // at one arbitrary vertex by counting edge crossings. Initially we do this by brute force, but
    // if there are too many connected components then we build an S2ShapeIndex to speed up the
    // process.
    //
    // Building an index takes about as long as 25 brute force queries. The classic competitive
    // algorithm technique would be to do 25 brute force queries and then build the index. However
    // in practice, inputs with one connected component are common while inputs with 2-25 connected
    // components are rare. Therefore we build the index as soon as we realize that there is more
    // than one component.
    //
    // Another idea is to count the connected components in advance, however it turns out that this
    // takes about 25% as long as building the index does.
    private int bruteForceWindingTestsLeft = 1;
    private final S2ShapeIndex index; // Built only if needed.

    public WindingOracle(
        int refInputEdgeId, int refWindingIn, S2Builder builder, S2BuilderGraph g) {
      this.g = g;
      assert g.options().edgeType() == EdgeType.DIRECTED;
      assert g.options().degenerateEdges() == DegenerateEdges.KEEP;
      assert g.options().duplicateEdges() == DuplicateEdges.KEEP;
      assert g.options().siblingPairs() == SiblingPairs.KEEP;

      // Compute the winding number at the reference point after snapping (see
      // GetSnappedWindingDelta).
      S2Point refIn = builder.inputEdgeSrcVertex(refInputEdgeId);
      int refVertexId = GetSnappedWindingDelta.findFirstVertexId(refInputEdgeId, g);
      assert refVertexId >= 0; // No errors are possible.
      refP = g.vertex(refVertexId);
      S2Error error = new S2Error();
      refWinding =
          refWindingIn
              + GetSnappedWindingDelta.getSnappedWindingDelta(
                  refIn,
                  refVertexId,
                  GetSnappedWindingDelta.InputEdgeFilter.noFilter(),
                  builder,
                  g,
                  error);
      assert error.ok() : error.text(); // No errors are possible.

      // Winding numbers at other points are computed by counting signed edge crossings. If we need
      // to do this many times, it is worthwhile to build an index. Note that although we initialize
      // the index here, it is only built the first time we use it (i.e., when
      // bruteForceWindingTestsLeft < 0).
      //
      // As a small optimization, we set maxEdgesPerCell() higher than its default value. This makes
      // it faster to build the index but increases the time per query. This is a good tradeoff
      // because the number of calls to getWindingNumber() is typically small relative to the number
      // of graph edges. It also saves memory, which is important when the graph is very large (e.g.
      // because the input loops have many self-intersections).
      final int kMaxEdgesPerCell = 40;
      S2ShapeIndex.Options options = new S2ShapeIndex.Options();
      options.setMaxEdgesPerCell(kMaxEdgesPerCell);
      this.index = new S2ShapeIndex(options);
      index.add(new S2BuilderUtil.GraphShape(g));
    }

    /** Returns the winding number at the given point after snapping. */
    public int getWindingNumber(S2Point p) {
      // Count signed edge crossings starting from the reference point, whose winding number is
      // known. If we need to do this many times then we build an S2ShapeIndex to speed up this
      // process.
      EdgeCrosser crosser = new EdgeCrosser(refP, p);
      int winding = refWinding;
      if (--bruteForceWindingTestsLeft >= 0) {
        for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
          winding += signedCrossingDelta(crosser, edgeId);
        }
      } else {
        S2CrossingEdgeQuery query = new S2CrossingEdgeQuery(index);
        S2CrossingEdgeQuery.Edges edges = query.getCandidates(refP, p, 0);
        while (!edges.isEmpty()) {
          winding += signedCrossingDelta(crosser, edges.nextEdge());
        }
      }
      // It turns out that getWindingNumber() is called with a sequence of points that are sorted in
      // approximate S2CellId order. This means that if we update our reference point as we go
      // along, the edges for which we need to count crossings are much shorter on average.
      refP = p;
      refWinding = winding;
      return winding;
    }

    /** Returns the winding number at the current reference point. */
    public int currentRefWinding() {
      return refWinding;
    }

    /** Returns the change in winding number due to crossing the given graph edge. */
    private int signedCrossingDelta(EdgeCrosser crosser, int edgeId) {
      return crosser.signedEdgeOrVertexCrossing(
          g.vertex(g.edgeSrcId(edgeId)), g.vertex(g.edgeDstId(edgeId)));
    }
  }
}
