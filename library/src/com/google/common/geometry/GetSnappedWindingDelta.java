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

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.MultimapBuilder;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IdSetLexicon.IdSet;
import com.google.common.geometry.primitives.Ints.OfInt;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.function.IntPredicate;

/**
 * Provides the static {@link #getSnappedWindingDelta} and {@link #findFirstVertexId} methods and
 * supporting code used by S2WindingOperation.
 */
@SuppressWarnings("Assertion")
final class GetSnappedWindingDelta {

  /**
   * An input edge may snap to zero, one, or two non-degenerate output edges incident to the
   * reference vertex, consisting of at most one incoming and one outgoing edge.
   */
  static class EdgeSnap {
    /** The input edge. */
    public MutableEdge inputEdge = new MutableEdge();

    /**
     * If vertexIdIn >= 0, the input edge snapped to an edge from vertexIdIn to the reference
     * vertex. Otherwise, if vertexIdIn == -1, then the snapped source vertex is the reference
     * vertex.
     */
    public int vertexIdIn = -1;

    /**
     * If vertexIdOut >= 0, the input edge snapped to an edge from the reference vertex to
     * vertexIdOut. Otherwise, if vertexIdOut == -1, then the snapped destination vertex is the
     * reference vertex.
     */
    public int vertexIdOut = -1;

    @Override
    public String toString() {
      return "EdgeSnap{"
          + "inputEdge="
          + inputEdge
          + ", vertexIdIn="
          + vertexIdIn
          + ", vertexIdOut="
          + vertexIdOut
          + '}';
    }
  }

  /**
   * A function that returns true if the given S2Builder input edge id should be ignored in the
   * winding number calculation. This means either that the edge is not a loop edge (e.g., a
   * non-closed polyline) or that this loop should not affect the winding number. This is useful for
   * two purposes:
   *
   * <ul>
   *   <li>To process graphs that contain polylines and points in addition to loops.
   *   <li>To process graphs where the winding number is computed with respect to only a subset of
   *       the input loops.
   * </ul>
   *
   * <p>The factory function taking no inputs ignores no edges.
   */
  static interface InputEdgeFilter extends IntPredicate {
    public static InputEdgeFilter noFilter() {
      return edgeId -> false;
    }
  }

  /**
   * Builds a chain of input edges (and a chain of corresponding snapped output edges), where the
   * input edges are either a loop where all vertices snapped to the reference vertex, or a path
   * from a vertex outside the Voronoi region of the reference vertex, to one or more vertices
   * snapped to the reference vertex, and then to a vertex outside the Voronoi region.
   *
   * <p>Fills in chainIn and chainOut from the given graph and inputVertexEdgeMap, which is a
   * multimap from input edge starting vertex to EdgeSnap, with an entry for every input edge that
   * snapped to a graph edge incident to the reference vertex.
   *
   * <p>Returns true if the chain was built successfully, or false if an error occurred, in which
   * case the given S2Error is initialized.
   */
  private static boolean buildChain(
      int refVertexId,
      S2BuilderGraph g,
      Multimap<S2Point, EdgeSnap> inputVertexEdgeMap,
      ListPlus<S2Point> chainIn,  // Output parameter: The chain of input edges
      ListPlus<S2Point> chainOut, // Output parameter: The chain of output edges
      S2Error error) {
    assert chainIn.isEmpty();
    assert chainOut.isEmpty();

    // First look for an incoming graph edge to the reference vertex. (This will be the start of a
    // chain that eventually leads to an outgoing edge.)
    {
      // Iterate over EdgeSnaps (snapped graph edges that are incident to the reference vertex)
      // keyed by the starting input edge vertex.
      Iterator<Entry<S2Point, EdgeSnap>> iter = inputVertexEdgeMap.entries().iterator();
      Entry<S2Point, EdgeSnap> entry = null;
      EdgeSnap snap = null;
      while (iter.hasNext()) {
        entry = iter.next();
        snap = entry.getValue();
        if (snap.vertexIdIn >= 0) {
          // The snapped start vertex isn't the reference vertex, so this is an incoming graph edge
          // to the reference vertex.
          chainOut.add(g.vertex(snap.vertexIdIn));
          break;
        } // Otherwise, the snapped start vertex is the reference vertex.
        // The snapped end vertex might also be the reference vertex, or it may be outgoing.
      }

      // If we didn't find an incoming graph edge, we're at an arbitrary incident edge where we can
      // start a closed loop. Either way, build the chains from this snapped edge.
      iter.remove();
      chainIn.add(snap.inputEdge.getStart());
      chainIn.add(snap.inputEdge.getEnd());
      chainOut.add(g.vertex(refVertexId));

      if (snap.vertexIdOut >= 0) {
        // This input edge exits the Voronoi region, we are done.
        chainOut.add(g.vertex(snap.vertexIdOut));
        return true;
      }
    }

    // The end of the last input edge snapped to the reference vertex.
    // Now repeatedly add edges until the chain or loop is finished.
    while (!chainIn.back().equalsPoint(chainIn.front())) {
      // We have not yet formed a loop. Get outgoing snapped edges from the end of the chain of
      // input edges.
      S2Point end = chainIn.back();
      Collection<EdgeSnap> range = inputVertexEdgeMap.get(end);
      if (range.isEmpty()) {
        // There are no snapped graph edges for input edges leaving the end of the chain.
        error.init(
            S2Error.Code.INVALID_ARGUMENT, "Input edges (after filtering) do not form loops");
        return false;
      }
      // Extend the chain with an arbitrary outgoing edge, both snapped and input.
      EdgeSnap snap = range.iterator().next();
      inputVertexEdgeMap.remove(end, snap);
      chainIn.add(snap.inputEdge.getEnd());

      if (snap.vertexIdOut >= 0) {
        // The chain has exited the Voronoi region.
        chainOut.add(g.vertex(snap.vertexIdOut));
        break;
      }
    }
    return true;
  }

  /**
   * Returns the change in winding number along the edge AB with respect to the given edge chain.
   * This is simply the sum of the signed edge crossings.
   */
  private static int getEdgeWindingDelta(S2Point a, S2Point b, List<S2Point> chain) {
    assert !chain.isEmpty();

    int delta = 0;
    EdgeCrosser crosser = new EdgeCrosser(a, b, chain.get(0));
    for (int i = 1; i < chain.size(); ++i) {
      delta += crosser.signedEdgeOrVertexCrossing(chain.get(i));
    }
    return delta;
  }

  /**
   * Given an input edge (B0, B1) that snaps to an edge chain (B0', B1', ...), returns a connecting
   * vertex "Bc" that can be used as a substitute for the remaining snapped vertices "..." when
   * computing winding numbers. This requires that (1) the edge (B1', Bc) does not intersect the
   * Voronoi region of B0', and (2) the edge chain (B0', B1', Bc, B1) always stays within the snap
   * radius of the input edge (B0, B1).
   */
  private static S2Point getConnector(S2Point b0, S2Point b1, S2Point b1Snapped) {
    // If B1' within 90 degrees of B1, no connecting vertex is necessary.
    if (b1Snapped.dotProd(b1) >= 0) {
      return b1;
    }

    // Otherwise we use the point on (B0, B1) that is 90 degrees away from B1'.
    // This is sufficient to ensure conditions (1) and (2).
    S2Point x = S2RobustCrossProd.robustCrossProd(b0, b1).crossProd(b1Snapped).normalize();
    return (x.dotProd(S2EdgeUtil.interpolate(b0, b1, 0.5)) >= 0) ? x : x.neg();
  }

  /**
   * Returns the set of incoming and outgoing edges incident to the given vertex id. This method
   * takes time linear in the size of the graph "g"; if you need to call this function many times
   * then it is more efficient to use Graph.VertexOutMap and Graph.VertexInMap instead.
   */
  private static IntArrayList getIncidentEdgesBruteForce(int vertexId, S2BuilderGraph g) {
    IntArrayList result = new IntArrayList();
    for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
      if (g.edgeSrcId(edgeId) == vertexId || g.edgeDstId(edgeId) == vertexId) {
        result.add(edgeId);
      }
    }
    return result;
  }

  /**
   * See the documentation of getSnappedWindingDelta() below. This version can be used when
   * getSnappedWindingDelta() needs to be called many times on the same graph. It is faster than the
   * function below, but less convenient to use because it requires the client to provide the set of
   * graph edges incident to the snapped reference vertex. It runs in time proportional to the size
   * of this set.
   *
   * <p>"incidentEdges" is the set of incoming and outgoing graph edges incident to refVertexId.
   * (These edges can be found efficiently using S2BuilderGraph.VertexOutMap and
   * S2BuilderGraph.VertexInMap.)
   *
   * <p>See the function below for the remaining parameters and requirements.
   */
  public static int getSnappedWindingDelta(
      S2Point refIn,
      int refVertexId,
      IntArrayList incidentEdges,
      InputEdgeFilter inputEdgeFilter,
      S2Builder builder,
      S2BuilderGraph g,
      S2Error error) {
    assert !builder.simplifyEdgeChains();
    assert g.options().edgeType() == S2Builder.EdgeType.DIRECTED;
    assert g.options().degenerateEdges() == GraphOptions.DegenerateEdges.KEEP;
    assert g.options().siblingPairs() == GraphOptions.SiblingPairs.KEEP
        || g.options().siblingPairs() == GraphOptions.SiblingPairs.REQUIRE
        || g.options().siblingPairs() == GraphOptions.SiblingPairs.CREATE;

    // First we group all the incident edges by input edge id, to handle the problem that input
    // edges can map to either one or two snapped edges. Specifically, we construct a map from input
    // edge ids to EdgeSnaps, for every input edge that snapped to every graph edge incident to the
    // reference vertex.
    Int2ObjectOpenHashMap<EdgeSnap> inputIdEdgeMap = new Int2ObjectOpenHashMap<>();
    for (int edgeId : incidentEdges) {
      // The vertex ids of the graph edge endpoints.
      int edgeSrcId = g.edgeSrcId(edgeId);
      int edgeDstId = g.edgeDstId(edgeId);

      // Iterate over the input edges snapped to this incident edge.
      IdSet inputEdgeIds = g.inputEdgeIds(edgeId);
      inputEdgeIds.forEach(
          inputEdgeId -> {
            if (inputEdgeFilter != null && inputEdgeFilter.test(inputEdgeId)) {
              return;
            }
            // If this input edge maps to two snapped edges, we may already have an EdgeSnap for it.
            // If there's an existing EdgeSnap, overwrite it. Otherwise create one.
            // TODO(torrey): Explain why it is ok to just overwrite the existing EdgeSnap.
            EdgeSnap snap = inputIdEdgeMap.computeIfAbsent(inputEdgeId, k -> new EdgeSnap());
            builder.getInputEdge(inputEdgeId, snap.inputEdge);
            if (edgeSrcId != refVertexId) {
              // Input edge snapped to an incoming edge to the reference vertex.
              snap.vertexIdIn = edgeSrcId;
            }
            if (edgeDstId != refVertexId) {
              // Input edge snapped to an outgoing edge from the reference vertex.
              snap.vertexIdOut = edgeDstId;
            }
          });
    }

    // Now we regroup the edges according to the starting vertex of the corresponding input edge.
    // This makes it easier to assemble these edges into (portions of) input edge loops.
    // A Multimap from input edge starting vertex to EdgeSnap is an inefficient data structure, but
    // we expect the number of edges here to be small, normally.
    ListMultimap<S2Point, EdgeSnap> inputVertexEdgeMap =
        MultimapBuilder.linkedHashKeys().arrayListValues().build();
    inputIdEdgeMap
        .values()
        .forEach(snap -> inputVertexEdgeMap.put(snap.inputEdge.getStart(), snap));

    // The position of the reference vertex after snapping. In comments we will refer to the
    // reference vertex before and after snapping as R and R'.
    S2Point refOut = g.vertex(refVertexId);

    // Now we repeatedly assemble edges into an edge chain and determine the change in winding
    // number due to snapping of that edge chain. These values are summed to compute the final
    // winding number delta.
    //
    // An edge chain is either a closed loop of input vertices where every vertex snapped to the
    // reference vertex R', or a partial loop such that all interior vertices snap to R' while the
    // first and last vertex do not. Note that the latter includes the case where the first and last
    // input vertices in the chain are identical but do not snap to R'.
    //
    // Essentially we compute the winding number of the unsnapped reference vertex R with respect to
    // the input chain and the winding number of the snapped reference vertex R' with respect to the
    // output chain, and then subtract them. In the case of open chains, we compute winding numbers
    // as if the chains had been closed in a way that preserves topology while snapping (i.e., in a
    // way that does not the cause chain to pass through the reference vertex as it continuously
    // deforms from the input to the output).
    //
    // Any changes to this code should be validated by running the RandomLoops unit test with at
    // least 10 million iterations.
    int windingDelta = 0;
    while (!inputVertexEdgeMap.isEmpty()) {
      // An edge chain of input vertices
      ListPlus<S2Point> chainIn = new ListPlus<>();
      ListPlus<S2Point> chainOut = new ListPlus<>();

      if (!buildChain(refVertexId, g, inputVertexEdgeMap, chainIn, chainOut, error)) {
        return Integer.MAX_VALUE;
      }

      if (chainOut.size() == 1) {
        // We have a closed chain C of input vertices such that every vertex snaps to the reference
        // vertex R'. Therefore we can easily find a point Z whose winding number is not affected by
        // the snapping of C; it just needs to be outside the Voronoi region of R'. Since the snap
        // radius is at most 70 degrees, we can use a point 90 degrees away such as S2.ortho(R').
        //
        // We then compute the winding numbers of R and R' relative to Z. We compute the winding
        // number of R by counting signed crossings of the edge ZR, while the winding number of R'
        // relative to Z is always zero because the snapped chain collapses to a single point.
        assert chainOut.get(0).equalsPoint(refOut); // Snaps to R'.
        assert chainIn.get(0).equalsPoint(chainIn.back()); // Chain is a loop.
        S2Point z = S2.ortho(refOut);
        windingDelta += 0 - getEdgeWindingDelta(z, refIn, chainIn);
      } else {
        // We have an input chain C = (A0, A1, ..., B0, B1) that snaps to a chain
        // C' = (A0', R', B1'), where A0 and B1 are outside the Voronoi region of R' and all other
        // input vertices snap to R'. This includes the case where A0 == B1 and also the case where
        // the input chain consists of only two vertices. Note that technically the chain C snaps to
        // a supersequence of C', since A0A1 snaps to a chain whose last two vertices are (A0', R')
        // and B0B1 snaps to a chain whose first two vertices are (R', B1'). This implies that A0
        // does not necessarily snap to A0', and similarly for B1 and B1'.
        //
        // Note that A0 and B1 can be arbitrarily far away from R'. This makes it difficult (on the
        // sphere) to construct a point Z whose winding number is guaranteed not to be affected by
        // snapping the edges A0A1 and B0B1. Instead we construct two points Za and Zb such that Za
        // is guaranteed not be affected by the snapping of A0A1, Zb is guaranteed not to be
        // affected by the snapping of B0B1, and both points are guaranteed not to be affected by
        // the snapping of any other input edges. We can then compute the change in winding number
        // of Zb by considering only the single edge A0A1 that snaps to A0'R'. Furthermore we can
        // compute the latter by using Za as the reference point, since its winding number is
        // guaranteed not be affected by this particular edge.
        //
        // Given the point Zb, whose change in winding number is now known, we can compute the
        // change in winding number of the reference vertex R. We essentially want to count the
        // signed crossings of ZbR' with respect to C' and subtract the signed crossings of ZbR with
        // respect to C, which we will write as s(ZbR', C') - s(ZbR, C).
        //
        // However to do this we need to close both chains in a way that is topologically consistent
        // and does not affect the winding number of Zb. This can be achieved by counting the signed
        // crossings of ZbR' by following the two-edge path (Zb, R, R'). In other words, we compute
        // this as s(ZbR, C') + s(RR', C') - s(ZbR, C). We can then compute s(ZbR, C') - s(ZbR, C)
        // by simply concatenating the vertices of C in reverse order to C' to form a single closed
        // loop. The remaining term s(RR', C') can be implemented as signed edge crossing tests, or
        // more directly by testing whether R is contained by the wedge C'.
        assert chainOut.size() == 3;
        assert chainOut.get(1).equalsPoint(refOut);

        // Compute two points Za and Zb such that Za is not affected by the snapping of any edge
        // except possibly B0B1, and Zb is not affected by the snapping of any edge except possibly
        // A0A1. Za and Zb are simply the normals to the edges A0A1 and B0B1 respectively, with
        // their sign chosen to point away from the Voronoi site R'. This ensures at least 20
        // degrees of separation from all edges except the ones mentioned.
        S2Point za = S2RobustCrossProd.robustCrossProd(chainIn.get(0), chainIn.get(1)).normalize();
        S2Point zb = S2RobustCrossProd.robustCrossProd(chainIn.get(-2), chainIn.back()).normalize();
        if (za.dotProd(refOut) > 0) {
          za = za.neg();
        }
        if (zb.dotProd(refOut) > 0) {
          zb = zb.neg();
        }

        // We now want to determine the change in winding number of Zb due to A0A1 snapping to
        // A0'R'. Conceptually we do this by closing these two single-edge chains into loops L and
        // L' and then computing s(ZaZb, L') - s(ZaZb, L). Recall that Za was constructed so as not
        // to be affected by the snapping of edge A0A1, however this is only true provided that L
        // can snap to L' without passing through Za.
        //
        // To achieve this we let L be the degenerate loop (A0, A1, A0), and L' be the loop
        // (A0', R', A1, A0, A0'). The only problem is that we need to ensure that the edge A0A0'
        // stays within 90 degrees of A0A1, since otherwise when the latter edge snaps to the former
        // it might pass through Za. (This problem arises because we only consider the last two
        // vertices (A0', R') that A0A1 snaps to. If we used the full chain of snapped vertices for
        // A0A1 then L' would always stay within the snap radius of this edge.)
        //
        // The simplest way to fix this problem is to insert a connecting vertex Ac between A0 and
        // A0'. THis vertex acts as a proxy for the missing snapped vertices, yielding the loop L'
        // = (A0', R', A1, A0, Ac, A0'). The vertex Ac is located on the edge A0A1 and is at most 90
        // degrees away from A0'. This ensures that the chain (A0, Ac, A0') always stays within the
        // snap radius of the input edge A0A1.
        //
        // Similarly we insert a connecting vertex Bc between B0 and B1 to ensure that the edge
        // B1'B1 passes on the correct side of Zb.
        S2Point a0Connector = getConnector(chainIn.get(1), chainIn.get(0), chainOut.get(0));
        S2Point b1Connector = getConnector(chainIn.get(-2), chainIn.back(), chainOut.get(2));

        // Compute the change in winding number for reference vertex Zb. Note that we must duplicate
        // the first/last vertex of the loop since the argument to getEdgeWindingDelta() is a
        // polyline.
        ImmutableList<S2Point> chainZ =
            ImmutableList.of(
                chainOut.get(0),
                chainOut.get(1),
                chainIn.get(1),
                chainIn.get(0),
                a0Connector,
                chainOut.get(0));
        windingDelta += getEdgeWindingDelta(za, zb, chainZ);

        // Compute the change in winding number of ZbR due to snapping C to C'. As before,
        // conceptually we do this by closing these chains into loops L and L' such that L snaps to
        // L' without passing through Zb. Again this can be done by concatenating the vertices of C'
        // with the reversed vertices of C, along with the two extra connecting vertices Ac and Bc
        // to ensure that L and L' pass on the same side of Za and Zb. This yields the loop
        // (A0', R', B1', Bc, B1, B0, ..., A1, A0, Ac, A0').
        List<S2Point> chainDiff = new ArrayList<>(chainOut);
        chainDiff.add(b1Connector);
        chainDiff.addAll(chainIn.reverse());
        chainDiff.add(a0Connector);
        chainDiff.add(chainOut.get(0)); // Close the loop.
        windingDelta += getEdgeWindingDelta(zb, refIn, chainDiff);

        // Compute the change in winding number of RR' with respect to C' only.
        // (This could also be computed using two calls to s2pred.OrderedCCW.)
        assert chainOut.get(1).equalsPoint(refOut);
        windingDelta += getEdgeWindingDelta(refIn, refOut, chainOut);
      }
    }
    return windingDelta;
  }

  /**
   * Given an S2BuilderGraph of output edges after snap rounding and a reference vertex R, computes
   * the change in winding number of R due to snapping. (See S2WindingOperation for an introduction
   * to winding numbers on the sphere.) The return value can be added to the original winding number
   * of R to obtain the winding number of the corresponding snapped vertex R'.
   *
   * <p>The algorithm requires that the S2Builder input edges consist entirely of (possibly self-
   * intersecting) closed loops. If you need to process inputs that include other types of geometry
   * (e.g., non-closed polylines), you will need to either (1) put them into a different S2Builder
   * layer, (2) close the polylines into loops (e.g. using GraphOptions.SiblingEdges.CREATE), or (3)
   * provide a suitable InputEdgeFilter (see above) so that the non-loop edges can be ignored.
   *
   * <p>The algorithm is designed to be robust for any input edge configuration and snapping result.
   * However note that it cannot be used in conjunction with edge chain simplification
   * (S2Builder.Options.simplifyEdgeChains). It also requires that S2Builder.GraphOptions be
   * configured to keep all snapped edges, even degenerate ones (see requirements below).
   *
   * <p>"refIn" is the reference vertex location before snapping. It *must* be an input vertex to
   * S2Builder, however this is not checked.
   *
   * <p>"refVertexId" is the vertex id of the reference vertex after snapping. (This can be found
   * using the findFirstVertexId() function below if desired.)
   *
   * <p>"inputEdgeFilter" can optionally be used to ignore input edges that should not affect the
   * winding number calculation (such as polyline edges). The value can {@link
   * InputEdgeFilter#noFilter()} to use all edges.
   *
   * <p>"builder" is the S2Builder that produced the given edge graph. It is used to map input edge
   * ids back to the original edge definitions, and also to verify that no incompatible
   * S2Builder.Options were used (see below).
   *
   * <p>"g" is the S2Builder output graph of snapped edges.
   *
   * <p>The only possible errors are usage errors, in which case "error" is initialized with an
   * appropriate error message and {@code Integer.MAX_VALUE} is returned.
   *
   * <p>Example usage:
   *
   * <p>This function is generally called from an S2Builder.Layer implementation. We assume here
   * that the reference vertex is the first vertex of the input edge identified by "refInputEdgeId",
   * and that its desired winding number with respect to the input loops is "refWinding".
   *
   * {@snippet :
   *   class SomeLayer : public S2Builder.Layer {
   *     private int refInputEdgeId;
   *     private int refWinding;
   *     private S2Builder builder;
   *
   *     public void build(S2BuilderGraph g, S2Error error) {
   *       // Find the positions of the reference vertex before and after snapping.
   *       S2Point refIn = builder.inputEdgeSrcVertex(refInputEdgeId);
   *       int refVertexId = findFirstVertexId(refInputEdgeId, g);
   *       S2Point refOut = g.vertex(refVertexId);
   *
   *       // Compute the change in winding number due to snapping.
   *       S2Error error = new S2Error();
   *       refWinding += getSnappedWindingDelta(
   *           refIn,
   *           refVertexId,
   *           InputEdgeFilter.noFilter(),
   *           builder,
   *           g,
   *           error);
   *       assert error.ok();  // All errors are usage errors.
   *
   *       // Winding numbers of others points can now be found by counting signed edge crossings
   *       // (S2EdgeCrosser.signedEdgeOrVertexCrossing) between "refOut" and the desired point.
   *       // Note that if DuplicateEdges.MERGE or SiblingPairs.CREATE was used, each crossing has
   *       // a multiplicity equal to the number of non-filtered input edges that snapped to that
   *       // output edge.
   *     }
   *   }
   * }
   *
   * <p>REQUIRES: The input edges after filtering consist entirely of closed loops. (If
   * DuplicateEdges.MERGE or SiblingPairs.CREATE was used, each graph edge has a multiplicity equal
   * to the number of non-filtered input edges that snapped to it.)
   *
   * <p>REQUIRES: g.options().edgeType() == DIRECTED
   *
   * <p>REQUIRES: g.options().degenerateEdges() == KEEP
   *
   * <p>REQUIRES: g.options().siblingPairs() == {KEEP, REQUIRE, CREATE}
   *
   * <p>REQUIRES: builder.options().simplifyEdgeChains() == false
   *
   * <p>CAVEAT: The running time is proportional to the S2BuilderGraph size. If you need to call
   * this function many times on the same graph then use the alternate version below. (Most clients
   * only need to call getSnappedWindingDelta() once per graph because the winding numbers of other
   * points can be computed by counting signed edge crossings.)
   */
  public static int getSnappedWindingDelta(
      S2Point refIn,
      int refVertexId,
      InputEdgeFilter inputEdgeFilter,
      S2Builder builder,
      S2BuilderGraph g,
      S2Error error) {
    return getSnappedWindingDelta(
        refIn,
        refVertexId,
        getIncidentEdgesBruteForce(refVertexId, g),
        inputEdgeFilter,
        builder,
        g,
        error);
  }

  /**
   * Returns the first vertex id of the snapped edge chain for the given input edge, or -1 if this
   * input edge does not exist in the graph "g".
   */
  public static int findFirstVertexId(int inputEdgeId, S2BuilderGraph g) {
    // A given input edge snaps to a chain of output edges. To determine which output vertex the
    // source of the given input edge snaps to, we must find the first edge in this chain.
    //
    // The search below takes linear time in the number of edges; it can be done more efficiently if
    // duplicate edges are not being merged and the mapping returned by Graph.GetInputEdgeOrder() is
    // available. The algorithm would also be much simpler if inputEdgeId were known to be
    // degenerate.

    Int2IntOpenHashMap excessDegreeMap = new Int2IntOpenHashMap();

    for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
      IdSetLexicon.IdSet idSet = g.inputEdgeIds(edgeId);
      for (OfInt iter = idSet.intIterator(); iter.hasNext(); ) {
        int id = iter.nextInt();
        if (id == inputEdgeId) {
          excessDegreeMap.compute(g.edgeSrcId(edgeId), (k, v) -> v == null ? 1 : v + 1);
          excessDegreeMap.compute(g.edgeDstId(edgeId), (k, v) -> v == null ? -1 : v - 1);
          break;
        }
      }
    }
    if (excessDegreeMap.isEmpty()) {
      return -1; // Does not exist.
    }

    // Look for the (unique) vertex whose excess degree is +1.
    ObjectIterator<Int2IntMap.Entry> iter = excessDegreeMap.int2IntEntrySet().fastIterator();
    while (iter.hasNext()) {
      Int2IntMap.Entry entry = iter.next();
      if (entry.getIntValue() == 1) {
        return entry.getIntKey();
      }
    }

    // Otherwise "inputEdgeId" must snap to a single degenerate edge.
    assert excessDegreeMap.size() == 1;
    iter = excessDegreeMap.int2IntEntrySet().fastIterator();
    return iter.next().getIntKey();
  }

  /**
   * Extends ArrayList with some minor conveniences to simplify porting from C++ vectors and spans.
   */
  private static class ListPlus<A> extends ArrayList<A> {

    /**
     * Returns a new list containing the contents of this list in reverse order. Note that it would
     * generally be better to use {@link java.util.List#reversed()}, but that is not yet available
     * in J2CL, or in Android SDK 21.
     */
    public List<A> reverse() {
      ArrayList<A> r = new ArrayList<>(this); // Copy this list.
      Collections.reverse(r); // Reverse in place.
      return r;
    }

    /** Returns the last entry in the list. */
    public A back() {
      return super.get(size() - 1);
    }

    /** Returns the first entry in the list. */
    public A front() {
      return super.get(0);
    }

    /**
     * If the index is greater than or equal to zero, sets the value at the index. If the index is
     * less than zero, sets the value counting backward from the end of the list, i.e. set(-1, X) is
     * the same as set(size() - 1, X).
     */
    @Override
    public A set(int index, A val) {
      if (index < 0) {
        return super.set(size() + index, val);
      }
      return super.set(index, val);
    }

    /**
     * If the index is greater than or equal to zero, gets the value at the index. If the index is
     * less than zero, gets the value counting backward from the end of the list, i.e. get(-1) is
     * the same as get(size() - 1).
     */
    @Override
    public A get(int index) {
      if (index < 0) {
        return super.get(size() + index);
      }
      return super.get(index);
    }
  }

  /** Not instantiable. */
  private GetSnappedWindingDelta() {}
}
