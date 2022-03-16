package com.google.common.geometry;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import java.util.List;

public class S2PointVectorShape implements S2Shape {

    private final List<S2Point> points;

    public S2PointVectorShape(List<S2Point> points) {
        this.points = points;
    }

    public S2PointVectorShape() {
        this(Lists.newArrayList());
    }

    public int numPoints() {
        return points.size();
    }

    /**
     * Returns the number of edges in this shape.
     */
    @Override
    public int numEdges() {
        return numPoints();
    }

    /**
     * Returns the edge for the given index in {@code result}. Must not return zero-length edges.
     *
     * @param index  which edge to set into {@code result}, from 0 to {@link #numEdges()} - 1
     * @param result
     */
    @Override
    public void getEdge(int index, MutableEdge result) {
        result.a = points.get(index);
        result.b = points.get(index);
    }

    /**
     * Returns true if this shape has an interior, i.e. the shape consists of one or more closed
     * non-intersecting loops.
     */
    @Override
    public boolean hasInterior() {
        return false;
    }

    /**
     * Returns true if this shape contains {@link S2#origin()}. Should return false for shapes that do
     * not have an interior.
     */
    @Override
    public boolean containsOrigin() {
        return false;
    }

    /**
     * Returns the number of contiguous edge chains in the shape. For example, a shape whose edges are
     * [AB, BC, CD, AE, EF] may consist of two chains [A, B, C, D] and [A, E, F]. Every chain is
     * assigned a chain id numbered sequentially starting from zero.
     *
     * <p>An empty shape has no chains. A full shape (which contains the entire globe) has one chain
     * with no edges. Other shapes should have at least one chain, and the sum of all valid {@link
     * #getChainLength(int) chain lengths} should equal {@link #numEdges()} (that is, edges may only
     * be used by a single chain).
     *
     * <p>Note that it is always acceptable to implement this method by returning {@link #numEdges()}
     * (i.e. every chain consists of a single edge), but this may reduce the efficiency of some
     * algorithms.
     */
    @Override
    public int numChains() {
        return numPoints();
    }

    /**
     * Returns the first edge id corresponding to the edge chain for the given chain id. The edge
     * chains must form contiguous, non-overlapping ranges that cover the entire range of edge ids.
     *
     * @param chainId which edge chain to return its start, from 0 to {@link #numChains()} - 1
     */
    @Override
    public int getChainStart(int chainId) {
        return chainId;
    }

    /**
     * Returns the number of edge ids corresponding to the edge chain for the given chain id. The edge
     * chains must form contiguous, non-overlapping ranges that cover the entire range of edge ids.
     *
     * @param chainId which edge chain to return its length, from 0 to {@link #numChains()} - 1
     */
    @Override
    public int getChainLength(int chainId) {
        return 1;
    }

    /**
     * Returns the edge for the given chain id and offset in {@code result}. Must not return
     * zero-length edges.
     *
     * @param chainId which chain contains the edge to return, from 0 to {@link #numChains()} - 1
     * @param offset  position from chain start for the edge to return, from 0 to {@link
     *                #getChainLength(int)} - 1
     * @param result
     */
    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge result) {
        Preconditions.checkArgument(offset == 0);
        result.a = points.get(chainId);
        result.b = result.a;
    }

    /**
     * Returns the start point of the edge that would be returned by {@link S2Shape#getChainEdge},
     * or the endpoint of the last edge if {@code edgeOffset==getChainLength(chainId)}.
     *
     * @param chainId
     * @param edgeOffset
     */
    @Override
    public S2Point getChainVertex(int chainId, int edgeOffset) {
        Preconditions.checkArgument(edgeOffset == 0);
        return points.get(chainId);
    }

    /**
     * Returns the dimension of the geometry represented by this shape.
     *
     * <ul>
     *   <li>0 - Point geometry. Each point is represented as a degenerate edge.
     *   <li>1 - Polyline geometry. Polyline edges may be degenerate. A shape may represent any number
     *       of polylines. Polylines edges may intersect.
     *   <li>2 - Polygon geometry. Edges should be oriented such that the polygon interior is always
     *       on the left. In theory the edges may be returned in any order, but typically the edges
     *       are organized as a collection of edge chains where each chain represents one polygon
     *       loop. Polygons may have degeneracies, e.g., degenerate edges or sibling pairs consisting
     *       of an edge and its corresponding reversed edge. A polygon loop may also be full
     *       (containing all points on the sphere); by convention this is represented as a chain with
     *       no edges.
     * </ul>
     *
     * <p>Note that this method allows degenerate geometry of different dimensions to be
     * distinguished, e.g., it allows a point to be distinguished from a polyline or polygon that has
     * been simplified to a single point.
     */
    @Override
    public int dimension() {
        return 0;
    }
}
