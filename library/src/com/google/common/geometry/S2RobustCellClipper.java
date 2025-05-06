/*
 * Copyright 2024 Google Inc.
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

import static com.google.common.geometry.R2EdgeClipper.BOTTOM;
import static com.google.common.geometry.R2EdgeClipper.INSIDE;
import static com.google.common.geometry.R2EdgeClipper.LEFT;
import static com.google.common.geometry.R2EdgeClipper.OUTSIDE;
import static com.google.common.geometry.R2EdgeClipper.RIGHT;
import static com.google.common.geometry.R2EdgeClipper.TOP;
import static com.google.common.geometry.S2Cell.Boundary.BOTTOM_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.LEFT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.RIGHT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.TOP_EDGE;
import static com.google.common.geometry.S2RobustCrossProd.robustCrossProd;
import static java.lang.Math.abs;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.geometry.S2Cell.Boundary;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2IndexCellData.EdgeAndIdChain;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.primitives.PooledList;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import org.jspecify.annotations.Nullable;

/**
 * This class clips edges of a single S2Shape to cell boundaries robustly, and provides information
 * on which vertices are contained in the cell, as well as where the cell boundary is crossed and
 * contained by the shape.
 *
 * <p>By robust, we don't mean that the intersection point between an edge and the boundary of a
 * cell is computed exactly. Rather, we determine consistently and robustly which boundaries of a
 * cell an edge crosses, and order the crossings along the boundary properly.
 *
 * <p>Importantly, no matter how close a crossing is to the boundary of a cell, we will never
 * "accidentally" move the crossing to the adjacent boundary.
 *
 * <p>To summarize, crossing intercept values may have error, but they're always the same exact
 * value in adjacent cells, and which boundary and the order they cross are always correct!
 *
 * <p>Specific guarantees:
 *
 * <ol>
 *   <li>Whether or not an edge crosses the boundary of a cell is determined using exact predicates
 *       when necessary. We use a consistent perturbation scheme so that any given edge either
 *       crosses or does not cross a given cell boundary.
 *   <li>If an edge does cross a cell boundary, then we compute an intersection point on the
 *       boundary for it.
 *   <li>The intercept point computed for a crossing will be exactly the same as the intercept
 *       computed in any cell(s) on the opposite side of the boundary. A particular edge crossing,
 *       no matter how close to a cell corner, will always cross consistently in either the U axis
 *       or the V axis.
 *       <p>The intersection point may have error on it, but that error is promised to be exactly
 *       the same when clipping in any cell on the other side of the boundary too.
 *       <p>A consequence of requiring that the intercept points compute to exactly the same value
 *       is that they may fall slightly outside of the cell proper due to rounding error, and that
 *       needs to be handled in user logic.
 *   <li>Boundary crossings are ordered correctly s.t. one can walk around the boundary of the cell,
 *       set an interior bit at each outgoing crossing and unset it at each incoming crossing to
 *       determine a consistent set of cell boundary segments that are contained by a polygon.
 *       <p>If crossings are closer than some maximum error, then we fall back to using an exact
 *       predicate to compute this ordering. It will always be correct.
 * </ol>
 *
 * <p>Example Usage
 *
 * <p>The clipper needs all the crossing points available so that it can resolve the ordering of
 * crossings if needed; So all the edges should be added first, and then the crossings can be used:
 *
 * {@snippet :
 * S2RobustCellClipper clipper = new S2RobustCellClipper();
 * clipper.startCell(cell);
 *
 * for (S2Shape.MutableEdge edge : edges) {
 *   if (clipper.clipEdge(edge)) {
 *     // edge intersected the cell
 *   }
 * }
 *
 * for (Crossing crossing : clipper.getCrossings()) {
 *   // <use crossings>
 * }
 * }
 *
 * <p>Note: We currently expect the clipped edges to belong to a single valid shape so that
 * assumptions about e.g. exactly equal boundary crossings hold.
 */
@SuppressWarnings("Assertion")
public final class S2RobustCellClipper {
  /** Maximum UV error we can have from conversion and clipping. */
  public static final double MAX_ERROR =
      S2Projections.MAX_XYZ_TO_UV_ERROR
          + S2EdgeUtil.EDGE_CLIP_ERROR_UV_COORD
          + S2EdgeUtil.FACE_CLIP_ERROR_UV_COORD;

  /** The options for this robust cell clipper. */
  private final Options options;

  /** The cell currently being clipped to. */
  private S2Cell cell = null;

  /** Cached boundary components of the cell. */
  private final MutableEdge[] boundaries = {
    new MutableEdge(), new MutableEdge(), new MutableEdge(), new MutableEdge()
  };

  /**
   * UV coordinates of the cell in the order [y.lo, x.hi, y.hi, x.lo]. Also known as [vMin, uMax,
   * vMax, uMin]
   */
  private final double[] uvcoords = new double[4];

  /** Normals of the cell edges. */
  private final S2Point[] normals = new S2Point[4];

  /** The center of the current cell. */
  private S2Point cellCenter;

  /** A point that is outside the current cell. */
  private S2Point outside;

  /** The clipper used to clip edges to the cell. */
  private final UVEdgeClipper clipper = new UVEdgeClipper();

  /**
   * Crossings that were computed. These are sorted and deduplicated before being returned. They are
   * stored in a pooled list for reuse for clipping edges in another cell.
   */
  private final PooledList<Crossing> crossings = new PooledList<>(Crossing::new);

  /** A flag indicating whether or not we've sorted "crossings" yet. */
  private boolean needSorting = false;

  /**
   * The edges that crossed the boundary of the cell. Crossings reference this list, and we need the
   * actual edge information to sort the crossings on demand when they're accessed.
   */
  private final List<MutableEdge> crossingEdges = new ArrayList<>();

  /**
   * The edges that were strictly contained by the cell. We need these edges to determine whether
   * the cell boundary is contained or not.
   */
  private final List<MutableEdge> containedEdges = new ArrayList<>();

  /** Constructs an S2RobustCellClipper with default options. */
  public S2RobustCellClipper() {
    this.options = new Options();
  }

  /** Constructs an S2RobustCellClipper with the given options. */
  public S2RobustCellClipper(Options options) {
    this.options = options;
  }

  /** Returns the options for this clipper. */
  public Options options() {
    return options;
  }

  /** Sets the cell to clip to and clears the crossing list. */
  public void startCell(S2Cell cell) {
    this.cell = cell;
    cellCenter = cell.getCenter();
    // The point outside the cell is arbitrary, so we just choose one that's 90 degrees away. Even
    // face cells are only 45 degrees from center to edge.
    outside = cellCenter.ortho();

    clipper.init(cell);
    reset();

    // Precompute edges and normals for the boundary of the cell.
    for (int k = 0; k < 4; ++k) {
      boundaries[k].set(cell.getVertex(k), cell.getVertex(k + 1));
      uvcoords[k] = cell.getUVCoordOfBoundary(Boundary.fromValue(k));
      normals[k] = cell.getEdgeRaw(k);
    }
  }

  /** Returns the current cell being clipped against. Returns null if no cell has been started. */
  public @Nullable S2Cell cell() {
    return cell;
  }

  /**
   * Clips the edge with the given vertices to the current cell. The starting vertex must be equal
   * to the end vertex of the previously clipped edge.
   */
  @CanIgnoreReturnValue
  public RobustClipResult clipEdge(S2Point v0, S2Point v1, boolean connected) {
    return clipEdge(MutableEdge.of(v0, v1), connected);
  }

  /** Clips the edge with the given vertices to the current cell. */
  @CanIgnoreReturnValue
  public RobustClipResult clipEdge(S2Point v0, S2Point v1) {
    return clipEdge(MutableEdge.of(v0, v1), false);
  }

  /**
   * Tests an edge for intersection with the current cell. When {@link Options#enableCrossings()} is
   * true (the default), boundary crossing information is added to the crossings list.
   *
   * <p>Returns a result object which can be examined to see whether the edge intersected the cell
   * and whether each vertex was contained or not.
   *
   * <p>If connected is true, then the clipper assumes the current edge follows the previous edge
   * passed to clipEdge() (i.e. curr.v0 == prev.v1) and will reuse previous computations.
   */
  @CanIgnoreReturnValue
  public RobustClipResult clipEdge(MutableEdge edge, boolean connected) {
    // Clip the edge first. This give us the UV coordinates for it. If we missed the cell face
    // entirely, then we can't intersect the cell so return.
    boolean hit = clipper.clipEdge(edge.a, edge.b, connected);
    if (!hit && clipper.missedFace()) {
      return RobustClipResult.MISS;
    }

    // Check if the vertices of the edge fall within the error region around the cell, and fall back
    // to exact intersection tests if so.
    if (withinUvErrorMargin(clipper.faceUvEdge(), clipper.uvError())) {
      return clipEdgeExactly(edge, clipper.faceUvEdge());
    }

    if (hit) {
      byte out0 = clipper.outcode(0);
      byte out1 = clipper.outcode(1);

      // If both vertices were inside the cell, there's no boundary crossings.
      if (out0 == INSIDE && out1 == INSIDE) {
        containedEdges.add(edge);
        return RobustClipResult.HIT_BOTH;
      }

      // Make sure intersection point(s) aren't too close to the corners. If they are then we have
      // to fall back to exact tests. This prevents false-hits where the edge actually misses the
      // cell, but we rounded to hit it in uv space:
      //
      //    xyz        uv
      //
      //        ╲       ╲
      //    ┌───┐╲ => ┌──╲┐
      //    │   │     │   │
      //    └───┘     └───┘
      //
      boolean tooClose = false;
      tooClose |= tooCloseToCorner(clipper.clippedUvEdge().v0, out0, clipper.clipError());
      tooClose |= tooCloseToCorner(clipper.clippedUvEdge().v1, out1, clipper.clipError());

      if (tooClose) {
        return clipEdgeExactly(edge, clipper.faceUvEdge());
      }

      // The intersection was far enough from the corner that it's reliable.
      if (options().enableCrossings()) {
        addCrossing(edge, clipper.clippedUvEdge().v0, out0);
        addCrossing(edge, clipper.clippedUvEdge().v1, out1);
      }
      return RobustClipResult.hit(out0, out1);
    }

    // Now we think the edge missed the cell completely. We have to check to make sure it wasn't a
    // false-miss where the edge grazed a corner of the cell within the error margin:
    //
    //      xyz       uv
    //
    //        ╲          ╲
    //      ┌──╲┐ => ┌───┐╲
    //      │   │    │   │
    //      └───┘    └───┘
    //
    // We need to do a little more testing to be sure that we can discard the edge. A lot of the
    // time an edge will be entirely to one side of a cell boundary, which we can detect very
    // quickly (two dot products in the best case, 8 in the worst), and reject the edge. In
    // practice, when clipping down from a larger cell to a smaller one, this tends to be most
    // edges, so falling back to exact tests isn't as bad as it sounds. A UV line-point distance
    // predicate might be faster still though.
    return clipEdgeExactly(edge, clipper.faceUvEdge());
  }

  /**
   * Tests an edge for intersection with the current cell. See {@link #clipEdge(MutableEdge,
   * boolean)} for details. Does not assume the edge is connected to the previous edge.
   */
  @CanIgnoreReturnValue
  public RobustClipResult clipEdge(MutableEdge edge) {
    return clipEdge(edge, false);
  }

  /**
   * Tests an edge for intersection with the current cell. See {@link #clipEdge(MutableEdge,
   * boolean)} for details. Does not assume the edge is connected to the previous edge.
   */
  public RobustClipResult clipEdge(EdgeAndIdChain edge) {
    MutableEdge e = new MutableEdge();
    e.set(edge.start(), edge.end());
    return clipEdge(e, false);
  }

  /** Clears the crossings list and other internal state. Call before passing a new set of edges. */
  public void reset() {
    crossings.clear();
    crossingEdges.clear();
    containedEdges.clear();
  }

  /**
   * Returns the list of crossings since the last call to {@link #reset()}. Crossings are returned
   * sorted (see ordering information on {@link Crossing}. Note that the Crossings and the
   * PooledList containing them are reused, so don't hold references to them beyond the next call to
   * {@link #reset()} or {@link #startCell(S2Cell)}.
   */
  public PooledList<Crossing> getCrossings() {
    if (needSorting) {
      sortCrossings();
    }

    return crossings;
  }

  /**
   * Determines whether or not the cell boundary is contained by the shape defined by the edges that
   * have been clipped to the cell (before they were clipped).
   *
   * <p>After clipping all the edges to the cell, if there are no crossings found, the shape must
   * either contain all of the cell boundary, or contain none of it. This method may be called to
   * determine which is which.
   *
   * <p>The tests used to determine the boundary containment are exact and correspond exactly to the
   * boundary used to clip edges.
   *
   * <p>REQUIRES: No crossings of the boundary were found.
   */
  public boolean isBoundaryContained(boolean containsCenter) {
    Preconditions.checkState(crossings.isEmpty());

    // If there are no contained edges in the cell, there's no edge that could toggle the interior
    // state between the center and the boundary, thus the boundary containment is the same as the
    // center containment.
    if (containedEdges.isEmpty()) {
      return containsCenter;
    }

    // Otherwise, we have one or more contained edges. We can draw a line from the cell center to a
    // point that we know is outside the cell and cross edges to determine its containment.
    //
    // If we cross any edges, they must be inside the boundary of the cell by construction, so we've
    // eliminated any edges that might flip the interior state between the boundary and the final
    // point we chose. Thus the final containment state of the point is the same as the boundary
    // containment.
    boolean inside = containsCenter;
    EdgeCrosser crosser = new EdgeCrosser(cellCenter, outside);
    for (MutableEdge edge : containedEdges) {
      if (crosser.edgeOrVertexCrossing(edge.a, edge.b)) {
        inside = !inside;
      }
    }
    return inside;
  }

  /**
   * Takes a point in XYZ coordinates and projects it onto boundary k of the cell in UV coordinates.
   * Since cells are rectangles in UV space, this is equivalent to just taking the U or V coordinate
   * depending on the boundary.
   */
  private double projectToBoundary(int k, S2Point p) {
    R2Vector uvp = new R2Vector();
    S2Projections.validFaceXyzToUv(clipper.clipFace(), p, uvp);
    return getCoordByBoundary(k, uvp);
  }

  /**
   * Adds a crossing of the given boundary, maintaining sorted order. Takes an edge, the boundary
   * that it crossed (0-3) and intercept which is the point on the U or V axis where the edge
   * crossed the boundary.
   *
   * <p>The crossing type is also specified, but if it is UNKNOWN, then it's computed using a Sign
   * test with the cell vertices.
   */
  private void addCrossing(
      MutableEdge edge, Boundary boundary, double intercept, CrossingType crossingType) {
    // If v0 is to the left of the edge (inside the cell) the edge must be exiting across this
    // boundary.
    if (crossingType == CrossingType.UNKNOWN) {
      MutableEdge cellEdge = boundaries[boundary.value];
      int sign = S2Predicates.sign(cellEdge.a, cellEdge.b, edge.a);
      assert sign != 0;
      crossingType = (sign > 0) ? CrossingType.OUTGOING : CrossingType.INCOMING;
    }

    int edgeIndex = crossingEdges.size();
    crossingEdges.add(edge);

    addCrossing(boundary, crossingType, uvcoords[boundary.value], intercept, edgeIndex);
  }

  /** As above, but with a default crossing type of UNKNOWN. */
  private void addCrossing(MutableEdge edge, Boundary boundary, double intercept) {
    addCrossing(edge, boundary, intercept, CrossingType.UNKNOWN);
  }

  /**
   * Adds a crossing as above, but computes the boundary and intercept points directly from a UV
   * vertex and the boundary it crossed when clipped.
   */
  private void addCrossing(MutableEdge edge, R2Vector uv, int outcode) {
    // If the point is inside or outside the clip region, ignore it.
    if (outcode == INSIDE || outcode == OUTSIDE) {
      return;
    }

    // TODO(torrey): Rather than switch here on OUTCODE only to pass a Boundary enum to addCrossing,
    // we should look at how 'boundary' is used, and pass the things addCrossing uses instead. We
    // can copy the required values to a record and pass that, so each function call passes a single
    // pointer. Update the C++ implementation to match.

    switch (outcode) {
      case BOTTOM:
        addCrossing(edge, BOTTOM_EDGE, uv.x);
        break;

      case RIGHT:
        addCrossing(edge, RIGHT_EDGE, uv.y);
        break;

      case TOP:
        addCrossing(edge, TOP_EDGE, uv.x);
        break;

      case LEFT:
        addCrossing(edge, LEFT_EDGE, uv.y);
        break;

      default:
        throw new IllegalStateException("Unknown outcode: " + outcode);
    }
  }

  /** Adds a crossing with the given values to the boundary set and marks it as needing sorting. */
  private void addCrossing(
      Boundary boundary,
      CrossingType crossingType,
      double uvcoord,
      double intercept,
      int edgeIndex) {
    crossings.add().set(boundary, crossingType, uvcoord, intercept, edgeIndex);
    needSorting = true;
  }

  /**
   * Returns true if a point is within maxError of corner along a boundary. The boundary segment is
   * specified by opcode (see {@link UVEdgeClipper}, and the position along the segment is taken
   * from the uv point given.
   *
   * <p>E.g. the top boundary of the cell is constant in the V axis, so the intercept position is
   * taken to be the u coordinate of the point which is checked against the end points of the top
   * cell boundary.
   */
  private boolean tooCloseToCorner(R2Vector uv, int outcode, double maxError) {
    // Point isn't on the boundary, so it wasn't clipped, ignore it.
    if (outcode == INSIDE || outcode == OUTSIDE) {
      return false;
    }

    // In both cases below we check that the intercept point is far enough in from the cell corners
    // that we can trust the intersection result. It can't lie within the forbidden regions:
    //
    //      maxError        maxError
    //        ├───┤           ├───┤
    //   lo() ┌────────•──────────┐ hi()
    //   ─────────┤           ├─────────
    //   forbidden             forbidden
    //
    // uvcoords are [y.lo, x.hi, y.hi, x.lo], also known as [vMin, uMax, vMax, uMin]
    boolean tooClose = false;
    if (outcode == BOTTOM || outcode == TOP) {
      // Top and bottom boundaries, compare against u bounds.
      double intercept = uv.x();
      tooClose |= (intercept - maxError <= uvcoords[3]); // x.lo or uMin
      tooClose |= (intercept + maxError >= uvcoords[1]); // x.hi or uMax
    } else {
      // Left and right boundaries, compare against v bounds.
      double intercept = uv.y();
      tooClose |= (intercept - maxError <= uvcoords[0]); // y.lo or vMin
      tooClose |= (intercept + maxError >= uvcoords[2]); // y.hi or vMax
    }
    return tooClose;
  }

  /** Returns true if either vertex of a UVEdge is within maxError of the current cell boundary. */
  @VisibleForTesting
  boolean withinUvErrorMargin(R2Edge uvEdge, double maxError) {
    // We have to test each x and y coordinate against two boundaries. Expand out into 4 element
    // arrays explicitly to help the vectorizer.
    // uvcoords are [y.lo, x.hi, y.hi, x.lo], also known as [vMin, uMax, vMax, uMin]
    double[] coord0 = {uvEdge.v0.y(), uvEdge.v0.x(), uvEdge.v0.y(), uvEdge.v0.x()};
    double[] coord1 = {uvEdge.v1.y(), uvEdge.v1.x(), uvEdge.v1.y(), uvEdge.v1.x()};

    boolean inMargin = false;
    for (int i = 0; i < 4; ++i) {
      inMargin |= (abs(coord0[i] - uvcoords[i]) <= maxError);
      inMargin |= (abs(coord1[i] - uvcoords[i]) <= maxError);
    }

    return inMargin;
  }

  /**
   * Performs exact crossing tests between an edge and the boundaries of the cell to determine
   * whether it actually crosses or not. This is done using the great circle normals that define the
   * boundary of the cell (as given by S2Cell.getEdgeRaw()).
   *
   * <p>We use these rather than testing using the cell vertices in order to ensure that results are
   * consistent even when we have different size cells across a boundary. The cell corners aren't
   * actually guaranteed to be on the boundary, rather they're the nearest representable S2Point to
   * the true cell corner.
   *
   * <p>Returns clip result as documented in clipEdge() indicating whether the edge hit the cell and
   * which vertices were contained.
   */
  private RobustClipResult clipEdgeExactly(MutableEdge edge, R2Edge uvEdge) {
    // If we're calling this function, then we were too close to the boundary of the cell to rely on
    // our intersection results (or we want to verify that an edge did indeed entirely miss the
    // cell.
    //
    // We can fall back to using exact predicates to determine which boundaries we crossed, and
    // re-use the UV intercepts to get consistent crossing points.
    //
    // BoundarySign() is just a dot product with the exact boundary normals with a perturbation to
    // break ties in the event that we land exactly on a boundary.
    //
    // We classify each of the two edge vertices with respect to the four cell boundaries, which
    // tells us which of the great circles we crossed.
    //
    // This is effectively a spherical version of Cohen-Sutherland, except that we use dot products
    // to determine axis crossings.
    //
    boolean allGreaterThan00 = true;
    boolean allGreaterThan01 = true;
    byte[] sign0 = new byte[4];
    byte[] sign1 = new byte[4];

    for (int k = 0; k < 4; ++k) {
      sign0[k] = (byte) boundarySign(k, edge.a);
      sign1[k] = (byte) boundarySign(k, edge.b);

      // If both vertices are on the negative side of the same cell boundary, the edge can't
      // intersect the cell so we're done.
      if (sign0[k] < 0 && sign1[k] < 0) {
        return RobustClipResult.MISS;
      }

      allGreaterThan00 &= (sign0[k] > 0);
      allGreaterThan01 &= (sign1[k] > 0);
    }

    // If the edge was properly inside all the boundaries, then it's contained so there's no
    // crossings to add, just indicate both vertices were inside.
    if (allGreaterThan00 && allGreaterThan01) {
      containedEdges.add(edge);
      return RobustClipResult.HIT_BOTH;
    }

    // We know we crossed at least one boundary, let's determine the details.
    RobustClipResult result = RobustClipResult.MISS;
    for (int k = 0; k < 4; ++k) {
      int kNext = (k + 1) % 4;
      int kPrev = (k + 3) % 4;

      // We have a crossing when the sign changed between the vertices.
      if (sign0[k] != sign1[k]) {
        // Get the normal for the crossed boundary, and its two neighbors, kPrev and kNext. Check
        // that the crossing of k occurred between the two neighbors boundaries, meaning the
        // intersection is in bounds of the cell:
        //                                 v1
        //                          ┊   k  ╱ ┊
        //                        ┄┄┌─────•──┐┄┄
        //                   kNext  │→   ╱  ←│ kPrev
        //                              v0
        //
        int signNext =
            S2Predicates.circleEdgeIntersectionSign(edge.a, edge.b, normals[k], normals[kNext]);

        int signPrev =
            S2Predicates.circleEdgeIntersectionSign(edge.a, edge.b, normals[k], normals[kPrev]);

        // If the signs are both >= 0 or <= 0 then we crossed within the lune or anti-lune formed by
        // kNext and kPrev, so add the crossing.
        if ((signNext >= 0 && signPrev >= 0) || (signNext <= 0 && signPrev <= 0)) {
          if (options().enableCrossings()) {
            // Now we know that:
            //   1. The edge crossed this cell boundary
            //   2. It did so within range of the neighboring boundary planes.
            //
            // Therefore we have a crossing.
            //
            // Compute the intercept point of the edge on this boundary. Since we're down in the
            // region where we're not sure whether the edge crossed in UV space, it's possible to
            // get a situation where the edge crosses the cell boundary in XYZ but becomes
            // horizontal or vertical in UV, which would cause the regular linear interpolation
            // to fail.
            //
            // Instead, when the edge is horizontal or vertical, we'll compute the interception
            // point manually with cross products and project it onto the boundary. As long as we're
            // consistent in how we do this so that neighboring cells see the same intercept, it
            // will give correct results.
            double intercept;
            if (uvEdge.v0.x() != uvEdge.v1.x() && uvEdge.v0.y() != uvEdge.v1.y()) {
              R2Vector uvEdgeClipResult = new R2Vector();
              clipper.r2Clipper.clip(uvEdge, (byte) (1 << k), uvEdgeClipResult);
              intercept = getCoordByBoundary(k, uvEdgeClipResult);
            } else {
              S2Point intersection =
                  robustCrossProd(
                          robustCrossProd(edge.a, edge.b).normalize(), normals[k].normalize())
                      .normalize();

              // We might get the antipodal intersection depending on which way we crossed, so flip
              // the sign if needed.
              if (sign0[k] < 0) {
                intersection = intersection.neg();
              }

              intercept = projectToBoundary(k, intersection);
            }
            Boundary boundary = Boundary.fromValue(k);
            if (sign0[k] > sign1[k]) {
              addCrossing(edge, boundary, intercept, CrossingType.OUTGOING);
            } else {
              addCrossing(edge, boundary, intercept, CrossingType.INCOMING);
            }
          }

          result = RobustClipResult.hit(allGreaterThan00, allGreaterThan01);
        }
      }

      assert (sign0[k] != 0 && sign1[k] != 0);
    }

    return result;
  }

  /**
   * Returns the U or V coordinate of P, whichever is not constant in the axis the given cell
   * boundary lies in.
   */
  private double getCoordByBoundary(int boundary, R2Vector p) {
    if (boundary % 2 == 0) {
      // Bottom (0) and Top (2) are constant in V, so take the U coordinate.
      return p.x();
    } else {
      // Right (1) and Left (3) are constant in U, so take the V coordinate.
      return p.y();
    }
  }

  /**
   * Determines which side of a cell boundary segment a given point is on.
   *
   * <p>The convention is that the positive side of a cell boundary points towards the cell center,
   * with the negative side pointing away. Thus, when a point is in the cell, it must be on the
   * non-negative side of all four edges:
   *
   * <pre>
   *
   *                             ┊   -   ┊
   *                           ┄┄┌───────┐┄┄
   *                             │   +   │
   *                            -│+  •  +│-
   *                             │   +   │
   *                           ┄┄└───────┘┄┄
   *                             ┊   -   ┊
   * </pre>
   *
   * The vertices of a cell are not guaranteed to be on the cell boundary (instead they're the
   * closest representable S2Point to the corner). Even when cells share a boundary, this
   * approximation can cause us to compute slightly different normal vectors on either side, giving
   * inconsistent results across the boundary.
   *
   * <p>So, instead of using the vertices, this function relies on the boundary edge normals
   * produced by getEdgeRaw(), which are computed exactly, and naturally face inward towards the
   * cell interior.
   *
   * <p>We use an exact predicate to determine the sign of the dot product, so the result is correct
   * for all S2Points.
   *
   * <p>If the sign is exactly 0, a sign is chosen to break the tie. The resulting sign is
   * arbitrary, but it's guaranteed to be consistent across all cell boundaries on the same great
   * circle.
   *
   * <p>For convenience, the argument k is reduced modulo 4 to the range [0..3].
   *
   * <p>REQUIRES: p.Norm2() <= 2
   *
   * <pre>
   * Returns:
   *   -1 if p is on the negative side of edge k
   *   +1 if p is on the positive side of edge k
   * </pre>
   */
  private int boundarySign(int k, S2Point p) {
    assert p.norm2() <= 2;
    S2Point normal = normals[k % 4];
    int sign = S2Predicates.signDotProd(normal, p);
    if (sign == 0) {
      // Zero sign means we landed exactly on the great circle, and the user requested perturbation
      // so we need to choose a side.
      //
      // Given the following:
      //   * getEdgeRaw(k) defines a great circle, which bisects the sphere.
      //   * For any cell with a boundary segment touching that circle:
      //       * Cells with a center on the positive side are on side 'A'
      //       * Cells with a center on the negative side are on side 'B', then:
      //
      // We want cells on side A to see one consistent sign, which must be opposite the cells on
      // side B.
      //
      // The normal vector itself can't be all zeros, and it is naturally equal-but-opposite across
      // the cell boundary, so we can take the sign of its first non-zero component as the sign to
      // break the tie.
      for (int i = 0; i < 3; ++i) {
        if (normal.get(i) != 0) {
          return normal.get(i) > 0 ? +1 : -1;
        }
      }
    }
    return sign;
  }

  /** Sort the crossings using exact predicates if necessary, and remove duplicates. */
  private void sortCrossings() {
    needSorting = false;
    Comparator<Crossing> crossingComparator = new ExactCrossingComparator();
    crossings.sort(crossingComparator);

    // We may have had duplicate crossings at the exact same point, which can only be caused by a
    // reverse duplicate edge or a collinear vertex. In either case the crossings won't cause a net
    // chance in the interior-ness of the boundary: they cancel each other out and can be ignored.
    // Scan the array and remove such duplicates.

    // Erasing in-place from a vector isn't terribly efficient but our vectors should be small, and
    // we'll almost never have equal crossings, so this almost always turns into a single pass
    // through the array which then does nothing.
    for (int i = 0; i < crossings.size() - 1; ++i) {
      if (crossingComparator.compare(crossings.get(i), crossings.get(i + 1)) == 0) {
        crossings.remove(i + 1);
        crossings.remove(i);
        --i;
      }
    }
  }

  /** The result of clipping an edge to a cell. */
  public enum RobustClipResult {
    MISS(0b000),
    HIT_NONE(0b100),
    HIT_V0(0b101),
    HIT_V1(0b110),
    HIT_BOTH(0b111);

    /**
     * Encoded status that tells us whether an edge hit a cell and containment information for each
     * of its vertices.
     */
    private final int clipCode;

    /** Creates a RobustClipResult with the given clip code. */
    RobustClipResult(int clipCode) {
      this.clipCode = clipCode;
    }

    /**
     * Returns a result where the edge intersected the cell and matching the given containment flags
     * for each vertex.
     */
    public static RobustClipResult hit(boolean v0Contained, boolean v1Contained) {
      return v1Contained //
          ? (v0Contained ? HIT_BOTH : HIT_V1) //
          : (v0Contained ? HIT_V0 : HIT_NONE);
    }

    /**
     * Returns a result where the edge intersected the cell and matching the given outcodes for each
     * vertex.
     */
    public static RobustClipResult hit(byte out0, byte out1) {
      return hit(out0 == INSIDE, out1 == INSIDE);
    }

    /** Return true if this result indicates the edge hit the cell, false otherwise. */
    public boolean hit() {
      return this != MISS;
    }

    /**
     * Returns true if this result indicates the start vertex of the clipped edge was contained by
     * the cell.
     */
    public boolean v0Inside() {
      return (clipCode & 0b001) != 0;
    }

    /**
     * Returns true if this result indicates the end vertex of the clipped edge was contained by the
     * cell.
     */
    public boolean v1Inside() {
      return (clipCode & 0b010) != 0;
    }

    @Override
    public String toString() {
      StringBuilder sb = new StringBuilder();
      if (hit()) {
        sb.append("Hit ");

        if (v0Inside() && v1Inside()) {
          sb.append("(Both)");
        } else {
          sb.append("(v0: ").append(v0Inside()).append(", v1: ").append(v1Inside()).append(")");
        }
      } else {
        sb.append("Miss");
      }
      return sb.toString();
    }
  }

  /** The direction of a boundary crossing. */
  public enum CrossingType {
    // Invalid crossing type only used as a sentinel value.
    UNKNOWN,
    INCOMING,
    OUTGOING;
  }

  /**
   * A mutable container to hold cell boundary crossing information. This stores the intercept point
   * along an implied axis (U or V), as well as whether the edge is entering or exiting the cell at
   * that point.
   *
   * <p>Crossings are sorted first by boundary, and then intercept point. Within each boundary
   * segment, the crossings are ordered so that they're ordered counter-clockwise around the cell:
   *
   * <pre>
   *
   *            ←
   *        ┌───────┐
   *        │   T   │
   *      ↓ │L     R│↑
   *        │   B   │
   *        └───────┘
   *            →
   * </pre>
   *
   * <p>This means that crossings on the bottom and right edges are sorted in ascending order, and
   * crossings on the top and left are sorted in descending order.
   *
   * <p>Crossings are stored in a PooledList so that they can be reused.
   */
  public static class Crossing {

    /** Value in the non-constant axis where edge crosses. */
    public double intercept;

    /** Value in the constant axis where edge crosses. */
    public double coord;

    /** Which cell boundary is crossed. */
    public Boundary boundary;

    /** Whether the edge is entering or exiting the cell. */
    public CrossingType crossingType;

    /** Index into the edges list for the edge that caused this crossing. */
    public int edgeIndex;

    /** Creates a new empty Crossing. */
    public Crossing() {}

    /** Sets the given parameters on this Crossing. */
    public void set(
        Boundary boundary,
        CrossingType crossingType,
        double coord,
        double intercept,
        int edgeIndex) {
      this.intercept = intercept;
      this.coord = coord;
      this.boundary = boundary;
      this.crossingType = crossingType;
      this.edgeIndex = edgeIndex;
    }

    @Override
    public String toString() {
      StringBuilder sb = new StringBuilder();
      if (crossingType == CrossingType.OUTGOING) {
        sb.append("Outgoing");
      }
      if (crossingType == CrossingType.INCOMING) {
        sb.append("Incoming");
      }

      sb.append(" on boundary ").append(boundary).append(" -- ");
      sb.append(intercept).append("@").append(coord);
      sb.append(" ").append("(edge ").append(edgeIndex).append(")");
      return sb.toString();
    }

    /** Returns true if the current values of the given crossing and this one are equal. */
    public boolean isEqualTo(Crossing other) {
      return this.crossingType == other.crossingType
          && this.intercept == other.intercept
          && this.boundary == other.boundary
          && this.coord == other.coord;
    }
  }

  /** Options for the S2RobustCellClipper. */
  public static class Options {
    private boolean enableCrossings = true;

    /**
     * When enableCrossings is true (the default), any cell boundary crossings are added to the
     * crossings list. Otherwise, only hit information is returned from clipEdge and we promise not
     * to allocate memory when called.
     */
    @CanIgnoreReturnValue
    public Options setEnableCrossings(boolean flag) {
      this.enableCrossings = flag;
      return this;
    }

    /** Returns true if crossings are enabled. */
    public boolean enableCrossings() {
      return enableCrossings;
    }
  }

  /** A Crossing comparator that uses exact tests to compare the ordering of A and B. */
  private class ExactCrossingComparator implements Comparator<Crossing> {

    @Override
    public int compare(Crossing a, Crossing b) {
      // If the crossings aren't on the same boundary, then we can just compare boundaries directly.
      if (a.boundary.value < b.boundary.value) {
        return -1;
      }
      if (a.boundary.value > b.boundary.value) {
        return +1;
      }

      // Each intercept can have MAX_ERROR on it, so we need to check that they're at least two
      // times that far apart before relying on the intercepts.
      if (abs(a.intercept - b.intercept) > 2 * MAX_ERROR) {
        switch (a.boundary) {
          case BOTTOM_EDGE: // fallthrough
          case RIGHT_EDGE: // Bottom and right boundaries sort ascending.
            return (a.intercept < b.intercept) ? -1 : +1;

          case TOP_EDGE: // fallthrough
          case LEFT_EDGE: // Top and left boundaries sort descending.
            return (a.intercept > b.intercept) ? -1 : +1;
        }
        throw new IllegalStateException("Invalid boundary: " + a.boundary);
      }

      // Intercept points are on the same boundary and too close, so we need to use exact predicates
      // to compare them. Grab the normal vectors for this boundary segment and the previous one.
      int k = a.boundary.value;
      S2Point norm = normals[k];
      S2Point prev = normals[(k + 3) % 4];

      MutableEdge edgeA = crossingEdges.get(a.edgeIndex);
      MutableEdge edgeB = crossingEdges.get(b.edgeIndex);

      // Incoming edges have v0 on the negative side of the boundary, and we want v0 to be positive,
      // v1 negative for consistency when calling the predicate. Flip the order if needed, without
      // actually modifying the edge.
      if (a.crossingType == CrossingType.INCOMING) {
        edgeA = MutableEdge.of(edgeA.b, edgeA.a);
      }

      if (b.crossingType == CrossingType.INCOMING) {
        edgeB = MutableEdge.of(edgeB.b, edgeB.a);
      }

      // Order the edges relative to the previous boundary. The corner between the current (k) and
      // previous (k-1) boundary segments is our zero point, so we'll naturally order crossings in
      // increasing intercept order along this boundary, regardless of which boundary segment it is:
      //
      //                        Increasing intercept
      //                              <----
      //                             B     D
      //                             •  k  •
      //                            ┌─╲───╱─┐
      //                        k+1 │  • •  │ k-1
      //                               A C
      //
      return S2Predicates.circleEdgeIntersectionOrdering(
          edgeA.a, edgeA.b, edgeB.a, edgeB.b, norm, prev);
    }
  }
}
