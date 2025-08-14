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

import static com.google.common.geometry.S2Cell.Boundary.BOTTOM_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.LEFT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.RIGHT_EDGE;
import static com.google.common.geometry.S2Cell.Boundary.TOP_EDGE;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.geometry.S2Cell.Boundary;
import com.google.common.geometry.S2RobustCellClipper.Crossing;
import com.google.common.geometry.S2RobustCellClipper.CrossingType;
import com.google.common.geometry.primitives.CountingBitset;
import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import java.util.List;
import java.util.function.BiFunction;
import org.jspecify.annotations.Nullable;

/**
 * A class for tracking a shape in an S2ShapeIndex and detecting when all the pieces of the shape
 * have been visited.
 *
 * <p>This works by tracking segments of S2Cell boundaries in IJ coordinates that are contained by a
 * shape. These segments must always match up across cell boundaries. By "toggling" segments on and
 * off as we visit them, we can arrange for them all to cancel once we have seen all the interior
 * pieces of a shape.
 *
 * <p>Once we have seen every chain at least once, and there are no outstanding boundary segments
 * contained by the shape, then we have visited all of the shape.
 *
 * <p>1. Usage
 *
 * <pre>
 * Instantiate an S2ShapeTracker instance and initialize it with an S2Shape:
 * {@code
 *   S2ShapeTracker tracker = new S2ShapeTracker(shape);
 * } or {@code
 *   tracker.reset(shape);
 * }
 *
 * Then, as you iterate through an index, mark individual chains as seen when you come across them:
 *
 * {@snippet :
 * tracker.markChain(chain);
 * }
 *
 * This method is cheap to call (O(1) complexity) and idempotent so it's safe to call multiple times
 * per chain.
 *
 * <p>Use S2RobustCellClipper to generate boundary crossings for the current cell and process them
 * with shape tracker:
 *
 * {@snippet :
 * S2RobustCellClipper clipper = new S2RobustCellClipper();
 *
 * clipper.startCell(cell);
 * for (each edge from shape in cell) {
 *   clipper.clipEdge(edge);
 * }
 *
 * tracker.processCrossings(cell, clipper.getCrossings());
 * }
 *
 * <p>If the cell has no edges but is contained by the shape, use addCellBoundary() to account for
 * the entire boundary of the cell in the tracker:
 *
 * {@snippet :
 * tracker.addCellBoundary(cell);
 * }
 *
 * The current of the shape can be queried via the finished() method:
 *
 * {@snippet :
 * if (tracker.finished()) {
 *   ... report results
 * }
 * }
 *
 * <p>The tracker can then be reused for a new shape by calling reset().
 *
 * <p>The rest of this documentation covers the theory and details behind interior tracking. If you
 * just want to get started tracking shapes, feel free to skip.
 *
 * <p>2. Theory
 *
 * <p>For bulk queries in general, and containment in particular, we need to know when we've seen
 * all of a shape so that we know when it's safe to make a final decision regarding the shape and
 * return results to the user.
 *
 * <p>Since shape edges can span multiple cells, it's not sufficient to just track that we've seen
 * all the edges in a shape. If we imagine a shape with one long edge that crosses many cells, the
 * first cell we visit we will mark the edge as seen and report results. But the edge continued
 * outside the cell and we may have missed information relevant to the query result.
 *
 * <p>The same reasoning applies to 2D shapes. Imagine a ring circling the Earth around the equator.
 * If we create another shape by reversing the order of the first, the two shapes won't be contained
 * even though they have the same exact vertices. The interior of each escapes the interior of the
 * other.
 *
 * <p>But, importantly, we have no way to know this directly, since all the vertices are coincident,
 * none of the vertices will test as outside the shapes nor will edges cross.
 *
 * <p>This interior "escapage" is only detectable by knowing that there is "more interior" pending
 * when we finish processing all the edges.
 *
 * <p>By tracking portions of the cell boundary that are contained by a shape and removing them when
 * we visit the opposite side of the boundary, we can track what pieces of a shape are still
 * unresolved.
 *
 * <p>Multi-points are trivial, they can't cross a cell boundary (or any edge) and are each their
 * own chain so we merely have to check that we've seen all the chains in the shape.
 *
 * <p>Polylines are slightly more complex because they can cross a boundary at a point. These
 * crossings are 0 dimensional however so we only have to match crossings and don't have to worry
 * about extents.
 *
 * <p>Finally, polygons can contain a 1 dimensional segment of the cell boundary, so we have to
 * track both the presence of a crossing and its extent.
 *
 * <p>It's entirely legitimate (and common) for us to have a large cell bounded by smaller ones, so
 * we may see a segment of a boundary that is then visited on its opposite side in smaller pieces
 * (or vice-versa):
 *
 * <pre>
 *             ┌─────┬─┐     ┌───────┐     ┌───────┐    ┌───────┐
 *             │     ├─┘  →  │     ┌─┤  →  │       │ →  │       │
 *             │     │       │     ├─┘     │     ┌─┤    │       │
 *             └─────┘       └─────┘       └─────┴─┘    └───────┘
 * </pre>
 *
 * <p>Thus, we need a way to track boundary intervals that also allows us to build and cancel in
 * small pieces. We can't merely store an interval (A,B) and expect to cancel it in a single piece
 * with another interval (B,A).
 *
 * <p>2.1 The Problem
 *
 * <p>Given a grid of non-overlapping cells (which any S2ShapeIndex is), we'd like to devise a
 * system where we can track boundary segments efficiently. We know a few things that we can use to
 * our advantage:
 *
 * <ol>
 *   <li>If we're destined to see all of a shape, then any time it crosses a cell boundary, we will
 *       always see both sides of the boundary due to the cell padding we use when building indices.
 *   <li>We will always see exactly the same crossing points in cells on either side of the boundary
 *       due to the exactness in UVEdgeClipper.
 *   <li>In IJ space, cells become axis-aligned squares with their boundary segments lying along
 *       constant I or J values.
 * </ol>
 *
 * <p>Since we can rely on seeing the same section of cell boundary twice the problem is one of
 * pairing and canceling matching segments.
 *
 * <p>2.2 The Closure Grid
 *
 * <p>Towards this goal, let's note that that we can build a consistent notion of sign along cell
 * boundaries. By consistent, we mean that we can assign every boundary segment of every S2 cell a
 * sign value (either +1 or -1) that is opposite across the boundary.
 *
 * <p>Cells have consistently numbered edges, ordered CCW around the boundary, numbered 0 to 3. If
 * we define edges 0 and 1 to be positive and edges 2 and 3 to be negative, then we can see that if
 * we tessellate the cells the sign values are equal-and-opposite across boundaries:
 *
 * <pre>
 *             ┌───────┐      ┌───────┐      ┌───────┬───────┐
 *             │   2   │      │   -   │      │   -   │   -   │
 *             │3     1│  →   │-     +│  →   │-     +│-     +│
 *             │   0   │      │   +   │      │   +   │   +   │
 *             └───────┘      └───────┘      ├───────┼───────┤
 *                                           │   -   │   -   │
 *                                           │-     +│-     +│
 *                                           │   +   │   +   │
 *                                           └───────┴───────┘
 * </pre>
 *
 * <p>Because cell orientations on a given face are constant, this holds true even when there are
 * cells of different sizes on each side of the boundary. One side will always have the positive
 * boundaries and the other side the negative.
 *
 * <p>This works automatically within a given face, but is it also true between faces? The naive
 * answer is "almost". Let's look at the face cells with sign values assigned as described above:
 *
 * <pre>
 *                 ┌─────────────────────┐
 *                 │                 ┌───────┬───────┐
 *                 │                 │J  -   │I  -   │
 *                 │              ┌╴ │-  4  +│-  5  +│ ┐
 *                 │                 │   +  I│   +  J│ │
 *                 │         ┌───────┼───────┼───────┘ │
 *                 │         │J  -   │I  -   │         │
 *                 │     ┌╴  │-  2  +│-  3  +│ ╶┘      │
 *                 │         │   +  I│   +  J│         │
 *                 │ ┌───────┼───────┼───────┘         │
 *                 │ │J  -   │I  -   │                 │
 *                 └ │-  0  +│-  1  +│ ╶┘              │
 *                   │   +  I│   +  J│                 │
 *                   └───────┴───────┘                 │
 *                               └─────────────────────┘
 * </pre>
 *
 * <p>Touching face edges that have the same sign (and are thus incompatible) are connected by
 * lines. Interestingly all the "interior" face boundaries in the diagram are consistent.
 *
 * <p>The edge of a face at a maximum IJ coordinate (S2Projections.LIMIT_IJ == 2^30) spatially
 * aliases the edge of the adjacent face at value 0. We need to break this tie so that intervals
 * along the boundaries between faces are canceled properly.
 *
 * <p>The exact choice in breaking ties doesn't matter as long as it's consistent, so we take the
 * convention that maximum values alias back down to the minimum value on the adjacent face and we
 * do our interval accounting there. This fits nicely with the mental model of faces being "half
 * open" and not containing their upper boundaries.
 *
 * <p>As an example, the maximum J coordinate on face 0 ends up being stored using the minimum I
 * coordinate of face 2, and we flip the values of the intervals so that they match the coordinate
 * system on face two.
 *
 * <p>Looking at the diagram above, all of our incompatible face edges are at a maximum J value.
 * Since we have to wrap those boundaries anyways, while we're doing that, we can negate the sign as
 * well and we now have a fully consistent signed boundary system.
 *
 * <p>2.3 Summing Interval Lengths
 *
 * <p>Now that we have a sign convention that naturally cancels across boundaries, we could track
 * polygon closure by finding segments of the grid that are contained (by finding cell boundary
 * crossings), and summing their lengths with appropriate signs. Since we'll see both sides of a
 * boundary, we'll always add a given length and subtract it out later. The lengths naturally
 * cancel.
 *
 * <p>To do this, take a section of boundary that contains one or more crossings:
 *
 * <pre>
 *                                    ←
 *        •────────•━━━━•───────────•━━━━•─────────────•━━━━━━━━━━━━•
 *        X        A    B           C    D             E            Y
 *                                    →
 * </pre>
 *
 * <p>When we're traversing left-to-right, the first contained interval we'll see is (A,B). (B,C) is
 * not contained so we ignore it, then we take (C,D) and (E,Y).
 *
 * <p>We'll maintain a sum, and for each contained interval we'll add the first point and subtract
 * the second point:
 *
 * <pre>
 *    sum += A-B + C-D + E-Y
 * </pre>
 *
 * <p>Now, when we traverse the opposite side, we'll proceed right-to-left, and we'll see the
 * intervals with their vertices swapped: (Y,E), (D,C), and (B,A).
 *
 * <p>We'll still add the first, and subtract the second:
 *
 * <pre>
 *   sum += Y-E + D-C + B-A
 * </pre>
 *
 * <p>If we look at this together with the first sum:
 *
 * <pre>
 *   sum = A-B + C-D + E-Y + Y-E + D-C + B-A = 0
 * </pre>
 *
 * <p>Each crossing point has been canceled out exactly, leaving a zero total.
 *
 * <p>Note that this left-to-right and right-to-left ordering is a natural consequence of cells all
 * being CCW oriented. Cells with adjacent boundaries will traverse that boundary in opposite
 * directions:
 *
 * <pre>
 *                             ┌───────┬───────┐
 *                             │   ←   │   ←   │
 *                             │↓     ↑│↓     ↑│
 *                             │   →   │   →   │
 *                             ├───────┼───────┤
 *                             │   ←   │   ←   │
 *                             │↓     ↑│↓     ↑│
 *                             │   →   │   →   │
 *                             └───────┴───────┘
 * </pre>
 *
 * <p>This summation also works naturally for splitting up larger intervals, imagine a large cell
 * with a piece of its boundary contained, bordered by smaller cells:
 *
 * <pre>
 *
 *      smaller cells │   0    │    1   │    2   │    3
 *                    │   →    B    →   C    →   D    →   │
 *                    ├─────•━━┷━━━━━━━━┷━━━━━━━━┷━━━━•───┤
 *      large cell    │     X           ←             Y   │
 * </pre>
 *
 * Moving right-to-left first, we'll add (Y,X) first:
 *
 * <pre>
 *   total = Y-X
 * </pre>
 *
 * Now, we process left-to-right (we can process the cells in any order, let's assume 1,2,0,3):
 *
 * <pre>
 *   cell1 - total = Y-X + B-C
 *   cell2 - total = Y-X + B-C + C-D = Y-X + B-D
 *   cell0 - total = Y-X + B-D + X-B = Y-D
 *   cell3 - total = Y-D + D-Y = 0
 * </pre>
 *
 * The values close, and the same is true in any order, let's do 3,2,1,0:
 *
 * <pre>
 *   cell3 - total = Y-X + D-Y = D-X
 *   cell2 - total = D-X + C-D = C-X
 *   cell1 - total = C-X + B-C = B-X
 *   cell0 - total = B-X + X-B = 0
 * </pre>
 *
 * <p>In this example we can see more clearly how we "cancel" pieces of the Y-X interval one at a
 * time, but we're always left with the remaining interval.
 *
 * <p>2.3.1 General Iteration Order Safety
 *
 * <p>The attentive reader might worry that we might accidentally have a crossing of, say, length 2,
 * and pair it with two unrelated length 1 crossings on the other side of the boundary, and get a
 * net zero without actually closing all the intervals. We need to ensure we can't get a false
 * positive for a given shape.
 *
 * <p>First, let's assume all of our points are distinct. This is easy enough to arrange in IJ
 * coordinates, since (face, i, j) specifies a unique point, and cells that are used in an
 * S2ShapeIndex do not overlap.
 *
 * <p>If we look at a closure sum:
 *
 * <pre>
 * (A-B) + (B-C) + (C-A) = 0
 * </pre>
 *
 * <p>Note the general pattern is that we add an element, and then eventually subtract its exact
 * value again to cancel it out. If our points are distinct, we don't have to rely on actually
 * summing them together, we can just track which points we've accumulated with a set, inserting and
 * removing them as we go.
 *
 * <p>We can implicitly track an IJ interval (A,B) by storing its endpoints.
 *
 * <p>It's possible to see multiple of the same point, e.g. if you had multiple crossings within the
 * same IJ coordinate. To support this, we'll use a hashmap of points mapping to a count.
 *
 * <p>Our final method for tracking closure becomes simply:
 *
 * <pre>
 *   hashmap[point]++
 * </pre>
 *
 * <p>And subtracting a point:
 *
 * <pre>
 *   hashmap[point]--
 * </pre>
 *
 * <p>With zero being the default value for a new entry. If the count goes to zero in either case,
 * we remove the point. When the set is empty, we have an exact closure.
 *
 * <p>2.4 Rounding intervals
 *
 * <p>We have to be careful with very small crossings. It's entirely possible for a shape to cross a
 * cell boundary with a width less than one unit in IJ space:
 *
 * <pre>
 *            ╲    ╱
 *     ┌─────┬─•  •┬─────┬─────┐
 *              ╲╱
 * </pre>
 *
 * <p>If we round these intervals naively, we can end up with a zero-length interval in IJ space
 * with equal endpoints. That would result in an interval (A,A) which we would add as A-A = 0, which
 * wouldn't contribute to the tracking information.
 *
 * <p>We can resolve this by rounding intervals _outwards_. We'll round the lower endpoint down, and
 * the higher endpoint up:
 *
 * <pre>
 *     ┌─────┬─• •─┬─────┬─────┐
 *            ←   →
 * </pre>
 *
 * <p>This will give always give us an interval that's at least one unit length, and it works
 * consistently across T points and corners, so that they will properly close:
 *
 * <pre>
 *         │   \   /   │       │   ←  →    │      │    ___    │
 *         ├────•╷•────┤   →   ├────•╷•────┤  →   ├─────┬─────┤
 *         │     │     │       │   ← │→    │      │    ‾│‾    │
 * </pre>
 */
@SuppressWarnings("Assertion")
public class S2ShapeTracker {

  // Tracking state.
  private int dimension = -1;
  private final CountingBitset chainsSeen = new CountingBitset();

  // These are the maps for tracking and pairing matching segments crossing IJ intervals, as
  // described above. For any IJ coordinate, we track the number of points seen at that coordinate.
  // We track points in the I and J axis separately to avoid collisions and reduce the hashmap
  // sizes.

  /**
   * Map from (face,i,j) encoded as a long, to the number of crossings at that that point, for IJ
   * axis 0. (Which of i or j is 0 or 1 depends on the face.)
   */
  private final Long2IntMap points0 = new Long2IntOpenHashMap();

  /**
   * Map from (face,i,j) encoded as a long, to the number of crossings at that that point, for IJ
   * axis 1. (Which of i or j is 0 or 1 depends on the face.)
   */
  private final Long2IntMap points1 = new Long2IntOpenHashMap();

  public S2ShapeTracker() {
    points0.defaultReturnValue(0);
    points1.defaultReturnValue(0);
  }

  public S2ShapeTracker(int dimension, int numChains) {
    points0.defaultReturnValue(0);
    points1.defaultReturnValue(0);
    init(dimension, numChains);
  }

  /** Constructs a new S2ShapeTracker for the given S2Shape. */
  public S2ShapeTracker(S2Shape shape) {
    points0.defaultReturnValue(0);
    points1.defaultReturnValue(0);
    init(shape);
  }

  /** Returns the map of points for the given I or J axis. */
  private Long2IntMap points(int axis) {
    return axis == 0 ? points0 : points1;
  }

  /** Resets state to track a new shape. */
  public void init(int dimension, int numChains) {
    this.dimension = dimension;

    chainsSeen.reset(numChains);
    points0.clear();
    points1.clear();
  }

  /** Resets state for a new shape as above, but takes the shape directly. */
  public void init(S2Shape shape) {
    // Don't track edges for the full and empty shapes, since we'll never actually see an edge
    // despite them having one chain.
    int numChains = shape.numChains();
    if (shape.isEmpty() || shape.isFull()) {
      // TODO(torrey): Add test cases for empty shapes of various dimensions, as well as the full
      // shape. Also in C++. Make corrections here as needed.
      assert numChains == 1;
      numChains = 0;
    }

    init(shape.dimension(), numChains);
  }

  /**
   * Marks a chain as having been seen. This operation is idempotent and may be called multiple
   * times for the same chain safely.
   */
  public void markChain(int chain) {
    chainsSeen.set(chain);
  }

  /**
   * Processes a list of cell boundary crossings produced by S2RobustCellClipper and either adds
   * interval crossings (for polygons), point crossings (for polylines), or does nothing (for
   * points).
   *
   * <p>For a 2D shape, if there are no crossings, then this operation does nothing. The user should
   * check whether the shape contains the cell boundary and add it via addCellBoundary() if
   * necessary.
   */
  public void processCrossings(S2Cell cell, List<Crossing> crossings) {
    // If there are no crossings, or all we have are points, do nothing.
    int ncrossing = crossings.size();
    if (ncrossing == 0 || dimension == 0) {
      return;
    }

    // Lookup cell boundary constants once.
    double[] cellUvCoords = new double[4];
    int[] cellIjCoords = new int[4];
    for (int k = 0; k < 4; ++k) {
      Boundary b = Boundary.fromValue(k);
      cellUvCoords[k] = cell.getUVCoordOfBoundary(b);
      cellIjCoords[k] = cell.getIJCoordOfBoundary(b);
    }

    if (dimension == 1) { //
      // Note that we convert from UV to IJ using std::round which always rounds ties away from
      // zero. This is needed to ensure that we round to correct values across face boundaries that
      // might be flipped.
      for (Crossing crossing : crossings) {
        int ij = uvToIjRound(crossing.intercept);
        int face = cell.face();
        int axis = constantBoundaryAxis(cell.face(), crossing.boundary);
        int coord = cellIjCoords[crossing.boundary.value];

        if (crossing.boundary.value < 2) {
          addPoint(face, axis, coord, ij);
        } else {
          delPoint(face, axis, coord, ij);
        }
      }
      return;
    }

    // We have one or more crossings for a polygon, we need to scan around the boundary and make
    // intervals where the shape contains the boundary.
    //
    // Crossings are ordered counter-clockwise only for the purposes of toggling shape-insideness,
    // the actual U or V coordinates of the intercepts aren't necessarily monotonic because we don't
    // toggle edges, we increment the winding number at coordinates, which is computed in the same
    // way on the other side of the cell boundary.
    boolean interior = crossings.get(0).crossingType == CrossingType.INCOMING;
    int i = 0;
    int[] bs = {0}; // Workaround for Java not supporting lambda capture. Just b.
    Boundary boundary;
    for (int b = 0; b < 4; ++b) {
      bs[0] = b;
      boundary = Boundary.fromValue(b);
      // Lookup the axis, constant coordinate, and start and end coordinates of this boundary
      // segment.
      int axis = constantBoundaryAxis(cell.face(), boundary);
      int coord = cellIjCoords[b];

      int bnext = (b + 1) % 4;
      int bprev = (b + 3) % 4;

      double uvbeg = cellUvCoords[bprev];
      double uvend = cellUvCoords[bnext];
      int ijbeg = cellIjCoords[bprev];
      int ijend = cellIjCoords[bnext];

      // Returns true if two intercepts are properly order along the boundary. For boundaries 0 and
      // 1, intercepts increase, and for 2 and 3, they decrease.
      XYPredicate ordered = (x, y) -> bs[0] < 2 ? x < y : x > y;

      // If there's a crossing on this boundary and we start in the interior, then add an interval
      // between the start corner and the crossing. We do this as a separate step so that we can use
      // the ij coordinate of the start point directly (which we know exactly). If we used the uv
      // coordinate of the start point and rounded, it's possible not to end up on the correct ij
      // value for the corner.
      double uvprev;
      if (i < ncrossing && crossings.get(i).boundary == boundary) {
        Crossing crossing = crossings.get(i);
        uvprev = crossing.intercept;
        ++i;

        // In very rare circumstances, it's possible to get a crossing on a boundary, but that has
        // an intercept value before (or after) the corner. This can happen when an edge crosses
        // very close to the corner.
        //
        // Fortunately, either:
        //   1. We're in the interior so we just added the b-1 segment.
        //   2. Or we're not and we're about to add this boundary segment.
        //
        // Either way, the negative crossing doesn't change the portion of the cell boundary in this
        // cell that's contained, so we can ignore it.
        if (interior && ordered.test(uvbeg, crossing.intercept)) {
          double uv = crossing.intercept;

          int ij0 = ijbeg;
          int ij1 = (b < 2) ? uvToIjCeil(uv) : uvToIjFloor(uv);

          if (ij0 != ij1) {
            addInterval(cell.face(), axis, coord, ij0, ij1);
          }
        }

        // Crossing an outgoing edge puts us in the interior.
        interior = (crossing.crossingType == CrossingType.OUTGOING);
      } else {
        // No crossings on this boundary. Just add the entire boundary if needed.
        if (interior) {
          addInterval(cell.face(), axis, coord, ijbeg, ijend);
        }
        continue;
      }

      // Proceed crossing by crossing, added each interval that covers part of the polygon interior,
      // and toggle interior state at each crossing.
      while (i < ncrossing && crossings.get(i).boundary == boundary) {
        Crossing crossing = crossings.get(i);

        if (interior) {
          // Get canonical ordering of UV coordinates where uv0 < uv1.
          double uv0 = min(uvprev, crossing.intercept);
          double uv1 = max(uvprev, crossing.intercept);

          // Convert to IJ coordinates but round _outwards_ from the interior.
          int ij0 = uvToIjFloor(uv0);
          int ij1 = uvToIjCeil(uv1);
          if (b >= 2) {
            int tmp = ij0;
            ij0 = ij1;
            ij1 = tmp;
          }

          // If the coordinates are equal despite rounding in different directions, then we must
          // have landed on an integer ij value to the limits of floating point precision.
          //
          // For an interior interval, this can only happen if:
          //   1. We had a reverse duplicate edge (not allowed)
          //   2. We crossed exactly at a corner.
          //
          // One case isn't allowed and the other doesn't affect the result, so we can ignore it.
          if (ij0 != ij1) {
            addInterval(cell.face(), axis, coord, ij0, ij1);
          }
        }

        // Crossing an outgoing edge puts us in the interior.
        interior = (crossing.crossingType == CrossingType.OUTGOING);
        uvprev = crossing.intercept;
        ++i;
      }

      // Either there's no more crossings or we started a new boundary. If we're still in the
      // interior, then we need to add an interval between the last crossing and the endpoint.
      //
      // We do this as a separate step so that we can use the ij coordinate of the start point
      // directly (which we know exactly). If we used the uv coordinate of the start point and
      // rounded, it's possible not to end up on the correct ij value for the corner.
      if (interior) {
        // See comment in the start point section above.
        if (ordered.test(uvprev, uvend)) {
          int ij0 = (b < 2) ? uvToIjFloor(uvprev) : uvToIjCeil(uvprev);
          int ij1 = ijend;

          // Safe to ignore, see comment above.
          if (ij0 != ij1) {
            addInterval(cell.face(), axis, coord, ij0, ij1);
          }
        }
      }
    }
    Preconditions.checkState(i == ncrossing);
  }

  /**
   * Adds the boundaries of the given cell to the tracker.
   *
   * <p>When a cell is contained by a shape but has no crossings, we still need to account for those
   * contained cell boundaries and ensure that we see them on both sides of the boundary.
   *
   * <p>This function only applies to 2D shapes. For points and polylines, this function is a noop.
   */
  public void addCellBoundary(S2Cell cell) {
    int face = cell.face();
    int bottomAxis = constantBoundaryAxis(face, BOTTOM_EDGE);
    int rightAxis = constantBoundaryAxis(face, RIGHT_EDGE);
    int topAxis = constantBoundaryAxis(face, TOP_EDGE);
    int leftAxis = constantBoundaryAxis(face, LEFT_EDGE);
    int bottomEdgeIJ = cell.getIJCoordOfBoundary(BOTTOM_EDGE);
    int rightEdgeIJ = cell.getIJCoordOfBoundary(RIGHT_EDGE);
    int topEdgeIJ = cell.getIJCoordOfBoundary(TOP_EDGE);
    int leftEdgeIJ = cell.getIJCoordOfBoundary(LEFT_EDGE);

    addInterval(face, bottomAxis, bottomEdgeIJ, leftEdgeIJ, rightEdgeIJ);
    addInterval(face, rightAxis, rightEdgeIJ, bottomEdgeIJ, topEdgeIJ);
    addInterval(face, topAxis, topEdgeIJ, rightEdgeIJ, leftEdgeIJ);
    addInterval(face, leftAxis, leftEdgeIJ, topEdgeIJ, bottomEdgeIJ);
  }

  /**
   * Returns true if the shape is finished. This means that all the chains in the shape have been
   * seen at least once and there are no outstanding interior pieces.
   */
  public boolean finished() {
    return chainsSeen.full() && points0.isEmpty() && points1.isEmpty();
  }

  /**
   * Creates a 64 bit integer (which is treated as unsigned) by appending face, i and j coordinates.
   */
  private long ijKey(int face, int i, int j) {
    assert face < 6;
    assert i <= 1 << 30;
    assert j <= 1 << 30;

    long ans = 0;
    ans |= ((long) face) << 60;
    ans |= ((long) i) << 30;
    ans |= (long) j;
    return ans;
  }

  /**
   * Adds the endpoints of an interval to the tracker. The sign of ij0 and ij1 are always the same,
   * this function always adds the point corresponding to ij0 and subtracts the point corresponding
   * to ij1.
   *
   * <p>Any given interval must be canceled out by an interval of the opposite sign, i.e. call this
   * function with the endpoints reversed: ij1, ij0.
   *
   * <p>This ordering appears naturally when working with cell boundary crossings. The order would
   * be ij0,ij1 on a positive boundary and ij1,ij0 on a negative boundary.
   *
   * <p>REQUIRES: Shape dimension is 2, and ij0 and ij1 aren't equal.
   */
  @VisibleForTesting
  void addInterval(int face, int axis, int ijcoord, int ij0, int ij1) {
    assert ij0 != ij1;
    assert dimension == 2;

    // If we're at the maximum coordinate value, wrap back to coordinate 0 on the
    // adjacent face. We always switch axes when we're transitioning across the
    // face boundary, but only when we go across the maximum J value do we negate
    // the value we were given.
    if (ijcoord == S2Projections.LIMIT_IJ) {
      ijcoord = 0;
      face = adjacentFace(face, axis);
      if (axis == 1) {
        ij0 = S2Projections.LIMIT_IJ - ij0;
        ij1 = S2Projections.LIMIT_IJ - ij1;
      }
      axis = 1 - axis;
    }

    // Increment the value at the ijKey for ij0, removing it from the map it then becomes zero.
    points(axis).compute(ijKey(face, ijcoord, ij0), incrementFn);
    // And decrement the value at the ijKey for ij1, removing it from the map if it becomes zero.
    points(axis).compute(ijKey(face, ijcoord, ij1), decrementFn);
  }

  /** Adds a point to the tracker. */
  @VisibleForTesting
  void addPoint(int face, int axis, int ijcoord, int ij) {
    // If we're at the maximum coordinate value, wrap back to coordinate 0 on the adjacent face. We
    // always switch axes when we're transitioning across the face boundary, but only when we go
    // across the maximum J value do we negate the value we were given.

    boolean flip = false;
    if (ijcoord == S2Projections.LIMIT_IJ) {
      ijcoord = 0;
      face = adjacentFace(face, axis);
      if (axis == 1) {
        flip = true;
        ij = S2Projections.LIMIT_IJ - ij;
      }
      axis = 1 - axis;
    }

    if (flip) {
      delPoint(face, axis, ijcoord, ij);
      return;
    }

    // Increment the value at the ijKey, removing it from the map it then becomes zero.
    points(axis).compute(ijKey(face, ijcoord, ij), incrementFn);
  }

  /** Removes a point from the tracker. */
  @VisibleForTesting
  void delPoint(int face, int axis, int ijcoord, int ij) {
    // If we're at the maximum coordinate value, wrap back to coordinate 0 on the adjacent face. We
    // always switch axes when we're transitioning across the face boundary, but only when we go
    // across the maximum J value do we negate the value we were given.
    boolean flip = false;
    if (ijcoord == S2Projections.LIMIT_IJ) {
      ijcoord = 0;
      face = adjacentFace(face, axis);
      if (axis == 1) {
        flip = true;
        ij = S2Projections.LIMIT_IJ - ij;
      }
      axis = 1 - axis;
    }

    if (flip) {
      addPoint(face, axis, ijcoord, ij);
      return;
    }

    // Decrement the value at the ijKey, removing it from the map if it becomes zero.
    points(axis).compute(ijKey(face, ijcoord, ij), decrementFn);
  }

  /**
   * Function for map.compute() that increments a map entry value. Creates the entry with starting
   * value 1 if the key is not present. Removes the entry if the value becomes zero.
   */
  private static final BiFunction<Long, Integer, Integer> incrementFn =
      (Long key, Integer value) -> {
        if (value == null) {
          return 1;
        }
        if (value == -1) {
          return null;
        }
        return value + 1;
      };

  /**
   * Function for map.compute() that decrements a map entry value. Creates the entry with starting
   * value -1 if the key is not present. Removes the entry if the value becomes zero.
   */
  private static final BiFunction<Long, Integer, Integer> decrementFn =
      new BiFunction<Long, Integer, Integer>() {
        @Override
        public @Nullable Integer apply(Long key, Integer value) {
          if (value == null) {
            return -1;
          }
          if (value == 1) {
            return null;
          }
          return value - 1;
        }
      };

  /** Converts from UV to IJ coordinates but always rounds down. */
  private static int uvToIjFloor(double uv) {
    return (int) floor(S2Projections.LIMIT_IJ * S2Projections.uvToST(uv));
  }

  /** Converts from UV to IJ coordinates but always rounds up. */
  private static int uvToIjCeil(double uv) {
    return (int) ceil(S2Projections.LIMIT_IJ * S2Projections.uvToST(uv));
  }

  /** Converts from UV to IJ coordinates but rounds ties away from zero. */
  private static int uvToIjRound(double uv) {
    return (int) round(S2Projections.LIMIT_IJ * S2Projections.uvToST(uv));
  }

  /** Returns the axis along which a given cell boundary is constant. */
  private static int constantBoundaryAxis(int face, Boundary boundary) {
    int axis = 0;
    if ((boundary == BOTTOM_EDGE) || (boundary == TOP_EDGE)) {
      axis = 1;
    }

    // Odd faces have axes flipped.
    return (face % 2 == 0) ? axis : 1 - axis;
  }

  /** Returns the face adjacent to the given face across the given axis. */
  private static int adjacentFace(int face, int axis) {
    return (face + (axis == 1 ? 2 : 1)) % 6;
  }

  private interface XYPredicate {
    boolean test(double x, double y);
  }
}
