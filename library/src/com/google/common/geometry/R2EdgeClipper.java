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

/**
 * Class to clip edges to rectangular regions in 2D space.
 *
 * <p>R2EdgeClipper does not clip exactly or use exact tests to determine boundary crossings. It's
 * possible for points very close to the boundary to falsely test as crossing.
 *
 * <p>Example usage (given some R2Rect region):
 *
 * {@snippet :
 *
 *   R2EdgeClipper clipper = new R2EdgeClipper(region);
 *   if (clipper.clipEdge(edge, false)) {
 *     if (clipper.outcode0 == LEFT) {
 *        System.out.printf("The first vertex landed on the left boundary!\n");
 *     }
 *   }
 * }
 *
 * <p>The clipped edge can be obtained via {@link #clippedEdge} (assuming it hit the clip region),
 * along with a {@link #outcode0} and {@link #outcode1} indicating which boundary of the clip region
 * the two vertices were clipped to.
 *
 * <p>We treat the clip region as a closed set (i.e. points on the boundary test as contained).
 * We're careful to ensure the following properties hold:
 *
 * <p>Given an Edge E, and a clip region R:
 *
 * <ol>
 *   <li>If a vertex of E is contained in R, then that vertex is left unmodified.
 *   <li>If E doesn't intersect R at all then it's left unmodified.
 *   <li>Otherwise E is clipped so that each vertex that falls outside of R is moved onto one of the
 *       our boundary segments _exactly_. I.e. at least one component (X or Y) of the clipped vertex
 *       will _equal_ a boundary coordinate. Intersection points near corners may fall outside the
 *       actual clip region due to numerical error.
 *   <li>If an edge crosses a boundary segment S that is the same between two different clip regions
 *       R and Q, then clipping the edge against R or against Q yields the _exact_ same result,
 *       which is the intersection point on S.
 * </ol>
 *
 * <p>We use the Cohen-Sutherland algorithm which classifies each endpoint of an edge by which
 * region it falls into relative to the clip region: top, bottom, left and right. Points that fall
 * into one of the corner regions are tagged as being in both of the regions that overlap it (e.g.
 * top-left).
 *
 * <p>We take the convention that left and right regions partition the X axis and the top and bottom
 * regions partition the Y axis:
 *
 * <pre>
 *
 *     Left       Right
 *
 *       TL │ T │ TR     Top
 *      ────┼───┼────  ^
 *        L │ I │ R    | Y
 *      ────┼───┼────  |
 *       BL │ B │ BR    Bottom
 *
 *          ---->
 *            X
 * </pre>
 *
 * <p>We then clip line segments (edges) based on which region its endpoints landed in. In the worst
 * case, we have to perform 4 clipping operations to fully clip the endpoints, though we take steps
 * to amortize this cost.
 *
 * <p>Cohen-Sutherland has a few advantages for clipping edges for our use case:
 *
 * <pre>
 *
 *   1. It can quickly accept and reject line segments that are entirely within (or outside) the
 *      clip region, which is the majority of our edges.
 *
 *   2. It can re-use information from a previous edge to speed up clipping.
 *
 *   3. The clipped results are guaranteed to exactly land on a clipping boundary.
 * </pre>
 */
@SuppressWarnings("Assertion")
public class R2EdgeClipper {
  /**
   * When the magnitudes of the clip region coordinates and the clip points are less than or equal
   * to 1, this is a bound on the absolute error in each coordinate of the clipped edge (i.e. each
   * of the x() and y() coordinates).
   */
  public static final double MAX_UNIT_CLIP_ERROR = 2 * S2EdgeUtil.EDGE_CLIP_ERROR_UV_COORD;

  /*
   * Constants for allowed outcode values which can be operated on with bitwise AND and OR. These
   * are designed to be consistent with {@link S2Cell.Boundary} naming.
   */

  /** An outcode indicating that a vertex is inside the clip region. */
  public static final byte INSIDE = 0x00;

  /** An outcode indicating that a vertex is on the bottom boundary of the clip region. */
  public static final byte BOTTOM = 0x01;

  /** An outcode indicating that a vertex is on the right boundary of the clip region. */
  public static final byte RIGHT = 0x02;

  /** An outcode indicating that a vertex is on the top boundary of the clip region. */
  public static final byte TOP = 0x04;

  /** An outcode indicating that a vertex is on the left boundary of the clip region. */
  public static final byte LEFT = 0x08;

  /** An outcode indicating that a vertex is outside the clip region. */
  public static final byte OUTSIDE = (byte) 0xFF;

  /** The minimum x boundary of the region being clipped to. */
  private double xMin;

  /** The maximum x boundary of the region being clipped to. */
  private double xMax;

  /** The minimum y boundary of the region being clipped to. */
  private double yMin;

  /** The maximum y boundary of the region being clipped to. */
  private double yMax;

  /** The edge being clipped. */
  private R2Edge edge;

  /** The outcode of the last vertex clipped. */
  private byte lastOutcode;

  /** The clipped edge. R2Edge is mutable, and is reused for each clipping operation. */
  public final R2Edge clippedEdge = new R2Edge();

  /**
   * If the last edge was not WRONG_FACE, then this is the outcode for the first vertex of the
   * clipped edge, indicating which boundary it landed on after clipping.
   */
  public byte outcode0;

  /**
   * If the last edge was not WRONG_FACE, then this is the outcode for the second vertex of the
   * clipped edge, indicating which boundary it landed on after clipping.
   */
  public byte outcode1;

  /** The default constructor does not set a clip rectangle. */
  public R2EdgeClipper() {}

  /** Constructor that sets the clip rectangle. */
  public R2EdgeClipper(R2Rect rectangle) {
    init(rectangle);
  }

  /** Sets the clip rectangle to the given {@code rectangle}. */
  public void init(R2Rect rectangle) {
    xMin = rectangle.x().lo();
    xMax = rectangle.x().hi();
    yMin = rectangle.y().lo();
    yMax = rectangle.y().hi();
  }

  /** Returns the current clipping rectangle. */
  public R2Rect clipRect() {
    return new R2Rect(new R1Interval(xMin, xMax), new R1Interval(yMin, yMax));
  }

  // TODO(torrey): Remove the connected parameter, and instead provide an overload of clipEdge that
  // takes the next connected vertex, for Java and C++.

  /**
   * Clips an edge to the current clip rectangle.
   *
   * <p>Returns true when the edge intersected the clip region, false otherwise. If connected is
   * true, then the clipper assumes that v1 of the last edge passed is equal to v0 of the current
   * edge and re-uses any previous calculations it can.
   */
  public boolean clipEdge(R2Edge edge, boolean connected) {
    // Assume we're outside the clip rectangle until we know otherwise.
    outcode0 = OUTSIDE;
    outcode1 = OUTSIDE;

    byte code0 = connected ? lastOutcode : outcode(edge.v0);
    byte code1 = outcode(edge.v1);
    lastOutcode = code1;

    // If both vertices are in the same outside region, the edge between them can't intersect the
    // clip region, so there's nothing to do.
    if ((code0 & code1) != INSIDE) {
      return false;
    }

    clippedEdge.init(edge);
    if (code0 != INSIDE) {
      code0 = clipVertex(clippedEdge.v0, edge, code0);
    }

    if (code1 != INSIDE) {
      code1 = clipVertex(clippedEdge.v1, edge, code1);
    }

    outcode0 = code0;
    outcode1 = code1;

    // If either clipped endpoint is outside, then the line was outside.
    return code0 != OUTSIDE && code1 != OUTSIDE;
  }

  /**
   * Unconditionally clips an edge to the boundary represented by an outcode and returns the result
   * in the given {@code result}.
   *
   * <p>REQUIRES: outcode is a power of two (represents a single region)
   */
  public void clip(R2Edge edge, byte outcode, R2Vector result) {
    this.edge = edge;
    assert outcode > 0;
    assert (outcode & (outcode - 1)) == 0;

    switch (outcode) {
      case BOTTOM:
        clipBottom(result);
        break;
      case RIGHT:
        clipRight(result);
        break;
      case TOP:
        clipTop(result);
        break;
      case LEFT:
        clipLeft(result);
        break;
      default:
        throw new IllegalArgumentException("Invalid outcode: " + outcode);
    }
  }

  /**
   * Returns a logical or of the outcodes indicating which clip region boundaries the point falls
   * outside of, or {@link Outcode#INSIDE} if the point is inside the clip region.
   */
  private byte outcode(R2Vector uv) {
    byte code = 0;
    if (uv.x() < xMin) {
      code |= LEFT;
    } else if (uv.x() > xMax) {
      code |= RIGHT;
    }
    if (uv.y() < yMin) {
      code |= BOTTOM;
    } else if (uv.y() > yMax) {
      code |= TOP;
    }
    return code;
  }

  /** Finds the intersection point between an edge and the bottom clip edge. */
  private void clipBottom(R2Vector intersection) {
    intersection.set(interpolateX(yMin), yMin);
  }

  /** Finds the intersection point between an edge and the right clip edge. */
  private void clipRight(R2Vector intersection) {
    intersection.set(xMax, interpolateY(xMax));
  }

  /** Finds the intersection point between an edge and the top clip edge. */
  private void clipTop(R2Vector intersection) {
    intersection.set(interpolateX(yMax), yMax);
  }

  /** Finds the intersection point between an edge and the left clip edge. */
  private void clipLeft(R2Vector intersection) {
    intersection.set(xMin, interpolateY(xMin));
  }

  /** Interpolates Y between two endpoints of a XY edge based on a X coordinate. */
  private double interpolateY(double x) {
    R2Vector v0 = edge.v0;
    R2Vector v1 = edge.v1;
    return S2EdgeUtil.interpolateDouble(x, v0.x(), v1.x(), v0.y(), v1.y());
  }

  /** Interpolates X between two endpoints of a XY edge based on a Y coordinate. */
  private double interpolateX(double y) {
    R2Vector v0 = edge.v0;
    R2Vector v1 = edge.v1;
    return S2EdgeUtil.interpolateDouble(y, v0.y(), v1.y(), v0.x(), v1.x());
  }

  /**
   * Clips a vertex of an edge based on its outcode and returns an outcode representing the boundary
   * that the vertex was clipped to. If the edge fell outside the clipping region, then OUTSIDE is
   * returned.
   */
  private byte clipVertex(R2Vector v0, R2Edge edge, byte code) {
    this.edge = edge;
    assert code != INSIDE;
    assert code != OUTSIDE;

    // Simple regions just refer to a single side of a boundary, so they only have to be clipped
    // once. Power of two outcodes are simple.
    if ((code & (code - 1)) == 0) {
      clip(edge, code, v0);

      // Return OUTSIDE if the vertex is still outside the clip region, otherwise return
      // the boundary we clipped to (which is just the outcode passed in).
      if (outcode(v0) != INSIDE) {
        return OUTSIDE;
      }
      return code;
    }

    // If the vertex is in one of the corners we may have to clip it twice to find the right axis
    // to clip in. We can clip in both at the same time and then afterwards pick the correct one.
    R2Vector va = new R2Vector();
    R2Vector vb = new R2Vector();
    byte outa;
    byte outb;

    switch (code) {
      case TOP | LEFT:
        outa = TOP;
        outb = LEFT;
        clipTop(va);
        clipLeft(vb);
        break;
      case TOP | RIGHT:
        outa = TOP;
        outb = RIGHT;
        clipTop(va);
        clipRight(vb);
        break;
      case BOTTOM | LEFT:
        outa = BOTTOM;
        outb = LEFT;
        clipBottom(va);
        clipLeft(vb);
        break;
      case BOTTOM | RIGHT:
        outa = BOTTOM;
        outb = RIGHT;
        clipBottom(va);
        clipRight(vb);
        break;
      default:
        throw new IllegalStateException("Invalid outcode: " + code);
    }

    if (outcode(va) == INSIDE) {
      v0.set(va);
      return outa;
    }

    if (outcode(vb) == INSIDE) {
      v0.set(vb);
      return outb;
    }

    // Neither clipped point landed inside the clip region, the edge must not actually intersect
    // it.
    return OUTSIDE;
  }
}
