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

import com.google.errorprone.annotations.CanIgnoreReturnValue;

/**
 * A clipper of shape edges to rectangular regions in UV space. Layered over R2EdgeClipper,
 * providing UV specific semantics and error bounds. UVEdgeClipper isn't cell specific. It can clip
 * to any rectangular region in UV space, set via calling {@link #init(int, R2Rect)} or the
 * equivalent constructor before clipping edges. Since clipping to {@link S2Cell} boundaries is the
 * most common use of this class, we provide a convenience method {@link #init(S2Cell)} and
 * equivalent constructor which set the face and clip rect from the cell.
 *
 * <p>UVEdgeClipper does not clip exactly or use exact tests to determine boundary crossings, so it
 * is possible for points very close to the boundary to falsely test as crossing. We include an
 * error bound when reporting UV edges so that care can be taken in those cases.
 *
 * <p>Example usage (given some S2Cell 'cell' and an {@code S2Point[2]} edge:
 *
 * {@snippet :
 * UVEdgeClipper clipper = new UVEdgeClipper();
 *
 * clipper.init(cell);
 * if (clipper.clipEdge(edge[0], edge[1])) {
 *   if (clipper.outcode(0) == LEFT) {
 *     System.out.println("The first vertex landed on the left boundary!");
 *   }
 * }
 * }
 *
 * <p>The clipEdge method returns true if the edge intersected the clip face and clip region, false
 * otherwise. When the edge misses, the {@link #missedFace()} method can be used to distinguish
 * whether the edge missed the clip face entirely, or hit the face but missed the clip region.
 *
 * <p>The clipper also provides the following information after clipping an edge:
 *
 * <pre>
 *
 *   faceUvEdge(): The UV coordinates of the edge after being clipped to the clip face, but before
 *     being clipped to the clip region.
 *
 *   clippedUvEdge(): the UV coordinates of the edge after being clipped to the clip face and then
 *     the clip region.
 *
 *    outcode(): A simple (single bit set) outcode for each vertex indicating which edge of the
 *      region it was clipped to (i.e. which boundary segment it landed on). If the vertex landed
 *      inside the clip region to begin with, this will be Outcode.INSIDE. If the edge doesn't
 *      intersect the clip region (including missing the face), this will be Outcode.OUTSIDE.
 *
 *   uvError(): The maximum absolute error incurred converting to UV coordinates and possibly
 *     clipping the edge to the clip face. If no face clipping is performed, this is just the
 *     XYZ => UV conversion error. This is the maximum absolute error in each of the U and V
 *     components of the vertices of faceUvEdge().
 *
 *   clipError(): The maximum absolute error in a clipped vertex. This necessarily includes
 *     uvError() since we always have to convert the input edge to UV coordinates and also includes
 *     any error incurred when interpolating to find the intersection point on a region boundary.
 *     This is the maximum absolute error in each of the U and V components of the vertices of
 *     clippedUvEdge().
 *
 * The information available depends on the result status:
 *
 *   hit?   missedFace?  faceUvEdge  clippedUvEdge  outcode  uvError  clipError
 *   ───────────────────────────────────────────────────────────────────────────
 *   false   true         NO          NO            YES      YES      YES
 *   false  false        YES          NO            YES      YES      YES
 *    true  -----        YES         YES            YES      YES      YES
 *
 * </pre>
 *
 * <p>We treat the clip region as a closed set (i.e. points on the boundary test as contained).
 * Given an Edge E, and a clip region R, we ensure the following properties hold:
 *
 * <ol>
 *   <li>If a vertex of E is contained in R, then that vertex is left unmodified.
 *   <li>If E doesn't intersect R at all then it's left unmodified.
 *   <li>Otherwise E is clipped so that each vertex that falls outside of R is moved onto one of the
 *       four boundary segments _exactly_. I.e. at least one component (U or V) of the clipped
 *       vertex will _equal_ a boundary coordinate. Intersection points near corners may fall
 *       outside the actual clip region due to numerical error.
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
 * <p>We take the convention that left and right regions partition the U axis and the top and bottom
 * regions partition the V axis:
 *
 * {@snippet :
 * Left       Right
 *
 *   TL │ T │ TR     Top
 *  ────┼───┼────  ^
 *    L │ I │ R    | v
 *  ────┼───┼────  |
 *   BL │ B │ BR    Bottom
 *
 *      ---->
 *        u
 * }
 *
 * <p>We then clip line segments (edges) based on which region its endpoints landed in. In the worst
 * case, we have to perform 4 clipping operations to fully clip the endpoints, though we take steps
 * to amortize this cost.
 *
 * <p>Cohen-Sutherland is ~3-4x faster than {@link S2EdgeUtil#clipEdge}, for several reasons:
 *
 * <ol>
 *   <li>It can quickly accept and reject edges that don't cross the clipping region
 *   <li>It can re-use information from a previous connected edge to speed up clipping
 *   <li>The clipped results are easier to use since they will exactly intersect the clipping
 *       boundary
 * </ol>
 */
@SuppressWarnings("Assertion")
public class UVEdgeClipper {
  final R2EdgeClipper r2Clipper = new R2EdgeClipper();

  /** The face being clipped to. */
  private int clipFace = 0;

  /** Face of the last vertex. */
  private byte lastFace = 0;

  private final R2Edge faceUvEdge = new R2Edge();
  private double uvError;
  private boolean missedFace;

  /** Constructor that does not set a face or region. */
  public UVEdgeClipper() {}

  /** Constructor that sets the face and clip region to the given values. */
  public UVEdgeClipper(int face, R2Rect region) {
    this();
    init(face, region);
  }

  /** Constructor that sets the face and clip region to the boundary of the given cell. */
  public UVEdgeClipper(S2Cell cell) {
    this();
    init(cell);
  }

  /** Initialize the clipper to clip to the given UV region on the given face. */
  public void init(int face, R2Rect region) {
    this.clipFace = face;
    r2Clipper.init(region);
  }

  /** Initialize the clipper to clip to the given S2Cell. */
  public void init(S2Cell cell) {
    init(cell.face(), cell.getBoundUV());
  }

  /** Returns the face being clipped to. */
  public int clipFace() {
    return clipFace;
  }

  /** Returns the current clipping rectangle. */
  public R2Rect clipRect() {
    return r2Clipper.clipRect();
  }

  /**
   * Clip an edge to the current face and clip region. After clipping, result details can be
   * obtained by calling {@link #missedFace()}, {@link #uvError()}, {@link #clipError()}, {@link
   * #faceUvEdge()}, {@link #clippedUvEdge()}, and {@link #outcode(int)} as needed.
   *
   * <p>If the edge intersects both the face and clip region, returns true, false otherwise. If the
   * edge misses the clip region, then it may have missed either because it missed the face entirely
   * or hit the face but missed the clip region. The {@link #missedFace()} method can be used to
   * distinguish between the two cases.
   *
   * <p>If connected is true, then the clipper assumes the current edge follows the previous edge
   * passed to clipEdge() (i.e. current.v0 == previous.v1) and will reuse previous computations.
   */
  @CanIgnoreReturnValue
  public boolean clipEdge(S2Point v0, S2Point v1, boolean connected) {
    // TODO(user): After refactoring, merge this with the clipping code in S2ShapeIndex.

    // Check that we didn't get the origin somehow, since we use it as a sentinel value, and that
    // none of the points are larger than we allowed for in our error bounds.
    assert !v0.equalsPoint(S2Point.ZERO) && !v1.equalsPoint(S2Point.ZERO);
    assert v0.largestAbsComponent() <= 2;
    assert v1.largestAbsComponent() <= 2;

    // We'll convert to UV directly into the faceUvEdge field.
    byte face0;
    byte face1;

    boolean needFaceClip = false;
    if (connected) {
      face0 = lastFace;
      face1 = (byte) S2Projections.xyzToFace(v1);
      if (face0 != face1 || face0 != clipFace) {
        // Vertices are on different faces, or the same face, but not this face.
        needFaceClip = true;
      } else {
        // Vertices are both on the clip face, just convert v1.
        faceUvEdge.v0.set(faceUvEdge.v1);
        S2Projections.validFaceXyzToUv(face1, v1, faceUvEdge.v1);
      }
    } else {
      // We can't re-use the values from the last vertex for whatever reason, so just convert and
      // clip both vertices.
      face0 = (byte) S2Projections.xyzToFace(v0);
      face1 = (byte) S2Projections.xyzToFace(v1);
      if (face0 != face1 || face0 != clipFace) {
        // Vertices are on different faces, or the same face, but not this face.
        needFaceClip = true;
      } else {
        // Both vertices are on the clip face, convert both to UV directly. Avoid having a
        // conditional between the conversion calls as the compiler can fuse them with SIMD giving
        // us one for ~free.
        S2Projections.validFaceXyzToUv(face0, v0, faceUvEdge.v0);
        S2Projections.validFaceXyzToUv(face0, v1, faceUvEdge.v1);
      }
    }
    lastFace = face1;

    // The vertices are on different faces, clip the edge to the clip face. It's possible that we
    // think we missed the face, but the edge is close enough to the boundary that it will hit the
    // padded face, so we need to clip in that case too to verify we truly missed.
    missedFace = false;
    if (needFaceClip) {
      if (!S2EdgeUtil.clipToPaddedFace(
          v0, v1, clipFace, S2EdgeUtil.FACE_CLIP_ERROR_UV_COORD, faceUvEdge.v0, faceUvEdge.v1)) {
        missedFace = true;
        return false;
      }
    }

    // Set error bounds for converting vertices to the face.
    uvError = S2Projections.MAX_XYZ_TO_UV_ERROR;
    if (needFaceClip) {
      uvError += S2EdgeUtil.FACE_CLIP_ERROR_UV_COORD;
    }

    // Clip the edge as though it's a regular R2 edge.
    return r2Clipper.clipEdge(faceUvEdge, connected);
  }

  /**
   * Like {@link #clipEdge(S2Point, S2Point, boolean)} but with a default false value for
   * "connected". See the documentation for {@link #clipEdge(S2Point, S2Point, boolean)} for more
   * details.
   */
  @CanIgnoreReturnValue
  public boolean clipEdge(S2Point v0, S2Point v1) {
    return clipEdge(v0, v1, false);
  }

  /**
   * When the edge misses the clip region, this indicates whether it was because the edge missed the
   * face entirely or not.
   */
  public boolean missedFace() {
    return missedFace;
  }

  /**
   * Returns the maximum absolute error incurred from converting the edge to UV coordinates and
   * clipping to the clip face. This is the error bound in each UV coordinate of the vertices of
   * faceUvEdge().
   *
   * <p>When both vertices are on the clip face, we only have to convert from XYZ to UV coordinates,
   * which is very accurate (DBL_EPSILON/4 absolute error, see {link S2#faceXYZToUV}.
   *
   * <p>If we have to clip the edge to the clip face too, more error is incurred,
   * (FACE_CLIP_ERROR_UV_COORD), which is added to the error bound.
   */
  public double uvError() {
    return uvError;
  }

  /**
   * Returns the maximum absolute error incurred computing an intersection point on the boundary of
   * the clip region.
   *
   * <p>When an edge crosses the clip region and we have to compute new vertices for it, this it the
   * error bound for each of the UV coordinates of the resulting vertex.
   *
   * <p>This necessarily includes uv_error() and any additional error incurred from linearly
   * interpolating the vertex onto the boundary. We may end up clipping an edge twice and we can
   * have uv conversion error on both vertices, so we multiply by two.
   */
  // TODO(eengle): Verify this is a safe upper bound in the tests
  public double clipError() {
    return 2 * (uvError() + R2EdgeClipper.MAX_UNIT_CLIP_ERROR);
  }

  /**
   * If the edge intersected the clip face, returns the full UV edge after it has been clipped to
   * the face but before it has been clipped to the clip region.
   */
  public R2Edge faceUvEdge() {
    return faceUvEdge;
  }

  /** If the edge intersected the clip region, returns the clipped edge. */
  public R2Edge clippedUvEdge() {
    return r2Clipper.clippedEdge;
  }

  /**
   * Returns the final outcode computed for a vertex. This value indicates which edge of the clip
   * region the vertex landed on. If the vertex was inside the clip region to begin with, it's not
   * modified so this returns {@code Outcode#INSIDE}.
   *
   * <p>If missedFace() is true or clipEdge() returned false (indicating a miss), this will return
   * {@code Outcode.OUTSIDE}.
   *
   * <p>If a vertex manages to land exactly on a corner of the region (touching two edges), this
   * will return one of the two edges. Which one is unspecified but will always be the same during a
   * particular run of the program.
   *
   * <p>REQUIRES: vertex is 0 or 1.
   */
  public byte outcode(int vertex) {
    if (vertex == 0) {
      return r2Clipper.outcode0;
    } else if (vertex == 1) {
      return r2Clipper.outcode1;
    }
    throw new IllegalArgumentException("vertex must be 0 or 1");
  }
}
