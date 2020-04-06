/*
 * Copyright 2005 Google Inc.
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

import com.google.common.annotations.GwtCompatible;
import com.google.common.geometry.S2.Metric;

/**
 * This class specifies the coordinate systems and transforms used to project points from the sphere
 * to the unit cube to an {@link S2CellId}.
 *
 * <p>In the process of converting a latitude-longitude pair to a 64-bit cell id, the following
 * coordinate systems are used:
 *
 * <ul>
 *   <li>(id): An S2CellId is a 64-bit encoding of a face and a Hilbert curve position on that face.
 *       The Hilbert curve position implicitly encodes both the position of a cell and its
 *       subdivision level (see s2cellid.h).
 *   <li>(face, i, j): Leaf-cell coordinates. "i" and "j" are integers in the range [0,(2**30)-1]
 *       that identify a particular leaf cell on the given face. The (i, j) coordinate system is
 *       right-handed on each face, and the faces are oriented such that Hilbert curves connect
 *       continuously from one face to the next.
 *   <li>(face, s, t): Cell-space coordinates. "s" and "t" are real numbers in the range [0,1] that
 *       identify a point on the given face. For example, the point (s, t) = (0.5, 0.5) corresponds
 *       to the center of the top-level face cell. This point is also a vertex of exactly four cells
 *       at each subdivision level greater than zero.
 *   <li>(face, si, ti): Discrete cell-space coordinates. These are obtained by multiplying "s" and
 *       "t" by 2**31 and rounding to the nearest unsigned integer. Discrete coordinates lie in the
 *       range [0,2**31]. This coordinate system can represent the edge and center positions of all
 *       cells with no loss of precision (including non-leaf cells). In binary, each coordinate of a
 *       level-k cell center ends with a 1 followed by (30 - k) 0s. The coordinates of its edges end
 *       with (at least) (31 - k) 0s.
 *   <li>(face, u, v): Cube-space coordinates. To make the cells at each level more uniform in size
 *       after they are projected onto the sphere, we apply a nonlinear transformation of the form
 *       u=f(s), v=f(t). The (u, v) coordinates after this transformation give the actual
 *       coordinates on the cube face (modulo some 90 degree rotations) before it is projected onto
 *       the unit sphere.
 *   <li>(face, u, v, w): Per-face coordinate frame. This is an extension of the (face, u, v)
 *       cube-space coordinates that adds a third axis "w" in the direction of the face normal. It
 *       is always a right-handed 3D coordinate system. Cube-space coordinates can be converted to
 *       this frame by setting w=1, while (u,v,w) coordinates can be projected onto the cube face by
 *       dividing by w, i.e. (face, u/w, v/w).
 *   <li>(x, y, z): Direction vector (S2Point). Direction vectors are not necessarily unit length,
 *       and are often chosen to be points on the biunit cube [-1,+1]x[-1,+1]x[-1,+1]. They can be
 *       normalized to obtain the corresponding point on the unit sphere.
 *   <li>(lat, lng): Latitude and longitude (S2LatLng). Latitudes must be between -90 and 90 degrees
 *       inclusive, and longitudes must be between -180 and 180 degrees inclusive.
 * </ul>
 *
 * <p>Note that the (i, j), (s, t), (si, ti), and (u, v) coordinate systems are right-handed on all
 * six faces.
 *
 * <p>We have implemented three different projections from cell-space (s,t) to cube-space (u,v):
 * {@link S2Projections#S2_LINEAR_PROJECTION}, {@link S2Projections#S2_TAN_PROJECTION}, and {@link
 * S2Projections#S2_QUADRATIC_PROJECTION}. The default is in {@link S2Projections#PROJ}, and uses
 * the quadratic projection since it has the best overall behavior.
 *
 * <p>Here is a table comparing the cell uniformity using each projection. "Area Ratio" is the
 * maximum ratio over all subdivision levels of the largest cell area to the smallest cell area at
 * that level, "Edge Ratio" is the maximum ratio of the longest edge of any cell to the shortest
 * edge of any cell at the same level, and "Diag Ratio" is the ratio of the longest diagonal of any
 * cell to the shortest diagonal of any cell at the same level. "ToPoint" and "FromPoint" are the
 * times in microseconds required to convert cell IDs to and from points (unit vectors)
 * respectively. "ToPointRaw" is the time to convert to a non-unit-length vector, which is all that
 * is needed for some purposes.
 *
 * <table>
 * <tr>
 * <th>Projection</th>
 * <th>Area Ratio</th>
 * <th>Edge Ratio</th>
 * <th>Diag Ratio</th>
 * <th>ToPointRaw (microseconds)</th>
 * <th>ToPoint (microseconds)</th>
 * <th>FromPoint (microseconds)</th>
 * </tr>
 * <tr>
 * <td>Linear</td>
 * <td>5.200</td>
 * <td>2.117</td>
 * <td>2.959</td>
 * <td>0.020</td>
 * <td>0.087</td>
 * <td>0.085</td>
 * </tr>
 * <tr>
 * <td>Tangent</td>
 * <td>1.414</td>
 * <td>1.414</td>
 * <td>1.704</td>
 * <td>0.237</td>
 * <td>0.299</td>
 * <td>0.258</td>
 * </tr>
 * <tr>
 * <td>Quadratic</td>
 * <td>2.082</td>
 * <td>1.802</td>
 * <td>1.932</td>
 * <td>0.033</td>
 * <td>0.096</td>
 * <td>0.108</td>
 * </tr>
 * </table>
 *
 * <p>The worst-case cell aspect ratios are about the same with all three projections. The maximum
 * ratio of the longest edge to the shortest edge within the same cell is about 1.4 and the maximum
 * ratio of the diagonals within the same cell is about 1.7.
 *
 * <p>This data was produced using {@code S2CellTest} and {@code S2CellIdTest}.
 *
 * @author eengle@google.com (Eric Engle) ported from util/geometry
 */
@GwtCompatible
public strictfp enum S2Projections {
  // All of the values below were obtained by a combination of hand analysis and
  // Mathematica. In general, S2_TAN_PROJECTION produces the most uniform
  // shapes and sizes of cells, S2_LINEAR_PROJECTION is considerably worse, and
  // S2_QUADRATIC_PROJECTION is somewhere in between (but generally closer to
  // the tangent projection than the linear one).

  /**
   * This is the fastest transformation, but also produces the least uniform cell sizes. Cell areas
   * vary by a factor of about 5.2, with the largest cells at the center of each face and the
   * smallest cells in the corners.
   */
  S2_LINEAR_PROJECTION(
      4 / (3 * Math.sqrt(3)), // minArea 0.770
      4, // maxArea 4.000
      1.0, // minAngleSpan 1.000
      2, // maxAngleSpan 2.000
      Math.sqrt(2. / 3), // minWidth 0.816
      1.411459345844456965, // avgWidth 1.411
      2 * Math.sqrt(2) / 3, // minEdge 0.943
      1.440034192955603643, // avgEdge 1.440
      2 * Math.sqrt(2) / 3, // minDiag 0.943
      2 * Math.sqrt(2), // maxDiag 2.828
      2.031817866418812674, // avgDiag 2.032
      Math.sqrt(2)) { // maxEdgeAspect 1.414
    @Override
    public double stToUV(double s) {
      return 2 * s - 1;
    }

    @Override
    public double uvToST(double u) {
      return 0.5 * (u + 1);
    }
  },
  /**
   * Transforming the coordinates via atan() makes the cell sizes more uniform. The areas vary by a
   * maximum ratio of 1.4 as opposed to a maximum ratio of 5.2. However, each call to atan() is
   * about as expensive as all of the other calculations combined when converting from points to
   * cell IDs, i.e. it reduces performance by a factor of 3.
   */
  S2_TAN_PROJECTION(
      (S2.M_PI * S2.M_PI) / (4 * Math.sqrt(2)), // minArea 1.745
      S2.M_PI * S2.M_PI / 4, // maxArea 2.467
      S2.M_PI / 2, // minAngleSpan 1.571
      S2.M_PI / 2, // maxAngleSpan 1.571
      S2.M_PI / (2 * Math.sqrt(2)), // minWidth 1.111
      1.437318638925160885, // avgWidth 1.437
      S2.M_PI / (2 * Math.sqrt(2)), // minEdge 1.111
      1.461667032546739266, // avgEdge 1.462
      S2.M_PI * Math.sqrt(2) / 3, // minDiag 1.481
      S2.M_PI * Math.sqrt(2. / 3), // maxDiag 2.565
      2.063623197195635753, // avgDiag 2.064
      Math.sqrt(2)) { // maxEdgeAspect 1.414
    @Override
    public double stToUV(double s) {
      // Unfortunately, tan(M_PI_4) is slightly less than 1.0. This isn't due to a flaw in the
      // implementation of tan(), it's because the derivative of tan(x) at x=pi/4 is 2, and it
      // happens that the two adjacent floating point numbers on either side of the infinite-
      // precision value of pi/4 have tangents that are slightly below and slightly above 1.0 when
      // rounded to the nearest double-precision result.
      s = Math.tan(S2.M_PI_2 * s - S2.M_PI_4);
      return s + (1.0 / (1L << 53)) * s;
    }

    @Override
    public double uvToST(double u) {
      return (2 * S2.M_1_PI) * (Math.atan(u) + S2.M_PI_4);
    }
  },
  /**
   * This is an approximation of the tangent projection that is much faster and produces cells that
   * are almost as uniform in size. It is about 3 times faster than the tangent projection for
   * converting cell IDs to points or vice versa. Cell areas vary by a maximum ratio of about 2.1.
   */
  S2_QUADRATIC_PROJECTION(
      8 * Math.sqrt(2) / 9, // minArea 1.257
      2.635799256963161491, // maxArea 2.636
      4. / 3, // minAngleSpan 1.333
      1.704897179199218452, // maxAngleSpan 1.705
      2 * Math.sqrt(2) / 3, // minWidth 0.943
      1.434523672886099389, // avgWidth 1.435
      2 * Math.sqrt(2) / 3, // minEdge 0.943
      1.459213746386106062, // avgEdge 1.459
      8 * Math.sqrt(2) / 9, // minDiag 1.257
      2.438654594434021032, // maxDiag 2.439
      2.060422738998471683, // avgDiag 2.060
      1.442615274452682920) { // maxEdgeAspect 1.443
    @Override
    public double stToUV(double s) {
      if (s >= 0.5) {
        return (1 / 3.) * (4 * s * s - 1);
      } else {
        return (1 / 3.) * (1 - 4 * (1 - s) * (1 - s));
      }
    }

    @Override
    public double uvToST(double u) {
      if (u >= 0) {
        return 0.5 * Math.sqrt(1 + 3 * u);
      } else {
        return 1 - 0.5 * Math.sqrt(1 - 3 * u);
      }
    }
  };

  /**
   * The maximum value of an si- or ti-coordinate. The range of valid (si,ti) values is
   * [0..MAX_SiTi].
   */
  public static final long MAX_SITI = 1L << (S2CellId.MAX_LEVEL + 1);

  /** The U,V,W axes for each face. */
  private static final S2Point[][] FACE_UVW_AXES = {
    {S2Point.Y_POS, S2Point.Z_POS, S2Point.X_POS},
    {S2Point.X_NEG, S2Point.Z_POS, S2Point.Y_POS},
    {S2Point.X_NEG, S2Point.Y_NEG, S2Point.Z_POS},
    {S2Point.Z_NEG, S2Point.Y_NEG, S2Point.X_NEG},
    {S2Point.Z_NEG, S2Point.X_POS, S2Point.Y_NEG},
    {S2Point.Y_POS, S2Point.X_POS, S2Point.Z_NEG}
  };

  /** The precomputed neighbors of each face. See {@link #getUVWFace}. */
  private static final int[][][] FACE_UVW_FACES = {
    {{4, 1}, {5, 2}, {3, 0}},
    {{0, 3}, {5, 2}, {4, 1}},
    {{0, 3}, {1, 4}, {5, 2}},
    {{2, 5}, {1, 4}, {0, 3}},
    {{2, 5}, {3, 0}, {1, 4}},
    {{4, 1}, {3, 0}, {2, 5}}
  };

  /**
   * A transform from 3D cartesian coordinates to the 2D coordinates of a face. For (x, y, z)
   * coordinates within the face, the resulting UV coordinates should each lie in the inclusive
   * range [-1,1], with the center of the face along that axis at 0.
   */
  public abstract static class UvTransform {
    /** Internal implementations only. */
    private UvTransform() {}

    /**
     * Returns the 'u' coordinate of the [u, v] point projected onto a cube face from the given [x,
     * y, z] position.
     */
    public abstract double xyzToU(double x, double y, double z);

    /**
     * Returns the 'u' coordinate of the [u, v] point projected onto a cube face from the given [x,
     * y, z] position.
     */
    public final double xyzToU(S2Point p) {
      return xyzToU(p.x, p.y, p.z);
    }

    /**
     * Returns the 'v' coordinate of the [u, v] point projected onto a cube face from the given [x,
     * y, z] position.
     */
    public abstract double xyzToV(double x, double y, double z);

    /**
     * Returns the 'v' coordinate of the [u, v] point projected onto a cube face from the given [x,
     * y, z] position.
     */
    public final double xyzToV(S2Point p) {
      return xyzToV(p.x, p.y, p.z);
    }
  }

  /**
   * The transforms to convert (x, y, z) coordinates to u and v coordinates on a specific face,
   * indexed by face.
   */
  private static final UvTransform[] UV_TRANSFORMS = {
    new UvTransform() {
      @Override
      public double xyzToU(double x, double y, double z) {
        return y / x;
      }

      @Override
      public double xyzToV(double x, double y, double z) {
        return z / x;
      }
    },
    new UvTransform() {
      @Override
      public double xyzToU(double x, double y, double z) {
        return -x / y;
      }

      @Override
      public double xyzToV(double x, double y, double z) {
        return z / y;
      }
    },
    new UvTransform() {
      @Override
      public double xyzToU(double x, double y, double z) {
        return -x / z;
      }

      @Override
      public double xyzToV(double x, double y, double z) {
        return -y / z;
      }
    },
    new UvTransform() {
      @Override
      public double xyzToU(double x, double y, double z) {
        return z / x;
      }

      @Override
      public double xyzToV(double x, double y, double z) {
        return y / x;
      }
    },
    new UvTransform() {
      @Override
      public double xyzToU(double x, double y, double z) {
        return z / y;
      }

      @Override
      public double xyzToV(double x, double y, double z) {
        return -x / y;
      }
    },
    new UvTransform() {
      @Override
      public double xyzToU(double x, double y, double z) {
        return -y / z;
      }

      @Override
      public double xyzToV(double x, double y, double z) {
        return -x / z;
      }
    }
  };

  /**
   * A transform from 2D cartesian coordinates of a face to 3D directional vectors. The resulting
   * vectors are not necessarily of unit length.
   */
  interface XyzTransform {
    /**
     * Returns the 'x' coordinate for the [x, y, z] point on the unit sphere that projects to the
     * given [u, v] point on a cube face.
     */
    public double uvToX(double u, double v);

    /**
     * Returns the 'y' coordinate for the [x, y, z] point on the unit sphere that projects to the
     * given [u, v] point on a cube face.
     */
    public double uvToY(double u, double v);

    /**
     * Returns the 'z' coordinate for the [x, y, z] point on the unit sphere that projects to the
     * given [u, v] point on a cube face.
     */
    public double uvToZ(double u, double v);
  }

  /**
   * The transforms to convert (u, v) coordinates on a specific face to x-, y-, and z- coordinates,
   * indexed by face.
   */
  private static final XyzTransform[] XYZ_TRANSFORMS = {
    new XyzTransform() {
      @Override
      public double uvToX(double u, double v) {
        return 1;
      }

      @Override
      public double uvToY(double u, double v) {
        return u;
      }

      @Override
      public double uvToZ(double u, double v) {
        return v;
      }
    },
    new XyzTransform() {
      @Override
      public double uvToX(double u, double v) {
        return -u;
      }

      @Override
      public double uvToY(double u, double v) {
        return 1;
      }

      @Override
      public double uvToZ(double u, double v) {
        return v;
      }
    },
    new XyzTransform() {
      @Override
      public double uvToX(double u, double v) {
        return -u;
      }

      @Override
      public double uvToY(double u, double v) {
        return -v;
      }

      @Override
      public double uvToZ(double u, double v) {
        return 1;
      }
    },
    new XyzTransform() {
      @Override
      public double uvToX(double u, double v) {
        return -1;
      }

      @Override
      public double uvToY(double u, double v) {
        return -v;
      }

      @Override
      public double uvToZ(double u, double v) {
        return -u;
      }
    },
    new XyzTransform() {
      @Override
      public double uvToX(double u, double v) {
        return v;
      }

      @Override
      public double uvToY(double u, double v) {
        return -1;
      }

      @Override
      public double uvToZ(double u, double v) {
        return -u;
      }
    },
    new XyzTransform() {
      @Override
      public double uvToX(double u, double v) {
        return v;
      }

      @Override
      public double uvToY(double u, double v) {
        return u;
      }

      @Override
      public double uvToZ(double u, double v) {
        return -1;
      }
    }
  };

  /** Minimum area of a cell at level k. */
  public final Metric minArea;

  /** Maximum area of a cell at level k. */
  public final Metric maxArea;

  /** Average area of a cell at level k. */
  public final Metric avgArea;

  /**
   * Minimum angular separation between opposite edges of a cell at level k. Each cell is bounded by
   * four planes passing through its four edges and the center of the sphere. The angle span metrics
   * relate to the angle between each pair of opposite bounding planes, or equivalently, between the
   * planes corresponding to two different s-values or two different t-values.
   */
  public final Metric minAngleSpan;

  /** Maximum angular separation between opposite edges of a cell at level k. */
  public final Metric maxAngleSpan;

  /** Average angular separation between opposite edges of a cell at level k. */
  public final Metric avgAngleSpan;

  /**
   * Minimum perpendicular angular separation between opposite edges of a cell at level k.
   *
   * <p>The width of a geometric figure is defined as the distance between two parallel bounding
   * lines in a given direction. For cells, the minimum width is always attained between two
   * opposite edges, and the maximum width is attained between two opposite vertices. However, for
   * our purposes we redefine the width of a cell as the perpendicular distance between a pair of
   * opposite edges. A cell therefore has two widths, one in each direction. The minimum width
   * according to this definition agrees with the classic geometric one, but the maximum width is
   * different. (The maximum geometric width corresponds to {@link #maxDiag}.)
   *
   * <p>This is useful for bounding the minimum or maximum distance from a point on one edge of a
   * cell to the closest point on the opposite edge. For example, this is useful when "growing"
   * regions by a fixed distance.
   */
  public final Metric minWidth;

  /** Maximum perpendicular angular separation between opposite edges of a cell at level k. */
  public final Metric maxWidth;

  /** Average perpendicular angular separation between opposite edges of a cell at level k. */
  public final Metric avgWidth;

  /**
   * Minimum angular length of any cell edge at level k. The edge length metrics can also be used to
   * bound the minimum, maximum, or average distance from the center of one cell to the center of
   * one of its edge neighbors. In particular, it can be used to bound the distance between adjacent
   * cell centers along the space-filling Hilbert curve for cells at any given level.
   */
  public final Metric minEdge;

  /** Maximum angular length of any cell edge at level k. */
  public final Metric maxEdge;

  /** Average angular length of any cell edge at level k. */
  public final Metric avgEdge;

  /** Minimum diagonal size of cells at level k. */
  public final Metric minDiag;

  /**
   * Maximum diagonal size of cells at level k. The maximum diagonal also happens to be the maximum
   * diameter of any cell, and also the maximum geometric width. So for example, the distance from
   * an arbitrary point to the closest cell center at a given level is at most half the maximum
   * diagonal length.
   */
  public final Metric maxDiag;

  /** Average diagonal size of cells at level k. */
  public final Metric avgDiag;

  /**
   * Maximum edge aspect ratio over all cells at any level, where the edge aspect ratio of a cell is
   * defined as the ratio of its longest edge length to its shortest edge length.
   */
  public final double maxEdgeAspect;

  /**
   * This is the maximum diagonal aspect ratio over all cells at any level, where the diagonal
   * aspect ratio of a cell is defined as the ratio of its longest diagonal length to its shortest
   * diagonal length.
   */
  public final double maxDiagAspect = Math.sqrt(3); // 1.732

  S2Projections(
      double minAreaDeriv,
      double maxAreaDeriv,
      double minAngleSpanDeriv,
      double maxAngleSpanDeriv,
      double minWidthDeriv,
      double avgWidthDeriv,
      double minEdgeDeriv,
      double avgEdgeDeriv,
      double minDiagDeriv,
      double maxDiagDeriv,
      double avgDiagDeriv,
      double maxEdgeAspect) {
    this.minArea = new Metric(2, minAreaDeriv);
    this.maxArea = new Metric(2, maxAreaDeriv);
    this.avgArea = new Metric(2, 4 * S2.M_PI / 6); // ~2.094
    this.minAngleSpan = new Metric(1, minAngleSpanDeriv);
    this.maxAngleSpan = new Metric(1, maxAngleSpanDeriv);
    this.avgAngleSpan = new Metric(1, S2.M_PI / 2); // ~1.571
    this.minWidth = new Metric(1, minWidthDeriv);
    this.maxWidth = new Metric(1, maxAngleSpanDeriv);
    this.avgWidth = new Metric(1, avgWidthDeriv);
    this.minEdge = new Metric(1, minEdgeDeriv);
    this.maxEdge = new Metric(1, maxAngleSpanDeriv);
    this.avgEdge = new Metric(1, avgEdgeDeriv);
    this.minDiag = new Metric(1, minDiagDeriv);
    this.maxDiag = new Metric(1, maxDiagDeriv);
    this.avgDiag = new Metric(1, avgDiagDeriv);
    this.maxEdgeAspect = maxEdgeAspect;
  }

  /**
   * Convert an s- or t-value to the corresponding u- or v-value. This is a non-linear
   * transformation from [-1,1] to [-1,1] that attempts to make the cell sizes more uniform.
   */
  public abstract double stToUV(double s);

  /**
   * Returns the i- or j-index of the leaf cell containing the given s- or t-value. If the argument
   * is outside the range spanned by valid leaf cell indices, return the index of the closest valid
   * leaf cell (i.e., return values are clamped to the range of valid leaf cell indices).
   */
  public static int stToIj(double s) {
    return Math.max(
        0, Math.min(S2CellId.MAX_SIZE - 1, (int) Math.round(S2CellId.MAX_SIZE * s - 0.5)));
  }

  /**
   * Converts the i- or j-index of a leaf cell to the minimum corresponding s- or t-value contained
   * by that cell. The argument must be in the range [0..2**30], i.e. up to one position beyond the
   * normal range of valid leaf cell indices.
   */
  public static double ijToStMin(int i) {
    // assert (i >= 0 && i <= S2CellId.MAX_SIZE);
    return (1.0 / S2CellId.MAX_SIZE) * i;
  }

  /**
   * Converts the specified i- or j-coordinate into its corresponding u- or v-coordinate,
   * respectively, for the given cell size.
   */
  public double ijToUV(int ij, int cellSize) {
    return stToUV(ijToStMin(ij & -cellSize));
  }

  /** Returns the s- or t-value corresponding to the given si- or ti-value. */
  public static double siTiToSt(long si) {
    // assert (si >= 0 && si <= MAX_SITI);
    return (1.0 / MAX_SITI) * si;
  }

  /**
   * Returns the si- or ti-coordinate that is nearest to the given s- or t-value. The result may be
   * outside the range of valid (si,ti)-values.
   */
  public static long stToSiTi(double s) {
    return Math.round(s * MAX_SITI);
  }

  /**
   * The inverse of {@link #stToUV(double)}. Note that it is not always true that {@code
   * uvToST(stToUV(x)) == x} due to numerical errors.
   */
  public abstract double uvToST(double u);

  /**
   * Convert (face, u, v) coordinates to a direction vector (not necessarily unit length).
   *
   * <p>Requires that the face is between 0 and 5, inclusive.
   */
  public static S2Point faceUvToXyz(int face, double u, double v) {
    XyzTransform t = faceToXyzTransform(face);
    return new S2Point(t.uvToX(u, v), t.uvToY(u, v), t.uvToZ(u, v));
  }

  /**
   * Convert (face, u, v) coordinates to a direction vector (not necessarily unit length).
   *
   * <p>Requires that the face is between 0 and 5, inclusive.
   */
  public static S2Point faceUvToXyz(int face, R2Vector uv) {
    return faceUvToXyz(face, uv.x(), uv.y());
  }

  /** Returns the {@link XyzTransform} for the specified face. */
  static XyzTransform faceToXyzTransform(int face) {
    // We map illegal face indices to the largest face index to preserve legacy behavior, i.e., we
    // do not (yet) want to throw an index out of bounds exception. Note that S2CellId.face() is
    // guaranteed to return a non-negative face index even for invalid S2 cells, so it is sufficient
    // to just map all face indices greater than 5 to a face index of 5.
    //
    // TODO(bjj): Remove this legacy behavior.
    return XYZ_TRANSFORMS[Math.min(5, face)];
  }

  /**
   * If the dot product of p with the given face normal is positive, set the corresponding u and v
   * values (which may lie outside the range [-1,1]) and return true. Otherwise return null.
   */
  public static R2Vector faceXyzToUv(int face, S2Point p) {
    if (face < 3) {
      if (p.get(face) <= 0) {
        return null;
      }
    } else {
      if (p.get(face - 3) >= 0) {
        return null;
      }
    }
    return validFaceXyzToUv(face, p);
  }

  /**
   * Given a *valid* face for the given point p (meaning that dot product of p with the face normal
   * is positive), return the corresponding u and v values (which may lie outside the range [-1,1]).
   *
   * <p>Requires that the face is between 0 and 5, inclusive.
   */
  public static R2Vector validFaceXyzToUv(int face, S2Point p) {
    R2Vector result = new R2Vector();
    validFaceXyzToUv(face, p, result);
    return result;
  }

  /**
   * As {@link #validFaceXyzToUv(int, S2Point)}, except {@code result} is updated, instead of a
   * being returned in a new instance. Package-private because non-S2 classes should not be mutating
   * R2Vectors.
   */
  static void validFaceXyzToUv(int face, S2Point p, R2Vector result) {
    UvTransform t = faceToUvTransform(face);
    result.set(t.xyzToU(p.x, p.y, p.z), t.xyzToV(p.x, p.y, p.z));
  }

  /** Returns the {@link UvTransform} for the specified face. */
  public static UvTransform faceToUvTransform(int face) {
    return UV_TRANSFORMS[face];
  }

  /**
   * Returns the given point P transformed to the (u,v,w) coordinate frame of the given face (where
   * the w-axis represents the face normal).
   */
  public static S2Point faceXyzToUvw(int face, S2Point p) {
    // The result coordinates are simply the dot products of P with the (u,v,w) axes for the given
    // face (see FACE_UVW_AXES).
    switch (face) {
      case 0:
        return new S2Point(p.y, p.z, p.x);
      case 1:
        return new S2Point(-p.x, p.z, p.y);
      case 2:
        return new S2Point(-p.x, -p.y, p.z);
      case 3:
        return new S2Point(-p.z, -p.y, -p.x);
      case 4:
        return new S2Point(-p.z, p.x, -p.y);
      default:
        return new S2Point(p.y, p.x, -p.z);
    }
  }

  /** Returns the level of the given si or ti coordinate. */
  private static final int siTiToLevel(long siTi) {
    return S2CellId.MAX_LEVEL - Long.numberOfTrailingZeros(siTi | MAX_SITI);
  }

  /**
   * A [face, si, ti] position. This is package private for now, since we may want to rework the
   * class to use 32-bit ints instead.
   */
  static final class FaceSiTi {
    /** The face on which the position exists. */
    public final int face;
    /** The si coordinate. See {@link S2Projections} for details. */
    public final long si;
    /** The ti coordinate. See {@link S2Projections} for details. */
    public final long ti;

    /** Package private constructor. Only S2 should create these for now. */
    FaceSiTi(int face, long si, long ti) {
      this.face = face;
      this.si = si;
      this.ti = ti;
    }
  }

  /** Convert (face, si, ti) coordinates to a direction vector (not necessarily unit length.) */
  public S2Point faceSiTiToXyz(int face, long si, long ti) {
    double u = stToUV(siTiToSt(si));
    double v = stToUV(siTiToSt(ti));
    return faceUvToXyz(face, u, v);
  }

  /** Convert a direction vector (not necessarily unit length) to (face, si, ti) coordinates. */
  FaceSiTi xyzToFaceSiTi(S2Point p) {
    int face = xyzToFace(p);
    R2Vector uv = validFaceXyzToUv(face, p);
    long si = stToSiTi(uvToST(uv.x()));
    long ti = stToSiTi(uvToST(uv.y()));
    return new FaceSiTi(face, si, ti);
  }

  /** If p is exactly a cell center, returns the level of the cell, -1 otherwise. */
  int levelIfCenter(FaceSiTi fst, S2Point p) {
    // If the levels corresponding to si,ti are not equal, then p is not a cell
    // center.  The si,ti values 0 and MAX_SITI need to be handled specially
    // because they do not correspond to cell centers at any valid level; they
    // are mapped to level -1 by the code below.
    int level = siTiToLevel(fst.si);
    if (level < 0 || level != siTiToLevel(fst.ti)) {
      return -1;
    } else {
      // assert (level <= S2CellId.MAX_LEVEL);
      // In infinite precision, this test could be changed to ST == SiTi. However,
      // due to rounding errors, UVtoST(XYZtoFaceUV(FaceUVtoXYZ(STtoUV(...)))) is
      // not idempotent. On the other hand, centerRaw is computed exactly the same
      // way p was originally computed (if it is indeed the center of an S2Cell):
      // the comparison can be exact.
      S2Point center = S2Point.normalize(faceSiTiToXyz(fst.face, fst.si, fst.ti));
      if (p.equals(center)) {
        return level;
      } else {
        return -1;
      }
    }
  }

  /**
   * Returns the face containing the given direction vector (for points on the boundary between
   * faces, the result is arbitrary but repeatable.)
   */
  public static int xyzToFace(S2Point p) {
    return xyzToFace(p.x, p.y, p.z);
  }

  /**
   * As {@link #xyzToFace(S2Point)}, but accepts the coordinates as primitive doubles instead.
   * Useful when the caller has coordinates and doesn't want to allocate an S2Point.
   */
  static int xyzToFace(double x, double y, double z) {
    switch (S2Point.largestAbsComponent(x, y, z)) {
      case 0:
        return (x < 0) ? 3 : 0;
      case 1:
        return (y < 0) ? 4 : 1;
      default:
        return (z < 0) ? 5 : 2;
    }
  }

  /**
   * Returns the right-handed normal (not necessarily unit length) for an edge in the direction of
   * the positive v-axis at the given u-value on the given face. (This vector is perpendicular to
   * the plane through the sphere origin that contains the given edge.)
   */
  public static S2Point getUNorm(int face, double u) {
    switch (face) {
      case 0:
        return new S2Point(u, -1, 0);
      case 1:
        return new S2Point(1, u, 0);
      case 2:
        return new S2Point(1, 0, u);
      case 3:
        return new S2Point(-u, 0, 1);
      case 4:
        return new S2Point(0, -u, 1);
      default:
        return new S2Point(0, -1, -u);
    }
  }

  /**
   * Returns the right-handed normal (not necessarily unit length) for an edge in the direction of
   * the positive u-axis at the given v-value on the given face.
   */
  public static S2Point getVNorm(int face, double v) {
    switch (face) {
      case 0:
        return new S2Point(-v, 0, 1);
      case 1:
        return new S2Point(0, -v, 1);
      case 2:
        return new S2Point(0, -1, -v);
      case 3:
        return new S2Point(v, -1, 0);
      case 4:
        return new S2Point(1, v, 0);
      default:
        return new S2Point(1, 0, v);
    }
  }

  /** Returns the u-axis for the given face. */
  public static S2Point getUAxis(int face) {
    return getUVWAxis(face, 0);
  }

  /** Returns the v-axis for the given face. */
  public static S2Point getVAxis(int face) {
    return getUVWAxis(face, 1);
  }

  /** Returns the unit-length normal for the given face. */
  public static S2Point getNorm(int face) {
    return getUVWAxis(face, 2);
  }

  /** Returns the given axis of the given face (u=0, v=1, w=2). */
  static S2Point getUVWAxis(int face, int axis) {
    return FACE_UVW_AXES[face][axis];
  }

  /**
   * Returns the face that lies in the given direction (negative=0, positive=1) of the given axis
   * (u=0, v=1, w=2) in the given face. For example, {@code getUVWFace(4, 0, 1)} returns the face
   * that is adjacent to face 4 in the positive u-axis direction.
   */
  static int getUVWFace(int face, int axis, int direction) {
    // assert (face >= 0 && face <= 5);
    // assert (axis >= 0 && axis <= 2);
    // assert (direction >= 0 && direction <= 1);
    return FACE_UVW_FACES[face][axis][direction];
  }

  /** The default transformation between ST and UV coordinates. */
  public static final S2Projections PROJ = S2Projections.S2_QUADRATIC_PROJECTION;
}
