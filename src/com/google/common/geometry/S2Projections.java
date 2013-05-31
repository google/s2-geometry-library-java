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
 * This class specifies the details of how the cube faces are projected onto the
 * unit sphere. This includes getting the face ordering and orientation correct
 * so that sequentially increasing cell ids follow a continuous space-filling
 * curve over the entire sphere, and defining the transformation from cell-space
 * to cube-space (see s2.h) in order to make the cells more uniform in size.
 *
 *
 *  We have implemented three different projections from cell-space (s,t) to
 * cube-space (u,v): linear, quadratic, and tangent. They have the following
 * tradeoffs:
 *
 *  Linear - This is the fastest transformation, but also produces the least
 * uniform cell sizes. Cell areas vary by a factor of about 5.2, with the
 * largest cells at the center of each face and the smallest cells in the
 * corners.
 *
 *  Tangent - Transforming the coordinates via atan() makes the cell sizes more
 * uniform. The areas vary by a maximum ratio of 1.4 as opposed to a maximum
 * ratio of 5.2. However, each call to atan() is about as expensive as all of
 * the other calculations combined when converting from points to cell ids, i.e.
 * it reduces performance by a factor of 3.
 *
 *  Quadratic - This is an approximation of the tangent projection that is much
 * faster and produces cells that are almost as uniform in size. It is about 3
 * times faster than the tangent projection for converting cell ids to points,
 * and 2 times faster for converting points to cell ids. Cell areas vary by a
 * maximum ratio of about 2.1.
 *
 *  Here is a table comparing the cell uniformity using each projection. "Area
 * ratio" is the maximum ratio over all subdivision levels of the largest cell
 * area to the smallest cell area at that level, "edge ratio" is the maximum
 * ratio of the longest edge of any cell to the shortest edge of any cell at the
 * same level, and "diag ratio" is the ratio of the longest diagonal of any cell
 * to the shortest diagonal of any cell at the same level. "ToPoint" and
 * "FromPoint" are the times in microseconds required to convert cell ids to and
 * from points (unit vectors) respectively.
 *
 *  Area Edge Diag ToPoint FromPoint Ratio Ratio Ratio (microseconds)
 * ------------------------------------------------------- Linear: 5.200 2.117
 * 2.959 0.103 0.123 Tangent: 1.414 1.414 1.704 0.290 0.306 Quadratic: 2.082
 * 1.802 1.932 0.116 0.161
 *
 *  The worst-case cell aspect ratios are about the same with all three
 * projections. The maximum ratio of the longest edge to the shortest edge
 * within the same cell is about 1.4 and the maximum ratio of the diagonals
 * within the same cell is about 1.7.
 *
 * This data was produced using s2cell_unittest and s2cellid_unittest.
 *
 */
@GwtCompatible
public strictfp enum S2Projections {
  // All of the values below were obtained by a combination of hand analysis and
  // Mathematica. In general, S2_TAN_PROJECTION produces the most uniform
  // shapes and sizes of cells, S2_LINEAR_PROJECTION is considerably worse, and
  // S2_QUADRATIC_PROJECTION is somewhere in between (but generally closer to
  // the tangent projection than the linear one).
  S2_LINEAR_PROJECTION(
      1 / (3 * Math.sqrt(3)), // minAreaDeriv ~0.192
      1, // maxAreaDeriv ~1.000
      0.5, // minAngleSpanDeriv ~0.500
      1, // maxAngleSpanDeriv ~1.000
      1 / Math.sqrt(6), // minWidthDeriv ~0.408
      0.70572967292222848, // avgWidthDeriv ~0.706
      S2.M_SQRT2 / 3, // minEdgeDeriv ~0.471
      0.72001709647780182, // avgEdgeDeriv ~0.720
      S2.M_SQRT2 / 3, // minDiagDeriv ~0.471
      S2.M_SQRT2, // maxDiagDeriv ~1.414
      1.0159089332094063, // avgDiagDeriv ~1.016
      S2.M_SQRT2) { // maxEdgeAspect ~1.414
    @Override
    public double stToUV(double s) {
      return s;
    }
    @Override
    public double uvToST(double u) {
      return u;
    }
  },
  S2_TAN_PROJECTION(
      (S2.M_PI * S2.M_PI) / (16 * S2.M_SQRT2), // minAreaDeriv ~0.436
      S2.M_PI * S2.M_PI / 16, // maxAreaDeriv ~0.617
      S2.M_PI / 4, // minAngleSpanDeriv ~0.785
      S2.M_PI / 4, // maxAngleSpanDeriv ~0.785
      S2.M_PI / (4 * S2.M_SQRT2), // minWidthDeriv ~0.555
      0.71865931946258044, // avgWidthDeriv ~0.719
      S2.M_PI / (4 * S2.M_SQRT2), // minEdgeDeriv ~0.555
      0.73083351627336963, // avgEdgeDeriv ~0.731
      S2.M_PI / (3 * S2.M_SQRT2), // minDiagDeriv ~0.740
      S2.M_PI / Math.sqrt(6), // maxDiagDeriv ~1.283
      1.0318115985978178, // avgDiagDeriv ~1.032
      S2.M_SQRT2) { // maxEdgeAspect ~1.414
    @Override
    public double stToUV(double s) {
      // Unfortunately, tan(M_PI_4) is slightly less than 1.0. This isn't due to a flaw in the
      // implementation of tan(), it's because the derivative of tan(x) at x=pi/4 is 2, and it
      // happens that the two adjacent floating point numbers on either side of the infinite-
      // precision value of pi/4 have tangents that are slightly below and slightly above 1.0 when
      // rounded to the nearest double-precision result.
      s = Math.tan(S2.M_PI_4 * s);
      return s + (1.0 / (1L << 53)) * s;
    }
    @Override
    public double uvToST(double u) {
      return (4 * S2.M_1_PI) * Math.atan(u);
    }
  },
  S2_QUADRATIC_PROJECTION(
      2 * S2.M_SQRT2 / 9, // minAreaDeriv ~0.314
      0.65894981424079037, // maxAreaDeriv ~0.659
      2. / 3, // minAngleSpanDeriv ~0.667
      0.85244858959960922, // maxAngleSpanDeriv ~0.852
      S2.M_SQRT2 / 3, // minWidthDeriv ~0.471
      0.71726183644304969, // avgWidthDeriv ~0.717
      S2.M_SQRT2 / 3, // minEdgeDeriv ~0.471
      0.72960687319305303, // avgEdgeDeriv ~0.730
      4 * S2.M_SQRT2 / 9, // minDiagDeriv ~0.629
      1.2193272972170106, // maxDiagDeriv ~1.219
      1.03021136949923584, // avgDiagDeriv ~1.030
      1.44261527445268292) { // maxEdgeAspect ~1.443
    @Override
    public double stToUV(double s) {
      if (s >= 0) {
        return (1 / 3.) * ((1 + s) * (1 + s) - 1);
      } else {
        return (1 / 3.) * (1 - (1 - s) * (1 - s));
      }
    }
    @Override
    public double uvToST(double u) {
      if (u >= 0) {
        return Math.sqrt(1 + 3 * u) - 1;
      } else {
        return 1 - Math.sqrt(1 - 3 * u);
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
   * Minimum angular separation between opposite edges of a cell at level k. Each cell is bounded
   * by four planes passing through its four edges and the center of the sphere. The angle span
   * metrics relate to the angle between each pair of opposite bounding planes, or equivalently,
   * between the planes corresponding to two different s-values or two different t-values.
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
   * Minimum angular length of any cell edge at level k. The edge length metrics can also be used
   * to bound the minimum, maximum, or average distance from the center of one cell to the center
   * of one of its edge neighbors. In particular, it can be used to bound the distance between
   * adjacent cell centers along the space-filling Hilbert curve for cells at any given level.
   */
  public final Metric minEdge;

  /** Maximum angular length of any cell edge at level k. */
  public final Metric maxEdge;

  /** Average angular length of any cell edge at level k. */
  public final Metric avgEdge;

  /** Minimum diagonal size of cells at level k. */
  public final Metric minDiag;

  /**
   * Maximum diagonal size of cells at level k. The maximum diagonal also happens to be the
   * maximum diameter of any cell, and also the maximum geometric width. So for example, the
   * distance from an arbitrary point to the closest cell center at a given level is at most half
   * the maximum diagonal length.
   */
  public final Metric maxDiag;

  /** Average diagonal size of cells at level k. */
  public final Metric avgDiag;

  /**
   * Maximum edge aspect ratio over all cells at any level, where the edge aspect ratio of a cell
   * is defined as the ratio of its longest edge length to its shortest edge length.
   */
  public final double maxEdgeAspect;

  /**
   * This is the maximum diagonal aspect ratio over all cells at any level, where the diagonal
   * aspect ratio of a cell is defined as the ratio of its longest diagonal length to its
   * shortest diagonal length.
   */
  public final double maxDiagAspect = Math.sqrt(3); // 1.732

  S2Projections(double minAreaDeriv, double maxAreaDeriv,
      double minAngleSpanDeriv, double maxAngleSpanDeriv,
      double minWidthDeriv, double avgWidthDeriv,
      double minEdgeDeriv, double avgEdgeDeriv,
      double minDiagDeriv, double maxDiagDeriv, double avgDiagDeriv,
      double maxEdgeAspect) {
    this.minArea = new Metric(2, minAreaDeriv);
    this.maxArea = new Metric(2, maxAreaDeriv);
    this.avgArea = new Metric(2, S2.M_PI / 6); // 0.524)
    this.minAngleSpan = new Metric(1, minAngleSpanDeriv);
    this.maxAngleSpan = new Metric(1, maxAngleSpanDeriv);
    this.avgAngleSpan = new Metric(1, S2.M_PI / 4); // 0.785
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
   * Convert an s- or t-value to the corresponding u- or v-value.  This is a non-linear
   * transformation from [-1,1] to [-1,1] that attempts to make the cell sizes more uniform.
   */
  public abstract double stToUV(double s);

  /**
   * The inverse of {@link #stToUV(double)}. Note that it is not always true that
   * {@code uvToST(stToUV(x)) == x} due to numerical errors.
   */
  public abstract double uvToST(double u);

  /**
   * Convert (face, u, v) coordinates to a direction vector (not necessarily
   * unit length).
   */
  public static S2Point faceUvToXyz(int face, double u, double v) {
    switch (face) {
      case 0:
        return new S2Point(1, u, v);
      case 1:
        return new S2Point(-u, 1, v);
      case 2:
        return new S2Point(-u, -v, 1);
      case 3:
        return new S2Point(-1, -v, -u);
      case 4:
        return new S2Point(v, -1, -u);
      default:
        return new S2Point(v, u, -1);
    }
  }

  public static R2Vector validFaceXyzToUv(int face, S2Point p) {
    // assert (p.dotProd(faceUvToXyz(face, 0, 0)) > 0);
    double pu;
    double pv;
    switch (face) {
      case 0:
        pu = p.y / p.x;
        pv = p.z / p.x;
        break;
      case 1:
        pu = -p.x / p.y;
        pv = p.z / p.y;
        break;
      case 2:
        pu = -p.x / p.z;
        pv = -p.y / p.z;
        break;
      case 3:
        pu = p.z / p.x;
        pv = p.y / p.x;
        break;
      case 4:
        pu = p.z / p.y;
        pv = -p.x / p.y;
        break;
      default:
        pu = -p.y / p.z;
        pv = -p.x / p.z;
        break;
    }
    return new R2Vector(pu, pv);
  }

  public static int xyzToFace(S2Point p) {
    int face = p.largestAbsComponent();
    if (p.get(face) < 0) {
      face += 3;
    }
    return face;
  }

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

  public static S2Point getNorm(int face) {
    return faceUvToXyz(face, 0, 0);
  }

  public static S2Point getUAxis(int face) {
    switch (face) {
      case 0:
        return S2Point.Y_POS;
      case 1:
        return S2Point.X_NEG;
      case 2:
        return S2Point.X_NEG;
      case 3:
        return S2Point.Z_NEG;
      case 4:
        return S2Point.Z_NEG;
      default:
        return S2Point.Y_POS;
    }
  }

  public static S2Point getVAxis(int face) {
    switch (face) {
      case 0:
        return S2Point.Z_POS;
      case 1:
        return S2Point.Z_POS;
      case 2:
        return S2Point.Y_NEG;
      case 3:
        return S2Point.Y_NEG;
      case 4:
        return S2Point.X_POS;
      default:
        return S2Point.X_POS;
    }
  }

  /** The default transformation between ST and UV coordinates. */
  public static final S2Projections PROJ = S2Projections.S2_QUADRATIC_PROJECTION;
}
