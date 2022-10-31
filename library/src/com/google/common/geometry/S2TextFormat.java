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

import static java.lang.Double.parseDouble;

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import javax.annotation.Nullable;

/**
 * S2TextFormat contains a collection of functions for converting geometry to and from a human-
 * readable format. It is mainly intended for testing and debugging. Be aware that the human-
 * readable format is *not* designed to preserve the full precision of the original object, so it
 * should not be used for data storage.
 */
public strictfp class S2TextFormat {

  /**
   * Returns an S2Point corresponding to the given a latitude-longitude coordinate in degrees.
   * Example of the input format: "-20:150"
   */
  public static S2Point makePointOrDie(String str) {
    S2Point point = makePoint(str);
    Preconditions.checkState(null != point, ": str == \"%s\"", str);
    return point;
  }

  /**
   * As above, but do not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2Point makePoint(String str) {
    List<S2Point> vertices = parsePoints(str);
    if (vertices == null || vertices.size() != 1) {
      return null;
    }
    return vertices.get(0);
  }

  /**
   * Parses a string of one or more latitude-longitude coordinates in degrees, and return the
   * corresponding List of S2LatLng points.
   *
   * <p>Examples of the input format:
   *
   * <pre>
   *     ""                          // no points
   *     "-20:150"                   // one point
   *     "-20:150, -20:151, -19:150" // three points
   * </pre>
   */
  public static List<S2LatLng> parseLatLngsOrDie(String str) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    Preconditions.checkState(latlngs != null, ": str == \"%s\"", str);
    return latlngs;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  static List<S2LatLng> parseLatLngs(String str) {
    List<ParseEntry> ps = dictionaryParse(str);
    if (ps == null) {
      return null;
    }
    List<S2LatLng> latlngs = new ArrayList<>();
    for (ParseEntry p : ps) {
      Double lat = parseDouble(p.key);
      if (lat == null) {
        return null;
      }
      Double lng = parseDouble(p.value);
      if (lng == null) {
        return null;
      }
      latlngs.add(S2LatLng.fromDegrees(lat, lng));
    }
    return latlngs;
  }

  /**
   * Parses a string in the same format as parseLatLngs, and return the corresponding List of
   * S2Point values.
   */
  public static List<S2Point> parsePointsOrDie(String str) {
    List<S2Point> vertices = parsePoints(str);
    Preconditions.checkState(vertices != null, ": str == \"%s\"", str);
    return vertices;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static List<S2Point> parsePoints(String str) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    if (latlngs == null) {
      return null;
    }
    List<S2Point> vertices = new ArrayList<>();
    for (S2LatLng latlng : latlngs) {
      vertices.add(latlng.toPoint());
    }
    return vertices;
  }

  /** Given a string in the same format as ParseLatLngs, returns a single S2LatLng. */
  public static S2LatLng makeLatLngOrDie(String str) {
    S2LatLng latlng = makeLatLng(str);
    Preconditions.checkState(null != latlng, ": str == \"%s\"", str);
    return latlng;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2LatLng makeLatLng(String str) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    if (null == latlngs || latlngs.size() != 1) {
      return null;
    }
    return latlngs.get(0);
  }

  /**
   * Given a string in the same format as ParseLatLngs, returns the minimal bounding S2LatLngRect
   * that contains the coordinates.
   */
  S2LatLngRect makeLatLngRectOrDie(String str) {
    S2LatLngRect rect = makeLatLngRect(str);
    Preconditions.checkState(null != rect, ": str == \"%s\"", str);
    return rect;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2LatLngRect makeLatLngRect(String str) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    if (null == latlngs || latlngs.isEmpty()) {
      return null;
    }
    S2LatLngRect rect = S2LatLngRect.fromPoint(latlngs.get(0));
    for (int i = 1; i < latlngs.size(); ++i) {
      rect = rect.addPoint(latlngs.get(i));
    }
    return rect;
  }

  /**
   * Parses an S2CellId in the format "f/dd..d" where "f" is a digit in the range [0-5] representing
   * the S2CellId face, and "dd..d" is a string of digits in the range [0-3] representing each
   * child's position with respect to its parent. (Note that the latter string may be empty.)
   *
   * <p>For example "4/" represents S2CellId.fromFace(4), and "3/02" represents
   * S2CellId.fromFace(3).child(0).child(2).
   *
   * <p>This function is a wrapper for S2CellId.fromDebugString().
   */
  public static S2CellId makeCellIdOrDie(String str) {
    S2CellId cellId = makeCellId(str);
    Preconditions.checkState(null != cellId, ": str == \"%s\"", str);
    return cellId;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2CellId makeCellId(String str) {
    S2CellId cellId = S2CellId.fromDebugString(str);
    if (cellId.equals(S2CellId.none())) {
      return null;
    }
    return cellId;
  }

  /**
   * Parses a comma-separated list of S2CellIds in the format above, and returns the corresponding
   * S2CellUnion. (Note that S2CellUnions are automatically normalized by sorting, removing
   * duplicates, and replacing groups of 4 child cells by their parent cell.)
   */
  public static S2CellUnion makeCellUnionOrDie(String str) {
    S2CellUnion cellUnion = makeCellUnion(str);
    Preconditions.checkState(null != cellUnion, ": str == \"%s\"", str);
    return cellUnion;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2CellUnion makeCellUnion(String str) {
    ArrayList<S2CellId> cellIds = new ArrayList<>();
    for (String cellStr : splitString(str, ",")) {
      cellStr = cellStr.trim();
      S2CellId cellId = makeCellId(cellStr);
      if (null == cellId) {
        return null;
      }
      cellIds.add(cellId);
    }
    S2CellUnion cellUnion = new S2CellUnion();
    cellUnion.initFromCellIds(cellIds);
    return cellUnion;
  }

  /**
   * Given a string of latitude-longitude coordinates in degrees, returns a newly allocated loop.
   * Example of the input format:
   *
   * <p>"-20:150, 10:-120, 0.123:-170.652"
   *
   * <p>The strings "empty" or "full" create an empty or full loop respectively.
   */
  public static S2Loop makeLoopOrDie(String str) {
    S2Loop loop = makeLoop(str);
    Preconditions.checkState(loop != null, ": str == \"%s\"", str);
    return loop;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2Loop makeLoop(String str) {
    if (str.equals("empty")) {
      return S2Loop.empty();
    }
    if (str.equals("full")) {
      return S2Loop.full();
    }
    List<S2Point> vertices = parsePoints(str);
    if (vertices == null) {
      return null;
    }
    return new S2Loop(vertices);
  }

  /** Similar to makeLoop(), but returns an S2Polyline rather than an S2Loop. */
  public static S2Polyline makePolylineOrDie(String str) {
    S2Polyline polyline = makePolyline(str);
    Preconditions.checkState(null != polyline, ": str == \"%s\"", str);
    return polyline;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2Polyline makePolyline(String str) {
    List<S2Point> vertices = parsePoints(str);
    if (null == vertices) {
      return null;
    }
    return new S2Polyline(vertices);
  }

  /** Like makePolyline, but returns an S2LaxPolylineShape instead. */
  public static S2LaxPolylineShape makeLaxPolylineOrDie(String str) {
    S2LaxPolylineShape laxPolyline = makeLaxPolyline(str);
    Preconditions.checkState(null != laxPolyline, ": str == \"%s\"", str);
    return laxPolyline;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2LaxPolylineShape makeLaxPolyline(String str) {
    List<S2Point> vertices = parsePoints(str);
    if (null == vertices) {
      return null;
    }
    return S2LaxPolylineShape.create(vertices);
  }

  /**
   * Given a sequence of loops separated by semicolons, returns a newly allocated polygon. Loops are
   * automatically normalized by inverting them if necessary so that they enclose at most half of
   * the unit sphere. (Historically this was once a requirement of polygon loops. It also hides the
   * problem that if the user thinks of the coordinates as X:Y rather than LAT:LNG, it yields a loop
   * with the opposite orientation.)
   *
   * <p>Examples of the input format:
   *
   * <pre>
   *     "10:20, 90:0, 20:30"                                 // one loop
   *     "10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2:20.3"  // two loops
   *     ""                   // the empty polygon (consisting of no loops)
   *     "empty"              // the empty polygon (consisting of no loops)
   *     "full"               // the full polygon (consisting of one full loop).
   * </pre>
   */
  public static S2Polygon makePolygonOrDie(String str) {
    S2Polygon polygon = makePolygon(str);
    Preconditions.checkState(polygon != null, ": str == \"%s\"", str);
    return polygon;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2Polygon makePolygon(String str) {
    return internalMakePolygon(str, true);
  }

  /**
   * Like MakePolygon(), except that it does not normalize loops (i.e., it gives you exactly what
   * you asked for).
   */
  public static S2Polygon makeVerbatimPolygonOrDie(String str) {
    S2Polygon polygon = makeVerbatimPolygon(str);
    Preconditions.checkState(polygon != null, ": str == \"%s\"", str);
    return polygon;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2Polygon makeVerbatimPolygon(String str) {
    return internalMakePolygon(str, false);
  }

  /**
   * Parses a string in the same format as MakePolygon, except that loops must be oriented so that
   * the interior of the loop is always on the left, and polygons with degeneracies are supported.
   * As with MakePolygon, "full" denotes the full polygon, and "" or "empty" denote the empty
   * polygon.
   */
  public static S2LaxPolygonShape makeLaxPolygonOrDie(String str) {
    S2LaxPolygonShape laxPolygon = makeLaxPolygon(str);
    Preconditions.checkState(null != laxPolygon, ": str == \"%s\"", str);
    return laxPolygon;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2LaxPolygonShape makeLaxPolygon(String str) {
    List<String> loopStrs = splitString(str, ";");
    List<List<S2Point>> loops = new ArrayList<>();
    for (String loopStr : loopStrs) {
      if (loopStr.equals("full")) {
        loops.add(new ArrayList<S2Point>());
      } else if (!loopStr.equals("empty")) {
        List<S2Point> points = parsePoints(loopStr);
        if (null == points) {
          return null;
        }
        loops.add(points);
      }
    }
    return S2LaxPolygonShape.create(loops);
  }

  /**
   * Returns a S2ShapeIndex containing the points, polylines, and loops (in the form of one polygon
   * for each group of loops) described by the following format:
   *
   * <p>point1|point2|... # line1|line2|... # polygon1|polygon2|...
   *
   * <p>Examples:
   *
   * <ul>
   * <li>1:2 | 2:3 # #                     // Two points (one S2Point.Shape)
   * <li># 0:0, 1:1, 2:2 | 3:3, 4:4 #      // Two polylines
   * <li># # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2  // Two nested loops (one polygon)
   * <li>5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0    // One of each point, line, and polygon
   * <li># # empty                         // One empty polygon
   * <li># # empty | full                  // One empty polygon, one full polygon
   * </ul>
   *
   * <p>All the points, if any, are stored as a single S2Point.Shape in the index. Polylines are
   * stored as individual S2LaxPolylineShapes. Polygons are separated by '|', with distinct loops
   * for a polygon separated by ';'. Each group of loops is stored as an individual
   * S2LaxPolygonShape.
   *
   * <p>Loops should be directed so that the region's interior is on the left. Loops can be
   * degenerate (they do not need to meet S2Loop requirements).
   *
   * <p>CAVEAT: Because whitespace is ignored, empty polygons must be specified as the string
   * "empty" rather than as the empty string ("").
   */
  public static S2ShapeIndex makeIndexOrDie(String str) {
    S2ShapeIndex index = makeIndex(str);
    Preconditions.checkState(index != null, ": str == \"%s\"", str);
    return index;
  }

  /**
   * As above, but does not CHECK-fail on invalid input. Returns null if conversion is unsuccessful.
   */
  @Nullable
  public static S2ShapeIndex makeIndex(String str) {
    String[] strs = str.split("#", -1); // Here, we want to include empty strings in split.
    if (strs.length != 3) {
      return null;
    }
    List<S2Point> points = new ArrayList<>();
    for (String pointStr : splitString(strs[0], "\\|")) {
      S2Point point = makePoint(pointStr);
      if (point == null) {
        return null;
      }
      points.add(point);
    }

    S2ShapeIndex index = new S2ShapeIndex();
    if (!points.isEmpty()) {
      index.add(S2Point.Shape.fromList(points));
    }

    for (String lineStr : splitString(strs[1], "\\|")) {
      S2LaxPolylineShape laxPolyline = makeLaxPolyline(lineStr);
      if (laxPolyline == null) {
        return null;
      }
      index.add(laxPolyline);
    }

    for (String polygonStr : splitString(strs[2], "\\|")) {
      S2LaxPolygonShape laxPolygon = makeLaxPolygon(polygonStr);
      if (laxPolygon == null) {
        return null;
      }
      index.add(laxPolygon);
    }
    return index;
  }

  /** Convert an S2Point to the S2TextFormat string representation documented above. */
  public static String toString(S2Point s2Point) {
    StringBuilder out = new StringBuilder();
    appendVertex(s2Point, out);
    return out.toString();
  }

  /** Convert an S2LatLng to the S2TextFormat string representation documented above. */
  public static String toString(S2LatLng latlng) {
    StringBuilder out = new StringBuilder();
    appendVertex(latlng, out);
    return out.toString();
  }

  /** Convert an S2LatLngRect to the S2TextFormat string representation documented above. */
  public static String toString(S2LatLngRect rect) {
    StringBuilder out = new StringBuilder();
    appendVertex(rect.lo(), out);
    out.append(", ");
    appendVertex(rect.hi(), out);
    return out.toString();
  }

  /** Convert an S2CellId to the S2TextFormat string representation documented above. */
  public static String toString(S2CellId cellId) {
    return cellId.toString();
  }

  /** Convert an S2CellUnion to the S2TextFormat string representation documented above. */
  public static String toString(S2CellUnion cellUnion) {
    StringBuilder out = new StringBuilder();
    for (S2CellId cellId : cellUnion) {
      if (out.length() > 0) {
        out.append(", ");
      }
      out.append(cellId.toString());
    }
    return out.toString();
  }

  /** Convert an S2Loop to the S2TextFormat string representation documented above. */
  public static String toString(S2Loop loop) {
    if (loop.isEmpty()) {
      return "empty";
    } else if (loop.isFull()) {
      return "full";
    }
    StringBuilder out = new StringBuilder();
    if (loop.numVertices() > 0) {
      appendVertices(loop.vertices(), out);
    }
    return out.toString();
  }

  /** Convert an S2Polyline to the S2TextFormat string representation documented above. */
  public static String toString(S2Polyline polyline) {
    StringBuilder out = new StringBuilder();
    if (polyline.numVertices() > 0) {
      appendVertices(polyline.vertices(), out);
    }
    return out.toString();
  }

  /** Convert an S2Polygon to the S2TextFormat string representation documented above. */
  public static String toString(S2Polygon polygon) {
    return toString(polygon, ";\n");
  }

  /**
   * Convert an S2 polygon to the S2TextFormat string representation documented above, using the
   * given loopSeparator between each loop. Empty and Full polygons are represented as "empty" and
   * "full" respectively.
   */
  public static String toString(S2Polygon polygon, String loopSeparator) {
    if (polygon.isEmpty()) {
      return "empty";
    } else if (polygon.isFull()) {
      return "full";
    }
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < polygon.numLoops(); ++i) {
      if (i > 0) {
        out.append(loopSeparator);
      }
      appendVertices(polygon.loop(i).vertices(), out);
    }
    return out.toString();
  }

  /** Convert a list of S2Points to the S2TextFormat string representation documented above. */
  public static String s2PointsToString(List<S2Point> points) {
    StringBuilder out = new StringBuilder();
    appendVertices(points, out);
    return out.toString();
  }

  /** Convert a list of S2LatLngs to the S2TextFormat string representation documented above. */
  public static String s2LatLngsToString(List<S2LatLng> latlngs) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < latlngs.size(); ++i) {
      if (i > 0) {
        out.append(", ");
      }
      appendVertex(latlngs.get(i), out);
    }
    return out.toString();
  }

  /** Convert an S2LaxPolylineShape to the S2TextFormat string representation documented above. */
  public static String toString(S2LaxPolylineShape polyline) {
    StringBuilder out = new StringBuilder();
    if (polyline.numVertices() > 0) {
      appendVertices(polyline.vertices(), out);
    }
    return out.toString();
  }

  /** Convert an S2LaxPolygonShape to the S2TextFormat string representation documented above. */
  public static String toString(S2LaxPolygonShape polygon) {
    return toString(polygon, ";\n");
  }

  /**
   * Convert an S2LaxPolygonShape to the S2TextFormat string representation documented above, using
   * the given loopSeparator.
   */
  public static String toString(S2LaxPolygonShape polygon, String loopSeparator) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < polygon.numChains(); ++i) {
      if (i > 0) {
        out.append(loopSeparator);
      }
      int chainLength = polygon.getChainLength(i);
      if (chainLength == 0) {
        out.append("full");
      } else {
        for (int edgeOffset = 0; edgeOffset < chainLength; edgeOffset++) {
          appendVertex(polygon.getChainVertex(i, edgeOffset), out);
          if (edgeOffset < chainLength - 1) {
            out.append(", ");
          }
        }
      }
    }
    return out.toString();
  }

  /**
   * Convert an S2CellUnion to the S2TextFormat string representation documented above. The index
   * may contain S2Shapes of any type. Shapes are reordered if necessary so that all point geometry
   * (shapes of dimension 0) are first, followed by all polyline geometry, followed by all polygon
   * geometry.
   */
  public static String toString(S2ShapeIndex index) {
    StringBuilder out = new StringBuilder();
    MutableEdge edge = new MutableEdge();

    for (int dim = 0; dim < 3; ++dim) {
      if (dim > 0) {
        out.append("#");
      }
      int count = 0;
      for (S2Shape shape : index.getShapes()) {
        if (shape == null || shape.dimension() != dim) {
          continue;
        }
        out.append((count > 0) ? " | " : (dim > 0) ? " " : "");
        for (int i = 0; i < shape.numChains(); ++i, ++count) {
          if (i > 0) {
            out.append((dim == 2) ? "; " : " | ");
          }
          if (shape.getChainLength(i) == 0) {
            out.append("full");
          } else {
            shape.getChainEdge(i, 0, edge);
            appendVertex(edge.getStart(), out);
          }
          int limit = shape.getChainLength(i);
          if (dim != 1) {
            --limit;
          }
          for (int edgeOffset = 0; edgeOffset < limit; ++edgeOffset) {
            out.append(", ");
            shape.getChainEdge(i, edgeOffset, edge);
            appendVertex(edge.getEnd(), out);
          }
        }
      }
      // Example output: "# #", "0:0 # #", "# # 0:0, 0:1, 1:0"
      if (dim == 1 || (dim == 0 && count > 0)) {
        out.append(" ");
      }
    }
    return out.toString();
  }

  // Split on the given regexp. Trim whitespace and skip empty strings to produce the result.
  private static List<String> splitString(String str, String regexp) {
    String[] parts = str.split(regexp);
    List<String> result = new ArrayList<>();
    for (String part : parts) {
      if (!part.trim().isEmpty()) {
        result.add(part.trim());
      }
    }
    return result;
  }

  private static class ParseEntry {
    public String key;
    public String value;

    public ParseEntry(String k, String v) {
      this.key = k;
      this.value = v;
    }
  }

  /** Modeled on the DictionaryParse method of strings/serialize.cc */
  @Nullable
  private static List<ParseEntry> dictionaryParse(String str) {
    List<ParseEntry> items = new ArrayList<>();
    String[] entries = str.split(",", -1);
    for (String entry : entries) {
      if (entry.trim().isEmpty()) { // skip empty
        continue;
      }
      String[] fields = entry.split(":", -1);
      if (fields.length != 2) { // parsing error
        return null;
      }
      items.add(new ParseEntry(fields[0], fields[1]));
    }
    return items;
  }

  @Nullable
  private static S2Polygon internalMakePolygon(String str, boolean normalizeLoops) {
    if (str.equals("empty")) {
      return new S2Polygon(new ArrayList<S2Loop>()); // Can't be an ImmutableList, it is clear()ed.
    }
    List<String> loopStrs = splitString(str, ";");
    List<S2Loop> loops = new ArrayList<>();
    for (String loopStr : loopStrs) {
      S2Loop loop = makeLoop(loopStr);
      if (loop == null) {
        return null;
      }
      // Don't normalize loops that were explicitly specified as "full".
      if (normalizeLoops && !loop.isFull()) {
        loop.normalize();
      }
      loops.add(loop);
    }
    return new S2Polygon(loops);
  }

  private static void appendVertex(S2LatLng ll, StringBuilder out) {
    out.append(Platform.formatDouble((ll.latDegrees())))
        .append(':')
        .append(Platform.formatDouble((ll.lngDegrees())));
  }

  private static void appendVertex(S2Point p, StringBuilder out) {
    appendVertex(new S2LatLng(p), out);
  }

  private static void appendVertices(Iterable<S2Point> points, StringBuilder out) {
    Iterator<S2Point> i = points.iterator();
    while (i.hasNext()) {
      appendVertex(i.next(), out);
      if (i.hasNext()) {
        out.append(", ");
      }
    }
  }
}
