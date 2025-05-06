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
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.jspecify.annotations.Nullable;

/**
 * S2TextFormat contains a collection of functions for converting geometry to and from a human-
 * readable format. It is mainly intended for testing and debugging. Be aware that the human-
 * readable format is *not* designed to preserve the full precision of the original object, so it
 * should not be used for data storage.
 */
public final class S2TextFormat {

  /**
   * Returns an S2Point corresponding to the given a latitude-longitude coordinate in degrees.
   * Example of the input format: "-20:150".
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Point makePointOrDie(String str) throws IllegalArgumentException {
    S2Point point = makePoint(str);
    Preconditions.checkArgument(point != null, ": str == \"%s\"", str);
    return point;
  }

  /**
   * As {@link #makePointOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2Point makePoint(String str) {
    List<S2Point> vertices = parsePoints(str);
    if (vertices == null || vertices.size() != 1) {
      return null;
    }
    return vertices.get(0);
  }

  /**
   * Parses a string of one or more comma-separated latitude-longitude coordinates in degrees, and
   * returns the corresponding List of S2LatLng points.
   *
   * <pre> Examples of the input format:
   *     ""                          // no points
   *     "-20:150"                   // one point
   *     "-20:150, -20:151, -19:150" // three points
   * </pre>
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static List<S2LatLng> parseLatLngsOrDie(String str) throws IllegalArgumentException {
    List<S2LatLng> latlngs = parseLatLngs(str);
    Preconditions.checkArgument(latlngs != null, ": str == \"%s\"", str);
    return latlngs;
  }

  /**
   * As {@link #parseLatLngsOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  static @Nullable List<S2LatLng> parseLatLngs(String str) {
    ImmutableList<ParseEntry> ps = dictionaryParse(str);
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
   * Parses multiple points from a string, where the points are separated by '|', using {@link
   * #makePointOrDie(String)} to parse each point. Note that this "|"-delimited format is different
   * than {@link #parsePointsOrDie(String)}, which uses a comma-delimited format.
   */
  public static List<S2Point> makePointsOrDie(String str) {
    ArrayList<S2Point> result = new ArrayList<>();
    for (String lineStr : splitString(str, "\\|")) {
      result.add(makePointOrDie(lineStr));
    }
    return result;
  }

  /**
   * Parses a string in the same format as {@link #parseLatLngs(String)}, and returns the
   * corresponding List of S2Point values.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static List<S2Point> parsePointsOrDie(String str) throws IllegalArgumentException {
    List<S2Point> vertices = parsePoints(str);
    Preconditions.checkArgument(vertices != null, ": str == \"%s\"", str);
    return vertices;
  }

  /**
   * As {@link #parsePointsOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable List<S2Point> parsePoints(String str) {
    return parsePoints(str, -1);
  }

  /**
   * As {@link #parsePoints(String)} above, but also snaps the vertices to the specified S2CellId
   * level if level >= 0. Returns null if conversion is unsuccessful.
   */
  public static @Nullable List<S2Point> parsePoints(String str, int level) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    if (latlngs == null) {
      return null;
    }
    List<S2Point> vertices = new ArrayList<>();
    for (S2LatLng latlng : latlngs) {
      S2Point vertex = latlng.toPoint();
      if (level >= 0) {
        vertex = snapPointToLevel(vertex, level);
      }
      vertices.add(vertex);
    }
    return vertices;
  }

  /**
   * Parses the given String into S2Points and appends them to the provided list of vertices. The
   * provided string must be a comma-separated sequence of lat:lng values in degrees, where each lat
   * and lng can be parsed as a Double, like "10:20, 90:0, 20:30". Returns the bounding rectangle of
   * the parsed points.
   */
  @CanIgnoreReturnValue
  public static S2LatLngRect parseVertices(String str, List<S2Point> vertices) {
    List<S2Point> points = parsePoints(str);
    S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();
    for (S2Point point : points) {
      vertices.add(point);
      builder.addPoint(point);
    }
    return builder.build();
  }

  /** Snaps the given S2Point to the nearest S2CellId center at the given level. */
  public static S2Point snapPointToLevel(S2Point point, int level) {
    return S2CellId.fromPoint(point).parent(level).toPoint();
  }

  /** Snaps the given list of S2Points to the nearest S2CellId centers at the given level. */
  public static List<S2Point> snapPointsToLevel(List<S2Point> points, int level) {
    List<S2Point> snapped = new ArrayList<>(points.size());
    for (S2Point p : points) {
      snapped.add(snapPointToLevel(p, level));
    }
    return snapped;
  }

  /**
   * Given a string in the same format as {@link #parseLatLngs(String)}, returns a single S2LatLng.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2LatLng makeLatLngOrDie(String str) throws IllegalArgumentException {
    S2LatLng latlng = makeLatLng(str);
    Preconditions.checkArgument(latlng != null, ": str == \"%s\"", str);
    return latlng;
  }

  /**
   * As {@link #makeLatLngOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2LatLng makeLatLng(String str) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    if (latlngs == null || latlngs.size() != 1) {
      return null;
    }
    return latlngs.get(0);
  }

  /**
   * Given a string in the same format as {@link #parseLatLngs(String)}, returns the minimal
   * bounding S2LatLngRect that contains the coordinates.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2LatLngRect makeLatLngRectOrDie(String str) throws IllegalArgumentException {
    S2LatLngRect rect = makeLatLngRect(str);
    Preconditions.checkArgument(rect != null, ": str == \"%s\"", str);
    return rect;
  }

  /**
   * As {@link #makeLatLngRectOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2LatLngRect makeLatLngRect(String str) {
    List<S2LatLng> latlngs = parseLatLngs(str);
    if (latlngs == null || latlngs.isEmpty()) {
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
   * <p>This function is a wrapper for {@link S2CellId#fromDebugString(String)}.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2CellId makeCellIdOrDie(String str) throws IllegalArgumentException {
    S2CellId cellId = makeCellId(str);
    Preconditions.checkArgument(cellId != null, ": str == \"%s\"", str);
    return cellId;
  }

  /**
   * As {@link #makeCellIdOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2CellId makeCellId(String str) {
    S2CellId cellId = S2CellId.fromDebugString(str);
    if (cellId.equals(S2CellId.none())) {
      return null;
    }
    return cellId;
  }

  /**
   * Parses a comma-separated list of S2CellIds as described in {@link #makeCellIdOrDie(String)},
   * and returns the corresponding S2CellUnion. (Note that S2CellUnions are automatically normalized
   * by sorting, removing duplicates, and replacing groups of 4 child cells by their parent cell.)
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2CellUnion makeCellUnionOrDie(String str) throws IllegalArgumentException {
    S2CellUnion cellUnion = makeCellUnion(str);
    Preconditions.checkArgument(cellUnion != null, ": str == \"%s\"", str);
    return cellUnion;
  }

  /**
   * As {@link #makeCellUnionOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2CellUnion makeCellUnion(String str) {
    ArrayList<S2CellId> cellIds = new ArrayList<>();
    for (String cellStr : splitString(str, ",")) {
      cellStr = cellStr.trim();
      S2CellId cellId = makeCellId(cellStr);
      if (cellId == null) {
        return null;
      }
      cellIds.add(cellId);
    }
    return new S2CellUnion().initFromCellIds(cellIds);
  }

  /**
   * Given a string of latitude-longitude coordinates in degrees, returns a newly allocated loop.
   * Example of the input format:
   *
   * <p>"-20:150, 10:-120, 0.123:-170.652"
   *
   * <p>The strings "empty" or "full" create an empty or full loop respectively.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Loop makeLoopOrDie(String str) throws IllegalArgumentException {
    return makeLoopOrDie(str, -1);
  }

  /**
   * Like {@link #makeLoopOrDie(String)} above, but also snaps the vertices of the resultant loop to
   * the given S2CellId level.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Loop makeLoopOrDie(String str, int level) throws IllegalArgumentException {
    S2Loop loop = makeLoop(str, level);
    Preconditions.checkArgument(loop != null, ": str == \"%s\"", str);
    return loop;
  }

  /**
   * Like {@link #makeLoopOrDie(String)} above, but does not throw IllegalArgumentException on
   * failure, instead returning null.
   */
  public static @Nullable S2Loop makeLoop(String str) {
    return makeLoop(str, -1);
  }

  /**
   * Like {@link #makeLoopOrDie(String, int)} above, but does not throw IllegalArgumentException on
   * failure, instead returning null.
   */
  public static @Nullable S2Loop makeLoop(String str, int level) {
    if ("empty".equalsIgnoreCase(str)) {
      return S2Loop.empty();
    } else if ("full".equalsIgnoreCase(str)) {
      return S2Loop.full();
    }
    List<S2Point> vertices = parsePoints(str, level);
    if (vertices == null) {
      return null;
    }
    return new S2Loop(vertices);
  }

  /** Wraps the S2Loop constructor to accept an array of vertices. */
  public static S2Loop makeLoop(S2Point... vertices) {
    return new S2Loop(vertices);
  }

  /**
   * Parses multiple polylines from a string, where the polylines are separated by '|', using {@link
   * #makePolylineOrDie(String)} to parse each polyline.
   */
  public static List<S2Polyline> makePolylinesOrDie(String str) {
    ArrayList<S2Polyline> result = new ArrayList<>();
    for (String lineStr : splitString(str, "\\|")) {
      result.add(makePolylineOrDie(lineStr));
    }
    return result;
  }

  /**
   * Parses an input in the same format as {@link #makeLoop(String)}, but returns an S2Polyline
   * rather than an S2Loop.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Polyline makePolylineOrDie(String str) throws IllegalArgumentException {
    S2Polyline polyline = makePolyline(str);
    Preconditions.checkArgument(polyline != null, ": str == \"%s\"", str);
    return polyline;
  }

  /**
   * As {@link #makePolylineOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2Polyline makePolyline(String str) {
    List<S2Point> vertices = parsePoints(str);
    if (vertices == null) {
      return null;
    }
    return new S2Polyline(vertices);
  }

  /**
   * Given a string of lat,lng points in degrees, returns an S2Polyline containing those points
   * after passing it through the encode/decode roundtrip.
   */
  public static S2Polyline makePolylineWithEncodeDecodeRoundtrip(String str) {
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    try {
      ByteArrayOutputStream bos = new ByteArrayOutputStream();
      new S2Polyline(vertices).encode(bos);
      return S2Polyline.decode(new ByteArrayInputStream(bos.toByteArray()));
    } catch (IOException e) {
      throw new IllegalArgumentException(e);
    }
  }

  /**
   * Like {@link #makePolyline(String)} above, but returns an S2LaxPolylineShape instead.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2LaxPolylineShape makeLaxPolylineOrDie(String str)
      throws IllegalArgumentException {
    S2LaxPolylineShape laxPolyline = makeLaxPolyline(str);
    Preconditions.checkArgument(laxPolyline != null, ": str == \"%s\"", str);
    return laxPolyline;
  }

  /**
   * As {@link #makeLaxPolylineOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2LaxPolylineShape makeLaxPolyline(String str) {
    List<S2Point> vertices = parsePoints(str);
    if (vertices == null) {
      return null;
    }
    return S2LaxPolylineShape.create(vertices);
  }

  /**
   * Parses multiple polygons from a string, where the polygons are separated by '|', using {@link
   * #makePolygonOrDie(String)} to parse each polygon.
   */
  public static List<S2Polygon> makePolygonsOrDie(String str) {
    ArrayList<S2Polygon> result = new ArrayList<>();
    for (String polyStr : splitString(str, "\\|")) {
      result.add(makePolygonOrDie(polyStr));
    }
    return result;
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
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Polygon makePolygonOrDie(String str) throws IllegalArgumentException {
    S2Polygon polygon = makePolygon(str);
    Preconditions.checkArgument(polygon != null, ": str == \"%s\"", str);
    return polygon;
  }

  /**
   * As {@link #makePolygonOrDie(String)} above, but also snaps the polygon to the specified level.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Polygon makePolygonOrDie(String str, int level) throws IllegalArgumentException {
    S2Polygon polygon = makePolygon(str, level);
    Preconditions.checkArgument(polygon != null, ": str == \"%s\"", str);
    return polygon;
  }

  /**
   * As {@link #makePolygonOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2Polygon makePolygon(String str) {
    return makePolygon(str, -1);
  }

  /**
   * As {@link #makePolygonOrDie(String, int)}, above but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2Polygon makePolygon(String str, int level) {
    return internalMakePolygon(str, true, level);
  }

  /**
   * Parses a string in the same format as {@link #makePolygonOrDie(String)}, except that loops must
   * be oriented so that the interior of the loop is always on the left, and polygons with
   * degeneracies are supported. As with {@link #makePolygon(String)}, "full" denotes the full
   * polygon, and "" or "empty" denote the empty polygon.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2LaxPolygonShape makeLaxPolygonOrDie(String str) throws IllegalArgumentException {
    S2LaxPolygonShape laxPolygon = makeLaxPolygon(str);
    Preconditions.checkArgument(laxPolygon != null, ": str == \"%s\"", str);
    return laxPolygon;
  }

  /**
   * As {@link #makeLaxPolygonOrDie(String)} above, but does not throw IllegalArgumentException on
   * invalid input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2LaxPolygonShape makeLaxPolygon(String str) {
    List<String> loopStrs = splitString(str, ";");
    List<List<S2Point>> loops = new ArrayList<>();
    for (String loopStr : loopStrs) {
      if (loopStr.equals("full")) {
        loops.add(new ArrayList<S2Point>());
      } else if (!loopStr.equals("empty")) {
        List<S2Point> points = parsePoints(loopStr);
        if (points == null) {
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
   *   <li>1:2 | 2:3 # # // Two points (one S2Point.Shape)
   *   <li># 0:0, 1:1, 2:2 | 3:3, 4:4 # // Two polylines
   *   <li># # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2 // Two nested loops (one polygon)
   *   <li>5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0 // One of each point, line, and polygon
   *   <li># # empty // One empty polygon
   *   <li># # empty | full // One empty polygon, one full polygon
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
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2ShapeIndex makeIndexOrDie(String str) throws IllegalArgumentException {
    S2ShapeIndex index = makeIndex(str);
    Preconditions.checkArgument(index != null, ": str == \"%s\"", str);
    return index;
  }

  /**
   * As {@link #makeIndexOrDie(String)} above, but does not throw IllegalStatException on invalid
   * input. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2ShapeIndex makeIndex(String str) {
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

  /**
   * Like {@link #makeIndex(String)} above, but returns an index containing S2Polylines instead of
   * S2LaxPolylines, and S2Polygon.Shapes instead of S2LaxPolygonShapes.
   *
   * <p>Returns an S2ShapeIndex containing the points, polylines, and loops (in the form of a single
   * polygon) described by the following format:
   *
   * <pre>point1|point2|... # line1|line2|... # polygon1|polygon2|...</pre>
   *
   * <table>
   * <caption>Examples</caption>
   * <tr><td>Two points</td><td><pre>1:2 | 2:3 # #</pre></td></tr>
   * <tr><td>Two polylines</td><td><pre># 0:0, 1:1, 2:2 | 3:3, 4:4 #</pre></td></tr>
   * <tr><td>Shell with hole</td><td><pre># # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2</pre></td></tr>
   * <tr><td>One of each</td><td><pre>5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0</pre></td></tr>
   * </table>
   *
   * <p>Loops should be directed so that the region's interior is on the left.
   */
  public static S2ShapeIndex makeIndexWithLegacyShapes(String str) {
    S2ShapeIndex index = new S2ShapeIndex();
    ImmutableList<String> strs = ImmutableList.copyOf(Splitter.on('#').split(str));
    Preconditions.checkArgument(strs.size() == 3, "Must contain two # characters: %s", str);
    List<S2Point> points = new ArrayList<>();
    for (String point : Splitter.on('|').omitEmptyStrings().split(strs.get(0).trim())) {
      points.add(makePoint(point));
    }
    if (!points.isEmpty()) {
      index.add(S2Point.Shape.fromList(points));
    }
    for (String line : Splitter.on('|').omitEmptyStrings().split(strs.get(1).trim())) {
      index.add(makePolyline(line));
    }
    for (String polygon : Splitter.on('|').omitEmptyStrings().split(strs.get(2).trim())) {
      index.add(makePolygon(polygon).shape());
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

  /**
   * Convert a list of S2Points to the S2TextFormat string representation documented above as if it
   * were a line.
   */
  public static String toString(List<S2Point> linePoints) {
    StringBuilder out = new StringBuilder();
    if (!linePoints.isEmpty()) {
      appendVertices(linePoints, out);
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
   * Convert an S2ShapeIndex to the S2TextFormat string representation documented above. The index
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

  /**
   * Returns the current content of the given S2BuilderGraph as a multi-line string, with no limits
   * on size. Intended for debugging unit tests.
   */
  public static String toString(S2BuilderGraph graph) {
    StringBuilder sb = new StringBuilder();
    sb.append("Vertices: \n");
    for (int vertexId = 0; vertexId < graph.numVertices(); ++vertexId) {
      sb.append("  VertexId ")
          .append(vertexId)
          .append(" is point ")
          .append(graph.vertex(vertexId))
          .append("\n");
    }
    sb.append("Edges: \n");
    for (int edgeId = 0; edgeId < graph.numEdges(); ++edgeId) {
      int setId = graph.inputEdgeIdSetIds().get(edgeId);
      sb.append("  EdgeId ")
          .append(edgeId)
          .append(" is edge (")
          .append(graph.edgeSrcId(edgeId))
          .append(" - ")
          .append(graph.edgeDstId(edgeId))
          .append(")\n")
          .append("      with input edge ids: ")
          .append(graph.inputEdgeIdSetLexicon().idSet(setId))
          .append("\n");
    }
    return sb.toString();
  }

  /**
   * Summarizes some key information about the given index: for each shape, the dimension, the
   * number of edges and chains, and for each chain the chainStart and length. Useful for debugging.
   */
  public static String summarizeIndex(S2ShapeIndex index) {
    StringBuilder sb = new StringBuilder();
    sb.append("Index has").append(index.getShapes().size()).append(" shapes.\n");
    for (int shapeId = 0; shapeId < index.getShapes().size(); ++shapeId) {
      S2Shape shape = index.getShapes().get(shapeId);
      sb.append("  Shape #")
          .append(shapeId)
          .append(" has dimension ")
          .append(shape.dimension())
          .append(", ")
          .append(shape.numEdges())
          .append(" edges, and ")
          .append(shape.numChains())
          .append(" chains.")
          .append("\n");
      for (int c = 0; c < shape.numChains(); ++c) {
        sb.append("    Chain #")
            .append(c)
            .append(" starts at shape edge #")
            .append(shape.getChainStart(c))
            .append(" and has length ")
            .append(shape.getChainLength(c))
            .append("\n");
      }
    }
    return sb.toString();
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

  // Split on the given regexp. Trim whitespace and skip empty strings to produce the result.
  public static List<String> splitString(String str, String regexp) {
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
  private static @Nullable ImmutableList<ParseEntry> dictionaryParse(String str) {
    if (str == null || str.isEmpty()) {
      return ImmutableList.of();
    }
    ImmutableList.Builder<ParseEntry> items = new ImmutableList.Builder<>();
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
    return items.build();
  }

  static @Nullable S2Polygon internalMakePolygon(String str, boolean normalizeLoops, int level) {
    if (str.equals("empty")) {
      return new S2Polygon();
    }
    List<String> loopStrs = splitString(str, ";");
    List<S2Loop> loops = new ArrayList<>();
    for (String loopStr : loopStrs) {
      S2Loop loop = makeLoop(loopStr, level);
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

  /** Formats the given S2LatLng as a colon-separated pair of values, like "12.345:-67.89". */
  private static void appendVertex(S2LatLng ll, StringBuilder out) {
    out.append(Platform.formatDouble(ll.latDegrees()))
        .append(':')
        .append(Platform.formatDouble(ll.lngDegrees()));
  }

  /**
   * Formats the given S2Point as an S2LatLng using {@link #appendVertex(S2LatLng, StringBuilder)}.
   */
  private static void appendVertex(S2Point p, StringBuilder out) {
    if (p == null) {
      out.append("null");
      return;
    }
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

  private S2TextFormat() {}
}
