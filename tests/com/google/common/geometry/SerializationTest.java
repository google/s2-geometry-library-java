package com.google.common.geometry;

import com.google.common.annotations.GwtCompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

import junit.framework.TestCase;

import java.util.Comparator;
import java.util.List;

/** Tests both GWT and Java serialization of supported class types. */
@GwtCompatible
public class SerializationTest extends TestCase {

  private void doTest(Object obj) {
    doTest(obj, null);
  }

  private <T> void doTest(T obj, Comparator<T> comparator) {
    try {
      TestPlatform.testSerialization(obj, comparator);
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public void testMutableInteger() {
    doTest(new MutableInteger(5));
  }

  public void testParametrizedS2Point() {
    doTest(new ParametrizedS2Point(0.123, new S2Point(1, 1e7, 1e-9)));
  }

  public void testR1Interval() {
    doTest(new R1Interval(1e5, 1e-6));
  }

  public void testR2Rect() {
    doTest(new R2Rect(new R1Interval(0, 0.1), new R1Interval(1.1, 1.2)));
  }

  public void testR2Vector() {
    doTest(new R2Vector(0.025, 1e6));
  }

  public void testS1Angle() {
    doTest(S1Angle.radians(0.234));
  }

  public void testS1Interval() {
    doTest(new S1Interval(0.5, 2.0));
  }

  public void testS2AreaCentroid() {
    doTest(new S2AreaCentroid(5.0, new S2Point(0.1, 10.3, 7.5)));
  }

  public void testS2Cap() {
    doTest(S2Cap.fromAxisArea(new S2Point(0.1, 0.2, 0.3), 5));
  }

  public void testS2Cell() {
    doTest(new S2Cell(S2LatLng.fromDegrees(0.1, 0.2)));
  }

  public void testS2CellId() {
    doTest(new S2CellId(1234567890123L));
  }

  public void testCellUnion() {
    S2CellUnion union = new S2CellUnion();
    union.initFromIds(Lists.newArrayList(ImmutableList.of(1111L, 2222L, 3333L)));
    doTest(union);
  }

  public void testS2Edge() {
    doTest(new S2Edge(new S2Point(0.1, 0.2, 0.3), new S2Point(0.5, 0.6, -100)));
  }

  public void testS2LatLng() {
    doTest(S2LatLng.fromRadians(0.1, 0.5));
  }

  public void testS2LatLngRect() {
    doTest(new S2LatLngRect(S2LatLng.fromRadians(0.1, 0.2), S2LatLng.fromRadians(1.1, 1.2)));
  }

  public void testS2Loop() {
    S2Loop loop = new S2Loop(ImmutableList.of(
        new S2Point(0.1, 0.2, 0.3),
        new S2Point(0.5, 0, 0),
        new S2Point(4, 4, 4),
        new S2Point(6, 7, 8)));
    doTest(loop);
  }

  public void testS2Point() {
    doTest(new S2Point(0.1, 0.2, 0.3));
  }

  public void testS2Polygon() {
    S2Loop loop1 = new S2Loop(new S2Cell(new S2Point(0.1, 0.2, 0.3)),
        new S2LatLngRect(S2LatLng.CENTER, S2LatLng.fromRadians(1, 2)));
    S2Loop loop2 = new S2Loop(ImmutableList.of(
        new S2Point(0.1, 0.2, 0.3),
        new S2Point(0.5, 0, 0),
        new S2Point(4, 4, 4),
        new S2Point(6, 7, 8)));
    doTest(new S2Polygon(Lists.newArrayList(loop1, loop2)), Ordering.<S2Polygon>natural());
  }

  public void testS2Polyline() {
    List<S2Point> vertices = Lists.newArrayList();
    vertices.add(new S2Point(0, 0, 0));
    vertices.add(new S2Point(0, 0.1, 0));
    vertices.add(new S2Point(0.5, 1e6, 7));
    vertices.add(new S2Point(8, 3, 5));
    doTest(new S2Polyline(vertices));
  }

  public void testS2RegionCoverer() {
    S2RegionCoverer coverer = new S2RegionCoverer();
    coverer.setLevelMod(2);
    coverer.setMaxCells(100);
    coverer.setMinLevel(50);
    coverer.setMaxLevel(500);
    doTest(coverer);
  }
}
