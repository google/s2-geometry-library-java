package com.google.common.geometry;

import com.google.caliper.Benchmark;

/** 
 * A benchmark for {@link S2LatLng}. Note that return values of benchmark functions are
 * ignored by Caliper.
 */
public class S2LatLngBenchmark {  
  @Benchmark double toPoint(int reps) {
    S2LatLng ll = new S2LatLng(S1Angle.e7(0x150bc888), S1Angle.e7(0x5099d63f));
    double xSum = 0;
    for (int r = reps; r > 0; --r) {
      xSum += ll.toPoint().getX();
    }
    return xSum;
  }
  
  @Benchmark double getDistance(int reps) {
    double sum = 0;
    S2LatLng x = S2LatLng.fromDegrees(25.0,  -78.0);
    for (int r = reps; r > 0; --r) {
      S2LatLng y = new S2LatLng(S1Angle.e7(r), S1Angle.degrees(56.0));
      sum += x.getDistance(y).radians();
    }
    return sum;
  }
}