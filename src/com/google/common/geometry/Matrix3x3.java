/*
 * Copyright 2013 Google Inc.
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
import com.google.common.base.Preconditions;

import java.util.List;

import javax.annotation.CheckReturnValue;

/** A simple 3x3 matrix. */
// TODO(eengle): Rename this to Matrix as it is not necessarily 3x3, and make Matrix3x3 a subclass.
@GwtCompatible
public final class Matrix3x3 {
  private final double[] values;
  private final int rows, cols;

  /** Constructs a matrix from a series of column vectors. */
  public static Matrix3x3 fromCols(S2Point ... columns) {
    Matrix3x3 result = new Matrix3x3(3, columns.length);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        result.set(row, col, columns[col].get(row));
      }
    }
    return result;
  }
  
  /** Constructs a matrix from a series of column vectors. */
  public static Matrix3x3 fromCols(List<S2Point> frame) {
    return fromCols(frame.toArray(new S2Point[frame.size()]));
  }

  /** Constructs a 2D matrix of the given width and values. */
  public Matrix3x3(int cols, double ... values) {
    Preconditions.checkArgument(cols >= 0, "Negative rows not allowed.");
    rows = values.length / cols;
    this.cols = cols;
    Preconditions.checkArgument(rows * cols == values.length,
        "Values not an even multiple of 'cols'");
    this.values = values;
  }

  /** Constructs a 2D matrix of a fixed size. */
  public Matrix3x3(int rows, int cols) {
    Preconditions.checkArgument(rows >= 0, "Negative rows not allowed.");
    Preconditions.checkArgument(cols >= 0, "Negative cols not allowed.");
    this.rows = rows;
    this.cols = cols;
    this.values = new double[rows * cols];
  }

  /** Returns the number of rows in this matrix. */
  public int rows() {
    return rows;
  }

  /** Returns the number of columns in this matrix. */
  public int cols() {
    return cols;
  }

  /** Sets a value. */
  public void set(int row, int col, double value) {
    values[row * cols + col] = value;
  }

  /** Gets a value. */
  public double get(int row, int col) {
    return values[row * cols + col];
  }

  /** Returns the transpose of this. */
  @CheckReturnValue
  public Matrix3x3 transpose() {
    Matrix3x3 result = new Matrix3x3(cols, rows);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        result.set(row, col, get(col, row));
      }
    }
    return result;
  }

  /** Returns the result of multiplying this x m. */
  @CheckReturnValue
  public Matrix3x3 mult(Matrix3x3 m) {
    Preconditions.checkArgument(cols == m.rows);
    Matrix3x3 result = new Matrix3x3(rows, m.cols);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        double sum = 0;
        for (int i = 0; i < cols; i++) {
          sum += get(row, i) * m.get(i, col);
        }
        result.set(row, col, sum);
      }
    }
    return result;
  }

  /**
   * Return the vector of the given column.
   */
  public S2Point getCol(int col) {
    Preconditions.checkState(rows == 3);
    Preconditions.checkArgument(0 <= col && col < cols);
    return new S2Point(values[col], values[cols + col], values[2 * cols + col]);
  }

  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Matrix3x3)) {
      return false;
    }
    Matrix3x3 m = (Matrix3x3) o;
    if (rows != m.rows || cols != m.cols) {
      return false;
    }
    for (int i = 0; i < values.length; i++) {
      if (values[i] != m.values[i]) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int hashCode() {
    long hash = 37L * cols;
    for (int i = 0; i < values.length; i++) {
      hash = 37L * hash + Platform.doubleHash(values[i]);
    }
    return (int) hash;
  }
}

