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

import com.google.common.base.Preconditions;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.List;

/** A simple dense matrix class. */
public final class Matrix {
  private final double[] values;
  private final int rows;
  private final int cols;

  /** Constructs a matrix from a series of column vectors. */
  public static Matrix fromCols(S2Point... columns) {
    Matrix result = new Matrix(3, columns.length);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        result.set(row, col, columns[col].get(row));
      }
    }
    return result;
  }

  /** Constructs a matrix from a series of column vectors. */
  public static Matrix fromCols(List<S2Point> frame) {
    return fromCols(frame.toArray(new S2Point[0]));
  }

  /** Constructs a 2D matrix of the given width and values. */
  public Matrix(int cols, double... values) {
    Preconditions.checkArgument(cols >= 0, "Negative cols not allowed.");
    rows = values.length / cols;
    this.cols = cols;
    Preconditions.checkArgument(
        rows * cols == values.length, "Values not an even multiple of 'cols'");
    this.values = values;
  }

  /** Constructs a 2D matrix of a fixed size. */
  public Matrix(int rows, int cols) {
    Preconditions.checkArgument(rows >= 0, "Negative rows not allowed.");
    Preconditions.checkArgument(cols >= 0, "Negative cols not allowed.");
    this.rows = rows;
    this.cols = cols;
    this.values = new double[rows * cols];
  }

  /** Returns the outer product of two S2Point vectors. */
  public static Matrix fromOuter(S2Point ma, S2Point mb) {
    return new Matrix(
        3,
        ma.x * mb.x, ma.x * mb.y, ma.x * mb.z,  // mb.mul(ma.x)
        ma.y * mb.x, ma.y * mb.y, ma.y * mb.z,  // mb.mul(ma.y)
        ma.z * mb.x, ma.z * mb.y, ma.z * mb.z   // mb.mul(ma.z)
      );
  }

  /** Returns the 3x3 identity matrix. */
  public static Matrix identity3x3() {
    return new Matrix(
        3,
        new double[] {
          1, 0, 0, //
          0, 1, 0, //
          0, 0, 1
        });
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
  public Matrix transpose() {
    Matrix result = new Matrix(cols, rows);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        result.set(row, col, get(col, row));
      }
    }
    return result;
  }

  /**
   * Return a matrix that reflects a point across the plane defined by the given normal, which
   * does not need to be unit-length.
   */
  public static Matrix householder(S2Point normal) {
    S2Point unit = normal.normalize();
    return identity3x3().sub(fromOuter(unit, unit).mult(2));
  }

  /**
   * Returns the result of adding the given matrix to this matrix. The two matrices must have equal
   * dimensions.
   */
  @CheckReturnValue
  public Matrix add(Matrix m) {
    Preconditions.checkArgument(rows == m.rows && cols == m.cols);
    Matrix result = new Matrix(rows, cols);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        result.set(row, col, get(row, col) + m.get(row, col));
      }
    }
    return result;
  }

  /**
   * Returns the result of subtracting the given matrix from this matrix. The two matrices must
   * have equal dimensions.
   */
  @CheckReturnValue
  public Matrix sub(Matrix m) {
    Preconditions.checkArgument(rows == m.rows && cols == m.cols);
    Matrix result = new Matrix(rows, cols);
    for (int row = 0; row < result.rows; row++) {
      for (int col = 0; col < result.cols; col++) {
        result.set(row, col, get(row, col) - m.get(row, col));
      }
    }
    return result;
  }

  /** Returns the result of multiplying "this" by the given matrix "m". */
  @CheckReturnValue
  public Matrix mult(Matrix m) {
    Preconditions.checkArgument(cols == m.rows);
    Matrix result = new Matrix(rows, m.cols);
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
   * Returns the result of multiplying "this" 3x3 matrix by the given S2Point column vector "v" to
   * produce an S2Point column vector.
  */
  @CheckReturnValue
  public S2Point mult(S2Point v) {
    Preconditions.checkArgument(rows == 3);
    Preconditions.checkArgument(cols == 3);
    return new S2Point(
        get(0, 0) * v.get(0) + get(0, 1) * v.get(1) + get(0, 2) * v.get(2),
        get(1, 0) * v.get(0) + get(1, 1) * v.get(1) + get(1, 2) * v.get(2),
        get(2, 0) * v.get(0) + get(2, 1) * v.get(1) + get(2, 2) * v.get(2));
  }

  /** Returns the result of multiplying "this" by the given scalar "k". */
  @CheckReturnValue
  public Matrix mult(double k) {
    Matrix result = new Matrix(rows, cols);
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        result.set(row, col, get(row, col) * k);
      }
    }
    return result;
  }

  /** Return the vector of the given column. The matrix must have exactly three rows. */
  public S2Point getCol(int col) {
    Preconditions.checkState(rows == 3);
    Preconditions.checkArgument(0 <= col && col < cols);
    return new S2Point(values[col], values[cols + col], values[2 * cols + col]);
  }

  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Matrix)) {
      return false;
    }
    Matrix m = (Matrix) o;
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

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("Matrix(");
    sb.append(rows);
    sb.append("x");
    sb.append(cols);
    sb.append("): ");
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        sb.append(get(row, col));
        sb.append(" ");
      }
      sb.append("\n");
    }
    return sb.toString();
  }
}
