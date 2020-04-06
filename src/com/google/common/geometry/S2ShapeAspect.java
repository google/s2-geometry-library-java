/*
 * Copyright 2019 Google LLC.
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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.AbstractList;
import java.util.Arrays;
import java.util.List;

/**
 * A set of partial {@link S2Shape shape} implementations, effectively breaking down the S2Shape API
 * into several aspects, each focused on a subset of the overall API:
 * 
 * <ul>
 * <li>{@link VertexAspect} provides a logical list of vertices, where the 'vertexId' is
 * at least 0 and less than {@link VertexAspect#numVertices}. Each implementation stores the list in
 * a different way, for example a {@link ChainAspect.Simple.Packed packed array}. This isn't part of
 * the S2Shape API, but is provided for use by the other aspects.
 * <li>{@link EdgeAspect} provides the 'vertexId' that starts and ends each edge or each
 * chain/offset, where the 'edgeId' is at least 0 and less than {@link EdgeAspect#numEdges}, the
 * 'chainId' is at least 0 and less than {@link ChainAspect#numChains}, and the 'edgeOffset' is at
 * least {@code edgeId(chainId)} and less than {@code edgeId(chainid+1)}. For example, the endpoint
 * of the last {@link EdgeAspect.Closed closed} edge wraps back to the first vertex of that chain.
 * <li>{@link ChainAspect} provides a mapping between chains and edge ranges, where 'chainId' is at
 * least 0 and less than {@link ChainAspect#numChains}, and the {@link ChainAspect#getChainStart}
 * and {@link ChainAspect#getChainLength} methods provide the 'edgeId' range of each chain.
 * <li>{@link TopoAspect} provides the methods to relate a point in the world to the interior,
 * exterior, or boundary of the shape.
 * 
 * <p>There may be fewer edges than vertices, e.g. 2 vertices can define 1 edge.
 */
@GwtIncompatible("Insufficient support for generics")
interface S2ShapeAspect {
  /** A provider of S2Point given a 'vertexId', allowing alternate storage options. */
  interface VertexAspect {
    /** Returns the number of vertices. May be different from {@link S2Shape#numEdges}. */
    int numVertices();
    
    /** Returns a vertex of this shape, from 0 (inclusive) to {@link #numVertices} (exclusive). */
    S2Point vertex(int vertexId);
    
    /** Returns the vertices in this shape. Less efficient but may be more convenient. */
    default List<S2Point> vertices() {
      return new AbstractList<S2Point>() {
        @Override public int size() {
          return numVertices();
        }
        @Override public S2Point get(int index) {
          return vertex(index);
        }
      };
    }
  }

  /**
   * A provider of the 'vertexId' for the start and end of each 'edgeId' or 'chainId'/'edgeOffset',
   * allowing alternate edge/vertex mappings.
   */
  interface EdgeAspect {
    /**
     * Returns the vertexId that starts 'edgeId', assuming
     * {@code edgeId(chainId) <= edgeId && edgeId < edgeId(chainId + 1)}.
     */
    int vertexId(int chainId, int edgeId);
    
    /**
     * Converts the given array of 'vertexId' values in place, yielding an array of 'edgeId' values
     * that start each chain. This requires knowledge of the edge/vertex mapping, and hence this
     * aspect of S2Shape construction is delegated here.
     */
    void adjustChains(int ... chainStarts);
    
    /**
     * Returns the start point of the edge that would be returned by {@link S2Shape#getChainEdge},
     * or the endpoint of the last edge if {@code edgeOffset==getChainLength(chainId)}.
     */
    S2Point getChainVertex(int chainId, int edgeOffset);
    
    /** Provides {@link S2Shape#numEdges}. */
    int numEdges();
    
    /** Provides {@link S2Shape#getEdge}. */
    void getEdge(int edgeId, MutableEdge result);
    
    /** Provides {@link S2Shape#getChainEdge}. */
    void getChainEdge(int chainId, int edgeOffset, MutableEdge result);
    
    /** Chains are closed, that is, there is an implicit edge between the ends of each chain. */
    interface Closed extends Mixed {
      @Override default void adjustChains(int ... chainStarts) {
      }
      
      @Override default int numEdges() {
        return numVertices();
      }
      
      @Override default void getEdge(int edgeId, MutableEdge result) {
        // Note edgeId=vertexId, since the last edge is implicit.
        result.set(vertex(edgeId), vertex(vertexId(chainId(edgeId), edgeId + 1)));
      }
      
      @Override default void getChainEdge(int chainId, int edgeOffset, MutableEdge result) {
        int edgeId = getChainStart(chainId) + edgeOffset;
        result.set(vertex(edgeId), vertex(vertexId(chainId, edgeId + 1)));
      }
      
      @Override default S2Point getChainVertex(int chainId, int edgeOffset) {
        return vertex(vertexId(chainId, getChainStart(chainId) + edgeOffset));
      }
      
      @Override default int vertexId(int chainId, int edgeId) {
        return edgeId < edgeId(chainId + 1) ? edgeId : getChainStart(chainId);
      }
    }
    
    /** Chains are open, that is, there is no implicit edge between the ends of each chain. */
    interface Open extends Mixed {
      @Override default void adjustChains(int ... chainStarts) {
        Preconditions.checkArgument(chainStarts.length > 0, "Must have at least 1 chain.");
        int last = chainStarts[0];
        for (int i = 1; i < chainStarts.length; i++) {
          int offset = chainStarts[i];
          chainStarts[i] -= i;
          Preconditions.checkArgument(last != offset, "Must have at least 1 edge.");
          last = offset;
        }
      }
      
      @Override default int numEdges() {
        return numVertices() - numChains();
      }
      
      @Override default void getEdge(int edgeId, MutableEdge result) {
        int vertexId = vertexId(chainId(edgeId), edgeId);
        result.set(vertex(vertexId), vertex(vertexId + 1));
      }
      
      @Override default void getChainEdge(int chainId, int edgeOffset, MutableEdge result) {
        int vertexId = vertexId(chainId, getChainStart(chainId) + edgeOffset);
        result.set(vertex(vertexId), vertex(vertexId + 1));
      }
      
      @Override default S2Point getChainVertex(int chainId, int edgeOffset) {
        return vertex(vertexId(chainId, getChainStart(chainId) + edgeOffset));
      }
      
      @Override default int vertexId(int chainId, int edgeId) {
        return chainId + edgeId;
      }
    }
  }
  
  /** A provider of the 'edgeId' ranges for each chain, allowing alternate chain representations. */
  interface ChainAspect {
    /** Returns the chain ID of a given edge. */
    int chainId(int edgeId);
    
    /** Returns start edge ID of a chain, or the number of edges if {@code chainId==numChains()}. */
    int edgeId(int chainId);

    /** Provides {@link S2Shape#numChains}. */
    int numChains();
    
    /** Provides {@link S2Shape#getChainStart}. */
    int getChainStart(int chainId);
    
    /** Provides {@link S2Shape#getChainLength}. */
    int getChainLength(int chainId);
    
    /** A single non-empty chain. */
    abstract class Simple implements Mixed {
      @Override public int numChains() {
        return 1;
      }
      
      @Override public int getChainStart(int chainId) {
        Preconditions.checkElementIndex(chainId, 1);
        return 0;
      }
      
      @Override public int getChainLength(int chainId) {
        Preconditions.checkElementIndex(chainId, 1);
        return numEdges();
      }
      
      @Override public int edgeId(int chainId) {
        switch (chainId) {
          case 0:
            return 0;
          case 1:
            return numEdges();
          default:
            throw new IndexOutOfBoundsException("Invalid chain " + chainId);
        }
      }
      
      @Override public int chainId(int edgeIndex) {
        Preconditions.checkElementIndex(edgeIndex, numEdges());
        return 0;
      }

      /** A simple chain of S2Point references. */
      abstract static class Array extends Simple {
        private final S2Point[] vertices;
        
        Array(Iterable<S2Point> vertices) {
          this.vertices = toArray(vertices);
        }
        
        @Override public int numVertices() {
          return vertices.length;
        }
        
        @Override public S2Point vertex(int index) {
          return vertices[index];
        }
        
        /** Returns an array of the given vertices. */
        // Note this implementation overcomes lack of GWT support for Iterables.toArray.
        private static S2Point[] toArray(Iterable<S2Point> vertices) {
          S2Point[] array = new S2Point[Iterables.size(vertices)];
          int offset = 0;
          for (S2Point v : vertices) {
            array[offset++] = v;
          }
          return array;
        }
      }
      
      /** A simple chain of packed coordinates. */
      abstract static class Packed extends Simple {
        private final double[] coordinates;
        
        public Packed(Iterable<S2Point> vertices) {
          this.coordinates = toArray(vertices);
        }
        
        @Override public int numVertices() {
          return coordinates.length / 3;
        }
        
        @Override public S2Point vertex(int index) {
          return vertex(coordinates, index);
        }
        
        private static double[] toArray(Iterable<S2Point> vertices) {
          double[] coordinates = new double[3 * Iterables.size(vertices)];
          int offset = 0;
          for (S2Point v : vertices) {
            coordinates[offset++] = v.x;
            coordinates[offset++] = v.y;
            coordinates[offset++] = v.z;
          }
          return coordinates;
        }
        
        private static S2Point vertex(double[] coordinates, int index) {
          int offset = 3 * index;
          return new S2Point(coordinates[offset], coordinates[offset + 1], coordinates[offset + 2]);
        }
      }
      
      /** A simple chain of packed cell centers. */
      abstract static class Snapped extends Simple {
        private final long[] vertices;
        
        public Snapped(Iterable<S2CellId> vertices) {
          this.vertices = toArray(vertices);
        }
        
        @Override public int numVertices() {
          return vertices.length;
        }
        
        @Override public S2Point vertex(int index) {
          return new S2CellId(vertices[index]).toPoint();
        }
        
        private static long[] toArray(Iterable<S2CellId> vertices) {
          long[] ids = new long[Iterables.size(vertices)];
          int offset = 0;
          for (S2CellId vertex : vertices) {
            ids[offset++] = vertex.id();
          }
          return ids;
        }
      }
    }
    
    /** A sequence of chains, represented as an array of the first 'edgeId' for each chain. */
    abstract class Multi implements Mixed {
      private final int[] cumulativeEdges;
      
      public Multi(Iterable<? extends Iterable<?>> chains) {
        this.cumulativeEdges = new int[Iterables.size(chains) + 1];
        int sum = 0;
        int offset = 0;
        for (Iterable<?> chain : chains) {
          cumulativeEdges[offset++] = sum;
          sum += Iterables.size(chain);
        }
        cumulativeEdges[offset] = sum;
        adjustChains(cumulativeEdges);
      }

      Multi(int[] cumulativeEdges) {
        this.cumulativeEdges = cumulativeEdges;
        adjustChains(cumulativeEdges);
      }

      @Override public final int numChains() {
        return cumulativeEdges.length - 1;
      }
      
      @Override public final int edgeId(int chainId) {
        return cumulativeEdges[chainId];
      }
      
      @Override public final int getChainStart(int chainId) {
        Preconditions.checkElementIndex(chainId, numChains());
        return edgeId(chainId);
      }
      
      @Override public final int getChainLength(int chainId) {
        return edgeId(chainId + 1) - edgeId(chainId);
      }
      
      @Override public final int chainId(int edgeId) {
        int chainId = Arrays.binarySearch(cumulativeEdges, edgeId);
        if (chainId < 0) {
          chainId = -chainId - 2;
        }
        // The binary search may have landed on an empty chain, which cannot match 'edgeId'.
        while (getChainLength(chainId) == 0) {
          chainId++;
        }
        return chainId;
      }

      /** An array of S2Point references for multiple chains. */
      abstract static class Array extends Multi {
        private final S2Point[] vertices;
        
        Array(Iterable<? extends Iterable<S2Point>> chains) {
          super(chains);
          this.vertices = Simple.Array.toArray(Iterables.concat(chains));
        }
        
        @Override public int numVertices() {
          return vertices.length;
        }
        
        @Override public S2Point vertex(int index) {
          return vertices[index];
        }
      }
      
      /** Packed coordinates for multiple chains. */
      abstract static class Packed extends Multi {
        private final double[] coordinates;
        
        public Packed(Iterable<? extends Iterable<S2Point>> chains) {
          super(chains);
          this.coordinates = Simple.Packed.toArray(Iterables.concat(chains));
        }
        
        @Override public int numVertices() {
          return coordinates.length / 3;
        }
        
        @Override public S2Point vertex(int index) {
          return Simple.Packed.vertex(coordinates, index);
        }
      }
      
      /** Snapped cell centers for multiple chains. */
      abstract static class Snapped extends Multi {
        private final long[] vertices;
        
        public Snapped(Iterable<? extends Iterable<S2CellId>> chains) {
          super(chains);
          this.vertices = Simple.Snapped.toArray(Iterables.concat(chains));
        }
        
        @Override public int numVertices() {
          return vertices.length;
        }
        
        @Override public S2Point vertex(int index) {
          return new S2CellId(vertices[index]).toPoint();
        }
      }
    }
  }
  
  /** How world positions are classified as exterior, interior, or on the boundary of the object. */
  interface TopoAspect {
    /** Provides {@link S2Shape#hasInterior}. */
    boolean hasInterior();

    /** Provides {@link S2Shape#containsOrigin}. */
    boolean containsOrigin();

    /** Provides {@link S2Shape#dimension}. */
    int dimension();
  }
  
  /** A full S2Shape that mixes together each aspect. */
  interface Mixed extends S2Shape, VertexAspect, EdgeAspect, ChainAspect, TopoAspect {}
}
