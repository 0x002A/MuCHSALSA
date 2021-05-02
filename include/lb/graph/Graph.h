// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of LazyBastardOnMate <https://github.com/0x002A/LazyBastardOnMate>.
//
// LazyBastardOnMate is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// LazyBastardOnMate is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with LazyBastardOnMate.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#ifndef INCLUDED_LAZYBASTARD_GRAPH
#define INCLUDED_LAZYBASTARD_GRAPH

#pragma once

#include <algorithm>
#include <cstddef>
#include <deque>
#include <gsl/pointers>
#include <iterator>
#include <list>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Lb.fwd.h"
#include "graph/Vertex.h"

namespace lazybastard {

// =====================================================================================================================
//                                                   UTILITY FUNCTIONS
// =====================================================================================================================

// ----------------
// struct GraphUtil
// ----------------

/**
 * Struct providing a namespace for utility functions regarding Graph and DiGraph.
 */
struct GraphUtil {
  /**
   * Calculates the shortest path between two Vertex instances within a given Graph or DiGraph.
   *
   * @param pGraph a pointer pointing to a Graph or DiGraph
   * @param pVertices a std::pair containing pointers pointing to the start and end Vertex
   * @return A std::deque containing pointers pointing to the Vertex instances representing the shortest path between
   *         the supplied Vertex instances
   */
  template <typename TYPE, typename std::enable_if_t<std::is_base_of<graph::GraphBase, TYPE>::value, bool> = true>
  static auto getShortestPath(
      gsl::not_null<TYPE const *>                                                                       pGraph,
      std::pair<gsl::not_null<graph::Vertex const *> const, gsl::not_null<graph::Vertex const *> const> pVertices);

private:
  /**
   * Returns all Vertex instances connected to the supplied Vertex.
   *
   * @param graph a const reference to the Graph containing the Vertex
   * @param pVertex a pointer pointing to the Vertex to return the reachable Vertex instances for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  [[maybe_unused]] static auto _getReachableVertices(graph::Graph const &                              graph,
                                                     gsl::not_null<lazybastard::graph::Vertex const *> pVertex);

  /**
   * Returns all Vertex instances connected to the supplied Vertex.
   *
   * @param diGraph a const reference to the DiGraph containing the Vertex
   * @param pVertex a pointer pointing to the Vertex to return the reachable Vertex instances for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  [[maybe_unused]] static auto _getReachableVertices(graph::DiGraph const &                            diGraph,
                                                     gsl::not_null<lazybastard::graph::Vertex const *> pVertex);
};

namespace graph {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ---------------
// class GraphBase
// ---------------

/**
 * Base class implementing common functionality required by concrete Graph implementations (Graph and DiGraph).
 *
 * In general a Graph holds a list of Vertex and Edge instances.
 * The methods of this class are **partially** thread-safe.
 */
class GraphBase {
public:
  /**
   * Destructor.
   */
  virtual ~GraphBase() = default;

  /**
   * Returns a bool indicating whether the supplied Vertex is present or not.
   *
   * @param nanoporeId a const reference to a std::string representing the id of the Vertex to be checked
   * @return A bool indicating whether the Vertex is present or not
   */
  bool hasVertex(std::string const &nanoporeId) const;

  /**
   * Returns a std::shared_ptr to the requested Vertex instance if found.
   * If no result is found, the std::shared_ptr will be initialized with nullptr.
   * This function is **thread-safe**.
   *
   * @param nanoporeId a const reference to a std::string representing the id of the Vertex to be returned
   * @return A std::shared_ptr to the Vertex (nullptr if not found)
   */
  std::shared_ptr<Vertex> getVertexAsSharedPtr(std::string const &nanoporeId) const;

  /**
   * Returns a pointer pointing to the requested Vertex instance if found.
   * If no result is found, nullptr will be returned.
   * This function is **thread-safe**.
   *
   * @param nanoporeId a const reference to a std::string representing the id of the Vertex to be returned
   * @return A pointer pointing to the Vertex if found (nullptr if not found)
   */
  Vertex *getVertex(std::string const &nanoporeId) const;

  /**
   * Returns a std::vector containing pointers to all Vertex instances assigned to this Graph.
   *
   * @return A std::vector containing pointers to all Vertex instances
   */
  std::vector<Vertex *> getVertices() const;

  /**
   * Checks whether an Edge between Vertex instances exists or not.
   *
   * @param vertices a reference to a std::pair containing pointers to the Vertex instances to be checked for an Edge
   * @return A bool indicating whether an Edge exists or not
   */
  bool hasEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &vertices) const;

  /**
   * Checks whether an Edge between Vertex instances exists or not.
   *
   * @param vertices an rvalue reference to a std::pair containing pointers to the Vertex instances to be checked for
   *                 an Edge
   * @return A bool indicating whether an Edge exists or not
   */
  bool hasEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &&vertices) const;

  /**
   * Getter returning a specific Edge instance.
   * This function returns nullptr if the requested Edge is not present.
   *
   * @param vertices a reference to a std::pair containing pointers to the Vertex instances connected by the Edge
   * @return A pointer pointing to the requested Edge if found (nullptr if not found)
   */
  Edge *getEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &vertices) const;

  /**
   * Getter returning a specific Edge.
   * This function returns nullptr if the requested Edge wasn't found.
   *
   * @param vertices an rvalue reference to a std::pair containing pointers to the Vertex instances connected by the
   *                 Edge
   * @return A pointer pointing to the requested Edge if found (nullptr if not found)
   */
  Edge *getEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &&vertices) const;

  /**
   * Getter returning all Edge instances attached to this Graph.
   *
   * @return A std::vector containing pointers to all Edge instances attached to this Graph
   */
  std::vector<Edge *> getEdges() const;

  /**
   * Replaces the Edge store of the Graph with the supplied one.
   *
   * @param edges a std::unordered_map containing the shared_ptr instances pointing to the Edge instances mapped to
   *              their respective Edge identifiers
   */
  void replaceEdges(std::unordered_map<std::string, std::shared_ptr<Edge>> &&edges);

  /**
   * Getter returning the number of Vertex instances attached to the Graph.
   *
   * @return The number of Vertex instances attached to the Graph
   */
  std::size_t getOrder() const;

  /**
   * Getter returning the number of Edge instances attached to the Graph.
   *
   * @return The number of Edge instances attached to the Graph
   */
  std::size_t getSize() const;

protected:
  /**
   * Constructor creating a new instance.
   */
  GraphBase() = default;

  /**
   * Copy constructor.
   */
  GraphBase(GraphBase const &other);

  /**
   * Move constructor.
   */
  GraphBase(GraphBase &&other) noexcept;

  /**
   * Copy assignment operator.
   */
  GraphBase &operator=(GraphBase other);

  /**
   * Move assignment operator.
   */
  GraphBase &operator=(GraphBase &&other) noexcept;

  /**
   * Swaps two instances of GraphBase.
   *
   * @param lhs a reference to an instance of GraphBase
   * @param rhs a reference to an instance of GraphBase
   */
  friend void swap(GraphBase &lhs, GraphBase &rhs) noexcept;

  /**
   * Constructor creating a new instance.
   *
   * @param vertices an rvalue reference to a std::unordered_map mapping the std::shared_ptr instances pointing to
   *                 all Vertex instances to their ids
   * @param edges an rvalue reference to a std::vector containing the std::shared_ptr instances pointing to all
   *                 Edge instances
   * @param hasBidirectionalEdges a bool indicating whether Edge instances should be inserted bidirectional
   */
  GraphBase(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
            std::vector<std::shared_ptr<Edge>> &&edges, bool hasBidirectionalEdges);

  /**
   * Adds a std::shared_ptr to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr pointing to the Vertex to be added to the Graph
   */
  void _addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertex a pointer pointing to the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void _deleteVertex(gsl::not_null<Vertex const *> pVertex, bool hasBidirectionalEdges);

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   * This method is **thread-safe**.
   *
   * @param vertexIds a const reference to a std::pair containing pointers pointing to the Vertex instances to be
   *                  connected by the Edge
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void _addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &vertexIds,
                bool isBidirectional);

  /**
   * Deletes an Edge from the Graph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer pointing to the Edge to become deleted
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void _deleteEdge(gsl::not_null<Edge const *> pEdge, bool isBidirectional);

  /**
   * Getter returning the successors of a particular Vertex instance.
   *
   * @param vertexId a pointer pointing to the Vertex to return the successors for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> _getSuccessors(gsl::not_null<Vertex const *> pVertex) const;

  /**
   * Getter returning the predecessors of a particular Vertex instance.
   *
   * @param vertexId a pointer pointing to the Vertex to return the predecessors for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> _getPredecessors(gsl::not_null<Vertex const *> pVertex) const;

  /**
   * Hook which is getting called every time an Edge is about to get added.
   * This method can be overridden by deriving classes.
   *
   * @param pVertices a std::pair of pointers to the Vertex instances connected by the Edge which is
   *                  about to get added
   */
  virtual void
  _onEdgeAdded(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices);

  /**
   * Hook which is getting called every time an Edge is about to get deleted.
   * This method can be overridden by deriving classes.
   *
   * @param pVertices a std::pair of pointers to the Vertex instances connected by the Edge which is
   *                  about to get deleted
   */
  virtual void
  _onEdgeDeleted(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices);

private:
  template <class KEY, class VALUE> using um_t = std::unordered_map<KEY, VALUE>;
  um_t<std::string, std::shared_ptr<Vertex>> m_vertices; /*!< std::unordered_map containing all the Vertex instances */
  um_t<std::string, std::shared_ptr<Edge>>   m_edges;    /*!< std::unordered_map  containing all the Edge instances */
  um_t<std::string, um_t<std::string, Edge *const>>
                            m_adjacencyList; /*!< std::unordered_map  containing all the Edge instances */
  mutable std::shared_mutex m_mutexVertex;   /*!< std::shared_mutex for securing the parallel use of the Vertex map */
  mutable std::shared_mutex m_mutexEdge;     /*!< std::shared_mutex for securing the parallel use of the Edge map */

  /**
   * Adds an Edge to this Graph.
   * This is the class-internal version which is **not thread-safe**.
   *
   * @param spEdge an rvalue reference to the std::shared_ptr pointing to the Edge instance to be added to the Graph
   *               (by moving)
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void _addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional);
};

// -----------
// class Graph
// -----------

/**
 * Class representing an undirected Graph.
 *
 * In contrast to a directed Graph the insertion order of Edges is not preserved and Edge instances are inserted
 * bidirectional.
 */
class Graph final : public GraphBase {
public:
  /**
   * Constructor creating a new instance.
   */
  Graph() = default;

  /**
   * Constructor creating a new instance.
   *
   * @param vertices an rvalue reference to a std::unordered_map mapping the std::shared_ptr instances pointing to
   *                 all Vertex instances to their ids
   * @param edges an rvalue reference to a std::vector containing the std::shared_ptr instances pointing to all
   *                 Edge instances
   */
  Graph(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
        std::vector<std::shared_ptr<Edge>> &&                      edges);

  /**
   * Destructor.
   */
  ~Graph() override = default;

  /**
   * Copy constructor.
   */
  Graph(Graph const &other) = default;

  /**
   * Moving is disallowed.
   */
  Graph(Graph &&) = delete;

  /**
   * Copy assignment operator.
   */
  Graph &operator=(Graph other);

  /**
   * Move assignment is disallowed.
   */
  Graph &operator=(Graph &&) = delete;

  /**
   * Adds a std::shared_ptr pointing to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr pointing to the Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Adds a std::shared_ptr pointing to Vertex to this Graph if the Vertex is not already attached to the Graph.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr pointing to the Vertex to be added to the Graph
   */
  [[maybe_unused]] void addMissingVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertexId a pointer pointing to the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void deleteVertex(gsl::not_null<Vertex const *> pVertex);

  /**
   * Adds an Edge to this Graph. Already existing Edge instances are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIds a const reference to a std::pair containing pointers pointing to the Vertex instances to be
   *                  connected by the Edge
   */
  void addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &vertexIds);

  /**
   * Adds an Edge to this Graph. Already existing Edge instances are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIds an rvalue reference to a std::pair containing pointers pointing to the Vertex instances to be
   *                  connected by the Edge
   */
  void addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&vertexIds);

  /**
   * Deletes an Edge from the Graph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer pointing to the Edge to become deleted
   */
  void deleteEdge(gsl::not_null<Edge const *> pEdge);

  /**
   * Getter returning the neighbors of a particular Vertex instance.
   *
   * @param vertexId a pointer pointing to the Vertex to return the neighbors for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> getNeighbors(gsl::not_null<Vertex const *> pVertex) const;

  /**
   * Returns a Graph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a const reference to a std::vector containing pointers pointing to the Vertex instances inducing
   *                 the requested subgraph
   * @return A std::unique_ptr to the Graph representing the induced subgraph
   */
  std::unique_ptr<Graph> getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices);
};

// -------------
// class DiGraph
// -------------

/**
 * Class representing a directed Graph.
 *
 * In contrast to an undirected Graph the insertion order of Edges is preserved and Edges are not inserted
 * bidirectional.
 */
class DiGraph final : public GraphBase {
public:
  /**
   * Constructor creating a new instance.
   */
  DiGraph() = default;

  /**
   * Constructor creating a new instance.
   *
   * @param vertices an rvalue reference to a std::unordered_map mapping the std::shared_ptr instances pointing to
   *                 all Vertex instances to their ids
   * @param edges an rvalue reference to a std::vector containing the std::shared_ptr instances pointing to all
   *                 Edge instances
   */
  DiGraph(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
          std::vector<std::shared_ptr<Edge>> &&                      edges);

  /**
   * Destructor.
   */
  ~DiGraph() override = default;

  /**
   * Copy constructor.
   */
  DiGraph(DiGraph const &other);

  /**
   * Moving is disallowed.
   */
  DiGraph(DiGraph &&) = delete;

  /**
   * Copy assignment operator.
   */
  DiGraph &operator=(DiGraph other);

  /**
   * Move assignment is disallowed.
   */
  DiGraph &operator=(DiGraph &&) = delete;

  /**
   * Swaps two instances of DiGraph.
   *
   * @param lhs a reference to an instance of DiGraph
   * @param rhs a reference to an instance of DiGraph
   */
  friend void swap(DiGraph &lhs, DiGraph &rhs) noexcept;

  /**
   * Adds a std::shared_ptr to Vertex to this DiGraph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr pointing to the Vertex to be added to the DiGraph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Adds a std::shared_ptr to Vertex to this DiGraph if the Vertex is not already attached to the DiGraph.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr pointing to the Vertex to be added to the DiGraph
   */
  [[maybe_unused]] void addMissingVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Deletes a Vertex from the DiGraph.
   * If the Vertex is not assigned to another (Di-)Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertex a pointer pointing to the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void deleteVertex(gsl::not_null<Vertex const *> pVertex);

  /**
   * Adds an Edge to this DiGraph. Already existing Edge instances are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIds a reference to a std::pair containing pointers pointing to the Vertex instances to be connected by
   *                  the Edge
   */
  void addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &vertexIds);

  /**
   * Adds an Edge to this DiGraph. Already existing edges are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIds an rvalue reference to a std::pair containing pointers pointing to the Vertex instances to be
   *                  connected by the Edge
   */
  void addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&vertexIds);

  /**
   * Deletes an Edge from the DiGraph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer pointing to the Edge to become deleted
   */
  void deleteEdge(gsl::not_null<Edge const *> pEdge);

  /**
   * Getter returning the successors of a particular Vertex instance.
   *
   * @param vertexId a pointer pointing to the Vertex to return the successors for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> getSuccessors(gsl::not_null<Vertex const *> pVertex) const;

  /**
   * Getter returning the predecessors of a particular Vertex instance.
   *
   * @param vertexId a pointer pointing to the Vertex to return the predecessors for
   * @return A std::unordered_map mapping the ids of the connected Vertex instances to a pointer pointing to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> getPredecessors(gsl::not_null<Vertex const *> pVertex) const;

  /**
   * Returns a DiGraph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a const reference to a std::vector containing pointers pointing to the Vertex instances inducing
   *                 the requested subgraph
   * @return A std::unique_ptr to the DiGraph representing the induced subgraph
   */
  std::unique_ptr<DiGraph> getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices);

  /**
   * Returns the std::unordered_map containing the mapping of all Vertex instances to their in-degrees.
   *
   * @return A const reference to the std::unordered_map containing the mapping of all Vertex instances to their
   *         in-degrees
   */
  auto const &getInDegrees() const;

  /**
   * Returns the std::unordered_map containing the mapping of all Vertex instances to their out-degrees.
   *
   * @return A const reference to the std::unordered_map containing the mapping of all Vertex instances to their
   *         out-degrees
   */
  auto const &getOutDegrees() const;

private:
  std::unordered_map<Vertex const *, std::size_t> m_inDegrees;  /*!< std::unordered_map containing the in-degrees */
  std::unordered_map<Vertex const *, std::size_t> m_outDegrees; /*!< std::unordered_map containing the out-degrees */
  std::mutex m_mutexDegrees; /*!< std::mutex for securing the parallel use of the in-degree and out-degree maps */

  /**
   * Increases the in-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer pointing to the Vertex instance to update the in-degree for
   */
  void _increaseInDegree(Vertex const *pVertex);

  /**
   * Decreases the in-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer pointing to the Vertex instance to update the in-degree for
   */
  void _decreaseInDegree(Vertex const *pVertex);

  /**
   * Increases the out-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer pointing to the Vertex instance to update the out-degree for
   */
  void _increaseOutDegree(Vertex const *pVertex);

  /**
   * Decreases the out-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer pointing to the Vertex instance to update the out-degree for
   */
  void _decreaseOutDegree(Vertex const *pVertex);

  /**
   * Hook which is getting called every time an Edge is about to get added.
   * This method updates the degrees of the Vertex instances involved.
   *
   * @param pVertices a std::pair of pointers pointing to the Vertex instances connected by the Edge which is about to
   *                  get added
   */
  void _onEdgeAdded(
      std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) override;

  /**
   * Hook which is getting called every time an Edge is about to get deleted.
   * This method updates the degrees of the Vertex instances involved.
   *
   * @param pVertices a std::pair of pointers pointing to the Vertex instances connected by the Edge which is about to
   *                  get deleted
   */
  void _onEdgeDeleted(
      std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) override;

  /**
   * Recalculates the in-degrees and out-degrees based on the Edge instances attached to the DiGraph.
   */
  void _updateDegrees();
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ---------------
// class GraphBase
// ---------------

// PUBLIC CLASS METHODS

inline Edge *
GraphBase::getEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &&vertices) const {
  auto temp = std::move(vertices);
  return getEdge(temp);
}

inline bool
GraphBase::hasEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &&vertices) const {
  auto temp = std::move(vertices);
  return hasEdge(temp);
}

inline void GraphBase::replaceEdges(std::unordered_map<std::string, std::shared_ptr<Edge>> &&edges) {
  m_edges = std::move(edges);
}

inline std::size_t GraphBase::getOrder() const { return m_vertices.size(); }

inline std::size_t GraphBase::getSize() const { return m_edges.size(); }

// PROTECTED CLASS METHODS

inline GraphBase::GraphBase(GraphBase const &other)
    : m_vertices(other.m_vertices), m_edges(other.m_edges), m_adjacencyList(other.m_adjacencyList) {}

inline GraphBase::GraphBase(GraphBase &&other) noexcept : GraphBase() { swap(*this, other); }

inline GraphBase &GraphBase::operator=(GraphBase other) {
  swap(*this, other);

  return *this;
}

inline GraphBase &GraphBase::operator=(GraphBase &&other) noexcept {
  GraphBase tmp(std::move(other));

  swap(*this, tmp);
  return *this;
}

inline void swap(GraphBase &lhs, GraphBase &rhs) noexcept {
  using std::swap;

  swap(lhs.m_vertices, rhs.m_vertices);
  swap(lhs.m_edges, rhs.m_edges);
  swap(lhs.m_adjacencyList, rhs.m_adjacencyList);
}

inline void GraphBase::_onEdgeAdded(
    std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> && /*pVertices*/) {}

inline void GraphBase::_onEdgeDeleted(
    std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> && /*pVertices*/) {}

// -----------
// class Graph
// -----------

// PUBLIC CLASS METHODS

inline Graph::Graph(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
                    std::vector<std::shared_ptr<Edge>> &&                      edges)
    : GraphBase(std::move(vertices), std::move(edges), true) {}

inline Graph &Graph::operator=(Graph other) {
  if (&other != this) {
    swap(static_cast<GraphBase &>(*this), static_cast<GraphBase &>(other));
  }
  return *this;
}

inline void Graph::addVertex(std::shared_ptr<Vertex> &&spVertex) { _addVertex(std::move(spVertex)); }

[[maybe_unused]] inline void Graph::addMissingVertex(std::shared_ptr<Vertex> &&spVertex) {
  if (!hasVertex(spVertex->getId())) {
    addVertex(std::move(spVertex));
  }
}

inline void Graph::deleteVertex(gsl::not_null<Vertex const *> pVertex) { _deleteVertex(pVertex, true); }

inline void
Graph::addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &vertexIds) {
  _addEdge(vertexIds, true);
}

inline void
Graph::addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&vertexIds) {
  auto temp = std::move(vertexIds);
  addEdge(temp);
}

inline void Graph::deleteEdge(gsl::not_null<Edge const *> pEdge) { _deleteEdge(pEdge, true); }

inline std::unordered_map<std::string, Edge *const>
Graph::getNeighbors(gsl::not_null<Vertex const *> const pVertex) const {
  return _getSuccessors(pVertex);
}

// -------------
// class DiGraph
// -------------

// PUBLIC CLASS METHODS

inline DiGraph::DiGraph(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
                        std::vector<std::shared_ptr<Edge>> &&                      edges)
    : GraphBase(std::move(vertices), std::move(edges), false) {
  _updateDegrees();
}

inline DiGraph::DiGraph(DiGraph const &other)
    : GraphBase(other), m_inDegrees(other.m_inDegrees), m_outDegrees(other.m_outDegrees) {}

inline DiGraph &DiGraph::operator=(DiGraph other) {
  if (&other != this) {
    swap(*this, other);
  }
  return *this;
}

inline void swap(DiGraph &lhs, DiGraph &rhs) noexcept {
  using std::swap;

  swap(static_cast<GraphBase &>(lhs), static_cast<GraphBase &>(rhs));
  swap(lhs.m_inDegrees, rhs.m_inDegrees);
  swap(lhs.m_outDegrees, rhs.m_outDegrees);
}

[[maybe_unused]] inline void DiGraph::addMissingVertex(std::shared_ptr<Vertex> &&spVertex) {
  if (!hasVertex(spVertex->getId())) {
    addVertex(std::move(spVertex));
  }
}

inline void DiGraph::deleteVertex(gsl::not_null<Vertex const *> pVertex) {
  _deleteVertex(pVertex, false);

  m_inDegrees.erase(pVertex);
  m_outDegrees.erase(pVertex);
}

inline void
DiGraph::addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &vertexIds) {
  _addEdge(vertexIds, false);
}

inline void
DiGraph::addEdge(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&vertexIds) {
  auto temp = std::move(vertexIds);
  return addEdge(temp);
}

inline void DiGraph::deleteEdge(gsl::not_null<const Edge *> pEdge) { _deleteEdge(pEdge, false); }

inline std::unordered_map<std::string, Edge *const>
DiGraph::getSuccessors(gsl::not_null<Vertex const *> const pVertex) const {
  return _getSuccessors(pVertex);
}

inline std::unordered_map<std::string, Edge *const>
DiGraph::getPredecessors(gsl::not_null<Vertex const *> const pVertex) const {
  return _getPredecessors(pVertex);
}

inline auto const &DiGraph::getInDegrees() const { return m_inDegrees; }

inline auto const &DiGraph::getOutDegrees() const { return m_outDegrees; }

// PRIVATE CLASS METHODS

inline void
DiGraph::_onEdgeAdded(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) {
  _increaseOutDegree(pVertices.first);
  _increaseInDegree(pVertices.second);
}

inline void DiGraph::_onEdgeDeleted(
    std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) {
  _decreaseOutDegree(pVertices.first);
  _decreaseInDegree(pVertices.second);
}

} // namespace graph

// ---------------
// class GraphUtil
// ---------------

// PUBLIC CLASS METHODS

template <typename TYPE, typename std::enable_if_t<std::is_base_of<graph::GraphBase, TYPE>::value, bool>>
auto GraphUtil::getShortestPath(
    gsl::not_null<TYPE const *> const                                                                       pGraph,
    std::pair<gsl::not_null<graph::Vertex const *> const, gsl::not_null<graph::Vertex const *> const> const pVertices) {
  std::deque<graph::Vertex const *const>                           result;
  std::unordered_map<graph::Vertex const *, std::size_t>           dist;
  std::unordered_map<graph::Vertex const *, graph::Vertex const *> prev;

  std::list<graph::Vertex const *> vertices;

  for (auto const *const pVertex : pGraph->getVertices()) {
    vertices.push_back(pVertex);

    if (pVertex == pVertices.first) {
      dist.insert({pVertex, 0});
    } else {
      dist.insert({pVertex, std::numeric_limits<std::size_t>::max() - 1});
    }
  }

  while (!vertices.empty()) {
    auto const iterMinDistVertex =
        std::min_element(vertices.begin(), vertices.end(),
                         [&](auto const *const pV1, auto const *const pV2) { return dist.at(pV1) < dist.at(pV2); });
    auto pMinDistVertex = *iterMinDistVertex;

    vertices.erase(iterMinDistVertex);

    if (pMinDistVertex == pVertices.second) {
      if (prev.contains(pMinDistVertex) || pMinDistVertex == pVertices.first) {
        while (true) {
          result.push_front(pMinDistVertex);

          if (!prev.contains(pMinDistVertex)) {
            break;
          }

          pMinDistVertex = prev.at(pMinDistVertex);
        }
      }

      return result;
    }

    for (auto const &neighbor : _getReachableVertices(*pGraph, pMinDistVertex)) {
      auto const        alt       = dist.at(pMinDistVertex) + 1;
      auto const *const pNeighbor = pGraph->getVertex(neighbor.first);
      if (alt < dist.at(pNeighbor)) {
        dist.at(pNeighbor) = alt;
        prev.insert_or_assign(pNeighbor, pMinDistVertex);
      }
    }
  }

  return result;
}

// PRIVATE CLASS METHODS

[[maybe_unused]] inline auto
GraphUtil::_getReachableVertices(graph::Graph const &                                    graph,
                                 gsl::not_null<lazybastard::graph::Vertex const *> const pVertex) {
  return graph.getNeighbors(pVertex);
}

[[maybe_unused]] inline auto
GraphUtil::_getReachableVertices(graph::DiGraph const &                                  diGraph,
                                 gsl::not_null<lazybastard::graph::Vertex const *> const pVertex) {
  return diGraph.getSuccessors(pVertex);
}

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_GRAPH

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------