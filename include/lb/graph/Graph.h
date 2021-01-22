#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <shared_mutex>
#include <unordered_map>

namespace lazybastard::graph {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Edge;
class Vertex;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * Class representing a graph.
 *
 * A Graph holds a list of vertices and edges.
 * Instances of this class are designed to be thread-safe.
 */
class Graph {
public:
  /**
   * Constructor.
   */
  Graph();

  /**
   * Destructor.
   */
  ~Graph();

  /**
   * Moving is disallowed.
   */
  Graph(Graph const &) = delete;

  /**
   * Copying is disallowed.
   */
  Graph(Graph &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Graph &operator=(Graph &&) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Graph &operator=(Graph const &) = delete;

  /**
   * Adds a shared pointer pointing to a Vertex to this Graph.
   *
   * @param spVertex the shared pointer to Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Returns a shared pointer to the requested Vertex instance if it could be found.
   * If no result could be found the shared pointer will be initialized with nullptr.
   *
   * @param nanoporeID the id of the Vertex to be returned
   * @return A shared pointer to the Vertex if found
   */
  auto getVertex(std::string const &nanoporeID) const;

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   *
   * @param vertexIDs the IDs of the vertices to be connected by an Edge
   */
  std::string addEdge(std::pair<std::string, std::string> const &vertexIDs);

  /**
   * Getter for the adjacency list.
   *
   * @return The adjacency list of the graph
   */
  auto const &getAdjacencyList() const { return m_adjacencyList; };

  /**
   * Getter for the number of Vertex instances attached to the Graph.
   *
   * @return The number of Vertex instances attached to the Graph
   */
  std::size_t getOrder() const;

  /**
   * Getter for the number of Edge instances attached to the Graph.
   *
   * @return The number of Edge instances attached to the Graph
   */
  std::size_t getSize() const;

private:
  /**
   * Adds an Edge to this Graph.
   * This is the class-internal, unsynchronised version.
   *
   * @param upEdge the unique pointer to Edge to be added to the Graph
   */
  void addEdgeInternal(std::unique_ptr<Edge> &&upEdge);

  std::unordered_map<std::string, std::shared_ptr<Vertex>> m_vertices; /*!< Map containing all the Vertex instances */
  std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<Edge>>>
      m_adjacencyList;                     /*!< Map containing all the Edge instances */
  std::size_t m_edgeCount{};               /*!< Number of Edge instances attached to the Graph */
  mutable std::shared_mutex m_mutexVertex; /*!< Mutex for securing the parallel use of the Vertex map */
  mutable std::shared_mutex m_mutexEdge;   /*!< Mutex for securing the parallel use of the Edge map */
};

// Inline definitions
inline std::size_t Graph::getOrder() const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  return m_vertices.size();
}

inline std::size_t Graph::getSize() const {
  std::shared_lock<std::shared_mutex> lck(m_mutexEdge);

  return m_edgeCount;
}

} // namespace lazybastard::graph
