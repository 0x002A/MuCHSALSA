#pragma once

#include <memory>
#include <mutex>
#include <unordered_map>
#include <string>

#include "graph/Edge.h"

namespace lazybastard {
namespace graph {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
   * Adds a shared pointer pointing to a Vertex to this Graph.
   *
   * @param spVertex the shared pointer to Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex>&& spVertex);

  /**
   * Returns a shared pointer to the requested Vertex instance if it could be found.
   * If no result could be found the shared pointer will be initialized with nullptr.
   *
   * @param nanoporeID the id of the Vertex to be returned
   */
  std::shared_ptr<Vertex> getVertex(const std::string& nanoporeID);

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   *
   * @param upEdge the unique pointer to Edge to be added to the Graph
   */
  void addEdge(const std::pair<std::string, std::string>& vertexIDs);

  /**
   * Returns the number of vertices attached to the graph.
   */
  const std::size_t getOrder();

  /**
   * Returns the number of edge attached to the graph.
   */
  const std::size_t getSize();
private:
  /**
   * Adds an Edge to this Graph.
   * This is the class-internal, unsynchronised version.
   *
   * @param upEdge the unique pointer to Edge to be added to the Graph
   */
  void addEdgeInternal(std::unique_ptr<Edge>&& upEdge);

  std::unordered_map<std::string, std::shared_ptr<Vertex>> m_vertices; /*!< Map containing all the vertices */
  std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<Edge>>> m_adjacencyList; /*!< Map containing all the edges */
  std::size_t m_edgeCount; /*!< Number of edges attached to the graph */
  std::mutex m_mutexVertex; /*!< Mutex for securing the parallel use of the vertex map */
  std::mutex m_mutexEdge; /*!< Mutex for securing the parallel use of the edge map */
};

// Inline definitions
inline const std::size_t
Graph::getOrder()
{
  std::scoped_lock<std::mutex> lck(m_mutexVertex);

  return m_vertices.size();
}

inline const std::size_t
Graph::getSize()
{
  std::scoped_lock<std::mutex> lck(m_mutexEdge);

  return m_edgeCount;
}

}
}
