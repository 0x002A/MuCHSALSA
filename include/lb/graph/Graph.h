#pragma once

#include <algorithm>
#include <cstddef>
#include <deque>
#include <gsl/pointers>
#include <iosfwd>
#include <iterator>
#include <memory>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Lb.fwd.h"

namespace lazybastard::graph {

//// FUNCTIONS ////

/**
 * Calculates the shortest path between two Vertex instances within a given Graph.
 *
 * @param pGraph a pointer to a Graph
 * @param vertexIDs a std::pair containing pointers to the start and end Vertex
 * @return A std::deque representing the shortest path between the supplied Vertex instances
 */
std::deque<gsl::not_null<std::string const *> const>
getShortestPath(gsl::not_null<Graph const *> pGraph,
                std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> vertexIDs);

//// TYPES ////

/**
 * Class representing a directed Graph.
 *
 * In contrast to an undirected Graph the insertion order of Edges is preserved.
 */
class DiGraph {};

/**
 * Class representing a Graph.
 *
 * A Graph holds a list of Vertex and Edge instances.
 * Instances of this class are designed to be **partially** thread-safe.
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
   * Adds a std::shared_ptr to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr to Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Returns a std::shared_ptr to the requested Vertex instance if found.
   * If no result is found, the std::shared_ptr will be initialized with nullptr.
   * This function is **thread-safe**.
   *
   * @param nanoporeID a constant reference to the ID of the Vertex to be returned
   * @return A std::shared_ptr to the Vertex if found
   */
  auto getVertex(std::string const &nanoporeID) const;

  /**
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertexID a pointer to a constant std::string representing the ID of the Vertex to be deleted
   */
  void deleteVertex(gsl::not_null<std::string const *> pVertexID);

  /**
   * Returns a std::vector containing pointers to all Vertex instances assigned to this Graph.
   *
   * @return A std::vector containing pointers to all Vertex instances
   */
  std::vector<Vertex const *> getVertices() const {
    std::vector<Vertex const *> vertices;

    std::transform(m_vertices.begin(), m_vertices.end(), std::back_inserter(vertices),
                   [](const std::unordered_map<std::string, std::shared_ptr<Vertex>>::value_type &pair) {
                     return pair.second.get();
                   });
    return vertices;
  }

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIDs a constant reference to a std::pair containing the IDs of the vertices to be connected by the Edge
   * @return The ID of the Edge
   */
  std::string addEdge(std::pair<std::string, std::string> const &vertexIDs);

  /**
   * Getter returning a specific Edge.
   * This function returns nullptr if the requested Edge wasn't found.
   *
   * @param vertexIDs a reference to a std::pair containing pointers to the IDs of the Vertex instances connected by the
   *                  Edge
   * @return A pointer to the requested Edge (constant) if found nullptr otherwise
   */
  Edge const *
  getEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const;

  /**
   * Getter returning a specific Edge.
   * This function returns nullptr if the requested Edge wasn't found.
   *
   * @param vertexIDs an rvalue reference to a std::pair containing pointers to the IDs of the Vertex instances
   *                  connected by the Edge
   * @return A pointer to the requested Edge (constant) if found nullptr otherwise
   */
  Edge const *
  getEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &&vertexIDs) const {
    auto temp = std::move(vertexIDs);
    return getEdge(temp);
  };

  /**
   * Deletes an Edge from the Graph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer to the Edge to become deleted
   */
  void deleteEdge(Edge const *pEdge);

  /**
   * Checks whether an Edge between Vertex instances exists or not.
   *
   * @param vertexIDs a reference to a std::pair containing pointers to the IDs of the Vertex instances to be
   *                  checked for an Edge
   * @return A bool indicating whether an Edge exists or not
   */
  bool hasEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const;

  /**
   * Checks whether an Edge between Vertex instances exists or not.
   *
   * @param vertexIDs an rvalue reference to a std::pair containing pointers to the IDs of the Vertex instances to be
   *                  checked for an Edge
   * @return A bool indicating whether an Edge exists or not
   */
  bool hasEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &&vertexIDs) const {
    auto temp = std::move(vertexIDs);
    return hasEdge(temp);
  };

  /**
   * Getter returning the Edge instances connected to a particular Vertex instance.
   *
   * @param vertexID a constant reference to the ID of the Vertex
   * @return A std::unordered_map containing the connected Vertex instances with a pointer to their IDs and a pointer to
   *         the corresponding Edge instance
   */
  std::unordered_map<std::string const *, Edge const *> getEdgesOfVertex(std::string const &vertexID) const;

  /**
   * Getter returning all Edge instances attached to this Graph.
   *
   * @return A std::vector containing pointers to all Edge instances attached to this Graph
   */
  std::vector<Edge *> getEdges() const;

  /**
   * Returns a Graph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a constant reference to a std::vector containing pointers to the Vertex instances inducing the
   *                 requested subgraph
   * @return A std::unique_ptr to the Graph representing the induced subgraph
   */
  std::unique_ptr<Graph> getSubgraph(std::vector<gsl::not_null<std::string const *>> const & /*vertices*/) {
    return std::make_unique<Graph>();
  }

  /**
   * Getter returning the number of Vertex instances attached to the Graph.
   * This function is **thread-safe**.
   *
   * @return The number of Vertex instances attached to the Graph
   */
  std::size_t getOrder() const;

  /**
   * Getter returning the number of Edge instances attached to the Graph.
   * This function is **thread-safe**.
   *
   * @return The number of Edge instances attached to the Graph
   */
  std::size_t getSize() const;

private:
  /**
   * Adds an Edge to this Graph.
   * This is the class-internal, **not-threadsafe** version.
   *
   * @param upEdge an rvalue reference to the unique pointer to the Edge instance to be added to the Graph (by moving)
   */
  void addEdgeInternal(std::unique_ptr<Edge> &&upEdge);

  std::unordered_map<std::string, std::shared_ptr<Vertex>> m_vertices; /*!< Map containing all the Vertex instances */
  std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<Edge>>>
      m_adjacencyList;                     /*!< Map containing all the Edge instances */
  std::size_t m_edgeCount{};               /*!< Number of Edge instances attached to the Graph */
  mutable std::shared_mutex m_mutexVertex; /*!< std::shared_mutex for securing the parallel use of the Vertex map */
  mutable std::shared_mutex m_mutexEdge;   /*!< std::shared_mutex for securing the parallel use of the Edge map */
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
