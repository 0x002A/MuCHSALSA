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
#include "Util.h"
#include "graph/Vertex.h"

namespace lazybastard {
namespace graph {

//// TYPES ////

/**
 * Class representing a Graph.
 *
 * A Graph holds a list of Vertex and Edge instances.
 * Instances of this class are designed to be **partially** thread-safe.
 */
class GraphBase {
public:
  /**
   * Constructor creating a new instance.
   */
  GraphBase() = default;

  /**
   * Constructor  creating a new instance.
   *
   * @param vertices an rvalue reference to a std::unordered_map containing the std::shared_ptr instances pointing to
   *                 all Vertex instances and their IDs
   * @param edges an rvalue reference to a std::vector containing the std::shared_ptr instances pointing to all
   *                 Edge instances
   * @param hasBidirectionalEdges a bool indicating whether Edge instances should be inserted bidirectional
   */
  GraphBase(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
            std::vector<std::shared_ptr<Edge>> &&edges, bool hasBidirectionalEdges);

  /**
   * Destructor.
   */
  ~GraphBase() = default;

  /**
   * Moving is disallowed.
   */
  GraphBase(GraphBase const &) = delete;

  /**
   * Copying is disallowed.
   */
  GraphBase(GraphBase &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  GraphBase &operator=(GraphBase &&) = delete;

  /**
   * Copy assignment is disallowed.
   */
  GraphBase &operator=(GraphBase const &) = delete;

  /**
   * Adds a std::shared_ptr to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr to Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Returns a bool indicating whether the supplied Vertex is present or not.
   *
   * @param nanoporeID a constant reference to a std::string representing the ID of the Vertex to be checked
   * @return A bool indicating whether the Vertex is present or not
   */
  bool hasVertex(std::string const &nanoporeID) const { return m_vertices.find(nanoporeID) != std::end(m_vertices); };

  /**
   * Returns a std::shared_ptr to the requested Vertex instance if found.
   * If no result is found, the std::shared_ptr will be initialized with nullptr.
   * This function is **thread-safe**.
   *
   * @param nanoporeID a constant reference to a std::string representing the ID of the Vertex to be returned
   * @return A std::shared_ptr to the Vertex if found
   */
  std::shared_ptr<Vertex> getVertex(std::string const &nanoporeID) const;

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
   * Getter returning all Edge instances attached to this Graph.
   *
   * @return A std::vector containing pointers to all Edge instances attached to this Graph
   */
  std::vector<Edge *> getEdges() const;

  /**
   * Fills the supplied std::vector with std::shared_ptr instances pointing to all Edge instances attached to this
   * Graph.
   *
   * @param result a reference to a std::vector which should receive the std::shared_ptr instances
   */
  void getEdges(std::vector<std::shared_ptr<Edge>> &result) const;

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

protected:
  /**
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertexID a pointer to a constant std::string representing the ID of the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void _deleteVertex(gsl::not_null<std::string const *> pVertexID, bool hasBidirectionalEdges); // NOLINT

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   *
   * @param vertexIDs a constant reference to a std::pair containing pointers to constant std::string instances
   *                  representing the IDs of the Vertex instances to be connected by the Edge
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void _addEdge(std::pair<std::string const, std::string const> const &vertexIDs, // NOLINT
                bool isBidirectional);

  /**
   * Deletes an Edge from the Graph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer to the Edge to become deleted
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void _deleteEdge(Edge const *pEdge, bool isBidirectional); // NOLINT

  /**
   * Getter returning the successors of a particular Vertex instance.
   *
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A pointer to the std::unordered_map containing the connected Vertex instances with their IDs and a constant
   *         pointer to the corresponding Edge instance, nullptr if the Vertex wasn't found
   */
  std::unordered_map<std::string, Edge *const> _getSuccessors(std::string const &vertexID) const; // NOLINT

  /**
   * Getter returning the predecessors of a particular Vertex instance.
   *
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A std::unordered_map containing the connected Vertex instances with their IDs and a constant pointer to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> _getPredecessors(std::string const &vertexID) const; // NOLINT

private:
  std::unordered_map<std::string, std::shared_ptr<Vertex>> m_vertices; /*!< Map containing all the Vertex instances */
  std::unordered_map<std::string, std::shared_ptr<Edge>> m_edges;      /*!< Map containing all the Edge instances */
  std::unordered_map<std::string, std::unordered_map<std::string, Edge *const>>
      m_adjacencyList;                     /*!< Map containing all the Edge instances */
  mutable std::shared_mutex m_mutexVertex; /*!< std::shared_mutex for securing the parallel use of the Vertex map */
  mutable std::shared_mutex m_mutexEdge;   /*!< std::shared_mutex for securing the parallel use of the Edge map */

  /**
   * Adds an Edge to this Graph.
   * This is the class-internal, **threadsafe** version.
   *
   * @param spEdge an rvalue reference to the std::shared_ptr to the Edge instance to be added to the Graph (by moving)
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional);
};

/**
 * Class representing an undirected Graph.
 *
 * In contrast to a directed Graph the insertion order of Edges is not preserved.
 */
class Graph final : public GraphBase {
public:
  /**
   * Constructor creating a new instance.
   */
  Graph() = default;

  /**
   * Constructor  creating a new instance.
   *
   * @param vertices an rvalue reference to a std::unordered_map containing the std::shared_ptr instances pointing to
   *                 all Vertex instances and their IDs
   * @param edges an rvalue reference to a std::vector containing the std::shared_ptr instances pointing to all
   *                 Edge instances
   */
  Graph(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices, std::vector<std::shared_ptr<Edge>> &&edges)
      : GraphBase(std::move(vertices), std::move(edges), true){};

  /**
   * Destructor.
   */
  ~Graph() = default;

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
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertexID a pointer to a constant std::string representing the ID of the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void deleteVertex(gsl::not_null<std::string const *> pVertexID) { _deleteVertex(pVertexID, true); };

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIDs a reference to a std::pair containing the IDs of the Vertex instances to be connected by the Edge
   */
  void addEdge(std::pair<std::string, std::string> &vertexIDs);

  /**
   * Adds an Edge to this DiGraph. Already existing edges are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIDs an rvalue reference to a std::pair containing the IDs of the Vertex instances to be connected by
   *                  the Edge
   */
  void addEdge(std::pair<std::string, std::string> &&vertexIDs) {
    auto temp = std::move(vertexIDs);
    addEdge(temp);
  }

  /**
   * Deletes an Edge from the Graph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer to the Edge to become deleted
   */
  void deleteEdge(Edge const *pEdge) { _deleteEdge(pEdge, true); };

  /**
   * Getter returning the neighbors of a particular Vertex instance.
   *
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A std::unordered_map containing the connected Vertex instances with their IDs and a constant pointer to the
   *         corresponding Edge instance, nullptr if the Vertex wasn't found
   */
  std::unordered_map<std::string, Edge *const> getNeighbors(std::string const &vertexID) const {
    return _getSuccessors(vertexID);
  };

  /**
   * Returns a Graph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a constant reference to a std::vector containing pointers to the Vertex instances inducing the
   *                 requested subgraph
   * @return A std::unique_ptr to the Graph representing the induced subgraph
   */
  std::unique_ptr<Graph> getSubgraph(std::vector<gsl::not_null<std::string const *>> const &vertices);
};

/**
 * Class representing a directed Graph.
 *
 * In contrast to an undirected Graph the insertion order of Edges is preserved.
 */
class DiGraph final : public GraphBase {
public:
  /**
   * Constructor creating a new instance.
   */
  DiGraph() = default;

  /**
   * Constructor  creating a new instance.
   *
   * @param vertices an rvalue reference to a std::unordered_map containing the std::shared_ptr instances pointing to
   *                 all Vertex instances and their IDs
   * @param edges an rvalue reference to a std::vector containing the std::shared_ptr instances pointing to all
   *                 Edge instances
   */
  DiGraph(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
          std::vector<std::shared_ptr<Edge>> &&edges)
      : GraphBase(std::move(vertices), std::move(edges), false){};

  /**
   * Destructor.
   */
  ~DiGraph() = default;

  /**
   * Moving is disallowed.
   */
  DiGraph(DiGraph const &) = delete;

  /**
   * Copying is disallowed.
   */
  DiGraph(DiGraph &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  DiGraph &operator=(DiGraph &&) = delete;

  /**
   * Copy assignment is disallowed.
   */
  DiGraph &operator=(DiGraph const &) = delete;

  /**
   * Deletes a Vertex from the DiGraph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertexID a pointer to a constant std::string representing the ID of the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void deleteVertex(gsl::not_null<std::string const *> pVertexID) { _deleteVertex(pVertexID, false); };

  /**
   * Adds an Edge to this DiGraph. Already existing edges are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIDs a reference to a std::pair containing the IDs of the Vertex instances to be connected by the Edge
   */
  void addEdge(std::pair<std::string, std::string> &vertexIDs);

  /**
   * Adds an Edge to this DiGraph. Already existing edges are omitted.
   * This function is **thread-safe**.
   *
   * @param vertexIDs an rvalue reference to a std::pair containing the IDs of the Vertex instances to be connected by
   *                  the Edge
   */
  void addEdge(std::pair<std::string, std::string> &&vertexIDs) {
    auto temp = std::move(vertexIDs);
    return addEdge(temp);
  }

  /**
   * Deletes an Edge from the Graph.
   * The memory will be cleaned up and all references will become invalid.
   *
   * @param pEdge a pointer to the Edge to become deleted
   */
  void deleteEdge(Edge const *pEdge) { _deleteEdge(pEdge, false); };

  /**
   * Getter returning the successors of a particular Vertex instance.
   *
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A std::unordered_map containing the connected Vertex instances with their IDs and a constant pointer to the
   *         corresponding Edge instance, nullptr if the Vertex wasn't found
   */
  std::unordered_map<std::string, Edge *const> getSuccessors(std::string const &vertexID) const {
    return _getSuccessors(vertexID);
  };

  /**
   * Getter returning the predecessors of a particular Vertex instance.
   *
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A std::unordered_map containing the connected Vertex instances with their IDs and a constant pointer to the
   *         corresponding Edge instance
   */
  std::unordered_map<std::string, Edge *const> getPredecessors(std::string const &vertexID) const {
    return _getPredecessors(vertexID);
  };

  /**
   * Returns a DiGraph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a constant reference to a std::vector containing pointers to the Vertex instances inducing the
   *                 requested subgraph
   * @return A std::unique_ptr to the DiGraph representing the induced subgraph
   */
  std::unique_ptr<DiGraph> getSubgraph(std::vector<gsl::not_null<std::string const *>> const &vertices);
};

//// INLINE DEFINITIONS ////

inline void Graph::addEdge(std::pair<std::string, std::string> &vertexIDs) {
  lazybastard::util::sortPair(vertexIDs);
  _addEdge(vertexIDs, true);
}

inline void DiGraph::addEdge(std::pair<std::string, std::string> &vertexIDs) { _addEdge(vertexIDs, false); }

inline std::size_t GraphBase::getOrder() const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  return m_vertices.size();
}

inline std::size_t GraphBase::getSize() const {
  std::shared_lock<std::shared_mutex> lck(m_mutexEdge);

  return m_edges.size();
}

} // namespace graph

//// UTILITY FUNCTIONS ////
struct GraphUtil {
  /**
   * Calculates the shortest path between two Vertex instances within a given Graph.
   *
   * @param pGraph a pointer to a Graph
   * @param vertexIDs a std::pair containing pointers to the start and end Vertex
   * @return A std::deque representing the shortest path between the supplied Vertex instances
   */
  template <typename T, typename std::enable_if_t<std::is_base_of<graph::GraphBase, T>::value, bool> = true>
  static std::deque<std::string const> getShortestPath(
      gsl::not_null<T const *> pGraph,
      std::pair<gsl::not_null<graph::Vertex const *> const, gsl::not_null<graph::Vertex const *> const> vertexIDs) {
    std::deque<std::string const> result;
    std::unordered_map<std::string, std::size_t> dist;
    std::unordered_map<std::string, std::string> prev;

    auto vertices = std::vector<std::string>(pGraph->getOrder());

    for (auto const *const pVertex : pGraph->getVertices()) {
      vertices.push_back(pVertex->getID());

      if (pVertex->getID() == vertexIDs.first->getID()) {
        dist[pVertex->getID()] = 0;
      } else {
        dist[pVertex->getID()] = std::numeric_limits<std::size_t>::max() - 1;
      }
    }

    while (vertices.size() > 0) {
      auto minDistVertex = *std::min_element(vertices.begin(), vertices.end(),
                                             [&](auto const &v1, auto const &v2) { return dist[v1] < dist[v2]; });

      vertices.erase(std::remove(vertices.begin(), vertices.end(), minDistVertex), vertices.end());

      if (minDistVertex == vertexIDs.second->getID()) {
        if (prev.contains(minDistVertex) || minDistVertex == vertexIDs.first->getID()) {
          while (true) {
            result.push_front(minDistVertex);

            if (!prev.contains(minDistVertex)) {
              break;
            }

            minDistVertex = prev[minDistVertex];
          }
        }

        return result;
      }

      for (auto const &pNeighbor : getReachableVertices(*pGraph, minDistVertex)) {
        auto const alt = dist[minDistVertex] + 1;
        if (alt < dist[pNeighbor.first]) {
          dist[pNeighbor.first] = alt;
          prev.insert_or_assign(pNeighbor.first, minDistVertex);
        }
      }
    }

    return result;
  };

private:
  static auto getReachableVertices(graph::Graph const &graph, std::string const &vertexID) {
    return graph.getNeighbors(vertexID);
  };
  static auto getReachableVertices(graph::DiGraph const &graph, std::string const &vertexID) {
    return graph.getSuccessors(vertexID);
  };
};
//// /UTILITY FUNCTIONS ////

} // namespace lazybastard
