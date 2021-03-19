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
  virtual ~GraphBase() = default;

  /**
   * Copy constructor.
   */
  GraphBase(GraphBase const &other)
      : m_vertices(other.m_vertices), m_edges(other.m_edges), m_adjacencyList(other.m_adjacencyList){};

  /**
   * Move constructor.
   */
  GraphBase(GraphBase &&other) noexcept : GraphBase() { swap(*this, other); };

  /**
   * Copy assignment operator.
   */
  GraphBase &operator=(GraphBase other) {
    swap(*this, other);

    return *this;
  };

  /**
   * Move assignment operator.
   */
  GraphBase &operator=(GraphBase &&other) noexcept {
    GraphBase tmp(std::move(other));

    swap(*this, tmp);
    return *this;
  };

  /**
   * Swaps two instances of GraphBase.
   *
   * @param first a reference to an instance of GraphBase
   * @param second a reference to an instance of GraphBase
   */
  friend void swap(GraphBase &first, GraphBase &second) noexcept {
    using std::swap;

    swap(first.m_vertices, second.m_vertices);
    swap(first.m_edges, second.m_edges);
    swap(first.m_adjacencyList, second.m_adjacencyList);
  }

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
  std::shared_ptr<Vertex> getVertexAsSharedPtr(std::string const &nanoporeID) const;

  /**
   * Returns a pointer to the requested Vertex instance if found.
   * If no result is found, nullptr will be returned.
   * This function is **thread-safe**.
   *
   * @param nanoporeID a constant reference to a std::string representing the ID of the Vertex to be returned
   * @return A pointer to the Vertex if found, nullptr otherwise
   */
  Vertex *getVertex(std::string const &nanoporeID) const;

  /**
   * Returns a std::vector containing pointers to all Vertex instances assigned to this Graph.
   *
   * @return A std::vector containing pointers to all Vertex instances
   */
  std::vector<Vertex *> getVertices() const {
    std::vector<Vertex *> vertices;

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
   * @return A pointer to the requested Edge if found nullptr otherwise
   */
  Edge *getEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const;

  /**
   * Getter returning a specific Edge.
   * This function returns nullptr if the requested Edge wasn't found.
   *
   * @param vertexIDs an rvalue reference to a std::pair containing pointers to the IDs of the Vertex instances
   *                  connected by the Edge
   * @return A pointer to the requested Edge if found nullptr otherwise
   */
  Edge *getEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &&vertexIDs) const {
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
   * Replaces the Edge store of the Graph with the supplied one.
   *
   * @param edges a std::unordered_map containing the shared_ptr instances pointing to the Edge instances mapped to
   *              their respective Edge identifiers
   */
  void replaceEdges(std::unordered_map<std::string, std::shared_ptr<Edge>> &&edges) { m_edges = std::move(edges); }

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
   * Adds a std::shared_ptr to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr to Vertex to be added to the Graph
   */
  void _addVertex(std::shared_ptr<Vertex> &&spVertex);

  /**
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertex a pointer to the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void _deleteVertex(gsl::not_null<Vertex const *> pVertex, bool hasBidirectionalEdges); // NOLINT

  /**
   * Adds an Edge to this Graph. Already existing edges are omitted.
   * This method is **thread-safe**.
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
   * @return A pointer to a constant std::unordered_map containing the connected Vertex instances with their IDs and a
   *         constant pointer to the corresponding Edge instance, nullptr if the Vertex wasn't found or hasn't got any
   *         successive Vertex instances
   */
  std::unordered_map<std::string, Edge *const> const *_getSuccessors(std::string const &vertexID) const; // NOLINT

  /**
   * Getter returning the predecessors of a particular Vertex instance.
   *
   * @param result a reference to a std::unordered_map receiving the connected Vertex instances with their IDs and a
   * constant pointer to the corresponding Edge instance
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A bool indicating whether the Vertex has been found
   */
  bool _getPredecessors(std::unordered_map<std::string, Edge *const> &result,
                        std::string const &vertexID) const; // NOLINT

  /**
   * Hook which is getting called every time an Edge is about to get added.
   * This method can be overridden by child classes.
   *
   * @param pVertices a std::pair of const pointers to the Vertex instances connected by the Edge which is about to get
   *                  added
   */
  virtual void
  _onEdgeAdded(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) {
    LB_UNUSED(pVertices);
  };

  /**
   * Hook which is getting called every time an Edge is about to get deleted.
   * This method can be overridden by child classes.
   *
   * @param pVertices a std::pair of const pointers to the Vertex instances connected by the Edge which is about to get
   *                  deleted
   */
  virtual void
  _onEdgeDeleted(std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) {
    LB_UNUSED(pVertices);
  };

private:
  std::unordered_map<std::string, std::shared_ptr<Vertex>> m_vertices; /*!< Map containing all the Vertex instances */
  std::unordered_map<std::string, std::shared_ptr<Edge>> m_edges;      /*!< Map containing all the Edge instances */
  std::unordered_map<std::string, std::unordered_map<std::string, Edge *const>>
      m_adjacencyList;                     /*!< Map containing all the Edge instances */
  mutable std::shared_mutex m_mutexVertex; /*!< std::shared_mutex for securing the parallel use of the Vertex map */
  mutable std::shared_mutex m_mutexEdge;   /*!< std::shared_mutex for securing the parallel use of the Edge map */

  /**
   * Adds an Edge to this Graph.
   * This is the class-internal, unsynchronized version.
   *
   * @param spEdge an rvalue reference to the std::shared_ptr to the Edge instance to be added to the Graph (by moving)
   * @param isBidirectional a bool indicating whether the Edge is bidirectional or not
   */
  void _addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional);
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
  Graph &operator=(Graph other) {
    if (&other != this) {
      swap(static_cast<GraphBase &>(*this), static_cast<GraphBase &>(other));
    }
    return *this;
  };

  /**
   * Move assignment is disallowed.
   */
  Graph &operator=(Graph &&) = delete;

  /**
   * Adds a std::shared_ptr to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr to Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex) { _addVertex(std::move(spVertex)); };

  /**
   * Deletes a Vertex from the Graph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertexID a pointer to the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void deleteVertex(gsl::not_null<Vertex const *> pVertex) { _deleteVertex(pVertex, true); };

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
   * @return A pointer to a constant std::unordered_map containing the connected Vertex instances with their IDs and a
   *         constant pointer to the corresponding Edge instance, nullptr if the Vertex wasn't found or hasn't got any
   *         adjacent Vertex instances
   */
  std::unordered_map<std::string, Edge *const> const *getNeighbors(std::string const &vertexID) const {
    return _getSuccessors(vertexID);
  };

  /**
   * Returns a Graph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a constant reference to a std::vector containing pointers to the Vertex instances inducing the
   *                 requested subgraph
   * @return A std::unique_ptr to the Graph representing the induced subgraph
   */
  std::unique_ptr<Graph> getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices);
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
      : GraphBase(std::move(vertices), std::move(edges), false) {
    _updateDegrees();
  };

  /**
   * Destructor.
   */
  ~DiGraph() override = default;

  /**
   * Copy constructor.
   */
  DiGraph(DiGraph const &other) : GraphBase(other){};

  /**
   * Moving is disallowed.
   */
  DiGraph(DiGraph &&) = delete;

  /**
   * Copy assignment operator.
   */
  DiGraph &operator=(DiGraph other) {
    if (&other != this) {
      swap(*this, other);
    }
    return *this;
  };

  /**
   * Move assignment is disallowed.
   */
  DiGraph &operator=(DiGraph &&) = delete;

  /**
   * Swaps two instances of GraphBase.
   *
   * @param first a reference to an instance of DiGraph
   * @param second a reference to an instance of DiGraph
   */
  friend void swap(DiGraph &first, DiGraph &second) noexcept {
    using std::swap;

    swap(static_cast<GraphBase &>(first), static_cast<GraphBase &>(second));
    swap(first.m_inDegrees, second.m_inDegrees);
    swap(first.m_outDegrees, second.m_outDegrees);
  }

  /**
   * Adds a std::shared_ptr to Vertex to this Graph.
   * This function is **thread-safe**.
   *
   * @param spVertex an rvalue reference to the std::shared_ptr to Vertex to be added to the Graph
   */
  void addVertex(std::shared_ptr<Vertex> &&spVertex) {
    auto const *const pVertex = spVertex.get();

    _addVertex(std::move(spVertex));

    {
      std::unique_lock<std::mutex> lck(m_mutexDegrees);
      m_inDegrees.insert({pVertex, 0});
      m_outDegrees.insert({pVertex, 0});
    }
  };

  /**
   * Deletes a Vertex from the DiGraph.
   * If the Vertex is not assigned to another Graph the memory will be cleaned up and all references
   * will become invalid.
   *
   * @param pVertex a pointer to the Vertex to be deleted
   * @param hasBidirectionalEdges a bool indicating whether the Vertex is connected by bidirectional Edge instances
   */
  void deleteVertex(gsl::not_null<Vertex const *> pVertex) {
    _deleteVertex(pVertex, false);

    m_inDegrees.erase(pVertex);
    m_outDegrees.erase(pVertex);
  };

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
   * @return A pointer to a constant std::unordered_map containing the connected Vertex instances with their IDs and a
   *         constant pointer to the corresponding Edge instance, nullptr if the Vertex wasn't found
   */
  std::unordered_map<std::string, Edge *const> const *getSuccessors(std::string const &vertexID) const {
    return _getSuccessors(vertexID);
  };

  /**
   * Getter returning the predecessors of a particular Vertex instance.
   *
   * @param result a reference to a std::unordered_map receiving the connected Vertex instances with their IDs and a
   * constant pointer to the corresponding Edge instance
   * @param vertexID a constant reference to a std::string representing the ID of the Vertex
   * @return A bool indicating whether the Vertex has been found
   */
  bool getPredecessors(std::unordered_map<std::string, Edge *const> &result, std::string const &vertexID) const {
    return _getPredecessors(result, vertexID);
  };

  /**
   * Returns a DiGraph representing the subgraph induced by the supplied Vertex instances.
   *
   * @param vertices a constant reference to a std::vector containing pointers to the Vertex instances inducing the
   *                 requested subgraph
   * @return A std::unique_ptr to the DiGraph representing the induced subgraph
   */
  std::unique_ptr<DiGraph> getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices);

  /**
   * Returns the std::unordered_map containing the mapping of all Vertex instances to their in-degrees via constant
   * reference.
   *
   * @return The std::unordered_map containing the mapping of all Vertex instances to their in-degrees via constant
   * reference
   */
  auto const &getInDegrees() const { return m_inDegrees; };

  /**
   * Returns the std::unordered_map containing the mapping of all Vertex instances to their out-degrees via constant
   * reference.
   *
   * @return The std::unordered_map containing the mapping of all Vertex instances to their out-degrees via constant
   * reference
   */
  auto const &getOutDegrees() const { return m_outDegrees; };

private:
  std::unordered_map<Vertex const *, std::size_t> m_inDegrees;  /*!< Map containing the in-degrees */
  std::unordered_map<Vertex const *, std::size_t> m_outDegrees; /*!< Map containing the out-degrees */
  std::mutex m_mutexDegrees; /*!< std::mutex for securing the parallel use of the in-degree and out-degree maps */

  /**
   * Increases the in-degree value for the given Verrtex instance.
   *
   * @param pVertex a pointer to the Vertex instance the in-degree should be updated for
   */
  void _increaseInDegree(Vertex const *const pVertex) {
    auto const iter = m_inDegrees.find(pVertex);
    if (iter == std::end(m_inDegrees)) {
      return;
    }

    ++iter->second;
  }

  /**
   * Decreases the in-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer to the Vertex instance the in-degree should be updated for
   */
  void _decreaseInDegree(Vertex const *const pVertex) {
    auto const iter = m_inDegrees.find(pVertex);
    if (iter == std::end(m_inDegrees)) {
      return;
    }

    --iter->second;
  }

  /**
   * Increases the out-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer to the Vertex instance the out-degree should be updated for
   */
  void _increaseOutDegree(Vertex const *const pVertex) {
    auto const iter = m_outDegrees.find(pVertex);
    if (iter == std::end(m_outDegrees)) {
      return;
    }

    ++iter->second;
  }

  /**
   * Decreases the out-degree value for the given Vertex instance.
   *
   * @param pVertex a pointer to the Vertex instance the out-degree should be updated for
   */
  void _decreaseOutDegree(Vertex const *const pVertex) {
    auto const iter = m_outDegrees.find(pVertex);
    if (iter == std::end(m_outDegrees)) {
      return;
    }

    --iter->second;
  }

  /**
   * Hook which is getting called every time an Edge is about to get added.
   * This method updates the degrees of the Vertex instances involved.
   *
   * @param pVertices a std::pair of const pointers to the Vertex instances connected by the Edge which is about to get
   *                  added
   */
  void _onEdgeAdded(
      std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) override {
    _increaseOutDegree(pVertices.first);
    _increaseInDegree(pVertices.second);
  };

  /**
   * Hook which is getting called every time an Edge is about to get deleted.
   *
   * @param pVertices a std::pair of const pointers to the Vertex instances connected by the Edge which is about to get
   *                  deleted
   */
  void _onEdgeDeleted(
      std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> &&pVertices) override {
    _decreaseOutDegree(pVertices.first);
    _decreaseInDegree(pVertices.second);
  };

  void _updateDegrees();
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
   * @return A std::deque containing pointers to the Vertex instances representing the shortest path between the
   *         supplied Vertex instances
   */
  template <typename T, typename std::enable_if_t<std::is_base_of<graph::GraphBase, T>::value, bool> = true>
  static std::deque<gsl::not_null<graph::Vertex const *> const>
  getShortestPath(gsl::not_null<T const *> const pGraph,
                  std::pair<gsl::not_null<graph::Vertex const *> const,
                            gsl::not_null<graph::Vertex const *> const> const vertexIDs) {
    std::deque<gsl::not_null<graph::Vertex const *> const> result;
    std::unordered_map<std::string const *, std::size_t> dist;
    std::unordered_map<graph::Vertex const *, graph::Vertex const *> prev;

    std::vector<graph::Vertex const *> vertices;

    for (auto const *const pVertex : pGraph->getVertices()) {
      vertices.push_back(pVertex);

      if (pVertex == vertexIDs.first) {
        dist[&pVertex->getID()] = 0;
      } else {
        dist[&pVertex->getID()] = std::numeric_limits<std::size_t>::max() - 1;
      }
    }

    while (!vertices.empty()) {
      auto pMinDistVertex =
          *std::min_element(vertices.begin(), vertices.end(), [&](auto const *const pV1, auto const *const pV2) {
            return dist[&pV1->getID()] < dist[&pV2->getID()];
          });

      vertices.erase(std::remove(vertices.begin(), vertices.end(), pMinDistVertex), vertices.end());

      if (pMinDistVertex == vertexIDs.second) {
        if (prev.contains(pMinDistVertex) || pMinDistVertex == vertexIDs.first) {
          while (true) {
            result.push_front(pMinDistVertex);

            if (!prev.contains(pMinDistVertex)) {
              break;
            }

            pMinDistVertex = prev[pMinDistVertex];
          }
        }

        return result;
      }

      for (auto const &neighbor : *getReachableVertices(*pGraph, pMinDistVertex->getID())) {
        auto const alt = dist[&pMinDistVertex->getID()] + 1;
        auto const *const pNeighbor = pGraph->getVertex(neighbor.first);
        if (alt < dist[&pNeighbor->getID()]) {
          dist[&pNeighbor->getID()] = alt;
          prev.insert_or_assign(pNeighbor, pMinDistVertex);
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
