#include "graph/Graph.h"

#include <functional>
#include <set>
#include <stdexcept>
#include <string>

#include "Util.h"
#include "graph/Edge.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

//// HELPER ////
namespace {

void getVertexMap(std::unordered_map<std::string, std::shared_ptr<Vertex>> &result,
                  std::vector<lazybastard::graph::Vertex *> const &vertices) {
  result.clear();

  std::for_each(std::begin(vertices), std::end(vertices), [&](auto *pVertex) {
    if (pVertex != nullptr) {
      result.insert(std::make_pair(pVertex->getID(), pVertex->getSharedPtr()));
    }
  });
}

} // unnamed namespace
//// /HELPER ////

GraphBase::GraphBase(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
                     std::vector<std::shared_ptr<Edge>> &&edges, bool hasBidirectionalEdges)
    : m_vertices(std::move(vertices)) {
  for (auto &edge : edges) {
    auto pV1 = gsl::make_not_null(edge->getVertices().first);
    auto pV2 = gsl::make_not_null(edge->getVertices().second);
    if (hasVertex(pV1->getID()) && hasVertex(pV2->getID())) {
      _addEdgeInternal(std::move(edge), hasBidirectionalEdges);
    }
  }
}

void GraphBase::_addVertex(std::shared_ptr<Vertex> &&spVertex) {
  lazybastard::util::check_pointers(spVertex.get());
  std::unique_lock<std::shared_mutex> lck(m_mutexVertex);

  m_vertices.emplace(spVertex->getID(), std::move(spVertex));
}

std::shared_ptr<Vertex> GraphBase::getVertexAsSharedPtr(std::string const &nanoporeID) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeID);

  return iter != std::end(m_vertices) ? iter->second->getSharedPtr() : nullptr;
}

Vertex *GraphBase::getVertex(std::string const &nanoporeID) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeID);

  return iter != std::end(m_vertices) ? iter->second.get() : nullptr;
}

void GraphBase::_deleteVertex(gsl::not_null<Vertex const *> const pVertex, bool hasBidirectionalEdges) {
  auto const iterSuccessors = m_adjacencyList.find(pVertex->getID());

  if (iterSuccessors == std::end(m_adjacencyList)) {
    return;
  }

  using successor_t = decltype(iterSuccessors->second)::value_type;
  auto const eraseSuccessors = hasBidirectionalEdges //
                                   ? std::function{[&](successor_t const &s) {
                                       auto const iterPredecessor = m_adjacencyList.find(s.first);
                                       if (iterPredecessor != std::end(m_adjacencyList)) {
                                         iterPredecessor->second.erase(pVertex->getID());
                                       }

                                       auto const pEdge = gsl::make_not_null(s.second);
                                       _onEdgeDeleted(pEdge->getVertices());
                                       m_edges.erase(pEdge->getID());
                                     }}
                                   : std::function{[&](successor_t const &s) {
                                       auto const pEdge = gsl::make_not_null(s.second);
                                       _onEdgeDeleted(pEdge->getVertices());
                                       m_edges.erase(pEdge->getID());
                                     }};
  std::for_each(std::begin(iterSuccessors->second), std::end(iterSuccessors->second), eraseSuccessors);
  m_adjacencyList.erase(iterSuccessors);

  if (!hasBidirectionalEdges) {
    for (auto &[targetID, connectedVertices] : m_adjacencyList) {
      auto const iterPredecessor = connectedVertices.find(pVertex->getID());
      if (iterPredecessor != std::end(connectedVertices)) {
        auto const pEdge = gsl::make_not_null(iterPredecessor->second);
        _onEdgeDeleted(pEdge->getVertices());
        m_edges.erase(pEdge->getID());
        connectedVertices.erase(iterPredecessor);
      }
    }
  }

  m_vertices.erase(pVertex->getID());
}

void GraphBase::_addEdge(std::pair<std::string const, std::string const> const &vertexIDs, bool isBidirectional) {
  std::unique_lock<std::shared_mutex> lck(m_mutexEdge);

  auto pV1 = getVertexAsSharedPtr(vertexIDs.first);
  auto pV2 = getVertexAsSharedPtr(vertexIDs.second);

  if (!(pV1 != nullptr && pV2 != nullptr)) {
    // Edges between unknown vertices are omitted
    return;
  }

  auto const pVertices = std::make_pair(pV1.get(), pV2.get());

  auto spEdge = std::make_shared<Edge>(std::make_pair(std::move(pV1), std::move(pV2)));
  _addEdgeInternal(std::move(spEdge), isBidirectional);

  _onEdgeAdded(pVertices);
}

Edge *
GraphBase::getEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const {
  auto const outerIter = m_adjacencyList.find(*vertexIDs.first);
  if (outerIter != std::end(m_adjacencyList)) {
    auto const innerIter = outerIter->second.find(*vertexIDs.second);
    if (innerIter != std::end(outerIter->second)) {
      return innerIter->second;
    }
  }

  return nullptr;
}

void GraphBase::_deleteEdge(const Edge *const pEdge, bool isBirectional) {
  auto const findAndDeleteEdge = [&](std::string const &source, std::string const &target) {
    auto const outerIter = m_adjacencyList.find(source);

    if (outerIter != std::end(m_adjacencyList)) {
      auto const innerIter = outerIter->second.find(target);
      if (innerIter != std::end(outerIter->second)) {
        outerIter->second.erase(innerIter->first);
        auto const pEdge = gsl::make_not_null(innerIter->second);
        m_edges.erase(pEdge->getID());
      }
    }
  };
  auto const vertices = pEdge->getVertices();

  findAndDeleteEdge(vertices.first->getID(), vertices.second->getID());
  if (isBirectional) {
    findAndDeleteEdge(vertices.second->getID(), vertices.first->getID());
  }

  _onEdgeDeleted(vertices);
}

bool GraphBase::hasEdge(
    std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const {
  auto const iter = m_adjacencyList.find(*vertexIDs.first);

  if (iter != std::end(m_adjacencyList)) {
    return iter->second.find(*vertexIDs.second) != std::end(iter->second);
  }

  return false;
}

std::unordered_map<std::string, Edge *const> const *
GraphBase::_getSuccessors(gsl::not_null<Vertex const *> const pVertex) const {
  auto const outerIter = m_adjacencyList.find(pVertex->getID());
  if (outerIter != std::end(m_adjacencyList)) {
    return &outerIter->second;
  }

  return nullptr;
}

bool GraphBase::_getPredecessors(std::unordered_map<std::string, Edge *const> &result,
                                 gsl::not_null<Vertex const *> const pVertex) const {
  auto const *const pVertexID = &pVertex->getID();
  if (!hasVertex(*pVertexID)) {
    return false;
  }

  result.clear();
  for (auto const &[currentVertexID, connectedVertices] : m_adjacencyList) {
    if (*pVertexID == currentVertexID) {
      continue;
    }

    auto const iter = connectedVertices.find(*pVertexID);
    if (iter != std::end(connectedVertices)) {
      result.insert(std::make_pair(currentVertexID, iter->second));
    }
  }

  return true;
}

std::vector<Edge *> GraphBase::getEdges() const {
  std::vector<Edge *> edges;

  std::transform(std::begin(m_edges), std::end(m_edges), std::back_inserter(edges),
                 [](auto const &pair) { return pair.second.get(); });

  return edges;
}

void GraphBase::getEdges(std::vector<std::shared_ptr<Edge>> &result) const {
  result.clear();

  std::transform(std::begin(m_edges), std::end(m_edges), std::back_inserter(result),
                 [](const auto &pair) { return pair.second; });
}

void GraphBase::_addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional) {
  auto const emplaceEdge = [&](std::string const &source, std::string const &target, Edge *pEdge) {
    auto const iterAdjacencyList =
        m_adjacencyList.emplace(source, std::unordered_map<std::string, Edge *const>()).first;
    iterAdjacencyList->second.emplace(target, pEdge);
  };

  auto const pVertices = spEdge->getVertices();
  auto const iterEdges = m_edges.emplace(spEdge->getID(), std::move(spEdge)).first;
  auto *const pEdge = iterEdges != std::end(m_edges) ? iterEdges->second.get() : nullptr;

  emplaceEdge(pVertices.first->getID(), pVertices.second->getID(), pEdge);

  if (isBidirectional) {
    emplaceEdge(pVertices.second->getID(), pVertices.first->getID(), pEdge);
  }
}

std::unique_ptr<Graph> Graph::getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::unordered_map<std::string, std::shared_ptr<Vertex>> v(vertices.size());
  std::vector<std::shared_ptr<Edge>> e(getSize());

  getVertexMap(v, vertices);
  getEdges(e);

  return std::make_unique<Graph>(std::move(v), std::move(e));
}

std::unique_ptr<DiGraph> DiGraph::getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::unordered_map<std::string, std::shared_ptr<Vertex>> v(vertices.size());
  std::vector<std::shared_ptr<Edge>> e(getSize());

  getVertexMap(v, vertices);
  getEdges(e);

  return std::make_unique<DiGraph>(std::move(v), std::move(e));
}

void DiGraph::_updateDegrees() {
  m_inDegrees.clear();
  m_outDegrees.clear();

  auto const vertices = getVertices();
  std::for_each(std::begin(vertices), std::end(vertices), [this](auto const *const pVertex) {
    m_inDegrees.insert({pVertex, 0});
    m_outDegrees.insert({pVertex, 0});
  });

  auto const edges = getEdges();
  std::for_each(std::begin(edges), std::end(edges), [this](auto const *const pEdge) {
    auto const pVertices = pEdge->getVertices();
    _increaseOutDegree(pVertices.first);
    _increaseInDegree(pVertices.second);
  });
}

} // namespace lazybastard::graph
