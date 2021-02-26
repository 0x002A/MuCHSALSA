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

void getVertexMap(std::unordered_map<std::string, std::shared_ptr<Vertex>> &result, GraphBase const &graph,
                  std::vector<gsl::not_null<std::string const *>> const &vertices) {
  result.clear();

  std::for_each(std::begin(vertices), std::end(vertices), [&](auto const pVertexID) {
    auto const spVertex = graph.getVertex(*pVertexID);

    if (spVertex != nullptr) {
      result.insert(std::make_pair(spVertex->getID(), std::move(spVertex)));
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
      addEdgeInternal(std::move(edge), hasBidirectionalEdges);
    }
  }
}

void GraphBase::addVertex(std::shared_ptr<Vertex> &&spVertex) {
  lazybastard::util::check_pointers(spVertex.get());
  std::unique_lock<std::shared_mutex> lk(m_mutexVertex);

  m_vertices.emplace(spVertex->getID(), std::move(spVertex));
}

std::shared_ptr<Vertex> GraphBase::getVertex(std::string const &nanoporeID) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeID);

  return iter != std::end(m_vertices) ? iter->second->getSharedPtr() : nullptr;
}

void GraphBase::_deleteVertex(gsl::not_null<std::string const *> const pVertexID, bool hasBidirectionalEdges) {
  auto const iterSuccessors = m_adjacencyList.find(*pVertexID);

  if (iterSuccessors == std::end(m_adjacencyList)) {
    return;
  }

  using successor_t = decltype(iterSuccessors->second)::value_type;
  auto const eraseSuccessors = hasBidirectionalEdges //
                                   ? std::function{[&](successor_t const &s) {
                                       auto const iterPredecessor = m_adjacencyList.find(s.first);
                                       if (iterPredecessor != std::end(m_adjacencyList)) {
                                         iterPredecessor->second.erase(*pVertexID);
                                       }

                                       auto const pEdge = gsl::make_not_null(s.second);
                                       m_edges.erase(pEdge->getID());
                                     }}
                                   : std::function{[&](successor_t const &s) {
                                       auto const pEdge = gsl::make_not_null(s.second);
                                       m_edges.erase(pEdge->getID());
                                     }};
  std::for_each(std::begin(iterSuccessors->second), std::end(iterSuccessors->second), eraseSuccessors);
  m_adjacencyList.erase(iterSuccessors);

  if (!hasBidirectionalEdges) {
    for (auto &[targetID, connectedVertices] : m_adjacencyList) {
      auto const iterPredecessor = connectedVertices.find(*pVertexID);
      if (iterPredecessor != std::end(connectedVertices)) {
        auto const pEdge = gsl::make_not_null(iterPredecessor->second);
        m_edges.erase(pEdge->getID());
        connectedVertices.erase(iterPredecessor);
      }
    }
  }

  m_vertices.erase(*pVertexID);
}

void GraphBase::_addEdge(std::pair<std::string const, std::string const> const &vertexIDs, bool isBidirectional) {
  auto pV1 = getVertex(vertexIDs.first);
  auto pV2 = getVertex(vertexIDs.second);

  if (!(pV1 != nullptr && pV2 != nullptr)) {
    // Edges between unknown vertices are omitted
    return;
  }

  auto vertexPair = std::make_pair(std::move(pV1), std::move(pV2));

  auto spEdge = std::make_shared<Edge>(std::move(vertexPair));
  addEdgeInternal(std::move(spEdge), isBidirectional);
}

Edge const *
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
}

bool GraphBase::hasEdge(
    std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const {
  auto const iter = m_adjacencyList.find(*vertexIDs.first);

  if (iter != std::end(m_adjacencyList)) {
    return iter->second.find(*vertexIDs.second) != std::end(iter->second);
  }

  return false;
}

std::unordered_map<std::string, Edge *const> GraphBase::_getSuccessors(std::string const &vertexID) const {
  auto const outerIter = m_adjacencyList.find(vertexID);
  if (outerIter != std::end(m_adjacencyList)) {
    return outerIter->second;
  }

  return decltype(outerIter->second)();
}

std::unordered_map<std::string, Edge *const> GraphBase::_getPredecessors(const std::string &vertexID) const {
  std::unordered_map<std::string, Edge *const> predecessors;

  for (auto const &[currentVertexID, connectedVertices] : m_adjacencyList) {
    if (vertexID == currentVertexID) {
      continue;
    }

    auto const iter = connectedVertices.find(vertexID);
    if (iter != std::end(connectedVertices)) {
      predecessors.insert(std::make_pair(currentVertexID, iter->second));
    }
  }

  return predecessors;
}

std::vector<Edge *> GraphBase::getEdges() const {
  std::vector<Edge *> edges;

  std::transform(std::begin(m_edges), std::end(m_edges), std::back_inserter(edges),
                 [](const auto &pair) { return pair.second.get(); });

  return edges;
}

void GraphBase::getEdges(std::vector<std::shared_ptr<Edge>> &result) const {
  result.clear();

  std::transform(std::begin(m_edges), std::end(m_edges), std::back_inserter(result),
                 [](const auto &pair) { return pair.second; });
}

void GraphBase::addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional) {
  std::unique_lock<std::shared_mutex> lck(m_mutexEdge);

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

std::unique_ptr<Graph> Graph::getSubgraph(std::vector<gsl::not_null<std::string const *>> const &vertices) {
  std::unordered_map<std::string, std::shared_ptr<Vertex>> v(vertices.size());
  std::vector<std::shared_ptr<Edge>> e(getSize());

  getVertexMap(v, *this, vertices);
  getEdges(e);

  return std::make_unique<Graph>(std::move(v), std::move(e));
}

std::unique_ptr<DiGraph> DiGraph::getSubgraph(std::vector<gsl::not_null<std::string const *>> const &vertices) {
  std::unordered_map<std::string, std::shared_ptr<Vertex>> v(vertices.size());
  std::vector<std::shared_ptr<Edge>> e(getSize());

  getVertexMap(v, *this, vertices);
  getEdges(e);

  return std::make_unique<DiGraph>(std::move(v), std::move(e));
}

} // namespace lazybastard::graph
