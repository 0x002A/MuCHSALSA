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

#include "graph/Graph.h"

#include <stdexcept>

#include "Util.h"
#include "graph/Edge.h"

namespace lazybastard::graph {

// =====================================================================================================================
//                                                        HELPER
// =====================================================================================================================

namespace {

auto getVertexMap(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::unordered_map<std::string, std::shared_ptr<Vertex>> result;

  std::for_each(std::begin(vertices), std::end(vertices), [&](auto *pVertex) {
    if (pVertex) {
      result.insert(std::make_pair(pVertex->getId(), pVertex->getSharedPtr()));
    }
  });

  return result;
}

} // unnamed namespace

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ---------------
// class GraphBase
// ---------------

// PUBLIC CLASS METHODS

bool GraphBase::hasVertex(std::string const &nanoporeId) const {
  return m_vertices.find(nanoporeId) != std::end(m_vertices);
}

std::shared_ptr<Vertex> GraphBase::getVertexAsSharedPtr(std::string const &nanoporeId) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeId);

  return iter != std::end(m_vertices) ? iter->second : nullptr;
}

Vertex *GraphBase::getVertex(std::string const &nanoporeId) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeId);

  return iter != std::end(m_vertices) ? iter->second.get() : nullptr;
}

std::vector<Vertex *> GraphBase::getVertices() const {
  std::vector<Vertex *> vertices;

  std::transform(m_vertices.begin(), m_vertices.end(), std::back_inserter(vertices),
                 [](const std::unordered_map<std::string, std::shared_ptr<Vertex>>::value_type &pair) {
                   return pair.second.get();
                 });
  return vertices;
}

bool GraphBase::hasEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &vertices) const {
  auto const iter = m_adjacencyList.find(vertices.first->getId());

  if (iter != std::end(m_adjacencyList)) {
    return iter->second.find(vertices.second->getId()) != std::end(iter->second);
  }

  return false;
}

Edge *GraphBase::getEdge(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &vertices) const {
  auto const outerIter = m_adjacencyList.find(vertices.first->getId());
  if (outerIter != std::end(m_adjacencyList)) {
    auto const innerIter = outerIter->second.find(vertices.second->getId());
    if (innerIter != std::end(outerIter->second)) {
      return innerIter->second;
    }
  }

  return nullptr;
}

std::vector<Edge *> GraphBase::getEdges() const {
  std::vector<Edge *> edges;

  std::transform(std::begin(m_edges), std::end(m_edges), std::back_inserter(edges),
                 [](auto const &pair) { return pair.second.get(); });

  return edges;
}

// PROTECTED CLASS METHODS

GraphBase::GraphBase(std::unordered_map<std::string, std::shared_ptr<Vertex>> &&vertices,
                     std::vector<std::shared_ptr<Edge>> &&edges, bool hasBidirectionalEdges)
    : m_vertices(std::move(vertices)) {
  for (auto &edge : edges) {
    auto pV1 = gsl::make_not_null(edge->getVertices().first);
    auto pV2 = gsl::make_not_null(edge->getVertices().second);
    if (hasVertex(pV1->getId()) && hasVertex(pV2->getId())) {
      _addEdgeInternal(std::move(edge), hasBidirectionalEdges);
    }
  }
}

void GraphBase::_addVertex(std::shared_ptr<Vertex> &&spVertex) {
  if (!spVertex) {
    throw std::runtime_error("Unexpected nullptr.");
  }

  std::unique_lock<std::shared_mutex> lck(m_mutexVertex);

  m_vertices.emplace(spVertex->getId(), std::move(spVertex));
}

void GraphBase::_deleteVertex(gsl::not_null<Vertex const *> const pVertex, bool hasBidirectionalEdges) {
  auto const iterSuccessors = m_adjacencyList.find(pVertex->getId());

  if (iterSuccessors == std::end(m_adjacencyList)) {
    return;
  }

  using successor_t          = decltype(iterSuccessors->second)::value_type;
  auto const eraseSuccessors = hasBidirectionalEdges //
                                   ? std::function{[&](successor_t const &s) {
                                       auto const iterPredecessor = m_adjacencyList.find(s.first);
                                       if (iterPredecessor != std::end(m_adjacencyList)) {
                                         iterPredecessor->second.erase(pVertex->getId());
                                       }

                                       auto const pEdge = gsl::make_not_null(s.second);
                                       _onEdgeDeleted(pEdge->getVertices());
                                       m_edges.erase(pEdge->getId());
                                     }}
                                   : std::function{[&](successor_t const &s) {
                                       auto const pEdge = gsl::make_not_null(s.second);
                                       _onEdgeDeleted(pEdge->getVertices());
                                       m_edges.erase(pEdge->getId());
                                     }};
  std::for_each(std::begin(iterSuccessors->second), std::end(iterSuccessors->second), eraseSuccessors);
  m_adjacencyList.erase(iterSuccessors);

  if (!hasBidirectionalEdges) {
    for (auto &[targetId, connectedVertices] : m_adjacencyList) {
      LB_UNUSED(targetId);

      auto const iterPredecessor = connectedVertices.find(pVertex->getId());
      if (iterPredecessor != std::end(connectedVertices)) {
        auto const pEdge = gsl::make_not_null(iterPredecessor->second);
        _onEdgeDeleted(pEdge->getVertices());
        m_edges.erase(pEdge->getId());
        connectedVertices.erase(iterPredecessor);
      }
    }
  }

  m_vertices.erase(pVertex->getId());
}

void GraphBase::_addEdge(
    std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &vertexIds,
    bool                                                                                       isBidirectional) {
  std::unique_lock<std::shared_mutex> lck(m_mutexEdge);

  auto pV1 = vertexIds.first->getSharedPtr();
  auto pV2 = vertexIds.second->getSharedPtr();

  if (!pV1 || !pV2) {
    // Edges between unknown vertices are omitted
    return;
  }

  auto const pVertices = std::make_pair(pV1.get(), pV2.get());

  auto spEdge = std::make_shared<Edge>(std::make_pair(std::move(pV1), std::move(pV2)));
  _addEdgeInternal(std::move(spEdge), isBidirectional);

  _onEdgeAdded(pVertices);
}

void GraphBase::_deleteEdge(gsl::not_null<Edge const *> const pEdge, bool isBidirectional) {
  auto const findAndDeleteEdge = [&](std::string const &source, std::string const &target) {
    auto const outerIter = m_adjacencyList.find(source);

    if (outerIter != std::end(m_adjacencyList)) {
      auto const innerIter = outerIter->second.find(target);
      if (innerIter != std::end(outerIter->second)) {
        outerIter->second.erase(innerIter->first);
      }
    }
  };
  auto const vertices = pEdge->getVertices();

  findAndDeleteEdge(vertices.first->getId(), vertices.second->getId());
  if (isBidirectional) {
    findAndDeleteEdge(vertices.second->getId(), vertices.first->getId());
  }

  m_edges.erase(pEdge->getId());

  _onEdgeDeleted(vertices);
}

std::unordered_map<std::string, Edge *const>
GraphBase::_getSuccessors(gsl::not_null<Vertex const *> const pVertex) const {
  auto const outerIter = m_adjacencyList.find(pVertex->getId());
  if (outerIter != std::end(m_adjacencyList)) {
    return outerIter->second;
  }

  return decltype(outerIter->second)();
}

std::unordered_map<std::string, Edge *const>
GraphBase::_getPredecessors(gsl::not_null<Vertex const *> const pVertex) const {
  std::unordered_map<std::string, Edge *const> result;

  auto const &vertexId = pVertex->getId();
  for (auto const &[currentVertexId, connectedVertices] : m_adjacencyList) {
    if (vertexId == currentVertexId) {
      continue;
    }

    auto const iter = connectedVertices.find(vertexId);
    if (iter != std::end(connectedVertices)) {
      result.insert(std::make_pair(currentVertexId, iter->second));
    }
  }

  return result;
}

// PRIVATE CLASS METHODS

void GraphBase::_addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional) {
  auto const emplaceEdge = [&](std::string const &source, std::string const &target, Edge *pEdge) {
    auto const iterAdjacencyList =
        m_adjacencyList.emplace(source, std::unordered_map<std::string, Edge *const>()).first;
    iterAdjacencyList->second.emplace(target, pEdge);
  };

  auto const  pVertices = spEdge->getVertices();
  auto const  iterEdges = m_edges.emplace(spEdge->getId(), std::move(spEdge)).first;
  auto *const pEdge     = iterEdges != std::end(m_edges) ? iterEdges->second.get() : nullptr;

  emplaceEdge(pVertices.first->getId(), pVertices.second->getId(), pEdge);

  if (isBidirectional) {
    emplaceEdge(pVertices.second->getId(), pVertices.first->getId(), pEdge);
  }
}

// -----------
// class Graph
// -----------

// PUBLIC CLASS METHODS

std::unique_ptr<Graph> Graph::getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::vector<std::shared_ptr<Edge>> e;
  e.reserve(getSize());

  auto const edgesUntransformed = getEdges();
  std::transform(std::begin(edgesUntransformed), std::end(edgesUntransformed), std::back_inserter(e),
                 [](auto *const pEdge) { return pEdge->getSharedPtr(); });

  return std::make_unique<Graph>(getVertexMap(vertices), std::move(e));
}

// -------------
// class DiGraph
// -------------

// PUBLIC CLASS METHODS

void DiGraph::addVertex(std::shared_ptr<Vertex> &&spVertex) {
  auto const *const pVertex = spVertex.get();

  _addVertex(std::move(spVertex));

  {
    std::scoped_lock<std::mutex> lck(m_mutexDegrees);
    m_inDegrees.insert({pVertex, 0});
    m_outDegrees.insert({pVertex, 0});
  }
}

std::unique_ptr<DiGraph> DiGraph::getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::vector<std::shared_ptr<Edge>> e;
  e.reserve(getSize());

  auto const edgesUntransformed = getEdges();
  std::transform(std::begin(edgesUntransformed), std::end(edgesUntransformed), std::back_inserter(e),
                 [](auto *const pEdge) { return pEdge->getSharedPtr(); });

  return std::make_unique<DiGraph>(getVertexMap(vertices), std::move(e));
}

// PRIVATE CLASS METHODS

void DiGraph::_increaseInDegree(Vertex const *const pVertex) {
  auto const iter = m_inDegrees.find(pVertex);
  if (iter == std::end(m_inDegrees)) {
    return;
  }

  ++iter->second;
}

void DiGraph::_decreaseInDegree(const Vertex *const pVertex) {
  auto const iter = m_inDegrees.find(pVertex);
  if (iter == std::end(m_inDegrees)) {
    return;
  }

  --iter->second;
}

void DiGraph::_increaseOutDegree(Vertex const *const pVertex) {
  auto const iter = m_outDegrees.find(pVertex);
  if (iter == std::end(m_outDegrees)) {
    return;
  }

  ++iter->second;
}

void DiGraph::_decreaseOutDegree(Vertex const *const pVertex) {
  auto const iter = m_outDegrees.find(pVertex);
  if (iter == std::end(m_outDegrees)) {
    return;
  }

  --iter->second;
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

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------