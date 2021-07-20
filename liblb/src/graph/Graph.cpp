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

#include <deque>
#include <stdexcept>

#include "Util.h"
#include "graph/Edge.h"

namespace lazybastard::graph {

// =====================================================================================================================
//                                                        HELPER
// =====================================================================================================================

namespace {

auto getVertexMap(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::unordered_map<unsigned int, std::shared_ptr<Vertex>> result;

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

bool GraphBase::hasVertex(unsigned int vertexId) const { return m_vertices.find(vertexId) != std::end(m_vertices); }

std::shared_ptr<Vertex> GraphBase::getVertexAsSharedPtr(unsigned int vertexId) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto const iter = m_vertices.find(vertexId);

  return iter != std::end(m_vertices) ? iter->second : nullptr;
}

Vertex *GraphBase::getVertex(unsigned int vertexId) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto const iter = m_vertices.find(vertexId);

  return iter != std::end(m_vertices) ? iter->second.get() : nullptr;
}

std::vector<Vertex *> GraphBase::getVertices() const {
  std::vector<Vertex *> vertices;

  std::transform(m_vertices.begin(), m_vertices.end(), std::back_inserter(vertices),
                 [](const auto &pair) { return pair.second.get(); });
  return vertices;
}

std::unordered_set<Vertex *> GraphBase::getVerticesAsUnorderedSet() const {
  std::unordered_set<Vertex *> vertices;

  std::transform(m_vertices.begin(), m_vertices.end(), std::inserter(vertices, std::begin(vertices)),
                 [](const auto &pair) { return pair.second.get(); });
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

GraphBase::GraphBase(std::unordered_map<unsigned int, std::shared_ptr<Vertex>> &&vertices,
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
  using successor_t          = decltype(m_adjacencyList)::mapped_type::value_type;
  auto const eraseSuccessors = hasBidirectionalEdges //
                                   ? std::function{[&](successor_t const &s) {
                                       auto const iterPredecessor = m_adjacencyList.find(s.first);
                                       if (iterPredecessor != std::end(m_adjacencyList)) {
                                         iterPredecessor->second.erase(pVertex->getId());
                                       }

                                       auto const pEdge = gsl::make_not_null(s.second);

                                       std::for_each(std::begin(m_observers), std::end(m_observers),
                                                     [=](auto *pObserver) { pObserver->onEdgeDeleted(pEdge); });

                                       auto vertices = pEdge->getVertices();
                                       m_edges.erase(pEdge);
                                       _onEdgeDeleted(std::move(vertices));
                                     }}
                                   : std::function{[&](successor_t const &s) {
                                       auto const pEdge = gsl::make_not_null(s.second);

                                       std::for_each(std::begin(m_observers), std::end(m_observers),
                                                     [=](auto *pObserver) { pObserver->onEdgeDeleted(pEdge); });

                                       auto vertices = pEdge->getVertices();
                                       m_edges.erase(pEdge);
                                       _onEdgeDeleted(std::move(vertices));
                                     }};

  auto const iterSuccessors = m_adjacencyList.find(pVertex->getId());
  if (iterSuccessors != std::end(m_adjacencyList)) {
    std::for_each(std::begin(iterSuccessors->second), std::end(iterSuccessors->second), eraseSuccessors);
    m_adjacencyList.erase(iterSuccessors);
  }

  if (!hasBidirectionalEdges) {
    for (auto &[targetId, connectedVertices] : m_adjacencyList) {
      LB_UNUSED(targetId);

      auto const iterPredecessor = connectedVertices.find(pVertex->getId());
      if (iterPredecessor != std::end(connectedVertices)) {
        auto const pEdge = gsl::make_not_null(iterPredecessor->second);

        std::for_each(std::begin(m_observers), std::end(m_observers),
                      [=](auto *pObserver) { pObserver->onEdgeDeleted(pEdge); });

        auto vertices = pEdge->getVertices();
        m_edges.erase(pEdge);
        _onEdgeDeleted(std::move(vertices));

        connectedVertices.erase(iterPredecessor);
      }
    }
  }

  std::for_each(std::begin(m_observers), std::end(m_observers),
                [=](auto *pObserver) { pObserver->onVertexDeleted(pVertex); });

  m_vertices.erase(pVertex->getId());
}

void GraphBase::_addEdge(
    std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const &pVertices,
    bool                                                                                       isBidirectional) {
  std::unique_lock<std::shared_mutex> lck(m_mutexEdge);

  auto spV1 = pVertices.first->getSharedPtr();
  auto spV2 = pVertices.second->getSharedPtr();

  if (!hasVertex(spV1->getId()) || !hasVertex(spV2->getId())) {
    // Edges between unknown vertices are omitted
    return;
  }

  auto spEdge = std::make_shared<Edge>(std::make_pair(std::move(spV1), std::move(spV2)));

  if (_addEdgeInternal(std::move(spEdge), isBidirectional)) {
    _onEdgeAdded(pVertices);
  }
}

void GraphBase::_deleteEdge(gsl::not_null<Edge const *> const pEdge, bool isBidirectional) {
  auto const findAndDeleteEdge = [&](unsigned int source, unsigned int target) {
    auto const outerIter = m_adjacencyList.find(source);

    if (outerIter != std::end(m_adjacencyList)) {
      auto const innerIter = outerIter->second.find(target);
      if (innerIter != std::end(outerIter->second)) {
        auto const isErased = outerIter->second.erase(innerIter->first);

        if (isErased) {
          _onEdgeDeleted(std::make_pair(getVertex(source), getVertex(target)));
        }
      }
    }
  };
  auto const vertices = pEdge->getVertices();

  findAndDeleteEdge(vertices.first->getId(), vertices.second->getId());
  if (isBidirectional) {
    findAndDeleteEdge(vertices.second->getId(), vertices.first->getId());
  }

  std::for_each(std::begin(m_observers), std::end(m_observers),
                [=](auto *pObserver) { pObserver->onEdgeDeleted(pEdge); });

  m_edges.erase(pEdge);
}

std::unordered_map<unsigned int, Edge *const>
GraphBase::_getSuccessors(gsl::not_null<Vertex const *> const pVertex) const {
  auto const outerIter = m_adjacencyList.find(pVertex->getId());
  if (outerIter != std::end(m_adjacencyList)) {
    return outerIter->second;
  }

  return decltype(outerIter->second)();
}

std::unordered_map<unsigned int, Edge *const>
GraphBase::_getPredecessors(gsl::not_null<Vertex const *> const pVertex) const {
  std::unordered_map<unsigned int, Edge *const> result;

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

bool GraphBase::_addEdgeInternal(std::shared_ptr<Edge> &&spEdge, bool isBidirectional) {
  auto const emplaceEdge = [&](unsigned int source, unsigned int target, Edge *pEdge) {
    auto const insertAdjacencyList = m_adjacencyList.emplace(source, decltype(m_adjacencyList)::mapped_type()).first;
    auto const inserted            = insertAdjacencyList->second.emplace(target, pEdge).second;

    return inserted;
  };

  auto const pVertices = spEdge->getVertices();
  auto const inserted  = emplaceEdge(pVertices.first->getId(), pVertices.second->getId(), spEdge.get());

  if (isBidirectional) {
    emplaceEdge(pVertices.second->getId(), pVertices.first->getId(), spEdge.get());
  }

  if (inserted) {
    m_edges.emplace(spEdge.get(), std::move(spEdge));
  }

  return inserted;
}

// -----------
// class Graph
// -----------

// PUBLIC CLASS METHODS

Graph Graph::getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::vector<std::shared_ptr<Edge>> e;
  e.reserve(getSize());

  auto const edgesUntransformed = getEdges();
  std::transform(std::begin(edgesUntransformed), std::end(edgesUntransformed), std::back_inserter(e),
                 [](auto *const pEdge) { return pEdge->getSharedPtr(); });

  return Graph(getVertexMap(vertices), std::move(e));
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

DiGraph DiGraph::getSubgraph(std::vector<lazybastard::graph::Vertex *> const &vertices) {
  std::vector<std::shared_ptr<Edge>> e;
  e.reserve(getSize());

  auto const edgesUntransformed = getEdges();
  std::transform(std::begin(edgesUntransformed), std::end(edgesUntransformed), std::back_inserter(e),
                 [](auto *const pEdge) { return pEdge->getSharedPtr(); });

  return DiGraph(getVertexMap(vertices), std::move(e));
}

std::vector<lazybastard::graph::Vertex const *> DiGraph::sortTopologically() const {
  std::vector<lazybastard::graph::Vertex const *> result;

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> verticesWithNonNullInDegree;
  std::deque<lazybastard::graph::Vertex const *>                      verticesWithNullInDegree;

  for (auto const &[pVertex, inDegree] : getInDegrees()) {
    if (inDegree > 0) {
      verticesWithNonNullInDegree[pVertex] = inDegree;
    } else {
      verticesWithNullInDegree.push_back(pVertex);
    }
  }

  while (!verticesWithNullInDegree.empty()) {
    auto const *const pVertex = verticesWithNullInDegree.back();
    verticesWithNullInDegree.pop_back();

    auto const successors = getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      LB_UNUSED(pEdge);

      auto const *pSuccessor = getVertex(targetId);

      verticesWithNonNullInDegree[pSuccessor] -= 1;

      if (verticesWithNonNullInDegree[pSuccessor] == 0) {
        verticesWithNullInDegree.push_back(pSuccessor);
        verticesWithNonNullInDegree.erase(pSuccessor);
      }
    }

    result.push_back(pVertex);
  }

  return result;
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
