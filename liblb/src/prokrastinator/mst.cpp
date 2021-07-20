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

#include "Prokrastinator.h"

#include <algorithm>
#include <memory>
#include <unordered_map>

#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "types/Direction.h"

namespace {
class UnionFind {
public:
  lazybastard::graph::Vertex const *operator[](lazybastard::graph::Vertex const *const pVertex) {
    auto const iter = m_parents.find(pVertex);
    if (iter == std::end(m_parents)) {
      m_parents.insert({pVertex, pVertex});
      m_weights.insert({pVertex, 1});

      return pVertex;
    }

    std::vector<lazybastard::graph::Vertex const *> path({pVertex});
    auto const *                                    pRoot = iter->second;
    while (pRoot != path.back()) {
      path.push_back(pRoot);
      pRoot = m_parents.find(pRoot)->second;
    }

    for (auto const *pAncestor : path) {
      m_parents[pAncestor] = pRoot;
    }

    return pRoot;
  }

  void unify(lazybastard::graph::Vertex const *const pV1, lazybastard::graph::Vertex const *const pV2) {
    auto weightedPair = std::make_pair(operator[](pV1), operator[](pV2));
    if (m_weights[pV2] > m_weights[pV1]) {
      std::swap(weightedPair.first, weightedPair.second);
    }

    m_weights[weightedPair.first] += m_weights[weightedPair.second];
    m_parents[weightedPair.second] = weightedPair.first;
  }

private:
  std::unordered_map<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *> m_parents;
  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t>                        m_weights;
};

std::vector<std::shared_ptr<lazybastard::graph::Edge>> kruskal(lazybastard::graph::Graph const &graph) {
  auto edges = graph.getEdges();
  edges.erase(std::remove_if(std::begin(edges), std::end(edges),
                             [](auto const *const pEdge) {
                               return pEdge->getConsensusDirection() == lazybastard::Direction::e_NONE;
                             }),
              edges.end());
  std::sort(std::begin(edges), std::end(edges),
            [](auto const *pEdge1, auto const *pEdge2) { return pEdge1->getWeight() > pEdge2->getWeight(); });

  std::vector<std::shared_ptr<lazybastard::graph::Edge>> result;
  auto                                                   uf = UnionFind();
  for (auto *pEdge : edges) {
    auto const vertices = pEdge->getVertices();
    if (uf[vertices.first] != uf[vertices.second]) {
      result.emplace_back(pEdge->getSharedPtr());
      uf.unify(vertices.first, vertices.second);
    }
  }

  return result;
}
} // unnamed namespace

lazybastard::graph::Graph lazybastard::getMaxSpanTree(lazybastard::graph::Graph const &graph) {
  auto vertices = [=]() {
    auto const                                                                    vertices = graph.getVertices();
    std::unordered_map<unsigned int, std::shared_ptr<lazybastard::graph::Vertex>> newVertices;

    std::transform(std::begin(vertices), std::end(vertices), std::inserter(newVertices, std::begin(newVertices)),
                   [](auto *const pVertex) { return std::make_pair(pVertex->getId(), pVertex->getSharedPtr()); });

    return newVertices;
  }();

  return lazybastard::graph::Graph(std::move(vertices), kruskal(graph));
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
