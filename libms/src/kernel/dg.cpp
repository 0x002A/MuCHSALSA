// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of MuCHSALSA <https://github.com/0x002A/MuCHSALSA>.
//
// MuCHSALSA is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// MuCHSALSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with MuCHSALSA.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#include "Kernel.h"

#include <stack>
#include <tuple>
#include <utility>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"
#include "types/Direction.h"
#include "types/Toggle.h"

muchsalsa::graph::DiGraph muchsalsa::getDirectionGraph(gsl::not_null<muchsalsa::matching::MatchMap *> pMatchMap,
                                                       muchsalsa::graph::Graph const &                graph,
                                                       muchsalsa::graph::Graph const & connectedComponent,
                                                       muchsalsa::graph::Vertex const &startNode) {
  std::stack<std::tuple<graph::Vertex const *, Toggle>> stack;
  stack.push(std::make_tuple(&startNode, true));

  auto diGraph = graph::DiGraph();
  while (!stack.empty()) {
    auto const currentNode = stack.top();
    stack.pop();

    auto const *const pCurrentNode = std::get<0>(currentNode);

    if (!diGraph.hasVertex(pCurrentNode->getId())) {
      diGraph.addVertex(graph.getVertexAsSharedPtr(pCurrentNode->getId()));
    }

    if (pCurrentNode->getVertexDirection() == Direction::e_NONE) {
      diGraph.getVertex(pCurrentNode->getId())->setVertexDirection(std::get<1>(currentNode));
    }

    auto const neighbors = connectedComponent.getNeighbors(pCurrentNode);
    for (auto const &[neighborId, pNeighborEdge] : neighbors) {
      auto const *const pOtherNode      = graph.getVertex(neighborId);
      auto const        vertices        = pNeighborEdge->getVertices();
      auto              otherNodeExists = diGraph.hasVertex(neighborId);

      util::exchange_if(otherNodeExists, pOtherNode->getVertexDirection() != Direction::e_NONE, otherNodeExists);

      if (!otherNodeExists) {
        diGraph.addVertex(connectedComponent.getVertexAsSharedPtr(neighborId));
      }

      if (diGraph.hasEdge(std::make_pair(vertices.first, vertices.second)) ||
          diGraph.hasEdge(std::make_pair(vertices.second, vertices.first))) {
        continue;
      }

      auto const *const pEdge = graph.getEdge(std::make_pair(vertices.first, vertices.second));
      for (auto const &order : pNeighborEdge->getEdgeOrders()) {
        auto flip = false;

        util::exchange_if(flip, !flip, !order.direction && order.baseVertex == pOtherNode);
        util::exchange_if(flip, !flip, !std::get<1>(currentNode));

        auto const *pStart = order.startVertex;
        auto const *pEnd   = order.endVertex;

        util::swap_if(pStart, pEnd, flip);

        auto *pNewEdge = diGraph.getEdge(std::make_pair(pStart, pEnd));
        if (!pNewEdge) {
          auto verticesNewEdge = std::make_pair(pStart, pEnd);

          diGraph.addEdge(verticesNewEdge);
          pNewEdge = diGraph.getEdge(verticesNewEdge);

          pNewEdge->setShadow(pEdge->isShadow());

          if (!pEdge->isShadow()) {
            pNewEdge->setWeight(pNeighborEdge->getWeight());
          }

          for (auto const &match : *pMatchMap->getEdgeMatches(graph.getEdge(verticesNewEdge))) {
            pMatchMap->addEdgeMatch(pNewEdge, match.first, match.second);
          }
        }

        pNewEdge->appendOrder(order);
      }

      if (pNeighborEdge->getConsensusDirection() == Direction::e_NONE) {
        continue;
      }

      auto const nextMod =
          std::get<1>(currentNode) * Toggle(pNeighborEdge->getConsensusDirection() == Direction::e_POS);

      if (!otherNodeExists) {
        stack.push(std::make_tuple(pOtherNode, nextMod));
      }
    }
  }

  return diGraph;
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
