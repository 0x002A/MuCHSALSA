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

#include <stack>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "types/Direction.h"
#include "types/Toggle.h"

std::unique_ptr<lazybastard::graph::DiGraph>
lazybastard::getDirectionGraph(gsl::not_null<graph::Graph const *> const               pGraph,
                               gsl::not_null<graph::Graph const *> const               pConnectedComponent,
                               gsl::not_null<lazybastard::graph::Vertex const *> const pStartNode) {
  std::stack<std::tuple<graph::Vertex const *, Toggle>> stack;
  stack.push(std::make_tuple(pStartNode, true));

  auto diGraph = std::make_unique<graph::DiGraph>();
  while (!stack.empty()) {
    auto const currentNode = stack.top();
    stack.pop();

    auto const *const pCurrentNode = std::get<0>(currentNode);

    if (!diGraph->hasVertex(pCurrentNode->getId())) {
      diGraph->addVertex(pGraph->getVertexAsSharedPtr(pCurrentNode->getId()));
    }

    if (pCurrentNode->getVertexDirection() == Direction::e_NONE) {
      diGraph->getVertex(pCurrentNode->getId())->setVertexDirection(std::get<1>(currentNode));
    }

    auto const neighbors = pConnectedComponent->getNeighbors(pCurrentNode);
    for (auto const &[neighborId, pNeighborEdge] : neighbors) {
      auto const *const pOtherNode      = pGraph->getVertex(neighborId);
      auto const        vertices        = pNeighborEdge->getVertices();
      auto              otherNodeExists = diGraph->hasVertex(neighborId);

      util::exchange_if(otherNodeExists, pOtherNode->getVertexDirection() != Direction::e_NONE, otherNodeExists);

      if (!otherNodeExists) {
        diGraph->addVertex(pConnectedComponent->getVertexAsSharedPtr(neighborId));
      }

      if (diGraph->hasEdge(std::make_pair(vertices.first, vertices.second)) ||
          diGraph->hasEdge(std::make_pair(vertices.second, vertices.first))) {
        continue;
      }

      auto const *const pEdge = pGraph->getEdge(std::make_pair(vertices.first, vertices.second));
      for (auto const &order : pNeighborEdge->getEdgeOrders()) {
        auto flip = false;

        util::exchange_if(flip, !flip, !order.direction && order.baseVertex == pOtherNode);
        util::exchange_if(flip, !flip, !std::get<1>(currentNode));

        auto const *pStart = order.startVertex;
        auto const *pEnd   = order.endVertex;

        util::swap_if(pStart, pEnd, flip);

        auto *pNewEdge = diGraph->getEdge(std::make_pair(pStart, pEnd));
        if (!pNewEdge) {
          diGraph->addMissingVertex(pGraph->getVertexAsSharedPtr(pStart->getId()));
          diGraph->addMissingVertex(pGraph->getVertexAsSharedPtr(pEnd->getId()));
          diGraph->addEdge(std::make_pair(pStart, pEnd));
          pNewEdge = diGraph->getEdge(std::make_pair(pStart, pEnd));

          pNewEdge->setShadow(pEdge->isShadow());

          if (!pEdge->isShadow()) {
            pNewEdge->setWeight(pNeighborEdge->getWeight());
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