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

#include <deque>
#include <functional>
#include <iterator>
#include <random>
#include <set>
#include <utility>

#include "graph/Graph.h"
#include "types/Direction.h"

std::vector<std::vector<muchsalsa::graph::Vertex *>>
muchsalsa::getConnectedComponents(muchsalsa::graph::Graph const &graph) {
  std::vector<std::vector<muchsalsa::graph::Vertex *>> result;
  std::set<muchsalsa::graph::Vertex const *>           visited;

  auto const vertices = graph.getVertices();
  for (auto *const pSourceVertex : vertices) {
    if (visited.contains(pSourceVertex)) {
      continue;
    }

    std::vector<muchsalsa::graph::Vertex *>           component({pSourceVertex});
    std::deque<muchsalsa::graph::Vertex const *const> queue({pSourceVertex});

    visited.insert(pSourceVertex);

    while (!queue.empty()) {
      auto const *const pCurrentVertex = queue.front();
      queue.pop_front();

      auto const currentNeighbors = graph.getNeighbors(pCurrentVertex);
      for (auto iterNeighbor = std::begin(currentNeighbors); iterNeighbor != std::end(currentNeighbors);
           ++iterNeighbor) {
        auto *pNeighbor = graph.getVertex(iterNeighbor->first);
        if (!visited.contains(pNeighbor) &&
            iterNeighbor->second->getConsensusDirection() != muchsalsa::Direction::e_NONE) {
          component.push_back(pNeighbor);
          queue.push_back(pNeighbor);
          visited.insert(pNeighbor);
        }
      }
    }

    result.push_back(std::move(component));
  }

  return result;
}


// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
