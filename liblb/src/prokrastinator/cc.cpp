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
#include <set>

#include "graph/Graph.h"
#include "types/Direction.h"

std::vector<std::vector<lazybastard::graph::Vertex *>>
lazybastard::getConnectedComponents(gsl::not_null<const lazybastard::graph::Graph *> pGraph) {
  std::vector<std::vector<lazybastard::graph::Vertex *>> result;
  std::set<lazybastard::graph::Vertex const *>           visited;

  auto const vertices = pGraph->getVertices();
  for (auto *const pVertex : vertices) {
    if (visited.contains(pVertex)) {
      continue;
    }

    auto const                                neighbors = pGraph->getNeighbors(pVertex);
    std::vector<lazybastard::graph::Vertex *> component({pVertex});
    visited.insert(pVertex);
    component.reserve(neighbors.size());
    std::transform(std::begin(neighbors), std::end(neighbors), std::back_inserter(component),
                   [&](auto const &p) { return pGraph->getVertex(p.first); });

    std::deque<lazybastard::graph::Vertex *> queue({pVertex});

    while (!queue.empty()) {
      auto *pCurrentVertex = queue.front();
      queue.pop_front();

      auto const currentNeighbors = pGraph->getNeighbors(pCurrentVertex);

      auto iterCurrentNeighbor = std::begin(currentNeighbors);
      auto biComponent         = std::back_inserter(component);
      while (iterCurrentNeighbor != std::end(currentNeighbors)) {

        pCurrentVertex = pGraph->getVertex(iterCurrentNeighbor->first);
        if (!visited.contains(pCurrentVertex)) {
          if (iterCurrentNeighbor->second->getConsensusDirection() != lazybastard::Direction::e_NONE) {
            *biComponent++ = pCurrentVertex;
            queue.push_back(pCurrentVertex);
          }
          visited.insert(pCurrentVertex);
        }

        ++iterCurrentNeighbor;
      }
    }

    result.push_back(std::move(component));
  }

  return result;
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------