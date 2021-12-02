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


std::vector<std::vector<muchsalsa::graph::Vertex *>>
muchsalsa::splitConnectedComponentsbyChineseWhispers(muchsalsa::graph::Graph const &graph, std::vector<muchsalsa::graph::Edge const *> * const outOfComponentEdges) {

    auto const vertices = graph.getVertices();
    std::size_t size = vertices.size();
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_int_distribution<std::size_t> unirand(0, size-1);
    
    std::unordered_map<graph::Vertex const *, std::size_t> mappingVertex2Idx;
    
    std::size_t groupCounter = 0;
    for (auto *const itVertex : vertices) {
        mappingVertex2Idx[itVertex] = groupCounter;
        groupCounter++;
    }

    for(std::size_t i = 0; i <= 50*size; i++){
    
         std::unordered_map<std::size_t, unsigned int> neighbourcount;
         auto *const randomVertex = vertices[ unirand(re) ];
         
         auto const currentNeighbors = graph.getNeighbors(randomVertex);
         for (auto iterNeighbor = std::begin(currentNeighbors); iterNeighbor != std::end(currentNeighbors); ++iterNeighbor) {
             auto *pNeighbor = graph.getVertex(iterNeighbor->first);
             
             auto *const edge = graph.getEdge(std::make_pair(randomVertex, pNeighbor));
             if (!edge->isShadow()) {
                 neighbourcount[ mappingVertex2Idx[pNeighbor] ]++;
             }
         }

         if (neighbourcount.empty()) {
             continue;
         }

         std::size_t const newGroup = std::max_element(neighbourcount.begin(), neighbourcount.end(),
                             [](const std::pair<std::size_t, unsigned int> &p1,
                                const std::pair<std::size_t, unsigned int> &p2)
                             {
                                 return p1.second < p2.second;
                             })->first;

         mappingVertex2Idx[randomVertex] = newGroup;
    }

    std::vector<std::vector<muchsalsa::graph::Vertex *>> result;
    std::unordered_map<std::size_t, std::size_t> groupToResultIdx;
    
    for (auto *const itVertex : vertices) {    
        const std::size_t group = mappingVertex2Idx[itVertex];
        if (groupToResultIdx.find(group) == groupToResultIdx.end()) {
            std::vector<muchsalsa::graph::Vertex *> newInner;
            result.push_back(std::move(newInner));
            groupToResultIdx[group] = result.size()-1;
        }

        result[groupToResultIdx[group]].push_back(itVertex);
    }
    
    
    auto const edges = graph.getEdges();
    for (auto const *const pEdge : edges) {
       if ( mappingVertex2Idx[ pEdge->getVertices().first ] != mappingVertex2Idx[ pEdge->getVertices().second ] ) {
          outOfComponentEdges->push_back(pEdge);
       }
    }
    
    return result;
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
