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

#include <algorithm>
#include <deque>
#include <initializer_list>
#include <iterator>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

namespace {

template <class IT_HAYSTACK, class IT_NEEDLE>
bool contains(IT_HAYSTACK first1, IT_HAYSTACK last1, IT_NEEDLE first2, IT_NEEDLE last2) {
  std::unordered_set<typename IT_HAYSTACK::value_type> haystack(first1, last1);

  for (; first2 != last2; ++first2) {
    if (!haystack.contains(*first2)) {
      return false;
    }
  }

  return true;
}

template <class INPUT_IT_NON_NULL_IN_DEGREES, class INPUT_IT_NULL_IN_DEGREES>
void initializeInDegreeMaps(INPUT_IT_NON_NULL_IN_DEGREES nonNullInDegrees, INPUT_IT_NULL_IN_DEGREES nullInDegrees,
                            gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraphCycle) {
  auto const inDegrees = pDiGraphCycle->getInDegrees();

  auto inDegreesIter = std::begin(inDegrees);
  while (inDegreesIter != std::end(inDegrees)) {
    if (inDegreesIter->second > 0) {
      *nonNullInDegrees++ = *inDegreesIter;
    } else {
      *nullInDegrees++ = inDegreesIter->first;
    }

    ++inDegreesIter;
  }
}


void findClusterWeightsHeuristic(gsl::not_null<std::unordered_map<muchsalsa::graph::Edge const *, std::size_t> *> const pResult,
                        gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraphCycle) {

  pResult->clear();
  auto const sortedVertices = pDiGraphCycle->sortTopologically();

  std::unordered_map<muchsalsa::graph::Vertex const *, std::size_t> mappingVertexToIdx;

  std::size_t idx = 0;
  for (auto const *const pVertex : sortedVertices) {
    mappingVertexToIdx.insert({pVertex, idx});
    ++idx;
  }

  auto const edges = pDiGraphCycle->getEdges();
  std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) { pResult->emplace(pEdge, 0); });

  for (auto const *const pVertex : sortedVertices) {
      std::set<std::size_t > sortedSuccessors;

      auto const successors    = pDiGraphCycle->getSuccessors(pVertex);
      for (auto const &[targetId, pEdge] : successors) {
         MS_UNUSED(pEdge);
         sortedSuccessors.insert( mappingVertexToIdx.at(pDiGraphCycle->getVertex(targetId)) );
      }
                           
      std::unordered_map< muchsalsa::graph::Vertex const *, std::vector<std::size_t > > candidates;
      candidates.insert( {pVertex,  { mappingVertexToIdx.at(pVertex) } } );
      for (const auto successorID: sortedSuccessors) {

            std::vector<std::size_t > bestPath;
            auto const *const v   = sortedVertices.at(successorID);
          
			auto const predecessors    = pDiGraphCycle->getPredecessors(v);
			for (auto const &[targetId, pEdge] : predecessors) {
				MS_UNUSED(pEdge);

                 auto const *const preV = pDiGraphCycle->getVertex(targetId);
                 auto const  preVCandidate = candidates.find(preV);
                 if ( preVCandidate != std::end(candidates) ) {
                      if (preVCandidate->second.size() > bestPath.size()) {
                           bestPath = preVCandidate->second;
                      }
                 }
			}

            bestPath.push_back(mappingVertexToIdx.at(v));
            candidates.insert( {v, bestPath} );
      }

      auto const bestPathofIDs = std::max_element(candidates.begin(), candidates.end(),
                             [](const std::pair<muchsalsa::graph::Vertex const *, std::vector<std::size_t > > &p1,
                                const std::pair<muchsalsa::graph::Vertex const *, std::vector<std::size_t > > &p2)
                             {
                                 return p1.second.size() < p2.second.size();
                             })->second;

      auto c = bestPathofIDs.size() - 1;
      auto const limit = std::max(bestPathofIDs.size(), 1UL) - 1;

      for (std::size_t i = 0; i < limit; ++i) {
        auto const *const pV1   = sortedVertices.at(bestPathofIDs.at(i));
        auto const *const pV2   = sortedVertices.at(bestPathofIDs.at(i + 1));
        auto const *const pEdge = pDiGraphCycle->getEdge(std::make_pair(pV1, pV2));

        (*pResult)[pEdge] += c;
        c -= 1;
      }
  }
}


void findClusterWeights(gsl::not_null<std::unordered_map<muchsalsa::graph::Edge const *, std::size_t> *> const pResult,
                        gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraphCycle) {
                        
  pResult->clear();

  auto const sortedVertices = pDiGraphCycle->sortTopologically();

  std::unordered_map<muchsalsa::graph::Vertex const *, std::size_t> mappingVertexToIdx;

  std::size_t idx = 0;
  for (auto const *const pVertex : sortedVertices) {
    mappingVertexToIdx.insert({pVertex, idx});

    ++idx;
  }

  auto const edges = pDiGraphCycle->getEdges();
  std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) { pResult->emplace(pEdge, 0); });

  std::unordered_map<muchsalsa::graph::Vertex const *, std::set<std::size_t>> mappingPredecessors;
  std::unordered_map<muchsalsa::graph::Vertex const *, std::set<std::size_t>> mappingSuccessors;

  std::for_each(std::begin(sortedVertices), std::end(sortedVertices), [&](auto const *const pVertex) {
    auto &     idxSuccessors = mappingSuccessors[pVertex]; // NOLINT
    auto const successors    = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      MS_UNUSED(pEdge);

      idxSuccessors.insert(mappingVertexToIdx.at(pDiGraphCycle->getVertex(targetId)));
    }

    auto &     idxPredecessors = mappingPredecessors[pVertex]; // NOLINT
    auto const predecessors    = pDiGraphCycle->getPredecessors(pVertex);
    for (auto const &[targetId, pEdge] : predecessors) {
      MS_UNUSED(pEdge);

      idxPredecessors.insert(mappingVertexToIdx.at(pDiGraphCycle->getVertex(targetId)));
    }
  });

  for (auto const *const pVertex : sortedVertices) {
    std::vector<
        std::tuple<decltype(mappingPredecessors)::mapped_type, std::vector<decltype(mappingVertexToIdx)::mapped_type>>>
        candidates({{mappingSuccessors.at(pVertex), {mappingVertexToIdx.at(pVertex)}}});

    for (auto const idxOut : mappingSuccessors.at(pVertex)) {
      auto const *const pActiveVertex = sortedVertices.at(idxOut);

      for (auto const idxIn : mappingPredecessors.at(pActiveVertex)) {
        for (std::size_t idxCandidate = 0; idxCandidate < candidates.size(); ++idxCandidate) {
          auto const &candidate = candidates.at(idxCandidate);
          auto const &open      = std::get<0>(candidate);
          auto const &visited   = std::get<1>(candidate);

          if (visited.back() == idxIn && open.contains(idxOut)) {
            decltype(mappingSuccessors)::mapped_type intersection;
            std::set_intersection(std::begin(open), std::end(open), std::begin(mappingSuccessors.at(pActiveVertex)),
                                  std::end(mappingSuccessors.at(pActiveVertex)),
                                  std::inserter(intersection, std::end(intersection)));
            auto newVisited = visited;
            newVisited.push_back(idxOut);

            candidates.emplace_back(std::move(intersection), std::move(newVisited));
          }
        }
      }

      decltype(candidates) filtered;
      filtered.reserve(candidates.size());
      for (std::size_t idxOuterCandidate = 0; idxOuterCandidate < candidates.size(); ++idxOuterCandidate) {
        auto const &outerCandidate = candidates.at(idxOuterCandidate);
        auto        isDominated    = false;

        for (std::size_t idxInnerCandidate = 0; idxInnerCandidate < candidates.size(); ++idxInnerCandidate) {
          auto const &innerCandidate = candidates.at(idxInnerCandidate);
          if (idxOuterCandidate != idxInnerCandidate &&
              contains(std::begin(std::get<0>(innerCandidate)), std::end(std::get<0>(innerCandidate)),
                       std::begin(std::get<0>(outerCandidate)), std::end(std::get<0>(outerCandidate))) &&
              contains(std::begin(std::get<1>(innerCandidate)), std::end(std::get<1>(innerCandidate)),
                       std::begin(std::get<1>(outerCandidate)), std::end(std::get<1>(outerCandidate)))) {
            isDominated = true;
            break;
          }
        }

        if (!isDominated) {
          filtered.push_back(outerCandidate);
        }
      }
      filtered.shrink_to_fit();
      candidates = std::move(filtered);
    }

    std::vector<std::vector<std::size_t>> maxVisited;
    std::size_t                           maxLength = 0;
    for (auto const &[open, visited] : candidates) {
      MS_UNUSED(open);

      if (visited.size() > maxLength) {
        maxVisited = {visited};
        maxLength  = visited.size();
      } else if (visited.size() == maxLength) {
        maxVisited.push_back(visited);
      }
    }

    for (auto const &mv : maxVisited) {
      auto c = mv.size() - 1;

      auto const limit = std::max(mv.size(), 1UL) - 1;
      for (std::size_t i = 0; i < limit; ++i) {
        auto const *const pV1   = sortedVertices.at(mv.at(i));
        auto const *const pV2   = sortedVertices.at(mv.at(i + 1));
        auto const *const pEdge = pDiGraphCycle->getEdge(std::make_pair(pV1, pV2));

        (*pResult)[pEdge] += c;
        c -= 1;
      }
    }
  }
}


std::vector<muchsalsa::graph::Vertex const *> findConservationPathAlt(
    gsl::not_null<muchsalsa::graph::DiGraph *> const                                             pDiGraphCycle,
    gsl::not_null<std::unordered_map<muchsalsa::graph::Edge const *, std::size_t> const *> const pClusterWeights) {
       
  auto const sortedVertices = pDiGraphCycle->sortTopologically();

  std::unordered_map<muchsalsa::graph::Vertex const *, std::size_t> mappingVertexToIdx;

  for (std::size_t idx = 0; idx < sortedVertices.size(); ++idx) {
    mappingVertexToIdx.insert({sortedVertices[idx], idx});
  }

  std::vector<muchsalsa::graph::Vertex const *> finalizedPath;
  
  std::unordered_map<muchsalsa::graph::Vertex const *, std::tuple<std::size_t, std::vector<muchsalsa::graph::Vertex const *> > > openPaths;

  auto const &outDegrees = pDiGraphCycle->getOutDegrees();
  for (auto const *const pVertex : sortedVertices) {
    if (outDegrees.at(pVertex) == 0) {
          if (openPaths.find(pVertex) == openPaths.end()) {
              if (finalizedPath.empty()) {
                  finalizedPath = { pVertex };
              }
          } else {   
              if (std::get<1>(openPaths[pVertex]).size() > finalizedPath.size() ) {
                  finalizedPath = std::move(std::get<1>(openPaths[pVertex]));
              } else {
                  std::get<1>(openPaths[pVertex]).clear();
              }       
          }     
      continue;
    }

    std::vector<std::pair<muchsalsa::graph::Vertex const *, muchsalsa::graph::Vertex const *>> maxOuts;
    std::size_t                                                                                maxOut = 0;
    auto const successors = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      auto edge = pEdge->getVertices();
      muchsalsa::util::swap_if(edge.first, edge.second, edge.second->getId() != targetId);

      auto const currentClusterWeight = pClusterWeights->at(pEdge);
      if (currentClusterWeight > maxOut) {
        maxOut  = currentClusterWeight;
        maxOuts = {std::move(edge)};
      } else if (currentClusterWeight == maxOut) {
        maxOuts.push_back(std::move(edge));
      }
    }

    for (auto const &edge : maxOuts) {
      const auto *const pNext = edge.second;
      if (openPaths.find(pNext) != openPaths.end()) {
      
        if (std::get<0>(openPaths[pNext]) < maxOut  ||  ( std::get<0>(openPaths[pNext]) == maxOut  &&  std::get<1>(openPaths[pNext]).size() <  std::get<1>(openPaths[pVertex]).size() + 1 ) ) {
            auto tmp = std::get<1>(openPaths[pVertex]);
            tmp.push_back(edge.second);
            openPaths[pNext] = std::make_tuple(maxOut, std::move(tmp));
        }
               
      } else {
      
         if (openPaths.find(pVertex) != openPaths.end()) {
            auto tmp = std::get<1>(openPaths[pVertex]);
            tmp.push_back(edge.second);
            openPaths[pNext] = std::make_tuple(maxOut, std::move(tmp));
         } else {
            openPaths[pNext] = std::make_tuple(maxOut,
                               std::initializer_list<muchsalsa::graph::Vertex const *>({edge.first, edge.second}));
         }
         
      }
    }
    std::get<1>(openPaths[pVertex]).clear();
    
  }

  return finalizedPath;
}


std::vector<std::vector<muchsalsa::graph::Vertex const *>>
extractPaths(gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraph) {

  auto diGraphCycle = *pDiGraph;

  auto const edges = diGraphCycle.getEdges();
  for (auto const *const pEdge : edges) {
    if (pEdge->isShadow()) {
      diGraphCycle.deleteEdge(pEdge);
    }
  }

  muchsalsa::sortReductionByWeight(&diGraphCycle);

  std::unordered_map<muchsalsa::graph::Edge const *, std::size_t> clusterWeights;
  if (diGraphCycle.getOrder() < 150000) {
     findClusterWeights(&clusterWeights, &diGraphCycle);
  } else {
     findClusterWeightsHeuristic(&clusterWeights, &diGraphCycle);
  }
  
  std::vector<std::vector<muchsalsa::graph::Vertex const *>> paths;
  std::unordered_set<muchsalsa::graph::Vertex const *>       visited;

  while (diGraphCycle.getSize() > 0) {
 
    auto const longestPath = findConservationPathAlt(&diGraphCycle, &clusterWeights);
    
    if (longestPath.size() < 10) {
      auto isInVisit = false;

      auto const predecessors = pDiGraph->getPredecessors(longestPath.front());

      for (auto const &[targetId, pEdge] : predecessors) {
        MS_UNUSED(pEdge);

        auto const *const pPredecessor = pDiGraph->getVertex(targetId);
        muchsalsa::util::exchange_if(isInVisit, true, visited.contains(pPredecessor));
      }

      auto       isOutVisit = false;
      auto const successors = pDiGraph->getSuccessors(longestPath.back());
      for (auto const &[targetId, pEdge] : successors) {
        MS_UNUSED(pEdge);

        auto const *const pSuccessor = pDiGraph->getVertex(targetId);
        muchsalsa::util::exchange_if(isOutVisit, true, visited.contains(pSuccessor));
      }

      if ((!isInVisit && !isOutVisit) || ((isInVisit || isOutVisit) && longestPath.size() > 5)) {
        paths.push_back(longestPath);
      }
    } else {
      paths.push_back(longestPath);
    }

    for (auto const *const pVertex : longestPath) {
      visited.insert(pVertex);
      diGraphCycle.deleteVertex(pVertex);
    }
  }

  for (auto const *const pVertex : diGraphCycle.getVertices()) {
    paths.push_back({pVertex});
  }

  return paths;
}

} // unnamed namespace

void muchsalsa::sortReductionByWeight(gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraphCycle) {
  std::map<muchsalsa::graph::Vertex const *, std::size_t> verticesWithNonNullInDegree;
  std::deque<muchsalsa::graph::Vertex const *>            verticesWithNullInDegree;
  std::unordered_set<muchsalsa::graph::Vertex const *>            resolvedVertices;

  initializeInDegreeMaps(std::inserter(verticesWithNonNullInDegree, std::begin(verticesWithNonNullInDegree)),
                         std::back_inserter(verticesWithNullInDegree), pDiGraphCycle);

  std::unordered_set<muchsalsa::graph::Vertex const *> neighbors;
  if (!verticesWithNonNullInDegree.empty()) {
    neighbors.insert(std::begin(verticesWithNonNullInDegree)->first);
  }

  while (true) {

    while (!verticesWithNullInDegree.empty()) {
 
      auto const *const pVertex = verticesWithNullInDegree.front();
      verticesWithNullInDegree.pop_front();

      resolvedVertices.insert(pVertex);

      auto const successors = pDiGraphCycle->getSuccessors(pVertex);
      for (auto const &[idSuccessor, pEdge] : successors) {
        MS_UNUSED(pEdge);

        auto const *const pSuccessor = pDiGraphCycle->getVertex(idSuccessor);
        verticesWithNonNullInDegree.at(pSuccessor) -= 1;

        if (verticesWithNonNullInDegree.at(pSuccessor) == 0) {
          verticesWithNullInDegree.push_back(pSuccessor);
          verticesWithNonNullInDegree.erase(pSuccessor);

          neighbors.erase(pSuccessor);
        } else {
          neighbors.insert(pSuccessor);
        }
      }
    }
    
    if (verticesWithNonNullInDegree.empty()) {
      break;
    }

    muchsalsa::graph::Edge *pMinEdge = nullptr;
    muchsalsa::graph::Vertex const *pMinVertex = nullptr;
    std::size_t minScore   = 0;    
    if (neighbors.empty()) { 
    
        for ( const auto& [key, value] : verticesWithNonNullInDegree) {

          const auto *const openVertex = key;
          auto const predecessors = pDiGraphCycle->getPredecessors(openVertex);

          std::for_each(std::begin(predecessors), std::end(predecessors), [&](auto const &p) {
 
             auto * const prev = pDiGraphCycle->getVertex(p.first);
             if (resolvedVertices.find(prev) == resolvedVertices.end() ) {

		         auto  score = p.second->getWeight();
		         if (!pMinEdge || score < minScore) {
		            pMinEdge = p.second;
		            pMinVertex = openVertex;
		            minScore   = score;
		         }
             }
          });
        }    
    } else {
        for (auto const *const pNeighbor : neighbors) {
          auto const predecessors = pDiGraphCycle->getPredecessors(pNeighbor);
          std::for_each(std::begin(predecessors), std::end(predecessors), [&](auto const &p) {

             auto * const prev = pDiGraphCycle->getVertex(p.first);
             if (resolvedVertices.find(prev) == resolvedVertices.end() ) {

		         auto  score = p.second->getWeight();
		         if (!pMinEdge || score < minScore) {
		            pMinEdge = p.second;
		            pMinVertex = pNeighbor;
		            minScore   = score;
		         }
              }
          });
        }
    }
    
    pMinEdge->setShadow(true);

    pDiGraphCycle->deleteEdge(pMinEdge);

    verticesWithNonNullInDegree.at(pMinVertex) -= 1;

    if ( verticesWithNonNullInDegree.at(pMinVertex) == 0 ) {

		verticesWithNonNullInDegree.erase(pMinVertex);
		verticesWithNullInDegree.push_back(pMinVertex);
		neighbors.erase(pMinVertex);
    }

  }
  
}

std::vector<std::vector<muchsalsa::graph::Vertex const *>>
muchsalsa::linearizeGraph(gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraph) {

  auto paths = extractPaths(pDiGraph);

  std::unordered_map<std::size_t, std::size_t>           colorCorrection;
  std::unordered_map<std::size_t, std::size_t>           mappingColor2Length;
  std::unordered_map<graph::Vertex const *, std::size_t> mappingVertex2Idx;
  for (auto iter = std::begin(paths); iter != std::end(paths); ++iter) {
    auto const idx = std::distance(std::begin(paths), iter);

    for (auto const *const pVertex : *iter) {
      mappingVertex2Idx.insert({pVertex, idx});
    }

    colorCorrection.insert({idx, idx});
    mappingColor2Length.insert({idx, iter->size()});
  }

  std::vector<std::tuple<std::size_t, graph::Edge const *>> potentialJoins;
  for (auto const *const pEdge : pDiGraph->getEdges()) {
    if (pEdge->isShadow()) {
      auto const vertices = pEdge->getVertices();

      if (!mappingVertex2Idx.contains(vertices.first) || !mappingVertex2Idx.contains(vertices.second)) {
        continue;
      }

      auto const idx1 = mappingVertex2Idx[vertices.first];
      auto const idx2 = mappingVertex2Idx[vertices.second];

      auto const idxL1Start = static_cast<unsigned>(std::abs(std::distance(
          std::begin(paths[idx1]), std::find(std::begin(paths[idx1]), std::end(paths[idx1]), vertices.first))));
      auto const idxL2Start = static_cast<unsigned>(std::abs(std::distance(
          std::begin(paths[idx2]), std::find(std::begin(paths[idx2]), std::end(paths[idx2]), vertices.second))));
      auto const l1End      = mappingColor2Length[idx1] - idxL1Start - 1;
      auto const l2End      = mappingColor2Length[idx2] - idxL2Start - 1;

      if (idx1 != idx2 && l1End < idxL1Start && idxL2Start < l2End) {
        potentialJoins.emplace_back(l1End + idxL2Start, pEdge);
      }
    }
  }
  std::sort(std::begin(potentialJoins), std::end(potentialJoins));

  for (auto const &potentialJoin : potentialJoins) {
    auto const distance = std::get<0>(potentialJoin);

    if (distance > 3) {
      break;
    }

    auto const vertices = std::get<1>(potentialJoin)->getVertices();
    auto const idx1     = mappingVertex2Idx[vertices.first];
    auto const idx2     = mappingVertex2Idx[vertices.second];

    auto const findColorCorrectionIdx = [&](std::size_t idx) {
      while (colorCorrection[idx] != idx) {
        idx = colorCorrection[idx];
      }

      return idx;
    };
    auto const color1 = findColorCorrectionIdx(idx1);
    auto const color2 = findColorCorrectionIdx(idx2);
    if (color1 == color2) {
      continue;
    }

    std::size_t idxL1 = 0;
    std::size_t idxL2 = 0;

    auto const findVertexInPath = [&](std::size_t *pIdx, std::size_t const color, graph::Vertex const *const pVertex) {
      auto const &haystack = paths[color];
      auto const  iter     = std::find(std::begin(haystack), std::end(haystack), pVertex);

      *pIdx = static_cast<unsigned>(std::abs(std::distance(std::begin(haystack), iter)));

      return iter != std::end(haystack);
    };

    if (!findVertexInPath(&idxL1, color1, vertices.first) || !findVertexInPath(&idxL2, color2, vertices.second)) {
      continue;
    }

    auto const l1End   = mappingColor2Length[color1] - idxL1 - 1;
    auto const l2Start = idxL2;

    if (l1End + l2Start != distance) {
      continue;
    }

    auto &pathColor1 = paths[color1];
    auto &pathColor2 = paths[color2];
    pathColor1.erase(std::next(std::begin(pathColor1), static_cast<signed>(idxL1 + 1)), std::end(pathColor1));
    pathColor1.insert(std::end(pathColor1), std::next(std::begin(pathColor2), static_cast<signed>(idxL2)),
                      std::end(pathColor2));
    pathColor2.clear();

    colorCorrection[color2]     = colorCorrection[color1];
    mappingColor2Length[color1] = paths[color1].size();
    mappingColor2Length[color2] = 0;
  }

  std::erase_if(paths, [](auto const &path) { return path.size() <= 1; });

  return paths;
}


std::vector<muchsalsa::graph::Vertex const *> findLongestPath(
    gsl::not_null<muchsalsa::graph::DiGraph *> const pDiGraph ) {
    
  auto const sortedVertices = pDiGraph->sortTopologically();
  std::unordered_map<muchsalsa::graph::Vertex const *, std::size_t> mappingVertexToIdx;
  for (std::size_t idx = 0; idx < sortedVertices.size(); ++idx) {
    mappingVertexToIdx.insert({sortedVertices[idx], idx});
  }

  std::vector<std::size_t> dist(sortedVertices.size());
    
  for (auto const *const pVertex : sortedVertices) {

      std::size_t longestDist = 0;  

      auto const predecessors = pDiGraph->getPredecessors(pVertex);
      for (auto const &[targetId, pEdge] : predecessors) {
        MS_UNUSED(pEdge);

        auto const *const pPredecessor = pDiGraph->getVertex(targetId);
        auto const cdist = dist.at( mappingVertexToIdx[pPredecessor] );
        if (cdist + 1 > longestDist) {
            longestDist = cdist + 1;
        }
      }

      dist[mappingVertexToIdx[pVertex]] = longestDist;
  }
  
  
  std::size_t maxDist = 0;
  std::size_t maxIndex = 0;
  for(std::size_t idx = 0; idx < dist.size(); ++idx){
      if ( dist[idx] > maxDist ) {
          maxDist = dist[idx];
          maxIndex = idx;
      } 
  }

  std::deque<muchsalsa::graph::Vertex const *> path;
  const auto *pVertex = sortedVertices[maxIndex];

  while(true) {
            
      path.push_front(pVertex);

      auto const predecessors = pDiGraph->getPredecessors(pVertex);

      if (predecessors.empty()) {
          break;
      }

      std::size_t prevDist = 0;
      muchsalsa::graph::Vertex const * prevVertex = nullptr;
      for (auto const &[targetId, pEdge] : predecessors) {
         MS_UNUSED(pEdge);

         auto const *const pPredecessor = pDiGraph->getVertex(targetId);
         auto const cdist = dist.at( mappingVertexToIdx[pPredecessor] );
         if (cdist > prevDist || prevVertex == nullptr) {
            prevVertex = pPredecessor;
            prevDist = cdist;
         }
      }

      pVertex = prevVertex;
  }
  
  std::vector<muchsalsa::graph::Vertex const *> res(std::make_move_iterator(path.begin()), std::make_move_iterator(path.end()) );
  return res;    
 
}



// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
