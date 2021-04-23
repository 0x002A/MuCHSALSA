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
#include <iterator>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "Util.h"
#include "graph/Graph.h"

namespace {

template <class INPUT_IT_NON_NULL_IN_DEGREES, class INPUT_IT_NULL_IN_DEGREES>
void initializeInDegreeMaps(INPUT_IT_NON_NULL_IN_DEGREES nonNullInDegrees, INPUT_IT_NULL_IN_DEGREES nullInDegrees,
                            gsl::not_null<lazybastard::graph::DiGraph *> const pDiGraphCycle) {
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

void sortReduction(gsl::not_null<lazybastard::graph::DiGraph *> const pDiGraphCycle) {
  std::map<lazybastard::graph::Vertex const *, std::size_t> verticesWithNonNullInDegree;
  std::deque<lazybastard::graph::Vertex const *>            verticesWithNullInDegree;

  initializeInDegreeMaps(std::inserter(verticesWithNonNullInDegree, std::begin(verticesWithNonNullInDegree)),
                         std::back_inserter(verticesWithNullInDegree), pDiGraphCycle);

  std::unordered_set<lazybastard::graph::Vertex const *> neighbors;
  if (!verticesWithNonNullInDegree.empty()) {
    neighbors.insert(std::begin(verticesWithNonNullInDegree)->first);
  }

  while (true) {
    while (!verticesWithNullInDegree.empty()) {
      auto const *const pVertex = verticesWithNullInDegree.front();
      verticesWithNullInDegree.pop_front();

      auto const successors = pDiGraphCycle->getSuccessors(pVertex);
      for (auto const &[idSuccessor, pEdge] : successors) {
        LB_UNUSED(pEdge);

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

    auto hasNoNeighbors = neighbors.empty();
    if (hasNoNeighbors) {
      std::for_each(std::begin(verticesWithNonNullInDegree), std::end(verticesWithNonNullInDegree),
                    [&](auto const &p) { neighbors.insert(p.first); });
    }

    lazybastard::graph::Vertex const *pMinVertex = nullptr;
    auto                              minScore   = 0L;
    for (auto const *const pNeighbor : neighbors) {
      auto       score        = 0L;
      auto const predecessors = pDiGraphCycle->getPredecessors(pNeighbor);

      score = std::count_if(std::begin(predecessors), std::end(predecessors), [&](auto const &p) {
        auto const *const pPredecessor = pDiGraphCycle->getVertex(p.first);
        return verticesWithNonNullInDegree.contains(pPredecessor);
      });

      if (!pMinVertex || score < minScore) {
        pMinVertex = pNeighbor;
        minScore   = score;
      }
    }

    auto const predecessors = pDiGraphCycle->getPredecessors(pMinVertex);
    std::for_each(std::begin(predecessors), std::end(predecessors), [&](auto const &p) {
      auto const *const pPredecessor = pDiGraphCycle->getVertex(p.first);
      if (verticesWithNonNullInDegree.contains(pPredecessor)) {
        pDiGraphCycle->deleteEdge(p.second);
        p.second->setShadow(true);
      }
    });

    verticesWithNonNullInDegree.erase(pMinVertex);
    verticesWithNullInDegree.push_back(pMinVertex);
    neighbors.erase(pMinVertex);

    if (hasNoNeighbors) {
      neighbors.clear();
    }
  }
}

void sortTopologically(gsl::not_null<std::vector<lazybastard::graph::Vertex const *> *> const pResult,
                       gsl::not_null<lazybastard::graph::DiGraph *> const                     pDiGraphCycle) {
  pResult->clear();

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> verticesWithNonNullInDegree;
  std::deque<lazybastard::graph::Vertex const *>                      verticesWithNullInDegree;

  initializeInDegreeMaps(std::inserter(verticesWithNonNullInDegree, std::begin(verticesWithNonNullInDegree)),
                         std::back_inserter(verticesWithNullInDegree), pDiGraphCycle);

  while (!verticesWithNullInDegree.empty()) {
    auto const *const pVertex = verticesWithNullInDegree.front();
    verticesWithNullInDegree.pop_front();

    auto const successors = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      LB_UNUSED(pEdge);

      auto const *pSuccessor = pDiGraphCycle->getVertex(targetId);

      verticesWithNonNullInDegree[pSuccessor] -= 1;

      if (verticesWithNonNullInDegree[pSuccessor] == 0) {
        verticesWithNullInDegree.push_back(pSuccessor);
        verticesWithNonNullInDegree.erase(pSuccessor);
      }
    }

    pResult->push_back(pVertex);
  }
}

void findClusterWeights(
    gsl::not_null<std::unordered_map<lazybastard::graph::Edge const *, std::size_t> *> const pResult,
    gsl::not_null<lazybastard::graph::DiGraph *> const                                       pDiGraphCycle) {
  pResult->clear();

  std::vector<lazybastard::graph::Vertex const *> sortedVertices;
  sortedVertices.reserve(pDiGraphCycle->getOrder());

  sortTopologically(&sortedVertices, pDiGraphCycle);
  sortedVertices.shrink_to_fit();

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> mappingVertexToIdx;
  std::unordered_map<std::size_t, lazybastard::graph::Vertex const *> mappingIdxToVertex;

  for (auto iter = std::begin(sortedVertices); iter != std::end(sortedVertices); ++iter) {
    auto const idx = std::distance(std::begin(sortedVertices), iter);

    mappingVertexToIdx.insert({*iter, idx});
    mappingIdxToVertex.insert({idx, *iter});
  }

  auto const edges = pDiGraphCycle->getEdges();
  std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) { pResult->insert({pEdge, 0}); });

  std::unordered_map<lazybastard::graph::Vertex const *, std::set<std::size_t>> mappingPredecessors;
  std::unordered_map<lazybastard::graph::Vertex const *, std::set<std::size_t>> mappingSuccessors;

  std::for_each(std::begin(sortedVertices), std::end(sortedVertices), [&](auto const *const pVertex) {
    auto insertIterMappingPredecessors =
        mappingPredecessors.insert({pVertex, decltype(mappingPredecessors)::mapped_type()});
    auto insertIterMappingSuccessors =
        mappingSuccessors.insert({pVertex, decltype(mappingPredecessors)::mapped_type()});

    auto const successors = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      LB_UNUSED(targetId);

      insertIterMappingSuccessors.first->second.insert(mappingVertexToIdx[pEdge->getVertices().second]);
    }

    auto const predecessors = pDiGraphCycle->getPredecessors(pVertex);
    for (auto const &[targetId, pEdge] : predecessors) {
      LB_UNUSED(targetId);

      insertIterMappingPredecessors.first->second.insert(mappingVertexToIdx[pEdge->getVertices().first]);
    }
  });

  for (auto const *const pVertex : sortedVertices) {
    std::vector<
        std::tuple<decltype(mappingPredecessors)::mapped_type, std::vector<decltype(mappingVertexToIdx)::mapped_type>>>
        candidates({{mappingSuccessors[pVertex], std::vector<std::size_t>({mappingVertexToIdx[pVertex]})}});

    for (auto const idx : mappingSuccessors[pVertex]) {
      auto const *const pActiveVertex = mappingIdxToVertex[idx];

      for (auto const idxIn : mappingPredecessors[pActiveVertex]) {
        for (auto iter = std::begin(candidates); iter != std::end(candidates); ++iter) {
          auto const &visited = std::get<1>(*iter);
          auto const &open    = std::get<0>(*iter);

          if (visited.back() == idxIn && open.contains(idx)) {
            decltype(mappingSuccessors)::mapped_type intersection;
            std::set_intersection(std::begin(open), std::end(open), std::begin(mappingSuccessors[pActiveVertex]),
                                  std::end(mappingSuccessors[pActiveVertex]),
                                  std::inserter(intersection, std::begin(intersection)));
            auto newVisited = visited;
            newVisited.push_back(idx);

            auto const currentIterIdx = std::distance(std::begin(candidates), iter);
            candidates.emplace_back(std::move(intersection), std::move(newVisited));
            iter = std::begin(candidates) + currentIterIdx;
          }
        }
      }

      decltype(candidates) filtered;
      filtered.reserve(candidates.size());
      for (auto outerIter = std::begin(candidates); outerIter != std::end(candidates); ++outerIter) {
        auto const i           = std::distance(std::begin(candidates), outerIter);
        auto       isDominated = false;

        for (auto innerIter = std::begin(candidates); innerIter != std::end(candidates); ++innerIter) {
          auto const j = std::distance(std::begin(candidates), innerIter);

          if (i != j && std::includes(std::begin(std::get<0>(*outerIter)), std::end(std::get<0>(*innerIter)),
                                      std::begin(std::get<1>(*outerIter)), std::end(std::get<1>(*innerIter)))) {
            isDominated = true;
            break;
          }
        }

        if (!isDominated) {
          filtered.push_back(*outerIter);
        }
      }
      filtered.shrink_to_fit();
      candidates = std::move(filtered);
    }

    std::vector<std::vector<std::size_t>> maxVisited;
    std::size_t                           maxLength = 0;
    for (auto const &[open, visited] : candidates) {
      LB_UNUSED(open);

      if (visited.size() > maxLength) {
        maxVisited = decltype(maxVisited)({visited});
        maxLength  = visited.size();
      } else if (visited.size() == maxLength) {
        maxVisited.push_back(visited);
      }
    }

    for (auto const &mv : maxVisited) {
      auto c = mv.size() - 1;

      auto const limit = std::max(mv.size() - 1, 0UL);
      for (std::size_t i = 0; i < limit; ++i) {
        auto const        pV1   = lazybastard::util::make_not_null_and_const(mappingIdxToVertex[i]);
        auto const        pV2   = lazybastard::util::make_not_null_and_const(mappingIdxToVertex[i + 1]);
        auto const *const pEdge = pDiGraphCycle->getEdge(std::make_pair(pV1, pV2));

        (*pResult)[pEdge] += c;
        c -= 1;
      }
    }
  }
}

std::vector<lazybastard::graph::Vertex const *> findConservationPath(
    gsl::not_null<lazybastard::graph::DiGraph *> const                                             pDiGraphCycle,
    gsl::not_null<std::unordered_map<lazybastard::graph::Edge const *, std::size_t> const *> const pClusterWeights) {
  std::vector<lazybastard::graph::Vertex const *> sortedVertices;
  sortedVertices.reserve(pDiGraphCycle->getOrder());

  sortTopologically(&sortedVertices, pDiGraphCycle);
  sortedVertices.shrink_to_fit();

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> mappingVertexToIdx;
  std::unordered_map<std::size_t, lazybastard::graph::Vertex const *> mappingIdxToVertex;

  for (auto iter = std::begin(sortedVertices); iter != std::end(sortedVertices); ++iter) {
    auto const idx = std::distance(std::begin(sortedVertices), iter);

    mappingVertexToIdx.insert({*iter, idx});
    mappingIdxToVertex.insert({idx, *iter});
  }

  std::vector<std::vector<lazybastard::graph::Vertex const *>>                          finalizedPaths;
  std::vector<std::tuple<std::size_t, std::vector<lazybastard::graph::Vertex const *>>> openPaths;

  auto const &outDegrees = pDiGraphCycle->getOutDegrees();
  for (auto const *const pVertex : sortedVertices) {
    if (outDegrees.at(pVertex) == 0) {
      for (auto const &[val, path] : openPaths) {
        LB_UNUSED(val);

        if (path.back() == pVertex) {
          finalizedPaths.emplace_back(decltype(finalizedPaths)::value_type(std::begin(path), std::end(path)));
          auto const iterMaxPath =
              std::max_element(std::begin(finalizedPaths), std::end(finalizedPaths),
                               [](auto const &p1, auto const &p2) { return p1.size() < p2.size(); });
          finalizedPaths = decltype(finalizedPaths)({std::move(*iterMaxPath)});
        }
      }
      continue;
    }

    std::vector<lazybastard::graph::Edge const *> maxOuts;
    std::size_t                                   maxOut     = 0;
    auto const                                    successors = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      LB_UNUSED(targetId);

      auto const currentClusterWeight = pClusterWeights->at(pEdge);
      if (currentClusterWeight > maxOut) {
        maxOut  = currentClusterWeight;
        maxOuts = decltype(maxOuts)({pEdge});
      } else if (currentClusterWeight == maxOut) {
        maxOuts.push_back(pEdge);
      }
    }

    std::vector<std::size_t> maxIns;
    std::size_t              maxIn = 0;
    for (auto iter = std::begin(openPaths); iter != std::end(openPaths); ++iter) {
      auto const idx = static_cast<size_t const>(std::distance(std::begin(openPaths), iter));

      if (std::get<1>(*iter).back() == pVertex) {

        auto const currentVal = std::get<0>(*iter);
        if (currentVal > maxIn) {
          maxIn  = currentVal;
          maxIns = decltype(maxIns)({idx});
        } else if (currentVal == maxIn) {
          maxIns.push_back(idx);
        }
      }
    }

    for (auto const *const pEdge : maxOuts) {
      if (!maxIns.empty()) {
        auto const iterMaxElem =
            std::max_element(std::begin(maxIns), std::end(maxIns), [&](auto const idx1, auto const idx2) {
              return std::get<1>(openPaths[idx1]).size() < std::get<1>(openPaths[idx2]).size();
            });

        auto tmp = std::get<1>(openPaths[*iterMaxElem]);
        tmp.push_back(pEdge->getVertices().second);

        openPaths.emplace_back(std::make_tuple(maxOut, std::move(tmp)));
      } else {
        auto const vertices = pEdge->getVertices();
        openPaths.emplace_back(std::make_tuple(
            maxOut, std::tuple_element_t<1, decltype(openPaths)::value_type>({vertices.first, vertices.second})));
      }
    }
    std::erase_if(openPaths, [&](auto const &t) {
      return mappingVertexToIdx[std::get<1>(t).back()] <= mappingVertexToIdx[pVertex];
    });
  }

  return finalizedPaths.front();
}

std::vector<std::vector<lazybastard::graph::Vertex const *>>
extractPaths(gsl::not_null<lazybastard::graph::DiGraph *> const pDiGraph) {
  auto diGraphCycle = *pDiGraph;

  auto const edges = diGraphCycle.getEdges();
  for (auto const *const pEdge : edges) {
    if (pEdge->isShadow()) {
      diGraphCycle.deleteEdge(pEdge);
    }
  }

  sortReduction(&diGraphCycle);

  std::unordered_map<lazybastard::graph::Edge const *, std::size_t> clusterWeights;
  findClusterWeights(&clusterWeights, &diGraphCycle);

  std::vector<std::vector<lazybastard::graph::Vertex const *>> paths;
  std::unordered_set<lazybastard::graph::Vertex const *>       visited;

  while (diGraphCycle.getSize() > 0) {
    auto const longestPath = findConservationPath(&diGraphCycle, &clusterWeights);

    if (longestPath.size() < 10) {
      auto isInVisit = false;

      auto const predecessors = pDiGraph->getPredecessors(longestPath.front());

      for (auto const &[targetId, pEdge] : predecessors) {
        LB_UNUSED(pEdge);

        auto const *const pPredecessor = pDiGraph->getVertex(targetId);
        lazybastard::util::exchange_if(isInVisit, true, visited.contains(pPredecessor));
      }

      auto       isOutVisit = false;
      auto const successors = pDiGraph->getSuccessors(longestPath.back());
      for (auto const &[targetId, pEdge] : successors) {
        LB_UNUSED(pEdge);

        auto const *const pSuccessor = pDiGraph->getVertex(targetId);
        lazybastard::util::exchange_if(isOutVisit, true, visited.contains(pSuccessor));
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
    paths.push_back(decltype(paths)::value_type({pVertex}));
  }

  return paths;
}

} // unnamed namespace

std::vector<std::vector<lazybastard::graph::Vertex const *>>
lazybastard::linearizeGraph(gsl::not_null<lazybastard::graph::DiGraph *> const pDiGraph) {
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
          std::begin(paths[idx2]), std::find(std::begin(paths[idx2]), std::end(paths[idx2]), vertices.first))));
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

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------