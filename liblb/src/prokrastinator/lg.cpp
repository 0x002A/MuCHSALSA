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
#include <deque>
#include <iterator>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "Util.h"
#include "graph/Graph.h"

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
        p.second->setShadow(true);
        pDiGraphCycle->deleteEdge(p.second, nullptr);
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

void findClusterWeights(
    gsl::not_null<std::unordered_map<lazybastard::graph::Edge const *, std::size_t> *> const pResult,
    gsl::not_null<lazybastard::graph::DiGraph *> const                                       pDiGraphCycle) {
  pResult->clear();

  auto const sortedVertices = pDiGraphCycle->sortTopologically();

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> mappingVertexToIdx;

  std::size_t idx = 0;
  for (auto const *const pVertex : sortedVertices) {
    mappingVertexToIdx.insert({pVertex, idx});

    ++idx;
  }

  auto const edges = pDiGraphCycle->getEdges();
  std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) { pResult->emplace(pEdge, 0); });

  std::unordered_map<lazybastard::graph::Vertex const *, std::set<std::size_t>> mappingPredecessors;
  std::unordered_map<lazybastard::graph::Vertex const *, std::set<std::size_t>> mappingSuccessors;

  std::for_each(std::begin(sortedVertices), std::end(sortedVertices), [&](auto const *const pVertex) {
    auto &     idxSuccessors = mappingSuccessors[pVertex]; // NOLINT
    auto const successors    = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      LB_UNUSED(pEdge);

      idxSuccessors.insert(mappingVertexToIdx.at(pDiGraphCycle->getVertex(targetId)));
    }

    auto &     idxPredecessors = mappingPredecessors[pVertex]; // NOLINT
    auto const predecessors    = pDiGraphCycle->getPredecessors(pVertex);
    for (auto const &[targetId, pEdge] : predecessors) {
      LB_UNUSED(pEdge);

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
      LB_UNUSED(open);

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

std::vector<lazybastard::graph::Vertex const *> findConservationPath(
    gsl::not_null<lazybastard::graph::DiGraph *> const                                             pDiGraphCycle,
    gsl::not_null<std::unordered_map<lazybastard::graph::Edge const *, std::size_t> const *> const pClusterWeights) {
  auto const sortedVertices = pDiGraphCycle->sortTopologically();

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> mappingVertexToIdx;

  for (std::size_t idx = 0; idx < sortedVertices.size(); ++idx) {
    mappingVertexToIdx.insert({sortedVertices[idx], idx});
  }

  std::vector<std::vector<lazybastard::graph::Vertex const *>>                          finalizedPaths;
  std::vector<std::tuple<std::size_t, std::vector<lazybastard::graph::Vertex const *>>> openPaths;

  auto const &outDegrees = pDiGraphCycle->getOutDegrees();
  for (auto const *const pVertex : sortedVertices) {
    if (outDegrees.at(pVertex) == 0) {
      for (auto const &[val, path] : openPaths) {
        LB_UNUSED(val);

        if (path.back() == pVertex) {
          finalizedPaths.emplace_back(std::begin(path), std::end(path));
          auto const iterMaxPath =
              std::max_element(std::begin(finalizedPaths), std::end(finalizedPaths),
                               [](auto const &pLhs, auto const &pRhs) { return pLhs.size() < pRhs.size(); });
          finalizedPaths = {std::move(*iterMaxPath)};
        }
      }
      continue;
    }

    std::vector<std::pair<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *>> maxOuts;
    std::size_t                                                                                    maxOut = 0;
    auto const successors = pDiGraphCycle->getSuccessors(pVertex);
    for (auto const &[targetId, pEdge] : successors) {
      auto edge = [&](auto const *const pEdge, auto const &targetId) {
        auto vertices = pEdge->getVertices();
        if (vertices.second->getId() == targetId) {
          return vertices;
        }

        return std::make_pair(vertices.second, vertices.first);
      }(pEdge, targetId);

      auto const currentClusterWeight = pClusterWeights->at(pEdge);
      if (currentClusterWeight > maxOut) {
        maxOut  = currentClusterWeight;
        maxOuts = {std::move(edge)};
      } else if (currentClusterWeight == maxOut) {
        maxOuts.push_back(std::move(edge));
      }
    }

    std::vector<std::size_t> maxIns;
    std::size_t              maxIn = 0;
    for (std::size_t idx = 0; idx < openPaths.size(); ++idx) {

      if (std::get<1>(openPaths.at(idx)).back() == pVertex) {

        auto const currentVal = std::get<0>(openPaths.at(idx));
        if (currentVal > maxIn) {
          maxIn  = currentVal;
          maxIns = decltype(maxIns)({idx});
        } else if (currentVal == maxIn) {
          maxIns.push_back(idx);
        }
      }
    }

    for (auto const &edge : maxOuts) {
      if (!maxIns.empty()) {
        auto const iterMaxElem =
            std::max_element(std::begin(maxIns), std::end(maxIns), [&](auto const idxLhs, auto const idxRhs) {
              return std::get<1>(openPaths[idxLhs]).size() < std::get<1>(openPaths[idxRhs]).size();
            });

        auto tmp = std::get<1>(openPaths[*iterMaxElem]);
        tmp.push_back(edge.second);

        openPaths.emplace_back(std::make_tuple(maxOut, std::move(tmp)));
      } else {
        openPaths.emplace_back(maxOut,
                               std::initializer_list<lazybastard::graph::Vertex const *>({edge.first, edge.second}));
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
      diGraphCycle.deleteEdge(pEdge, nullptr);
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
      diGraphCycle.deleteVertex(pVertex, nullptr);
    }
  }

  for (auto const *const pVertex : diGraphCycle.getVertices()) {
    paths.push_back({pVertex});
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

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
