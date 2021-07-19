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

#include "matching/MatchMap.h"

#include <any>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "Util.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "threading/Job.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr std::size_t TH_OVERLAP = 100;

namespace lazybastard::matching {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// --------------
// class MatchMap
// --------------

void MatchMap::addVertexMatch(unsigned int nanoporeId, unsigned int illuminaId,
                              std::shared_ptr<VertexMatch> const &spMatch) {
  if (!spMatch) {
    throw std::runtime_error("Unexpected nullptr.");
  }

  std::scoped_lock<std::mutex> lck(m_mutexVertexMatches);

  // Actual map containing the matches
  auto &illuminaIds = m_vertexMatches[nanoporeId];

  auto insertMatch = false;
  if (illuminaIds.contains(illuminaId)) {
    auto const lineNumber = illuminaIds[illuminaId]->lineNumber;

    if (lineNumber > spMatch->lineNumber) {
      insertMatch = true;
    }
  } else {
    insertMatch = true;
  }

  if (insertMatch) {
    illuminaIds[illuminaId] = spMatch;
    auto &scaffold          = m_scaffolds[illuminaId];

    scaffold[m_pGraph->getVertex(nanoporeId)] = spMatch;
  }
}

[[nodiscard]] VertexMatch const *MatchMap::getVertexMatch(unsigned int vertexId, unsigned int illuminaId) const {
  auto const vertexIter = m_vertexMatches.find(vertexId);
  if (vertexIter != std::end(m_vertexMatches)) {
    auto const illuminaIter = vertexIter->second.find(illuminaId);
    if (illuminaIter != std::end(vertexIter->second)) {
      return illuminaIter->second.get();
    }
  }

  return nullptr;
}

[[nodiscard]] std::unordered_map<unsigned int, std::shared_ptr<VertexMatch>> const *
MatchMap::getVertexMatches(unsigned int vertexId) const {
  auto iter = m_vertexMatches.find(vertexId);
  if (iter != std::end(m_vertexMatches)) {
    return &(iter->second);
  }

  return nullptr;
}

void MatchMap::addEdgeMatch(lazybastard::graph::Edge const *const pEdge, unsigned int illuminaId,
                            std::shared_ptr<EdgeMatch> const &spMatch) {
  if (!spMatch) {
    throw std::runtime_error("Unexpected nullptr.");
  }

  std::scoped_lock<std::mutex> lck(m_mutexEdgeMatches);

  // Actual map containing the matches
  auto &illuminaIds = m_edgeMatches[pEdge];

  auto insertMatch = false;
  if (illuminaIds.contains(illuminaId)) {
    auto const lineNumber = illuminaIds[illuminaId]->lineNumber;

    if (lineNumber > spMatch->lineNumber) {
      insertMatch = true;
    }
  } else {
    insertMatch = true;
  }

  if (insertMatch) {
    illuminaIds[illuminaId] = spMatch;
  }
}

[[nodiscard]] EdgeMatch const *MatchMap::getEdgeMatch(lazybastard::graph::Edge const *const pEdge,
                                                      unsigned int                          illuminaId) const {
  auto const edgeIter = m_edgeMatches.find(pEdge);
  if (edgeIter != std::end(m_edgeMatches)) {
    auto const illuminaIter = edgeIter->second.find(illuminaId);
    if (illuminaIter != std::end(edgeIter->second)) {
      return illuminaIter->second.get();
    }
  }

  return nullptr;
}

[[nodiscard]] std::unordered_map<unsigned int, std::shared_ptr<EdgeMatch>> const *
MatchMap::getEdgeMatches(lazybastard::graph::Edge const *const pEdge) const {
  auto const iter = m_edgeMatches.find(pEdge);
  if (iter != std::end(m_edgeMatches)) {
    return &iter->second;
  }

  return nullptr;
}

void MatchMap::calculateEdges() {
  threading::WaitGroup wg;
  auto const           jobFn = [this](threading::Job const *const pJob) { processScaffold(pJob); };

  for (auto const &[illuminaId, scaffold] : m_scaffolds) {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, illuminaId, scaffold);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void MatchMap::processScaffold(gsl::not_null<threading::Job const *> const pJob) {
  auto scaffold = std::any_cast<decltype(m_scaffolds)::mapped_type>(pJob->getParam(2));

  std::vector<std::pair<std::size_t, lazybastard::graph::Vertex *>> idx;
  for (auto iter = std::begin(scaffold); iter != std::end(scaffold); ++iter) {
    idx.emplace_back(iter->second->lineNumber, iter->first);
  }

  std::sort(std::begin(idx), std::end(idx), [](auto const &lhs, auto const &rhs) { return lhs.first < rhs.first; });

  for (auto outerIter = std::next(std::begin(idx)); outerIter != std::end(idx); ++outerIter) {
    auto const *const outerMatch = scaffold[outerIter->second].get();
    for (auto innerIter = std::begin(idx); innerIter != outerIter; ++innerIter) {
      auto const *const innerMatch = scaffold[innerIter->second].get();
      auto const overlap = std::make_pair(std::max(outerMatch->illuminaRange.first, innerMatch->illuminaRange.first),
                                          std::min(outerMatch->illuminaRange.second, innerMatch->illuminaRange.second));

      if (overlap.first <= overlap.second && overlap.second - overlap.first > static_cast<int>(TH_OVERLAP)) {
        auto const direction = outerMatch->direction == innerMatch->direction;
        auto const isPrimary = outerMatch->isPrimary && innerMatch->isPrimary;
        auto const outerLength =
            static_cast<double>(outerMatch->illuminaRange.second - outerMatch->illuminaRange.first + 1);
        auto const innerLength =
            static_cast<double>(innerMatch->illuminaRange.second - innerMatch->illuminaRange.first + 1);
        auto const commonLength = static_cast<double>(overlap.second - overlap.first + 1);
        auto const outerScore   = static_cast<double>(outerMatch->score) * commonLength / outerLength;
        auto const innerScore   = static_cast<double>(innerMatch->score) * commonLength / innerLength;
        auto const sumScore     = outerScore + innerScore;

        auto pVertices = [=]() {
          auto const lineIdxInnerVertex = innerIter->second->getMetaDatum<std::size_t>(0);
          auto const lineIdxOuterVertex = outerIter->second->getMetaDatum<std::size_t>(0);

          if (lineIdxOuterVertex < lineIdxInnerVertex) {
            return std::make_pair(outerIter->second, innerIter->second);
          }

          return std::make_pair(innerIter->second, outerIter->second);
        }();
        m_pGraph->addEdge(pVertices);

        addEdgeMatch(m_pGraph->getEdge(std::move(pVertices)), std::any_cast<unsigned int>(pJob->getParam(1)),
                     lazybastard::util::make_shared_aggregate<EdgeMatch>(overlap, direction, sumScore, isPrimary,
                                                                         outerMatch->lineNumber));
      }
    }
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

} // namespace lazybastard::matching

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
