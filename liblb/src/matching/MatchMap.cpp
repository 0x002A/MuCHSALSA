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

void MatchMap::addVertexMatch(const std::string &nanoporeId, const std::string &illuminaId,
                              std::shared_ptr<VertexMatch> &&spMatch) {
  if (!spMatch) {
    throw std::runtime_error("Unexpected nullptr.");
  }

  std::scoped_lock<std::mutex> lck(m_mutexVertexMatches);

  // Map representing the scaffolds
  auto nanoporeIDs =
      m_scaffolds
          .insert(std::make_pair(illuminaId,
                                 decltype(m_scaffolds)::mapped_type(lazybastard::MatchingUtil::scaffoldLineIdxCmp)))
          .first;
  nanoporeIDs->second.insert(std::make_pair(m_pGraph->getVertex(nanoporeId), spMatch));

  // Actual map containing the matches
  auto illuminaIDs =
      m_vertexMatches.insert(std::make_pair(nanoporeId, um_t<std::string, std::shared_ptr<VertexMatch>>())).first;
  illuminaIDs->second.insert(std::make_pair(illuminaId, std::move(spMatch)));
}

[[nodiscard]] VertexMatch const *MatchMap::getVertexMatch(std::string const &vertexId,
                                                          std::string const &illuminaId) const {
  auto const vertexIter = m_vertexMatches.find(vertexId);
  if (vertexIter != std::end(m_vertexMatches)) {
    auto const illuminaIter = vertexIter->second.find(illuminaId);
    if (illuminaIter != std::end(vertexIter->second)) {
      return illuminaIter->second.get();
    }
  }

  return nullptr;
}

[[nodiscard]] std::unordered_map<std::string, std::shared_ptr<VertexMatch>> const *
MatchMap::getVertexMatches(std::string const &vertexId) const {
  auto iter = m_vertexMatches.find(vertexId);
  if (iter != std::end(m_vertexMatches)) {
    return &(iter->second);
  }

  return nullptr;
}

void MatchMap::addEdgeMatch(std::string &&edgeId, std::string const &illuminaId, std::shared_ptr<EdgeMatch> &&spMatch) {
  if (!spMatch) {
    throw std::runtime_error("Unexpected nullptr.");
  }

  std::scoped_lock<std::mutex> lck(m_mutexEdgeMatches);

  // Actual map containing the matches
  auto illuminaIDs =
      m_edgeMatches.insert(std::make_pair(std::move(edgeId), um_t<std::string, std::shared_ptr<EdgeMatch>>())).first;
  illuminaIDs->second.insert(std::make_pair(illuminaId, std::move(spMatch)));
}

[[nodiscard]] EdgeMatch const *MatchMap::getEdgeMatch(std::string const &edgeId, std::string const &illuminaId) const {
  auto const edgeIter = m_edgeMatches.find(edgeId);
  if (edgeIter != std::end(m_edgeMatches)) {
    auto const illuminaIter = edgeIter->second.find(illuminaId);
    if (illuminaIter != std::end(edgeIter->second)) {
      return illuminaIter->second.get();
    }
  }

  return nullptr;
}

[[nodiscard]] std::unordered_map<std::string, std::shared_ptr<EdgeMatch>> const *
MatchMap::getEdgeMatches(std::string const &edgeId) const {
  auto const iter = m_edgeMatches.find(edgeId);
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

  // Clear line numbers
  for (auto *const pVertex : m_pGraph->getVertices()) {
    pVertex->clearMetaData();
  }
}

void MatchMap::processScaffold(gsl::not_null<threading::Job const *> const pJob) {
  auto scaffold = std::any_cast<decltype(m_scaffolds)::mapped_type>(pJob->getParam(2));

  for (auto outerIter = std::begin(scaffold); outerIter != std::end(scaffold); ++outerIter) {
    auto const *const outerMatch = outerIter->second.get();
    for (auto innerIter = std::begin(scaffold); innerIter != outerIter; ++innerIter) {
      auto const *const innerMatch = innerIter->second.get();
      auto const overlap = std::make_pair(std::max(outerMatch->illuminaRange.first, innerMatch->illuminaRange.first),
                                          std::min(outerMatch->illuminaRange.second, innerMatch->illuminaRange.second));

      if (overlap.first <= overlap.second && overlap.second - overlap.first > static_cast<int>(TH_OVERLAP)) {
        auto const direction = outerMatch->direction == innerMatch->direction;
        auto const isPrimary = outerMatch->isPrimary == innerMatch->isPrimary;
        auto const outerLength =
            static_cast<double>(outerMatch->illuminaRange.second - outerMatch->illuminaRange.first + 1);
        auto const innerLength =
            static_cast<double>(innerMatch->illuminaRange.second - innerMatch->illuminaRange.first + 1);
        auto const commonLength = static_cast<double>(overlap.second - overlap.first + 1);
        auto const outerScore   = static_cast<double>(outerMatch->score) * commonLength / outerLength;
        auto const innerScore   = static_cast<double>(innerMatch->score) * commonLength / innerLength;
        auto const sumScore     = outerScore + innerScore;

        auto vertexIDs = std::make_pair(innerIter->first, outerIter->first);
        m_pGraph->addEdge(vertexIDs);

        addEdgeMatch(lazybastard::graph::Edge::getEdgeId(std::move(vertexIDs)),
                     std::any_cast<std::string>(pJob->getParam(1)),
                     lazybastard::util::make_shared_aggregate<EdgeMatch>(overlap, direction, sumScore, isPrimary,
                                                                         std::make_pair(0, 0)));
      }
    }
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

} // namespace lazybastard::matching

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
