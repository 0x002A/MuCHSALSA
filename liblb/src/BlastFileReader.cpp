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

#include "BlastFileReader.h"

#include <algorithm>
#include <any>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "BlastFileAccessor.h"
#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"
#include "threading/Job.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr auto        MINIMUM_MATCHES = 400;
constexpr auto        TH_LENGTH       = 500;
constexpr std::size_t TH_MATCHES      = 500;

constexpr std::size_t POS_IID = 0;
constexpr std::size_t POS_NID = 5;
constexpr std::size_t POS_IRS = 2;
constexpr std::size_t POS_IRE = 3;
constexpr std::size_t POS_NOM = 9;
constexpr std::size_t POS_NLE = 6;
constexpr std::size_t POS_NRS = 7;
constexpr std::size_t POS_NRE = 8;
constexpr std::size_t POS_DIR = 4;

namespace lazybastard {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ---------------------
// class BlastFileReader
// ---------------------

void BlastFileReader::read() {
  threading::WaitGroup wg;
  auto                 jobFn = [this](threading::Job const *const pJob) { parseLine(pJob); };

  for (std::size_t lineIdx = 0; lineIdx < m_pBlastFileAccessor->getLineCount() - 1; ++lineIdx) {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, lineIdx);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void BlastFileReader::parseLine(gsl::not_null<threading::Job const *> const pJob) {
  std::vector<std::string> tokens;

  auto const         lineIdx = std::any_cast<std::size_t>(pJob->getParam(1));
  std::istringstream iss(m_pBlastFileAccessor->getLine(m_pBlastFileAccessor->getLineOffsets().at(lineIdx)),
                         std::ios_base::in);
  std::string        token;
  while (std::getline(iss, token, '\t')) {
    tokens.push_back(token);
  }

  if (tokens.size() < std::max({POS_IID, POS_NID, POS_IRS, POS_IRE, POS_NOM, POS_NLE, POS_NRS, POS_NRE, POS_DIR})) {
    throw std::runtime_error("Invalid BLAST file.");
  }

  auto const illuminaRange = std::make_pair(std::stoi(tokens[POS_IRS]), std::stoi(tokens[POS_IRE]) - 1);
  auto const matches       = static_cast<std::size_t>(std::stoi(tokens[POS_NOM]));

  auto const nanoporeLength = std::stoi(tokens[POS_NLE]);

  auto addNode = matches >= MINIMUM_MATCHES;
  addNode &= illuminaRange.second - illuminaRange.first + 1 >= MINIMUM_MATCHES;

  if (addNode) {
    auto spVertex = std::make_shared<graph::Vertex>(tokens[POS_NID], nanoporeLength, lineIdx);
    m_pGraph->addVertex(std::move(spVertex));

    auto const &nanoporeId = tokens[POS_NID];
    auto const &illuminaId = tokens[POS_IID];

    auto const nanoporeRange = std::make_pair(std::stoi(tokens[POS_NRS]), std::stoi(tokens[POS_NRE]) - 1);
    auto const direction     = tokens[POS_DIR] == "+";
    auto const rRatio        = static_cast<double>(illuminaRange.second - illuminaRange.first + 1) /
                        static_cast<double>(nanoporeRange.second - nanoporeRange.first + 1);

    auto isPrimary = illuminaRange.second - illuminaRange.first + 1 >= TH_LENGTH;
    isPrimary &= matches >= TH_MATCHES;

    auto spVertexMatch = lazybastard::util::make_shared_aggregate<lazybastard::matching::VertexMatch>(
        nanoporeRange, illuminaRange, rRatio, direction, matches, isPrimary, lineIdx, std::make_pair(0, 0));
    m_pMatchMap->addVertexMatch(nanoporeId, illuminaId, spVertexMatch);
  }

  std::any_cast<threading::WaitGroup *const>(pJob->getParam(0))->done();
}

} // namespace lazybastard

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
