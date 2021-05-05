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

#include "SequenceAccessor.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

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

constexpr char FASTA_IDENTIFIER_DESCLINE  = '>';
constexpr char FASTQ_IDENTIFIER_DESCLINE  = '@';
constexpr char FASTQ_IDENTIFIER_SPLITLINE = '+';

namespace lazybastard {

// =====================================================================================================================
//                                                        HELPER
// =====================================================================================================================

namespace {

std::string getSequenceFromFile(std::FILE *pFile, std::pair<std::size_t, std::size_t> const &seqPos) {
  std::fseek(pFile, static_cast<int64_t>(seqPos.first), SEEK_SET);

  std::vector<char> buffer;
  buffer.reserve(seqPos.second);

  size_t ret = fread(buffer.data(), sizeof(decltype(buffer)::value_type), seqPos.second, pFile);
  if (ret != seqPos.second) {
    throw std::runtime_error("Failed to read sequence.");
  }

  std::string sequence(buffer.data(), seqPos.second);
  sequence.erase(std::remove_if(std::begin(sequence), std::end(sequence), [](auto const c) { return std::isspace(c); }),
                 std::end(sequence));

  return sequence;
}

bool isFastQ(std::string const &filename) {
  auto fileExtension = filename.substr(filename.find_last_of('.') + 1);
  std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(),
                 static_cast<int (*)(int)>(&std::tolower));

  return fileExtension == "fastq";
}

void cleanSequenceId(std::string &sequenceId) {
  auto const iterSpace =
      std::find_if(std::begin(sequenceId), std::end(sequenceId), [](auto const c) { return std::isspace(c); });

  sequenceId.erase(iterSpace, std::end(sequenceId));
}

} // unnamed namespace

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ----------------------
// class SequenceAccessor
// ----------------------

SequenceAccessor::SequenceAccessor(gsl::not_null<threading::ThreadPool *> pThreadPool,
                                   gsl::not_null<graph::Graph *> pGraph, gsl::not_null<matching::MatchMap *> pMatchMap,
                                   std::string const &fpNanopore, std::string const &fpIllumina)
    : m_pThreadPool(pThreadPool), m_pGraph(pGraph), m_pMatchMap(pMatchMap),
      m_pNanoporeSequenceFile(decltype(m_pNanoporeSequenceFile)(
          fopen(fpNanopore.data(), "rbe"), [](gsl::owner<std::FILE *> pFile) { std::fclose(pFile); })),
      m_pIlluminaSequenceFile(decltype(m_pIlluminaSequenceFile)(
          fopen(fpIllumina.data(), "rbe"), [](gsl::owner<std::FILE *> pFile) { std::fclose(pFile); })),
      m_nanoporeFileIsFastQ(isFastQ(fpNanopore)) {
  if (!m_pNanoporeSequenceFile || !m_pIlluminaSequenceFile) {
    throw std::runtime_error("Can't open sequence file(s).");
  }
}

void SequenceAccessor::buildIndex() {
  threading::WaitGroup wg;

  auto jobFnNanoporeIdx = [this](threading::Job const *const pJob) { _buildNanoporeIdx(pJob); };
  auto jobFnIlluminaIdx = [this](threading::Job const *const pJob) { _buildIlluminaIdx(pJob); };

  wg.add(1);
  m_pThreadPool->addJob(threading::Job(jobFnNanoporeIdx, &wg));

  wg.add(1);
  m_pThreadPool->addJob(threading::Job(jobFnIlluminaIdx, &wg));

  wg.wait();
}

[[maybe_unused]] std::string SequenceAccessor::getNanoporeSequence(gsl::not_null<graph::Vertex const *> const pVertex) {
  std::scoped_lock<std::mutex> lck(m_mutexNanoporeSequenceFile);

  return getSequenceFromFile(m_pNanoporeSequenceFile.get(), pVertex->getSeqPos());
}

[[maybe_unused]] std::string
SequenceAccessor::getIlluminaSequence(gsl::not_null<matching::VertexMatch const *> const pVertexMatch) {
  std::scoped_lock<std::mutex> lck(m_mutexIlluminaSequenceFile);

  return getSequenceFromFile(m_pIlluminaSequenceFile.get(), pVertexMatch->seqPos);
}

[[maybe_unused]] std::string
SequenceAccessor::getIlluminaSequence(gsl::not_null<matching::EdgeMatch const *> const pEdgeMatch) {
  std::scoped_lock<std::mutex> lck(m_mutexIlluminaSequenceFile);

  return getSequenceFromFile(m_pIlluminaSequenceFile.get(), pEdgeMatch->seqPos);
}

void SequenceAccessor::_buildNanoporeIdx(gsl::not_null<threading::Job const *> pJob) {
  auto const identifierDescline  = m_nanoporeFileIsFastQ ? FASTQ_IDENTIFIER_DESCLINE : FASTA_IDENTIFIER_DESCLINE;
  auto const identifierSplitline = m_nanoporeFileIsFastQ ? FASTQ_IDENTIFIER_SPLITLINE : FASTA_IDENTIFIER_DESCLINE;

  char *      pLine       = nullptr;
  std::size_t bufferSize  = 0;
  auto        ret         = getline(&pLine, &bufferSize, m_pNanoporeSequenceFile.get());
  auto        offsetStart = std::ftell(m_pNanoporeSequenceFile.get());
  while (ret != -1) {
    if (*pLine == identifierDescline) {
      break;
    }

    ret         = getline(&pLine, &bufferSize, m_pNanoporeSequenceFile.get());
    offsetStart = std::ftell(m_pNanoporeSequenceFile.get());
  }

  while (*pLine == identifierDescline) {
    std::span const spanLine(pLine, bufferSize);
    auto            sequenceId = std::string(spanLine.subspan(1).data());
    cleanSequenceId(sequenceId);

    auto lengthCurrentSequence = 0L;

    while (true) {
      ret            = getline(&pLine, &bufferSize, m_pNanoporeSequenceFile.get());
      auto offsetEnd = std::ftell(m_pNanoporeSequenceFile.get());

      if (ret == -1 || *pLine == identifierSplitline) {
        auto *const pVertex = m_pGraph->getVertex(sequenceId);
        if (pVertex) {
          pVertex->setSeqPos(std::make_pair(offsetStart, lengthCurrentSequence));
        }

        offsetStart = offsetEnd;
        break;
      }

      lengthCurrentSequence += ret;
    }
  }

  if (pLine) {
    free(pLine); // NOLINT
  }

  std::any_cast<threading::WaitGroup *>(pJob->getParam(0))->done();
}

void SequenceAccessor::_buildIlluminaIdx(gsl::not_null<threading::Job const *> pJob) {
  std::unordered_map<std::string, std::vector<matching::VertexMatch *>> mappingIlluminaId2VertexMatches;
  for (auto const &[nanoporeId, vertexMatches] : m_pMatchMap->getVertexMatches()) {
    LB_UNUSED(nanoporeId);

    std::for_each(std::begin(vertexMatches), std::end(vertexMatches), [&](auto const &vertexMatch) {
      auto iterVertexMatches =
          mappingIlluminaId2VertexMatches
              .insert({vertexMatch.first, decltype(mappingIlluminaId2VertexMatches)::mapped_type()})
              .first;
      iterVertexMatches->second.push_back(vertexMatch.second.get());
    });
  }

  std::unordered_map<std::string, std::vector<matching::EdgeMatch *>> mappingIlluminaId2EdgeMatches;
  for (auto const &[nanoporeId, edgeMatches] : m_pMatchMap->getEdgeMatches()) {
    LB_UNUSED(nanoporeId);

    std::for_each(std::begin(edgeMatches), std::end(edgeMatches), [&](auto const &edgeMatch) {
      auto iterEdgeMatches = mappingIlluminaId2EdgeMatches
                                 .insert({edgeMatch.first, decltype(mappingIlluminaId2EdgeMatches)::mapped_type()})
                                 .first;
      iterEdgeMatches->second.push_back(edgeMatch.second.get());
    });
  }

  char *      pLine       = nullptr;
  std::size_t bufferSize  = 0;
  auto        ret         = getline(&pLine, &bufferSize, m_pIlluminaSequenceFile.get());
  auto        offsetStart = std::ftell(m_pIlluminaSequenceFile.get());
  while (ret != -1) {
    if (*pLine == FASTA_IDENTIFIER_DESCLINE) {
      break;
    }

    ret         = getline(&pLine, &bufferSize, m_pIlluminaSequenceFile.get());
    offsetStart = std::ftell(m_pIlluminaSequenceFile.get());
  }

  while (*pLine == FASTA_IDENTIFIER_DESCLINE) {
    std::span const spanLine(pLine, bufferSize);
    auto            sequenceId = std::string(spanLine.subspan(1).data());
    cleanSequenceId(sequenceId);

    auto lengthCurrentSequence = 0L;

    while (true) {
      ret            = getline(&pLine, &bufferSize, m_pIlluminaSequenceFile.get());
      auto offsetEnd = std::ftell(m_pIlluminaSequenceFile.get());

      if (ret == -1 || *pLine == FASTA_IDENTIFIER_DESCLINE) {
        if (mappingIlluminaId2VertexMatches.contains(sequenceId)) {
          auto const &vertexMatches = mappingIlluminaId2VertexMatches[sequenceId];
          std::for_each(std::begin(vertexMatches), std::end(vertexMatches), [&](auto *const pVertexMatch) {
            pVertexMatch->seqPos = std::make_pair(offsetStart, lengthCurrentSequence);
          });
        }

        if (mappingIlluminaId2EdgeMatches.contains(sequenceId)) {
          auto const &edgeMatches = mappingIlluminaId2EdgeMatches[sequenceId];
          std::for_each(std::begin(edgeMatches), std::end(edgeMatches), [&](auto *const pEdgeMatch) {
            pEdgeMatch->seqPos = std::make_pair(offsetStart, lengthCurrentSequence);
          });
        }

        offsetStart = offsetEnd;
        break;
      }

      lengthCurrentSequence += ret;
    }
  }

  if (pLine) {
    free(pLine); // NOLINT
  }

  std::any_cast<threading::WaitGroup *>(pJob->getParam(0))->done();
}

} // namespace lazybastard

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------