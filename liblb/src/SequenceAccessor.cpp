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
#include <cstdlib>
#include <cstring>
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

  std::vector<char> buffer(seqPos.second + sizeof(char), '\0');

  size_t ret = fread(buffer.data(), sizeof(decltype(buffer)::value_type), seqPos.second, pFile);
  if (ret != seqPos.second) {
    throw std::runtime_error("Failed to read sequence.");
  }

  std::string sequence(buffer.data());
  sequence.erase(std::remove_if(std::begin(sequence), std::end(sequence), [](auto const c) { return std::isspace(c); }),
                 std::end(sequence));

  return sequence;
}

bool isFastQ(std::string_view filename) {
  auto fileExtension = filename.substr(filename.find_last_of('.') + 1);

  std::vector<char> buffer(fileExtension.size() + 1, '\0');
  std::transform(fileExtension.begin(), fileExtension.end(), buffer.begin(), static_cast<int (*)(int)>(&std::tolower));

  return std::strcmp(fileExtension.data(), "fastq") == 0;
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

SequenceAccessor::SequenceAccessor(gsl::not_null<threading::ThreadPool *> pThreadPool, std::string_view fpNanopore,
                                   std::string_view fpIllumina)
    : m_pThreadPool(pThreadPool),
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

[[maybe_unused]] std::string SequenceAccessor::getNanoporeSequence(std::string const &nanoporeId) {
  std::scoped_lock<std::mutex> lck(m_mutexNanoporeSequenceFile);

  return getSequenceFromFile(m_pNanoporeSequenceFile.get(), m_idxNanopore[nanoporeId]);
}

[[maybe_unused]] std::string SequenceAccessor::getIlluminaSequence(std::string const &illuminaId) {
  std::scoped_lock<std::mutex> lck(m_mutexIlluminaSequenceFile);

  return getSequenceFromFile(m_pIlluminaSequenceFile.get(), m_idxIllumina[illuminaId]);
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
        m_idxNanopore.emplace(sequenceId, std::make_pair(offsetStart, lengthCurrentSequence));

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
        m_idxIllumina.emplace(sequenceId, std::make_pair(offsetStart, lengthCurrentSequence));

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