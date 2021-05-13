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

#ifndef INCLUDED_LAZYBASTARD_SEQUENCEACCESSOR
#define INCLUDED_LAZYBASTARD_SEQUENCEACCESSOR

#pragma once

#include <cstdio>
#include <gsl/pointers>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>

#include "Lb.fwd.h"

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ----------------------
// class SequenceAccessor
// ----------------------

/**
 * Class providing access to the nanopore and illumina sequences.
 *
 * The methods of this class are **partially** thread-safe.
 */
class SequenceAccessor {
public:
  /**
   * Constructor creating a new instance.
   *
   * @param pThreadPool a pointer pointing to the ThreadPool to be used for parallelization
   * @param fpNanopore a std::string_view representing the path to the file containing the nanopore sequences
   * @param fpIllumina a std::string_view representing the path to the file containing the illumina sequences
   */
  SequenceAccessor(gsl::not_null<threading::ThreadPool *> pThreadPool, std::string_view fpNanopore,
                   std::string_view fpIllumina);

  /**
   * Builds the sequence index.
   */
  void buildIndex();

  /**
   * Returns the nanopore sequence of the supplied nanopore id.
   * This function is **thread-safe**.
   *
   * @param nanoporeId a const reference to the std::string representing the nanopore id to return the sequence for
   * @return A std::string representing the nanopore sequence
   */
  std::string getNanoporeSequence(std::string const &nanoporeId);

  /**
   * Returns the illumina sequence of the supplied illumina id.
   * This function is **thread-safe**.
   *
   * @param illuminaId a const reference to the std::string representing the illumina id to return the sequence for
   * @return A std::string representing the illumina sequence
   */
  std::string getIlluminaSequence(std::string const &illuminaId);

private:
  threading::ThreadPool *const m_pThreadPool; /*!< Pointer to the ThreadPool used for parallelization */
  std::unique_ptr<std::FILE, void (*)(std::FILE *)>
      m_pNanoporeSequenceFile; /*!< Pointer pointing to the FILE handle required for accessing the nanopore sequences
                                */
  std::unique_ptr<std::FILE, void (*)(std::FILE *)>
       m_pIlluminaSequenceFile; /*!< Pointer pointing to the FILE handle required for accessing illumina sequences */
  bool m_nanoporeFileIsFastQ;   /*!< bool stating whether the file containing the illumina sequences is a FastQ file */
  std::unordered_map<std::string, std::pair<long, long>>
      m_idxNanopore; /*!< std::unordered_map containing the mapping of nanopore ids to the position of the sequences */
  std::unordered_map<std::string, std::pair<long, long>>
      m_idxIllumina; /*!< std::unordered_map containing the mapping of illumina ids to the position of the sequences */
  std::mutex m_mutexNanoporeSequenceFile; /*!< std::mutex for securing the parallel use of the file containing the
                                             nanopore sequences */
  std::mutex m_mutexIlluminaSequenceFile; /*!< std::mutex for securing the parallel use of the file containing the
                                             illumina sequences */

  /**
   * Assigns the offset and length data to the Vertex instances of the attached Graph.
   *
   * @param pJob a pointer to the Job containing the parameters
   */
  void _buildNanoporeIdx(gsl::not_null<threading::Job const *> pJob);

  /**
   * Assigns the offset and length data to the VertexMatch and EdgeMatch instances of the attached MatchMap.
   *
   * @param pJob a pointer to the Job containing the parameters
   */
  void _buildIlluminaIdx(gsl::not_null<threading::Job const *> pJob);
};

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_SEQUENCEACCESSOR

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------