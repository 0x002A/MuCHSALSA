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
  SequenceAccessor(gsl::not_null<threading::ThreadPool *> pThreadPool, gsl::not_null<graph::Graph *> pGraph,
                   gsl::not_null<matching::MatchMap *> pMatchMap, std::string const &fpNanopore,
                   std::string const &fpIllumina);

  /**
   * Builds the sequence index.
   */
  void buildIndex();

  /**
   * Returns the nanopore sequence of the supplied Vertex instance.
   * This function is **thread-safe**.
   *
   * @param pVertex a pointer pointing to the Vertex instance to return the nanopore sequence for
   * @return A std::string representing the nanopore sequence
   */
  [[maybe_unused]] std::string getNanoporeSequence(gsl::not_null<graph::Vertex const *> pVertex);

  /**
   * Returns the illumina sequence of the supplied VertexMatch instance.
   * This function is **thread-safe**.
   *
   * @param pVertexMatch a pointer pointing to the VertexMatch instance to return the illumina sequence for
   * @return A std::string representing the illumina sequence
   */
  [[maybe_unused]] std::string getIlluminaSequence(gsl::not_null<matching::VertexMatch const *> pVertexMatch);

  /**
   * Returns the illumina sequence of the supplied EdgeMatch instance.
   * This function is **thread-safe**.
   *
   * @param pEdgeMatch a pointer pointing to the EdgeMatch instance to return the illumina sequence for
   * @return A std::string representing the illumina sequence
   */
  [[maybe_unused]] std::string getIlluminaSequence(gsl::not_null<matching::EdgeMatch const *> pEdgeMatch);

private:
  threading::ThreadPool *const m_pThreadPool; /*!< Pointer to the ThreadPool used for parallelization */
  graph::Graph *const          m_pGraph;      /*!< Pointer to the Graph */
  matching::MatchMap *const    m_pMatchMap;   /*!< Pointer to the MatchMap */
  std::unique_ptr<std::FILE, void (*)(std::FILE *)>
      m_pNanoporeSequenceFile; /*!< Pointer pointing to the FILE object required for accessing the nanopore sequences
                                */
  std::unique_ptr<std::FILE, void (*)(std::FILE *)>
       m_pIlluminaSequenceFile; /*!< Pointer pointing to the FILE object required for accessing illumina sequences */
  bool m_nanoporeFileIsFastQ;   /*!< bool stating whether the file containing the illumina sequences is a FastQ file */
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