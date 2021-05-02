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

#ifndef INCLUDED_LAZYBASTARD_BLASTFILEREADER
#define INCLUDED_LAZYBASTARD_BLASTFILEREADER

#pragma once

#include <gsl/pointers>
#include <iosfwd>

#include "Lb.fwd.h"

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ---------------------
// class BlastFileReader
// ---------------------

/**
 * Class representing a reader for reading files formatted according to the so called BLAST format.
 *
 * The reader creates Vertex instances and adds them to the supplied Graph.
 */
class BlastFileReader {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param pThreadPool a pointer to the ThreadPool to be used for parallelization
   * @param inputStream the input stream of the file
   * @param pGraph a pointer to the Graph receiving the Vertex instances
   * @param pMatchMap a pointer to the MatchMap
   */
  BlastFileReader(gsl::not_null<threading::ThreadPool *> pThreadPool, std::ifstream &inputStream,
                  gsl::not_null<graph::Graph *> pGraph, gsl::not_null<matching::MatchMap *> pMatchMap);

  /**
   * Destructor.
   */
  ~BlastFileReader() = default;

  /**
   * Copying is disallowed.
   */
  BlastFileReader(BlastFileReader const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  BlastFileReader &operator=(BlastFileReader const &) = delete;

  /**
   * Moving is disallowed.
   */
  BlastFileReader(BlastFileReader &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  BlastFileReader &operator=(BlastFileReader &&) = delete;

  /**
   * Reads the file.
   */
  void read();

  /**
   * Parses a line of the file.
   *
   * @param pJob a pointer to the Job containing the parameters
   */
  void parseLine(gsl::not_null<threading::Job const *> pJob);

private:
  threading::ThreadPool *const m_pThreadPool; /*!< Pointer to the ThreadPool used for parallelization */
  std::ifstream &              m_inputStream; /*!< Input stream of the file */
  graph::Graph *const          m_pGraph;      /*!< Pointer to the Graph receiving the Vertex instances */
  matching::MatchMap *const    m_pMatchMap;   /*!< Pointer to the MatchMap */
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ---------------------
// class BlastFileReader
// ---------------------

// PUBLIC CLASS METHODS

inline BlastFileReader::BlastFileReader(gsl::not_null<threading::ThreadPool *> const pThreadPool,
                                        std::ifstream &inputStream, gsl::not_null<graph::Graph *> const pGraph,
                                        gsl::not_null<matching::MatchMap *> const pMatchMap)
    : m_pThreadPool(pThreadPool), m_inputStream(inputStream), m_pGraph(pGraph), m_pMatchMap(pMatchMap) {}

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_BLASTFILEREADER

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------