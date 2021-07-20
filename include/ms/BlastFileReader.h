// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of MuCHSALSA <https://github.com/0x002A/MuCHSALSA>.
//
// MuCHSALSA is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// MuCHSALSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with MuCHSALSA.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#ifndef INCLUDED_MUCHSALSA_BLASTFILEREADER
#define INCLUDED_MUCHSALSA_BLASTFILEREADER

#pragma once

#include <gsl/pointers>

#include "Lb.fwd.h"

namespace muchsalsa {

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
   * @param pBlastFileAccessor a pointer to the BlastFileAccessor to be used for accessing the BLAST file
   * @param pGraph a pointer to the Graph receiving the Vertex instances
   * @param pMatchMap a pointer to the MatchMap
   * @param pRegistryNanopore a pointer to the Registry used to register nanopore ids
   * @param pRegistryIllumina a pointer to the Registry used to register nanopore ids
   */
  BlastFileReader(gsl::not_null<threading::ThreadPool *> pThreadPool,
                  gsl::not_null<BlastFileAccessor *> pBlastFileAccessor, gsl::not_null<graph::Graph *> pGraph,
                  gsl::not_null<matching::MatchMap *> pMatchMap, gsl::not_null<Registry *> pRegistryNanopore,
                  gsl::not_null<Registry *> pRegistryIllumina);

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
  BlastFileAccessor
      *const m_pBlastFileAccessor;    /*!< Pointer to the BlastFileAccessor to be used for accessing the BLAST file */
  graph::Graph *const       m_pGraph; /*!< Pointer to the Graph receiving the Vertex instances */
  matching::MatchMap *const m_pMatchMap;         /*!< Pointer to the MatchMap */
  Registry *const           m_pRegistryNanopore; /*!< Pointer to the Registry used to register nanopore ids */
  Registry *const           m_pRegistryIllumina; /*!< Pointer to the Registry used to register illumina ids */
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ---------------------
// class BlastFileReader
// ---------------------

// PUBLIC CLASS METHODS

inline BlastFileReader::BlastFileReader(gsl::not_null<threading::ThreadPool *> const pThreadPool,
                                        gsl::not_null<BlastFileAccessor *>           pBlastFileAccessor,
                                        gsl::not_null<graph::Graph *> const          pGraph,
                                        gsl::not_null<matching::MatchMap *> const    pMatchMap,
                                        gsl::not_null<Registry *> const              pRegistryNanopore,
                                        gsl::not_null<Registry *> const              pRegistryIllumina)
    : m_pThreadPool(pThreadPool), m_pBlastFileAccessor(pBlastFileAccessor), m_pGraph(pGraph), m_pMatchMap(pMatchMap),
      m_pRegistryNanopore(pRegistryNanopore), m_pRegistryIllumina(pRegistryIllumina) {}

} // namespace muchsalsa

#endif // INCLUDED_MUCHSALSA_BLASTFILEREADER

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
