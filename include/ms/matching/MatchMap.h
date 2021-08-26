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

#ifndef INCLUDED_MUCHSALSA_MATCHMAP
#define INCLUDED_MUCHSALSA_MATCHMAP

#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <memory>
#include <shared_mutex>
#include <unordered_map>
#include <utility>

#include "Lb.fwd.h"
#include "graph/Graph.h"
#include "types/Toggle.h"

namespace muchsalsa::matching {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ------------------
// struct VertexMatch
// ------------------

/**
 * Struct representing a match attached to a Vertex.
 */
struct VertexMatch {
  std::pair<int, int> const nanoporeRange; /*!< Nanopore range */
  std::pair<int, int> const illuminaRange; /*!< Illumina range */
  double const              rRatio;        /*!< Read ratio */
  Toggle const              direction;     /*!< Read direction */
  std::size_t const         score;         /*!< Score (number of matches) */
  Toggle const              isPrimary;     /*!< Is a primary */
  std::size_t               lineNumber;    /*!< Position in BLAST file */
};

// ----------------
// struct EdgeMatch
// ----------------

/**
 * Struct representing a match attached to an Edge.
 */
struct EdgeMatch {
  std::pair<int, int> const overlap;    /*!< Overlap */
  Toggle const              direction;  /*!< Edge direction */
  double const              score;      /*!< Score */
  Toggle const              isPrimary;  /*!< Is a primary */
  std::size_t               lineNumber; /*!< Position in BLAST file */
};

// ---------------------
// struct ContainElement
// ---------------------

struct ContainElement {
  std::unordered_map<unsigned int, VertexMatch const *> const matches;
  unsigned int                                                nano;
  std::size_t const                                           nanoporeLength;
  std::size_t const                                           score;
  muchsalsa::Toggle const                                     direction;
  bool const                                                  isPrimary;
};

// --------------
// class MatchMap
// --------------

/**
 * Class representing the mapping of illumina ids to their matching nanopore ids.
 *
 * Instances of this class are designed to be thread-safe.
 */
class MatchMap : public muchsalsa::graph::IGraphObserver {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param pThreadPool a pointer pointing to the ThreadPool to be used for parallelization
   * @param pGraph a pointer pointing to the Graph receiving the Vertex instances
   */
  MatchMap(gsl::not_null<threading::ThreadPool *> pThreadPool, gsl::not_null<graph::Graph *> pGraph);

  /**
   * Adds a VertexMatch to the MatchMap.
   *
   * @param nanoporeId an unsigned int representing the nanopore id
   * @param illuminaId an unsigned int representing the illumina id
   * @param spMatch a const reference to the std::shared_ptr pointing to the VertexMatch to be copied into the MatchMap
   */
  void addVertexMatch(unsigned int nanoporeId, unsigned int illuminaId, std::shared_ptr<VertexMatch> const &spMatch);

  /**
   * Getter returning a pointer pointing to a specific VertexMatch of the Vertex having the supplied id.
   * This functions returns nullptr if the requested VertexMatch wasn't found.
   *
   * @param vertexId an unsigned int representing the id of the Vertex
   * @param illuminaId an unsigned int representing the illumina id associated with the VertexMatch
   * @return A pointer to the VertexMatch (nullptr if not found)
   */
  [[nodiscard]] VertexMatch const *getVertexMatch(unsigned int vertexId, unsigned int illuminaId) const;

  /**
   * Getter returning a pointer pointing to the std::unordered_map containing all VertexMatch instances for a specific
   * Vertex stored within this MatchMap.
   *
   * @param vertexId an unsigned int representing the id of the Vertex
   * @return A pointer to the std::unordered_map containing all VertexMatch instances associated with the requested
   *         Vertex (nullptr if not found)
   */
  [[nodiscard]] std::unordered_map<unsigned int, std::shared_ptr<VertexMatch>> const *
  getVertexMatches(unsigned int vertexId) const;

  /**
   * Deletes all VertexMatch instances for a specific Vertex stored within this MatchMap.
   *
   * @param vertexId an unsigned int representing the id of the Vertex
   */
  void deleteVertexMatches(unsigned int vertexId);

  /**
   * Adds an EdgeMatch to the MatchMap.
   *
   * @param pEdge a pointer pointing to the Edge to add the matches for
   * @param illuminaId an unsigned int representing the illumina id
   * @param spMatch a const reference to the std::shared_ptr pointing to the EdgeMatch to be copied into the MatchMap
   */
  void addEdgeMatch(muchsalsa::graph::Edge const *pEdge, unsigned int illuminaId,
                    std::shared_ptr<EdgeMatch> const &spMatch);

  /**
   * Getter returning a specific EdgeMatch of the Edge having the supplied id.
   * This functions returns nullptr if the requested EdgeMatch wasn't found.
   *
   * @param pEdge a pointer pointing to the Edge to return the matches for
   * @param illuminaId an unsigned int representing the illumina id
   * @return A pointer pointing to the EdgeMatch (nullptr if not found)
   */
  [[nodiscard]] EdgeMatch const *getEdgeMatch(muchsalsa::graph::Edge const *pEdge, unsigned int illuminaId) const;

  /**
   * Getter returning a std::unordered_map containing all EdgeMatch instances for a specific Edge stored within this
   * MatchMap.
   *
   * @param pEdge a pointer pointing to the Edge to return the matches for
   * @return A pointer pointing to the std::unordered_map containing all EdgeMatch instances associated with the
   *         requested Edge
   */
  [[nodiscard]] std::unordered_map<unsigned int, std::shared_ptr<EdgeMatch>> const *
  getEdgeMatches(muchsalsa::graph::Edge const *pEdge) const;

  /**
   * Deletes all EdgeMatch instances for a specific Edge stored within this MatchMap.
   *
   * @param pEdge a pointer pointing to the Edge to delete the matches for
   */
  void deleteEdgeMatches(muchsalsa::graph::Edge const *pEdge);

  /**
   * Creates or updates Edge instances according to the scaffolds.
   */
  void calculateEdges();

  /**
   * Processes a single scaffold.
   *
   * @param pJob a pointer pointing to the Job containing the parameters
   */
  void processScaffold(gsl::not_null<threading::Job const *> pJob);

  /**
   * Hook which is getting called every time a Vertex is about to get deleted.
   *
   * @param pVertex a const pointer to the Vertex instance which is about to get deleted
   */
  void onVertexDeleted(muchsalsa::graph::Vertex const *pVertex);

  /**
   * Hook which is getting called every time an Edge is about to get deleted.
   *
   * @param pEdge a const pointer to the Edge instance which is about to get deleted
   */
  void onEdgeDeleted(muchsalsa::graph::Edge const *pEdge);

private:
  template <class KEY, class VALUE> using um_t = std::unordered_map<KEY, VALUE>;

  um_t<unsigned int, um_t<unsigned int, std::shared_ptr<VertexMatch>>>
      m_vertexMatches; /*!< std::unordered_map containing the VertexMatch instances */
  um_t<muchsalsa::graph::Edge const *, um_t<unsigned int, std::shared_ptr<EdgeMatch>>>
      m_edgeMatches; /*!< std::unordered_map containing the EdgeMatch instances */
  um_t<unsigned int, um_t<graph::Vertex *, std::shared_ptr<VertexMatch>>>
                            m_scaffolds;          /*!< std::unordered_map containing the scaffolds */
  mutable std::shared_mutex m_mutexVertexMatches; /*!< std::shared_mutex for securing the parallel use of the
                               std::unordered_map containing the VertexMatches */
  mutable std::shared_mutex m_mutexEdgeMatches; /*!< std::shared_mutex for securing the parallel use of the std::unordered_map
                                    containing the EdgeMatches */
  threading::ThreadPool *const m_pThreadPool;   /*!< Pointer pointing to the ThreadPool used for parallelization */
  graph::Graph *const          m_pGraph;        /*!< Pointer pointing to the Graph receiving the Vertex instances */
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// --------------
// class MatchMap
// --------------

// PUBLIC CLASS METHODS
inline MatchMap::MatchMap(gsl::not_null<threading::ThreadPool *> const pThreadPool,
                          gsl::not_null<graph::Graph *> const          pGraph)
    : m_pThreadPool(pThreadPool), m_pGraph(pGraph) {
  pGraph->attachObserver(this);
}

inline void MatchMap::deleteVertexMatches(unsigned int vertexId) {
  std::lock_guard<std::shared_mutex> lck(m_mutexVertexMatches);
  m_vertexMatches.erase(vertexId);
}

inline void MatchMap::deleteEdgeMatches(muchsalsa::graph::Edge const *pEdge) {
  std::lock_guard<std::shared_mutex> lck(m_mutexEdgeMatches);
  m_edgeMatches.erase(pEdge);
}

} // namespace muchsalsa::matching

#endif // INCLUDED_MUCHSALSA_MATCHMAP

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
