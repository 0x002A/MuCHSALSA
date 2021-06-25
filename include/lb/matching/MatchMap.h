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

#ifndef INCLUDED_LAZYBASTARD_MATCHMAP
#define INCLUDED_LAZYBASTARD_MATCHMAP

#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>

#include "Lb.fwd.h"
#include "graph/Vertex.h"
#include "types/Toggle.h"

namespace lazybastard {

// =====================================================================================================================
//                                                   UTILITY FUNCTIONS
// =====================================================================================================================

// -------------------
// struct MatchingUtil
// -------------------

/**
 * Struct providing a namespace for utility functions regarding the matching.
 */
struct MatchingUtil {
  /**
   * Compares the line index of two Vertex instances.
   *
   * @param pLhs a pointer pointing to the left hand side Vertex
   * @param pRhs a pointer pointing to the right hand side Vertex
   * @return A bool indicating whether the left hand side Vertex has a smaller line index.
   */
  static bool scaffoldLineIdxCmp(graph::Vertex const *pLhs, graph::Vertex const *pRhs);
};

namespace matching {

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
  std::pair<int, int> const           nanoporeRange; /*!< Nanopore range */
  std::pair<int, int> const           illuminaRange; /*!< Illumina range */
  double const                        rRatio;        /*!< Read ratio */
  Toggle const                        direction;     /*!< Read direction */
  std::size_t const                   score;         /*!< Score (number of matches) */
  Toggle const                        isPrimary;     /*!< Is a primary */
  std::size_t                         lineNumber;    /*!< Position in BLAST file */
  std::pair<std::size_t, std::size_t> seqPos;        /*!< Sequence offset and length */
};

// ----------------
// struct EdgeMatch
// ----------------

/**
 * Struct representing a match attached to an Edge.
 */
struct EdgeMatch {
  std::pair<int, int> const           overlap;    /*!< Overlap */
  Toggle const                        direction;  /*!< Edge direction */
  double const                        score;      /*!< Score */
  Toggle const                        isPrimary;  /*!< Is a primary */
  std::size_t                         lineNumber; /*!< Position in BLAST file */
  std::pair<std::size_t, std::size_t> seqPos;     /*!< Sequence offset and length */
};

// ---------------------
// struct ContainElement
// ---------------------

struct ContainElement {
  std::unordered_map<std::string, VertexMatch const *> const matches;
  std::string                                                nano;
  std::size_t const                                          nanoporeLength;
  std::size_t const                                          score;
  lazybastard::Toggle const                                  direction;
  bool const                                                 isPrimary;
};

// --------------
// class MatchMap
// --------------

/**
 * Class representing the mapping of illumina ids to their matching nanopore ids.
 *
 * Instances of this class are designed to be thread-safe.
 */
class MatchMap {
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
   * @param nanoporeId a const reference to the nanopore id
   * @param illuminaId a const reference to the illumina id
   * @param spMatch a const reference to the std::shared_ptr pointing to the VertexMatch to be copied into the MatchMap
   */
  void addVertexMatch(std::string const &nanoporeId, std::string const &illuminaId,
                      std::shared_ptr<VertexMatch> const &spMatch);

  /**
   * Getter returning a pointer pointing to a specific VertexMatch of the Vertex having the supplied id.
   * This functions returns nullptr if the requested VertexMatch wasn't found.
   *
   * @param vertexId a const reference to a std::string representing the id of the Vertex
   * @param illuminaId a const reference to the illumina id associated with the VertexMatch
   * @return A pointer to the VertexMatch (nullptr if not found)
   */
  [[nodiscard]] VertexMatch const *getVertexMatch(std::string const &vertexId, std::string const &illuminaId) const;

  /**
   * Getter returning a pointer pointing to the std::unordered_map containing all VertexMatch instances for a specific
   * Vertex stored within this MatchMap.
   *
   * @param vertexId a const reference to a std::string representing the id of the Vertex
   * @return A pointer to the std::unordered_map containing all VertexMatch instances associated with the requested
   *         Vertex (nullptr if not found)
   */
  [[nodiscard]] std::unordered_map<std::string, std::shared_ptr<VertexMatch>> const *
  getVertexMatches(std::string const &vertexId) const;

  /**
   * Getter returning all VertexMatch instances stored within this MatchMap.
   *
   * @return A const reference to the std::unordered_map containing all VertexMatch instances
   */
  [[nodiscard]] auto const &getVertexMatches() const;

  /**
   * Deletes all VertexMatch instances for a specific Vertex stored within this MatchMap.
   *
   * @param vertexId a const reference to a std::string representing the id of the Vertex
   */
  void deleteVertexMatches(std::string const &vertexId);

  /**
   * Adds an EdgeMatch to the MatchMap.
   *
   * @param pEdge a pointer pointing to the Edge to add the matches for
   * @param illuminaId a const reference to a std::string representing the illumina id
   * @param spMatch a const reference to the std::shared_ptr pointing to the EdgeMatch to be copied into the MatchMap
   */
  void addEdgeMatch(lazybastard::graph::Edge const *pEdge, std::string const &illuminaId,
                    std::shared_ptr<EdgeMatch> const &spMatch);

  /**
   * Getter returning a specific EdgeMatch of the Edge having the supplied id.
   * This functions returns nullptr if the requested EdgeMatch wasn't found.
   *
   * @param pEdge a pointer pointing to the Edge to return the matches for
   * @param illuminaId a const reference to a std::string representing the illumina id
   * @return A pointer pointing to the EdgeMatch (nullptr if not found)
   */
  [[nodiscard]] EdgeMatch const *getEdgeMatch(lazybastard::graph::Edge const *pEdge,
                                              std::string const &             illuminaId) const;

  /**
   * Getter returning a std::unordered_map containing all EdgeMatch instances for a specific Edge stored within this
   * MatchMap.
   *
   * @param pEdge a pointer pointing to the Edge to return the matches for
   * @return A pointer pointing to the std::unordered_map containing all EdgeMatch instances associated with the
   *         requested Edge
   */
  [[nodiscard]] std::unordered_map<std::string, std::shared_ptr<EdgeMatch>> const *
  getEdgeMatches(lazybastard::graph::Edge const *pEdge) const;

  /**
   * Getter returning all EdgeMatch instances stored within this MatchMap.
   *
   * @return A const reference to the std::unordered_map containing all EdgeMatch instances
   */
  [[maybe_unused]] [[nodiscard]] auto const &getEdgeMatches() const;

  /**
   * Deletes all EdgeMatch instances for a specific Edge stored within this MatchMap.
   *
   * @param pEdge a pointer pointing to the Edge to delete the matches for
   */
  void deleteEdgeMatches(lazybastard::graph::Edge const *pEdge);

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

private:
  template <class KEY, class VALUE> using um_t = std::unordered_map<KEY, VALUE>;

  um_t<std::string, um_t<std::string, std::shared_ptr<VertexMatch>>>
      m_vertexMatches; /*!< std::unordered_map containing the VertexMatch instances */
  um_t<lazybastard::graph::Edge const *, um_t<std::string, std::shared_ptr<EdgeMatch>>>
      m_edgeMatches; /*!< std::unordered_map containing the EdgeMatch instances */
  um_t<std::string, um_t<graph::Vertex *, std::shared_ptr<VertexMatch>>>
             m_scaffolds;          /*!< std::unordered_map containing the scaffolds */
  std::mutex m_mutexVertexMatches; /*!< std::mutex for securing the parallel use of the std::unordered_map containing
                               the VertexMatches */
  std::mutex m_mutexEdgeMatches; /*!< std::mutex for securing the parallel use of the std::unordered_map containing the
                                    EdgeMatches */
  threading::ThreadPool *const m_pThreadPool; /*!< Pointer pointing to the ThreadPool used for parallelization */
  graph::Graph *const          m_pGraph;      /*!< Pointer pointing to the Graph receiving the Vertex instances */
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
    : m_pThreadPool(pThreadPool), m_pGraph(pGraph) {}

[[nodiscard]] inline auto const &MatchMap::getVertexMatches() const { return m_vertexMatches; }

inline void MatchMap::deleteVertexMatches(std::string const &vertexId) { m_vertexMatches.erase(vertexId); }

[[maybe_unused]] [[nodiscard]] inline auto const &MatchMap::getEdgeMatches() const { return m_edgeMatches; }

inline void MatchMap::deleteEdgeMatches(lazybastard::graph::Edge const *pEdge) { m_edgeMatches.erase(pEdge); }

} // namespace matching

// -------------------
// struct MatchingUtil
// -------------------

// PUBLIC CLASS METHODS

inline bool MatchingUtil::scaffoldLineIdxCmp(graph::Vertex const *pLhs, graph::Vertex const *pRhs) {
  return pLhs->getMetaDatum<std::size_t>(0) < pRhs->getMetaDatum<std::size_t>(0);
}

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_MATCHMAP

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
