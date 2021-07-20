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

#ifndef INCLUDED_LAZYBASTARD_PROKRASTINATOR
#define INCLUDED_LAZYBASTARD_PROKRASTINATOR

#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "Lb.fwd.h"
#include "graph/Edge.h"

namespace lazybastard {

// =====================================================================================================================
//                                                        KERNEL
// =====================================================================================================================

std::optional<graph::EdgeOrder> computeOverlap(gsl::not_null<matching::MatchMap const *> matches,
                                               std::vector<unsigned int> &ids, gsl::not_null<graph::Edge const *> pEdge,
                                               bool direction, std::size_t score, bool isPrimary);

std::vector<std::tuple<std::vector<unsigned int>, std::size_t, bool>>
getMaxPairwisePaths(gsl::not_null<matching::MatchMap const *> pMatches, gsl::not_null<graph::Edge const *> pEdge,
                    std::vector<unsigned int> const &illuminaIds, Toggle direction, std::size_t wiggleRoom);

bool sanityCheck(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<graph::Vertex const *> pSubnode,
                 gsl::not_null<graph::Vertex const *> pNode, gsl::not_null<graph::Vertex const *> pTarget,
                 gsl::not_null<graph::EdgeOrder const *> pOrder, std::size_t wiggleRoom);

graph::Graph getMaxSpanTree(gsl::not_null<graph::Graph const *> pGraph);

std::vector<std::vector<lazybastard::graph::Vertex *>>
getConnectedComponents(gsl::not_null<graph::Graph const *> pGraph);

graph::DiGraph getDirectionGraph(gsl::not_null<graph::Graph const *>  pGraph,
                                 gsl::not_null<matching::MatchMap *>  pMatchMap,
                                 gsl::not_null<graph::Graph const *>  pConnectedComponent,
                                 gsl::not_null<graph::Vertex const *> pStartNode);

std::vector<std::vector<lazybastard::graph::Vertex const *>> linearizeGraph(gsl::not_null<graph::DiGraph *> pDiGraph);

void assemblePath(
    gsl::not_null<matching::MatchMap const *> pMatchMap,
    gsl::not_null<std::unordered_map<graph::Vertex const *, std::vector<matching::ContainElement>> const *>
                                      pContainElements,
    gsl::not_null<SequenceAccessor *> pSequenceAccessor, gsl::not_null<matching::Id2OverlapMap *> pId2OverlapMap,
    gsl::not_null<std::vector<lazybastard::graph::Vertex const *> const *> pPath,
    gsl::not_null<graph::DiGraph const *> pDiGraph, int asmIdx, OutputWriter &writer);

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_PROKRASTINATOR

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
