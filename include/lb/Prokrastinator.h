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
#include <tuple>
#include <unordered_map>
#include <vector>

#include "Lb.fwd.h"
#include "graph/Edge.h"

namespace lazybastard {

// =====================================================================================================================
//                                                        KERNEL
// =====================================================================================================================

std::optional<graph::EdgeOrder> computeOverlap(matching::MatchMap const &matches, std::vector<unsigned int> &ids,
                                               graph::Edge const &edge, bool direction, std::size_t score,
                                               bool isPrimary);

std::vector<std::tuple<std::vector<unsigned int>, std::size_t, bool>>
getMaxPairwisePaths(matching::MatchMap const &matches, graph::Edge const &edge,
                    std::vector<unsigned int> const &illuminaIds, Toggle direction, std::size_t wiggleRoom);

bool sanityCheck(graph::Graph const &graph, graph::Vertex const &subnode, graph::Vertex const &node,
                 graph::Vertex const &target, graph::EdgeOrder const &order, std::size_t wiggleRoom);

graph::Graph getMaxSpanTree(graph::Graph const &graph);

std::vector<std::vector<lazybastard::graph::Vertex *>> getConnectedComponents(graph::Graph const &graph);

graph::DiGraph getDirectionGraph(gsl::not_null<matching::MatchMap *> pMatchMap, graph::Graph const &graph,
                                 graph::Graph const &connectedComponent, graph::Vertex const &startNode);

std::vector<std::vector<lazybastard::graph::Vertex const *>> linearizeGraph(gsl::not_null<graph::DiGraph *> pDiGraph);

void assemblePath(
    gsl::not_null<matching::Id2OverlapMap *> pId2OverlapMap, matching::MatchMap const &matchMap,
    std::unordered_map<graph::Vertex const *, std::vector<matching::ContainElement>> const &containElements,
    SequenceAccessor &sequenceAccessor, std::vector<lazybastard::graph::Vertex const *> const &path,
    graph::DiGraph const &diGraph, int asmIdx, OutputWriter &writer);

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_PROKRASTINATOR

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
