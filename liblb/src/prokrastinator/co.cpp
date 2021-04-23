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

#include "Prokrastinator.h"

#include <utility>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"

namespace {

std::pair<double, double> computeOverhangs(gsl::not_null<lazybastard::matching::MatchMap const *> const pMatches,
                                           gsl::not_null<lazybastard::graph::Vertex const *> const      pVertex,
                                           gsl::not_null<lazybastard::graph::Edge const *> const        pEdge,
                                           std::string const &                                          illuminaId) {
  auto const pVertexMatch = gsl::make_not_null(pMatches->getVertexMatch(pVertex->getId(), illuminaId));
  auto const pEdgeMatch   = gsl::make_not_null(pMatches->getEdgeMatch(pEdge->getId(), illuminaId));

  auto nanoCorrectionLeft = static_cast<double>(pEdgeMatch->overlap.first - pVertexMatch->illuminaRange.first) //
                            / pVertexMatch->rRatio;
  auto nanoCorrectionRight = static_cast<double>(pVertexMatch->illuminaRange.second - pEdgeMatch->overlap.second) //
                             / pVertexMatch->rRatio;

  lazybastard::util::swap_if(nanoCorrectionLeft, nanoCorrectionRight, !pVertexMatch->direction);

  auto const overhangLeft   = static_cast<double>(pVertexMatch->nanoporeRange.first) + nanoCorrectionLeft;
  auto const nanoporeLength = static_cast<int>(pVertex->getNanoporeLength());
  auto const overhangRight  = static_cast<double>(nanoporeLength - pVertexMatch->nanoporeRange.second) //
                             + nanoCorrectionRight;

  return std::make_pair(overhangLeft, overhangRight);
}

} // unnamed namespace

std::optional<lazybastard::graph::EdgeOrder>
lazybastard::computeOverlap(gsl::not_null<lazybastard::matching::MatchMap const *> const pMatches,
                            std::deque<std::string> &ids, gsl::not_null<lazybastard::graph::Edge const *> const pEdge,
                            bool direction, std::size_t score, bool isPrimary) {
  auto const &firstId = ids.front();
  auto const &lastId  = ids.back();

  auto const vertices = pEdge->getVertices();

  auto const overhangsFirstIdFirstVertex  = computeOverhangs(pMatches, vertices.first, pEdge, firstId);
  auto const overhangsLastIdFirstVertex   = computeOverhangs(pMatches, vertices.first, pEdge, lastId);
  auto const overhangsFirstIdSecondVertex = computeOverhangs(pMatches, vertices.second, pEdge, firstId);
  auto const overhangsLastIdSecondVertex  = computeOverhangs(pMatches, vertices.second, pEdge, lastId);

  auto const leftOverhangFirstVertex  = overhangsFirstIdFirstVertex.first;
  auto const rightOverhangFirstVertex = overhangsLastIdFirstVertex.second;

  auto leftOverhangSecondVertex  = overhangsFirstIdSecondVertex.first;
  auto rightOverhangSecondVertex = overhangsLastIdSecondVertex.second;

  if (!direction) {
    leftOverhangSecondVertex  = overhangsFirstIdSecondVertex.second;
    rightOverhangSecondVertex = overhangsLastIdSecondVertex.first;
  }

  std::optional<graph::EdgeOrder> result;
  if (leftOverhangFirstVertex <= leftOverhangSecondVertex && rightOverhangFirstVertex <= rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.first, vertices.second, leftOverhangSecondVertex - leftOverhangFirstVertex,
                                    rightOverhangSecondVertex - rightOverhangFirstVertex, true, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  } else if (leftOverhangFirstVertex >= leftOverhangSecondVertex &&
             rightOverhangFirstVertex >= rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.second, vertices.first, leftOverhangFirstVertex - leftOverhangSecondVertex,
                                    rightOverhangFirstVertex - rightOverhangSecondVertex, true, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  } else if (leftOverhangFirstVertex > leftOverhangSecondVertex &&
             rightOverhangFirstVertex < rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.first, vertices.second, leftOverhangFirstVertex - leftOverhangSecondVertex,
                                    rightOverhangSecondVertex - rightOverhangFirstVertex, false, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  } else if (leftOverhangFirstVertex < leftOverhangSecondVertex &&
             rightOverhangFirstVertex > rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.second, vertices.first, leftOverhangSecondVertex - leftOverhangFirstVertex,
                                    rightOverhangFirstVertex - rightOverhangSecondVertex, false, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  }

  return result;
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------