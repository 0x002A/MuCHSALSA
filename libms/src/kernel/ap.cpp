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

#include "Kernel.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "OutputWriter.h"
#include "Registry.h"
#include "SequenceAccessor.h"
#include "SequenceUtils.h"
#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/Id2OverlapMap.h"
#include "matching/MatchMap.h"
#include "types/Direction.h"
#include "types/Toggle.h"

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr auto SEQUENCE_LINE_LENGTH = 60;
constexpr auto TH_SEQUENCE_LENGTH   = 200;

// =====================================================================================================================
//                                                        METHODS
// =====================================================================================================================

namespace {

std::string limitLength(std::string_view original) {
  std::string formatted = std::string(original);
  for (auto iter = std::begin(formatted); iter != std::end(formatted);) {
    if (iter != std::begin(formatted)) {
      iter = std::next(formatted.insert(iter, '\n'));
    }

    if (std::distance(iter, std::end(formatted)) > SEQUENCE_LINE_LENGTH) {
      std::advance(iter, SEQUENCE_LINE_LENGTH);
    } else {
      iter = std::end(formatted);
    }
  }

  return formatted;
}

std::string tupleTuple2Id(std::tuple<std::tuple<unsigned int, std::size_t>, std::size_t> const &tupleTuple) {
  std::string id;

  auto const &innerTuple = std::get<0>(tupleTuple);
  id.append(std::to_string(std::get<0>(innerTuple)));
  id.append(",");
  id.append(std::to_string(std::get<1>(innerTuple)));
  id.append(",");
  id.append(std::to_string(std::get<1>(tupleTuple)));

  return id;
}

std::vector<muchsalsa::graph::Vertex *> ramseyR2(muchsalsa::graph::Graph const &         graph,
                                                 std::vector<muchsalsa::graph::Vertex *> vertices) {
  if (vertices.empty()) {
    return std::vector<muchsalsa::graph::Vertex *>();
  }

  auto        neighbors    = decltype(vertices)();
  auto        nonNeighbors = decltype(vertices)();
  auto *const firstVertex  = vertices[0];
  std::for_each(std::next(std::begin(vertices)), std::end(vertices),
                [=, &graph, &neighbors, &nonNeighbors](auto *const pVertex) {
                  if (graph.hasEdge(std::make_pair(firstVertex, pVertex))) {
                    neighbors.push_back(pVertex);
                  } else {
                    nonNeighbors.push_back(pVertex);
                  }
                });

  auto       cliqueNeighbors    = ramseyR2(graph, neighbors);
  auto const cliqueNonNeighbors = ramseyR2(graph, nonNeighbors);

  cliqueNeighbors.push_back(firstVertex);

  return cliqueNeighbors.size() >= cliqueNonNeighbors.size() ? cliqueNeighbors : cliqueNonNeighbors;
}

std::vector<std::vector<muchsalsa::graph::Vertex *>> getAnchorCliques(muchsalsa::graph::Graph const &graph) {
  auto const toVector = [](auto const &source) {
    return std::vector<typename std::remove_reference_t<decltype(source)>::value_type>(std::begin(source),
                                                                                       std::end(source));
  };

  auto                                 vertices      = graph.getVerticesAsUnorderedSet();
  auto                                 currentClique = ramseyR2(graph, toVector(vertices));
  std::vector<decltype(currentClique)> cliques{currentClique};

  while (!vertices.empty()) {
    std::for_each(std::begin(currentClique), std::end(currentClique),
                  [&](auto *const pVertex) { vertices.erase(pVertex); });

    currentClique = ramseyR2(graph, toVector(vertices));
    if (!currentClique.empty()) {
      cliques.push_back(currentClique);
    }
  }

  return cliques;
}

void getClusterAnchors(
    gsl::not_null<std::vector<std::unordered_map<unsigned int, std::size_t>> *> const pClusterModifier,
    gsl::not_null<muchsalsa::matching::Id2OverlapMap *> const                         pId2OverlapMap,
    muchsalsa::matching::MatchMap const &matchMap, unsigned int illuminaIdBase,
    std::vector<unsigned int> const &edgeIdx, std::vector<muchsalsa::graph::Edge const *> const &edges) {

  muchsalsa::graph::Graph g;
  for (auto const idxEdge1 : edgeIdx) {
    g.addVertex(std::make_shared<muchsalsa::graph::Vertex>(idxEdge1, 0));

    for (auto const idxEdge2 : edgeIdx) {
      if (idxEdge1 == idxEdge2) {
        break;
      }

      auto const overlapEdge1 = matchMap.getEdgeMatch(edges.at(idxEdge1), illuminaIdBase)->overlap;
      auto const overlapEdge2 = matchMap.getEdgeMatch(edges.at(idxEdge2), illuminaIdBase)->overlap;

      auto const overlap = std::make_pair(std::max(overlapEdge1.first, overlapEdge2.first),
                                          std::min(overlapEdge1.second, overlapEdge2.second));

      if (overlap.first <= overlap.second) {
        g.addEdge(std::make_pair(g.getVertex(idxEdge2), g.getVertex(idxEdge1)));
      }
    }
  }

  auto const anchorCliques = getAnchorCliques(g);

  for (std::size_t idx = 0; idx < anchorCliques.size(); ++idx) {
    auto const correctedId = std::make_tuple(illuminaIdBase, idx);

    std::optional<std::pair<int, int>> commonOverlap;
    for (auto const *const pVertex : anchorCliques.at(idx)) {
      (*pClusterModifier)[pVertex->getId()][illuminaIdBase] = idx;

      auto const *const pMatch = matchMap.getEdgeMatch(edges.at(pVertex->getId()), illuminaIdBase);
      if (!commonOverlap) {
        commonOverlap = pMatch->overlap;
      } else {
        auto const otherOverlap = pMatch->overlap;

        commonOverlap = std::make_pair(std::max(commonOverlap->first, otherOverlap.first),
                                       std::min(commonOverlap->second, otherOverlap.second));
      }
    }

    (*pId2OverlapMap)[correctedId] = commonOverlap.value();
  }
}

std::pair<double, double> getCorrectedNanoporeRange(muchsalsa::matching::MatchMap const &matchMap,
                                                    unsigned int nanoporeId, unsigned int illuminaId,
                                                    std::pair<int, int> const &overlap) {
  auto const *const pMatch = matchMap.getVertexMatch(nanoporeId, illuminaId);

  auto nanoCorrectionLeft  = (overlap.first - pMatch->illuminaRange.first) / pMatch->rRatio;
  auto nanoCorrectionRight = (pMatch->illuminaRange.second - overlap.second) / pMatch->rRatio;

  muchsalsa::util::swap_if(nanoCorrectionLeft, nanoCorrectionRight, !pMatch->direction);

  return std::make_pair(pMatch->nanoporeRange.first + nanoCorrectionLeft,
                        pMatch->nanoporeRange.second - nanoCorrectionRight);
}

std::tuple<std::string, int, int> updateConsensusBase(std::optional<std::string_view> oldSequence,
                                                      std::pair<int, int> oldBorders, std::string_view newSequence,
                                                      std::pair<int, int> newBorders) {
  if (!oldSequence.has_value()) {
    return std::make_tuple(std::string(newSequence), newBorders.first, newBorders.second);
  }

  std::string updatedSequence;
  if (newBorders.first < oldBorders.first) {
    auto const borderRight = oldBorders.first - newBorders.first;

    updatedSequence.append([=]() { return muchsalsa::strSlice(newSequence, 0, borderRight); }());
    updatedSequence.append(oldSequence.value());
  } else if (newBorders.second > oldBorders.second) {
    updatedSequence.append(oldSequence.value());

    auto const borderLeft = -(newBorders.second - oldBorders.second);
    updatedSequence.append(muchsalsa::strSlice(newSequence, borderLeft, static_cast<int>(newSequence.size())));
  } else {
    updatedSequence = oldSequence.value();
  }

  return std::make_tuple(std::move(updatedSequence), std::min(oldBorders.first, newBorders.first),
                         std::max(oldBorders.second, newBorders.second));
}

std::tuple<std::optional<std::string>, int, int>
visitOrdered(gsl::not_null<std::unordered_map<muchsalsa::graph::Vertex const *, bool> *> const pVisitedVertices,
             gsl::not_null<std::unordered_map<muchsalsa::graph::Vertex const *, std::tuple<int, int>> *> pTap,
             muchsalsa::graph::DiGraph const &                                                           adg,
             std::unordered_map<unsigned int, std::tuple<unsigned int, std::size_t>> const &             regIdx2Id,
             std::unordered_map<muchsalsa::graph::Vertex const *, std::size_t> const &           mappingVertex2OrderIdx,
             std::vector<muchsalsa::graph::Vertex const *> const &                               order,
             std::unordered_map<muchsalsa::graph::Edge const *, int> const &                     distances,
             std::unordered_map<muchsalsa::graph::Edge const *, std::vector<std::string>> const &sequences,
             std::unordered_map<muchsalsa::graph::Vertex const *, std::string> const &           anchorSequences,
             muchsalsa::matching::Id2OverlapMap const &id2OverlapMap, muchsalsa::graph::Vertex const &vertexStart) {
  std::optional<std::string> sequence;

  std::set<std::pair<std::size_t, int>, bool (*)(std::pair<int, int>, std::pair<int, int>)> queueEdges(
      [](std::pair<int, int> lhs, std::pair<int, int> rhs) {
        lhs.second = -lhs.second;
        rhs.second = -rhs.second;

        return lhs < rhs;
      });
  std::set<std::size_t> queueVertices;

  int borderLeft  = 0;
  int borderRight = 0;

  queueVertices.insert(mappingVertex2OrderIdx.at(&vertexStart));
  while (!queueVertices.empty()) {
    auto const        idx     = *std::begin(queueVertices);
    auto const *const pVertex = order[idx];

    queueVertices.erase(std::begin(queueVertices));

    if (!(*pVisitedVertices)[pVertex]) {
      (*pVisitedVertices)[pVertex] = true;

      for (auto const &[targetId, pEdge] : adg.getSuccessors(pVertex)) {

        auto const orderSuccessor = mappingVertex2OrderIdx.at(adg.getVertex(targetId));
        queueEdges.emplace(orderSuccessor, idx);
        queueVertices.insert(orderSuccessor);
      }

      while (!queueEdges.empty() && (*std::begin(queueEdges)).first == idx) {
        auto const getAnchor = [&](bool first) {
          auto const idx = first ? static_cast<std::size_t>((*std::begin(queueEdges)).first)
                                 : static_cast<std::size_t>((*std::begin(queueEdges)).second);
          return order[idx];
        };

        auto const *const pAnchorLeft  = getAnchor(false);
        auto const *const pAnchorRight = getAnchor(true);

        auto const hasAnchorLeft  = pTap->contains(pAnchorLeft);
        auto const hasAnchorRight = pTap->contains(pAnchorRight);

        auto const overlapLeft  = id2OverlapMap.at(regIdx2Id.at(pAnchorLeft->getId()));
        auto const overlapRight = id2OverlapMap.at(regIdx2Id.at(pAnchorRight->getId()));

        auto const *const pEdge  = adg.getEdge(std::make_pair(pAnchorLeft, pAnchorRight));
        auto const        offset = distances.at(pEdge);

        auto const lengthLeft  = overlapLeft.second - overlapLeft.first + 1;
        auto const lengthRight = overlapRight.second - overlapRight.first + 1;

        if (hasAnchorLeft and !hasAnchorRight) {
          auto const posRight   = std::get<1>(pTap->at(pAnchorLeft));
          (*pTap)[pAnchorRight] = std::make_tuple(posRight + offset + 1, posRight + offset + lengthRight);

          if (offset > 0) {
            std::tie(sequence, borderLeft, borderRight) =
                updateConsensusBase(sequence, std::make_pair(borderLeft, borderRight), sequences.at(pEdge).front(),
                                    std::make_pair(posRight + 1, posRight + offset));
          }

          std::tie(sequence, borderLeft, borderRight) = updateConsensusBase(
              sequence, std::make_pair(borderLeft, borderRight), anchorSequences.at(pAnchorRight),
              std::make_pair(std::get<0>(pTap->at(pAnchorRight)), std::get<1>(pTap->at(pAnchorRight))));
        } else if (!hasAnchorLeft and hasAnchorRight) {
          auto const posRight  = std::get<0>(pTap->at(pAnchorRight));
          (*pTap)[pAnchorLeft] = std::make_tuple(posRight - offset - lengthLeft, posRight - offset - 1);

          if (offset > 0) {
            std::tie(sequence, borderLeft, borderRight) =
                updateConsensusBase(sequence, std::make_pair(borderLeft, borderRight), sequences.at(pEdge).front(),
                                    std::make_pair(posRight - offset, posRight));
          }

          std::tie(sequence, borderLeft, borderRight) = updateConsensusBase(
              sequence, std::make_pair(borderLeft, borderRight), anchorSequences.at(pAnchorLeft),
              std::make_pair(std::get<0>(pTap->at(pAnchorLeft)), std::get<1>(pTap->at(pAnchorLeft))));
        } else if (!hasAnchorLeft and !hasAnchorRight) {
          (*pTap)[pAnchorLeft]  = std::make_tuple(0, lengthLeft - 1);
          (*pTap)[pAnchorRight] = std::make_tuple(lengthLeft + offset, lengthLeft + offset + lengthRight - 1);

          if (offset > 0) {
            std::tie(sequence, borderLeft, borderRight) =
                updateConsensusBase(sequence, std::make_pair(borderLeft, borderRight), sequences.at(pEdge).front(),
                                    std::make_pair(lengthLeft, lengthLeft + offset - 1));
          }

          std::tie(sequence, borderLeft, borderRight) = updateConsensusBase(
              sequence, std::make_pair(borderLeft, borderRight), anchorSequences.at(pAnchorLeft),
              std::make_pair(std::get<0>(pTap->at(pAnchorLeft)), std::get<1>(pTap->at(pAnchorLeft))));
          std::tie(sequence, borderLeft, borderRight) = updateConsensusBase(
              sequence, std::make_pair(borderLeft, borderRight), anchorSequences.at(pAnchorRight),
              std::make_pair(std::get<0>(pTap->at(pAnchorRight)), std::get<1>(pTap->at(pAnchorRight))));
        }

        queueEdges.erase(std::begin(queueEdges));
      }
    } else {
      while (!queueEdges.empty() && (*std::begin(queueEdges)).first == idx) {
        queueEdges.erase(std::begin(queueEdges));
      }
    }
  }

  return std::make_tuple(sequence, borderLeft, borderRight);
}


std::string getSequenceLeftOfAnchor(muchsalsa::matching::MatchMap const &matchMap,
                                    muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int nanoporeId,
                                    std::size_t nanoporeLength, unsigned int illuminaId,
                                    std::pair<int, int> illuminaOverlap, muchsalsa::Toggle const direction) {
  auto const *const pMatch = matchMap.getVertexMatch(nanoporeId, illuminaId);

  if (!direction) {
    std::string illuminaSequence;

    if (!pMatch->direction) {
      illuminaSequence =
          muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, false);
    } else {
      illuminaSequence =
          muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second, pMatch->illuminaRange.second, true);
    }

    illuminaSequence.append(muchsalsa::getNanoporeSequence(sequenceAccessor, nanoporeId, pMatch->nanoporeRange.second,
                                                static_cast<int>(nanoporeLength) - 1, true));

    return muchsalsa::getReverseComplement(illuminaSequence);
  }

  auto nanoporeSequence = muchsalsa::getNanoporeSequence(sequenceAccessor, nanoporeId, 0, pMatch->nanoporeRange.first, true);

  if (!pMatch->direction) {
    nanoporeSequence.append(
        muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second, pMatch->illuminaRange.second, false));
  } else {
    nanoporeSequence.append(
        muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, true));
  }

  return nanoporeSequence;
}

std::string getSequenceRightOfAnchor(muchsalsa::matching::MatchMap const &matchMap,
                                     muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int nanoporeId,
                                     std::size_t nanoporeLength, unsigned int illuminaId,
                                     std::pair<int, int> illuminaOverlap, muchsalsa::Toggle const direction) {
  auto const *const pMatch = matchMap.getVertexMatch(nanoporeId, illuminaId);

  if (!direction) {
    auto nanoporeSequence = muchsalsa::getNanoporeSequence(sequenceAccessor, nanoporeId, 0, pMatch->nanoporeRange.first, true);

    if (!pMatch->direction) {
      nanoporeSequence.append(muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second,
                                                  pMatch->illuminaRange.second, false));
    } else {
      nanoporeSequence.append(
          muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, true));
    }

    return muchsalsa::getReverseComplement(nanoporeSequence);
  }

  std::string illuminaSequence;

  if (!pMatch->direction) {
    illuminaSequence =
        muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, false);
  } else {
    illuminaSequence =
        muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second, pMatch->illuminaRange.second, true);
  }

  illuminaSequence.append(muchsalsa::getNanoporeSequence(sequenceAccessor, nanoporeId, pMatch->nanoporeRange.second,
                                              static_cast<int>(nanoporeLength) - 1, true));

  return illuminaSequence;
}

std::string getAnchorSequence(muchsalsa::matching::MatchMap const &matchMap,
                              muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int nanoporeId,
                              unsigned int illuminaId, std::pair<int, int> illuminaOverlap,
                              muchsalsa::Toggle const direction) {
  auto const combinedDirection = matchMap.getVertexMatch(nanoporeId, illuminaId)->direction * direction;
  auto       illuminaSequence  = muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.first,
                                              illuminaOverlap.second, combinedDirection);

  return illuminaSequence;
}

std::tuple<int, std::optional<std::string>>
getSequenceBetweenAnchors(muchsalsa::matching::MatchMap const &matchMap, muchsalsa::SequenceAccessor &sequenceAccessor,
                          unsigned int nanoporeId, unsigned int illuminaIdLeft, unsigned int illuminaIdRight,
                          std::pair<int, int> const &overlapLeft, std::pair<int, int> const &overlapRight,
                          muchsalsa::Toggle const direction) {
  auto const *const matchLeft  = matchMap.getVertexMatch(nanoporeId, illuminaIdLeft);
  auto const *const matchRight = matchMap.getVertexMatch(nanoporeId, illuminaIdRight);

  auto const &illuminaRangeLeft  = matchLeft->illuminaRange;
  auto const &illuminaRangeRight = matchRight->illuminaRange;

  auto const rRatioLeft  = matchLeft->rRatio;
  auto const rRatioRight = matchRight->rRatio;

  auto const &nanoporeRangeLeft  = matchLeft->nanoporeRange;
  auto const &nanoporeRangeRight = matchRight->nanoporeRange;

  auto correctionLeft  = 0;
  auto correctionRight = 0;

  if (!direction) {
    double errorOffset = nanoporeRangeRight.second - nanoporeRangeLeft.first;
    if (errorOffset > 0) {
      auto const correctedNanoporeLeft = getCorrectedNanoporeRange(matchMap, nanoporeId, illuminaIdLeft, overlapLeft);
      auto const correctedNanoporeRight =
          getCorrectedNanoporeRange(matchMap, nanoporeId, illuminaIdRight, overlapRight);

      if (correctedNanoporeLeft.first < correctedNanoporeRight.second) {
        return std::make_tuple(std::floor(correctedNanoporeLeft.first - correctedNanoporeRight.second), std::nullopt);
      }

      double availableLeft  = 0;
      double availableRight = 0;
      if (!matchLeft->direction) {
        availableLeft  = (illuminaRangeLeft.second - overlapLeft.second) / rRatioLeft;
        correctionLeft = illuminaRangeLeft.second - overlapLeft.second;
      } else {
        availableLeft  = (overlapLeft.first - illuminaRangeLeft.first) / rRatioLeft;
        correctionLeft = overlapLeft.first - illuminaRangeLeft.first;
      }

      if (availableLeft > errorOffset) {
        correctionLeft = static_cast<int>(std::floor(errorOffset * rRatioLeft));
        errorOffset    = 0;
      } else {
        errorOffset -= availableLeft;
      }

      if (!matchRight->direction) {
        availableRight  = (overlapRight.first - illuminaRangeRight.first) / rRatioRight;
        correctionRight = overlapRight.first - illuminaRangeRight.first;
      } else {
        availableRight  = (illuminaRangeRight.second - overlapRight.second) / rRatioRight;
        correctionRight = illuminaRangeRight.second - overlapRight.second;
      }

      if (availableRight > errorOffset) {
        correctionRight = static_cast<int>(std::floor(errorOffset * rRatioRight));
      }
    }

    std::string sequence;
    if (!matchRight->direction) {
      sequence = muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdRight, illuminaRangeRight.first + correctionRight,
                                     overlapRight.first, false);
    } else {
      sequence = muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdRight, overlapRight.second,
                                     illuminaRangeRight.second - correctionRight, true);
    }

    sequence.append(
        muchsalsa::getNanoporeSequence(sequenceAccessor, nanoporeId, nanoporeRangeRight.second, nanoporeRangeLeft.first, true));

    if (!matchLeft->direction) {
      sequence.append(muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdLeft, overlapLeft.second,
                                          illuminaRangeLeft.second - correctionLeft, false));
    } else {
      sequence.append(muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdLeft, illuminaRangeLeft.first + correctionLeft,
                                          overlapLeft.first, true));
    }

    return std::make_tuple(sequence.size(), muchsalsa::getReverseComplement(sequence));
  }

  double errorOffset = nanoporeRangeLeft.second - nanoporeRangeRight.first;
  if (errorOffset > 0) {
    auto const correctedNanoporeLeft  = getCorrectedNanoporeRange(matchMap, nanoporeId, illuminaIdLeft, overlapLeft);
    auto const correctedNanoporeRight = getCorrectedNanoporeRange(matchMap, nanoporeId, illuminaIdRight, overlapRight);

    if (correctedNanoporeLeft.second > correctedNanoporeRight.first) {
      return std::make_tuple(std::floor(correctedNanoporeRight.first - correctedNanoporeLeft.second), std::nullopt);
    }

    double availableLeft  = 0;
    double availableRight = 0;
    if (!matchLeft->direction) {
      availableLeft  = (overlapLeft.first - illuminaRangeLeft.first) / rRatioLeft;
      correctionLeft = overlapLeft.first - illuminaRangeLeft.first;
    } else {
      availableLeft  = (illuminaRangeLeft.second - overlapLeft.second) / rRatioLeft;
      correctionLeft = illuminaRangeLeft.second - overlapLeft.second;
    }

    if (availableLeft > errorOffset) {
      correctionLeft = static_cast<int>(std::floor(errorOffset * rRatioLeft));
      errorOffset    = 0;
    } else {
      errorOffset -= availableLeft;
    }

    if (!matchRight->direction) {
      availableRight  = (illuminaRangeRight.second - overlapRight.second) / rRatioRight;
      correctionRight = illuminaRangeRight.second - overlapRight.second;
    } else {
      availableRight  = (overlapRight.first - illuminaRangeRight.first) / rRatioRight;
      correctionRight = overlapRight.first - illuminaRangeRight.first;
    }

    if (availableRight > errorOffset) {
      correctionRight = static_cast<int>(std::floor(errorOffset * rRatioRight));
    }
  }

  std::string sequence;
  if (!matchLeft->direction) {
    sequence = muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdLeft, illuminaRangeLeft.first + correctionLeft,
                                   overlapLeft.first, false);
  } else {
    sequence = muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdLeft, overlapLeft.second,
                                   illuminaRangeLeft.second - correctionLeft, true);
  }

  sequence.append(
      muchsalsa::getNanoporeSequence(sequenceAccessor, nanoporeId, nanoporeRangeLeft.second, nanoporeRangeRight.first, true));

  if (!matchRight->direction) {
    sequence.append(muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdRight, overlapRight.second,
                                        illuminaRangeRight.second - correctionRight, false));
  } else {
    sequence.append(muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaIdRight, illuminaRangeRight.first + correctionRight,
                                        overlapRight.first, true));
  }

  return std::make_tuple(sequence.size(), sequence);
}

void alignAnchorRegion(
    gsl::not_null<std::unordered_map<muchsalsa::graph::Edge const *, std::vector<std::string>> *> pSequences,
    gsl::not_null<std::unordered_map<muchsalsa::graph::Edge const *, int> *>                      pDistances,
    muchsalsa::graph::DiGraph const &adg, muchsalsa::matching::MatchMap const &matchMap,
    std::unordered_map<unsigned int, std::tuple<unsigned int, std::size_t>> const &regIdx2Id,
    muchsalsa::SequenceAccessor &                                                  sequenceAccessor,
    std::unordered_map<muchsalsa::graph::Edge const *, std::vector<muchsalsa::graph::Vertex const *>> const &nanopores,
    muchsalsa::graph::Vertex const &anchorLeft, muchsalsa::graph::Vertex const &anchorRight,
    std::pair<int, int> const &overlapLeft, std::pair<int, int> const &overlapRight) {
  auto const *const pEdge = adg.getEdge(std::make_pair(&anchorLeft, &anchorRight));

  std::optional<int>       globalDistance;
  std::vector<std::string> sequences;
  for (auto const *const pVertex : nanopores.at(pEdge)) {
    auto const sequence = getSequenceBetweenAnchors(
        matchMap, sequenceAccessor, pVertex->getId(), std::get<0>(regIdx2Id.at(anchorLeft.getId())),
        std::get<0>(regIdx2Id.at(anchorRight.getId())), overlapLeft, overlapRight,
        pVertex->getVertexDirection() == muchsalsa::Direction::e_POS);

    if (std::get<1>(sequence).has_value()) {
      sequences.push_back(std::get<1>(sequence).value());
    }

    if (!globalDistance.has_value()) {
      globalDistance = std::get<0>(sequence);
    }
  }

  (*pDistances)[pEdge] = globalDistance.value();
  (*pSequences)[pEdge] = std::move(sequences);
}

} // unnamed namespace

void muchsalsa::assemblePath(
    gsl::not_null<muchsalsa::matching::Id2OverlapMap *> const                               pId2OverlapMap,
    muchsalsa::matching::MatchMap const &                                                   matchMap,
    std::unordered_map<graph::Vertex const *, std::vector<matching::ContainElement>> const &containElements,
    muchsalsa::SequenceAccessor &sequenceAccessor, std::vector<muchsalsa::graph::Vertex const *> const &path,
    muchsalsa::graph::DiGraph const &diGraph, int asmIdx, muchsalsa::OutputWriter &writer) {
  struct Candidate {
    std::set<unsigned int>                 openIds;
    std::set<unsigned int>                 visitedIds;
    std::size_t                            score{0};
    std::size_t                            kinks{0};
    std::vector<graph::Edge const *>       edges;
    std::vector<graph::EdgeOrder const *>  orders;
    std::vector<std::vector<unsigned int>> modifiers;
  };

  std::vector<Candidate> candidates(1);

  auto const findBestCandidate = [](std::size_t const **ppMinKinks, std::size_t const **ppMaxScore,
                                    decltype(candidates) const &candidates) {
    for (auto const &candidate : candidates) {
      if (*ppMinKinks == nullptr || candidate.kinks < **ppMinKinks ||
          (candidate.kinks == **ppMinKinks && (!*ppMaxScore || candidate.score > **ppMaxScore))) {
        *ppMinKinks = &candidate.kinks;
        *ppMaxScore = &candidate.score;
      }
    }
  };

  muchsalsa::Registry                                                     registryAdg;
  muchsalsa::graph::DiGraph                                               adg;
  std::unordered_map<unsigned int, std::tuple<unsigned int, std::size_t>> regIdx2Id;

  for (auto it = std::begin(path); it != std::prev(std::end(path)); ++it) {
    auto pPathEdge = util::make_not_null_and_const(diGraph.getEdge(std::make_pair(*it, *std::next(it))));

    std::vector<Candidate> nextCandidates;
    for (auto const &order : pPathEdge->getEdgeOrders()) {

      std::vector<Candidate> subCandidates;
      for (auto const &candidate : candidates) {
        auto const baseScore = candidate.score + order.score;

        auto ids = std::vector<decltype(order.ids)::value_type>(std::begin(order.ids), std::end(order.ids));
        if (order.baseVertex->getVertexDirection() == Direction::e_NEG) {
          std::reverse(std::begin(ids), std::end(ids));
        }

        auto edgeModifiers = decltype(candidate.modifiers)::value_type();
        edgeModifiers.reserve(ids.size());
        std::copy_if(std::begin(ids), std::end(ids), std::back_inserter(edgeModifiers), [&](auto const id) {
          return candidate.openIds.find(id) == std::end(candidate.openIds) &&
                 candidate.visitedIds.find(id) != std::end(candidate.visitedIds);
        });
        edgeModifiers.shrink_to_fit();

        auto const newKinks = candidate.kinks + edgeModifiers.size();

        auto newVisitedIds = candidate.visitedIds;
        std::copy(std::begin(ids), std::end(ids), std::inserter(newVisitedIds, std::end(newVisitedIds)));

        auto const copyAndAppend = []<template <class, class...> class CONTAINER, class TYPE, class... PARAMS>(
            CONTAINER<TYPE, PARAMS...> const &container, TYPE const &value) {
          auto copy = container;
          copy.insert(std::end(copy), value);

          return copy;
        };
        subCandidates.emplace_back(
            Candidate{std::set<decltype(ids)::value_type>(std::begin(ids), std::end(ids)), std::move(newVisitedIds),
                      baseScore, newKinks, copyAndAppend(candidate.edges, pPathEdge.get()),
                      copyAndAppend(candidate.orders, &order), copyAndAppend(candidate.modifiers, edgeModifiers)});
      }

      std::size_t const *pMinKinks = nullptr;
      std::size_t const *pMaxScore = nullptr;
      findBestCandidate(&pMinKinks, &pMaxScore, subCandidates);
      std::copy_if(std::begin(subCandidates), std::end(subCandidates), std::back_inserter(nextCandidates),
                   [=](auto const &candidate) {
                     return pMinKinks && pMaxScore && candidate.kinks == *pMinKinks && candidate.score == *pMaxScore;
                   });
    }
    candidates = std::move(nextCandidates);
  }

  std::size_t const *pMinKinks = nullptr;
  std::size_t const *pMaxScore = nullptr;
  findBestCandidate(&pMinKinks, &pMaxScore, candidates);

  auto const bestCandidate = *std::find_if(std::begin(candidates), std::end(candidates), [=](auto const &candidate) {
    return pMinKinks && pMaxScore && candidate.kinks == *pMinKinks && candidate.score == *pMaxScore;
  });

  std::unordered_map<unsigned int, std::vector<unsigned int>> clusters;
  for (unsigned int idx = 0; idx < bestCandidate.edges.size(); ++idx) {
    for (auto const &match : bestCandidate.orders[idx]->ids) {
      auto &cluster = clusters[match];
      cluster.push_back(idx);
    }
  }

  std::vector<std::unordered_map<unsigned int, std::size_t>> clusterModifier(bestCandidate.edges.size());
  std::for_each(std::begin(clusters), std::end(clusters), [&](auto const &cluster) {
    getClusterAnchors(&clusterModifier, pId2OverlapMap, matchMap, cluster.first, cluster.second, bestCandidate.edges);
  });

  std::vector<
      std::vector<std::pair<std::pair<int, int>, std::tuple<std::tuple<unsigned int, std::size_t>, std::size_t>>>>
                                                vertexInfo(bestCandidate.edges.size() + 1);
  std::vector<muchsalsa::graph::Vertex const *> vertices(bestCandidate.edges.size() + 1);
  std::unordered_map<unsigned int, std::size_t> matchModifiers;
  for (std::size_t idx = 0; idx < bestCandidate.edges.size(); ++idx) {
    for (auto const &modifier : bestCandidate.modifiers[idx]) {
      ++matchModifiers[modifier];
    }

    auto              ids        = bestCandidate.orders[idx]->ids;
    auto const *const baseVertex = bestCandidate.orders[idx]->baseVertex;
    if (baseVertex->getVertexDirection() == Direction::e_NEG) {
      std::reverse(std::begin(ids), std::end(ids));
    }

    auto const verticesOfEdge = bestCandidate.edges[idx]->getVertices();

    for (auto const &id : ids) {
      auto matches = std::make_tuple(std::make_tuple(id, clusterModifier[idx][id]),
                                     matchModifiers.contains(id) ? matchModifiers[id] : 0);

      auto nanoporeRangeA = matchMap.getVertexMatch(verticesOfEdge.first->getId(), id)->nanoporeRange;
      vertexInfo[idx].push_back(std::make_pair(std::move(nanoporeRangeA), matches));

      auto nanoporeRangeB = matchMap.getVertexMatch(verticesOfEdge.second->getId(), id)->nanoporeRange;
      vertexInfo[idx + 1].push_back(std::make_pair(std::move(nanoporeRangeB), std::move(matches)));
    }

    vertices[idx]     = verticesOfEdge.first;
    vertices[idx + 1] = verticesOfEdge.second;
  }

  std::unordered_map<muchsalsa::graph::Vertex const *, std::string>                                 anchorSequences;
  std::unordered_map<muchsalsa::graph::Edge const *, std::vector<muchsalsa::graph::Vertex const *>> nanopores;
  std::unordered_map<muchsalsa::graph::Vertex const *, std::vector<std::string>>                    preSequences;
  std::unordered_map<muchsalsa::graph::Vertex const *, std::vector<std::string>>                    postSequences;

  for (std::size_t idx = 0; idx < vertices.size(); ++idx) {
    std::sort(std::begin(vertexInfo[idx]), std::end(vertexInfo[idx]), [&](auto const &lhs, auto const &rhs) {
      if (lhs.first == rhs.first) {
        if (!matchMap.getVertexMatch(vertices.at(idx)->getId(), std::get<0>(std::get<0>(lhs.second)))->direction) {
          return pId2OverlapMap->at(std::get<0>(rhs.second)) < pId2OverlapMap->at(std::get<0>(lhs.second));
        }

        return pId2OverlapMap->at(std::get<0>(lhs.second)) < pId2OverlapMap->at(std::get<0>(rhs.second));
      }

      return lhs.first < rhs.first;
    });

    auto &vertexInfoOfInterest = vertexInfo[idx];
    if (vertices.at(idx)->getVertexDirection() == Direction::e_NEG) {
      std::reverse(std::begin(vertexInfoOfInterest), std::end(vertexInfoOfInterest));
    }

    vertexInfoOfInterest = vertexInfo[idx];
    if (vertexInfoOfInterest.empty()) {
      continue;
    }

    auto lastNr    = vertexInfoOfInterest.front().first;
    auto lastMatch = vertexInfoOfInterest.front().second;
    for (auto const &[nr, match] : vertexInfoOfInterest) {

      auto idMatch = tupleTuple2Id(match);
      if (!adg.hasVertex(registryAdg[idMatch])) {
        adg.addVertex(std::make_shared<muchsalsa::graph::Vertex>(registryAdg[idMatch], 0));
        anchorSequences[adg.getVertex(registryAdg[idMatch])] = getAnchorSequence(
            matchMap, sequenceAccessor, vertices.at(idx)->getId(), std::get<0>(std::get<0>(match)),
            pId2OverlapMap->at(std::get<0>(match)), vertices.at(idx)->getVertexDirection() == Direction::e_POS);
        regIdx2Id[registryAdg[idMatch]] = std::get<0>(match);
      }

      if (match == lastMatch) {
        continue;
      }

      auto idLastMatch = tupleTuple2Id(lastMatch);
      if (!adg.hasVertex(registryAdg[idLastMatch])) {
        adg.addVertex(std::make_shared<muchsalsa::graph::Vertex>(registryAdg[idLastMatch], 0));
        anchorSequences[adg.getVertex(registryAdg[idLastMatch])] = getAnchorSequence(
            matchMap, sequenceAccessor, vertices.at(idx)->getId(), std::get<0>(std::get<0>(lastMatch)),
            pId2OverlapMap->at(std::get<0>(lastMatch)), vertices.at(idx)->getVertexDirection() == Direction::e_POS);
        regIdx2Id[registryAdg[idLastMatch]] = std::get<0>(lastMatch);
      }

      auto flip = false;
      if ((lastNr.second > nr.second && lastNr.first < nr.first) ||
          (lastNr.second < nr.second && lastNr.first > nr.first)) {
        auto const cnLeft =
            getCorrectedNanoporeRange(matchMap, vertices.at(idx)->getId(), std::get<0>(std::get<0>(lastMatch)),
                                      pId2OverlapMap->at(std::get<0>(lastMatch)));
        auto const cnRight =
            getCorrectedNanoporeRange(matchMap, vertices.at(idx)->getId(), std::get<0>(std::get<0>(match)),
                                      pId2OverlapMap->at(std::get<0>(match)));

        flip = (vertices.at(idx)->getVertexDirection() == Direction::e_POS &&
                (cnLeft.first > cnRight.first || (cnLeft.first == cnRight.first && cnLeft.second > cnRight.second))) ||
               (vertices.at(idx)->getVertexDirection() == Direction::e_NEG &&
                (cnLeft.first < cnRight.first || (cnLeft.first == cnRight.first && cnLeft.second < cnRight.second)));
      }

      std::pair<muchsalsa::graph::Vertex const *, muchsalsa::graph::Vertex const *> edge;
      if (flip) {
        edge = std::make_pair(adg.getVertex(registryAdg[idMatch]), adg.getVertex(registryAdg[idLastMatch]));
      } else {
        edge = std::make_pair(adg.getVertex(registryAdg[idLastMatch]), adg.getVertex(registryAdg[idMatch]));
      }

      adg.addEdge(edge);

      nanopores[adg.getEdge(std::move(edge))].push_back(vertices.at(idx));

      lastMatch = match;
      lastNr    = nr;
    }

    auto const firstId  = std::get<1>(vertexInfoOfInterest.front());
    auto const secondId = std::get<1>(vertexInfoOfInterest.back());

    auto const *const pFirstVertex = adg.getVertex(registryAdg[tupleTuple2Id(vertexInfoOfInterest.front().second)]);
    preSequences[pFirstVertex].push_back(getSequenceLeftOfAnchor(
        matchMap, sequenceAccessor, vertices.at(idx)->getId(), vertices.at(idx)->getNanoporeLength(),
        std::get<0>(std::get<0>(firstId)), pId2OverlapMap->at(std::get<0>(firstId)),
        vertices.at(idx)->getVertexDirection() == Direction::e_POS));

    auto const *const pSecondVertex = adg.getVertex(registryAdg[tupleTuple2Id(vertexInfoOfInterest.back().second)]);
    postSequences[pSecondVertex].push_back(getSequenceRightOfAnchor(
        matchMap, sequenceAccessor, vertices.at(idx)->getId(), vertices.at(idx)->getNanoporeLength(),
        std::get<0>(std::get<0>(secondId)), pId2OverlapMap->at(std::get<0>(secondId)),
        vertices.at(idx)->getVertexDirection() == Direction::e_POS));
  }

  std::unordered_map<muchsalsa::graph::Edge const *, int>                      distances;
  std::unordered_map<muchsalsa::graph::Edge const *, std::vector<std::string>> sequences;
  for (auto const *const pEdge : adg.getEdges()) {
    auto const verticesOfEdge = pEdge->getVertices();
    alignAnchorRegion(&sequences, &distances, adg, matchMap, regIdx2Id, sequenceAccessor, nanopores,
                      *verticesOfEdge.first, *verticesOfEdge.second,
                      pId2OverlapMap->at(regIdx2Id.at(verticesOfEdge.first->getId())),
                      pId2OverlapMap->at(regIdx2Id.at(verticesOfEdge.second->getId())));
  }

  auto const sortedAdg = adg.sortTopologically();

  std::unordered_map<muchsalsa::graph::Vertex const *, std::size_t> mappingVertex2OrderIdx;
  std::vector<muchsalsa::graph::Vertex const *>                     order;
  order.reserve(sortedAdg.size());

  for (std::size_t idx = 0; idx < sortedAdg.size(); ++idx) {
    mappingVertex2OrderIdx[sortedAdg.at(idx)] = idx;
    order.push_back(sortedAdg.at(idx));
  }

  std::unordered_map<muchsalsa::graph::Vertex const *, bool>                 visitedVertices;
  std::unordered_map<muchsalsa::graph::Vertex const *, std::tuple<int, int>> tap;

  std::optional<std::string> globalSequence;
  int                        globalPos1 = 0;
  int                        globalPos2 = 0;
  std::tie(globalSequence, globalPos1, globalPos2) =
      visitOrdered(&visitedVertices, &tap, adg, regIdx2Id, mappingVertex2OrderIdx, order, distances, sequences,
                   anchorSequences, *pId2OverlapMap, *order.front());

  auto const adgVertices = adg.getVertices();
  if (adgVertices.size() == 1) {
    auto const *const pAnchor = adgVertices.front();
    auto const        overlap = pId2OverlapMap->at(regIdx2Id.at(pAnchor->getId()));

    tap[pAnchor]   = std::make_tuple(0, overlap.second - overlap.first);
    globalSequence = anchorSequences.at(pAnchor);
    globalPos1     = 0;
    globalPos2     = overlap.second - overlap.first;
  }

  std::vector<std::tuple<std::optional<std::string>, int, int, decltype(tap)>> additionalPaths;
  std::vector<bool>                                                            isPathAdded;
  for (auto iter = std::next(std::begin(order)); iter != std::end(order); ++iter) {
    auto const *const pVertex = *iter;

    if (visitedVertices.contains(pVertex)) {
      continue;
    }

    decltype(tap)              localTap;
    std::optional<std::string> localSequence;
    int                        localPos1 = 0;
    int                        localPos2 = 0;
    std::tie(localSequence, localPos1, localPos2) =
        visitOrdered(&visitedVertices, &localTap, adg, regIdx2Id, mappingVertex2OrderIdx, order, distances, sequences,
                     anchorSequences, *pId2OverlapMap, *pVertex);

    if (localTap.empty()) {
      auto const overlap = pId2OverlapMap->at(regIdx2Id.at(pVertex->getId()));

      localTap[pVertex] = std::make_tuple(0, overlap.second - overlap.first);
      localSequence     = anchorSequences.at(pVertex);
      localPos1         = 0;
      localPos2         = overlap.second - overlap.first;
    }

    additionalPaths.emplace_back(std::move(localSequence), localPos1, localPos2, std::move(localTap));
    isPathAdded.push_back(false);
  }

  auto loop = true;
  while (loop) {
    loop = false;

    for (std::size_t idx = 0; idx < additionalPaths.size(); ++idx) {
      auto isFound = false;

      if (isPathAdded[idx]) {
        continue;
      }

      auto localSequence = std::get<0>(additionalPaths.at(idx));
      auto localPos1     = std::get<1>(additionalPaths.at(idx));
      auto localPos2     = std::get<2>(additionalPaths.at(idx));

      int groupOffset = 0;

      auto const &localTap = std::get<3>(additionalPaths.at(idx));
      for (auto const &[pMatch, overlap] : localTap) {
        isFound = false;

        for (auto const &[targetId, pEdge] : adg.getSuccessors(pMatch)) {
          auto const *const pTargetVertex = adg.getVertex(targetId);

          if (tap.contains(pTargetVertex)) {
            groupOffset =
                std::get<0>(tap.at(pTargetVertex)) - distances.at(pEdge) - std::get<1>(localTap.at(pMatch)) - 1;

            if (!sequences.at(pEdge).empty()) {
              std::tie(localSequence, localPos1, localPos2) =
                  updateConsensusBase(localSequence, std::make_pair(localPos1, localPos2), sequences.at(pEdge).front(),
                                      std::make_pair(std::get<1>(localTap.at(pMatch)) + 1,
                                                     std::get<1>(localTap.at(pMatch)) + distances.at(pEdge)));
            }

            isFound = true;
            break;
          }
        }

        if (isFound) {
          break;
        }

        for (auto const &[targetId, pEdge] : adg.getPredecessors(pMatch)) {
          auto const *const pTargetVertex = adg.getVertex(targetId);

          if (tap.contains(pTargetVertex)) {
            groupOffset =
                std::get<1>(tap.at(pTargetVertex)) + distances.at(pEdge) + 1 - std::get<0>(localTap.at(pMatch)) + 1;

            if (!sequences.at(pEdge).empty()) {
              std::tie(localSequence, localPos1, localPos2) =
                  updateConsensusBase(localSequence, std::make_pair(localPos1, localPos2), sequences.at(pEdge).front(),
                                      std::make_pair(std::get<0>(localTap.at(pMatch)) - distances.at(pEdge),
                                                     std::get<0>(localTap.at(pMatch)) - 1));
            }

            isFound = true;
            break;
          }
        }

        if (isFound) {
          break;
        }
      }

      if (!isFound) {
        loop = true;

        continue;
      }

      isPathAdded[idx] = true;
      for (auto const &[pMatch, overlap] : localTap) {
        tap[pMatch] = std::make_tuple(std::get<0>(overlap) + groupOffset, std::get<1>(overlap) + groupOffset);
      }

      std::tie(globalSequence, globalPos1, globalPos2) =
          updateConsensusBase(globalSequence, std::make_pair(globalPos1, globalPos2), localSequence.value(),
                              std::make_pair(localPos1 + groupOffset, localPos2 + groupOffset));
    }
  }

  for (auto const *const pVertex : adg.getVertices()) {
    if (preSequences.contains(pVertex)) {
      auto const iterMaxSeq =
          std::max_element(std::begin(preSequences[pVertex]), std::end(preSequences[pVertex]),
                           [](auto const &lhs, auto const &rhs) { return lhs.size() < rhs.size(); });
      auto const lenMaxSeq                             = static_cast<int>((*iterMaxSeq).size());
      std::tie(globalSequence, globalPos1, globalPos2) = updateConsensusBase(
          globalSequence, std::make_pair(globalPos1, globalPos2), *iterMaxSeq,
          std::make_pair(std::get<0>(tap.at(pVertex)) - lenMaxSeq, std::get<0>(tap.at(pVertex)) - 1));
    }

    if (postSequences.contains(pVertex)) {
      auto const iterMaxSeq =
          std::max_element(std::begin(postSequences[pVertex]), std::end(postSequences[pVertex]),
                           [](auto const &lhs, auto const &rhs) { return lhs.size() < rhs.size(); });
      auto const lenMaxSeq                             = static_cast<int>((*iterMaxSeq).size());
      std::tie(globalSequence, globalPos1, globalPos2) = updateConsensusBase(
          globalSequence, std::make_pair(globalPos1, globalPos2), *iterMaxSeq,
          std::make_pair(std::get<1>(tap.at(pVertex)) + 1, std::get<1>(tap.at(pVertex)) + lenMaxSeq));
    }
  }

  auto       globalLeftMostPosition = globalPos1 * -1;
  auto const targetName             = [&]() {
    std::string targetName = ">muchsalsa_";
    targetName.append(std::to_string(asmIdx));

    return targetName;
  }();

  writer.writeTarget([&]() {
    auto targetSequence = targetName;
    targetSequence.append("\n");
    targetSequence.append(limitLength(globalSequence.value()));
    targetSequence.append("\n");

    return targetSequence;
  }());

  std::size_t queryIdx = 0;
  for (auto const *const pEdge : adg.getEdges()) {

    for (auto const &sequence : sequences.at(pEdge)) {
      if (sequence.empty()) {
        continue;
      }

      auto const querySequenceName = [&]() {
        std::string querySequenceName = ">Middle.";
        querySequenceName.append(std::to_string(asmIdx));
        querySequenceName.append(".");
        querySequenceName.append(std::to_string(queryIdx));

        return querySequenceName;
      }();

      writer.writeQuery([&]() {
        auto querySequence = querySequenceName;
        querySequence.append("\n");
        querySequence.append(limitLength(sequence));
        querySequence.append("\n");

        return querySequence;
      }());

      writer.writePaf([&]() {
        int const sequenceLength = static_cast<int>(sequence.size());

        auto const anchors = pEdge->getVertices();
        auto const lb      = std::get<1>(tap.at(anchors.first)) + 1 + globalLeftMostPosition;
        auto const rb      = std::get<0>(tap.at(anchors.second)) - 1 + globalLeftMostPosition;

        std::string paf = querySequenceName.substr(1);
        paf.append("\t");
        paf.append(std::to_string(sequenceLength));
        paf.append("\t0\t");
        paf.append(std::to_string(sequenceLength));
        paf.append("\t+\t");
        paf.append(targetName.substr(1));
        paf.append("\t");
        paf.append(std::to_string(globalSequence.value().size()));
        paf.append("\t");
        paf.append(std::to_string(lb));
        paf.append("\t");
        paf.append(std::to_string(rb));
        paf.append("\t");
        paf.append(std::to_string(rb - lb + 1));
        paf.append("\t");
        paf.append(std::to_string(rb - lb + 1));
        paf.append("\t255");
        paf.append("\n");

        return paf;
      }());

      ++queryIdx;
    }
  }

  for (auto const *const pVertex : adg.getVertices()) {
    if (preSequences.contains(pVertex)) {
      for (auto const &sequence : preSequences[pVertex]) {
        if (sequence.size() < TH_SEQUENCE_LENGTH) {
          continue;
        }

        auto const querySequenceName = [&]() {
          std::string querySequenceName = ">Left.";
          querySequenceName.append(std::to_string(asmIdx));
          querySequenceName.append(".");
          querySequenceName.append(std::to_string(queryIdx));

          return querySequenceName;
        }();

        writer.writeQuery([&]() {
          auto querySequence = querySequenceName;
          querySequence.append("\n");
          querySequence.append(limitLength(sequence));
          querySequence.append("\n");

          return querySequence;
        }());

        writer.writePaf([&]() {
          int const sequenceLength = static_cast<int>(sequence.size());

          auto const rb = std::get<0>(tap.at(pVertex)) - 1 + globalLeftMostPosition;
          auto const lb = rb - sequenceLength + 1;

          std::string paf = querySequenceName.substr(1);
          paf.append("\t");
          paf.append(std::to_string(sequenceLength));
          paf.append("\t0\t");
          paf.append(std::to_string(sequenceLength));
          paf.append("\t+\t");
          paf.append(targetName.substr(1));
          paf.append("\t");
          paf.append(std::to_string(globalSequence.value().size()));
          paf.append("\t");
          paf.append(std::to_string(lb));
          paf.append("\t");
          paf.append(std::to_string(rb));
          paf.append("\t");
          paf.append(std::to_string(rb - lb + 1));
          paf.append("\t");
          paf.append(std::to_string(rb - lb + 1));
          paf.append("\t255");
          paf.append("\n");

          return paf;
        }());

        ++queryIdx;
      }
    }

    if (postSequences.contains(pVertex)) {
      for (auto const &sequence : postSequences[pVertex]) {
        if (sequence.size() < TH_SEQUENCE_LENGTH) {
          continue;
        }

        auto const querySequenceName = [&]() {
          std::string querySequenceName = ">Right.";
          querySequenceName.append(std::to_string(asmIdx));
          querySequenceName.append(".");
          querySequenceName.append(std::to_string(queryIdx));

          return querySequenceName;
        }();

        writer.writeQuery([&]() {
          auto querySequence = querySequenceName;
          querySequence.append("\n");
          querySequence.append(limitLength(sequence));
          querySequence.append("\n");

          return querySequence;
        }());

        writer.writePaf([&]() {
          int const sequenceLength = static_cast<int>(sequence.size());

          auto const lb = std::get<1>(tap.at(pVertex)) + 1 + globalLeftMostPosition;
          auto const rb = lb + sequenceLength - 1;

          std::string paf = querySequenceName.substr(1);
          paf.append("\t");
          paf.append(std::to_string(sequenceLength));
          paf.append("\t0\t");
          paf.append(std::to_string(sequenceLength));
          paf.append("\t+\t");
          paf.append(targetName.substr(1));
          paf.append("\t");
          paf.append(std::to_string(globalSequence.value().size()));
          paf.append("\t");
          paf.append(std::to_string(lb));
          paf.append("\t");
          paf.append(std::to_string(rb));
          paf.append("\t");
          paf.append(std::to_string(rb - lb + 1));
          paf.append("\t");
          paf.append(std::to_string(rb - lb + 1));
          paf.append("\t255");
          paf.append("\n");

          return paf;
        }());

        ++queryIdx;
      }
    }
  }

  for (std::size_t idx = 0; idx < vertices.size(); ++idx) {
    std::unordered_map<unsigned int, std::tuple<std::tuple<unsigned int, std::size_t>, std::size_t>> mappingId2Anchor;
    for (auto const &info : vertexInfo.at(idx)) {
      mappingId2Anchor[std::get<0>(std::get<0>(info.second))] = info.second;
    }

    if (!containElements.contains(vertices.at(idx))) {
      continue;
    }

    for (auto const &containElement : containElements.at(vertices.at(idx))) {
      std::vector<std::tuple<std::pair<int, int>, unsigned int>> containInfo;
      containInfo.reserve(containElement.matches.size());

      for (auto const &match : containElement.matches) {
        if (mappingId2Anchor.contains(match.first)) {
          containInfo.emplace_back(match.second->nanoporeRange, match.first);
        }
      }

      if (containInfo.empty()) {
        continue;
      }

      std::sort(std::begin(containInfo), std::end(containInfo));

      auto const direction = containElement.direction * (vertices.at(idx)->getVertexDirection() == Direction::e_POS);
      if (!direction) {
        std::reverse(std::begin(containInfo), std::end(containInfo));
      }

      std::vector<std::tuple<int, int>> globalRanges;
      globalRanges.reserve(containInfo.size());

      for (auto const &info : containInfo) {
        auto const tapId  = mappingId2Anchor.at(std::get<1>(info));
        auto const tapDir = matchMap.getVertexMatch(vertices.at(idx)->getId(), std::get<1>(info))->direction *
                            (vertices.at(idx)->getVertexDirection() == Direction::e_POS);
        auto const illuminaRef =
            tapDir ? pId2OverlapMap->at(std::get<0>(tapId)).second : pId2OverlapMap->at(std::get<0>(tapId)).first;

        auto const totalRef =
            std::get<1>(tap.at(adg.getVertex(registryAdg[tupleTuple2Id(tapId)]))) + globalLeftMostPosition;

        auto const contDir       = containElement.matches.at(std::get<1>(info))->direction * direction;
        auto const illuminaRange = containElement.matches.at(std::get<1>(info))->illuminaRange;
        if (!contDir) {
          auto const offset = illuminaRange.first - illuminaRef;
          globalRanges.emplace_back(totalRef - offset - (illuminaRange.second - illuminaRange.first),
                                    totalRef - offset);
        } else {
          auto const offset = illuminaRange.second - illuminaRef;
          globalRanges.emplace_back(totalRef + offset - (illuminaRange.second - illuminaRange.first),
                                    totalRef + offset);
        }
      }

      std::vector<std::tuple<std::string, int, int, const char *>> sequences2Write;
      for (std::size_t idxGlobalRange = 0; idxGlobalRange < globalRanges.size(); ++idxGlobalRange) {
        auto const &illuminaId = std::get<1>(containInfo[idxGlobalRange]);
        auto const &match      = containElement.matches.at(illuminaId);
        sequences2Write.emplace_back(muchsalsa::getIlluminaSequence(sequenceAccessor, illuminaId, match->illuminaRange.first,
                                                         match->illuminaRange.second, match->direction * direction),
                                     std::get<0>(globalRanges.at(idxGlobalRange)),
                                     std::get<1>(globalRanges.at(idxGlobalRange)), "Illumina_Match");

        if (idxGlobalRange == 0) {
          continue;
        }

        auto const preNanopore = containElement.matches.at(std::get<1>(containInfo[idxGlobalRange - 1]))->nanoporeRange;
        sequences2Write.emplace_back(muchsalsa::getNanoporeSequence(sequenceAccessor, containElement.nano, preNanopore.second + 1,
                                                         match->nanoporeRange.first - 1, direction),
                                     std::get<1>(globalRanges.at(idxGlobalRange - 1)) + 1,
                                     std::get<0>(globalRanges.at(idxGlobalRange)) - 1, "Nano_Middle");
      }

      for (auto const &sequence : sequences2Write) {
        if (std::get<0>(sequence).size() < TH_SEQUENCE_LENGTH) {
          continue;
        }

        auto const querySequenceName = [&]() {
          std::string querySequenceName = ">Contain_";
          querySequenceName.append(std::get<3>(sequence));
          querySequenceName.append(".");
          querySequenceName.append(std::to_string(asmIdx));
          querySequenceName.append(".");
          querySequenceName.append(std::to_string(queryIdx));

          return querySequenceName;
        }();

        writer.writeQuery([&]() {
          auto querySequence = querySequenceName;
          querySequence.append("\n");
          querySequence.append(limitLength(std::get<0>(sequence)));
          querySequence.append("\n");

          return querySequence;
        }());

        writer.writePaf([&]() {
          auto const sequenceLength = std::get<0>(sequence).size();

          auto const lb = std::get<1>(sequence);
          auto const rb = std::get<2>(sequence);

          std::string paf = querySequenceName.substr(1);
          paf.append("\t");
          paf.append(std::to_string(sequenceLength));
          paf.append("\t0\t");
          paf.append(std::to_string(sequenceLength));
          paf.append("\t+\t");
          paf.append(targetName.substr(1));
          paf.append("\t");
          paf.append(std::to_string(globalSequence.value().size()));
          paf.append("\t");
          paf.append(std::to_string(lb));
          paf.append("\t");
          paf.append(std::to_string(rb));
          paf.append("\t");
          paf.append(std::to_string(rb - lb + 1));
          paf.append("\t");
          paf.append(std::to_string(rb - lb + 1));
          paf.append("\t255");
          paf.append("\n");

          return paf;
        }());

        ++queryIdx;
      }
    }
  }
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
