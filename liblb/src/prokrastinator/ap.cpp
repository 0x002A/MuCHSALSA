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

#include <algorithm>
#include <cmath>
#include <iterator>
#include <memory>
#include <set>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "OutputWriter.h"
#include "SequenceAccessor.h"
#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/Id2OverlapMap.h"
#include "matching/MatchMap.h"
#include "types/Direction.h"

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr auto SEQUENCE_LINE_LENGTH = 60;
constexpr auto TH_SEQUENCE_LENGTH   = 200;

// =====================================================================================================================
//                                                        METHODS
// =====================================================================================================================

namespace {

std::string_view strSlice(std::string_view original, int intStart, int intEnd) {
  auto const doSlicing = [=](int i, int j) {
    std::size_t intStart = static_cast<std::size_t>(std::max(0, i));
    std::size_t intEnd =
        std::max(std::min(original.size(), static_cast<std::size_t>(std::max(0, j))), static_cast<std::size_t>(i));

    return original.substr(intStart, intEnd - intStart);
  };

  int const size = static_cast<int>(original.size());
  return doSlicing(intStart >= 0 ? intStart : size + intStart, intEnd >= 0 ? intEnd : size + intEnd);
}

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

std::tuple<std::tuple<std::string, std::size_t>, std::size_t> id2TupleTuple(std::string const &id) {
  auto copy = id;

  std::string const delimiter = ",";
  auto              pos       = copy.find(delimiter);

  auto first = copy.substr(0, pos);
  copy.erase(0, pos + 1);
  pos = copy.find(delimiter);

  auto second = copy.substr(0, pos);
  copy.erase(0, pos + 1);

  return std::make_tuple(std::make_tuple(first, std::stoul(second)), std::stoul(copy));
}

std::string tupleTuple2Id(std::tuple<std::tuple<std::string, std::size_t>, std::size_t> const &tupleTuple) {
  std::string id;

  auto const &innerTuple = std::get<0>(tupleTuple);
  id.append(std::get<0>(innerTuple));
  id.append(",");
  id.append(std::to_string(std::get<1>(innerTuple)));
  id.append(",");
  id.append(std::to_string(std::get<1>(tupleTuple)));

  return id;
}

std::vector<lazybastard::graph::Vertex *> ramseyR2(lazybastard::graph::Graph const &         graph,
                                                   std::vector<lazybastard::graph::Vertex *> vertices) {
  if (vertices.empty()) {
    return std::vector<lazybastard::graph::Vertex *>();
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

std::vector<std::vector<lazybastard::graph::Vertex *>> getAnchorCliques(lazybastard::graph::Graph const &graph) {
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
    gsl::not_null<std::vector<std::unordered_map<std::string, std::size_t>> *> const pClusterModifier,
    gsl::not_null<lazybastard::matching::MatchMap const *> pMatchMap, std::string const &illuminaIdBase,
    std::vector<std::size_t> const &                                           edgeIdx,
    gsl::not_null<std::vector<lazybastard::graph::Edge const *> const *> const pEdges,
    gsl::not_null<lazybastard::matching::Id2OverlapMap *> const                pId2OverlapMap) {

  lazybastard::graph::Graph g;
  for (auto const idxEdge1 : edgeIdx) {
    g.addVertex(std::make_shared<lazybastard::graph::Vertex>(std::to_string(idxEdge1), 0));

    for (auto const idxEdge2 : edgeIdx) {
      if (idxEdge1 == idxEdge2) {
        break;
      }

      auto const overlapEdge1 = pMatchMap->getEdgeMatch((*pEdges)[idxEdge1], illuminaIdBase)->overlap;
      auto const overlapEdge2 = pMatchMap->getEdgeMatch((*pEdges)[idxEdge2], illuminaIdBase)->overlap;

      auto const overlap = std::make_pair(std::max(overlapEdge1.first, overlapEdge2.first),
                                          std::min(overlapEdge1.second, overlapEdge2.second));

      if (overlap.first <= overlap.second) {
        g.addEdge(std::make_pair(g.getVertex(std::to_string(idxEdge2)), g.getVertex(std::to_string(idxEdge1))));
      }
    }
  }

  auto const anchorCliques = getAnchorCliques(g);

  for (std::size_t idx = 0; idx < anchorCliques.size(); ++idx) {
    auto const correctedId = std::make_tuple(illuminaIdBase, idx);

    std::optional<std::pair<int, int>> commonOverlap;
    for (auto const *const pVertex : anchorCliques.at(idx)) {
      (*pClusterModifier)[std::stoul(pVertex->getId())][illuminaIdBase] = idx;

      if (!commonOverlap) {
        commonOverlap = pMatchMap->getEdgeMatch((*pEdges)[std::stoul(pVertex->getId())], illuminaIdBase)->overlap;
      } else {
        auto const otherOverlap =
            pMatchMap->getEdgeMatch((*pEdges)[std::stoul(pVertex->getId())], illuminaIdBase)->overlap;

        commonOverlap = std::make_pair(std::max(commonOverlap->first, otherOverlap.first),
                                       std::min(commonOverlap->second, otherOverlap.second));
      }
    }

    (*pId2OverlapMap)[correctedId] = commonOverlap.value();
  }
}

std::pair<double, double> getCorrectedNanoporeRange(lazybastard::matching::MatchMap const &matchMap,
                                                    std::string const &nanoporeId, std::string const &illuminaId,
                                                    std::pair<int, int> const &overlap) {
  auto const *const pMatch = matchMap.getVertexMatch(nanoporeId, illuminaId);

  auto nanoCorrectionLeft  = (overlap.first - pMatch->illuminaRange.first) / pMatch->rRatio;
  auto nanoCorrectionRight = (pMatch->illuminaRange.second - overlap.second) / pMatch->rRatio;

  lazybastard::util::swap_if(nanoCorrectionLeft, nanoCorrectionRight, !pMatch->direction);

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

    updatedSequence.append([=]() { return strSlice(newSequence, 0, borderRight); }());
    updatedSequence.append(oldSequence.value());
  } else if (newBorders.second > oldBorders.second) {
    updatedSequence.append(oldSequence.value());

    auto const borderLeft = -(newBorders.second - oldBorders.second);
    updatedSequence.append(strSlice(newSequence, borderLeft, static_cast<int>(newSequence.size())));
  } else {
    updatedSequence = oldSequence.value();
  }

  return std::make_tuple(std::move(updatedSequence), std::min(oldBorders.first, newBorders.first),
                         std::max(oldBorders.second, newBorders.second));
}

std::tuple<std::optional<std::string>, int, int>
visitOrdered(gsl::not_null<std::unordered_map<lazybastard::graph::Vertex const *, bool> *> const pVisitedVertices,
             gsl::not_null<std::unordered_map<lazybastard::graph::Vertex const *, std::tuple<int, int>> *> pTap,
             lazybastard::graph::DiGraph const &                                                           adg,
             std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> const &mappingVertex2OrderIdx,
             std::vector<lazybastard::graph::Vertex const *> const &                    order,
             std::unordered_map<lazybastard::graph::Edge const *, int> const &          distances,
             std::unordered_map<lazybastard::graph::Edge const *, std::vector<std::string>> const &sequences,
             std::unordered_map<lazybastard::graph::Vertex const *, std::string> const &           anchorSequences,
             lazybastard::matching::Id2OverlapMap const &id2OverlapMap, lazybastard::graph::Vertex const &vertexStart) {
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

        auto const overlapLeft  = id2OverlapMap.at(std::get<0>(id2TupleTuple(pAnchorLeft->getId())));
        auto const overlapRight = id2OverlapMap.at(std::get<0>(id2TupleTuple(pAnchorRight->getId())));

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
          (*pTap)[pAnchorLeft] = std::make_tuple(posRight - offset - 2 - lengthLeft, posRight - offset - 1);

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

std::string getReverseComplement(std::string const &sequence) {
  std::string reverseComplement;
  reverseComplement.reserve(sequence.size());

  std::transform(std::rbegin(sequence), std::rend(sequence), std::back_inserter(reverseComplement), [](auto const &c) {
    switch (c) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    default:
      return c;
    }
  });

  return reverseComplement;
}

std::string getIlluminaSequence(lazybastard::SequenceAccessor &sequenceAccessor, std::string const &illuminaId,
                                int illuminaOverlapLeft, int illuminaOverlapRight,
                                lazybastard::Toggle const direction) {
  auto illuminaSequence = std::string(
      strSlice(sequenceAccessor.getIlluminaSequence(illuminaId), illuminaOverlapLeft, illuminaOverlapRight + 1));

  if (!direction) {
    return getReverseComplement(illuminaSequence);
  }

  return illuminaSequence;
}

std::string getNanoporeSequence(lazybastard::SequenceAccessor &sequenceAccessor, std::string const &nanoporeId,
                                int nanoporeRegionLeft, int nanoporeRegionRight, lazybastard::Toggle const direction) {
  auto nanoporeSequence = std::string(
      strSlice(sequenceAccessor.getNanoporeSequence(nanoporeId), nanoporeRegionLeft, nanoporeRegionRight + 1));

  if (!direction) {
    return getReverseComplement(nanoporeSequence);
  }

  return nanoporeSequence;
}

std::string getSequenceLeftOfAnchor(lazybastard::matching::MatchMap const &matchMap,
                                    lazybastard::SequenceAccessor &sequenceAccessor, const std::string &nanoporeId,
                                    std::size_t nanoporeLength, const std::string &illuminaId,
                                    std::pair<int, int> illuminaOverlap, lazybastard::Toggle const direction) {
  auto const *const pMatch = matchMap.getVertexMatch(nanoporeId, illuminaId);

  if (!direction) {
    std::string illuminaSequence;

    if (!pMatch->direction) {
      illuminaSequence =
          getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, false);
    } else {
      illuminaSequence =
          getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second, pMatch->illuminaRange.second, true);
    }

    illuminaSequence.append(getNanoporeSequence(sequenceAccessor, nanoporeId, pMatch->nanoporeRange.second,
                                                static_cast<int>(nanoporeLength) - 1, true));

    return getReverseComplement(illuminaSequence);
  }

  auto nanoporeSequence = getNanoporeSequence(sequenceAccessor, nanoporeId, 0, pMatch->nanoporeRange.first, true);

  if (!pMatch->direction) {
    nanoporeSequence.append(
        getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second, pMatch->illuminaRange.second, false));
  } else {
    nanoporeSequence.append(
        getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, true));
  }

  return nanoporeSequence;
}

std::string getSequenceRightOfAnchor(lazybastard::matching::MatchMap const &matchMap,
                                     lazybastard::SequenceAccessor &sequenceAccessor, const std::string &nanoporeId,
                                     std::size_t nanoporeLength, const std::string &illuminaId,
                                     std::pair<int, int> illuminaOverlap, lazybastard::Toggle const direction) {
  auto const *const pMatch = matchMap.getVertexMatch(nanoporeId, illuminaId);

  if (!direction) {
    auto nanoporeSequence = getNanoporeSequence(sequenceAccessor, nanoporeId, 0, pMatch->nanoporeRange.first, true);

    if (!pMatch->direction) {
      nanoporeSequence.append(getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second,
                                                  pMatch->illuminaRange.second, false));
    } else {
      nanoporeSequence.append(
          getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, true));
    }

    return getReverseComplement(nanoporeSequence);
  }

  std::string illuminaSequence;

  if (!pMatch->direction) {
    illuminaSequence =
        getIlluminaSequence(sequenceAccessor, illuminaId, pMatch->illuminaRange.first, illuminaOverlap.first, false);
  } else {
    illuminaSequence =
        getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.second, pMatch->illuminaRange.second, true);
  }

  illuminaSequence.append(getNanoporeSequence(sequenceAccessor, nanoporeId, pMatch->nanoporeRange.second,
                                              static_cast<int>(nanoporeLength) - 1, true));

  return illuminaSequence;
}

std::string getAnchorSequence(lazybastard::matching::MatchMap const &matchMap,
                              lazybastard::SequenceAccessor &sequenceAccessor, std::string const &nanoporeId,
                              std::string const &illuminaId, std::pair<int, int> illuminaOverlap,
                              lazybastard::Toggle const direction) {
  auto const combinedDirection = matchMap.getVertexMatch(nanoporeId, illuminaId)->direction * direction;
  auto       illuminaSequence  = getIlluminaSequence(sequenceAccessor, illuminaId, illuminaOverlap.first,
                                              illuminaOverlap.second, combinedDirection);

  return illuminaSequence;
}

std::tuple<int, std::optional<std::string>>
getSequenceBetweenAnchors(lazybastard::matching::MatchMap const &matchMap,
                          lazybastard::SequenceAccessor &sequenceAccessor, std::string const &nanoporeId,
                          std::string const &illuminaIdLeft, std::string const &illuminaIdRight,
                          std::pair<int, int> const &overlapLeft, std::pair<int, int> const &overlapRight,
                          lazybastard::Toggle const direction) {
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
      sequence = getIlluminaSequence(sequenceAccessor, illuminaIdRight, illuminaRangeRight.first + correctionRight,
                                     overlapRight.first, false);
    } else {
      sequence = getIlluminaSequence(sequenceAccessor, illuminaIdRight, overlapRight.second,
                                     illuminaRangeRight.second - correctionRight, true);
    }

    sequence.append(
        getNanoporeSequence(sequenceAccessor, nanoporeId, nanoporeRangeRight.second, nanoporeRangeLeft.first, true));

    if (!matchLeft->direction) {
      sequence.append(getIlluminaSequence(sequenceAccessor, illuminaIdLeft, overlapLeft.second,
                                          illuminaRangeLeft.second - correctionLeft, false));
    } else {
      sequence.append(getIlluminaSequence(sequenceAccessor, illuminaIdLeft, illuminaRangeLeft.first + correctionLeft,
                                          overlapLeft.first, true));
    }

    return std::make_tuple(sequence.size(), getReverseComplement(sequence));
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
    sequence = getIlluminaSequence(sequenceAccessor, illuminaIdLeft, illuminaRangeLeft.first + correctionLeft,
                                   overlapLeft.first, false);
  } else {
    sequence = getIlluminaSequence(sequenceAccessor, illuminaIdLeft, overlapLeft.second,
                                   illuminaRangeLeft.second - correctionLeft, true);
  }

  sequence.append(
      getNanoporeSequence(sequenceAccessor, nanoporeId, nanoporeRangeLeft.second, nanoporeRangeRight.first, true));

  if (!matchRight->direction) {
    sequence.append(getIlluminaSequence(sequenceAccessor, illuminaIdRight, overlapRight.second,
                                        illuminaRangeRight.second - correctionRight, false));
  } else {
    sequence.append(getIlluminaSequence(sequenceAccessor, illuminaIdRight, illuminaRangeRight.first + correctionRight,
                                        overlapRight.first, true));
  }

  return std::make_tuple(sequence.size(), sequence);
}

void alignAnchorRegion(
    gsl::not_null<std::unordered_map<lazybastard::graph::Edge const *, std::vector<std::string>> *> pSequences,
    gsl::not_null<std::unordered_map<lazybastard::graph::Edge const *, int> *>                      pDistances,
    lazybastard::graph::DiGraph const &adg, lazybastard::matching::MatchMap const &matchMap,
    lazybastard::SequenceAccessor &sequenceAccessor,
    std::unordered_map<lazybastard::graph::Edge const *, std::vector<lazybastard::graph::Vertex const *>> const
        &                             nanopores,
    lazybastard::graph::Vertex const &anchorLeft, lazybastard::graph::Vertex const &anchorRight,
    std::pair<int, int> const &overlapLeft, std::pair<int, int> const &overlapRight) {
  auto const *const pEdge = adg.getEdge(std::make_pair(&anchorLeft, &anchorRight));

  std::optional<int>       globalDistance;
  std::vector<std::string> sequences;
  for (auto const *const pVertex : nanopores.at(pEdge)) {
    auto const sequence = getSequenceBetweenAnchors(
        matchMap, sequenceAccessor, pVertex->getId(), std::get<0>(std::get<0>(id2TupleTuple(anchorLeft.getId()))),
        std::get<0>(std::get<0>(id2TupleTuple(anchorRight.getId()))), overlapLeft, overlapRight,
        pVertex->getVertexDirection() == lazybastard::Direction::e_POS);

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

void lazybastard::assemblePath(
    gsl::not_null<lazybastard::matching::MatchMap const *> const pMatchMap,
    gsl::not_null<std::unordered_map<graph::Vertex const *, std::vector<matching::ContainElement>> const *> const
                                                                                 pContainElements,
    gsl::not_null<lazybastard::SequenceAccessor *> const                         pSequenceAccessor,
    gsl::not_null<lazybastard::matching::Id2OverlapMap *> const                  pId2OverlapMap,
    gsl::not_null<std::vector<lazybastard::graph::Vertex const *> const *> const pPath,
    gsl::not_null<lazybastard::graph::DiGraph const *> const pDiGraph, int asmIdx, lazybastard::OutputWriter &writer) {
  struct Candidate {
    std::set<std::string>                 openIds;
    std::set<std::string>                 visitedIds;
    std::size_t                           score{0};
    std::size_t                           kinks{0};
    std::vector<graph::Edge const *>      edges;
    std::vector<graph::EdgeOrder const *> orders;
    std::vector<std::vector<std::string>> modifiers;
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

  lazybastard::graph::DiGraph adg;

  for (auto it = std::begin(*pPath); it != std::prev(std::end(*pPath)); ++it) {
    auto pPathEdge = util::make_not_null_and_const(pDiGraph->getEdge(std::make_pair(*it, *std::next(it))));

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

  std::unordered_map<std::string, std::vector<std::size_t>> clusters;
  for (std::size_t idx = 0; idx < bestCandidate.edges.size(); ++idx) {
    for (auto const &match : bestCandidate.orders[idx]->ids) {
      auto &cluster = clusters[match];
      cluster.push_back(idx);
    }
  }

  std::vector<std::unordered_map<std::string, std::size_t>> clusterModifier(bestCandidate.edges.size());
  std::for_each(std::begin(clusters), std::end(clusters), [&](auto const &cluster) {
    getClusterAnchors(&clusterModifier, pMatchMap, cluster.first, cluster.second, &bestCandidate.edges, pId2OverlapMap);
  });

  std::vector<
      std::vector<std::pair<std::pair<int, int>, std::tuple<std::tuple<std::string, std::size_t>, std::size_t>>>>
                                                  vertexInfo(bestCandidate.edges.size() + 1);
  std::vector<lazybastard::graph::Vertex const *> vertices(bestCandidate.edges.size() + 1);
  std::unordered_map<std::string, std::size_t>    matchModifiers;
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

      auto nanoporeRangeA = pMatchMap->getVertexMatch(verticesOfEdge.first->getId(), id)->nanoporeRange;
      vertexInfo[idx].push_back(std::make_pair(std::move(nanoporeRangeA), matches));

      auto nanoporeRangeB = pMatchMap->getVertexMatch(verticesOfEdge.second->getId(), id)->nanoporeRange;
      vertexInfo[idx + 1].push_back(std::make_pair(std::move(nanoporeRangeB), std::move(matches)));
    }

    vertices[idx]     = verticesOfEdge.first;
    vertices[idx + 1] = verticesOfEdge.second;
  }

  std::unordered_map<lazybastard::graph::Vertex const *, std::string>                                   anchorSequences;
  std::unordered_map<lazybastard::graph::Edge const *, std::vector<lazybastard::graph::Vertex const *>> nanopores;
  std::unordered_map<lazybastard::graph::Vertex const *, std::vector<std::string>>                      preSequences;
  std::unordered_map<lazybastard::graph::Vertex const *, std::vector<std::string>>                      postSequences;

  for (std::size_t idx = 0; idx < vertices.size(); ++idx) {
    std::sort(std::begin(vertexInfo[idx]), std::end(vertexInfo[idx]), [=](auto const &lhs, auto const &rhs) {
      if (lhs.first == rhs.first) {
        if (!pMatchMap->getVertexMatch(vertices.at(idx)->getId(), std::get<0>(std::get<0>(lhs.second)))->direction) {
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

    // std::pair<std::pair<int, int>, std::tuple<std::tuple<std::string, std::size_t>, std::size_t>
    std::string values;
    for (auto const &value : vertexInfo[idx]) {
      values.append(std::to_string(value.first.first) + "," + std::to_string(value.first.second) + ",");
      values.append(std::get<0>(std::get<0>(value.second)) + "," +
                    std::to_string(std::get<1>(std::get<0>(value.second))));
      values.append(";");
    }
    // std::cout << vertices.at(idx)->getId() << ";" << (vertices.at(idx)->getVertexDirection() == Direction::e_NEG ?
    // "-1" : "1") << ";" << values << std::endl;

    vertexInfoOfInterest = vertexInfo[idx];
    if (vertexInfoOfInterest.empty()) {
      continue;
    }

    auto lastNr    = vertexInfoOfInterest.front().first;
    auto lastMatch = vertexInfoOfInterest.front().second;
    for (auto const &[nr, match] : vertexInfoOfInterest) {

      auto idMatch = tupleTuple2Id(match);
      if (!adg.hasVertex(idMatch)) {
        adg.addVertex(std::make_shared<lazybastard::graph::Vertex>(idMatch, 0));
        anchorSequences[adg.getVertex(idMatch)] = getAnchorSequence(
            *pMatchMap, *pSequenceAccessor, vertices.at(idx)->getId(), std::get<0>(std::get<0>(match)),
            pId2OverlapMap->at(std::get<0>(match)), vertices.at(idx)->getVertexDirection() == Direction::e_POS);
      }

      if (match == lastMatch) {
        continue;
      }

      auto idLastMatch = tupleTuple2Id(lastMatch);
      if (!adg.hasVertex(idLastMatch)) {
        adg.addVertex(std::make_shared<lazybastard::graph::Vertex>(idLastMatch, 0));
        anchorSequences[adg.getVertex(idLastMatch)] = getAnchorSequence(
            *pMatchMap, *pSequenceAccessor, vertices.at(idx)->getId(), std::get<0>(std::get<0>(lastMatch)),
            pId2OverlapMap->at(std::get<0>(lastMatch)), vertices.at(idx)->getVertexDirection() == Direction::e_POS);
      }

      auto flip = false;
      if ((lastNr.second > nr.second && lastNr.first < nr.first) ||
          (lastNr.second < nr.second && lastNr.first > nr.first)) {
        auto const cnLeft =
            getCorrectedNanoporeRange(*pMatchMap, vertices.at(idx)->getId(), std::get<0>(std::get<0>(lastMatch)),
                                      pId2OverlapMap->at(std::get<0>(lastMatch)));
        auto const cnRight =
            getCorrectedNanoporeRange(*pMatchMap, vertices.at(idx)->getId(), std::get<0>(std::get<0>(match)),
                                      pId2OverlapMap->at(std::get<0>(match)));

        flip = (vertices.at(idx)->getVertexDirection() == Direction::e_POS &&
                (cnLeft.first > cnRight.first || (cnLeft.first == cnRight.first && cnLeft.second > cnRight.second))) ||
               (vertices.at(idx)->getVertexDirection() == Direction::e_NEG &&
                (cnLeft.first < cnRight.first || (cnLeft.first == cnRight.first && cnLeft.second < cnRight.second)));
      }

      std::pair<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *> edge;
      if (flip) {
        edge = std::make_pair(adg.getVertex(tupleTuple2Id(match)), adg.getVertex(tupleTuple2Id(lastMatch)));
      } else {
        edge = std::make_pair(adg.getVertex(tupleTuple2Id(lastMatch)), adg.getVertex(tupleTuple2Id(match)));
      }

      adg.addEdge(edge);

      nanopores[adg.getEdge(std::move(edge))].push_back(vertices.at(idx));

      lastMatch = match;
      lastNr    = nr;
    }

    auto const firstId  = std::get<1>(vertexInfoOfInterest.front());
    auto const secondId = std::get<1>(vertexInfoOfInterest.back());

    auto const *const pFirstVertex = adg.getVertex(tupleTuple2Id(vertexInfoOfInterest.front().second));
    preSequences[pFirstVertex].push_back(getSequenceLeftOfAnchor(
        *pMatchMap, *pSequenceAccessor, vertices.at(idx)->getId(), vertices.at(idx)->getNanoporeLength(),
        std::get<0>(std::get<0>(firstId)), pId2OverlapMap->at(std::get<0>(firstId)),
        vertices.at(idx)->getVertexDirection() == Direction::e_POS));

    auto const *const pSecondVertex = adg.getVertex(tupleTuple2Id(vertexInfoOfInterest.back().second));
    postSequences[pSecondVertex].push_back(getSequenceRightOfAnchor(
        *pMatchMap, *pSequenceAccessor, vertices.at(idx)->getId(), vertices.at(idx)->getNanoporeLength(),
        std::get<0>(std::get<0>(secondId)), pId2OverlapMap->at(std::get<0>(secondId)),
        vertices.at(idx)->getVertexDirection() == Direction::e_POS));
  }

  std::unordered_map<lazybastard::graph::Edge const *, int>                      distances;
  std::unordered_map<lazybastard::graph::Edge const *, std::vector<std::string>> sequences;
  for (auto const *const pEdge : adg.getEdges()) {
    auto const verticesOfEdge = pEdge->getVertices();
    alignAnchorRegion(&sequences, &distances, adg, *pMatchMap, *pSequenceAccessor, nanopores, *verticesOfEdge.first,
                      *verticesOfEdge.second,
                      pId2OverlapMap->at(std::get<0>(id2TupleTuple(verticesOfEdge.first->getId()))),
                      pId2OverlapMap->at(std::get<0>(id2TupleTuple(verticesOfEdge.second->getId()))));
  }

  auto const sortedAdg = adg.sortTopologically();

  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> mappingVertex2OrderIdx;
  std::vector<lazybastard::graph::Vertex const *>                     order;
  order.reserve(sortedAdg.size());

  for (std::size_t idx = 0; idx < sortedAdg.size(); ++idx) {
    mappingVertex2OrderIdx[sortedAdg.at(idx)] = idx;
    order.push_back(sortedAdg.at(idx));
  }

  std::unordered_map<lazybastard::graph::Vertex const *, bool>                 visitedVertices;
  std::unordered_map<lazybastard::graph::Vertex const *, std::tuple<int, int>> tap;

  std::optional<std::string> globalSequence;
  int                        globalPos1 = 0;
  int                        globalPos2 = 0;
  std::tie(globalSequence, globalPos1, globalPos2) =
      visitOrdered(&visitedVertices, &tap, adg, mappingVertex2OrderIdx, order, distances, sequences, anchorSequences,
                   *pId2OverlapMap, *order.front());

  auto const adgVertices = adg.getVertices();
  if (adgVertices.size() == 1) {
    auto const *const pAnchor = adgVertices.front();
    auto const        overlap = pId2OverlapMap->at(std::get<0>(id2TupleTuple(pAnchor->getId())));

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
        visitOrdered(&visitedVertices, &localTap, adg, mappingVertex2OrderIdx, order, distances, sequences,
                     anchorSequences, *pId2OverlapMap, *pVertex);

    if (localTap.empty()) {
      auto const overlap = pId2OverlapMap->at(std::get<0>(id2TupleTuple(pVertex->getId())));

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
    std::string targetName = ">Prokrastinator_";
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
    std::unordered_map<std::string, std::tuple<std::tuple<std::string, std::size_t>, std::size_t>> mappingId2Anchor;
    for (auto const &info : vertexInfo.at(idx)) {
      mappingId2Anchor[std::get<0>(std::get<0>(info.second))] = info.second;
    }

    if (!pContainElements->contains(vertices.at(idx))) {
      continue;
    }

    for (auto const &containElement : pContainElements->at(vertices.at(idx))) {
      std::vector<std::tuple<std::pair<int, int>, std::string>> containInfo;
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
        auto const tapDir = pMatchMap->getVertexMatch(vertices.at(idx)->getId(), std::get<1>(info))->direction *
                            (vertices.at(idx)->getVertexDirection() == Direction::e_POS);
        auto const illuminaRef =
            tapDir ? pId2OverlapMap->at(std::get<0>(tapId)).second : pId2OverlapMap->at(std::get<0>(tapId)).first;

        auto const totalRef = std::get<1>(tap.at(adg.getVertex(tupleTuple2Id(tapId)))) + globalLeftMostPosition;

        auto const contDir       = containElement.matches.at(std::get<1>(info))->direction * direction;
        int        offset        = 0;
        auto const illuminaRange = containElement.matches.at(std::get<1>(info))->illuminaRange;
        if (!contDir) {
          offset = illuminaRange.first - illuminaRef;
          globalRanges.emplace_back(totalRef - offset - illuminaRange.first - illuminaRange.second, totalRef - offset);
        } else {
          offset = illuminaRange.second - illuminaRef;
          globalRanges.emplace_back(totalRef + offset - illuminaRange.first - illuminaRange.second, totalRef + offset);
        }
      }

      std::vector<std::tuple<std::string, int, int, const char *>> sequences2Write;
      for (std::size_t idxGlobalRange = 0; idxGlobalRange < globalRanges.size(); ++idxGlobalRange) {
        auto const &illuminaId = std::get<1>(containInfo[idxGlobalRange]);
        auto const &match      = containElement.matches.at(illuminaId);
        sequences2Write.emplace_back(getIlluminaSequence(*pSequenceAccessor, illuminaId, match->illuminaRange.first,
                                                         match->illuminaRange.second, match->direction * direction),
                                     std::get<0>(globalRanges.at(idxGlobalRange)),
                                     std::get<1>(globalRanges[idxGlobalRange]), "Illumina_Match");

        if (idxGlobalRange == 0) {
          continue;
        }

        auto const preNanopore = containElement.matches.at(std::get<1>(containInfo[idxGlobalRange - 1]))->nanoporeRange;
        sequences2Write.emplace_back(getNanoporeSequence(*pSequenceAccessor, containElement.nano,
                                                         preNanopore.second + 1, match->nanoporeRange.first - 1,
                                                         direction),
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
