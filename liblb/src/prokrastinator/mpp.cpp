#include "Prokrastinator.h"

#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>

#include "Util.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "matching/MatchMap.h"

using path_t = std::tuple<std::deque<gsl::not_null<std::string const *> const>, std::size_t, bool>;
using match_t = std::tuple<std::pair<int, int>, gsl::not_null<std::string const *>>;

namespace {
bool checkCompatibility(gsl::not_null<lazybastard::matching::MatchMap const *> const pMatches,
                        gsl::not_null<lazybastard::graph::Edge const *> const pEdge, std::string const &illumindaID1,
                        std::string const &illumindaID2, std::size_t wiggleRoom) {
  auto const nanoCheck = [&](int *pOrientation, float *pDiff,
                             gsl::not_null<lazybastard::graph::Vertex const *> const pVertex) -> bool {
    auto const edgeMatch1 = gsl::make_not_null(pMatches->getEdgeMatch(pEdge->getID(), illumindaID1));
    auto const edgeMatch2 = gsl::make_not_null(pMatches->getEdgeMatch(pEdge->getID(), illumindaID2));

    auto const vertexMatch1 = gsl::make_not_null(pMatches->getVertexMatch(pVertex->getID(), illumindaID1));
    auto const vertexMatch2 = gsl::make_not_null(pMatches->getVertexMatch(pVertex->getID(), illumindaID2));

    auto ncl1 = static_cast<float>(edgeMatch1->overlap.first - vertexMatch1->illuminaRange.first) //
                / vertexMatch1->rRatio;
    auto ncr1 = static_cast<float>(vertexMatch1->illuminaRange.second - edgeMatch1->overlap.second) //
                / vertexMatch1->rRatio;

    lazybastard::util::swap_if(ncl1, ncr1, !vertexMatch1->direction);

    auto ncl2 = static_cast<float>(edgeMatch2->overlap.first - vertexMatch2->illuminaRange.first) //
                / vertexMatch2->rRatio;
    auto ncr2 = static_cast<float>(vertexMatch2->illuminaRange.second - edgeMatch2->overlap.second) //
                / vertexMatch2->rRatio;

    lazybastard::util::swap_if(ncl2, ncr2, !vertexMatch2->direction);

    auto const correctedNano1 = std::make_pair(static_cast<float>(vertexMatch1->nanoporeRange.first) + ncl1,
                                               static_cast<float>(vertexMatch1->nanoporeRange.second) + ncr1);
    auto const correctedNano2 = std::make_pair(static_cast<float>(vertexMatch2->nanoporeRange.first) + ncl2,
                                               static_cast<float>(vertexMatch2->nanoporeRange.second) + ncr2);

    *pOrientation = 0, *pDiff = 0;
    if (correctedNano1.first <= correctedNano2.second && correctedNano2.first <= correctedNano1.second) {

      // First Case
      lazybastard::util::exchange_if(*pOrientation, 2,
                                     correctedNano1.first < correctedNano2.first &&
                                         correctedNano1.second < correctedNano2.second);
      lazybastard::util::exchange_if(*pDiff, correctedNano1.second - correctedNano2.first + 1.0F,
                                     correctedNano1.first < correctedNano2.first &&
                                         correctedNano1.second < correctedNano2.second);

      // Second Case
      lazybastard::util::exchange_if(*pOrientation, -2,
                                     correctedNano1.first > correctedNano2.first &&
                                         correctedNano1.second > correctedNano2.second);
      lazybastard::util::exchange_if(*pDiff, correctedNano2.second - correctedNano1.first + 1.0F,
                                     correctedNano1.first > correctedNano2.first &&
                                         correctedNano1.second > correctedNano2.second);
    } else if (correctedNano1.first < correctedNano2.first) {
      *pOrientation = 1;
      *pDiff = correctedNano2.first - correctedNano1.second + 1.0F;
    } else {
      *pOrientation = -1;
      *pDiff = correctedNano1.first - correctedNano2.second + 1.0F;
    }

    auto uco = 0;
    if (vertexMatch1->nanoporeRange.first <= vertexMatch2->nanoporeRange.second &&
        vertexMatch2->nanoporeRange.first <= vertexMatch1->nanoporeRange.second) {

      // First Case
      lazybastard::util::exchange_if(uco, 2,
                                     (vertexMatch1->nanoporeRange.first < vertexMatch2->nanoporeRange.first) &&
                                         (vertexMatch1->nanoporeRange.second < vertexMatch2->nanoporeRange.second));

      // Second Case
      lazybastard::util::exchange_if(uco, -2,
                                     vertexMatch1->nanoporeRange.first > vertexMatch2->nanoporeRange.first &&
                                         vertexMatch1->nanoporeRange.second > vertexMatch2->nanoporeRange.second);
    }

    auto abort = (*pOrientation < 0 && uco >= 0);
    abort |= (*pOrientation > 0 && uco <= 0);

    return abort;
  };

  auto const vertices = pEdge->getVertices();

  auto orientation1 = 0;
  auto orientation2 = 0;
  auto diff1 = 0.0F;
  auto diff2 = 0.0F;

  auto abort = false;
  abort &= nanoCheck(&orientation1, &diff1, vertices.first);
  abort &= nanoCheck(&orientation2, &diff2, vertices.second);

  if (abort) {
    return false;
  }

  auto pEdgeMatch1 = gsl::make_not_null(pMatches->getEdgeMatch(pEdge->getID(), illumindaID1));
  if (!pEdgeMatch1->direction) {
    orientation2 *= -1;
  }

  auto matching = false;
  if (orientation1 == orientation2 && orientation1 != 0) {
    auto const diff = std::max(diff1, diff2) - std::min(diff1, diff2);
    matching = (diff <= static_cast<float>(wiggleRoom)) || (diff * 100 / std::max(diff1, diff2) <= 15);
  } else if ((orientation1 < 0 and orientation2 < 0) || (orientation1 > 0 && orientation2 > 0)) {
    matching = diff1 + diff2 <= static_cast<float>(wiggleRoom);
  }

  return matching;
}
} // unnamed namespace

std::vector<path_t> lazybastard::getMaxPairwisePaths(
    gsl::not_null<lazybastard::matching::MatchMap const *> const pMatches,
    gsl::not_null<lazybastard::graph::Edge const *> const pEdge,
    std::set<gsl::not_null<std::string const *> const,
             lazybastard::util::LTCmp<gsl::not_null<std::string const *> const>> const &illuminaIDs,
    bool direction, std::size_t wiggleRoom) {
  std::vector<path_t> result;

  if (illuminaIDs.empty()) {
    return result;
  }

  std::vector<match_t> vStart;
  std::vector<match_t> vEnd;

  auto const idCount = illuminaIDs.size();
  vStart.reserve(idCount);
  vEnd.reserve(idCount);

  auto const vertexIDs = pEdge->getVertices();
  auto const storeIlluminaID = [&](auto const pIlluminaID) {
    auto const *const start = pMatches->getVertexMatch(vertexIDs.first->getID(), *pIlluminaID);
    auto const *const end = pMatches->getVertexMatch(vertexIDs.second->getID(), *pIlluminaID);

    vStart.emplace_back(start->nanoporeRange, pIlluminaID);
    vEnd.emplace_back(end->nanoporeRange, pIlluminaID);
  };
  std::for_each(std::begin(illuminaIDs), std::end(illuminaIDs), storeIlluminaID);

  std::sort(std::begin(vStart), std::end(vStart));
  std::sort(std::begin(vEnd), std::end(vEnd));

  util::reverse_if(std::begin(vEnd), std::end(vEnd), !direction);

  auto const length = vStart.size();
  std::vector<std::tuple<std::vector<std::size_t>, float>> population;
  population.reserve(length);

  std::transform(std::begin(vStart), std::end(vStart), std::back_inserter(population), [&](auto const &match) {
    return std::make_tuple(std::vector<std::size_t>(),
                           static_cast<float>(pMatches->getEdgeMatch(pEdge->getID(), *std::get<1>(match))->score));
  });

  auto const limit = std::max(length - 1, 0UL);
  for (std::size_t k = 0; k < limit; ++k) {
    for (std::size_t l = k + 1; l < limit; ++l) {
      auto isCompatible =
          checkCompatibility(pMatches, pEdge, *std::get<1>(vStart[k]), *std::get<1>(vStart[l]), wiggleRoom);
      auto const score =
          std::get<1>(population[k]) + pMatches->getEdgeMatch(pEdge->getID(), *std::get<1>(vStart[l]))->score;
      isCompatible &= score > std::get<1>(population[l]);

      if (isCompatible) {
        auto vTmp = std::get<0>(population[k]);
        vTmp.push_back(k);

        population[l] = std::make_tuple(std::move(vTmp), score);
      }
    }
  }

  float maxPathVal = 0.0F;
  auto maxPathIter = std::begin(population);
  for (auto iter = std::begin(population); iter != std::end(population); ++iter) {
    std::get<0>(*iter).push_back(static_cast<std::size_t>(iter - std::begin(population)));

    auto const currentPathVal = std::get<1>(*iter);
    util::exchange_if(maxPathVal, currentPathVal, currentPathVal > maxPathVal);
    util::exchange_if(maxPathIter, iter, currentPathVal > maxPathVal);
  }

  auto const &vMaxPath = std::get<0>(*maxPathIter);
  std::deque<gsl::not_null<std::string const *> const> dTmp;
  std::transform(std::begin(vMaxPath), std::end(vMaxPath), std::back_inserter(dTmp),
                 [&](auto const &idx) { return std::get<1>(vStart[idx]); });

  auto hasPrimary = std::any_of(std::begin(vMaxPath), std::end(vMaxPath), [&](auto idx) {
    return pMatches->getEdgeMatch(pEdge->getID(), *std::get<1>(vStart[idx]))->isPrimary;
  });
  hasPrimary |= vMaxPath.size() > 2;
  result.emplace_back(std::move(dTmp), maxPathVal, hasPrimary);

  for (auto const &populationEntry : population) {

    auto const score = std::get<1>(populationEntry);
    auto const pathThreshold = maxPathVal * 0.75F;
    if (score > pathThreshold) {
      auto disjoint = true;

      auto const &vIdx = std::get<0>(populationEntry);
      auto const isDisjoint = [&](auto const &viable) {
        disjoint &= !std::any_of(std::begin(vIdx), std::end(vIdx), [&](auto const &idx) {
          auto const &vSearch = std::get<0>(viable);
          auto const &illuminaID = std::get<1>(vStart[idx]);
          return std::find(std::begin(vSearch), std::end(vSearch), illuminaID) != std::end(vSearch);
        });
      };
      std::for_each(std::begin(result), std::end(result), isDisjoint);

      if (disjoint) {
        dTmp.clear(); // Return object to known state
        std::transform(std::begin(vIdx), std::end(vIdx), std::back_inserter(dTmp),
                       [&](auto const &idx) { return std::get<1>(vStart[idx]); });

        result.emplace_back(std::move(dTmp), score, std::any_of(std::begin(vIdx), std::end(vIdx), [&](auto idx) {
                              return pMatches->getEdgeMatch(pEdge->getID(), *std::get<1>(vStart[idx]))->isPrimary;
                            }));
      }
    }
  }

  auto alreadyAShadow = result.size() == 1;
  alreadyAShadow &= std::get<2>(result.front());
  if (alreadyAShadow) {
    using match_id_t = std::tuple<std::pair<int, int>, std::string const *>;
    std::vector<match_id_t> vIDsStart;
    std::vector<match_id_t> vIDsEnd;

    auto startVertexMatches = gsl::make_not_null(pMatches->getVertexMatches(vertexIDs.first->getID()));
    auto endVertexMatches = gsl::make_not_null(pMatches->getVertexMatches(vertexIDs.first->getID()));

    vIDsStart.reserve(startVertexMatches->size());
    vIDsEnd.reserve(endVertexMatches->size());

    std::transform(std::begin(*startVertexMatches), std::end(*startVertexMatches), std::back_inserter(vIDsStart),
                   [](auto const &match) { return std::make_tuple(match.second->nanoporeRange, &match.first); });
    std::transform(std::begin(*endVertexMatches), std::end(*endVertexMatches), std::back_inserter(vIDsEnd),
                   [](auto const &match) { return std::make_tuple(match.second->nanoporeRange, &match.first); });

    std::sort(std::begin(vIDsStart), std::end(vIDsStart));
    std::sort(std::begin(vIDsEnd), std::end(vIDsEnd));

    util::reverse_if(std::begin(vIDsEnd), std::end(vIDsEnd), !direction);

    auto const &pIDs = std::get<0>(result.front());
    if ((std::get<1>(vIDsStart.front()) != pIDs.front() && std::get<1>(vIDsEnd.front()) != pIDs.front()) ||
        (std::get<1>(vIDsStart.back()) != pIDs.back() && std::get<1>(vIDsEnd.back()) != pIDs.back())) {
      result = decltype(result)({std::make_tuple(std::get<0>(result.front()), std::get<1>(result.front()), false)});
    } else {
      auto i = 0L;
      auto j = 0L;
      auto isShadow = false;

      for (auto iter = std::begin(pIDs); !isShadow && iter != std::end(pIDs); ++iter) {
        auto const pred = [&](auto const &t) { return *iter == std::get<1>(t); };

        auto const searchIterStart = std::begin(vIDsStart) + i;
        auto const resultIterStart = std::find_if(searchIterStart, std::end(vIDsStart), pred);
        i = std::distance(searchIterStart, resultIterStart);
        auto isInter = i > 0;

        auto const searchIterEnd = std::begin(vIDsEnd) + j;
        auto const resultIterEnd = std::find_if(searchIterEnd, std::end(vIDsEnd), pred);
        j = std::distance(searchIterEnd, resultIterEnd);
        isInter &= j > 0;

        isShadow = isInter;
      }

      util::exchange_if(result, {std::make_tuple(std::get<0>(result.front()), std::get<1>(result.front()), false)},
                        isShadow);
    }
  }

  return result;
}