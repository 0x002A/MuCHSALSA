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

#include <algorithm>
#include <any>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <gsl/pointers>
#include <gsl/span>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <ms/BlastFileAccessor.h>
#include <ms/BlastFileReader.h>
#include <ms/Debug.h>
#include <ms/Kernel.h>
#include <ms/OutputWriter.h>
#include <ms/Registry.h>
#include <ms/SequenceAccessor.h>
#include <ms/Util.h>
#include <ms/graph/Edge.h>
#include <ms/graph/Graph.h>
#include <ms/graph/Vertex.h>
#include <ms/matching/Id2OverlapMap.h>
#include <ms/matching/MatchMap.h>
#include <ms/threading/Job.h>
#include <ms/threading/ThreadPool.h>
#include <ms/threading/WaitGroup.h>
#include <ms/types/Direction.h>
#include <ms/types/Toggle.h>

#include "Application.h"

// =====================================================================================================================
//                                                  USING DECLARATIONS
// =====================================================================================================================

using muchsalsa::BlastFileAccessor;
using muchsalsa::BlastFileReader;
using muchsalsa::Direction;
using muchsalsa::GraphUtil;
using muchsalsa::OutputWriter;
using muchsalsa::Registry;
using muchsalsa::SequenceAccessor;

using muchsalsa::graph::Edge;
using muchsalsa::graph::EdgeOrder;
using muchsalsa::graph::Graph;
using muchsalsa::graph::Vertex;

using muchsalsa::matching::ContainElement;
using muchsalsa::matching::Id2OverlapMap;
using muchsalsa::matching::MatchMap;
using muchsalsa::matching::VertexMatch;

using muchsalsa::threading::Job;
using muchsalsa::threading::ThreadPool;
using muchsalsa::threading::WaitGroup;

using muchsalsa::util::make_not_null_and_const;

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr static auto BASE_WEIGHT_MULTIPLICATOR = 1.1;
constexpr static auto MAX_WEIGHT_MULTIPLICATOR  = 0.8;

// =====================================================================================================================
//                                                   JOB DECLARATIONS
// =====================================================================================================================

void chainingAndOverlaps(gsl::not_null<Job const *> pJob);

void findContractionEdges(gsl::not_null<Job const *> pJob);

void findContractionTargets(gsl::not_null<Job const *> pJob);

void findDeletableVertices(gsl::not_null<Job const *> pJob);

void contract(gsl::not_null<Job const *> pJob);

void findDeletableEdges(gsl::not_null<Job const *> pJob);

void computeBitweight(gsl::not_null<Job const *> pJob);

void decycle(gsl::not_null<Job const *> pJob);

void assemblePaths(gsl::not_null<Job const *> pJob);

// =====================================================================================================================
//                                                         Main
// =====================================================================================================================

auto main(int const argc, char const *argv[]) -> int {
  try {
    gsl::span<char const *> const args = {argv, static_cast<std::size_t>(argc)};
    Application                   app(args);

    if (!app.checkIntegrity()) {
      std::cerr << "Paths are pointing to invalid/unusable locations" << std::endl;

      return -1;
    }

    // Initialize threading
    auto const threadCount = app.getThreadCount();
    auto       threadPool  = ThreadPool(threadCount);

    // Initialize structures
    auto graph            = Graph();
    auto matchMap         = MatchMap(&threadPool, &graph);
    auto registryNanopore = Registry();
    auto registryIllumina = Registry();

    // Read BLAST file
    auto bfAccessor = BlastFileAccessor(app.getContigsFilePath());
    auto blastFileReader =
        BlastFileReader(&threadPool, &bfAccessor, &graph, &matchMap, &registryNanopore, &registryIllumina);
    blastFileReader.read();
    matchMap.calculateEdges();

    TRACE("Order: %lu, Size: %lu\n", graph.getOrder(), graph.getSize());

    SequenceAccessor sequenceAccessor(&threadPool, app.getNanoporeFilePath(), app.getUnitigsFilePath(),
                                      &registryNanopore, &registryIllumina);
    sequenceAccessor.buildIndex();

    registryNanopore.clear();
    registryIllumina.clear();

    WaitGroup wg;

    auto chainingJob = [](Job const *const pJob) { chainingAndOverlaps(pJob); };
    auto edges       = graph.getEdges();
    std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
      wg.add(1);

      auto job = Job(chainingJob, &wg, &matchMap, pEdge, app.getWiggleRoom());
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    std::mutex                                          mutex;
    std::unordered_map<Edge const *, EdgeOrder const *> contractionEdges;
    auto contractionEdgesJob = [](Job const *const pJob) { findContractionEdges(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) {
      wg.add(1);

      auto job = Job(contractionEdgesJob, &wg, &graph, pEdge, &contractionEdges, std::ref(mutex), app.getWiggleRoom());
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    TRACE("Number of contraction edges: %lu\n", contractionEdges.size());

    std::unordered_map<Vertex const *, Vertex const *> contractionTargets;
    for (auto const *const pVertex : graph.getVertices()) {
      contractionTargets.insert({pVertex, pVertex});
    }

    auto contractionTargetsJob = [](Job const *const pJob) { findContractionTargets(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = Job(contractionTargetsJob, &wg, pOrder, &contractionTargets, std::ref(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    TRACE("Number of contraction target: %lu\n", contractionTargets.size());

    std::set<Vertex const *const> deletableVertices;
    std::set<Vertex const *const> contractionRoots;
    auto                          deletableVerticesJob = [](Job const *const pJob) { findDeletableVertices(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = Job(deletableVerticesJob, &wg, pOrder, &deletableVertices, &contractionRoots, &contractionTargets,
                     std::ref(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    TRACE("Number of contraction roots: %lu\n", contractionRoots.size());
    TRACE("Vertices to become deleted: %lu\n", deletableVertices.size());

    std::unordered_map<Vertex const *, std::vector<ContainElement>> containElements;
    auto contractionJob = [](Job const *const pJob) { contract(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      if (!contractionRoots.contains(pOrder->endVertex)) {
        continue;
      }

      wg.add(1);
      auto job = Job(contractionJob, &wg, pOrder, &containElements, &matchMap, std::ref(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::for_each(std::begin(deletableVertices), std::end(deletableVertices),
                  [&](auto const *const pTarget) { graph.deleteVertex(pTarget); });
    deletableVertices.clear();
    edges = graph.getEdges();

    std::set<Edge const *const> deletableEdges;
    auto                        deletableEdgesJob = [](Job const *const pJob) { findDeletableEdges(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
      wg.add(1);
      auto job = Job(deletableEdgesJob, &wg, pEdge, &deletableEdges, std::ref(mutex));
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    TRACE("Edges to become deleted: %lu\n", deletableEdges.size());

    std::for_each(std::begin(deletableEdges), std::end(deletableEdges),
                  [&](auto const *const pTarget) { graph.deleteEdge(pTarget); });
    deletableEdges.clear();

    TRACE("Order: %lu, Size: %lu\n", graph.getOrder(), graph.getSize());

    edges                    = graph.getEdges();
    auto computeBitweightJob = [](Job const *const pJob) { computeBitweight(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
      wg.add(1);
      auto job = Job(computeBitweightJob, &wg, pEdge);
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    auto const maxSpanTree = muchsalsa::getMaxSpanTree(graph);

    auto decycleJob = [](Job const *const pJob) { decycle(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) {
      wg.add(1);
      auto job = Job(decycleJob, &wg, &graph, pEdge, &maxSpanTree, &deletableEdges, std::ref<std::mutex>(mutex));
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    TRACE("Edges to become deleted: %lu\n", deletableEdges.size());

    std::for_each(std::begin(deletableEdges), std::end(deletableEdges),
                  [&](auto const *const pTarget) { graph.deleteEdge(pTarget); });
    deletableEdges.clear();
    edges = graph.getEdges();

    TRACE("Order: %lu, Size: %lu\n", graph.getOrder(), graph.getSize());
    TRACE("Starting assembly\n");

    //// OUTPUT FILES ////
    OutputWriter outputWriter(app.getOutputFilePath() + "/temp.query.fa", app.getOutputFilePath() + "/temp.align.paf",
                              app.getOutputFilePath() + "/temp.target.fa");
    //// /OUTPUT FILES ////

    std::atomic<int> assemblyIdx{-1};
    auto             connectedComponents = muchsalsa::getConnectedComponents(graph);
    auto             assemblyJob         = [](Job const *const pJob) { assemblePaths(pJob); };
    for (auto const &connectedComponent : connectedComponents) {

      wg.add(1);
      auto job = Job(assemblyJob, &wg, &graph, &matchMap, &containElements, &sequenceAccessor, &connectedComponent,
                     &assemblyIdx, std::ref(outputWriter));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::cout << "Finished assembly\n";
  } catch (std::exception const &e) {
    std::cerr << "Exception occurred: " << '\n' << e.what();
  }

  return 0;
}

// =====================================================================================================================
//                                                    JOB DEFINITIONS
// =====================================================================================================================

void chainingAndOverlaps(gsl::not_null<Job const *> const pJob) {
  std::vector<unsigned int> plusIDs;
  std::vector<unsigned int> minusIDs;

  auto const        pMatchMap    = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(1)));
  auto const        pEdge        = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(2)));
  auto const *const pEdgeMatches = pMatchMap->getEdgeMatches(pEdge);

  if (!pEdgeMatches) {
    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
    return;
  }

  for (auto const &[illuminaId, edgeMatch] : *pEdgeMatches) {
    LB_UNUSED(edgeMatch);

    if (edgeMatch->direction) {
      plusIDs.push_back(illuminaId);
    } else {
      minusIDs.push_back(illuminaId);
    }
  }

  auto const wiggleRoom = std::any_cast<std::size_t>(pJob->getParam(3));
  auto       minusPaths = muchsalsa::getMaxPairwisePaths(*pMatchMap, *pEdge, minusIDs, false, wiggleRoom);
  auto       plusPaths  = muchsalsa::getMaxPairwisePaths(*pMatchMap, *pEdge, plusIDs, true, wiggleRoom);

  auto hasPrimary =
      std::any_of(plusPaths.begin(), plusPaths.end(), [](auto const &plusPath) { return std::get<2>(plusPath); });

  if (!hasPrimary) {
    hasPrimary =
        std::any_of(minusPaths.begin(), minusPaths.end(), [](auto const &minusPath) { return std::get<2>(minusPath); });
  }

  if (hasPrimary) {
    plusPaths.erase(
        std::remove_if(plusPaths.begin(), plusPaths.end(), [](auto const &plusPath) { return !std::get<2>(plusPath); }),
        plusPaths.end());
    minusPaths.erase(std::remove_if(minusPaths.begin(), minusPaths.end(),
                                    [](auto const &minusPath) { return !std::get<2>(minusPath); }),
                     minusPaths.end());
  }

  auto hasMulti = std::any_of(plusPaths.begin(), plusPaths.end(),
                              [](auto const &plusPath) { return std::get<0>(plusPath).size() > 1; });

  if (!hasMulti) {
    hasMulti = std::any_of(minusPaths.begin(), minusPaths.end(),
                           [](auto const &minusPath) { return std::get<0>(minusPath).size() > 1; });
  }

  if (hasMulti) {
    plusPaths.erase(std::remove_if(plusPaths.begin(), plusPaths.end(),
                                   [](auto const &plusPath) { return std::get<0>(plusPath).size() <= 1; }),
                    plusPaths.end());
    minusPaths.erase(std::remove_if(minusPaths.begin(), minusPaths.end(),
                                    [](auto const &minusPath) { return std::get<0>(minusPath).size() <= 1; }),
                     minusPaths.end());
  }

  auto const combinedSize = plusPaths.size() + minusPaths.size();
  if (combinedSize > 1) {
    pEdge->setShadow(true);
  } else {
    auto const &path = !minusPaths.empty() ? minusPaths[0] : plusPaths[0];
    pEdge->setShadow(!std::get<2>(path));
  }

  std::for_each(std::begin(minusPaths), std::end(minusPaths), [=](auto &minusPath) {
    auto const overlap = muchsalsa::computeOverlap(*pMatchMap, std::get<0>(minusPath), *pEdge, false,
                                                   std::get<1>(minusPath), std::get<2>(minusPath));
    if (overlap.has_value()) {
      pEdge->appendOrder(std::move(overlap.value()));
    }
  });

  std::for_each(std::begin(plusPaths), std::end(plusPaths), [=](auto &plusPath) {
    auto const order = muchsalsa::computeOverlap(*pMatchMap, std::get<0>(plusPath), *pEdge, true, std::get<1>(plusPath),
                                                 std::get<2>(plusPath));
    if (order.has_value()) {
      pEdge->appendOrder(std::move(order.value()));
    }
  });

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findContractionEdges(gsl::not_null<Job const *> const pJob) {
  auto const pGraph = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(1)));
  auto const pEdge  = make_not_null_and_const(std::any_cast<Edge const *const>(pJob->getParam(2)));

  auto const &edgeOrders = pEdge->getEdgeOrders();
  for (auto const &order : edgeOrders) {
    if (order.isContained && order.isPrimary) {
      auto       isSane    = true;
      auto const neighbors = pGraph->getNeighbors(order.startVertex);
      auto const edges     = std::map<unsigned int, Edge *>(std::begin(neighbors), std::end(neighbors));
      for (auto const &[targetId, pSubedge] : edges) {
        LB_UNUSED(pSubedge);

        auto const *const pTargetVertex = pGraph->getVertex(targetId);

        if (targetId == order.endVertex->getId() || pSubedge->isShadow()) {
          continue;
        }

        isSane &= pGraph->hasEdge(std::make_pair(order.endVertex, pTargetVertex));
        if (!isSane) {
          break;
        }

        auto const wiggleRoom = std::any_cast<std::size_t>(pJob->getParam(5));
        isSane &=
            muchsalsa::sanityCheck(*pGraph, *order.startVertex, *order.endVertex, *pTargetVertex, order, wiggleRoom);

        if (!isSane) {
          break;
        }
      }

      if (isSane) {
        std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(4)).get());
        auto const                  pContractionEdges =
            gsl::make_not_null(std::any_cast<std::unordered_map<Edge const *, EdgeOrder const *> *>(pJob->getParam(3)));
        pContractionEdges->insert({pEdge, &order});

        break;
      }
    }
  }
  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findContractionTargets(gsl::not_null<Job const *> const pJob) {
  auto const pOrder = gsl::make_not_null(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
  auto const pContractionTargets =
      gsl::make_not_null(std::any_cast<std::unordered_map<Vertex const *, Vertex const *> *>(pJob->getParam(2)));
  {
    std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get());

    auto const *const contractTo = (*pContractionTargets)[pOrder->endVertex];

    if ((*pContractionTargets)[pOrder->startVertex] == pOrder->startVertex ||
        (*pContractionTargets)[pOrder->startVertex]->getMetaDatum<std::size_t>(0) >
            contractTo->getMetaDatum<std::size_t>(0)) {
      (*pContractionTargets)[pOrder->startVertex] = contractTo;
    }
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findDeletableVertices(gsl::not_null<const Job *> pJob) {
  auto const pOrder             = gsl::make_not_null(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
  auto const pDeletableVertices = gsl::make_not_null(std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(2)));
  {
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
    pDeletableVertices->insert(pOrder->startVertex);
  }

  auto const pContractionTargets =
      gsl::make_not_null(std::any_cast<std::unordered_map<Vertex const *, Vertex const *> *>(pJob->getParam(4)));
  auto const pContractionRoots = gsl::make_not_null(std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(3)));
  {
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT

    auto const *const contractTo = (*pContractionTargets)[pOrder->startVertex];

    pContractionRoots->insert(contractTo);
    pContractionRoots->erase(pOrder->startVertex);
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void contract(gsl::not_null<Job const *> const pJob) {
  auto const pOrder    = make_not_null_and_const(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
  auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(3)));

  std::unordered_map<unsigned int, VertexMatch const *> matches;
  std::for_each(std::begin(pOrder->ids), std::end(pOrder->ids), [&](auto const id) {
    auto const *pVertexMatches = pMatchMap->getVertexMatch(pOrder->startVertex->getId(), id);
    if (pVertexMatches != nullptr) {
      matches.insert({id, pVertexMatches});
    }
  });

  std::lock_guard<std::mutex> guard(
      std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(4)).get()); // NOLINT
  auto const pContainElements = gsl::make_not_null(
      std::any_cast<std::unordered_map<Vertex const *, std::vector<ContainElement>> *>(pJob->getParam(2)));

  (*pContainElements)[pOrder->endVertex].push_back(ContainElement{std::move(matches), pOrder->startVertex->getId(),
                                                                  pOrder->startVertex->getNanoporeLength(),
                                                                  pOrder->score, pOrder->direction, pOrder->isPrimary});

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findDeletableEdges(gsl::not_null<Job const *> const pJob) {
  auto const             pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(1)));
  std::vector<EdgeOrder> filteredOrders;
  std::remove_copy_if(std::begin(pEdge->getEdgeOrders()), std::end(pEdge->getEdgeOrders()),
                      std::back_inserter(filteredOrders), [](auto const &order) { return order.isContained; });

  if (filteredOrders.empty()) {
    auto const pDeletableEdges = gsl::make_not_null(std::any_cast<std::set<Edge const *const> *>(pJob->getParam(2)));
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get()); // NOLINT
    pDeletableEdges->insert(pEdge);
  }

  pEdge->replaceOrders(std::move(filteredOrders));
  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void computeBitweight(gsl::not_null<Job const *> const pJob) {
  auto const  pEdge  = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(1)));
  auto const &orders = pEdge->getEdgeOrders();

  if (orders.empty()) {
    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
    return;
  }

  if (pEdge->isShadow()) {
    auto const initial   = orders[0].direction;
    auto const findOther = std::find_if(orders.begin(), orders.end(),
                                        [&](auto const &order) { return order.direction != initial; }) != orders.end();
    if (!findOther) {
      pEdge->setConsensusDirection(orders[0].direction);
    }
  } else {
    pEdge->setWeight(orders[0].score);
    pEdge->setConsensusDirection(orders[0].direction);
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void decycle(gsl::not_null<Job const *> const pJob) {
  auto const pEdge        = gsl::make_not_null(std::any_cast<Edge const *>(pJob->getParam(2)));
  auto const pMaxSpanTree = make_not_null_and_const(std::any_cast<Graph const *>(pJob->getParam(3)));
  if (pEdge->getConsensusDirection() != Direction::e_NONE &&
      !pMaxSpanTree->hasEdge(std::make_pair(pEdge->getVertices().first, pEdge->getVertices().second))) {
    auto const        shortestPath = GraphUtil::getShortestPath(pMaxSpanTree, pEdge->getVertices());
    muchsalsa::Toggle direction    = pEdge->getConsensusDirection() == Direction::e_POS;
    auto const        baseWeight   = static_cast<double>(pEdge->getWeight());

    std::vector<double> weights;
    auto const          pGraph = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(1)));
    for (auto it = shortestPath.begin(); it != std::prev(shortestPath.end()); ++it) {
      auto const pPathEdge = make_not_null_and_const(pGraph->getEdge(std::make_pair(*it, *std::next(it))));

      direction *= pPathEdge->getConsensusDirection() == Direction::e_POS;
      weights.push_back(static_cast<double>(pPathEdge->getWeight()));
    }

    if (!direction && !weights.empty()) {
      auto const iterMinWeight = std::min_element(weights.begin(), weights.end());
      auto const iterMaxWeight = std::max_element(weights.begin(), weights.end());

      auto const pDeletableEdges = gsl::make_not_null(std::any_cast<std::set<Edge const *const> *>(pJob->getParam(4)));
      if (*iterMinWeight < baseWeight || (baseWeight * BASE_WEIGHT_MULTIPLICATOR >= *iterMinWeight &&
                                          *iterMinWeight < *iterMaxWeight * MAX_WEIGHT_MULTIPLICATOR)) {
        auto const        minWeightIdx   = std::distance(weights.begin(), iterMinWeight);
        auto const *const pDeletableEdge = pGraph->getEdge(std::make_pair(
            *std::next(shortestPath.begin(), minWeightIdx), *std::next(shortestPath.begin(), minWeightIdx + 1)));

        std::lock_guard<std::mutex> lck(
            std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
        pDeletableEdges->insert(pDeletableEdge);
      }

      {
        std::lock_guard<std::mutex> lck(
            std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
        pDeletableEdges->insert(pEdge);
      }
    }
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void assemblePaths(gsl::not_null<Job const *> const pJob) {
  auto const pGraph = gsl::make_not_null(std::any_cast<Graph *>(pJob->getParam(1)));
  auto const pConnectedComponent =
      make_not_null_and_const(std::any_cast<std::vector<muchsalsa::graph::Vertex *> const *>(pJob->getParam(5)));

  auto const subGraph         = pGraph->getSubgraph(*pConnectedComponent);
  auto const subGraphVertices = subGraph.getVertices();
  auto const maxNplVertexIter = std::max_element(
      std::begin(subGraphVertices), std::end(subGraphVertices),
      [](Vertex const *v1, Vertex const *v2) { return v1->getNanoporeLength() < v2->getNanoporeLength(); });
  auto const *const pMaxNplVertex = maxNplVertexIter == std::end(subGraph.getVertices()) ? nullptr : *maxNplVertexIter;
  auto const        pMatchMap     = gsl::make_not_null(std::any_cast<MatchMap *>(pJob->getParam(2)));

  if (pMaxNplVertex != nullptr) {
    auto       diGraph = muchsalsa::getDirectionGraph(pMatchMap, *pGraph, subGraph, *pMaxNplVertex);
    auto const paths   = muchsalsa::linearizeGraph(&diGraph);

    Id2OverlapMap id2OverlapMap;
    auto *const   pAssemblyIdx = std::any_cast<std::atomic<int> *>(pJob->getParam(6));
    std::for_each(std::begin(paths), std::end(paths), [&](auto const &path) {
      muchsalsa::assemblePath(
          &id2OverlapMap, *pMatchMap,
          *std::any_cast<std::unordered_map<Vertex const *, std::vector<ContainElement>> *>(pJob->getParam(3)),
          *std::any_cast<SequenceAccessor *>(pJob->getParam(4)), path, diGraph, ++(*pAssemblyIdx),
          std::any_cast<std::reference_wrapper<OutputWriter>>(pJob->getParam(7)).get());
    });
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
