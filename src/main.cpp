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

#include <algorithm>
#include <any>
#include <cstddef>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <gsl/pointers>
#include <gsl/span>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <lb/BlastFileReader.h>
#include <lb/OutputWriter.h>
#include <lb/Prokrastinator.h>
#include <lb/Util.h>
#include <lb/graph/Edge.h>
#include <lb/graph/Graph.h>
#include <lb/graph/Vertex.h>
#include <lb/matching/Id2OverlapMap.h>
#include <lb/matching/MatchMap.h>
#include <lb/threading/Job.h>
#include <lb/threading/ThreadPool.h>
#include <lb/threading/WaitGroup.h>
#include <lb/types/Toggle.h>

#include "Application.h"

// =====================================================================================================================
//                                                  USING DECLARATIONS
// =====================================================================================================================

using lazybastard::Direction;
using lazybastard::GraphUtil;
using lazybastard::OutputWriter;

using lazybastard::graph::Edge;
using lazybastard::graph::EdgeOrder;
using lazybastard::graph::Graph;
using lazybastard::graph::Vertex;

using lazybastard::matching::Id2OverlapMap;
using lazybastard::matching::MatchMap;
using lazybastard::matching::VertexMatch;

using lazybastard::threading::Job;
using lazybastard::threading::ThreadPool;
using lazybastard::threading::WaitGroup;

using lazybastard::util::make_not_null_and_const;

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr static auto BASE_WEIGHT_MULTIPLICATOR = 1.1;
constexpr static auto MAX_WEIGHT_MULTIPLICATOR = 0.8;

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ---------------------
// struct ContainElement
// ---------------------

struct ContainElement {
    std::unordered_map<std::string const *, VertexMatch const *> const matches;
    Vertex const *const nano;
    std::size_t const nanoporeLength;
    std::size_t const score;
    bool const direction;
    bool const isPrimary;
};

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
    gsl::span<char const *> const args = {argv, static_cast<std::size_t>(argc)};
    Application app(args);

    if (!app.checkIntegrity()) {
        std::cerr << "Paths are pointing to invalid/unusable locations" << std::endl;

        return -1;
    }

    // Keep one free for the main program
    auto const threadCount = app.getThreadCount() - 1;
    auto threadPool = ThreadPool(threadCount);

    // Read BLAST file
    auto graph = Graph();
    auto matchMap = MatchMap(&threadPool, &graph);
    if (std::ifstream inputStream{app.getContigsFilePath(), std::ios::binary | std::ios::in}) {
        auto blastReader = lazybastard::BlastFileReader(&threadPool, inputStream, &graph, &matchMap);
        blastReader.read();
    } else {
        std::cerr << "Can't open BLAST file for reading" << std::endl;

        return -1;
    }

    try {
        matchMap.calculateEdges();

        std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

        WaitGroup wg;

        auto chainingJob = [](Job const *const pJob) { chainingAndOverlaps(pJob); };
        auto edges = graph.getEdges();
        std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
            wg.add(1);

            auto job = Job(chainingJob, &wg, &matchMap, pEdge, app.getWiggleRoom());
            threadPool.addJob(std::move(job));
        });
        wg.wait();

        std::mutex mutex;
        std::unordered_map<Edge const *, EdgeOrder const *> contractionEdges;
        auto contractionEdgesJob = [](Job const *const pJob) { findContractionEdges(pJob); };
        std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) {
            wg.add(1);

            auto job = Job(contractionEdgesJob, &wg, &graph, pEdge, &contractionEdges, std::ref(mutex),
                           app.getWiggleRoom());
            threadPool.addJob(std::move(job));
        });
        wg.wait();

        std::cout << "Number of contraction edges: " << contractionEdges.size() << std::endl;

        std::unordered_map<Vertex const *, Vertex const *> contractionTargets;
        auto contractionTargetsJob = [](Job const *const pJob) { findContractionTargets(pJob); };
        for (auto const &[pEdge, pOrder] : contractionEdges) {
            LB_UNUSED(pEdge);

            wg.add(1);
            auto job = lazybastard::threading::Job(contractionTargetsJob, &wg, pOrder, &contractionTargets,
                                                   std::ref(mutex));
            threadPool.addJob(std::move(job));
        }
        wg.wait();

        std::cout << "Number of contraction target: " << contractionTargets.size() << std::endl;

        std::set<Vertex const *const> deletableVertices;
        std::set<Vertex const *const> contractionRoots;
        auto deletableVerticesJob = [](Job const *const pJob) { findDeletableVertices(pJob); };
        for (auto const &[pEdge, pOrder] : contractionEdges) {
            LB_UNUSED(pEdge);

            wg.add(1);
            auto job = lazybastard::threading::Job(deletableVerticesJob, &wg, pOrder, &deletableVertices,
                                                   &contractionRoots,
                                                   &contractionTargets, std::ref(mutex));
            threadPool.addJob(std::move(job));
        }
        wg.wait();

        std::cout << "Vertices to become deleted: " << deletableVertices.size() << std::endl;

        std::unordered_map<Vertex const *, std::unique_ptr<ContainElement const> const> containElements;
        auto contractionJob = [](Job const *const pJob) { contract(pJob); };
        for (auto const &[pEdge, pOrder] : contractionEdges) {
            LB_UNUSED(pEdge);

            wg.add(1);
            auto job = lazybastard::threading::Job(contractionJob, &wg, pOrder, &containElements, &contractionRoots,
                                                   &matchMap, std::ref(mutex));
            threadPool.addJob(std::move(job));
        }
        wg.wait();

        std::for_each(std::begin(deletableVertices), std::end(deletableVertices), [&](auto const *const pTarget) {
            matchMap.deleteVertexMatches(pTarget->getId());
            graph.deleteVertex(pTarget);
        });
        deletableVertices.clear();
        edges = graph.getEdges();

        std::set<Edge const *const> deletableEdges;
        auto deletableEdgesJob = [](Job const *const pJob) { findDeletableEdges(pJob); };
        std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
            wg.add(1);
            auto job = Job(deletableEdgesJob, &wg, pEdge, &deletableEdges, std::ref(mutex));
            threadPool.addJob(std::move(job));
        });
        wg.wait();

        std::cout << "Edges to become deleted: " << deletableEdges.size() << std::endl;

        std::for_each(std::begin(deletableEdges), std::end(deletableEdges), [&](auto const *const pTarget) {
            matchMap.deleteEdgeMatches(pTarget->getId());
            graph.deleteEdge(pTarget);
        });
        deletableEdges.clear();

        std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

        edges = graph.getEdges();
        auto computeBitweightJob = [](Job const *const pJob) { computeBitweight(pJob); };
        std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
            wg.add(1);
            auto job = Job(computeBitweightJob, &wg, pEdge);
            threadPool.addJob(std::move(job));
        });
        wg.wait();

        auto const maxSpanTree = lazybastard::getMaxSpanTree(&graph);

        auto decycleJob = [](Job const *const pJob) { decycle(pJob); };
        std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) {
            wg.add(1);
            auto job = Job(decycleJob, &wg, &graph, pEdge, maxSpanTree.get(), &deletableEdges,
                           std::ref<std::mutex>(mutex));
            threadPool.addJob(std::move(job));
        });
        wg.wait();

        std::cout << "Edges to become deleted: " << deletableEdges.size() << std::endl;

        std::for_each(std::begin(deletableEdges), std::end(deletableEdges), [&](auto const *const pTarget) {
            matchMap.deleteEdgeMatches(pTarget->getId());
            graph.deleteEdge(pTarget);
        });
        deletableEdges.clear();
        edges = graph.getEdges();

        std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl
                  << "Starting assembly" << std::endl;

        //// OUTPUT FILES ////
        std::ofstream ofQuery(app.getOutputFilePath() + "/temp.query.fa", std::ios::binary);
        std::ofstream ofPaf(app.getOutputFilePath() + "/temp.align.paf", std::ios::binary);
        std::ofstream ofTarget(app.getOutputFilePath() + "/temp.target.fa", std::ios::binary);

        OutputWriter outputWriter(ofQuery, ofPaf, ofTarget);
        //// /OUTPUT FILES ////

        std::size_t assemblyIdx = 0;
        auto connectedComponents = lazybastard::getConnectedComponents(&graph);
        auto assemblyJob = [](Job const *const pJob) { assemblePaths(pJob); };
        for (auto const &connectedComponent : connectedComponents) {

            wg.add(1);
            auto job = Job(assemblyJob, &wg, &graph, &matchMap, &connectedComponent, assemblyIdx,
                           std::ref(outputWriter));
            threadPool.addJob(std::move(job));

            ++assemblyIdx;
        }
        wg.wait();

        std::cout << "Finished assembly " << std::endl;
    } catch (std::exception const &e) {
        std::cerr << "Exception occurred: " << '\n' << e.what();
    }

    return 0;
}

// =====================================================================================================================
//                                                    JOB DEFINITIONS
// =====================================================================================================================

void chainingAndOverlaps(gsl::not_null<Job const *> const pJob) {
    std::vector<std::string> plusIDs;
    std::vector<std::string> minusIDs;

    auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(1)));
    auto const pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(2)));
    auto const *const pEdgeMatches = pMatchMap->getEdgeMatches(pEdge->getId());

    if (!pEdgeMatches) {
        std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
        return;
    }

    for (auto const &[illuminaId, edgeMatch] : *pEdgeMatches) {
        if (edgeMatch->direction) {
            plusIDs.push_back(illuminaId);
        } else {
            minusIDs.push_back(illuminaId);
        }
    }

    auto const wiggleRoom = std::any_cast<std::size_t>(pJob->getParam(3));
    auto plusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pEdge, plusIDs, true, wiggleRoom);
    auto minusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pEdge, minusIDs, false, wiggleRoom);

    auto hasPrimary =
            std::any_of(plusPaths.begin(), plusPaths.end(), [](auto const &plusPath) { return std::get<2>(plusPath); });

    if (!hasPrimary) {
        hasPrimary =
                std::any_of(minusPaths.begin(), minusPaths.end(),
                            [](auto const &minusPath) { return std::get<2>(minusPath); });
    }

    if (hasPrimary) {
        plusPaths.erase(
                std::remove_if(plusPaths.begin(), plusPaths.end(),
                               [](auto const &plusPath) { return !std::get<2>(plusPath); }),
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
        auto const overlap = lazybastard::computeOverlap(pMatchMap, std::get<0>(minusPath), pEdge, false,
                                                         std::get<1>(minusPath), std::get<2>(minusPath));
        if (overlap.has_value()) {
            pEdge->appendOrder(std::move(overlap.value()));
        }
    });

    std::for_each(std::begin(plusPaths), std::end(plusPaths), [=](auto &plusPath) {
        auto const order = lazybastard::computeOverlap(pMatchMap, std::get<0>(plusPath), pEdge, true,
                                                       std::get<1>(plusPath),
                                                       std::get<2>(plusPath));
        if (order.has_value()) {
            pEdge->appendOrder(std::move(order.value()));
        }
    });

    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findContractionEdges(gsl::not_null<Job const *> const pJob) {
    auto const pGraph = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(1)));
    auto const pEdge = make_not_null_and_const(std::any_cast<Edge const *const>(pJob->getParam(2)));

    auto const &edgeOrders = pEdge->getEdgeOrders();
    for (auto const &order : edgeOrders) {
        if (order.isContained && order.isPrimary) {
            auto isSane = true;
            auto const neighbors = pGraph->getNeighbors(order.startVertex);
            auto const edges = std::map<std::string, Edge *>(std::begin(neighbors), std::end(neighbors));
            for (auto const &[targetId, pSubedge] : edges) {
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
                        lazybastard::sanityCheck(pGraph, order.startVertex, order.endVertex, pTargetVertex, &order,
                                                 wiggleRoom);

                if (!isSane) {
                    break;
                }
            }

            if (isSane) {
                std::lock_guard<std::mutex> guard(
                        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(4)).get());
                auto const pContractionEdges =
                        gsl::make_not_null(std::any_cast<std::unordered_map<Edge const *, EdgeOrder const *> *>(
                                pJob->getParam(3)));
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

        auto const iter = pContractionTargets->find(pOrder->endVertex);
        auto const *const contractTo = iter == pContractionTargets->end() ? pOrder->endVertex : iter->second;

        if (pContractionTargets->find(pOrder->startVertex) == pContractionTargets->end()) {
            pContractionTargets->insert({pOrder->startVertex, contractTo});
        }
    }

    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findDeletableVertices(gsl::not_null<const Job *> pJob) {
    auto const pOrder = gsl::make_not_null(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
    auto const pDeletableVertices = gsl::make_not_null(
            std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(2)));
    {
        std::lock_guard<std::mutex> guard(
                std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
        pDeletableVertices->insert(pOrder->startVertex);
    }

    auto const pContractionTargets =
            gsl::make_not_null(std::any_cast<std::unordered_map<Vertex const *, Vertex const *> *>(pJob->getParam(4)));
    auto const pContractionRoots = gsl::make_not_null(
            std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(3)));
    {
        std::lock_guard<std::mutex> guard(
                std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT

        auto const iter = pContractionTargets->find(pOrder->startVertex);
        auto const *const contractTo = iter == pContractionTargets->end() ? pOrder->startVertex : iter->second;

        pContractionRoots->insert(contractTo);
        pContractionRoots->erase(pOrder->startVertex);
    }

    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void contract(gsl::not_null<Job const *> const pJob) {
    auto const pOrder = make_not_null_and_const(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
    auto const pContractionRoots = gsl::make_not_null(
            std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(3)));
    auto const iter = pContractionRoots->find(pOrder->endVertex);

    if (iter == pContractionRoots->end()) {
        std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
        return;
    }

    auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(4)));
    std::unordered_map<std::string const *, VertexMatch const *> matches;
    std::for_each(std::begin(pOrder->ids), std::end(pOrder->ids), [&](auto const id) {
        auto const *pVertexMatches = pMatchMap->getVertexMatch(pOrder->startVertex->getId(), id);
        if (pVertexMatches != nullptr) {
            matches.insert({&id, pVertexMatches});
        }
    });

    std::lock_guard<std::mutex> guard(
            std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
    auto const pContainElements = gsl::make_not_null(
            std::any_cast<std::unordered_map<Vertex const *, std::unique_ptr<ContainElement const> const> *>(
                    pJob->getParam(2)));
    pContainElements->insert(std::make_pair(
            pOrder->endVertex, std::make_unique<ContainElement>(ContainElement{
                    std::move(matches), pOrder->startVertex, pOrder->startVertex->getNanoporeLength(),
                    pOrder->score, pOrder->direction, pOrder->isPrimary})));

    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findDeletableEdges(gsl::not_null<Job const *> const pJob) {
    auto const pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(1)));
    std::vector<EdgeOrder> filteredOrders;
    std::remove_copy_if(std::begin(pEdge->getEdgeOrders()), std::end(pEdge->getEdgeOrders()),
                        std::back_inserter(filteredOrders), [](auto const &order) { return order.isContained; });

    if (filteredOrders.empty()) {
        auto const pDeletableEdges = gsl::make_not_null(
                std::any_cast<std::set<Edge const *const> *>(pJob->getParam(2)));
        std::lock_guard<std::mutex> guard(
                std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get()); // NOLINT
        pDeletableEdges->insert(pEdge);
    }

    pEdge->replaceOrders(std::move(filteredOrders));
    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void computeBitweight(gsl::not_null<Job const *> const pJob) {
    auto const pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(1)));
    auto const &orders = pEdge->getEdgeOrders();

    if (orders.empty()) {
        std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
        return;
    }

    if (pEdge->isShadow()) {
        auto const initial = orders[0].direction;
        auto const findOther = std::find_if(orders.begin(), orders.end(),
                                            [&](auto const &order) { return order.direction != initial; }) !=
                               orders.end();
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
    auto const pEdge = gsl::make_not_null(std::any_cast<Edge const *>(pJob->getParam(2)));
    auto const pMaxSpanTree = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(3)));
    if (pEdge->getConsensusDirection() != Direction::e_NONE &&
        !pMaxSpanTree->hasEdge(std::make_pair(pEdge->getVertices().first, pEdge->getVertices().second))) {
        auto const shortestPath = GraphUtil::getShortestPath(pMaxSpanTree, pEdge->getVertices());
        lazybastard::Toggle direction = pEdge->getConsensusDirection();
        auto const baseWeight = static_cast<double>(pEdge->getWeight());

        std::vector<double> weights;
        auto const pGraph = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(1)));
        for (auto it = shortestPath.begin(); it != std::prev(shortestPath.end()); ++it) {
            auto const pPathEdge = make_not_null_and_const(pGraph->getEdge(std::make_pair(*it, *std::next(it))));

            direction *= pPathEdge->getConsensusDirection();
            weights.push_back(static_cast<double>(pPathEdge->getWeight()));
        }

        if (!direction && !weights.empty()) {
            auto const iterMinWeight = std::min_element(weights.begin(), weights.end());
            auto const iterMaxWeight = std::max_element(weights.begin(), weights.end());

            auto const pDeletableEdges = gsl::make_not_null(
                    std::any_cast<std::set<Edge const *const> *>(pJob->getParam(4)));
            if (*iterMinWeight < baseWeight || (baseWeight * BASE_WEIGHT_MULTIPLICATOR >= *iterMinWeight &&
                                                *iterMinWeight < *iterMaxWeight * MAX_WEIGHT_MULTIPLICATOR)) {
                auto const minWeightIdx = std::distance(weights.begin(), iterMinWeight);
                auto const *const pDeletableEdge = pGraph->getEdge(std::make_pair(
                        *std::next(shortestPath.begin(), minWeightIdx),
                        *std::next(shortestPath.begin(), minWeightIdx + 1)));

                std::unique_lock<std::mutex> lck(
                        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
                pDeletableEdges->insert(pDeletableEdge);
            }

            {
                std::unique_lock<std::mutex> lck(
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
            make_not_null_and_const(
                    std::any_cast<std::vector<lazybastard::graph::Vertex *> const *>(pJob->getParam(3)));
    auto const pSubGraph = pGraph->getSubgraph(*pConnectedComponent);
    auto const subGraphVertices = pSubGraph->getVertices();
    auto const maxNplVertexIter = std::max_element(
            std::begin(subGraphVertices), std::end(subGraphVertices),
            [](Vertex const *v1, Vertex const *v2) { return v1->getNanoporeLength() < v2->getNanoporeLength(); });
    auto *const pMaxNplVertex = maxNplVertexIter == pSubGraph->getVertices().end() ? nullptr : *maxNplVertexIter;
    auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(2)));

    if (pMaxNplVertex != nullptr) {
        auto const pDiGraph = lazybastard::getDirectionGraph(pGraph, pSubGraph.get(), pMaxNplVertex);
        auto const paths = lazybastard::linearizeGraph(pDiGraph.get());

        Id2OverlapMap id2OverlapMap;
        auto const assemblyIdx = std::any_cast<std::size_t>(pJob->getParam(4));
        std::for_each(paths.begin(), paths.end(), [&](auto const &path) {
            lazybastard::assemblePath(pGraph, pMatchMap, &id2OverlapMap, &path, pDiGraph.get(), assemblyIdx,
                                      std::any_cast<std::reference_wrapper<OutputWriter>>(pJob->getParam(5)).get());
        });
    }

    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------