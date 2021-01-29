#include <algorithm>
#include <any>
#include <cstddef>
#include <fstream>
#include <gsl/pointers>
#include <gsl/span>
#include <iostream>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>

#include <lb/BlastFileReader.h>
#include <lb/Prokrastinator.h>
#include <lb/Util.h>
#include <lb/graph/Edge.h>
#include <lb/graph/Graph.h>
#include <lb/graph/Vertex.h>
#include <lb/matching/MatchMap.h>
#include <lb/threading/Job.h>
#include <lb/threading/ThreadPool.h>
#include <lb/threading/WaitGroup.h>

#include "Application.h"

struct ContainElement { // NOLINT
  std::unordered_map<std::string const *, lazybastard::matching::VertexMatch const *> const matches;
  lazybastard::graph::Vertex const *const nano;
  std::size_t const nanoporeLength;
  std::size_t const score;
  bool const direction;
  bool const isPrimary;
};

void chainingAndOverlaps(gsl::not_null<lazybastard::threading::Job const *> pJob);
void findContractionEdges(gsl::not_null<lazybastard::threading::Job const *> pJob);
void findContractionTargets(gsl::not_null<lazybastard::threading::Job const *> pJob);
void findDeletableVertices(gsl::not_null<lazybastard::threading::Job const *> pJob);
void contract(gsl::not_null<lazybastard::threading::Job const *> pJob);
void findDeletableEdges(gsl::not_null<lazybastard::threading::Job const *> pJob);
void computeBitweight(gsl::not_null<lazybastard::threading::Job const *> pJob);

auto main(int const argc, char const *argv[]) -> int {
  gsl::span<char const *> const args = {argv, static_cast<std::size_t>(argc)};
  Application app(args);

  if (!app.checkIntegrity()) {
    std::cerr << "Paths are pointing to invalid/unusable locations" << std::endl;

    return -1;
  }

  // Keep one free for the main program
  auto threadCount = app.getThreadCount() - 1;
  auto threadPool = lazybastard::threading::ThreadPool(threadCount);

  // Read BLAST file
  auto graph = lazybastard::graph::Graph();
  auto matchMap = lazybastard::matching::MatchMap(&threadPool, &graph);
  if (std::ifstream inputStream{app.getContigsFilePath(), std::ios::binary | std::ios::in}) {
    auto blastReader = lazybastard::BlastFileReader(&threadPool, inputStream, &graph, &matchMap);
    blastReader.read();
  } else {
    std::cerr << "Can't open BLAST file for reading" << std::endl;

    return -1;
  }

  matchMap.calculateEdges();

  std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

  try {
    lazybastard::threading::WaitGroup wg;

    auto chainingJob = [](lazybastard::threading::Job const *const pJob) { chainingAndOverlaps(pJob); };
    for (auto const &[vertexID, edges] : graph.getAdjacencyList()) {
      for (auto const &[targetID, edge] : edges) {
        wg.add(1);

        auto job = lazybastard::threading::Job(chainingJob, &wg,
                                               static_cast<lazybastard::matching::MatchMap const *const>(&matchMap),
                                               static_cast<lazybastard::graph::Graph const *const>(&graph), edge.get());
        threadPool.addJob(std::move(job));
      }
    }
    wg.wait();

    std::mutex mutex;
    std::unordered_map<lazybastard::graph::Edge const *, lazybastard::graph::EdgeOrder const *> contractionEdges;
    auto contractionEdgesJob = [](lazybastard::threading::Job const *const pJob) { findContractionEdges(pJob); };
    for (auto const &[edgeID, edges] : graph.getAdjacencyList()) {
      LB_UNUSED(edgeID);

      for (auto const &[targetEdgeID, pEdge] : edges) {
        LB_UNUSED(targetEdgeID);

        for (auto const &order : pEdge->getEdgeOrders()) {
          wg.add(1);

          auto job = lazybastard::threading::Job(contractionEdgesJob, &wg,
                                                 static_cast<lazybastard::graph::Graph const *const>(&graph), order,
                                                 &contractionEdges, std::reference_wrapper<std::mutex>(mutex));
          threadPool.addJob(std::move(job));
        }
      }
    }
    wg.wait();

    std::cout << "Number of contraction edges: " << contractionEdges.size() << std::endl;

    std::unordered_map<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *> contractionTargets;
    auto contractionTargetsJob = [](lazybastard::threading::Job const *const pJob) { findContractionTargets(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = lazybastard::threading::Job(contractionTargetsJob, &wg, pOrder, &contractionTargets,
                                             std::reference_wrapper<std::mutex>(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::set<std::string const *const> deletableVertices;
    std::set<lazybastard::graph::Vertex const *const> contractionRoots;
    auto deletableVerticesJob = [](lazybastard::threading::Job const *const pJob) { findDeletableVertices(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = lazybastard::threading::Job(deletableVerticesJob, &wg, pOrder, &deletableVertices, &contractionRoots,
                                             &contractionTargets, std::reference_wrapper<std::mutex>(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::unordered_map<lazybastard::graph::Vertex const *, std::unique_ptr<ContainElement const> const> containElements;
    auto contractionJob = [](lazybastard::threading::Job const *const pJob) { contract(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = lazybastard::threading::Job(contractionJob, &wg, pOrder, &containElements, &contractionRoots,
                                             &matchMap, std::reference_wrapper<std::mutex>(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::cout << "Vertices to become deleted: " << deletableVertices.size() << std::endl;

    for (auto const *const target : deletableVertices) {
      graph.deleteVertex(target);
    }
    deletableVertices.clear();

    std::set<lazybastard::graph::Edge const *const> deletableEdges;
    auto deletableEdgesJob = [](lazybastard::threading::Job const *const pJob) { findDeletableEdges(pJob); };
    for (auto const &[edgeID, edges] : graph.getAdjacencyList()) {
      LB_UNUSED(edgeID);

      for (auto const &[targetEdgeID, pEdge] : edges) {
        LB_UNUSED(targetEdgeID);

        wg.add(1);
        auto job = lazybastard::threading::Job(deletableEdgesJob, &wg, pEdge.get(), &deletableEdges,
                                               std::reference_wrapper<std::mutex>(mutex));
        threadPool.addJob(std::move(job));
      }
    }
    wg.wait();

    std::cout << "Edges to become deleted: " << deletableEdges.size() << std::endl;

    for (auto const *const target : deletableEdges) {
      graph.deleteEdge(target);
    }
    deletableEdges.clear();

    std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

    auto computeBitweightJob = [](lazybastard::threading::Job const *const pJob) { computeBitweight(pJob); };
    for (auto const &[edgeID, edges] : graph.getAdjacencyList()) {
      LB_UNUSED(edgeID);

      for (auto const &[targetEdgeID, pEdge] : edges) {
        LB_UNUSED(targetEdgeID);

        wg.add(1);
        auto job = lazybastard::threading::Job(computeBitweightJob, &wg, pEdge.get());
        threadPool.addJob(std::move(job));
      }
    }
    wg.wait();
  } catch (std::exception const &e) {
    std::cerr << "Exception occurred: " << '\n' << e.what();
  }
  return 0;
}

void chainingAndOverlaps(gsl::not_null<lazybastard::threading::Job const *> const pJob) {
  std::set<std::string const *const, lazybastard::util::LTCmp<std::string const *const>> plusIDs;
  std::set<const std::string *const, lazybastard::util::LTCmp<std::string const *const>> minusIDs;

  auto const *const pMatchMap = std::any_cast<lazybastard::matching::MatchMap const *const>(pJob->getParam(1));
  auto *const pEdge = std::any_cast<lazybastard::graph::Edge *const>(pJob->getParam(3));
  auto const edgeMatches = pMatchMap->getEdgeMatches().find(pEdge->getID());

  if (edgeMatches == pMatchMap->getEdgeMatches().end()) {
    std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
    return;
  }

  for (auto const &[illuminaID, edgeMatch] : edgeMatches->second) {

    if (edgeMatch->direction) {
      plusIDs.insert(&illuminaID);
    } else {
      minusIDs.insert(&illuminaID);
    }
  }

  auto const *const pGraph = std::any_cast<lazybastard::graph::Graph const *const>(pJob->getParam(2));
  auto plusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pGraph, pEdge, plusIDs, true);
  auto minusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pGraph, pEdge, minusIDs, false);

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
  } else if (combinedSize == 1) {
    auto const &path = !minusPaths.empty() ? minusPaths[0] : plusPaths[0];
    pEdge->setShadow(std::get<2>(path));
  }

  for (const auto &plusPath : plusPaths) {
    pEdge->appendOrder(lazybastard::computeOverlap(pMatchMap, pGraph, std::get<0>(plusPath), pEdge, false,
                                                   std::get<1>(plusPath), std::get<2>(plusPath)));
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

void findContractionEdges(gsl::not_null<lazybastard::threading::Job const *> const pJob) {
  auto const *const pGraph = std::any_cast<lazybastard::graph::Graph const *const>(pJob->getParam(2));
  auto const *const pOrder = std::any_cast<lazybastard::graph::EdgeOrder const *const>(pJob->getParam(2));
  if (pOrder->isContained && pOrder->isPrimary) {
    auto isSane = true;
    auto const edges = pGraph->getEdgesOfVertex(pOrder->startVertex->getID());
    for (auto const &[targetID, edge] : edges) {
      if (*targetID == pOrder->endVertex->getID() || edge->isShadow()) {
        continue;
      }

      auto vertices = std::make_pair(pOrder->endVertex->getID(), *targetID);
      isSane &= pGraph->hasEdge(vertices);
      if (!isSane) {
        break;
      }

      auto const *const pTargetVertex =
          edge->getVertices().first != pOrder->startVertex ? edge->getVertices().first : edge->getVertices().second;

      isSane &= lazybastard::sanityCheck(pGraph, pOrder->startVertex, pOrder->endVertex, pTargetVertex, pOrder);

      if (!isSane) {
        break;
      }

      auto *const contractionEdges = std::any_cast<
          std::unordered_map<lazybastard::graph::Edge const *, lazybastard::graph::EdgeOrder const *> *const>(
          pJob->getParam(3));
      std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(4)).get());
      contractionEdges->insert({edge, pOrder});
    }
  }
  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

void findContractionTargets(gsl::not_null<lazybastard::threading::Job const *> const pJob) {
  auto const *const pOrder = std::any_cast<lazybastard::graph::EdgeOrder const *const>(pJob->getParam(1));
  auto *const pContractionTargets =
      std::any_cast<std::unordered_map<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *> *const>(
          pJob->getParam(2));
  {
    std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get());

    auto const iter = pContractionTargets->find(pOrder->endVertex);
    auto const *const contractTo = iter == pContractionTargets->end() ? pOrder->endVertex : iter->second;

    if (pContractionTargets->find(pOrder->startVertex) == pContractionTargets->end()) {
      pContractionTargets->insert({pOrder->startVertex, contractTo});
    }
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

void findDeletableVertices(gsl::not_null<const lazybastard::threading::Job *> pJob) {
  auto const *const pOrder = std::any_cast<lazybastard::graph::EdgeOrder const *const>(pJob->getParam(1));
  auto *const pDeletableVertices = std::any_cast<std::set<std::string const *const> *const>(pJob->getParam(2));
  {
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
    pDeletableVertices->insert(&pOrder->startVertex->getID());
  }

  auto *const pContractionTargets =
      std::any_cast<std::unordered_map<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *> *const>(
          pJob->getParam(4));
  auto *const pContractionRoots =
      std::any_cast<std::set<lazybastard::graph::Vertex const *const> *const>(pJob->getParam(3));
  {
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT

    auto const iter = pContractionTargets->find(pOrder->startVertex);
    auto const *const contractTo = iter == pContractionTargets->end() ? pOrder->startVertex : iter->second;

    pContractionRoots->insert(contractTo);
    pContractionRoots->erase(pOrder->startVertex);
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

void contract(gsl::not_null<lazybastard::threading::Job const *> const pJob) {
  auto const *const pOrder = std::any_cast<lazybastard::graph::EdgeOrder const *const>(pJob->getParam(1));
  auto *const pContractionRoots =
      std::any_cast<std::set<lazybastard::graph::Vertex const *const> *const>(pJob->getParam(3));
  auto const iter = pContractionRoots->find(pOrder->endVertex);

  if (iter == pContractionRoots->end()) {
    std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
    return;
  }

  auto const *const pMatchMap = std::any_cast<lazybastard::matching::MatchMap const *const>(pJob->getParam(4));
  std::unordered_map<std::string const *, lazybastard::matching::VertexMatch const *> matches;
  for (auto const *id : pOrder->ids) {
    auto const *pVertexMatches = pMatchMap->getVertexMatch(pOrder->startVertex->getID(), *id);
    if (pVertexMatches != nullptr) {
      matches.insert({id, pVertexMatches});
    }
  }

  std::lock_guard<std::mutex> guard(
      std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
  auto *const pContainElements = std::any_cast<
      std::unordered_map<lazybastard::graph::Vertex const *, std::unique_ptr<ContainElement const> const> *const>(
      pJob->getParam(3));
  pContainElements->insert(std::move(std::make_pair(
      pOrder->endVertex, std::make_unique<ContainElement>(ContainElement{
                             std::move(matches), pOrder->startVertex, pOrder->startVertex->getNanoporeLength(),
                             pOrder->score, pOrder->direction, pOrder->isPrimary}))));

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

void findDeletableEdges(gsl::not_null<lazybastard::threading::Job const *> const pJob) {
  auto *const pEdge = std::any_cast<lazybastard::graph::Edge *const>(pJob->getParam(1));
  std::vector<lazybastard::graph::EdgeOrder> filteredOrders;
  std::remove_copy_if(pEdge->getEdgeOrders().begin(), pEdge->getEdgeOrders().end(), std::back_inserter(filteredOrders),
                      [](auto const &order) { return order.isContained; });

  if (filteredOrders.empty()) {
    auto *const pDeletableEdges =
        std::any_cast<std::set<lazybastard::graph::Edge const *const> *const>(pJob->getParam(2));
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get()); // NOLINT
    pDeletableEdges->insert(pEdge);
  }

  pEdge->replaceOrders(std::move(filteredOrders));
  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

void computeBitweight(gsl::not_null<lazybastard::threading::Job const *> const pJob) {
  auto *const pEdge = std::any_cast<lazybastard::graph::Edge *const>(pJob->getParam(1));
  auto const &orders = pEdge->getEdgeOrders();

  if (orders.empty()) {
    std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
    return;
  }

  if (pEdge->isShadow()) {
    auto const &initial = orders[0].direction;
    auto const findOther = std::find_if(orders.begin(), orders.end(),
                                        [&](auto const &order) { return order.direction != initial; }) != orders.end();
    if (!findOther) {
      pEdge->setConsensusDirection(orders[0].direction);
    }
  } else {
    pEdge->setWeight(orders[0].score);
    pEdge->setConsensusDirection(orders[0].direction);
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}