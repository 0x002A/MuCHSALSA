#include <algorithm>
#include <any>
#include <cstddef>
#include <fstream>
#include <functional>
#include <gsl/pointers>
#include <gsl/span>
#include <iostream>
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

void chainingAndOverlaps(gsl::not_null<lazybastard::threading::Job const *> pJob);
void findContractionEdges(gsl::not_null<lazybastard::threading::Job const *> pJob);

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

    std::mutex mutexContractionEdges;
    std::unordered_map<std::string const *, lazybastard::graph::EdgeOrder const *> contractionEdges;
    auto contractionJob = [](lazybastard::threading::Job const *const pJob) { findContractionEdges(pJob); };
    for (auto const &[vertexID, edges] : graph.getAdjacencyList()) {
      for (auto const &[targetID, edge] : edges) {
        for (auto const &order : edge->getEdgeOrders()) {
          wg.add(1);

          auto job = lazybastard::threading::Job(
              contractionJob, &wg, static_cast<lazybastard::graph::Graph const *const>(&graph), order,
              &contractionEdges, std::reference_wrapper<std::mutex>(mutexContractionEdges));
          threadPool.addJob(std::move(job));
        }
      }
    }
    wg.wait();

    std::cout << "Number of contraction edges: " << contractionEdges.size() << std::endl;
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
  auto const &edgeMatches = pMatchMap->getEdgeMatches().find(pEdge->getID());

  if (edgeMatches == pMatchMap->getEdgeMatches().end()) {
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
  if (pOrder->contained && pOrder->isPrimary) {
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

      auto const *const pTargetVertex = edge->getVertices().first.get() != pOrder->startVertex
                                            ? edge->getVertices().first.get()
                                            : edge->getVertices().second.get();

      isSane &= lazybastard::sanityCheck(pGraph, pOrder->startVertex, pOrder->endVertex, pTargetVertex, pOrder);

      if (!isSane) {
        break;
      }

      auto *const contractionEdges =
          std::any_cast<std::unordered_map<std::string const *, lazybastard::graph::EdgeOrder const *> *const>(
              pJob->getParam(3));
      std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(4)).get());
      contractionEdges->insert({&edge->getID(), pOrder});
    }
  }
  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}