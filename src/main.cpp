#include <fstream>
#include <gsl/pointers>
#include <gsl/span>
#include <iostream>
#include <set>

#include <lb/BlastFileReader.h>
#include <lb/Prokrastinator.h>
#include <lb/Util.h>
#include <lb/graph/Edge.h>
#include <lb/graph/Graph.h>
#include <lb/matching/MatchMap.h>
#include <lb/threading/Job.h>
#include <lb/threading/ThreadPool.h>
#include <lb/threading/WaitGroup.h>

#include "Application.h"

void chainingAndOverlaps(gsl::not_null<const lazybastard::threading::Job *> pJob);

auto main(int argc, char *argv[]) -> int {
  gsl::span<char *> args = {argv, static_cast<std::size_t>(argc)};
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
    auto jobFn = [](const lazybastard::threading::Job *pJob) { chainingAndOverlaps(pJob); };
    for (const auto &[vertexID, edges] : graph.getAdjacencyList()) {
      for (const auto &[targetID, edge] : edges) {
        wg.add(1);

        auto job = lazybastard::threading::Job(jobFn, &wg, &matchMap, &graph, edge.get());
        threadPool.addJob(std::move(job));
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "Exception occurred: " << '\n' << e.what();
  }
  return 0;
}

void chainingAndOverlaps(gsl::not_null<const lazybastard::threading::Job *> pJob) {
  std::set<const std::string *, lazybastard::util::LTCmp<const std::string *>> plusIDs;
  std::set<const std::string *, lazybastard::util::LTCmp<const std::string *>> minusIDs;

  auto *const pMatchMap = std::any_cast<lazybastard::matching::MatchMap *>(pJob->getParam(1));
  auto *const pEdge = std::any_cast<lazybastard::graph::Edge *>(pJob->getParam(3));
  const auto &edgeMatches = pMatchMap->getEdgeMatches().find(pEdge->getID());

  if (edgeMatches == pMatchMap->getEdgeMatches().end()) {
    return;
  }

  for (const auto &[illuminaID, edgeMatch] : edgeMatches->second) {

    if (edgeMatch->direction) {
      plusIDs.insert(&illuminaID);
    } else {
      minusIDs.insert(&illuminaID);
    }
  }

  auto *const pGraph = std::any_cast<lazybastard::graph::Graph *>(pJob->getParam(2));
  auto plusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pGraph, pEdge, plusIDs, true);
  auto minusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pGraph, pEdge, minusIDs, false);

  bool hasPrimary = true;
  for_each(plusPaths.begin(), plusPaths.end(), [&](const auto &plusPath) {
    if (!std::get<2>(plusPath)) {
      hasPrimary = false;
    }
  });

  if (!hasPrimary) {
    bool hasMinusPrimary = true;
    for_each(minusPaths.begin(), minusPaths.end(), [&](const auto &minusPath) {
      if (!std::get<2>(minusPath)) {
        hasMinusPrimary = false;
      }
    });

    hasPrimary = hasMinusPrimary;
  }

  if (hasPrimary) {
    plusPaths.erase(
        std::remove_if(plusPaths.begin(), plusPaths.end(), [](const auto &plusPath) { return !std::get<2>(plusPath); }),
        plusPaths.end());
    minusPaths.erase(std::remove_if(minusPaths.begin(), minusPaths.end(),
                                    [](const auto &minusPath) { return !std::get<2>(minusPath); }),
                     minusPaths.end());
  }

  bool hasMulti = true;
  for_each(plusPaths.begin(), plusPaths.end(), [&](const auto &plusPath) {
    if (std::get<0>(plusPath).size() <= 1) {
      hasMulti = false;
    }
  });

  if (!hasMulti) {
    bool hasMinusMulti = true;
    for_each(minusPaths.begin(), minusPaths.end(), [&](const auto &minusPath) {
      if (std::get<0>(minusPath).size() <= 1) {
        hasMinusMulti = false;
      }
    });

    hasMulti = hasMinusMulti;
  }

  if (hasMulti) {
    plusPaths.erase(std::remove_if(plusPaths.begin(), plusPaths.end(),
                                   [](const auto &plusPath) { return std::get<0>(plusPath).size() <= 1; }),
                    plusPaths.end());
    minusPaths.erase(std::remove_if(minusPaths.begin(), minusPaths.end(),
                                    [](const auto &minusPath) { return std::get<0>(minusPath).size() <= 1; }),
                     minusPaths.end());
  }

  const auto combinedSize = plusPaths.size() + minusPaths.size();
  if (combinedSize > 1) {
    pEdge->setShadow(true);
  } else if (combinedSize == 1) {
    const auto &path = !minusPaths.empty() ? minusPaths[0] : plusPaths[0];
    pEdge->setShadow(std::get<2>(path));
  }

  for (const auto &plusPath : plusPaths) {
    pEdge->appendOrder(lazybastard::computeOverlap(pMatchMap, pGraph, std::get<0>(plusPath), pEdge, false,
                                                   std::get<1>(plusPath), std::get<2>(plusPath)));
  }
}