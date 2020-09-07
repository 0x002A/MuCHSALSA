#include <fstream>
#include <iostream>

#include <lb/BlastFileReader.h>
#include <lb/graph/Graph.h>
#include <lb/matching/MatchMap.h>
#include <lb/threading/ThreadPool.h>

#include "Application.h"

auto main(int argc, char *argv[]) -> int {
  Application app(argc, argv);

  if (!app.checkIntegrity()) {
    std::cerr << "Paths are pointing to invalid/unuseable locations" << std::endl;

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

  matchMap.calculateEgdes();

  std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

  return 0;
}
