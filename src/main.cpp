#include <iostream>
#include <fstream>

#include <lb/BlastFileReader.h>
#include <lb/threading/ThreadPool.h>

#include "Application.h"

auto
main(int argc, char *argv[]) -> int
{
  Application app(argc, argv);

  if (!app.checkIntegrity()) {
    std::cerr << "Paths are pointing to invalid/unuseable locations" << std::endl;

    return -1;
  }

  // Keep one free for the main program
  auto threadCount = app.getThreadCount() - 1;
  auto threadPool = lazybastard::threading::ThreadPool(threadCount);

  // Read BLAST file
  if (std::ifstream inputStream{app.getContigsFilePath(), std::ios::binary | std::ios::in}) {
    lazybastard::BlastFileReader blastReader(inputStream, &threadPool);
  } else {
    std::cerr << "Can't open BLAST file for reading" << std::endl;

    return -1;
  }

  return 0;
}
