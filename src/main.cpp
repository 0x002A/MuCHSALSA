#include <iostream>

#include "Application.h"

auto
main(int argc, char *argv[]) -> int
{
  Application app(argc, argv);

  if (!app.checkIntegrity()) {
    std::cerr << "Paths are pointing to invalid/unuseable locations" << std::endl;

    return -1;
  }

  return 0;
}
