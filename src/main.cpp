#include <iostream>

#include <lb/LazyBastard.h>

auto
main(int argc, char *argv[]) -> int
{
  LazyBastard lb;
  std::cout << lb.getMagicValue() << std::endl;
  return 0;
}
