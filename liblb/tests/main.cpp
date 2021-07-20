#include <gtest/gtest.h>

char *pTestDataPath = nullptr;

auto main(int argc, char **argv) -> int {
  if (argc == 2) {
    pTestDataPath = argv[1];
  }

  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
