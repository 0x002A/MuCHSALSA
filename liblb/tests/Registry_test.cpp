#include <gtest/gtest.h>

#include "Registry.h"

TEST(RegistryTest, Test) {
  auto registry = lazybastard::Registry();

  ASSERT_EQ(registry["abc"], 0);
  ASSERT_EQ(registry["def"], 1);

  registry.clear();

  ASSERT_EQ(registry["def"], 0);
}
