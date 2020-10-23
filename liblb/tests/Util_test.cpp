#include <stdexcept>

#include <gtest/gtest.h>

#include "Util.h"

using lazybastard::util::check_pointers;
using lazybastard::util::less_than;
using lazybastard::util::sortPair;

TEST(UtilTest, CorrectlyChecksPointers) {
  std::size_t s = 42;
  auto vp = &s;
  std::size_t *np = nullptr;

  ASSERT_NO_THROW(check_pointers(&s));
  ASSERT_NO_THROW(check_pointers(&s, vp));
  ASSERT_THROW(check_pointers(np), std::invalid_argument);
  ASSERT_THROW(check_pointers(&s, np), std::invalid_argument);

  // Can't be tested as it forces compilation to fail
  // ASSERT_THROW(check_pointers(nullptr), std::invalid_argument);
  // ASSERT_THROW(check_pointers(s), std::invalid_argument);
}

TEST(UtilTest, CorrectlyComparesPointers) {
  std::size_t s1 = 42;
  std::size_t s2 = 84;

  ASSERT_TRUE(less_than(s1, s2));
  ASSERT_TRUE(less_than(&s1, &s2));
  ASSERT_FALSE(less_than(s2, s1));
  ASSERT_FALSE(less_than(&s2, &s1));
}

TEST(UtilTest, CorrectlySortsPairs) {
  std::size_t s1 = 42;
  std::size_t s2 = 84;

  auto p1 = std::make_pair(s1, s2);
  auto p2 = std::make_pair(&s1, &s2);
  auto p3 = std::make_pair(s2, s1);
  auto p4 = std::make_pair(&s2, &s1);

  sortPair(p1);
  sortPair(p2);
  sortPair(p3);
  sortPair(p4);

  ASSERT_EQ(p1.first, s1);
  ASSERT_EQ(p1.second, s2);

  ASSERT_EQ(*p2.first, s1);
  ASSERT_EQ(*p2.second, s2);

  ASSERT_EQ(p3.first, s1);
  ASSERT_EQ(p3.second, s2);

  ASSERT_EQ(*p4.first, s1);
  ASSERT_EQ(*p4.second, s2);
}
