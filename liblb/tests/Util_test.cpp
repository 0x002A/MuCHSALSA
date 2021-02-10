#include <memory>
#include <stdexcept>

#include <gtest/gtest.h>

#include "Util.h"

using lazybastard::util::check_pointers;
using lazybastard::util::less_than;
using lazybastard::util::sortPair;

TEST(UtilTest, CorrectlyChecksPointers) {
  std::size_t const s = 42;
  auto const vp = &s;
  std::size_t const *const np = nullptr;

  ASSERT_NO_THROW(check_pointers(&s));
  ASSERT_NO_THROW(check_pointers(&s, vp));
  ASSERT_THROW(check_pointers(np), std::invalid_argument);
  ASSERT_THROW(check_pointers(&s, np), std::invalid_argument);

  // Can't be tested as it forces compilation to fail
  // ASSERT_THROW(check_pointers(nullptr), std::invalid_argument);
  // ASSERT_THROW(check_pointers(s), std::invalid_argument);
}

TEST(UtilTest, CorrectlyComparesPointers) {
  std::size_t const s1 = 42;
  std::size_t const s2 = 84;

  auto const sp1 = std::make_shared<std::size_t>(42);
  auto const sp2 = std::make_shared<std::size_t>(84);

  ASSERT_TRUE(less_than(s1, s2));
  ASSERT_TRUE(less_than(&s1, &s2));
  ASSERT_TRUE(less_than(sp1, sp2));
  ASSERT_FALSE(less_than(s2, s1));
  ASSERT_FALSE(less_than(&s2, &s1));
  ASSERT_FALSE(less_than(sp2, sp1));
}

TEST(UtilTest, CorrectlySortsPairs) {
  std::size_t const s1 = 42;
  std::size_t const s2 = 84;

  auto const sp1 = std::make_shared<std::size_t>(42);
  auto const sp2 = std::make_shared<std::size_t>(84);

  auto p1 = std::make_pair(s1, s2);
  auto p2 = std::make_pair(&s1, &s2);
  auto p3 = std::make_pair(s2, s1);
  auto p4 = std::make_pair(&s2, &s1);
  auto p5 = std::make_pair(std::move(sp2), std::move(sp1));

  sortPair(p1);
  sortPair(p2);
  sortPair(p3);
  sortPair(p4);
  sortPair(p5);

  ASSERT_EQ(p1.first, s1);
  ASSERT_EQ(p1.second, s2);

  ASSERT_EQ(*p2.first, s1);
  ASSERT_EQ(*p2.second, s2);

  ASSERT_EQ(p3.first, s1);
  ASSERT_EQ(p3.second, s2);

  ASSERT_EQ(*p4.first, s1);
  ASSERT_EQ(*p4.second, s2);

  ASSERT_EQ(*p5.first, s1);
  ASSERT_EQ(*p5.second, s2);
}
