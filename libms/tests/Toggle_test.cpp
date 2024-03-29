#include <gtest/gtest.h>

#include "types/Toggle.h"

TEST(ToggleTest, BasicTest) {
  muchsalsa::Toggle const tTrue  = true;
  muchsalsa::Toggle const tFalse = false;

  muchsalsa::Toggle tShouldBeTrue = false;
  tShouldBeTrue *= false;

  ASSERT_EQ(((bool)tTrue), true);
  ASSERT_EQ(((bool)tFalse), false);
  ASSERT_EQ(((bool)!tTrue), false);
  ASSERT_EQ(((bool)!tFalse), true);
  ASSERT_EQ(((bool)tShouldBeTrue), true);
  ASSERT_EQ(((bool)!tShouldBeTrue), false);
  ASSERT_EQ(((bool)(tTrue && tTrue)), true);
  ASSERT_EQ(((bool)(tTrue && tFalse)), false);
  ASSERT_EQ(((bool)(tFalse && tTrue)), false);
  ASSERT_EQ(((bool)(tFalse && tFalse)), false);
  ASSERT_EQ(((bool)(tTrue == tTrue)), true);
  ASSERT_EQ(((bool)(tTrue != tTrue)), false);
  ASSERT_EQ(((bool)(tFalse == tFalse)), true);
  ASSERT_EQ(((bool)(tFalse != tFalse)), false);
}