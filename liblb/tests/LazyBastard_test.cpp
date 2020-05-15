#include <gtest/gtest.h>

#include <LazyBastard.h>

TEST(LazyBastardTest, ReturnsMagicValue) {
    LazyBastard lb;
    EXPECT_EQ(lb.getMagicValue(),  42);
}
