#include <gtest/gtest.h>

#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "IO.h"

extern char *pTestDataPath;

TEST(IOTest, ReadlineTest) {
  auto textInputFile = std::string(pTestDataPath);
  textInputFile.append("/text.txt");

  std::vector<std::string> result;

  auto *const pFile      = fopen(textInputFile.c_str(), "rbe");
  char       *pLine      = nullptr;
  std::size_t sizeBuffer = 0;
  auto        ret        = muchsalsa::readline(&pLine, &sizeBuffer, pFile);
  while (ret != -1) {
    result.push_back(std::string(pLine));

    ret = muchsalsa::readline(&pLine, &sizeBuffer, pFile);
  }

  ::operator delete(pLine);
  std::fclose(pFile);

  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0], "MuCHSALSA\n");
  ASSERT_EQ(result[1], "is\n");
  ASSERT_EQ(result[2], "GREAT");
}
