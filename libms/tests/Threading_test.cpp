#include <chrono>
#include <vector>

#include <gtest/gtest.h>

#include "threading/Job.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"

std::mutex m;

TEST(ThreadingTest, ReturnsCorrectResult) {
  const std::vector<std::size_t> nums{10, 20, 30, 40, 50, 60, 70, 80};
  std::vector<std::size_t>       res{};

  auto threadPool = muchsalsa::threading::ThreadPool(2);

  muchsalsa::threading::WaitGroup wg;
  auto                            jobFn = [&res](muchsalsa::threading::Job const *const pJob) {
    const auto number = std::any_cast<size_t>(pJob->getParam(1));

    // push value
    {
      std::scoped_lock<std::mutex> lock(m);
      res.push_back(number * 10);
    }

    // simulate a long running job to test wait group behaviour
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    std::any_cast<muchsalsa::threading::WaitGroup *const>(pJob->getParam(0))->done();
  };

  for (std::size_t num : nums) {
    wg.add(1);

    auto job = muchsalsa::threading::Job(jobFn, &wg, num);
    threadPool.addJob(std::move(job));
  }

  wg.wait();

  ASSERT_EQ(res.size(), 8);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 100) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 200) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 300) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 400) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 500) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 600) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 700) == 1);
  ASSERT_TRUE(std::count(res.begin(), res.end(), 800) == 1);
}

TEST(ThreadingTest, ThrowsNoException) {
  auto threadPool = muchsalsa::threading::ThreadPool(2);

  muchsalsa::threading::WaitGroup wg;
  auto                            jobFn = [](muchsalsa::threading::Job const *const pJob) {
    // simulate a long running job to test wait group behaviour
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    std::any_cast<muchsalsa::threading::WaitGroup *const>(pJob->getParam(0))->done();
  };

  for (std::size_t i = 0; i <= 8; ++i) {
    wg.add(1);

    auto job = muchsalsa::threading::Job(jobFn, &wg);
    threadPool.addJob(std::move(job));
  }

  wg.wait();
  ASSERT_NO_THROW(wg.add(1));
}
