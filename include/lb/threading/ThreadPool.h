#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "threading/Job.h"

namespace lazybastard::threading {

/**
 * Class for thread management.
 *
 * A ThreadPool holds all the threads available for work to be done.
 * Instances of this class are designed to be thread-safe.
 */
class ThreadPool {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param threadCount the number of threads to start
   */
  explicit ThreadPool(std::size_t threadCount);

  /**
   * Copying is disallowed.
   */
  ThreadPool(ThreadPool const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  ThreadPool &operator=(ThreadPool const &) = delete;

  /**
   * Moving is disallowed.
   */
  ThreadPool(ThreadPool &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  ThreadPool &operator=(ThreadPool &&) = delete;

  /**
   * Class destructor.
   */
  ~ThreadPool();

  /**
   * Adds a Job to the queue of the thread pool.
   *
   * @param job an rvalue reference to the Job to be added.
   */
  void addJob(Job &&job);

private:
  std::vector<std::thread> m_threads;       /*!< Store containing the available threads */
  std::queue<Job> m_jobs;                   /*!< std::queue of Jobs to be executed */
  std::mutex m_mutex;                       /*!< std::mutex for securing the parallel use of the std::queue */
  std::atomic<bool> m_terminatePool{false}; /*!< Bool indicating whether the ThreadPool is going to be terminated */
  std::condition_variable m_condition;      /*!< std::conditional_variable used to notify threads about new Jobs */
};

} // namespace lazybastard::threading
