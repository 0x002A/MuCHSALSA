#pragma once

#include <atomic>
#include <condition_variable>
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
   * Class constructor which creates a new instance.
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
   * Class destructor which cleans up.
   */
  ~ThreadPool();

  /**
   * Adds a Job to the queue of the thread pool.
   *
   * @param job The Job to be added.
   */
  void addJob(Job &&job);

private:
  std::vector<std::thread> m_threads;       /*!< Vector of threads available */
  std::queue<Job> m_jobs;                   /*!< Queue of Jobs to be executed */
  std::mutex m_mutex;                       /*!< Mutex for securing the parallel use of the queue */
  std::atomic<bool> m_terminatePool{false}; /*!< Indicator whether the ThreadPool is going to be terminated */
  std::condition_variable m_condition;      /*!< Conditional variable used to notify thread about new Jobs */
};

} // namespace lazybastard::threading
