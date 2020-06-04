#pragma once

#include <thread>
#include <atomic>
#include <condition_variable>

namespace lazybastard {
namespace threading {

/**
 * Class enabling a thread to block until a bunch of jobs have completed.
 *
 * A WaitGroup holds a count of jobs which have to complete until it returns and thread-safe methods to increment and
 * decrement this counter.
 * This concept is borrowed from the go programming language.
 */
class WaitGroup {
public:
  /**
   * Class constructor which creates a new instance.
   */
  WaitGroup() = default;

  /**
   * Copying is disallowed.
   */
   WaitGroup(const WaitGroup& ) = delete;

   /**
    * Copy assignment is disallowed.
    */
    WaitGroup& operator=(const WaitGroup& ) = delete;

  /**
   * Moving is disallowed.
   */
   WaitGroup(WaitGroup&& ) = delete;

   /**
    * Move assignment is disallowed.
    */
   WaitGroup& operator=(WaitGroup&& ) = delete;

   /**
    * Raises the job counter by the number supplied as argument.
    *
    * @param newJobs The number of jobs to raise the job counter by.
    */
   void add(std::size_t newJobs);

   /**
    * Notifies the WaitGroup that a job has completed.
    */
   void done();

   /**
    * Blocks the caller until all jobs assigned to the WaitGroup completed.
    */
   void wait();
private:
  std::atomic<bool> m_waitLock{false}; /*!< Toggle marking the WaitGroup as waiting and locking the add function */
  std::atomic<std::size_t> m_jobCount{0}; /*!< Number of jobs assigned */
  std::mutex m_mutex; /*!< Mutex used for locking */
  std::condition_variable m_cv; /*!< Conditional variable used to notify threads */
};

}
}
