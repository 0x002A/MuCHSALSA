#include "threading/WaitGroup.h"

#include <stdexcept>

namespace lazybastard::threading {

void WaitGroup::add(std::size_t newJobs) {
  if (m_waitLock) {
    throw std::logic_error("Adding jobs to an already waiting wait group is illegal behaviour.");
  }

  m_jobCount += newJobs;
  m_cv.notify_one();
};

void WaitGroup::done() {
  m_jobCount -= 1;
  m_cv.notify_one();
}

void WaitGroup::wait() {
  // prevent future calls to add
  m_waitLock = true;

  std::unique_lock<std::mutex> lck(m_mutex);
  m_cv.wait(lck, [this] { return m_jobCount.load() == 0; });
}

} // namespace lazybastard::threading
