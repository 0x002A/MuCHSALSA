#include "threading/ThreadPool.h"

#include <cstddef>
#include <utility>

#include "threading/Job.h"

namespace lazybastard {
namespace threading {

ThreadPool::ThreadPool(std::size_t threadCount)
  : m_threads(std::vector<std::thread>(threadCount))
{
  auto threadLoop = [this]()
  {
      while(!m_terminatePool)
      {
          Job job;
          {
              std::unique_lock<std::mutex> lock(m_mutex);

              m_condition.wait(lock, [this]{return !m_jobs.empty() || m_terminatePool; });

              if (!m_terminatePool) {
                job = std::move(m_jobs.front());
                m_jobs.pop();
              }
          }

          if (job) {
            job();
          }
      }
  };

  for(int i = 0; i < threadCount; ++i)
  {
    m_threads[i]  = std::thread(threadLoop);
  }
}

void
ThreadPool::addJob(Job&& job)
{
  {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_jobs.push(std::move(job));
  }
  m_condition.notify_one();
}

ThreadPool::~ThreadPool()
{
  m_terminatePool = true;
  m_condition.notify_all();

  for(auto& thread: m_threads)
  {
    thread.join();
  }
}

}
}
