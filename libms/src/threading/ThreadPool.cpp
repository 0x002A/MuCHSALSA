// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of MuCHSALSA <https://github.com/0x002A/MuCHSALSA>.
//
// MuCHSALSA is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// MuCHSALSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with MuCHSALSA.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#include "threading/ThreadPool.h"

#include <utility>

#include "threading/Job.h"

namespace muchsalsa::threading {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ----------------
// class ThreadPool
// ----------------

ThreadPool::ThreadPool(std::size_t threadCount) : m_threads(std::vector<std::thread>(threadCount)), s_threads(std::vector<std::thread>(threadCount)) {

  auto threadLoop = [this]() {
    while (!ms_terminatePool) {
      Job job;
      {
        std::unique_lock<std::mutex> lock(ms_mutex);

        m_condition.wait(lock, [this] { return !m_jobs.empty() || ms_terminatePool; });

        if (!ms_terminatePool) {
          job = std::move(m_jobs.front());
          m_jobs.pop();
        }
        ++ms_runningJobs;
      }

      if (job) {        
        job();
      }
      --ms_runningJobs;
    }
  };

  auto subThreadLoop = [this, threadCount]() {
    while (!ms_terminatePool) {
      Job job;
      {
        std::unique_lock<std::mutex> lock(ms_mutex);

        s_condition.wait(lock, [this, threadCount] { return (!s_jobs.empty() && static_cast<int>(threadCount) + s_passiveJobs > ms_runningJobs) || ms_terminatePool; });

        if (!ms_terminatePool) {
          job = std::move(s_jobs.front());
          s_jobs.pop();
        }
        ++ms_runningJobs;
      }

      if (job) {        
        job();
      }
      --ms_runningJobs;
    }
  };

  for (std::size_t i = 0; i < threadCount; ++i) {
      m_threads[i] = std::thread(threadLoop);
  }
  for (std::size_t i = 0; i < threadCount; ++i) {
      s_threads[i] = std::thread(subThreadLoop);
  }

}

void ThreadPool::addJob(Job &&job) {
  {
    std::scoped_lock<std::mutex> lock(ms_mutex);
    m_jobs.push(std::move(job));
  }
  m_condition.notify_one();
}

void ThreadPool::addSubJob(Job &&job) {
  {
    std::scoped_lock<std::mutex> lock(ms_mutex);
    s_jobs.push(std::move(job));
  }
  s_condition.notify_one();
}

ThreadPool::~ThreadPool() {
  ms_terminatePool = true;
  m_condition.notify_all();
  s_condition.notify_all();

  for (auto &thread : m_threads) {
    thread.join();
  }
  for (auto &thread : s_threads) {
    thread.join();
  }
}

 void ThreadPool::incPassive() {
      ++s_passiveJobs;
      s_condition.notify_one();
 } 

 void ThreadPool::decPassive() {
      --s_passiveJobs;
 } 

} // namespace muchsalsa::threading

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
