// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of LazyBastardOnMate <https://github.com/0x002A/LazyBastardOnMate>.
//
// LazyBastardOnMate is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// LazyBastardOnMate is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with LazyBastardOnMate.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#include "threading/ThreadPool.h"

#include <utility>

#include "threading/Job.h"

namespace lazybastard::threading {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ----------------
// class ThreadPool
// ----------------

ThreadPool::ThreadPool(std::size_t threadCount) : m_threads(std::vector<std::thread>(threadCount)) {
  auto threadLoop = [this]() {
    while (!m_terminatePool) {
      Job job;
      {
        std::unique_lock<std::mutex> lock(m_mutex);

        m_condition.wait(lock, [this] { return !m_jobs.empty() || m_terminatePool; });

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

  for (std::size_t i = 0; i < threadCount; ++i) {
    m_threads[i] = std::thread(threadLoop);
  }
}

void ThreadPool::addJob(Job &&job) {
  {
    std::scoped_lock<std::mutex> lock(m_mutex);
    m_jobs.push(std::move(job));
  }
  m_condition.notify_one();
}

ThreadPool::~ThreadPool() {
  m_terminatePool = true;
  m_condition.notify_all();

  for (auto &thread : m_threads) {
    thread.join();
  }
}

} // namespace lazybastard::threading

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
