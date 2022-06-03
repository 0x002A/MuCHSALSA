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

#ifndef INCLUDED_MUCHSALSA_THREADPOOL
#define INCLUDED_MUCHSALSA_THREADPOOL

#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "threading/Job.h"

namespace muchsalsa::threading {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ----------------
// class ThreadPool
// ----------------

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
   * @param threadCount the number of std::thread instances to start
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
   * Adds a Job to the queue of the ThreadPool.
   *
   * @param job an rvalue reference to the Job to be added.
   */
  void addJob(Job &&job);

  /**
   * Adds a Job to the inner queue of the ThreadPool.
   *
   * @param job an rvalue reference to the Job to be added.
   */
  void addSubJob(Job &&job);

  /**
   * Increase count of parents waiting for children.
   */
  void incPassive();

  /**
   * Decrease count of parents waiting for children.
   */
  void decPassive();

private:
  std::vector<std::thread> m_threads;       /*!< std::vector containing the available threads */
  std::queue<Job>          m_jobs;          /*!< std::queue of Jobs to be executed */
  std::condition_variable m_condition;      /*!< std::conditional_variable used to notify threads about new Jobs */
  
  std::vector<std::thread> s_threads;       /*!< std::vector containing the available sub-threads */
  std::queue<Job>          s_jobs;          /*!< std::queue of Jobs to be executed as sub-threads */
  std::condition_variable s_condition;      /*!< std::conditional_variable used to notify sub-threads about new Jobs */
  std::atomic<int> s_passiveJobs{0};        /*!< int indicating how many parent-threads have started sub-threads  */

  std::mutex               ms_mutex;         /*!< std::mutex for securing the parallel use of the std::queue and counters*/
  std::atomic<int> ms_runningJobs{0};        /*!< int indicating how many jobs are actually running  */

  std::atomic<bool> ms_terminatePool{false}; /*!< bool indicating whether the ThreadPool is going to be terminated */

};

} // namespace muchsalsa::threading

#endif // INCLUDED_MUCHSALSA_THREADPOOL

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
