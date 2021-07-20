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

#ifndef INCLUDED_MUCHSALSA_WAITGROUP
#define INCLUDED_MUCHSALSA_WAITGROUP

#pragma once

#include <condition_variable>
#include <cstddef>
#include <mutex>

namespace muchsalsa::threading {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ---------------
// class WaitGroup
// ---------------

/**
 * Class enabling a thread to block until a bunch of Job instances have completed.
 *
 * A WaitGroup holds a count of Job instances which have to complete until it returns and thread-safe methods to
 * increment and decrement this counter. This concept is borrowed from the go programming language.
 */
class WaitGroup {
public:
  /**
   * Class constructor creating a new instance.
   */
  WaitGroup() = default;

  /**
   * Destructor.
   */
  ~WaitGroup() = default;

  /**
   * Copying is disallowed.
   */
  WaitGroup(WaitGroup const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  WaitGroup &operator=(WaitGroup const &) = delete;

  /**
   * Moving is disallowed.
   */
  WaitGroup(WaitGroup &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  WaitGroup &operator=(WaitGroup &&) = delete;

  /**
   * Raises the Job counter by the number supplied as argument.
   *
   * @param newJobs the number to raise the Job counter by
   */
  void add(std::size_t newJobs);

  /**
   * Notifies the WaitGroup that a Job has completed.
   */
  void done();

  /**
   * Blocks the caller until all Job instances assigned to the WaitGroup have completed.
   * On completion the WaitGroup is resetted so it can be reused.
   */
  void wait();

private:
  bool                    m_waitLock{false}; /*!< bool marking the WaitGroup as waiting and locking the add function */
  std::size_t             m_jobCount{0};     /*!< Number of Job instances assigned */
  std::mutex              m_mutex;           /*!< std::mutex used for locking */
  std::condition_variable m_cv;              /*!< std::conditional_variable used to notify threads */
};

} // namespace muchsalsa::threading

#endif // INCLUDED_MUCHSALSA_WAITGROUP

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
