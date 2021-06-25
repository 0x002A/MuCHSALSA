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

#include "threading/WaitGroup.h"

#include <stdexcept>

namespace lazybastard::threading {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ---------------
// class WaitGroup
// ---------------

void WaitGroup::add(std::size_t newJobs) {
  {
    std::lock_guard<std::mutex> lck(m_mutex);

    if (m_waitLock) {
      throw std::logic_error("Adding jobs to an already waiting wait group is illegal behaviour.");
    }

    m_jobCount += newJobs;
  }

  m_cv.notify_one();
}

void WaitGroup::done() {
  {
    std::lock_guard<std::mutex> lck(m_mutex);

    if (m_jobCount > 0) {
      m_jobCount -= 1;
    }
  }

  m_cv.notify_one();
}

void WaitGroup::wait() {
  {
    std::lock_guard<std::mutex> lck(m_mutex);

    m_waitLock = true;
  }

  std::unique_lock<std::mutex> lck(m_mutex);
  m_cv.wait(lck, [this] { return m_jobCount == 0; });
  m_waitLock = false;
}

} // namespace lazybastard::threading

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
