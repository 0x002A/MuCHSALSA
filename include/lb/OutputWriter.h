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

#ifndef INCLUDED_LAZYBASTARD_OUTPUTWRITER
#define INCLUDED_LAZYBASTARD_OUTPUTWRITER

#pragma once

#include <iosfwd>
#include <mutex>

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// -------------------
// class OutputWriter
// -------------------

/**
 * Class offering **thread-safe** writing to output streams.
 */
class OutputWriter {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param osQuery a reference to the query output stream
   * @param osPaf a reference to the paf output stream
   * @param osTarget a reference to the target output stream
   */
  OutputWriter(std::ostream &osQuery, std::ostream &osPaf, std::ostream &osTarget);

private:
  //@{
  //** Reference to a specific output stream */
  std::ostream &m_osQuery, &m_osPaf, &m_osTarget;
  //@}
  //@{
  //** Mutex for securing the corresponding output stream */
  std::mutex m_mutexQuery, m_mutexPaf, m_mutexTarget;
  //@}
};

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_OUTPUTWRITER

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------