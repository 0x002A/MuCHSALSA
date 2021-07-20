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

#include "Registry.h"

namespace lazybastard {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// --------------
// class Registry
// --------------

// PUBLIC CLASS METHODS

const unsigned int &Registry::operator[](std::string const &str) {
  std::scoped_lock<std::mutex> lck(m_mutex);

  if (!m_identifiers.contains(str)) {
    m_identifiers[str] = m_nextValue;
    ++m_nextValue;
  }

  return m_identifiers.at(str);
}

} // namespace lazybastard

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
