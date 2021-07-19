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

#ifndef INCLUDED_LAZYBASTARD_REGISTRY
#define INCLUDED_LAZYBASTARD_REGISTRY

#pragma once

#include <mutex>
#include <string>
#include <unordered_map>

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// --------------
// class Registry
// --------------

/**
 * Class enabling the replacement of large identifiers with unsigned ints.
 *
 * The methods of this class are **partially** thread-safe.
 */
class Registry {
public:
  /**
   * Returns the number assigned to the supplied identifier.
   * Assigns the next free number beforehand if identifier is unknown.
   * This function is **thread-safe**.
   *
   * @param str a const reference to a std::string representing the identifier
   * @return An unsigned int assigned to the suppplied identifier
   */
  const unsigned int &operator[](std::string const &str);

  /**
   * Clears the identifier mapping.
   */
  void clear() noexcept;

private:
  unsigned int m_nextValue{0}; /*!< unsigned int representing the next free number */
  std::unordered_map<std::string, unsigned int>
             m_identifiers; /*!< std::unordered_map mapping the identifiers to numbers */
  std::mutex m_mutex;       /*!< std::mutex for securing the parallel use of the std::unordered_map */
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// --------------
// class Registry
// --------------

// PUBLIC CLASS METHODS

inline void Registry::clear() noexcept {
  m_identifiers.clear();
  m_nextValue = 0;
}

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_REGISTRY

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
