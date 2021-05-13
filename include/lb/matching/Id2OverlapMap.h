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

#ifndef INCLUDED_LAZYBASTARD_ID2OVERLAPMAP
#define INCLUDED_LAZYBASTARD_ID2OVERLAPMAP

#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>

namespace lazybastard::matching {

// =====================================================================================================================
//                                                        HELPER
// =====================================================================================================================

namespace detail {
using key_t = std::tuple<std::string, std::size_t>;

struct KeyHash : public std::unary_function<key_t, std::size_t> {
  std::size_t operator()(const key_t &key) const {
    return std::hash<std::string>{}(std::get<0>(key)) ^ std::get<1>(key);
  }
};

struct KeyEqual : public std::binary_function<key_t, key_t, bool> {
  bool operator()(const key_t &lhs, const key_t &rhs) const { return lhs == rhs; }
};
} // namespace detail

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// -------------------
// class Id2OverlapMap
// -------------------

class Id2OverlapMap {
public:
  /**
   * Constructor.
   */
  Id2OverlapMap() = default;

  /**
   * Destructor.
   */
  ~Id2OverlapMap() = default;

  /**
   * Moving is disallowed.
   */
  Id2OverlapMap(Id2OverlapMap const &) = delete;

  /**
   * Copying is disallowed.
   */
  Id2OverlapMap(Id2OverlapMap &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Id2OverlapMap &operator=(Id2OverlapMap &&) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Id2OverlapMap &operator=(Id2OverlapMap const &) = delete;

  /**
   * Assigment operator.
   *
   * @param key
   * @return
   */
  std::pair<int, int> &operator[](detail::key_t const &key);

  /**
   * Access operator.
   *
   * @param key
   * @return
   */
  std::pair<int, int> const &at(detail::key_t const &key) const;

  /**
   * Access operator.
   *
   * @param key
   * @return
   */
  std::pair<int, int> &at(detail::key_t const &key);

  std::size_t getSize() const;

private:
  std::unordered_map<detail::key_t, std::pair<int, int>, detail::KeyHash, detail::KeyEqual>
      m_map; /*!< std::unordered_map storing the mapping */
};

inline std::pair<int, int> &Id2OverlapMap::operator[](detail::key_t const &key) { return m_map[key]; }

inline std::pair<int, int> &Id2OverlapMap::at(detail::key_t const &key) { return m_map.at(key); }

inline std::pair<int, int> const &Id2OverlapMap::at(detail::key_t const &key) const { return m_map.at(key); }

inline std::size_t Id2OverlapMap::getSize() const { return m_map.size(); }

} // namespace lazybastard::matching

#endif // INCLUDED_LAZYBASTARD_ID2OVERLAPMAP

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------