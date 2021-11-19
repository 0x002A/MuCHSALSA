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

#ifndef INCLUDED_MUCHSALSA_ID2OVERLAPMAP
#define INCLUDED_MUCHSALSA_ID2OVERLAPMAP

#pragma once

#include <cstddef>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace muchsalsa::matching {

// =====================================================================================================================
//                                                        HELPER
// =====================================================================================================================

namespace detail {
using key_t = std::tuple<unsigned int, std::size_t>;

struct KeyHash : public std::unary_function<key_t, std::size_t> {
  std::size_t operator()(const key_t &key) const {
    return std::hash<unsigned int>{}(std::get<0>(key)) ^ std::get<1>(key);
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

/**
 * Class representing the mapping of pairs of a illumina id and a clique index to overlaps.
 */
class Id2OverlapMap {
public:
  /**
   * Assigment operator.
   *
   * @param key a const reference to a std::pair of an unsigned int and a std::size_t representing the illumina id and
   *            the clique index which the overlap is or should become mapped to
   * @return A reference to the mapped overlap, performing an insertion if no element with a key equivalent to key exist
   */
  std::pair<int, int> &operator[](detail::key_t const &key);

  /**
   * Access operator.
   *
   * @param key a const reference to a std::pair of an unsigned int and a std::size_t representing the illumina id and
   *            the clique index which the overlap is mapped to
   * @return A const reference to the mapped overlap
   * @throws std::out_of_range
   */
  std::pair<int, int> const &at(detail::key_t const &key) const;

  /**
   * Access operator.
   *
   * @param key a const reference to a std::pair of an unsigned int and a std::size_t representing the illumina id and
   *            the clique index which the overlap is mapped to
   * @return A reference to the mapped overlap
   * @throws std::out_of_range
   */
  std::pair<int, int> &at(detail::key_t const &key);

private:
  std::unordered_map<detail::key_t, std::pair<int, int>, detail::KeyHash, detail::KeyEqual>
      m_map; /*!< std::unordered_map storing the mapping */
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// -------------------
// class Id2OverlapMap
// -------------------

// PUBLIC CLASS METHODS

inline std::pair<int, int> &Id2OverlapMap::operator[](detail::key_t const &key) { return m_map[key]; }

inline std::pair<int, int> &Id2OverlapMap::at(detail::key_t const &key) { return m_map.at(key); }

inline std::pair<int, int> const &Id2OverlapMap::at(detail::key_t const &key) const { return m_map.at(key); }

} // namespace muchsalsa::matching

#endif // INCLUDED_MUCHSALSA_ID2OVERLAPMAP

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------