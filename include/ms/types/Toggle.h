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

#ifndef INCLUDED_MUCHSALSA_TOGGLE
#define INCLUDED_MUCHSALSA_TOGGLE

#pragma once

namespace muchsalsa {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ------------
// class Toggle
// ------------

/**
 * Class representing a toggle which can reach two possible states.
 *
 * A toggle has the size of a bool, but supports operations like an integer.
 *
 * Example:
 *
 * Operation: Multiplication
 *
 *    O1    O2    R
 * ___________________
 * |  1  |  1  |  1  |
 * |  1  | -1  | -1  |
 * | -1  | -1  |  1  |
 * | -1  |  1  | -1  |
 * -------------------
 *
 * All Toggle related operations:
 *
 *  Toggle(true)    *= true     ->  true
 *  Toggle(true)    *= false    ->  false
 *  Toggle(false)   *= false    ->  true
 *  Toggle(false)   *= true     ->  false
 *
 *  Toggle(true)    *= Toggle(true)     ->  true
 *  Toggle(true)    *= Toggle(false)    ->  false
 *  Toggle(false)   *= Toggle(false)    ->  true
 *  Toggle(false)   *= Toggle(true)     ->  false
 *
 *  Creation:
 *
 *  Toggle t = true;
 *
 *  Implicit conversion:
 *
 *  auto t = Toggle(false);
 *  t *= false; -> becomes true
 *
 *  if (t) {
 *      // t became true
 *  }
 */
struct Toggle {
  /**
   * Conversion constructor.
   *
   * @param b a bool to convert from
   */
  Toggle(bool b) : m_state(b){}; // NOLINT

  /**
   * Implicit conversion operator.
   *
   * @return The boolean representation of the Toggle
   */
  operator bool() const { return m_state; }; // NOLINT

  Toggle &operator*=(Toggle rhs);
  Toggle &operator*=(bool rhs);
  Toggle &operator&=(Toggle rhs);
  Toggle &operator&=(bool rhs);

  bool operator&&(Toggle rhs) const;
  bool operator&&(bool rhs) const;
  bool operator==(Toggle rhs) const;
  bool operator==(bool rhs) const;
  bool operator!=(Toggle rhs) const;
  bool operator!=(bool rhs) const;

  Toggle operator*(Toggle const &rhs) const;
  Toggle operator*(bool rhs) const;

private:
  bool m_state{false};
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// -------------
// struct Toggle
// -------------

inline Toggle &Toggle::operator*=(Toggle rhs) {
  *this *= rhs.m_state;

  return *this;
}

inline Toggle &Toggle::operator*=(bool rhs) {
  m_state = m_state == rhs;

  return *this;
}

inline Toggle &Toggle::operator&=(Toggle rhs) {
  *this &= rhs.m_state;

  return *this;
}

inline Toggle &Toggle::operator&=(bool rhs) {
  m_state &= rhs;

  return *this;
}

inline bool Toggle::operator&&(Toggle rhs) const { return *this && rhs.m_state; }
inline bool Toggle::operator&&(bool rhs) const { return m_state && rhs; }
inline bool Toggle::operator==(Toggle rhs) const { return *this == rhs.m_state; }
inline bool Toggle::operator==(bool rhs) const { return m_state == rhs; }
inline bool Toggle::operator!=(Toggle rhs) const { return !(*this == rhs); }
inline bool Toggle::operator!=(bool rhs) const { return !(*this == rhs); }

inline Toggle Toggle::operator*(Toggle const &rhs) const { return *this * rhs.m_state; }
inline Toggle Toggle::operator*(bool rhs) const { return Toggle(m_state == rhs); }

// PUBLIC CLASS METHODS

} // namespace muchsalsa

#endif // INCLUDED_MUCHSALSA_TOGGLE

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------