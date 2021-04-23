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

#ifndef INCLUDED_LAZYBASTARD_DIRECTION
#define INCLUDED_LAZYBASTARD_DIRECTION

#pragma once

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ----------------
// struct Direction
// ----------------

/**
 * Scoped enum representing a Direction.
 */
struct Direction {
  enum Enum : char { e_POS = 'a', e_NEG = 'b', e_NONE = 'c' };
};

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_DIRECTION

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
