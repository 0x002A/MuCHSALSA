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

#ifndef INCLUDED_MUCHSALSA_TRACKINGALLOCATOR
#define INCLUDED_MUCHSALSA_TRACKINGALLOCATOR

#pragma once

#include <cstddef>

// =====================================================================================================================
//                                                    FREE FUNCTIONS
// =====================================================================================================================

/**
 * Getter returning the peak memory usage in bytes.
 *
 * @return The peak memory usage in bytes
 */
std::size_t getMemoryUsagePeak();

#endif // INCLUDED_MUCHSALSA_TRACKINGALLOCATOR

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
