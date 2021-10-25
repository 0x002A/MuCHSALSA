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

#ifndef INCLUDED_MUCHSALSA_IO
#define INCLUDED_MUCHSALSA_IO

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>

namespace muchsalsa {

// =====================================================================================================================
//                                                    FREE FUNCTIONS
// =====================================================================================================================

/**
 * Reads an entire line from a file, storing the address of the buffer containing the text using the supplied pointer.
 * The buffer will be null-terminated and will contain the newline character, if one was found.
 *
 * @param pPtrBuffer a pointer to a pointer receiving the buffer address
 * @param pSizeBuffer a pointer to a std::size_t representing the size of the supplied buffer
 * @param pFile a pointer to a FILE handle required for accessing the file
 * @return The number of characters read or -1 in the event of a failure
 */
int64_t readline(char **pPtrBuffer, std::size_t *pSizeBuffer, std::FILE *pFile);

} // namespace muchsalsa

#endif // INCLUDED_MUCHSALSA_IO
