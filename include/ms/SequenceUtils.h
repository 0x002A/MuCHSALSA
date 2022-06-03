// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Thomas Gatter
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

#ifndef INCLUDED_MUCHSALSA_SEQUENCEUTILS
#define INCLUDED_MUCHSALSA_SEQUENCEUTILS

#pragma once

#include <string>
#include <string_view>


namespace muchsalsa {

// =====================================================================================================================
//                                                        SequenceUtils
// =====================================================================================================================

std::string_view strSlice(std::string_view original, int intStart, int intEnd);

std::string getReverseComplement(std::string const &sequence);

std::string getIlluminaSequence(muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int illuminaId,
                                int illuminaOverlapLeft, int illuminaOverlapRight, muchsalsa::Toggle const direction);

std::string getNanoporeSequence(muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int nanoporeId,
                                int nanoporeRegionLeft, int nanoporeRegionRight, muchsalsa::Toggle const direction);

}
#endif // INCLUDED_MUCHSALSA_SEQUENCEUTILS

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
