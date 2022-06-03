// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2022 Thomas Gatter
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

#include "SequenceAccessor.h"
#include "SequenceUtils.h"
#include "types/Direction.h"
#include "types/Toggle.h"

std::string_view muchsalsa::strSlice(std::string_view original, int intStart, int intEnd) {
  auto const doSlicing = [=](int i, int j) {
    std::size_t intStart = static_cast<std::size_t>(std::max(0, i));
    std::size_t intEnd =
        std::max(std::min(original.size(), static_cast<std::size_t>(std::max(0, j))), static_cast<std::size_t>(i));

    return original.substr(intStart, intEnd - intStart + 1);
  };

  int const size = static_cast<int>(original.size());
  return doSlicing(intStart >= 0 ? intStart : size + intStart, intEnd >= 0 ? intEnd : size + intEnd);
}


std::string muchsalsa::getReverseComplement(std::string const &sequence) {
  std::string reverseComplement;
  reverseComplement.reserve(sequence.size());

  std::transform(std::rbegin(sequence), std::rend(sequence), std::back_inserter(reverseComplement), [](auto const &c) {
    switch (c) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    default:
      return c;
    }
  });

  return reverseComplement;
}

std::string muchsalsa::getIlluminaSequence(muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int illuminaId,
                                int illuminaOverlapLeft, int illuminaOverlapRight, muchsalsa::Toggle const direction) {
  auto illuminaSequence = std::string(
      strSlice(sequenceAccessor.getIlluminaSequence(illuminaId), illuminaOverlapLeft, illuminaOverlapRight + 1));

  if (!direction) {
    return getReverseComplement(illuminaSequence);
  }

  return illuminaSequence;
}

std::string muchsalsa::getNanoporeSequence(muchsalsa::SequenceAccessor &sequenceAccessor, unsigned int nanoporeId,
                                int nanoporeRegionLeft, int nanoporeRegionRight, muchsalsa::Toggle const direction) {
  auto nanoporeSequence = std::string(
      strSlice(sequenceAccessor.getNanoporeSequence(nanoporeId), nanoporeRegionLeft, nanoporeRegionRight + 1));

  if (!direction) {
    return getReverseComplement(nanoporeSequence);
  }

  return nanoporeSequence;
}


// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
