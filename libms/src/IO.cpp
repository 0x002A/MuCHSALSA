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

#include "IO.h"

#include <cstring>

namespace muchsalsa {

// =====================================================================================================================
//                                                        HELPER
// =====================================================================================================================

namespace {

void *realloc(void *pOldMem, std::size_t sizeOld, std::size_t sizeNew) {
  auto *pNewMem = ::operator new(sizeNew);

  if (!pOldMem || !pNewMem) {
    return nullptr;
  }

  std::memcpy(pNewMem, pOldMem, sizeOld);

  ::operator delete(pOldMem);

  return pNewMem;
}

} // unnamed namespace

// =====================================================================================================================
//                                                    FREE FUNCTIONS
// =====================================================================================================================

int64_t readline(char **pPtrBuffer, std::size_t *pSizeBuffer, std::FILE *pFile) {
  if (*pPtrBuffer == nullptr || *pSizeBuffer == 0) {
    *pSizeBuffer = 120;

    *pPtrBuffer = static_cast<char *>(::operator new(*pSizeBuffer));
    if (*pPtrBuffer == nullptr) {
      return -1;
    }
  }

  auto *pCurrentBufferPos = *pPtrBuffer;
  auto *pLastBufferPos    = *pPtrBuffer + *pSizeBuffer - 1; // NOLINT
  auto  c                 = std::fgetc(pFile);
  while (c != -1) {
    *pCurrentBufferPos++ = static_cast<char>(c); // NOLINT
    if (c == '\n') {
      *pCurrentBufferPos = '\0';
      return pCurrentBufferPos - *pPtrBuffer;
    }
    if (pCurrentBufferPos + 1 > pLastBufferPos) { // NOLINT
      std::size_t sizeNewBuffer = *pSizeBuffer * 2;
      auto        charsInBuffer = static_cast<std::size_t>(pCurrentBufferPos - *pPtrBuffer);

      auto *const pNewBuffer = static_cast<char *>(realloc(*pPtrBuffer, charsInBuffer, sizeNewBuffer));
      if (!pNewBuffer) {
        return -1;
      }

      *pPtrBuffer       = pNewBuffer;
      *pSizeBuffer      = sizeNewBuffer;
      pLastBufferPos    = pNewBuffer + sizeNewBuffer; // NOLINT
      pCurrentBufferPos = pNewBuffer + charsInBuffer; // NOLINT
    }

    c = std::fgetc(pFile);
  }

  if (std::feof(pFile) && pCurrentBufferPos != *pPtrBuffer) {
    *pCurrentBufferPos = '\0';
    return pCurrentBufferPos - *pPtrBuffer;
  }

  return -1;
}

} // namespace muchsalsa

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
