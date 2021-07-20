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

#include "BlastFileAccessor.h"

#include <cstdlib>
#include <gsl/pointers>

namespace lazybastard {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// -----------------------
// class BlastFileAccessor
// -----------------------

// PUBLIC CLASS METHODS

BlastFileAccessor::BlastFileAccessor(std::string_view fpBlastFile)
    : m_pBlastFile(decltype(m_pBlastFile)(fopen(fpBlastFile.data(), "rbe"),
                                          [](gsl::owner<std::FILE *> pFile) { std::fclose(pFile); })) {
  if (!m_pBlastFile) {
    throw std::runtime_error("Can't open blast file.");
  }

  _buildIndex();
}

std::string BlastFileAccessor::getLine(int64_t const &offset) {
  std::scoped_lock<std::mutex> lck(m_mutexBlastFile);

  std::fseek(m_pBlastFile.get(), offset, SEEK_SET);

  char *      pLine      = nullptr;
  std::size_t bufferSize = 0;
  auto        ret        = getline(&pLine, &bufferSize, m_pBlastFile.get());

  auto result = [=]() {
    if (ret != -1 && bufferSize > 0) {
      auto line = std::string(pLine);
      line.pop_back();

      return line;
    }

    return std::string();
  }();

  if (pLine) {
    std::free(pLine); // NOLINT
  }

  return result;
}

// PRIVATE CLASS METHODS

void BlastFileAccessor::_buildIndex() {
  char *      pLine      = nullptr;
  std::size_t bufferSize = 0;
  auto        offset     = std::ftell(m_pBlastFile.get());
  auto        ret        = getline(&pLine, &bufferSize, m_pBlastFile.get());

  while (ret != -1) {
    m_offsets.push_back(offset);

    offset = std::ftell(m_pBlastFile.get());
    ret    = getline(&pLine, &bufferSize, m_pBlastFile.get());
  }

  if (pLine) {
    std::free(pLine); // NOLINT
  }
}

} // namespace lazybastard

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
