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

#include "OutputWriter.h"

#include <gsl/pointers>
#include <stdexcept>

namespace muchsalsa {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

// ------------------
// class OutputWriter
// ------------------

OutputWriter::OutputWriter(std::string_view fpQuery, std::string_view fpPaf, std::string_view fpTarget)
    : m_pQueryFile(decltype(m_pQueryFile)(fopen(fpQuery.data(), "we"),
                                          [](gsl::owner<std::FILE *> pFile) { std::fclose(pFile); })),
      m_pPafFile(
          decltype(m_pPafFile)(fopen(fpPaf.data(), "we"), [](gsl::owner<std::FILE *> pFile) { std::fclose(pFile); })),
      m_pTargetFile(decltype(m_pTargetFile)(fopen(fpTarget.data(), "we"),
                                            [](gsl::owner<std::FILE *> pFile) { std::fclose(pFile); })) {
  if (!m_pQueryFile || !m_pPafFile || !m_pTargetFile) {
    throw std::runtime_error("Can't open output file(s).");
  }
}

void OutputWriter::writeQuery(std::string_view data) {
  std::lock_guard<std::mutex> lck(m_mutexQuery);
  std::fwrite(data.data(), sizeof(data[0]), data.size(), m_pQueryFile.get());
}

void OutputWriter::writePaf(std::string_view data) {
  std::lock_guard<std::mutex> lck(m_mutexPaf);
  std::fwrite(data.data(), sizeof(data[0]), data.size(), m_pPafFile.get());
}

void OutputWriter::writeTarget(std::string_view data) {
  std::lock_guard<std::mutex> lck(m_mutexTarget);
  std::fwrite(data.data(), sizeof(data[0]), data.size(), m_pTargetFile.get());
}

} // namespace muchsalsa

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------