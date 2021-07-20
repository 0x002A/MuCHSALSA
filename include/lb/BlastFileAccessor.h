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

#ifndef INCLUDED_LAZYBASTARD_BLASTFILEACCESSOR
#define INCLUDED_LAZYBASTARD_BLASTFILEACCESSOR

#pragma once

#include <cstdint>
#include <cstdio>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// -----------------------
// class BlastFileAccessor
// -----------------------

/**
 * Class providing access to the BLAST file.
 *
 * The methods of this class are **partially** thread-safe.
 */
class BlastFileAccessor {
public:
  /**
   * Constructor creating a new instance.
   *
   * @param fpBlastFile a std::string_view representing the path to the BLAST file
   */
  BlastFileAccessor(std::string_view fpBlastFile);

  /**
   * Getter returning the std::vector containing the line offsets.
   *
   * @return A std::vector containing the line offsets
   */
  auto const &getLineOffsets() const;

  /**
   * Getter returning the number of lines within the BLAST file.
   *
   * @return The number of lines within the BLAST file
   */
  auto getLineCount() const;

  /**
   * Returns the line starting at the supplied offset.
   * This function is **thread-safe**.
   *
   * @param offset a const reference to an int64_t representing the offset of the line to return
   * @return A std::string representing the line
   */
  std::string getLine(int64_t const &offset);

private:
  std::unique_ptr<std::FILE, void (*)(std::FILE *)>
      m_pBlastFile;                      /*!< Pointer pointing to the FILE handle required for accessing the BLAST file
                                          */
  std::vector<int64_t> m_offsets;        /*!< std::vector containing the line offsets */
  std::mutex           m_mutexBlastFile; /*!< std::mutex for securing the parallel use of the BLAST file */

  /**
   * Builds the line index.
   */
  void _buildIndex();
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// -----------------------
// class BlastFileAccessor
// -----------------------

// PUBLIC CLASS METHODS

inline auto const &BlastFileAccessor::getLineOffsets() const { return m_offsets; }

inline auto BlastFileAccessor::getLineCount() const { return m_offsets.size(); }

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_SEQUENCEACCESSOR

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
