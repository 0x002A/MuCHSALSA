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

#ifndef INCLUDED_LAZYBASTARD_OUTPUTWRITER
#define INCLUDED_LAZYBASTARD_OUTPUTWRITER

#pragma once

#include <cstdio>
#include <memory>
#include <mutex>
#include <string_view>

namespace lazybastard {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// -------------------
// class OutputWriter
// -------------------

/**
 * Class offering **thread-safe** writing to output files.
 */
class OutputWriter {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param fpQuery a std::string_view representing the filepath to the query output file
   * @param fpPaf a std::string_view representing the filepath to the paf output file
   * @param fpTarget a std::string_view representing the filepath to the target output file
   */
  OutputWriter(std::string_view fpQuery, std::string_view fpPaf, std::string_view fpTarget);

  /**
   * Writes data to the query output file.
   * This function is **thread-safe**.
   *
   * @param data a std::string_view representing the data to write
   */
  void writeQuery(std::string_view data);

  /**
   * Writes data to the paf output file.
   * This function is **thread-safe**.
   *
   * @param data a std::string_view representing the data to write
   */
  void writePaf(std::string_view data);

  /**
   * Writes data to the target output file.
   * This function is **thread-safe**.
   *
   * @param data a std::string_view representing the data to write
   */
  void writeTarget(std::string_view data);

private:
  std::unique_ptr<std::FILE, void (*)(std::FILE *)> m_pQueryFile;  /*!< Pointer pointing to the FILE handle required
                                                                    * for accessing the query file
                                                                    */
  std::unique_ptr<std::FILE, void (*)(std::FILE *)> m_pPafFile;    /*!< Pointer pointing to the FILE handle required
                                                                    * for accessing the paf file
                                                                    */
  std::unique_ptr<std::FILE, void (*)(std::FILE *)> m_pTargetFile; /*!< Pointer pointing to the FILE handle required
                                                                    * for accessing the target file
                                                                    */
  //@{
  //** std::mutex for securing the corresponding file handle */
  std::mutex m_mutexQuery, m_mutexPaf, m_mutexTarget;
  //@}
};

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_OUTPUTWRITER

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------