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

#ifndef INCLUDED_MUCHSALSA_APPLICATION
#define INCLUDED_MUCHSALSA_APPLICATION

#pragma once

#include <cstddef>
#include <gsl/span>
#include <string>

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// -----------------
// class Application
// -----------------

/**
 * Class for parameters.
 *
 * This class holds the input parameters and provides application-specific functions.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Application {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param args the application arguments
   */
  explicit Application(gsl::span<char const *> const &args);

  /**
   * Destructor.
   */
  ~Application() = default;

  /**
   * Copying is disallowed.
   */
  Application(Application const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Application &operator=(Application const &) = delete;

  /**
   * Moving is disallowed.
   */
  Application(Application &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Application &operator=(Application &&) = delete;

  /**
   * Checks if files and directories exist and whether they are readable and writeable.
   */
  [[nodiscard]] auto checkIntegrity() const -> bool;

  /**
   * Getter for the supplied path to the file containing the contigs.
   *
   * @return The supplied path to the file containing the contigs
   */
  [[nodiscard]] std::string const &getContigsFilePath() const noexcept;

  /**
   * Getter for the supplied path to the file containing the unitigs.
   *
   * @return The supplied path to the file containing the unitigs
   */
  [[nodiscard]] std::string const &getUnitigsFilePath() const noexcept;

  /**
   * Getter for the supplied path to the file containing the nanopore data.
   *
   * @return The supplied path to the file containing the nanopore data
   */
  [[nodiscard]] std::string const &getNanoporeFilePath() const noexcept;

  /**
   * Getter for the supplied path to the output directory.
   *
   * @return The supplied path to the output directory
   */
  [[nodiscard]] std::string const &getOutputFilePath() const noexcept;

  /**
   * Getter for the supplied level of parallelism.
   *
   * @return The supplied level of parallelism
   */
  [[nodiscard]] std::size_t const &getThreadCount() const noexcept;

  /**
   * Getter for the supplied wiggle room.
   *
   * @return The supplied wiggle room
   */
  [[nodiscard]] std::size_t const &getWiggleRoom() const noexcept;

private:
  //@{
  //** Filepath */
  std::string m_contigsFilePath, m_unitigsFilePath, m_nanoporeFilePath, m_outputPath;
  //@}
  std::size_t m_threadCount;     /*!< Number of threads to use */
  std::size_t m_wiggleRoom{300}; /*!< Wiggle room */

  /**
   * Parses the parameters from the command line.
   *
   * @param argc the number of parameters.
   * @param argv the array of parameters.
   */
  void _readParameters(gsl::span<char const *> const &args);
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ---------------
// class GraphBase
// ---------------

// PUBLIC CLASS METHODS

[[nodiscard]] inline std::string const &Application::getContigsFilePath() const noexcept { return m_contigsFilePath; }

[[nodiscard]] inline std::string const &Application::getUnitigsFilePath() const noexcept { return m_unitigsFilePath; }

[[nodiscard]] inline std::string const &Application::getNanoporeFilePath() const noexcept { return m_nanoporeFilePath; }

[[nodiscard]] inline std::string const &Application::getOutputFilePath() const noexcept { return m_outputPath; }

[[nodiscard]] inline std::size_t const &Application::getThreadCount() const noexcept { return m_threadCount; }

[[nodiscard]] inline std::size_t const &Application::getWiggleRoom() const noexcept { return m_wiggleRoom; }

#endif // INCLUDED_MUCHSALSA_APPLICATION

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------