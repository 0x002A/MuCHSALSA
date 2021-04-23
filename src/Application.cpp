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

#include "Application.h"

#include <filesystem>
#include <stdexcept>
#include <thread>

// =====================================================================================================================
//                                                       CONSTANTS
// =====================================================================================================================

constexpr std::size_t MIN_PAR = 4;

constexpr std::size_t POS_CFP = 1;
constexpr std::size_t POS_UFP = 2;
constexpr std::size_t POS_NFP = 3;
constexpr std::size_t POS_OFP = 4;
constexpr std::size_t POS_NOT = 5;
constexpr std::size_t POS_WGR = 6;

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

Application::Application(gsl::span<char const *> const &args) : m_threadCount(std::thread::hardware_concurrency()) {
  if (args.size() >= MIN_PAR + 1) {
    _readParameters(args);
  } else {
    throw std::runtime_error("Invalid parameter count");
  }
}

auto Application::checkIntegrity() const -> bool {
  std::filesystem::path contigsPath = m_contigsFilePath;
  std::filesystem::path unitigsPath = m_unitigsFilePath;
  std::filesystem::path nanoporePath = m_nanoporeFilePath;
  std::filesystem::path outputPath = m_outputPath;

  auto const exists = std::filesystem::exists(contigsPath) && std::filesystem::exists(unitigsPath) &&
                      std::filesystem::exists(nanoporePath) && std::filesystem::exists(outputPath);

  return exists;
}

void Application::_readParameters(gsl::span<char const *> const &args) {

  m_contigsFilePath  = args[POS_CFP];
  m_unitigsFilePath  = args[POS_UFP];
  m_nanoporeFilePath = args[POS_NFP];
  m_outputPath       = args[POS_OFP];

  // Check for optional parameter
  auto argIdx = args.size() > POS_WGR ? POS_WGR : args.size() - 1;
  switch (argIdx) {
  case POS_WGR:
    m_wiggleRoom = static_cast<std::size_t>(std::stoi(args[POS_WGR]));
    [[fallthrough]];
  case POS_NOT:
    m_threadCount = static_cast<std::size_t>(std::stoi(args[POS_NOT]));
    break;
  }
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
