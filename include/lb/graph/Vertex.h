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

#ifndef INCLUDED_LAZYBASTARD_VERTEX
#define INCLUDED_LAZYBASTARD_VERTEX

#pragma once

#include <any>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "types/Direction.h"

namespace lazybastard::graph {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ------------
// class Vertex
// ------------

/**
 * Class representing a Vertex.
 *
 * A Vertex holds a bunch of meta data and can be assigned to a Graph or DiGraph.
 * It can also be connected to instances of Edge.
 */
class Vertex : public std::enable_shared_from_this<Vertex> {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @tparam TYPES the list of meta datum types
   * @param id a std::string representing the unique Id of the Vertex
   * @param nanoporeLength the nanopore length
   * @param metaData the meta data of the Vertex
   */
  template <class... TYPES> explicit Vertex(std::string id, std::size_t nanoporeLength, TYPES... metaData);

  /**
   * Destructor.
   */
  ~Vertex() = default;

  /**
   * Copying is disallowed.
   */
  Vertex(Vertex const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Vertex &operator=(Vertex const &) = delete;

  /**
   * Moving is disallowed.
   */
  Vertex(Vertex &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Vertex &operator=(Vertex &&) = delete;

  /**
   * Less than comparison operator.
   *
   * @param rhs a constant reference to the Vertex instance to compare
   * @return A bool indicating whether the supplied Vertex instance is greater or not
   */
  auto operator<(Vertex const &rhs) const;

  /**
   * Getter returning a std::shared_ptr to this instance of Vertex.
   *
   * @return A std::shared_ptr pointing to this instance of Vertex
   */
  auto getSharedPtr();

  /**
   * Getter returning a std::shared_ptr to this instance of Vertex.
   *
   * @return A std::shared_ptr pointing to this instance of Vertex
   */
  auto getSharedPtr() const;

  /**
   * Getter returning the unique Id of this Vertex.
   *
   * @return The unique Id of this Vertex
   */
  auto const &getId() const;

  /**
   * Getter returning the nanopore length.
   *
   * @return The nanopore length
   */
  auto getNanoporeLength() const;

  /**
   * Getter returning the Vertex Direction.
   *
   * @return The Vertex Direction
   */
  [[nodiscard]] Direction::Enum getVertexDirection() const;

  /**
   * Setter setting the Vertex Direction.
   *
   * @param vertexDirection the Vertex Direction represented by a bool
   */
  void setVertexDirection(bool vertexDirection);

  /**
   * Getter returning the meta datum at the specified index.
   *
   * @tparam TYPE the type of the meta datum
   * @param idx the index of the meta datum to return
   * @return The meta datum
   * @throws std::out_of_range
   */
  template <class TYPE> TYPE getMetaDatum(std::size_t idx) const;

  /**
   * Clears the internal store holding the meta data.
   */
  void clearMetaData();

private:
  std::string const     m_id;             /*!< Unique Vertex Id */
  std::size_t const     m_nanoporeLength; /*!< Nanopore length */
  Direction::Enum       m_direction;      /*!< Vertex Direction */
  std::vector<std::any> m_metaData;       /*!< Meta data */

  /**
   * Adds a meta datum to the Vertex.
   * This function recursively works through the parameter pack.
   *
   * @tparam TYPE the type of the current meta datum
   * @tparam TYPES the list of the other meta datum types
   * @param val the current meta datum
   * @param values the other meta data to be supplied to the next function call
   */
  template <class TYPE, class... TYPES> void _addMetaDatum(TYPE val, TYPES... values);

  /**
   * Empty function required to end the recursion.
   */
  void _addMetaDatum();
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ------------
// class Vertex
// ------------

// PUBLIC CLASS METHODS

template <class... TYPES> Vertex::Vertex(std::string id, std::size_t nanoporeLength, TYPES... metaData)
    : m_id(std::move(id)), m_nanoporeLength(nanoporeLength), m_direction(Direction::e_NONE) {
  _addMetaDatum(metaData...);
}

inline auto Vertex::operator<(Vertex const &rhs) const { return m_id < rhs.m_id; }

inline auto Vertex::getSharedPtr() { return shared_from_this(); }

inline auto Vertex::getSharedPtr() const { return shared_from_this(); }

inline auto const &Vertex::getId() const { return m_id; }

inline auto Vertex::getNanoporeLength() const { return m_nanoporeLength; }

[[nodiscard]] inline Direction::Enum Vertex::getVertexDirection() const { return m_direction; }

inline void Vertex::setVertexDirection(bool vertexDirection) {
  m_direction = vertexDirection ? Direction::e_POS : Direction::e_NEG;
}

template <class TYPE> TYPE Vertex::getMetaDatum(std::size_t idx) const {
  return std::any_cast<TYPE>(m_metaData.at(idx));
}

inline void Vertex::clearMetaData() { m_metaData.clear(); }

// PRIVATE CLASS METHODS

template <class TYPE, class... TYPES> void Vertex::_addMetaDatum(TYPE val, TYPES... values) {
  m_metaData.push_back(val);
  _addMetaDatum(values...);
}

inline void Vertex::_addMetaDatum() {}

} // namespace lazybastard::graph

#endif // INCLUDED_LAZYBASTARD_VERTEX

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------