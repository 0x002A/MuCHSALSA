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

#ifndef INCLUDED_LAZYBASTARD_EDGE
#define INCLUDED_LAZYBASTARD_EDGE

#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "Lb.fwd.h"
#include "types/Direction.h"
#include "types/Toggle.h"

namespace lazybastard::graph {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ----------------
// struct EdgeOrder
// ----------------

/**
 * Struct representing an order attached to an Edge.
 */
struct EdgeOrder {
  Vertex const *            startVertex; /*!< Start Vertex */
  Vertex const *            endVertex;   /*!< End Vertex */
  double                    leftOffset;  /*!< Left offset */
  double                    rightOffset; /*!< Right offset */
  Toggle                    isContained; /*!< Bool indicating containment */
  Vertex const *            baseVertex;  /*!< Base Vertex */
  std::size_t               score;       /*!< Score */
  std::vector<unsigned int> ids;         /*!< Ids */
  Toggle                    direction;   /*!< Bool indicating direction */
  Toggle                    isPrimary;   /*!< Bool indicating a primary */
};

// ----------
// class Edge
// ----------

/**
 * Class representing an Edge.
 *
 * An Edge is assigned to two instances of Vertex.
 * Instances of this class are **not** thread-safe.
 */
class Edge : public std::enable_shared_from_this<Edge> {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param vertices an rvalue reference to the std::pair of shared pointers pointing to the Vertex instances connected
   *                 by the Edge
   */
  explicit Edge(std::pair<std::shared_ptr<Vertex const> const, std::shared_ptr<Vertex const> const> &&vertices);

  /**
   * Default destructor.
   */
  ~Edge() = default;

  /**
   * Copying is disallowed.
   */
  Edge(Edge const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Edge &operator=(Edge) = delete;

  /**
   * Moving is disallowed.
   */
  Edge(Edge &&) = delete;

  /**
   * Getter returning a std::shared_ptr to this instance of Edge.
   *
   * @return A std::shared_ptr pointing to this instance of Edge
   */
  std::shared_ptr<Edge> getSharedPtr();

  /**
   * Getter returning a std::shared_ptr to this instance of Edge.
   *
   * @return A std::shared_ptr pointing to this instance of Edge
   */
  std::shared_ptr<Edge const> getSharedPtr() const;

  /**
   * Move assignment is disallowed.
   */
  Edge &operator=(Edge &&) = delete;

  /**
   * Getter returning the Vertex instances assigned to this Edge.
   *
   * @return The assigned Vertex instances
   */
  [[nodiscard]] std::pair<Vertex const *, Vertex const *> getVertices() const;

  /**
   * Getter returning the shadow Edge indicator.
   *
   * @return The shadow Edge indicator
   */
  [[nodiscard]] bool isShadow() const;

  /**
   * Setter setting the shadow Edge indicator.
   *
   * @param shadow a bool indicating whether this Edge is a shadow Edge
   */
  void setShadow(bool shadow);

  /**
   * Getter returning the weight of the Edge.
   *
   * @return The weight of the Edge
   */
  [[nodiscard]] std::size_t getWeight() const;

  /**
   * Setter setting the weight of the Edge.
   *
   * @param weight the weight of the Edge
   */
  void setWeight(std::size_t weight);

  /**
   * Getter returning the consensus Direction.
   *
   * @return The consensus Direction
   */
  [[nodiscard]] Direction::Enum getConsensusDirection() const;

  /**
   * Setter setting the consensus Direction.
   *
   * @param consensusDirection the consensus Direction represented by a bool
   */
  void setConsensusDirection(bool consensusDirection);

  /**
   * Getter returning all EdgeOrder elements.
   *
   * @return A constant reference to the std::vector holding the EdgeOrder elements assigned to this Edge
   */
  [[nodiscard]] auto const &getEdgeOrders() const;

  /**
   * Getter returning the EdgeOrder element at the specified index.
   *
   * @param idx the index of the EdgeOrder element to return
   * @return A reference to the requested EdgeOrder element
   * @throws std::out_of_range
   */
  EdgeOrder &orderAt(std::size_t idx);

  /**
   * Appends an EdgeOrder instance to the internal EdgeOrder store.
   *
   * @param edgeOrder a constant reference to an instance of EdgeOrder
   */
  void appendOrder(EdgeOrder const &edgeOrder);

  /**
   * Appends an EdgeOrder instance to the internal EdgeOrder store.
   *
   * @param edgeOrder an rvalue reference to an instance of EdgeOrder
   */
  void appendOrder(EdgeOrder &&edgeOrder);

  /**
   * Replaces the internal EdgeOrder store with the supplied one.
   *
   * @param edgeOrders an rvalue reference to a std::vector of EdgeOrder instances to replace the internal store with
   */
  void replaceOrders(std::vector<EdgeOrder> &&edgeOrders);

  /**
   * Clears the internal store of EdgeOrder instances.
   */
  void clearOrders();

private:
  std::pair<std::shared_ptr<Vertex const> const, std::shared_ptr<Vertex const> const> const
                         m_vertices;           /*!< Assigned Vertex instances */
  std::vector<EdgeOrder> m_orders;             /*!< Assigned EdgeOrder instances */
  bool                   m_shadow;             /*!< Shadow Edge indicator */
  std::size_t            m_weight;             /*!< Weight of the Edge */
  Direction::Enum        m_consensusDirection; /*!< Consensus Direction */
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ----------
// class Edge
// ----------

// PUBLIC CLASS METHODS

inline std::shared_ptr<Edge> Edge::getSharedPtr() { return shared_from_this(); }

inline std::shared_ptr<Edge const> Edge::getSharedPtr() const { return shared_from_this(); }

[[nodiscard]] inline std::pair<Vertex const *, Vertex const *> Edge::getVertices() const {
  return std::make_pair(m_vertices.first.get(), m_vertices.second.get());
}

[[nodiscard]] inline bool Edge::isShadow() const { return m_shadow; }

inline void Edge::setShadow(bool shadow) { m_shadow = shadow; }

[[nodiscard]] inline std::size_t Edge::getWeight() const { return m_weight; }

inline void Edge::setWeight(std::size_t weight) { m_weight = weight; }

[[nodiscard]] inline Direction::Enum Edge::getConsensusDirection() const { return m_consensusDirection; }

inline void Edge::setConsensusDirection(bool consensusDirection) {
  m_consensusDirection = consensusDirection ? Direction::e_POS : Direction::e_NEG;
}

[[nodiscard]] inline auto const &Edge::getEdgeOrders() const { return m_orders; }

inline EdgeOrder &Edge::orderAt(std::size_t idx) { return m_orders.at(idx); }

inline void Edge::appendOrder(EdgeOrder const &edgeOrder) { m_orders.push_back(edgeOrder); }

inline void Edge::appendOrder(EdgeOrder &&edgeOrder) { m_orders.push_back(std::move(edgeOrder)); }

inline void Edge::replaceOrders(std::vector<EdgeOrder> &&edgeOrders) { m_orders = std::move(edgeOrders); }

inline void Edge::clearOrders() { m_orders.clear(); }

} // namespace lazybastard::graph

#endif // INCLUDED_LAZYBASTARD_EDGE

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
