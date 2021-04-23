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

#include "graph/Edge.h"
#include "Util.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

// =====================================================================================================================
//                                                     CLASS METHODS
// =====================================================================================================================

Edge::Edge(std::pair<std::shared_ptr<Vertex const> const, std::shared_ptr<Vertex const> const> &&vertices)
    : m_id(Edge::getEdgeId(std::make_pair(vertices.first.get(), vertices.second.get()))),
      m_vertices(std::move(vertices)), m_shadow(false), m_weight(0), m_consensusDirection(Direction::e_NONE) {}

std::string Edge::getEdgeId(std::pair<gsl::not_null<Vertex const *>, gsl::not_null<Vertex const *>> &&vertices) {
  auto id = vertices.first->getId();
  id.append(",");
  id.append(vertices.second->getId());

  return id;
}

} // namespace lazybastard::graph

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------