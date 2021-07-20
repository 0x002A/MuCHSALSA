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

#ifndef INCLUDED_LAZYBASTARD_LB_FWD
#define INCLUDED_LAZYBASTARD_LB_FWD

#pragma once

namespace lazybastard {

// =====================================================================================================================
//                                                 FORWARD DECLARATIONS
// =====================================================================================================================

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace graph {

class DiGraph;
class Edge;
struct EdgeOrder;
class IGraphObserver;
class Graph;
class GraphBase;
class Vertex;

} // namespace graph

namespace matching {

struct VertexMatch;
struct EdgeMatch;
struct ContainElement;
class MatchMap;
class Id2OverlapMap;

} // namespace matching

namespace threading {

class Job;
class ThreadPool;

} // namespace threading

class BlastFileAccessor;
class BlastFileReader;
class OutputWriter;
class Registry;
class SequenceAccessor;
struct Toggle;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} // namespace lazybastard

#endif // INCLUDED_LAZYBASTARD_LB_FWD
