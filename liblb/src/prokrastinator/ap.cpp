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

#include "Prokrastinator.h"

#include "OutputWriter.h"
#include "SequenceAccessor.h"
#include "graph/Graph.h"
#include "matching/Id2OverlapMap.h"
#include "matching/MatchMap.h"

void lazybastard::assemblePath(gsl::not_null<graph::Graph const *> const /*pGraph*/,
                               gsl::not_null<lazybastard::matching::MatchMap const *> const /*pMatches*/,
                               gsl::not_null<lazybastard::SequenceAccessor *> const /*pSequenceAccessor*/,
                               gsl::not_null<lazybastard::matching::Id2OverlapMap *> const /*pId2OverlapMap*/,
                               gsl::not_null<std::vector<lazybastard::graph::Vertex const *> const *> const /*pPath*/,
                               gsl::not_null<lazybastard::graph::DiGraph const *> const /*pDiGraph*/,
                               std::size_t /*idx*/, lazybastard::OutputWriter & /*writer*/) {}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------