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

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

bool lazybastard::sanityCheck(gsl::not_null<lazybastard::graph::Graph const *> const     pGraph,
                              gsl::not_null<lazybastard::graph::Vertex const *> const    pSubnode,
                              gsl::not_null<lazybastard::graph::Vertex const *> const    pNode,
                              gsl::not_null<lazybastard::graph::Vertex const *> const    pTarget,
                              gsl::not_null<lazybastard::graph::EdgeOrder const *> const pOrder,
                              std::size_t                                                wiggleRoom) {
  auto const checkOnEdge  = util::make_not_null_and_const(pGraph->getEdge(std::make_pair(pNode, pTarget)));
  auto const checkForEdge = util::make_not_null_and_const(pGraph->getEdge(std::make_pair(pSubnode, pTarget)));
  for (auto const &checkOnOrder : checkOnEdge->getEdgeOrders()) {
    for (auto const &checkForOrder : checkForEdge->getEdgeOrders()) {
      auto isSane = pOrder->direction * checkOnOrder.direction == checkForOrder.direction;

      auto const isContainedCheckFor = checkForOrder.isContained;
      auto const isContainedCheckOn  = checkOnOrder.isContained;

      if (isContainedCheckFor && isContainedCheckOn) {
        isSane &= (checkForOrder.startVertex == pTarget || checkForOrder.endVertex == pTarget) &&
                  checkOnOrder.startVertex == pTarget;
      } else if (isContainedCheckFor && !isContainedCheckOn) {
        if (checkForOrder.endVertex != pTarget) {
          auto l1 = false;
          auto l2 = false;
          auto l3 = false;

          if ((!checkOnOrder.direction && ((pNode == checkOnOrder.baseVertex && checkOnOrder.endVertex == pTarget) ||
                                           (pNode != checkOnOrder.baseVertex && checkOnOrder.endVertex == pTarget))) ||
              (checkOnOrder.direction && checkOnOrder.endVertex == pTarget)) {
            util::exchange_if(l2, true, !pOrder->direction);
          } else {
            l1 = true;
            l3 = true;

            util::exchange_if(l2, true, pOrder->direction);
          }

          util::exchange_if(l1, !l1, !pOrder->direction && pOrder->baseVertex != pOrder->endVertex);
          util::exchange_if(l2, !l2, !checkForOrder.direction && checkForOrder.baseVertex != checkForOrder.endVertex);

          auto const d1 = l1 ? pOrder->leftOffset : pOrder->rightOffset;
          auto const d2 = l2 ? checkForOrder.leftOffset : checkForOrder.rightOffset;
          auto const d3 = l3 ? checkOnOrder.leftOffset : checkOnOrder.rightOffset;

          isSane &= (d1 + d2 + d3) < static_cast<double>(wiggleRoom);
        }
      } else if (!isContainedCheckFor && isContainedCheckOn) {
        isSane &= checkOnOrder.startVertex == pTarget;
      } else {
        auto d1 = checkForOrder.startVertex == pTarget;
        auto d2 = checkOnOrder.startVertex == pTarget;

        util::exchange_if(d1, !d1, !checkForOrder.direction && checkForOrder.baseVertex == pTarget);
        util::exchange_if(d2, !d2, !checkOnOrder.direction && checkOnOrder.baseVertex == pTarget);
        util::exchange_if(d1, !d1, !pOrder->direction);

        isSane &= d1 == d2;
      }

      if (isSane) {
        return true;
      }
    }
  }

  return false;
}

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------