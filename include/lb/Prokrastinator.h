#pragma once

#include <cstddef>
#include <deque>
#include <gsl/pointers>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Lb.fwd.h"
#include "graph/Edge.h"

namespace lazybastard {

graph::EdgeOrder computeOverlap(gsl::not_null<graph::Graph const *> pGraph,
                                gsl::not_null<matching::MatchMap const *> matches,
                                std::deque<gsl::not_null<std::string const *> const> const &ids,
                                gsl::not_null<graph::Edge const *> pEdge, bool direction, std::size_t score,
                                bool isPrimary);

std::vector<std::tuple<std::deque<gsl::not_null<std::string const *> const>, std::size_t, bool>> getMaxPairwisePaths(
    gsl::not_null<matching::MatchMap const *> matches, gsl::not_null<graph::Edge const *> pEdge,
    std::set<gsl::not_null<std::string const *> const, util::LTCmp<gsl::not_null<std::string const *> const>> const
        &illuminaIDs,
    bool direction, std::size_t wiggleRoom);

bool sanityCheck(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<graph::Vertex const *> pSubnode,
                 gsl::not_null<graph::Vertex const *> pNode, gsl::not_null<graph::Vertex const *> pTarget,
                 gsl::not_null<graph::EdgeOrder const *> pOrder);

std::unique_ptr<graph::Graph> getMaxSpanTree(gsl::not_null<graph::Graph const *> pGraph);

coroutine::generator<std::vector<gsl::not_null<std::string const *>>>
getConnectedComponents(gsl::not_null<graph::Graph const *> pGraph);

std::unique_ptr<graph::DiGraph>
getDirectionGraph(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<matching::MatchMap const *> matches,
                  gsl::not_null<std::vector<gsl::not_null<std::string const *>> const *> connectedComponent,
                  gsl::not_null<graph::Vertex const *> pStartNode);

std::vector<std::vector<gsl::not_null<std::string const *>>>
linearizeGraph(gsl::not_null<graph::DiGraph const *> pDiGraph);

void assemblePath(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<matching::MatchMap const *> matches,
                  gsl::not_null<matching::ID2OverlapMap *> pID2OverlapMap,
                  gsl::not_null<std::vector<gsl::not_null<std::string const *>> const *> pPath,
                  gsl::not_null<graph::DiGraph const *> pDiGraph, std::size_t idx, OutputWriter &writer);

} // namespace lazybastard
