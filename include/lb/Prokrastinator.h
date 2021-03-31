#pragma once

#include <cstddef>
#include <deque>
#include <gsl/pointers>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "Lb.fwd.h"
#include "graph/Edge.h"

namespace lazybastard {

std::optional<graph::EdgeOrder> computeOverlap(gsl::not_null<matching::MatchMap const *> matches,
                                               std::deque<std::string> &ids, gsl::not_null<graph::Edge const *> pEdge,
                                               bool direction, std::size_t score, bool isPrimary);

std::vector<std::tuple<std::deque<std::string>, std::size_t, bool>>
getMaxPairwisePaths(gsl::not_null<matching::MatchMap const *> matches, gsl::not_null<graph::Edge const *> pEdge,
                    std::vector<std::string> const &illuminaIDs, Toggle direction, std::size_t wiggleRoom);

bool sanityCheck(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<graph::Vertex const *> pSubnode,
                 gsl::not_null<graph::Vertex const *> pNode, gsl::not_null<graph::Vertex const *> pTarget,
                 gsl::not_null<graph::EdgeOrder const *> pOrder, std::size_t wiggleRoom);

std::unique_ptr<graph::Graph> getMaxSpanTree(gsl::not_null<graph::Graph const *> pGraph);

std::vector<std::vector<lazybastard::graph::Vertex *>>
getConnectedComponents(gsl::not_null<graph::Graph const *> pGraph);

std::unique_ptr<graph::DiGraph> getDirectionGraph(gsl::not_null<graph::Graph const *> pGraph,
                                                  gsl::not_null<graph::Graph const *> pConnectedComponent,
                                                  gsl::not_null<graph::Vertex const *> pStartNode);

std::vector<std::vector<lazybastard::graph::Vertex const *>> linearizeGraph(gsl::not_null<graph::DiGraph *> pDiGraph);

void assemblePath(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<matching::MatchMap const *> matches,
                  gsl::not_null<matching::ID2OverlapMap *> pID2OverlapMap,
                  gsl::not_null<std::vector<lazybastard::graph::Vertex const *> const *> pPath,
                  gsl::not_null<graph::DiGraph const *> pDiGraph, std::size_t idx, OutputWriter &writer);

} // namespace lazybastard
