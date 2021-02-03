#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <iosfwd>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace graph {
class DiGraph;
class Edge;
class EdgeOrder;
class Graph;
class Vertex;
} // namespace graph

namespace coroutine {
template <typename T> class generator;
} // namespace coroutine

namespace matching {
class MatchMap;
class ID2OverlapMap;
} // namespace matching

namespace util {
template <typename T> struct LTCmp;
} // namespace util
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

graph::EdgeOrder computeOverlap(gsl::not_null<matching::MatchMap const *> matches,
                                gsl::not_null<graph::Graph const *> pGraph,
                                std::set<std::string const *const, util::LTCmp<std::string const *const>> const &ids,
                                gsl::not_null<graph::Edge const *> pEdge, bool direction, std::size_t score,
                                bool isPrimary);

std::vector<std::tuple<std::set<std::string const *const, util::LTCmp<std::string const *const>>, std::size_t, bool>>
getMaxPairwisePaths(gsl::not_null<matching::MatchMap const *> matches, gsl::not_null<graph::Graph const *> pGraph,
                    gsl::not_null<graph::Edge const *> pEdge,
                    std::set<std::string const *const, util::LTCmp<std::string const *const>> const &illuminaIDs,
                    bool direction);

bool sanityCheck(gsl::not_null<graph::Graph const *> pGraph, gsl::not_null<graph::Vertex const *> pSubnode,
                 gsl::not_null<graph::Vertex const *> pNode, gsl::not_null<graph::Vertex const *> pTarget,
                 gsl::not_null<graph::EdgeOrder const *> pOrder);

std::unique_ptr<graph::Graph> getMaxSpanTree(gsl::not_null<graph::Graph const *> pGraph);

coroutine::generator<std::vector<std::string const *>>
getConnectedComponents(gsl::not_null<graph::Graph const *> pGraph);

std::unique_ptr<graph::DiGraph>
getDirectionGraph(gsl::not_null<matching::MatchMap const *> matches, gsl::not_null<graph::Graph const *> pGraph,
                  gsl::not_null<std::vector<std::string const *> const *> connectedComponent,
                  gsl::not_null<graph::Vertex const *> pStartNode);

std::vector<std::vector<std::string const *>> linearizeGraph(gsl::not_null<graph::DiGraph const *> pDiGraph);

void assemblePath(gsl::not_null<matching::MatchMap const *> matches, gsl::not_null<graph::Graph const *> pGraph,
                  gsl::not_null<matching::ID2OverlapMap *> pID2OverlapMap,
                  gsl::not_null<std::vector<std::string const *> const *> pPath,
                  gsl::not_null<graph::DiGraph const *> pDiGraph, std::size_t idx, std::ostream &osQuery,
                  std::ostream &osPAF, std::ostream &osTarget);

} // namespace lazybastard
