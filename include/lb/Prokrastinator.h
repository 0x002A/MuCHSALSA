#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <iosfwd>
#include <memory>
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
} // namespace matching

namespace detail {
using Key = std::tuple<std::string, std::size_t>;

struct KeyHash : public std::unary_function<Key, std::size_t> {
  std::size_t operator()(const Key &k) const { return std::hash<std::string>{}(std::get<0>(k)) ^ std::get<1>(k); }
};

struct KeyEqual : public std::binary_function<Key, Key, bool> {
  bool operator()(const Key &v0, const Key &v1) const { return v0 == v1; }
};
} // namespace detail
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

using Id2OverlapMap =
    std::unordered_map<detail::Key, std::tuple<std::size_t, std::size_t>, detail::KeyHash, detail::KeyEqual>;

std::unique_ptr<graph::EdgeOrder> computeOverlap(const lazybastard::matching::MatchMap &matches,
                                                 gsl::not_null<graph::Graph *> pGraph,
                                                 const std::vector<std::string> &ids,
                                                 gsl::not_null<graph::Edge *> pEdge, bool direction, std::size_t score,
                                                 bool isPrimary);

std::vector<std::tuple<std::vector<std::string>, bool>>
getMaxPairwisePaths(const lazybastard::matching::MatchMap &matches, gsl::not_null<graph::Graph *> pGraph,
                    gsl::not_null<graph::Edge *> pEdge, const std::vector<std::string> &illuminaIDs, bool direction);

bool sanityCheck(gsl::not_null<graph::Graph *> pGraph, gsl::not_null<graph::Vertex *> pSubnode,
                 gsl::not_null<graph::Vertex *> pNode, gsl::not_null<graph::Vertex *> pTarget,
                 gsl::not_null<graph::EdgeOrder *> pOrder);

std::unique_ptr<graph::Graph> getMaxSpanTree(gsl::not_null<graph::Graph *> pGraph);

lazybastard::coroutine::generator<std::vector<std::string>>
getConnectedComponents(gsl::not_null<graph::Graph *> pGraph);

std::unique_ptr<graph::DiGraph> getDirectionGraph(const lazybastard::matching::MatchMap &matches,
                                                  gsl::not_null<graph::Graph *> pGraph,
                                                  const std::vector<std::string> &connectedComponent,
                                                  gsl::not_null<graph::Vertex *> pStartNode);

std::vector<std::vector<std::string>> linearizeGraph(gsl::not_null<graph::DiGraph *> pDiGraph);

void assemblePath(const lazybastard::matching::MatchMap &matches, Id2OverlapMap &id2OverlapMap,
                  gsl::not_null<graph::Graph *> pGraph, gsl::not_null<graph::DiGraph *> pDiGraph, std::size_t idx,
                  std::ostream &osQuery, std::ostream &osPAF, std::ostream &osTarget);

} // namespace lazybastard
