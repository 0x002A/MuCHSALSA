#pragma once

#include <gsl/pointers>
#include <iosfwd>
#include <memory>
#include <tuple>
#include <unordered_map>

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
template <typename T> class Generator;
}

namespace matching {
class MatchMap;
}

namespace detail {
typedef std::tuple<std::string, std::size_t> Key;

struct key_hash : public std::unary_function<Key, std::size_t> {
  std::size_t operator()(const Key &k) const { return std::hash<std::string>{}(std::get<0>(k)) ^ std::get<1>(k); }
};

struct key_equal : public std::binary_function<Key, Key, bool> {
  bool operator()(const Key &v0, const Key &v1) const { return v0 == v1; }
};
} // namespace detail
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

typedef std::unordered_map<detail::Key, std::tuple<std::size_t, std::size_t>, detail::key_hash, detail::key_equal>
    Id2OverlapMap;

std::unique_ptr<graph::EdgeOrder> computeOverlap(const lazybastard::matching::MatchMap &matches,
                                                 const gsl::not_null<graph::Graph *> pGraph,
                                                 const std::vector<std::string> &ids,
                                                 const gsl::not_null<graph::Edge *> pEdge, const bool direction,
                                                 const std::size_t score, const bool isPrimary);

std::vector<std::tuple<std::vector<std::string>, bool>>
getMaxPairwisePaths(const lazybastard::matching::MatchMap &matches, const gsl::not_null<graph::Graph *> pGraph,
                    const gsl::not_null<graph::Edge *> pEdge, const std::vector<std::string> &illuminaIDs,
                    const bool direction);

bool sanityCheck(const gsl::not_null<graph::Graph *> pGraph, const gsl::not_null<graph::Vertex *> pSubnode,
                 const gsl::not_null<graph::Vertex *> pNode, gsl::not_null<graph::Vertex *> pTarget,
                 const gsl::not_null<graph::EdgeOrder *> pOrder);

std::unique_ptr<graph::Graph> getMaxSpanTree(const gsl::not_null<graph::Graph *> pGraph);

lazybastard::coroutine::Generator<std::vector<std::string>>
getConnectedComponents(const gsl::not_null<graph::Graph *> pGraph);

std::unique_ptr<graph::DiGraph> getDirectionGraph(const lazybastard::matching::MatchMap &matches,
                                                  const gsl::not_null<graph::Graph *> pGraph,
                                                  const std::vector<std::string> &connectedComponent,
                                                  const gsl::not_null<graph::Vertex *> pStartNode);

std::vector<std::vector<std::string>> linearizeGraph(const gsl::not_null<graph::DiGraph *> pDiGraph);

void assemblePath(const lazybastard::matching::MatchMap &matches, Id2OverlapMap &id2OverlapMap,
                  const gsl::not_null<graph::Graph *> pGraph, const gsl::not_null<graph::DiGraph *> pDiGraph,
                  std::size_t idx, std::ostream &osQuery, std::ostream &osPAF, std::ostream &osTarget);

} // namespace lazybastard
