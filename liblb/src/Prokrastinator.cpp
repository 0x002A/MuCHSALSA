#include "Prokrastinator.h"

#include <matching/ID2OverlapMap.h>
#include <matching/MatchMap.h>
#include <vector>

#include "Util.h"
#include "coroutine/generator.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

using lazybastard::coroutine::generator;

namespace lazybastard {

graph::EdgeOrder computeOverlap(
    gsl::not_null<matching::MatchMap const *> const /*matches*/, gsl::not_null<graph::Graph const *> const /*pGraph*/,
    std::set<std::string const *const, util::LTCmp<std::string const *const>> const & /*ids*/,
    gsl::not_null<graph::Edge const *> const /*pEdge*/, bool /*direction*/, std::size_t /*score*/, bool /*isPrimary*/) {
  return graph::EdgeOrder{
      nullptr, nullptr, {0.0, 0.0}, {0.0, 0.0},
      false,   nullptr, 0,          std::set<std::string const *const, util::LTCmp<std::string const *const>>(),
      false,   false};
}

std::vector<std::tuple<std::set<std::string const *const, util::LTCmp<std::string const *const>>, std::size_t, bool>>
getMaxPairwisePaths(gsl::not_null<matching::MatchMap const *> const /*matches*/,
                    gsl::not_null<graph::Graph const *> const /*pGraph*/,
                    gsl::not_null<graph::Edge const *> const /*pEdge*/,
                    std::set<std::string const *const, util::LTCmp<std::string const *const>> const & /*illuminaIDs*/,
                    bool /*direction*/) {
  return std::vector<
      std::tuple<std::set<std::string const *const, util::LTCmp<std::string const *const>>, std::size_t, bool>>();
}

bool sanityCheck(gsl::not_null<graph::Graph const *> const /*pGraph*/,
                 gsl::not_null<graph::Vertex const *> const /*pSubnode*/,
                 gsl::not_null<graph::Vertex const *> const /*pNode*/,
                 gsl::not_null<graph::Vertex const *> const /*pTarget*/,
                 gsl::not_null<graph::EdgeOrder const *> const /*pOrder*/) {
  return true;
}

std::unique_ptr<graph::Graph> getMaxSpanTree(gsl::not_null<graph::Graph const *> const /*pGraph*/) {
  return std::make_unique<graph::Graph>();
}

coroutine::generator<std::vector<std::string const *>>
getConnectedComponents(gsl::not_null<const graph::Graph *> /*pGraph*/) {
  auto emptyVec = std::vector<std::string const *>();
  co_yield emptyVec;
}

std::unique_ptr<graph::DiGraph>
getDirectionGraph(gsl::not_null<matching::MatchMap const *> const /*matches*/,
                  gsl::not_null<graph::Graph const *> const /*pGraph*/,
                  gsl::not_null<std::vector<std::string const *> const *> const /*connectedComponent*/,
                  gsl::not_null<graph::Vertex const *> const /*pStartNode*/) {
  return std::make_unique<graph::DiGraph>();
}

std::vector<std::vector<std::string const *>> linearizeGraph(gsl::not_null<graph::DiGraph const *> const /*pDiGraph*/) {
  return std::vector<std::vector<std::string const *>>();
}

void assemblePath(gsl::not_null<matching::MatchMap const *> const /*matches*/,
                  gsl::not_null<graph::Graph const *> const /*pGraph*/,
                  gsl::not_null<matching::ID2OverlapMap *> const /*pID2OverlapMap*/,
                  gsl::not_null<std::vector<std::string const *> const *> const /*pPath*/,
                  gsl::not_null<graph::DiGraph const *> const /*pDiGraph*/, std::size_t /*idx*/,
                  std::ostream & /*osQuery*/, std::ostream & /*osPAF*/, std::ostream & /*osTarget*/) {}

} // namespace lazybastard
