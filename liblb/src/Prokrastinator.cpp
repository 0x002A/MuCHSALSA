#include "Prokrastinator.h"

#include <vector>

#include "OutputWriter.h"
#include "Util.h"
#include "coroutine/generator.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/ID2OverlapMap.h"
#include "matching/MatchMap.h"

using lazybastard::coroutine::generator;

bool lazybastard::sanityCheck(gsl::not_null<lazybastard::graph::Graph const *> const /*pGraph*/,
                              gsl::not_null<lazybastard::graph::Vertex const *> const /*pSubnode*/,
                              gsl::not_null<lazybastard::graph::Vertex const *> const /*pNode*/,
                              gsl::not_null<lazybastard::graph::Vertex const *> const /*pTarget*/,
                              gsl::not_null<lazybastard::graph::EdgeOrder const *> const /*pOrder*/) {
  return false;
}

std::unique_ptr<lazybastard::graph::Graph>
lazybastard::getMaxSpanTree(gsl::not_null<lazybastard::graph::Graph const *> const /*pGraph*/) {
  return std::make_unique<lazybastard::graph::Graph>();
}

lazybastard::coroutine::generator<std::vector<gsl::not_null<std::string const *>>>
lazybastard::getConnectedComponents(gsl::not_null<const lazybastard::graph::Graph *> /*pGraph*/) {
  auto emptyVec = std::vector<gsl::not_null<std::string const *>>();
  co_yield emptyVec;
}

std::unique_ptr<lazybastard::graph::DiGraph> lazybastard::getDirectionGraph(
    gsl::not_null<graph::Graph const *> const /*pGraph*/,
    gsl::not_null<lazybastard::matching::MatchMap const *> const /*matches*/,
    gsl::not_null<std::vector<gsl::not_null<std::string const *>> const *> const /*connectedComponent*/,
    gsl::not_null<lazybastard::graph::Vertex const *> const /*pStartNode*/) {
  return std::make_unique<lazybastard::graph::DiGraph>();
}

std::vector<std::vector<gsl::not_null<std::string const *>>>
lazybastard::linearizeGraph(gsl::not_null<lazybastard::graph::DiGraph const *> const /*pDiGraph*/) {
  return std::vector<std::vector<gsl::not_null<std::string const *>>>();
}

void lazybastard::assemblePath(gsl::not_null<graph::Graph const *> const /*pGraph*/,
                               gsl::not_null<lazybastard::matching::MatchMap const *> const /*matches*/,
                               gsl::not_null<lazybastard::matching::ID2OverlapMap *> const /*pID2OverlapMap*/,
                               gsl::not_null<std::vector<gsl::not_null<std::string const *>> const *> const /*pPath*/,
                               gsl::not_null<lazybastard::graph::DiGraph const *> const /*pDiGraph*/,
                               std::size_t /*idx*/, lazybastard::OutputWriter & /*writer*/) {}
