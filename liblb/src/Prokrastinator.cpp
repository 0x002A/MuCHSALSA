#include "Prokrastinator.h"

#include <matching/MatchMap.h>
#include <vector>

#include "coroutine/generator.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

using lazybastard::coroutine::generator;

namespace lazybastard {

namespace graph {
class DiGraph {};
} // namespace graph

std::unique_ptr<graph::EdgeOrder> computeOverlap(const lazybastard::matching::MatchMap /*&matches*/,
                                                 const gsl::not_null<graph::Graph *> /*pGraph*/,
                                                 const std::vector<std::string> & /*ids*/,
                                                 const gsl::not_null<graph::Edge *> /*pEdge*/, const bool /*direction*/,
                                                 const std::size_t /*score*/, const bool /*isPrimary*/) {
  return std::make_unique<graph::EdgeOrder>();
}

std::vector<std::tuple<std::vector<std::string>, bool>>
getMaxPairwisePaths(const lazybastard::matching::MatchMap /*&matches*/, const gsl::not_null<graph::Graph *> /*pGraph*/,
                    const gsl::not_null<graph::Edge *> /*pEdge*/, const std::vector<std::string> & /*illuminaIDs*/,
                    const bool /*direction*/) {
  return std::vector<std::tuple<std::vector<std::string>, bool>>();
}

bool sanityCheck(const gsl::not_null<graph::Graph *> /*pGraph*/, const gsl::not_null<graph::Vertex *> /*pSubnode*/,
                 const gsl::not_null<graph::Vertex *> /*pNode*/, const gsl::not_null<graph::Vertex *> /*pTarget*/,
                 const gsl::not_null<graph::EdgeOrder *> /*pOrder*/) {
  return true;
}

std::unique_ptr<graph::Graph> getMaxSpanTree(const gsl::not_null<graph::Graph *> /*pGraph*/) {
  return std::make_unique<graph::Graph>();
}

generator<std::vector<std::string>> getConnectedComponents(const gsl::not_null<graph::Graph *> /*pGraph*/) {
  auto emptyVec = std::vector<std::string>();
  co_yield emptyVec;
}

std::unique_ptr<graph::DiGraph> getDirectionGraph(const lazybastard::matching::MatchMap /*&matches*/,
                                                  const gsl::not_null<graph::Graph *> /*pGraph*/,
                                                  const std::vector<std::string> & /*connectedComponent*/,
                                                  const gsl::not_null<graph::Vertex *> /*pStartNode*/) {
  return std::make_unique<graph::DiGraph>();
}

std::vector<std::vector<std::string>> linearizeGraph(const gsl::not_null<graph::DiGraph *> /*pDiGraph*/) {
  return std::vector<std::vector<std::string>>();
}

void assemblePath(const lazybastard::matching::MatchMap & /*matches*/, const gsl::not_null<graph::Graph *> /*pGraph*/,
                  Id2OverlapMap & /*id2OverlapMap*/, const gsl::not_null<graph::DiGraph *> /*pDiGraph*/,
                  std::size_t /*idx*/, std::ostream & /*osQuery*/, std::ostream & /*osPAF*/,
                  std::ostream & /*osTarget*/) {}

} // namespace lazybastard
