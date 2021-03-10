#include "Prokrastinator.h"

#include "graph/Graph.h"

std::unique_ptr<lazybastard::graph::Graph>
lazybastard::getMaxSpanTree(gsl::not_null<lazybastard::graph::Graph const *> const /*pGraph*/) {
  return std::make_unique<lazybastard::graph::Graph>();
}