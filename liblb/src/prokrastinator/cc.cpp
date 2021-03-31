#include "Prokrastinator.h"

#include <algorithm>
#include <deque>
#include <set>

#include "graph/Graph.h"

std::vector<std::vector<lazybastard::graph::Vertex *>>
lazybastard::getConnectedComponents(gsl::not_null<const lazybastard::graph::Graph *> pGraph) {
  std::vector<std::vector<lazybastard::graph::Vertex *>> result;
  std::set<lazybastard::graph::Vertex const *> visited;

  std::deque<lazybastard::graph::Vertex *> neighbors;
  auto const vertices = pGraph->getVertices();
  for (auto *const pVertex : vertices) {
    if (visited.contains(pVertex)) {
      continue;
    }

    auto const *const pNeighbors = pGraph->getNeighbors(pVertex);
    std::vector<lazybastard::graph::Vertex *> component({pVertex});
    visited.insert(pVertex);
    if (pNeighbors) {
      component.reserve(pNeighbors->size());
      std::transform(std::begin(*pNeighbors), std::end(*pNeighbors), std::back_inserter(component),
                     [&](auto const &p) { return pGraph->getVertex(p.first); });

      neighbors.push_back(pVertex);

      while (!neighbors.empty()) {
        auto *pCurrentNeighbor = neighbors.front();
        auto const *const pNeighborsOfCurrentNeighbor = pGraph->getNeighbors(pCurrentNeighbor);

        auto iterNeighbor = std::begin(*pNeighborsOfCurrentNeighbor);
        auto biComponent = std::back_inserter(component);
        while (iterNeighbor != std::end(*pNeighborsOfCurrentNeighbor)) {

          pCurrentNeighbor = pGraph->getVertex(iterNeighbor->first);
          if (!visited.contains(pCurrentNeighbor)) {
            if (iterNeighbor->second->getConsensusDirection() != lazybastard::graph::ConsensusDirection::e_NONE) {
              *biComponent++ = pCurrentNeighbor;
              neighbors.push_back(pCurrentNeighbor);
            }
            visited.insert(pCurrentNeighbor);
          }

          ++iterNeighbor;
        }
        neighbors.pop_front();
      }
    }

    result.push_back(std::move(component));
  }

  return result;
}