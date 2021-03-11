#include "Prokrastinator.h"

#include <algorithm>
#include <unordered_map>

#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

namespace {
class UnionFind {
public:
  lazybastard::graph::Vertex const *operator[](lazybastard::graph::Vertex const *const pVertex) {
    auto const iter = m_parents.find(pVertex);
    if (iter == std::end(m_parents)) {
      m_parents.insert({pVertex, pVertex});
      m_weights.insert({pVertex, 1});

      return pVertex;
    }

    std::vector<lazybastard::graph::Vertex const *> path({pVertex});
    auto const *pRoot = iter->second;
    while (pRoot != path.back()) {
      path.push_back(pRoot);
      pRoot = m_parents.find(pRoot)->second;
    }

    for (auto const *pAncestor : path) {
      m_parents[pAncestor] = pRoot;
    }

    return pRoot;
  }

  void unify(lazybastard::graph::Vertex const *const pV1, lazybastard::graph::Vertex const *const pV2) {
    auto weightedPair = std::make_pair(operator[](pV1), operator[](pV2));
    if (m_weights[pV2] > m_weights[pV1]) {
      std::swap(weightedPair.first, weightedPair.second);
    }

    m_weights[weightedPair.first] += m_weights[weightedPair.second];
    m_parents[weightedPair.second] = weightedPair.first;
  }

private:
  std::unordered_map<lazybastard::graph::Vertex const *, lazybastard::graph::Vertex const *> m_parents;
  std::unordered_map<lazybastard::graph::Vertex const *, std::size_t> m_weights;
};

std::unordered_map<std::string, std::shared_ptr<lazybastard::graph::Edge>>
kruskal(gsl::not_null<lazybastard::graph::Graph const *> const pGraph) {
  auto edges = pGraph->getEdges();
  edges.erase(std::remove_if(std::begin(edges), std::end(edges),
                             [](auto const *const pEdge) {
                               return pEdge->getConsensusDirection() == lazybastard::graph::ConsensusDirection::e_NONE;
                             }),
              edges.end());
  std::sort(std::begin(edges), std::end(edges),
            [](auto const *pEdge1, auto const *pEdge2) { return pEdge1->getWeight() > pEdge2->getWeight(); });

  std::unordered_map<std::string, std::shared_ptr<lazybastard::graph::Edge>> result;
  auto uf = UnionFind();
  for (auto *pEdge : edges) {
    auto const vertices = pEdge->getVertices();
    if (uf[vertices.first] != uf[vertices.second]) {
      result.insert({pEdge->getID(), pEdge->getSharedPtr()});
      uf.unify(vertices.first, vertices.second);
    }
  }

  return result;
}
} // unnamed namespace

std::unique_ptr<lazybastard::graph::Graph>
lazybastard::getMaxSpanTree(gsl::not_null<lazybastard::graph::Graph const *> const pGraph) {
  auto pNewGraph = std::make_unique<lazybastard::graph::Graph>(*pGraph);
  pNewGraph->replaceEdges(kruskal(pGraph));

  return pNewGraph;
}