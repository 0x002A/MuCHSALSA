#include "Prokrastinator.h"

#include <stack>

#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"

std::unique_ptr<lazybastard::graph::DiGraph>
lazybastard::getDirectionGraph(gsl::not_null<graph::Graph const *> const pGraph,
                               gsl::not_null<graph::Graph const *> const pConnectedComponent,
                               gsl::not_null<lazybastard::graph::Vertex *> const pStartNode) {
  std::stack<std::tuple<lazybastard::graph::Vertex const *, bool>> stack;
  stack.push(std::make_tuple(pStartNode, true));

  auto diGraph = std::make_unique<lazybastard::graph::DiGraph>();
  while (!stack.empty()) {
    auto const currentNode = stack.top();
    stack.pop();

    if (!diGraph->hasVertex(std::get<0>(currentNode)->getID())) {
      diGraph->addVertex(pGraph->getVertexAsSharedPtr(std::get<0>(currentNode)->getID()));
    }

    auto const edges = pConnectedComponent->getEdges();
    for (auto const *const pCCEdge : edges) {
      auto const vertices = pCCEdge->getVertices();
      auto const *const otherNode = vertices.second;
      auto otherNodeExists = diGraph->hasVertex(otherNode->getID());

      if (!otherNodeExists) {
        diGraph->addVertex(pConnectedComponent->getVertexAsSharedPtr(otherNode->getID()));
      }

      if (diGraph->hasEdge(std::make_pair(&vertices.first->getID(), &vertices.second->getID())) ||
          diGraph->hasEdge(std::make_pair(&vertices.second->getID(), &vertices.first->getID()))) {
        continue;
      }

      auto const *const pEdge = pGraph->getEdge(std::make_pair(&vertices.first->getID(), &vertices.second->getID()));
      for (auto const &order : pCCEdge->getEdgeOrders()) {
        auto flip = false;
        if (!order.direction && order.baseVertex == otherNode) {
          flip = true;
        }

        if (!std::get<1>(currentNode)) {
          flip = false;
        }

        auto const *pStart = order.startVertex;
        auto const *pEnd = order.endVertex;
        if (flip) {
          std::swap(pStart, pEnd);
        }

        auto *pNewEdge = diGraph->getEdge(std::make_pair(&pStart->getID(), &pEnd->getID()));
        if (!pNewEdge) {
          diGraph->addEdge(std::make_pair(pStart->getID(), pEnd->getID()));
          pNewEdge = diGraph->getEdge(std::make_pair(&pStart->getID(), &pEnd->getID()));

          pNewEdge->setShadow(pEdge->isShadow());

          if (!pEdge->isShadow()) {
            pNewEdge->setWeight(pCCEdge->getWeight());
          }
        }

        pNewEdge->appendOrder(order);
      }

      if (pCCEdge->getConsensusDirection() == lazybastard::graph::ConsensusDirection::e_NONE) {
        continue;
      }

      auto const nextMod = static_cast<bool>(
          lazybastard::Toggle(std::get<1>(currentNode)) &&
          lazybastard::Toggle(pCCEdge->getConsensusDirection() == lazybastard::graph::ConsensusDirection::e_POS));

      if (!otherNodeExists) {
        stack.push(std::make_tuple(otherNode, nextMod));
      }
    }
  }

  return diGraph;
}