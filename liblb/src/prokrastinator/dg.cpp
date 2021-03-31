#include "Prokrastinator.h"

#include <cassert>
#include <stack>

#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "types/Toggle.h"

std::unique_ptr<lazybastard::graph::DiGraph>
lazybastard::getDirectionGraph(gsl::not_null<graph::Graph const *> const pGraph,
                               gsl::not_null<graph::Graph const *> const pConnectedComponent,
                               gsl::not_null<lazybastard::graph::Vertex const *> const pStartNode) {
  std::stack<std::tuple<lazybastard::graph::Vertex const *, Toggle>> stack;
  stack.push(std::make_tuple(pStartNode, true));

  auto diGraph = std::make_unique<lazybastard::graph::DiGraph>();
  while (!stack.empty()) {
    auto const currentNode = stack.top();
    stack.pop();

    auto const *const pCurrentNode = std::get<0>(currentNode);

    if (!diGraph->hasVertex(pCurrentNode->getID())) {
      diGraph->addVertex(pGraph->getVertexAsSharedPtr(pCurrentNode->getID()));
    }

    if (pCurrentNode->getVertexDirection() == graph::VertexDirection::e_NONE) {
      diGraph->getVertex(pCurrentNode->getID())->setVertexDirection(std::get<1>(currentNode));
    }

    auto const *const pNeighbors = pConnectedComponent->getNeighbors(pCurrentNode);

    if (!pNeighbors) {
      continue;
    }

    for (auto const &[neighborID, pNeighborEdge] : *pNeighbors) {
      auto const *const pOtherNode = pGraph->getVertex(neighborID);
      auto const vertices = pNeighborEdge->getVertices();
      auto otherNodeExists = diGraph->hasVertex(neighborID);

      if (otherNodeExists) {
        otherNodeExists = pOtherNode->getVertexDirection() != graph::VertexDirection::e_NONE;
      }

      if (!otherNodeExists) {
        diGraph->addVertex(pConnectedComponent->getVertexAsSharedPtr(neighborID));
      }

      if (diGraph->hasEdge(std::make_pair(&vertices.first->getID(), &vertices.second->getID())) ||
          diGraph->hasEdge(std::make_pair(&vertices.second->getID(), &vertices.first->getID()))) {
        continue;
      }

      auto const *const pEdge = pGraph->getEdge(std::make_pair(&vertices.first->getID(), &vertices.second->getID()));
      for (auto const &order : pNeighborEdge->getEdgeOrders()) {
        auto flip = false;
        if (!order.direction && order.baseVertex == pOtherNode) {
          flip = !flip;
        }

        if (!std::get<1>(currentNode)) {
          flip = !flip;
        }

        auto const *pStart = order.startVertex;
        auto const *pEnd = order.endVertex;
        if (flip) {
          std::swap(pStart, pEnd);
        }

        auto *pNewEdge = diGraph->getEdge(std::make_pair(&pStart->getID(), &pEnd->getID()));
        if (!pNewEdge) {
          diGraph->addMissingVertex(pGraph->getVertexAsSharedPtr(pStart->getID()));
          diGraph->addMissingVertex(pGraph->getVertexAsSharedPtr(pEnd->getID()));
          diGraph->addEdge(std::make_pair(pStart->getID(), pEnd->getID()));
          pNewEdge = diGraph->getEdge(std::make_pair(&pStart->getID(), &pEnd->getID()));

          pNewEdge->setShadow(pEdge->isShadow());

          if (!pEdge->isShadow()) {
            pNewEdge->setWeight(pNeighborEdge->getWeight());
          }
        }

        pNewEdge->appendOrder(order);
      }

      if (pNeighborEdge->getConsensusDirection() == lazybastard::graph::ConsensusDirection::e_NONE) {
        continue;
      }

      auto const nextMod =
          std::get<1>(currentNode) *
          lazybastard::Toggle(pNeighborEdge->getConsensusDirection() == lazybastard::graph::ConsensusDirection::e_POS);

      if (!otherNodeExists) {
        stack.push(std::make_tuple(pOtherNode, nextMod));
      }
    }
  }

  return diGraph;
}