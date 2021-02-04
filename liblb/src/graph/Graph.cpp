#include "graph/Graph.h"

#include <stdexcept>
#include <string>

#include "Util.h"
#include "graph/Edge.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

std::deque<gsl::not_null<std::string const *> const> getShortestPath(
    gsl::not_null<Graph const *> const /*pGraph*/,
    std::pair<gsl::not_null<Vertex const *> const, gsl::not_null<Vertex const *> const> const /*vertexIDs*/) {
  return std::deque<gsl::not_null<std::string const *> const>();
}

Graph::Graph() = default;
Graph::~Graph() = default;

void Graph::addVertex(std::shared_ptr<Vertex> &&spVertex) {
  lazybastard::util::check_pointers(spVertex.get());
  std::unique_lock<std::shared_mutex> lk(m_mutexVertex);

  m_vertices.emplace(spVertex->getID(), std::move(spVertex));
}

auto Graph::getVertex(std::string const &nanoporeID) const {
  std::shared_lock<std::shared_mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeID);

  return iter != std::end(m_vertices) ? iter->second->getSharedPtr() : nullptr;
}

void Graph::deleteVertex(gsl::not_null<std::string const *> const pVertexID) {
  for (auto adjacencyIter = std::begin(m_adjacencyList); adjacencyIter != std::end(m_adjacencyList);) {
    if (adjacencyIter->first == *pVertexID) {
      adjacencyIter = m_adjacencyList.erase(adjacencyIter);
    } else {
      auto targetVerticesIter = adjacencyIter->second.find(*pVertexID);
      if (targetVerticesIter != std::end(adjacencyIter->second)) {
        adjacencyIter->second.erase(targetVerticesIter);
      }
      ++adjacencyIter;
    }
  }

  /*
   * This MUST be done AFTERWARDS as it deletes the map entry of the Vertex and therefore
   * releases the shared_ptr. This probably invalidates all pointers to the ID of the Vertex.
   */
  auto const vertexIter = m_vertices.find(*pVertexID);
  if (vertexIter != std::end(m_vertices)) {
    m_vertices.erase(vertexIter);
  }
}

std::string Graph::addEdge(std::pair<std::string, std::string> const &vertexIDs) {
  std::unique_lock<std::shared_mutex> lck(m_mutexEdge);
  auto pV1 = getVertex(vertexIDs.first);
  auto pV2 = getVertex(vertexIDs.second);

  if (!(pV1 && pV2)) {
    throw std::invalid_argument("Vertices supplied for edge creation aren't allowed to be null.");
  }

  auto edgeID = Edge::getEdgeID(std::make_pair(pV1.get(), pV2.get()));
  auto vertexPair = std::make_pair(std::move(pV1), std::move(pV2));

  auto upEdge = std::make_unique<Edge>(std::move(vertexPair));
  addEdgeInternal(std::move(upEdge));

  return edgeID;
}

Edge const *
Graph::getEdge(std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const {
  lazybastard::util::sortPair(vertexIDs);

  auto const outerIter = m_adjacencyList.find(*vertexIDs.first);
  if (outerIter != std::end(m_adjacencyList)) {
    auto const innerIter = outerIter->second.find(*vertexIDs.second);
    if (innerIter != std::end(outerIter->second)) {
      return innerIter->second.get();
    }
  }

  return nullptr;
}

void Graph::deleteEdge(const Edge *const pEdge) {
  auto const vertices = pEdge->getVertices();
  auto const outerIter = m_adjacencyList.find(vertices.first->getID());

  if (outerIter != std::end(m_adjacencyList)) {
    auto const innerIter = outerIter->second.find(vertices.second->getID());
    if (innerIter != std::end(outerIter->second)) {
      outerIter->second.erase(innerIter);
      --m_edgeCount;
    }
  }
}

bool Graph::hasEdge(
    std::pair<gsl::not_null<std::string const *>, gsl::not_null<std::string const *>> &vertexIDs) const {
  lazybastard::util::sortPair(vertexIDs);

  auto const iter = m_adjacencyList.find(*vertexIDs.first);

  if (iter != std::end(m_adjacencyList)) {
    return iter->second.find(*vertexIDs.second) != std::end(iter->second);
  }

  return false;
}

std::unordered_map<std::string const *, Edge const *> Graph::getEdgesOfVertex(std::string const &vertexID) const {
  std::unordered_map<std::string const *, Edge const *> edgeMap;

  auto const outerIter = m_adjacencyList.find(vertexID);
  if (outerIter != std::end(m_adjacencyList)) {
    for (auto const &[targetID, edge] : outerIter->second) {
      edgeMap.insert({&targetID, edge.get()});
    }
  }

  for (auto const &[targetID, edges] : m_adjacencyList) {
    if (targetID == vertexID) {
      continue;
    }

    auto const interIter = edges.find(vertexID);
    if (interIter != std::end(edges)) {
      edgeMap.insert({&targetID, interIter->second.get()});
    }
  }

  return edgeMap;
}

std::vector<Edge *> Graph::getEdges() const {
  std::vector<Edge *> edges;

  for (const auto &[vertexID, connectedEdges] : m_adjacencyList) {
    LB_UNUSED(vertexID);

    std::transform(std::begin(connectedEdges), std::end(connectedEdges), std::back_inserter(edges),
                   [](const auto &pair) { return pair.second.get(); });
  }

  return edges;
}

void Graph::addEdgeInternal(std::unique_ptr<Edge> &&upEdge) {
  auto assignedVertices = upEdge->getVertices();
  auto iter =
      m_adjacencyList.emplace(assignedVertices.first->getID(), std::unordered_map<std::string, std::unique_ptr<Edge>>())
          .first;
  auto inserted = iter->second.emplace(assignedVertices.second->getID(), std::move(upEdge)).second;

  if (inserted) {
    ++m_edgeCount;
  }
}

} // namespace lazybastard::graph
