#include "graph/Graph.h"

#include <stdexcept>
#include <string>

#include "Util.h"
#include "graph/Edge.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

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

  return iter != m_vertices.end() ? iter->second->getSharedPtr() : nullptr;
}

void Graph::deleteVertex(const std::string *const pVertexID) {
  for (auto adjacencyIter = m_adjacencyList.begin(); adjacencyIter != m_adjacencyList.end();) {
    if (adjacencyIter->first == *pVertexID) {
      adjacencyIter = m_adjacencyList.erase(adjacencyIter);
    } else {
      auto targetVerticesIter = adjacencyIter->second.find(*pVertexID);
      if (targetVerticesIter != adjacencyIter->second.end()) {
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
  if (vertexIter != m_vertices.end()) {
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

void Graph::deleteEdge(const Edge *const pEdge) {
  auto const vertices = pEdge->getVertices();
  auto const outerIter = m_adjacencyList.find(vertices.first->getID());

  if (outerIter != m_adjacencyList.end()) {
    auto const innerIter = outerIter->second.find(vertices.second->getID());
    if (innerIter != outerIter->second.end()) {
      outerIter->second.erase(innerIter);
    }
  }
}

bool Graph::hasEdge(std::pair<std::string, std::string> &vertexIDs) const {
  lazybastard::util::sortPair(vertexIDs);

  auto const iter = m_adjacencyList.find(vertexIDs.first);

  if (iter != m_adjacencyList.end()) {
    return iter->second.find(vertexIDs.second) != iter->second.end();
  }

  return false;
}

std::unordered_map<std::string const *, Edge const *> Graph::getEdgesOfVertex(std::string const &vertexID) const {
  std::unordered_map<std::string const *, Edge const *> edgeMap;

  auto const outerIter = m_adjacencyList.find(vertexID);
  if (outerIter != m_adjacencyList.end()) {
    for (auto const &[targetID, edge] : outerIter->second) {
      edgeMap.insert({&targetID, edge.get()});
    }
  }

  for (auto const &[targetID, edges] : m_adjacencyList) {
    if (targetID == vertexID) {
      continue;
    }

    auto const interIter = edges.find(vertexID);
    if (interIter != edges.end()) {
      edgeMap.insert({&targetID, interIter->second.get()});
    }
  }

  return edgeMap;
}

void Graph::addEdgeInternal(std::unique_ptr<Edge> &&upEdge) {
  auto assignedVertices = upEdge->getVertices();
  auto iter =
      m_adjacencyList.emplace(assignedVertices.first->getID(), std::unordered_map<std::string, std::unique_ptr<Edge>>())
          .first;
  auto inserted = iter->second.emplace(assignedVertices.second->getID(), std::move(upEdge)).second;

  if (inserted) {
    m_edgeCount += 1;
  }
}

} // namespace lazybastard::graph
