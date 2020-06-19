#include "graph/Graph.h"

#include <stdexcept>

#include "graph/Vertex.h"
#include "Util.h"

namespace lazybastard {
namespace graph {

void
Graph::addVertex(std::shared_ptr<Vertex>&& spVertex)
{
  std::scoped_lock<std::mutex> lck(m_mutexVertex);
  m_vertices.emplace(spVertex.get()->getID(), std::move(spVertex));
}

std::shared_ptr<Vertex>
Graph::getVertex(const std::string& nanoporeID)
{
  std::scoped_lock<std::mutex> lck(m_mutexVertex);

  auto iter = m_vertices.find(nanoporeID);

  return iter != m_vertices.end() ? iter->second.get()->getSharedPtr() : nullptr;
}

void
Graph::addEdge(const std::pair<std::string, std::string>& vertexIDs)
{
  std::scoped_lock<std::mutex> lck(m_mutexEdge);
  auto pV1 = getVertex(vertexIDs.first);
  auto pV2 = getVertex(vertexIDs.second);

  if (!(pV1 && pV2)) {
    throw std::invalid_argument("Vertices supplied for edge creation aren't allowed to be null.");
  }

  auto vertexPair = std::make_pair(std::move(pV1), std::move(pV2));

  auto upEdge = std::make_unique<Edge>(std::move(lazybastard::util::sortPairOfSharedPointers(vertexPair)));
  addEdgeInternal(std::move(upEdge));
}

void
Graph::addEdgeInternal(std::unique_ptr<Edge>&& upEdge)
{
  auto assignedVertices = upEdge->getVertices();
  auto iter = m_adjacencyList.emplace(
                assignedVertices.first->getID(),
                std::unordered_map<std::string, std::unique_ptr<Edge>>())
                .first;
  auto inserted = iter->second.emplace(assignedVertices.second->getID(), std::move(upEdge)).second;

  if (inserted) {
    m_edgeCount+= 1;
  }
}

}
}
