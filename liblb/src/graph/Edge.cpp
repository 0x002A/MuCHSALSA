#include "graph/Edge.h"
#include "Util.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

Edge::Edge(std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> &&vertices)
    : m_id(Edge::getEdgeID(std::make_pair(vertices.first->getID(), vertices.second->getID()))),
      m_vertices(std::move(vertices)), m_shadow(false), m_weight(0),
      m_consensusDirection(ConsensusDirection::Enum::e_NONE) {}

std::string Edge::getEdgeID(std::pair<std::string, std::string> &&vertexIDs) {
  return std::move(vertexIDs.first) + "," + std::move(vertexIDs.second);
}

} // namespace lazybastard::graph