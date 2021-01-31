#include "graph/Edge.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

Edge::Edge(std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> &&vertices)
    : m_id(Edge::getEdgeID(std::make_pair(vertices.first.get(), vertices.second.get()))),
      m_vertices(std::move(lazybastard::util::sortPair(vertices))), m_shadow(false), m_weight(0),
      m_consensusDirection(ConsensusDirection::Enum::e_NONE) {}

std::string Edge::getEdgeID(std::pair<const Vertex *, const Vertex *> &&vertices) {
  lazybastard::util::sortPair(vertices);
  return vertices.first->getID() + "," + vertices.second->getID();
}

} // namespace lazybastard::graph