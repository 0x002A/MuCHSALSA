#pragma once

#include <memory>
#include <string>

#include "Util.h"
#include "graph/Vertex.h"

namespace lazybastard {
namespace graph {

/**
 * Class representing a edge.
 *
 * An Edge is assigned to two instances of Vertex.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Edge  {
public:
  /**
   * Class constructor which creates a new instance.
   *
   * @param vertices pair of shared pointers to the Vertex instaces connected by the Edge
   */
  Edge(std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>>&& vertices)
    : m_vertices(std::move(lazybastard::util::sortPairOfSharedPointers(vertices))) {};

  /**
   * Copying is disallowed.
   */
  Edge(const Edge& ) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Edge& operator=(const Edge& ) = delete;

  /**
   * Moving is disallowed.
   */
  Edge(Edge&& ) = delete;

  /**
   * Move assignment is disallowed.
   */
  Edge& operator=(Edge&& ) = delete;

  /**
   * Returns the vertices assigned to the edge.
   */
  std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> getVertices() const { return m_vertices; };

  /**
   * Generates an edge identifier.
   *
   * @param idV1 identifier of the first vertex
   * @param idV2 identifier of the second vertex
   */
   static std::string getEdgeID(std::string idV1, std::string idV2) { return idV1 + "," + idV1; };
private:
  std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> m_vertices; /*!< Assigned vertices */
};

}
}
