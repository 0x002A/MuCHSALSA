#pragma once

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Util.h"
#include "graph/Vertex.h"

namespace lazybastard::graph {

struct EdgeOrder {
  Vertex *startVertex;
  Vertex *endVertex;
  std::pair<float, float> leftOffset;
  std::pair<float, float> rightOffset;
  bool contained;
  Vertex *baseVertex;
  size_t score;
  std::tuple<std::string, std::string> ids;
  bool direction;
  bool isPrimary;
};

/**
 * Class representing an Edge.
 *
 * An Edge is assigned to two instances of Vertex.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Edge {
public:
  /**
   * Class constructor which creates a new instance.
   *
   * @param vertices pair of shared pointers to the Vertex instances connected by
   * the Edge
   */
  Edge(std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> &&vertices)
      : m_vertices(std::move(lazybastard::util::sortPair(vertices))){};

  /**
   * Copying is disallowed.
   */
  Edge(const Edge &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Edge &operator=(const Edge &) = delete;

  /**
   * Moving is disallowed.
   */
  Edge(Edge &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Edge &operator=(Edge &&) = delete;

  /**
   * Getter for the assigned vertices.
   *
   * @return The assigned vertices
   */
  std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> getVertices() const { return m_vertices; };

  /**
   * Getter for the specified EdgeOrder element with bounds checking.
   *
   * @param idx index of the EdgeOrder element
   * @return Reference to the requested EdgeOrder element
   * @throws std::out_of_range
   */
  EdgeOrder& orderAt(std::size_t idx)
  {
    return m_orders.at(idx);
  }

  /**
   * Appends an EdgeOrder instance to the privately held vector.
   *
   * @param edgeOrder instance of EdgeOrder to copy
   */
  void appendOrder(const EdgeOrder &edgeOrder)
  {
    m_orders.push_back(edgeOrder);
  };

  /**
   * Appends an EdgeOrder instance to the privately held vector.
   *
   * @param edgeOrder instance of EdgeOrder to move
   */
  void appendOrder(EdgeOrder &&edgeOrder)
  {
    m_orders.push_back(std::move(edgeOrder));
  };

  /**
   * Replaces the privately held vector with the supplied one.
   *
   * @param vEdgeOrders the vector of EdgeOrder instances to replace (by moving) the internal vector with
   */
  void replaceOrders(std::vector<EdgeOrder> &&vEdgeOrders)
  {
    m_orders = std::move(vEdgeOrders);
  };

  /**
   * Clears the privately held vector of EdgeOrder instances.
   */
  void clearOrders()
  {
    m_orders.clear();
  };

  /**
   * Generates identifier based on two identifiers of a Vertex.
   *
   * @param idV1 identifier of the first Vertex
   * @param idV2 identifier of the second Vertex
   * @return The identifier
   */
  static std::string getEdgeID(std::string idV1, std::string idV2) { return idV1 + "," + idV1; };

private:
  std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> m_vertices; /*!< Assigned Vertex instances */
  std::vector<EdgeOrder> m_orders; /*!< Assigned EdgeOrder instances */
};

} // namespace lazybastard::graph
