#pragma once

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Util.h"

namespace lazybastard::graph {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Vertex;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

struct EdgeOrder {
  Vertex const *const startVertex;
  Vertex const *const endVertex;
  std::pair<float const, float const> const leftOffset;
  std::pair<float const, float const> const rightOffset;
  bool const isContained;
  Vertex const *const baseVertex;
  std::size_t const score;
  std::set<std::string const *const, util::LTCmp<std::string const *const>> const ids;
  bool const direction;
  bool const isPrimary;
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
  explicit Edge(std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> &&vertices);

  /**
   * Default destructor.
   */
  ~Edge() = default;

  /**
   * Copying is disallowed.
   */
  Edge(Edge const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Edge &operator=(Edge const &) = delete;

  /**
   * Moving is disallowed.
   */
  Edge(Edge &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Edge &operator=(Edge &&) = delete;

  /**
   * Getter for the unique id of the Edge.
   *
   * @return The unique id of the Edge
   */
  [[nodiscard]] auto const &getID() const { return m_id; };

  /**
   * Getter for the assigned vertices.
   *
   * @return The assigned vertices
   */
  [[nodiscard]] std::pair<Vertex const *const, Vertex const *const> getVertices() const {
    return std::make_pair(m_vertices.first.get(), m_vertices.second.get());
  };

  /**
   * Getter for the shadow Edge indicator.
   *
   * @return Shadow Edge indicator
   */
  [[nodiscard]] bool isShadow() const { return m_shadow; };

  /**
   * Setter for the shadow Edge indicator.
   *
   * @param shadow shadow Edge indicator
   */
  void setShadow(bool shadow) { m_shadow = shadow; };

  /**
   * Getter for the Edge weight.
   *
   * @return Edge weight
   */
  [[nodiscard]] std::size_t getWeight() const { return m_weight; };

  /**
   * Setter for the Edge weight.
   *
   * @param weight Edge weight
   */
  void setWeight(std::size_t weight) { m_weight = weight; };

  /**
   * Getter for the consensus direction.
   *
   * @return Consensus direction
   */
  [[nodiscard]] bool getConsensusDirection() const { return m_consensusDirection; };

  /**
   * Setter for the consensus direction.
   *
   * @param direction consensus direction
   */
  void setConsensusDirection(bool consensusDirection) { m_consensusDirection = consensusDirection; };

  /**
   * Getter for all EdgeOrder elements.
   *
   * @return Reference to the vector holding the EdgeOrder elements
   */
  [[nodiscard]] auto const &getEdgeOrders() const { return m_orders; };

  /**
   * Getter for the specified EdgeOrder element with bounds checking.
   *
   * @param idx index of the EdgeOrder element
   * @return Reference to the requested EdgeOrder element
   * @throws std::out_of_range
   */
  EdgeOrder &orderAt(std::size_t idx) { return m_orders.at(idx); }

  /**
   * Appends an EdgeOrder instance to the privately held vector.
   *
   * @param edgeOrder instance of EdgeOrder to copy
   */
  void appendOrder(EdgeOrder const &edgeOrder) { m_orders.push_back(edgeOrder); };

  /**
   * Appends an EdgeOrder instance to the privately held vector.
   *
   * @param edgeOrder instance of EdgeOrder to move
   */
  void appendOrder(EdgeOrder &&edgeOrder) { m_orders.push_back(std::move(edgeOrder)); };

  /**
   * Replaces the privately held vector with the supplied one.
   *
   * @param vEdgeOrders the vector of EdgeOrder instances to replace (by moving) the internal vector with
   */
  void replaceOrders(std::vector<EdgeOrder> &&vEdgeOrders) { m_orders = std::move(vEdgeOrders); };

  /**
   * Clears the privately held vector of EdgeOrder instances.
   */
  void clearOrders() { m_orders.clear(); };

  /**
   * Generates identifier based on two identifiers of a Vertex.
   *
   * @param vertices pair of shared_ptr to Vertex
   * @return The identifier
   */
  static std::string getEdgeID(std::pair<Vertex *, Vertex *> &&vertices);

private:
  std::string const m_id; /*!< ID */
  std::pair<std::shared_ptr<Vertex const> const, std::shared_ptr<Vertex const> const> const
      m_vertices;                  /*!< Assigned Vertex instances */
  std::vector<EdgeOrder> m_orders; /*!< Assigned EdgeOrder instances */
  bool m_shadow{false};            /*!< Is shadow Edge */
  std::size_t m_weight;            /*!< Edge weight */
  bool m_consensusDirection;       /*!< Consensus direction */
};

} // namespace lazybastard::graph
