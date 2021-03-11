#pragma once

#include <cstddef>
#include <deque>
#include <gsl/pointers>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Lb.fwd.h"

namespace lazybastard::graph {

/**
 * Struct representing an order attached to an Edge.
 */
struct EdgeOrder {
  Vertex const *startVertex;                                /*!< Start Vertex */
  Vertex const *endVertex;                                  /*!< End Vertex */
  float leftOffset;                                         /*!< Left offset */
  float rightOffset;                                        /*!< Right offset */
  bool isContained;                                         /*!< Bool indicating containment */
  Vertex const *baseVertex;                                 /*!< Base Vertex */
  std::size_t score;                                        /*!< Score */
  std::deque<gsl::not_null<std::string const *> const> ids; /*!< IDs */
  bool direction;                                           /*!< Bool indicating direction */
  bool isPrimary;                                           /*!< Bool indicating a primary */
};

/**
 * Scoped enum representing the consensus direction.
 */
struct ConsensusDirection {
  enum Enum : char { e_POS = 'a', e_NEG = 'b', e_NONE = 'c' };
};

/**
 * Class representing an Edge.
 *
 * An Edge is assigned to two instances of Vertex.
 * Instances of this class are **not** thread-safe.
 */
class Edge : public std::enable_shared_from_this<Edge> {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param vertices an rvalue reference to the std::pair of shared pointers to the Vertex instances connected by
   *                 the Edge
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
   * Getter returning a std::shared_ptr to this instance of Edge.
   *
   * @return A std::shared_ptr to this instance of Edge
   */
  std::shared_ptr<Edge> getSharedPtr() { return shared_from_this(); };

  /**
   * Returns a std::weak_ptr to this instance of Edge.
   *
   * @return A std::weak_ptr to this instance of Edge
   */
  std::weak_ptr<Edge> getWeakPtr() { return weak_from_this(); };

  /**
   * Move assignment is disallowed.
   */
  Edge &operator=(Edge &&) = delete;

  /**
   * Getter returning the unique Edge id.
   *
   * @return The unique Edge id
   */
  [[nodiscard]] auto const &getID() const { return m_id; };

  /**
   * Getter returning the vertices assigned to this Edge.
   *
   * @return The assigned vertices
   */
  [[nodiscard]] std::pair<Vertex const *const, Vertex const *const> getVertices() const {
    return std::make_pair(m_vertices.first.get(), m_vertices.second.get());
  };

  /**
   * Getter returning the shadow Edge indicator.
   *
   * @return The shadow Edge indicator
   */
  [[nodiscard]] bool isShadow() const { return m_shadow; };

  /**
   * Setter setting the shadow Edge indicator.
   *
   * @param shadow the shadow Edge indicator
   */
  void setShadow(bool shadow) { m_shadow = shadow; };

  /**
   * Getter returning the Edge weight.
   *
   * @return The Edge weight
   */
  [[nodiscard]] std::size_t getWeight() const { return m_weight; };

  /**
   * Setter setting the Edge weight.
   *
   * @param weight the Edge weight
   */
  void setWeight(std::size_t weight) { m_weight = weight; };

  /**
   * Getter returning the ConsensusDirection.
   *
   * @return The ConsensusDirection
   */
  [[nodiscard]] ConsensusDirection::Enum getConsensusDirection() const { return m_consensusDirection; };

  /**
   * Setter setting the ConsensusDirection.
   *
   * @param consensusDirection the ConsensusDirection
   */
  void setConsensusDirection(bool consensusDirection) {
    m_consensusDirection = consensusDirection ? ConsensusDirection::Enum::e_POS : ConsensusDirection::Enum::e_NEG;
  };

  /**
   * Getter returning all EdgeOrder elements.
   *
   * @return A constant reference to the std::vector holding the EdgeOrder elements assigned to this Edge
   */
  [[nodiscard]] auto const &getEdgeOrders() const { return m_orders; };

  /**
   * Getter returning the requested EdgeOrder element if it exists.
   *
   * @param idx the index of the EdgeOrder element
   * @return A reference to the requested EdgeOrder element
   * @throws std::out_of_range
   */
  EdgeOrder &orderAt(std::size_t idx) { return m_orders.at(idx); }

  /**
   * Appends an EdgeOrder instance (by copy) to the internal EdgeOrder store.
   *
   * @param edgeOrder a constant reference to an instance of EdgeOrder
   */
  void appendOrder(EdgeOrder const &edgeOrder) { m_orders.push_back(edgeOrder); };

  /**
   * Appends an EdgeOrder instance (by moving) to the internal EdgeOrder store.
   *
   * @param edgeOrder an rvalue reference to an instance of EdgeOrder
   */
  void appendOrder(EdgeOrder &&edgeOrder) { m_orders.push_back(std::move(edgeOrder)); };

  /**
   * Replaces the internal EdgeOrder store with the supplied one.
   *
   * @param edgeOrders an rvalue reference to a std::vector of EdgeOrder instances to replace the internal store with
   *                   (by moving)
   */
  void replaceOrders(std::vector<EdgeOrder> &&edgeOrders) { m_orders = std::move(edgeOrders); };

  /**
   * Clears the internal store of EdgeOrder instances.
   */
  void clearOrders() { m_orders.clear(); };

  /**
   * Generates an Edge identifier based on two Vertex identifiers.
   *
   * @param vertexIDs an rvalue reference to a std::pair of std::string instances representing the Vertex IDs
   * @return The Edge identifier
   */
  static std::string getEdgeID(std::pair<std::string, std::string> &&vertexIDs);

private:
  std::string const m_id; /*!< ID */
  std::pair<std::shared_ptr<Vertex const> const, std::shared_ptr<Vertex const> const> const
      m_vertices;                                /*!< Assigned Vertex instances */
  std::vector<EdgeOrder> m_orders;               /*!< Assigned EdgeOrder instances */
  bool m_shadow;                                 /*!< Is shadow Edge */
  std::size_t m_weight;                          /*!< Edge weight */
  ConsensusDirection::Enum m_consensusDirection; /*!< Consensus direction */
};

} // namespace lazybastard::graph
