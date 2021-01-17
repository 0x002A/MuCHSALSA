#pragma once

#include <any>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lazybastard::graph {

/**
 * Class representing a vertex.
 *
 * A Vertex holds a bunch of meta data and can be assigned to Graph.
 * It can also be connected to instances of Edge.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Vertex : public std::enable_shared_from_this<Vertex> {
public:
  /**
   * Class constructor which creates a new instance.
   *
   * @tparam Types the list of meta datum types
   * @param id The unique id of the vertex
   * @param metaData The meta data of the vertex
   */
  template <typename... Types> explicit Vertex(std::string id, Types... metaData) : m_id(std::move(id)) {
    addMetaDatum(metaData...);
  };

  /**
   * Destructor.
   */
  ~Vertex() = default;

  /**
   * Copying is disallowed.
   */
  Vertex(const Vertex &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Vertex &operator=(const Vertex &) = delete;

  /**
   * Moving is disallowed.
   */
  Vertex(Vertex &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Vertex &operator=(Vertex &&) = delete;

  /**
   * Less than comparison operator.
   *
   * @return Whether the supplied Vertex instance is greater or not
   */
  bool operator<(const Vertex &v) const { return m_id < v.m_id; };

  /**
   * Getter for a shared pointer pointing to this instance of Vertex.
   *
   * @return A new shared pointer pointing to this instance of Vertex
   */
  std::shared_ptr<Vertex> getSharedPtr() { return shared_from_this(); };

  /**
   * Returns a weak pointer pointing to this instance of Vertex.
   *
   * @return A new weak pointer pointing to this instance of Vertex
   */
  std::weak_ptr<Vertex> getWeakPtr() { return weak_from_this(); };

  /**
   * Getter for the unique ID of the Vertex.
   *
   * @return The unique ID of the Vertex
   */
  const std::string &getID() const { return m_id; };

private:
  std::string m_id;                 /*!< Unique vertex id */
  std::vector<std::any> m_metaData; /*!< Vertex's meta data */

  /**
   * Adds a meta datum to the Vertex.
   * This function recursively works through the parameter pack.
   *
   * @tparam T the type of the meta datum to be added to the Vertex
   * @tparam Ts the list of the other meta datum types
   * @param val The meta datum to be added to the Vertex
   * @param vals The other meta data to be supplied to the next function call
   */
  template <typename T, typename... Ts> void addMetaDatum(T val, Ts... vals) {
    m_metaData.push_back(val);
    addMetaDatum(vals...);
  };
  /*
   * Empty function required to end the recursion.
   */
  void addMetaDatum(){};
};

} // namespace lazybastard::graph
