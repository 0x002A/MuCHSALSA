#pragma once

#include <string>
#include <any>
#include <vector>
#include <memory>

namespace lazybastard {
namespace graph {

/**
 * Class representing a vertex.
 *
 * A Vertex holds a bunch of meta data and can be assigned to Graph.
 * It can also be part of edges.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Vertex : public std::enable_shared_from_this<Vertex>  {
public:
  /**
   * Class constructor which creates a new instance.
   *
   * @tparam Types the list of meta datum types
   * @param id The unique id of the vertex
   * @param metaData The meta data of the vertex
   */
  template<typename... Types>
  Vertex(std::string id, Types... metaData)
    : m_id(id)
  {
    addMetaDatum(metaData...);
  };

  /**
   * Copying is disallowed.
   */
  Vertex(const Vertex& ) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Vertex& operator=(const Vertex& ) = delete;

  /**
   * Moving is disallowed.
   */
  Vertex(Vertex&& ) = delete;

  /**
   * Move assignment is disallowed.
   */
  Vertex& operator=(Vertex&& ) = delete;

  /**
   * Less than comparison operator.
   */
  bool operator <(const Vertex &v) const { return m_id < v.m_id; };

  /**
   * Returns a new shared pointer pointing to this vertex.
   */
  std::shared_ptr<Vertex> getSharedPtr()
  {
    return shared_from_this();
  };

  /**
   * Returns a new weak pointer pointing to this vertex.
   */
  std::weak_ptr<Vertex> getWeakPtr()
  {
    return weak_from_this();
  };

  /**
   * Returns the unique id of the vertex.
   */
  const std::string& getID() const { return m_id; };
private:
  std::string m_id; /*!< Unique vertex id */
  std::vector<std::any> m_metaData; /*!< Vertex's meta data */

  /**
   * Adds a meta datum to the vertex.
   * This function recursively works through the parameter pack.
   *
   * @tparam T the type of the meta datum to be added to the vertex
   * @tparam Ts the list of the other meta datum types
   * @param val The meta datum to be added to the vertex
   * @param vals The other meta data to be supplied to the next function call
   */
  template<typename T, typename... Ts>
  void addMetaDatum(T val, Ts... vals)
  {
    m_metaData.push_back(val);
    addMetaDatum(vals...);
  };
  /*
   * Empty function required to end the recursion.
   */
  void addMetaDatum() {};
};

}
}
