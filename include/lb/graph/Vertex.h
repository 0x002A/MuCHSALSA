#pragma once

#include <any>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace lazybastard::graph {

/**
 * Scoped enum representing the Vertex direction.
 */
struct VertexDirection {
  enum Enum : char { e_POS = 'a', e_NEG = 'b', e_NONE = 'c' };
};

/**
 * Class representing a Vertex.
 *
 * A Vertex holds a bunch of meta data and can be assigned to Graph.
 * It can also be connected to instances of Edge.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Vertex : public std::enable_shared_from_this<Vertex> {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @tparam Types the list of meta datum types
   * @param id a std::string representing the unique ID of the Vertex
   * @param nanoporeLength the nanopore length
   * @param metaData the meta data of the Vertex
   */
  template <typename... Types>
  explicit Vertex(std::string id, std::size_t nanoporeLength, Types... metaData)
      : m_id(std::move(id)), m_nanoporeLength(nanoporeLength) {
    addMetaDatum(metaData...);
  };

  /**
   * Destructor.
   */
  ~Vertex() = default;

  /**
   * Copying is disallowed.
   */
  Vertex(Vertex const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Vertex &operator=(Vertex const &) = delete;

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
   * @param v a constant reference to the Vertex instance to compare
   * @return A bool indicating whether the supplied Vertex instance is greater or not
   */
  bool operator<(Vertex const &v) const { return m_id < v.m_id; };

  /**
   * Getter returning a std::shared_ptr to this instance of Vertex.
   *
   * @return A std::shared_ptr to this instance of Vertex
   */
  std::shared_ptr<Vertex> getSharedPtr() { return shared_from_this(); };

  /**
   * Returns a std::weak_ptr to this instance of Vertex.
   *
   * @return A std::weak_ptr to this instance of Vertex
   */
  std::weak_ptr<Vertex> getWeakPtr() { return weak_from_this(); };

  /**
   * Getter returning the unique ID of this Vertex.
   *
   * @return The unique ID of this Vertex
   */
  auto const &getID() const { return m_id; };

  /**
   * Getter returning the nanopore length.
   *
   * @return The nanopore length
   */
  std::size_t getNanoporeLength() const { return m_nanoporeLength; }

  /**
   * Getter returning the VertexDirection.
   *
   * @return The VertexDirection
   */
  [[nodiscard]] VertexDirection::Enum getVertexDirection() const { return m_direction; };

  /**
   * Setter setting the VertexDirection.
   *
   * @param vertexDirection the VertexDirection
   */
  void setVertexDirection(bool vertexDirection) {
    m_direction = vertexDirection ? VertexDirection::Enum::e_POS : VertexDirection::Enum::e_NEG;
  };

  /**
   * Getter returning the requested meta datum.
   *
   * @tparam T the type of the meta datum
   * @param idx the index of the meta datum to return
   * @return The meta datum
   */
  template <typename T> T getMetaDatum(std::size_t idx) const { return std::any_cast<T>(m_metaData.at(idx)); }

  /**
   * Clears the internal store holding the meta data.
   */
  void clearMetaData() { m_metaData.clear(); };

private:
  std::string const m_id;             /*!< Unique Vertex ID */
  std::size_t const m_nanoporeLength; /*!< Nanopore length*/
  VertexDirection::Enum m_direction;  /*!< Vertex direction */
  std::vector<std::any> m_metaData;   /*!< Vertex's meta data */

  /**
   * Adds a meta datum to the Vertex.
   * This function recursively works through the parameter pack.
   *
   * @tparam T the type of the meta datum to be added to the Vertex
   * @tparam Ts the list of the other meta datum types
   * @param val the meta datum to be added to the Vertex
   * @param vals the other meta data to be supplied to the next function call
   */
  template <typename T, typename... Ts> void addMetaDatum(T val, Ts... vals) {
    m_metaData.push_back(val);
    addMetaDatum(vals...);
  };
  /**
   * Empty function required to end the recursion.
   */
  void addMetaDatum(){};
};

} // namespace lazybastard::graph
