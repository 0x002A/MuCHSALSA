#pragma once

#include <cstddef>
#include <functional>
#include <shared_mutex>
#include <string>
#include <tuple>
#include <unordered_map>

namespace lazybastard::matching {

namespace detail {
using Key = std::tuple<std::string, std::size_t>;

struct KeyHash : public std::unary_function<Key, std::size_t> {
  std::size_t operator()(const Key &key) const { return std::hash<std::string>{}(std::get<0>(key)) ^ std::get<1>(key); }
};

struct KeyEqual : public std::binary_function<Key, Key, bool> {
  bool operator()(const Key &_1, const Key &_2) const { return _1 == _2; }
};
} // namespace detail

class ID2OverlapMap {
public:
  /**
   * Constructor.
   */
  ID2OverlapMap() = default;

  /**
   * Destructor.
   */
  ~ID2OverlapMap() = default;

  /**
   * Moving is disallowed.
   */
  ID2OverlapMap(ID2OverlapMap const &) = delete;

  /**
   * Copying is disallowed.
   */
  ID2OverlapMap(ID2OverlapMap &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  ID2OverlapMap &operator=(ID2OverlapMap &&) = delete;

  /**
   * Copy assignment is disallowed.
   */
  ID2OverlapMap &operator=(ID2OverlapMap const &) = delete;

private:
  std::unordered_map<detail::Key, std::tuple<std::size_t, std::size_t>, detail::KeyHash, detail::KeyEqual> m_map;
  mutable std::shared_mutex mutex;
};
} // namespace lazybastard::matching