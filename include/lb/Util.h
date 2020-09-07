#pragma once

#include <cassert>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace lazybastard {
namespace util {

template <typename T, typename Enable = void> struct is_valid_pointer : std::false_type {};

template <typename T>
struct is_valid_pointer<T, typename std::enable_if_t<std::is_pointer_v<T> && !std::is_null_pointer_v<T>>>
    : std::true_type {};

template <typename T> inline constexpr bool is_valid_pointer_v = is_valid_pointer<T>::value;

template <typename T, typename std::enable_if_t<!is_valid_pointer_v<T>, int> E = 0> void check_pointers(const T &) {
  static_assert(E != 0, "Unexpected nullptr or value of non-pointer type.");
}

template <typename T, typename std::enable_if_t<is_valid_pointer_v<T>, int> = 0> void check_pointers(const T &p) {
  if (p == nullptr) {
    throw std::invalid_argument("Unexpected nullptr.");
  }
}

template <typename T, typename... Ts,
          typename std::enable_if_t<!std::conjunction_v<is_valid_pointer<T>, is_valid_pointer<Ts>...>, int> E = 0>
void check_pointers(const T &, const Ts &...) {
  static_assert(E != 0, "Unexpected nullptr or value of non-pointer type.");
}

template <typename T, typename... Ts,
          typename std::enable_if_t<std::conjunction_v<is_valid_pointer<T>, is_valid_pointer<Ts>...>, int> = 0>
void check_pointers(const T &p, const Ts &... ps) {
  check_pointers(p);
  check_pointers(ps...);
}

template <typename T, std::enable_if_t<std::is_pointer_v<T>, int> = 0> bool less_than(const T &_1, const T &_2) {
  check_pointers(_1, _2);

  return *_1 < *_2;
}

template <typename T, std::enable_if_t<!std::is_pointer_v<T>, int> = 0> bool less_than(const T &_1, const T &_2) {
  return _1 < _2;
}

template <typename T> std::pair<T, T> &sortPair(std::pair<T, T> &p) {
  if (less_than(p.second, p.first)) {
    auto _ = std::move(p.first);
    p.first = std::move(p.second);
    p.second = std::move(_);
  }

  return p;
}

template <typename T, typename... Args> std::shared_ptr<T> make_shared_aggregate(Args &&... args) {
  return std::make_shared<T>(T{std::forward<Args>(args)...});
}

} // namespace util
} // namespace lazybastard
