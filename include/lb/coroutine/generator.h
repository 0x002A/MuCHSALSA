#pragma once

#include <experimental/coroutine>
#include <memory>

namespace lazybastard::coroutine {

// Based on code found in:
// Rainer Grimm, Concurrency with Modern C++ (Leanpub, 2017 - 2019), 207-209.
template <typename T> class generator {
public:
  struct promise_type;
  using handle_type = std::experimental::coroutine_handle<promise_type>;

  explicit generator(handle_type ht) : coroutine(ht) {}
  generator(const generator &) = delete;
  generator(generator &&other) noexcept : coroutine(other.coroutine) { other.coroutine = nullptr; }
  generator &operator=(const generator &) = delete;
  generator &operator=(generator &&other) noexcept {
    coroutine = other.coroutine;
    other.coroutine = nullptr;
    return *this;
  }
  ~generator() {
    if (coroutine) {
      coroutine.destroy();
    }
  }

  bool next() {
    coroutine.resume();
    return not coroutine.done();
  }

  T getValue() { return coroutine.promise().current_value; }

  struct promise_type {
  public:
    promise_type() = default;
    ~promise_type() = default;
    promise_type(const promise_type &) = delete;
    promise_type(promise_type &&) = delete;
    promise_type &operator=(const promise_type &) = delete;
    promise_type &operator=(promise_type &&) = delete;

    auto initial_suspend() { return std::experimental::suspend_always{}; }

    auto final_suspend() noexcept { return std::experimental::suspend_always{}; }

    auto get_return_object() { return generator{handle_type::from_promise(*this)}; }

    auto return_void() { return std::experimental::suspend_never{}; }

    auto yield_value(T some_value) {
      current_value = some_value;
      return std::experimental::suspend_always{};
    }

    void unhandled_exception() { std::exit(1); }

  private:
    T current_value{};
    friend class generator;
  };

private:
  handle_type coroutine;
};

} // namespace lazybastard::coroutine