#pragma once

#include <any>
#include <functional>
#include <gsl/pointers>
#include <utility>
#include <vector>

namespace lazybastard::threading {

/**
 * Class representing a Job.
 *
 * A Job holds a function to be executed by a thread.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Job {
public:
  /**
   * Class constructor creating a new instance.
   */
  Job() = default;

  /**
   * Destructor.
   */
  ~Job() = default;

  /**
   * Class constructor creating a new instance.
   *
   * @tparam Types the list of function parameter types
   * @param fn the std::function to be executed
   * @param params the parameters
   */
  template <typename... Types>
  explicit Job(std::function<void(gsl::not_null<Job const *> const)> fn, Types... params) : m_fn(std::move(fn)) {
    addParam(params...);
  };

  /**
   * Copying is disallowed.
   */
  Job(Job const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Job &operator=(Job const &) = delete;

  /**
   * Move constructor.
   */
  Job(Job &&) = default;

  /**
   * Move assignment operator.
   */
  Job &operator=(Job &&) = default;

  /**
   * Boolean conversion operator.
   *
   * @return A bool indicating whether the Job holds a function
   */
  explicit operator bool() const { return m_fn.operator bool(); };

  /**
   * Function call operator.
   */
  void operator()() const { return m_fn(this); };

  /**
   * Getter for Job parameters.
   * This function does not perform any range checking.
   *
   * @param idx the index of the parameter to be returned.
   * @return The parameter value
   */
  [[nodiscard]] std::any getParam(std::size_t idx) const { return m_params[idx]; };

private:
  std::function<void(gsl::not_null<Job const *> const)> m_fn; /*!< std::function to be executed */
  std::vector<std::any> m_params;                             /*!< Parameters */

  /**
   * Adds a parameter to the internal parameter store.
   * This function recursively works through the parameter pack.
   *
   * @tparam T the type of the parameter to be added to the store
   * @tparam Ts the list of the other parameter types
   * @param val the parameter to be added to the store
   * @param vals the other parameters to be supplied to the next function call
   */
  template <typename T, typename... Ts> void addParam(T val, Ts... vals) {
    m_params.push_back(val);
    addParam(vals...);
  };
  /*
   * Empty function required to end the recursion.
   */
  void addParam(){};
};

} // namespace lazybastard::threading
