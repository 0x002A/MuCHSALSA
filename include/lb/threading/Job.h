#pragma once

#include <any>
#include <cstddef>
#include <functional>
#include <vector>

namespace lazybastard::threading {

/**
 * Class representing a job.
 *
 * A Job holds a function to be executed by a thread.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Job {
public:
  /**
   * Class constructor which creates a new instance.
   */
  Job() = default;

  /**
   * Class constructor which creates a new instance.
   *
   * @tparam Types the list of function parameter types
   * @param fn The function to be executed
   * @param params The function parameters
   */
  template <typename... Types> Job(std::function<void(const Job *)> fn, Types... params) : m_fn(fn) {
    addParam(params...);
  };

  /**
   * Copying is disallowed.
   */
  Job(const Job &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Job &operator=(const Job &) = delete;

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
   * @return Whether the Job holds a function
   */
  operator bool() const { return m_fn.operator bool(); };

  /**
   * Function call operator.
   */
  void operator()() const { return m_fn(this); };

  /**
   * Getter for job parameters.
   * This function does not perform any range checking.
   *
   * @param idx The index of the parameter to be returned.
   * @return The parameter value
   */
  std::any getParam(std::size_t idx) const { return m_params[idx]; };

private:
  std::function<void(const Job *)> m_fn; /*!< Function to be executed */
  std::vector<std::any> m_params;        /*!< Function parameters */

  /**
   * Adds a parameter to the internal parameter vector.
   * This function recursively works through the parameter pack.
   *
   * @tparam T the type of the parameter to be added to the vector
   * @tparam Ts the list of the other parameter types
   * @param val The parameter to be added to the vector
   * @param vals The other parameters to be supplied to the next function call
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
