// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of LazyBastardOnMate <https://github.com/0x002A/LazyBastardOnMate>.
//
// LazyBastardOnMate is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// LazyBastardOnMate is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with LazyBastardOnMate.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#ifndef INCLUDED_LAZYBASTARD_JOB
#define INCLUDED_LAZYBASTARD_JOB

#pragma once

#include <any>
#include <cstddef>
#include <functional>
#include <gsl/pointers>
#include <utility>
#include <vector>

namespace lazybastard::threading {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// ---------
// class Job
// ---------

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
   * @tparam TYPES the list of function parameter types
   * @param fn the std::function to be executed
   * @param params the parameters
   */
  template <class... TYPES> explicit Job(std::function<void(gsl::not_null<Job const *> const)> fn, TYPES... params);

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
  explicit operator bool() const;

  /**
   * Function call operator.
   */
  void operator()() const;

  /**
   * Getter for Job parameters.
   * This function does not perform any range checking.
   *
   * @param idx the index of the parameter to be returned
   * @return The parameter value
   */
  [[nodiscard]] std::any getParam(std::size_t idx) const;

private:
  std::function<void(gsl::not_null<Job const *> const)> m_fn;     /*!< std::function to be executed */
  std::vector<std::any>                                 m_params; /*!< std::vector containing the parameters */

  /**
   * Adds a parameter to the internal parameter store.
   * This function recursively works through the parameter pack.
   *
   * @tparam TYPE the type of the parameter to be added to the store
   * @tparam TYPES the list of the other parameter types
   * @param val the parameter to be added to the store
   * @param values the other parameters to be supplied to the next function call
   */
  template <class TYPE, class... TYPES> void _addParam(TYPE val, TYPES... values);

  /**
   * Empty function required to end the recursion.
   */
  void _addParam();
};

// =====================================================================================================================
//                                                  INLINE DEFINITIONS
// =====================================================================================================================

// ---------
// class Job
// ---------

// PUBLIC CLASS METHODS

template<class... TYPES>
Job::Job(std::function<void(gsl::not_null<Job const *> const)> fn, TYPES... params)
        : m_fn(std::move(fn)) {
    _addParam(params...);
}

    inline Job::operator bool() const { return m_fn.operator bool(); }

    inline void Job::operator()() const { return m_fn(this); }

    [[nodiscard]] inline std::any Job::getParam(std::size_t idx) const { return m_params[idx]; }

// PRIVATE CLASS METHODS

    template<class TYPE, class... TYPES>
    void Job::_addParam(TYPE val, TYPES... values) {
        m_params.push_back(val);
        _addParam(values...);
    }

    inline void Job::_addParam() {};

} // namespace lazybastard::threading

#endif // INCLUDED_LAZYBASTARD_JOB

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------