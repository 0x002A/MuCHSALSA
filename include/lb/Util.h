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

#ifndef INCLUDED_LAZYBASTARD_UTIL
#define INCLUDED_LAZYBASTARD_UTIL

#pragma once

#include <algorithm>
#include <gsl/pointers>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#define LB_UNUSED(x) (void)(x)

namespace lazybastard::util {

// =====================================================================================================================
//                                                       UTILITIES
// =====================================================================================================================

template <class TYPE, class... ARGS> std::shared_ptr<TYPE> make_shared_aggregate(ARGS &&...args) {
  return std::make_shared<TYPE>(TYPE{std::forward<ARGS>(args)...});
}

template <class TYPE> TYPE const *make_const(TYPE *pT) { return static_cast<TYPE const *>(pT); }

template <class TYPE> TYPE const *make_const(TYPE const *pT) { return pT; }

template <class TYPE> gsl::not_null<TYPE const *> make_not_null_and_const(TYPE *pT) {
  return gsl::make_not_null(make_const(pT));
}

template <class TYPE> constexpr void swap_if(TYPE &first, TYPE &second, bool p) {
  if (p) {
    std::swap(first, second);
  }
}

template <class BI_DIR_ITER> constexpr void reverse_if(BI_DIR_ITER first, BI_DIR_ITER last, bool p) {
  if (p) {
    std::reverse(first, last);
  }
}

template <class TYPE, class U = TYPE> constexpr void exchange_if(TYPE &obj, U &&new_value, bool p) {
  if (p) {
    std::exchange(obj, std::forward<U>(new_value));
  }
}

} // namespace lazybastard::util

#endif // INCLUDED_LAZYBASTARD_UTIL

// ---------------------------------------------------- END-OF-FILE ----------------------------------------------------
