#pragma once

#include <memory>
#include <utility>
#include <stdexcept>

namespace lazybastard {
namespace util {

template<typename T>
std::pair<std::shared_ptr<T>, std::shared_ptr<T>>&
sortPairOfSharedPointers(std::pair<std::shared_ptr<T>, std::shared_ptr<T>>& p)
{
  if (!(p.first && p.second)) {
      throw std::invalid_argument("Sorting requires two valid shared_ptr!");
  }

  if (*p.second.get() < *p.first.get()) {
      auto _ = p.first;
      p.first = std::move(p.second);
      p.second = std::move(_);
  }

  return p;
}

}
}
