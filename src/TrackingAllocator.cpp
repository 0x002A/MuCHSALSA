// -*- C++ -*-
//===---------------------------------------------------------------------------------------------------------------==//
//
// Copyright (C) 2021 Kevin Klein
// This file is part of MuCHSALSA <https://github.com/0x002A/MuCHSALSA>.
//
// MuCHSALSA is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// MuCHSALSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with MuCHSALSA.
// If not, see <http://www.gnu.org/licenses/>.
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
//===---------------------------------------------------------------------------------------------------------------==//

#include "TrackingAllocator.h"

#include <cstdlib>
#include <mutex>
#include <new>

namespace {

// =====================================================================================================================
//                                                         TYPES
// =====================================================================================================================

// -----------------------
// struct AllocationHeader
// -----------------------

/**
 * Struct representing the header placed at the beginning of a reserved zone.
 */
struct AllocationHeader {
  std::size_t size; /*!< Size of the memory chunk */
};

// -----------------------
// class TrackingAllocator
// -----------------------

/**
 * Class implementing an allocator capable of tracking allocations and deallocations.
 * This allocator can be queried for the max number of bytes in use at the same time.
 *
 * The methods of this class are thread-safe.
 */
class TrackingAllocator {
public:
  /**
   * Allocates a block of memory of the requested size that is suitably aligned for any scalar type.
   * The underlying memory chunk looks like this:
   * H H R R U U U U U U
   *   R -- reserved zone (size = __STDCPP_DEFAULT_NEW_ALIGNMENT__).
   *   H -- AllocationHeader (size = sizeof(std::size_t)), placed at the beginning of the reserved zone.
   *   U -- user memory.
   *
   * @param count a std::size_t representing the number of bytes to be allocated
   * @return A pointer to the beginning of the user memory (nullptr if the allocation fails)
   */
  void *allocate(std::size_t count) noexcept {
    std::lock_guard<std::mutex> lck{m_mutex};

    if (count == 0) {
      count = 1;
    }

    auto const sizeReservedZone = __STDCPP_DEFAULT_NEW_ALIGNMENT__;
    auto const sizeToAllocate   = sizeReservedZone + count;

    void *const pMem = std::malloc(sizeToAllocate); // NOLINT

    if (!pMem) {
      return nullptr;
    }

    m_memUsage += count;
    if (m_memUsage > m_memUsagePeak) {
      m_memUsagePeak = m_memUsage;
    }

    auto *const pMemStart = static_cast<std::byte *>(pMem);
    auto *const pMemUser  = static_cast<void *>(pMemStart + sizeReservedZone); // NOLINT

    auto *const pHeader = static_cast<AllocationHeader *>(pMem);
    pHeader->size       = count;

    return pMemUser;
  }

  /**
   * Frees a block of memory previously allocated using this allocator.
   *
   * @param ptr a pointer to the beginning of the user memory of the memory block to be freed
   */
  void deallocate(void *ptr) noexcept {
    std::lock_guard<std::mutex> lck{m_mutex};

    if (!ptr) {
      return;
    }

    auto *const       pMemUser  = static_cast<std::byte *>(ptr);
    auto *const       pMemStart = static_cast<void *>(pMemUser - __STDCPP_DEFAULT_NEW_ALIGNMENT__); // NOLINT
    auto const *const pHeader   = static_cast<AllocationHeader *>(pMemStart);

    m_memUsage -= pHeader->size;

    std::free(pMemStart); // NOLINT
  }

  /**
   * Getter returning the peak memory usage in bytes.
   *
   * @return The peak memory usage in bytes
   */
  [[nodiscard]] std::size_t getMemoryUsagePeak() const { return m_memUsagePeak; }

private:
  std::size_t m_memUsage{0};     /*!< Current memory in use */
  std::size_t m_memUsagePeak{0}; /*!< Peak of memory in use */
  std::mutex  m_mutex;           /*!< std::mutex for securing the parallel use of the allocator */
} allocator;                     // NOLINT
} // unnamed namespace

// =====================================================================================================================
//                                                    FREE FUNCTIONS
// =====================================================================================================================

std::size_t getMemoryUsagePeak() { return allocator.getMemoryUsagePeak(); }

// =====================================================================================================================
//                                                 OVERLOADED OPERATORS
// =====================================================================================================================

#ifdef TRACK_MEMORY_USAGE
void *operator new(std::size_t count) { return allocator.allocate(count); }
void  operator delete(void *ptr) noexcept { allocator.deallocate(ptr); }
#endif
