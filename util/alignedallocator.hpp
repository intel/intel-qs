// Copyright (C) 2012 Andreas Hehn  <hehn@phys.ethz.ch>.

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef IQS_ALIGNED_ALLOCATOR
#define IQS_ALIGNED_ALLOCATOR

// FIXME: namespace for AlignedAllocator changed from openqu to qhipster

/// \addtogroup sim
/// @{

/// @file alignedallocator.hpp
/// This header defines the class @c AlignedAllocator to provide aligned memory allocation

#include "bitops.hpp"

#ifdef _WIN32
#include <malloc.h>
#else
#include <cstdlib>
#endif

#include <cstddef>
#include <memory>
#include <new>

#if __cplusplus < 201103L
#define noexcept
#endif

namespace qhipster {

/// \class AlignedAllocator
/// @brief An allocator returning aligned memory
///
/// This class provides an aligned C++98 and C++11 conforming allocator.
/// \pre The alignment must be a power of 2.

template <typename T, unsigned int Alignment>
class AlignedAllocator
{
 public:

  typedef T* pointer;
  typedef T const* const_pointer;
  typedef T& reference;
  typedef T const& const_reference;
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  template <typename U>
  struct rebind
  {
    typedef AlignedAllocator<U, Alignment> other;
  };

  AlignedAllocator() noexcept {}
  AlignedAllocator(AlignedAllocator const&) noexcept {}

  template <typename U>
  AlignedAllocator(AlignedAllocator<U, Alignment> const&) noexcept
  {
  }

  pointer allocate(size_type n)
  {
    pointer p;

    static_assert(isPowerOf2(Alignment), "Alignment not a power of 2");

#ifdef _WIN32
    p = reinterpret_cast<pointer>(_aligned_malloc(n * sizeof(T), Alignment));
    if (p == 0) throw std::bad_alloc();
#else
    if (posix_memalign(reinterpret_cast<void**>(&p), Alignment, n * sizeof(T)))
      throw std::bad_alloc();
#endif
    return p;
  }

  void deallocate(pointer p, size_type) noexcept
  {
#ifdef _WIN32
    _aligned_free(p);
#else
    std::free(p);
#endif
  }

  size_type max_size() const noexcept
  {
    std::allocator<T> a;
    return a.max_size();
  }

#if __cplusplus >= 201103L
  template <typename C, class... Args>
  void construct(C* c, Args&&... args)
  {
    new ((void*)c) C(std::forward<Args>(args)...);
  }
#else
  void construct(pointer p, const_reference t) { new ((void*)p) T(t); }
#endif

  template <typename C>
  void destroy(C* c)
  {
    c->~C();
  }

  bool operator==(AlignedAllocator const&) const noexcept { return true; }
  bool operator!=(AlignedAllocator const&) const noexcept { return false; }
  template <typename U, unsigned int UAlignment>
  bool operator==(AlignedAllocator<U, UAlignment> const&) const noexcept
  {
    return false;
  }

  template <typename U, unsigned int UAlignment>
  bool operator!=(AlignedAllocator<U, UAlignment> const&) const noexcept
  {
    return true;
  }
};


}	// namespace qhipster

#if __cplusplus < 201103L
#undef noexcept
#endif

/** @}*/

#endif	// header guard IQS_ALIGNED_ALLOCATOR
