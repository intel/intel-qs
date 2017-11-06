// Copyright (C) 2015 Theoretical Physics, ETH Zurich

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

/** @file alignedallocator_test.cpp
 *
 *  Tests for alignedallocator.hpp
 */

#include "alignedallocator.hpp"
#include <cassert>

int main()
{
  // check for a number of repetitions whether alignment is okay
  const unsigned repetitions = 1000;

  using Alloc = openqu::AlignedAllocator<double, 64>;

  Alloc alloc;

  Alloc::pointer pointers[repetitions];

  for (Alloc::pointer& p : pointers) {
    p = alloc.allocate(1);
    assert((std::size_t)(p) % 64 == 0);
  }

  for (Alloc::pointer p : pointers) alloc.deallocate(p, 1);

  return 0;
}
