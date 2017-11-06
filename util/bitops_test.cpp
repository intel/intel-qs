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

/** @file bitops_test.cpp
 *
 *  Tests for bitops.hpp
 */

#include "bitops.hpp"
#include <cassert>
#include <limits>

int main()
{
  assert(openqu::highestBit(1) == 0);
  assert(openqu::highestBit(2) == 1);
  assert(openqu::highestBit(3) == 1);
  assert(openqu::highestBit(4) == 2);
  assert(openqu::highestBit(std::numeric_limits<int>::max()) == 30);
  assert(openqu::highestBit(std::numeric_limits<unsigned>::max()) == 31);

  assert(openqu::ilog2(1) == 0);
  assert(openqu::ilog2(2) == 1);
  assert(openqu::ilog2(4) == 2);

  assert(!openqu::isPowerOf2(0));
  assert(openqu::isPowerOf2(1));
  assert(openqu::isPowerOf2(2));
  assert(!openqu::isPowerOf2(3));
  assert(openqu::isPowerOf2(4));
  assert(!openqu::isPowerOf2(std::numeric_limits<int>::max()));

  return 0;
}
