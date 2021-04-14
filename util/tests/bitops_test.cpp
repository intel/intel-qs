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

// FIXME: namespace changed from openqu to iqs
// FIXME: most of the code is not necessary

/// @file bitops_test.cpp
/// Tests for bitops.hpp

#include "bitops.hpp"

#include <cassert>
#include <limits>

int main()
{
  assert(iqs::highestBit(1) == 0);
  assert(iqs::highestBit(2) == 1);
  assert(iqs::highestBit(3) == 1);
  assert(iqs::highestBit(4) == 2);
  assert(iqs::highestBit(std::numeric_limits<int>::max()) == 30);
  assert(iqs::highestBit(std::numeric_limits<unsigned>::max()) == 31);

  assert(iqs::ilog2(1) == 0);
  assert(iqs::ilog2(2) == 1);
  assert(iqs::ilog2(4) == 2);

  assert(!iqs::isPowerOf2(0));
  assert(iqs::isPowerOf2(1));
  assert(iqs::isPowerOf2(2));
  assert(!iqs::isPowerOf2(3));
  assert(iqs::isPowerOf2(4));
  assert(!iqs::isPowerOf2(std::numeric_limits<int>::max()));

  return 0;
}
