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
 *  Tests for bitops.hpp using Boost.Test
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SImUtilBitops
#include <boost/test/unit_test.hpp>

#include <limits>
#include "bitops.hpp"

BOOST_AUTO_TEST_CASE(highestBit)
{
  BOOST_CHECK(openqu::highestBit(1) == 0);
  BOOST_CHECK(openqu::highestBit(2) == 1);
  BOOST_CHECK(openqu::highestBit(3) == 1);
  BOOST_CHECK(openqu::highestBit(4) == 2);
  BOOST_CHECK(openqu::highestBit(std::numeric_limits<int>::max()) == 30);
  BOOST_CHECK(openqu::highestBit(std::numeric_limits<unsigned>::max()) == 31);
}

BOOST_AUTO_TEST_CASE(isPowerOf2)
{
  BOOST_CHECK(!openqu::isPowerOf2(0));
  BOOST_CHECK(openqu::isPowerOf2(1));
  BOOST_CHECK(openqu::isPowerOf2(2));
  BOOST_CHECK(!openqu::isPowerOf2(3));
  BOOST_CHECK(openqu::isPowerOf2(4));
  BOOST_CHECK(!openqu::isPowerOf2(std::numeric_limits<int>::max()));
}
