// Copyright (C) 2015 Theoretical Physics, ETH Zurich
//
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

#ifndef IQS_BITOPS_HPP
#define IQS_BITOPS_HPP

/// \addtogroup util
/// @{

/// @file bitops.hpp
///
/// This header defines useful bit manipilation functions, especially determining the
/// highest bit set to 1 and checking whether a number is a power of 2.

#include <cassert>
#include <cstdint>

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////

namespace detail {

template <class Integral>
constexpr unsigned highestBitImpl(Integral i, unsigned pos)
{
  return i == 1 ? pos : highestBitImpl(i >> 1, pos + 1);
}

}	// namespace detail

/////////////////////////////////////////////////////////////////////////////////////////

template <class Integral>
unsigned floor_power_of_two(Integral x)
{
  int power = 1;
  while (x >>= 1) power <<= 1;

  return power;
}

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief returns the highest bit set in a non-zero integer
///
/// This function returns the highest bit set in a non-zero integer.
/// \param[in] i the non-zero integer of which the highest bit is returned
/// \pre The integer @c i is non-zero

template <class Integral>
constexpr unsigned highestBit(Integral i)
{
#if __cplusplus > 201103
  assert(i != 0);  // needs C++14 to remain constexpr
#endif
  return detail::highestBitImpl(i, 0u);
}

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Returns the logarithm base 2 of a non-zero integer.
///
/// This function returns the the logarithm base 2 of a non-zero integer
/// \param[in] i the non-zero integer of which the logarithm base 2 is returned
/// \pre The integer @c i is a power of 2

template <class Integral>
unsigned int ilog2(Integral n)
{
  for (unsigned width = 0; width < 8 * sizeof(Integral); ++width)
    if ((static_cast<Integral>(1) << width) == n) return width;
  // not a power of 2
  assert(false);
  return 0;  // dummy return
}

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief checks whether an integer is a power of 2
///
/// \param[in] i an integer to be checked

template <class Integral>
inline constexpr bool isPowerOf2(Integral i)
{
  return i == static_cast<Integral>(0) ? false
                                       : static_cast<unsigned long>(i) == (1UL << highestBit(i));
}

/////////////////////////////////////////////////////////////////////////////////////////

namespace detail {

inline int BX(long x)
{
  return ((x) - (((x) >> 1) & 0x77777777) - (((x) >> 2) & 0x33333333) -
          (((x) >> 3) & 0x11111111));
}

}	// namespace detail

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Determines the population count (number of set bits) of a 32bit integer
///
/// \param[in] x Integer of which to determine the population count
inline long popcnt(uint32_t x)
{
  return (((detail::BX(x) + (detail::BX(x) >> 4)) & 0x0F0F0F0F) % 255);
}

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Determines the population count (number of set bits) of a 64bit integer
/// using alps::popcnt
///
/// \param[in] x Integer of which to determine the population count
inline long popcnt(uint64_t x)
{
  uint64_t mask32 = (1UL << 32) - 1;
  return popcnt((uint32_t)(x & mask32)) + popcnt((uint32_t)((x >> 32) & mask32));
}

/////////////////////////////////////////////////////////////////////////////////////////

}	// namespace iqs

/// @}*/

#endif	// header guard IQS_BITOPS_HPP
