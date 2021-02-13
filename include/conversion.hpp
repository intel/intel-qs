// Copyright (C) 2016 Theoretical Physics, ETHZ Zurich
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

#ifndef IQS_CONVERSION_HPP
#define IQS_CONVERSION_HPP

/// \addtogroup util
/// @{

/// @file conversion.hpp
///
/// This header defines type conversion functions.

#include <cassert>
#include <sstream>
#include <string>

namespace iqs {

/// @brief convert to a string
///
/// This function converts any value to a string, by writing it into a string stream.
/// \param[in] val the value to be converted to a string
/// \pre Writing into a std::istream using operator<< needs to be implemented for the type

template <class T>
std::string toString(T const& val)
{
  std::ostringstream os;
  os << val;
  return os.str();
}

}	// namespace iqs

/// @}*/

#endif	// header guard IQS_CONVERSION_HPP
