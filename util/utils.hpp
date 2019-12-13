//------------------------------------------------------------------------------
// Copyright (C) 2017 Intel Corporation 
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
//------------------------------------------------------------------------------

#ifndef IQS_UTILS_HPP
#define IQS_UTILS_HPP

#include <complex>

// Helpful defines, if not already provided.
#define DO_PRAGMA(x) _Pragma(#x)

#ifndef TODO
#define TODO(x) DO_PRAGMA(message("\033[30;43mTODO\033[0m - " #x))
#endif

#ifndef INFO
#define INFO(x) DO_PRAGMA(message("\033[30;46mINFO\033[0m - " #x))
#endif

#define D(x) ((double)(x))
#define UL(x) ((std::size_t)(x))
#define sec() time_in_seconds()
#define xstr(s) __str__(s)
#define __str__(s) #s

/////////////////////////////////////////////////////////////////////////////////////////

using ComplexSP = std::complex<float>;
using ComplexDP = std::complex<double>;

double time_in_seconds(void);

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

/// Utility method to inform on the currently set compiler flags.
void WhatCompileDefinitions();

}

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard IQS_UTILS_HPP
