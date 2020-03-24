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

/// @file mpi_env.cpp
///
///  This header implements MPI support functionality

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

#include <cassert>
#include <stdexcept>

#include "mpi_exception.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

namespace mpi {

/////////////////////////////////////////////////////////////////////////////////////////

Exception::Exception(const char* routine, int error_code)
  : routine_(routine), error_code_(error_code)
{
#ifndef INTELQS_HAS_MPI
  assert(0);
#else
  // Get the error message from the MPI implementation.
  char buffer[MPI_MAX_ERROR_STRING];
  int len;
  MPI_Error_string(error_code, buffer, &len);

  if (routine!=nullptr)
    description_ = std::string(routine) + ": " + std::string(buffer);
  else
    description_ = "Error in exception handling since 'routine' or 'buffer' are NULL";
#endif
}

Exception::~Exception() throw() {}

/////////////////////////////////////////////////////////////////////////////////////////

}	// end namespace mpi
}	// end namespace qhipster

/////////////////////////////////////////////////////////////////////////////////////////
