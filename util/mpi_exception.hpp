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

#ifndef IQS_MPI_EXCEPTION_HPP
#define IQS_MPI_EXCEPTION_HPP

/// \addtogroup util
/// @{

/// @file mpi_exception.hpp
///
/// This header includes the MPI library if available and provides some
/// convenience functions that are implemented even if the MPI library is not provided.

#include <exception>
#include <string>
#include <unistd.h>

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

namespace mpi {

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Catch-all exception class for MPI errors.
///
/// Similar to the MPI exception class in the Boost libraries.
/// Instances of this class will be thrown when an MPI error occurs.
class Exception : public std::exception
{
 public:

  /// Build a new @c Exception exception.
  ///
  /// @param routine The MPI routine in which the error occurred. This should be a pointer
  /// to a string constant: it will not be copied.
  ///
  /// @param error_code The result code returned from the MPI routine that aborted with
  /// an error.
  Exception(const char* routine, int error_code);

  virtual ~Exception() throw();

  /// A description of the error that occurred.
  virtual const char* what() const throw() { return this->description_.c_str(); }

  /// Retrieve the name of the MPI routine that reported the error.
  const char* routine() const { return routine_; }

  /// @brief Obtain the result code returned from the MPI routine that caused an error.
  int error_code() const { return error_code_; }

 private:
  /// The MPI routine that triggered the error.
  const char * routine_;
  /// The failure code reported by the MPI implementation.
  int error_code_;
  /// The formatted error description.
  std::string description_;
};

/////////////////////////////////////////////////////////////////////////////////////////

/// Call the MPI routine MPIFunc with arguments Args (surrounded by parentheses).
/// Checks the return value of MPIFunc call and throws a qhipster::mpi::exception
/// if the result is not MPI_SUCCESS.
#define QHIPSTER_MPI_CHECK_RESULT(MPIFunc, Args)                                           \
  {                                                                                        \
    int check_result = MPIFunc Args;                                                       \
    if (check_result != MPI_SUCCESS) throw qhipster::mpi::Exception(#MPIFunc,check_result);\
  }

/////////////////////////////////////////////////////////////////////////////////////////

}	// end namespace mpi
}	// end namespace qhipster

/// @}*/

#endif	// header guard IQS_MPI_EXCEPTION_HPP
