// Copyright (C) 2016 Theoretical Physics, ETHZ Zurich

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

#pragma once

/** \addtogroup util
 *  @{
 */

/** @file mpi.hpp
 *
 *  This header includes the MPI library if available and provides some
 *  convenience functions that are implemented even if the MPI library is not provided
 */

#ifdef OPENQU_HAVE_MPI
#include "mpi.h"
#endif

#ifdef _OPENMP
#include "openmp_affinity_corei7.hpp"
extern qhipster::openmp::AffinityCoreI7 glb_affinity;
#else
#include "openmp_affinity_noomp.hpp"
extern qhipster::openmp::AffinityNoOmp glb_affinity;
#endif

#include <exception>
#include <string>
#include <unistd.h>

namespace openqu {

namespace mpi {

/**
 * A trimmed down version of the BOOST::MPI environment.  It's purpose is to 
 * initialize the MPI library and partition the cluster or single threaded
 * environment for parallel operations.
 */
class Environment
{

  public:
  /** Intialize the MPI Environment
   *
   * passing in the argc and argv arguments passed to the main function.
   * If MPI is present but has not been initialized then MPI_Init will be called
   *
   */

  Environment(int& argc, char**& argv);

  /** Finalize the MPI Environment
   *
   * If MPI is present and has been initizlized in the constructor then
   * MPI_Finalize will be called here.
   *
   */

  ~Environment();
  Environment(Environment const&) = delete;
  Environment& operator=(Environment const&) = delete;

  bool is_usefull_rank() {return useful_rank;}

  /** the rank of the current MPI process
   *
   * returns the global rank of the current MPI process in the MPI_COMM_WORLD
   * communicator or 0 if MPI is not present.
   *
   * @pre If MPI is present, this can only be called after intializing MPI.
   */

  static int rank();

  /** the total nunber of MPI processes
   *
   * returns the size of the MPI_COMM_WORLD communicator, indicating the number
   * of MPI processes available or 1 if MPI is not present.
   *
   * @pre If MPI is present, this can only be called after intializing MPI.
   */

  static int size();
  static int get_nrankspernode();
  static int get_nnodes();
  static int get_nodeid();
#ifdef OPENQU_HAVE_MPI
  static MPI_Comm comm();
#endif
  static void remaprank(int newme);

#ifdef OPENQU_HAVE_MPI
  static MPI_Comm communicator;
#endif

 private:
  bool inited_;
  bool useful_rank;
  unsigned nrankspernode;
  static unsigned nnodes;
  static unsigned mynodeid;
  MPI_Request synch_request;

};

/** @brief An MPI barrier
 *
 * waits until all MPI processes have reached this call. This function does nothing if
 * MPI is not present.
 */
void barrier();

/** @brief Print from all MPI processes
 *
 * prints a string from all processes if @c all is true or just from the master process
 * with rank 0 if @c all is false. If MPI is not present it prints the string.
 *
 * If @c all is set, the string is prefixed by the number of the MPI process.
 *
 * @param s the string to be printed
 * @param all a flag to specify if all processes should print or just the master process
 */
void print(std::string s, bool all);
void print(std::string s, MPI_Comm comm);

/** @brief Catch-all exception class for MPI errors.
 *
 * Similar to the MPI exeption class in the Boost libraries.
 * Instances of this class will be thrown when an MPI error
 * occurs.
 */
class Exception : public std::exception
{
 public:
  /**
   * Build a new @c exception exception.
   *
   *   @param routine The MPI routine in which the error
   *   occurred. This should be a pointer to a string constant: it
   *   will not be copied.
   *
   *   @param result_code The result code returned from the MPI
   *   routine that aborted with an error.
   */
  Exception(int result);

  virtual ~Exception() throw();

  /**
   * A description of the error that occurred.
   */
  virtual const char* what() const throw() { return this->text_.c_str(); }

  /**
   * @brief Obtain the result code returned from the MPI routine
   * that caused an error.
   */
  int errorCode() const { return result_; }

 private:
  int result_;
  std::string text_;
};

/**
 * Checks the return value of an MPI call MPICALL and throws an
 * openqu::mpi::exception if the result is not MPI_SUCCESS.
 */
#define OPENQU_MPI_CHECK_RESULT(MPICALL)                                           \
  {                                                                                \
    int _check_result = MPICALL;                                                   \
    if (_check_result != MPI_SUCCESS) throw openqu::mpi::Exception(_check_result); \
  }
}
}

/** @}*/
