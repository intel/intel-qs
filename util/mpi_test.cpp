// Copyright (C) 2015 Theoretical Physics, ETHZ Zurich

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

/** @file mpi_test.cpp
 *
 *  Tests for mpi.hpp
 */

#include "openqu/util/mpi.hpp"
#include <cassert>
#include <iostream>

int main(int argc, char** argv)
{
  // initialize the MPI environment if MPI exists
  openqu::mpi::Environment env(argc, argv);

  // these should work even without MPI
  assert(openqu::mpi::Environment::size() > 0);
  assert(openqu::mpi::Environment::rank() < openqu::mpi::Environment::size());

#ifdef OPENQU_HAVE_MPI
  OPENQU_MPI_CHECK_RESULT((MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN)))
  // test if exceptions work
  try {
    int data;
    // start an invalid send
    OPENQU_MPI_CHECK_RESULT(
        (MPI_Send(&data, 1, MPI_INT, openqu::mpi::Environment::size(), 0, MPI_COMM_WORLD)))
    assert(false);
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
  }
#endif

  openqu::mpi::print("Hello World!", false);
  openqu::mpi::print("Hello World!", true);

  return 0;
}
