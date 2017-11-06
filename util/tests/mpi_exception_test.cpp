//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

/**
 * @file mpi_exception_test.cpp
 *
 * Test the MPI Wrapper exception handling code.
 */
#include "../mpi_wrapper.hpp"
#include <iostream>
#include <cassert>


using namespace std;


int main(int argc, char** argv)
{
    // Initialize MPI framework.
    qhipster::MpiWrapper mpi(argc,argv);

    try {

        // Set the error handler for C++ framework calls through the MPI_COMM_WORLD communicator.
        QH_MPI_STATUS_CHECK(MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN));

        // do something dumb.
        QH_MPI_STATUS_CHECK((MPI_Bsend(&argc, 1, MPI_INT, 600, 0, MPI_COMM_WORLD)));

        // No more MPI calls allowed after this.
        QH_MPI_STATUS_CHECK(MPI_Barrier(MPI_COMM_WORLD));
        QH_MPI_STATUS_CHECK(MPI_Finalize());
    }
    catch (exception const& e) {
        cerr << e.what() << endl;
    }

  return 0;
}
