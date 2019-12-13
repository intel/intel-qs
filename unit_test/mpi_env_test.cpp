//------------------------------------------------------------------------------
// Copyright 2019 Intel Corporation
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

/// @file mpi_env_test.cpp
///
/// Tests for the MPI environment created in ../util/mpi_env.hpp

#include <cassert>
#include <iostream>
#include <string>
#include <sys/time.h>

#include "../util/mpi_env.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
/// Extension of the assert function that prints to screen a description of the failure.

#ifndef ASSERT
# define ASSERT(condition, message) \
do { \
    if (! (condition)) { \
        std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                  << " line " << __LINE__ << std::endl \
                  << "\033[30;41mExplanation:\033[0m " << message << "\n\n"; \
        std::terminate(); \
    } \
} while (false)
#else
#pragma message "ASSERT() already defined"
#endif

/////////////////////////////////////////////////////////////////////////////////////////
/// Shortcut for the MPI_COMM_WORLD barrier.

#ifdef INTELQS_HAS_MPI
  #define WORLD_BARRIER \
  { MPI_Barrier(MPI_COMM_WORLD); }
#else 
  #define WORLD_BARRIER \
  {}
#endif

/////////////////////////////////////////////////////////////////////////////////////////
// main code
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  // Initialize the MPI environmentk, if MPI exists.
  qhipster::mpi::Environment env(argc, argv);
  // These should work even without MPI.
  int pool_rank = qhipster::mpi::Environment::GetPoolRank();
  int pool_size  = qhipster::mpi::Environment::GetPoolSize();

  // When MPI is used, obtain rank and size of MPI_COMM_WORLD.
  int world_rank = 0, world_size = 1;
#ifdef INTELQS_HAS_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

  // Trivial tests.
  ASSERT(pool_size > 0,"At least one process needed.");
  ASSERT(pool_rank < pool_size, "rank id >= num processes.");

  // Verify that a single state has been created.
  if (env.IsUsefulRank()==true)
  {
      ASSERT(pool_rank==qhipster::mpi::Environment::GetStateRank(),"Consistency check."); 
      ASSERT( pool_size>=qhipster::mpi::Environment::GetStateSize(),"Consistency check."); 
  }
  else
      std::cout << "\nDummy pool: world rank = " << world_rank << "\n"
                <<   "            pool  rank = " << pool_rank  << "\n\n";

  WORLD_BARRIER;

/////////////////////////////////////////////////////////////////////////////////////////

  std::string buffer;
  if (env.IsUsefulRank()==true)
    buffer = "Pool rank " + std::to_string(qhipster::mpi::Environment::GetPoolRank())
             + "/" + std::to_string(qhipster::mpi::Environment::GetPoolSize()) + " , "
             + "State rank " + std::to_string(qhipster::mpi::Environment::GetStateRank())
             + "/" + std::to_string(qhipster::mpi::Environment::GetStateSize()) + " , "
             + "Node id " + std::to_string(qhipster::mpi::Environment::GetNodeId())
             + "/" + std::to_string(qhipster::mpi::Environment::GetNumNodes());
  else
    buffer = "Pool rank " + std::to_string(qhipster::mpi::Environment::GetPoolRank())
             + "/" + std::to_string(qhipster::mpi::Environment::GetPoolSize()) + " , "
             + "Node id " + std::to_string(qhipster::mpi::Environment::GetNodeId())
             + "/" + std::to_string(qhipster::mpi::Environment::GetNumNodes());


/////////////////////////////////////////////////////////////////////////////////////////

  qhipster::mpi::PoolPrint("Hello Pool (printed only from pool_rank=0)!", false);
  
  WORLD_BARRIER;

/////////////////////////////////////////////////////////////////////////////////////////

  if (env.IsUsefulRank()==true)
  {
      qhipster::mpi::PoolPrint("\n --------- print once at a time, useful ranks.\n");
      qhipster::mpi::PoolPrint(buffer, true);
  }

  WORLD_BARRIER;

  if (env.IsUsefulRank()==false)
  {
      qhipster::mpi::PoolPrint("\n --------- print once at a time, dummy ranks.\n");
      qhipster::mpi::PoolPrint(buffer, true);
  }

  WORLD_BARRIER;

/////////////////////////////////////////////////////////////////////////////////////////

  if (world_rank==0)
      std::cout << "\n --------- test of mpi_exception (they should fail, one per rank).\n\n";

/////////////////////////////////////////////////////////////////////////////////////////
// test of the exception handling for the MPI environment

#ifdef INTELQS_HAS_MPI
  QHIPSTER_MPI_CHECK_RESULT(MPI_Errhandler_set,(MPI_COMM_WORLD, MPI_ERRORS_RETURN))
  // test if exceptions work
  try {
    int data;
    // start an invalid send
    QHIPSTER_MPI_CHECK_RESULT(
        MPI_Send,(&data, 1, MPI_INT, world_size+2, 0, MPI_COMM_WORLD))
    assert(false);
  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

/////////////////////////////////////////////////////////////////////////////////////////

  if (world_rank==0)
      std::cout << "\n --------- update the communicators by using 3 states.\n\n";

  int num_states = 3;

  if (world_size<num_states)
  {
      if (world_rank==0)
          std::cout << "There are not enoug MPI ranks for all states.\n";
  }
  else
  {
      env.UpdateStateComm(num_states);

      WORLD_BARRIER;

      buffer = "Pool rank " + std::to_string(qhipster::mpi::Environment::GetPoolRank())
               + "/" + std::to_string(qhipster::mpi::Environment::GetPoolSize()) + " , "
               + "State rank " + std::to_string(qhipster::mpi::Environment::GetStateRank())
               + "/" + std::to_string(qhipster::mpi::Environment::GetStateSize()) + " , "
               + "Node id " + std::to_string(qhipster::mpi::Environment::GetNodeId())
               + "/" + std::to_string(qhipster::mpi::Environment::GetNumNodes());

      if (env.IsUsefulRank()==true)
          qhipster::mpi::PoolPrint(buffer, true);
      else
          std::cout << "[dummy world rank = " << world_rank << "] pool size = "
                    << env.GetPoolSize() << "\n";

      if (env.IsUsefulRank()==true)
          ASSERT(env.GetNumStates()==3,"Check that the states in the pool are 3.");
      else
          ASSERT(env.GetNumStates()==3,"Check that the states in the pool are 3.");

  }

  WORLD_BARRIER;

/////////////////////////////////////////////////////////////////////////////////////////

  if (world_rank==0)
      std::cout << "\n --------- terminate the dummy processes.\n\n";

  if (env.IsUsefulRank()==false)
      return 0;

  qhipster::mpi::PoolBarrier();

  // Verify that the PoolBarrier works.
  qhipster::mpi::PoolPrint(" -- before the PoolBarrier()\n");
  qhipster::mpi::PoolBarrier();
  qhipster::mpi::PoolPrint(" --  after the PoolBarrier()\n\n");
 
/////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}
