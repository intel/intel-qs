/// @file mpi_speedup_test.cpp
///
/// Tests for the speedup of the MPI environment of IQS.

#include <cmath>
#include <iostream>
#include <string>
#include <sys/time.h>

#include "../include/mpi_env.hpp"
#include "../include/qureg.hpp"

//#define MPI_WITHOUT_IQS_ENV

/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  int pool_rank;
  int pool_size;
#ifdef MPI_WITHOUT_IQS_ENV
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &pool_rank);
  if (!pool_rank)
      std::cout << "\n-- Without using the MPI environment of IQS --\n";
#else
  // Initialize the MPI environmentk, if MPI exists.
  iqs::mpi::Environment env(argc, argv);
  // These should work even without MPI.
  pool_rank = iqs::mpi::Environment::GetPoolRank();
  pool_size  = iqs::mpi::Environment::GetPoolSize();
  if (!pool_rank)
      std::cout << "\n-- Using the MPI environment of IQS --\n";
#endif

  std::size_t power_of_2 = 32;
  std::size_t global_value = 1L << power_of_2;
  std::size_t local_value = global_value/std::size_t(pool_size);

  if (!pool_rank)
      std::cout << "\nQuick test of MPI parallel speedup:\n"
                << "Count to 2^" << std::log2(global_value)
                << " in groups of 2^" << std::log2(local_value) << "\n";

  struct timeval time;
  double start, end;
  std::size_t dumb_counter=0;

  // Start the timer.
  gettimeofday(&time, (struct timezone*)0);
  start =  time.tv_sec + time.tv_usec * 1.0e-6;

  for (std::size_t j=0; j<local_value; ++j)
  {
     if (j%2==0)
         dumb_counter += 3;
     else
         dumb_counter -= 1;
  }

  // Stop the timer.
  gettimeofday(&time, (struct timezone*)0);
  end =  time.tv_sec + time.tv_usec * 1.0e-6;

  std::cout << "[" << pool_rank << "] Time elapsed = " << (end-start)
            << " to count to 2^" << std::log2(dumb_counter) << "\n\n";

#ifndef MPI_WITHOUT_IQS_ENV
#if 10
  iqs::mpi::StateBarrier();

  int num_qubits =26;
  if (!pool_rank)
      std::cout << "\nLet us initialize " << num_qubits << " qubits in |0>, "
                << "and perform a single 1-qubit-gate per qubit.\n";

  iqs::mpi::StateBarrier();

  iqs::QubitRegister<ComplexDP> psi(num_qubits, "base", 0);

  iqs::mpi::StateBarrier();
  if (!pool_rank) std::cout << "\nInitialization of the qubit register.\n";
  iqs::mpi::StateBarrier();

  psi.TurnOnSpecialize();
  psi.EnableStatistics();
  for(int qubit = 0; qubit < num_qubits; qubit++)
      psi.ApplyHadamard(qubit);
  psi.GetStatistics();
  psi.DisableStatistics();
#endif
#else
  MPI_Finalize();
#endif

  return 0;
}
