//------------------------------------------------------------------------------
// Copyright (C) 2020 Intel Corporation 
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

#include <sys/time.h>

#include "../include/qureg.hpp"

#ifdef USE_MKL
#include <mkl.h>
#endif

/// @file communication_reduction_via_qubit_reordering.cpp
/// Simple program to illustrate how the communication overhead between MPI
/// processes can be reduced by reordering the qubits.


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  qhipster::mpi::Environment env(argc, argv);
  if (env.IsUsefulRank() == false) return 0;
  int myrank = env.GetStateRank();
  if (myrank==0)
      std::cout << " ------------------------------------------------------------ \n"
                << " ---- Reducing communication overhead via qubit reordering -- \n"
                << " ------------------------------------------------------------ \n";

  int num_qubits = 22;
  if (argc != 2)
  {
    if (myrank==0)
        fprintf(stderr, "usage: %s <num_qubits>\n", argv[0]);
    exit(1);
  }
  else
  {
    num_qubits = atoi(argv[1]);
  }

  if (num_qubits<22)
  {
    if (myrank==0)
        fprintf(stderr, "simulate circuits of at least 22 qubits\n");
    exit(1);
  }


#pragma omp parallel
#pragma omp master
  {
    int nthreads = omp_get_num_threads();
    if (myrank==0)
      fprintf(stdout, "OMP number of threads = %d \n", nthreads);
  }

/////////////////////////////////////////////////////////////////////////////////////////

  // Variables for the timing.
  struct timeval time;
  double start, end;
  double duration_trivial_order;
  double duration_inverse_order;

  // We consider a simple circuit, with Hadamard, X, Y, Z acting on the last 10 qubits.
  QubitRegister<ComplexDP> psi_trivial_order(num_qubits, "base", 0);
  psi_trivial_order.ApplyPauliZ(num_qubits-1); // dummy operation to avoid the first MPI communication during timing

  std::vector<std::size_t> original_map = psi_trivial_order.qubit_permutation->map;
  std::vector<std::size_t> new_map(num_qubits);
  for (std::size_t q = 0; q < num_qubits; q++)
  {
      // Verify that qubit order is trivial.
      assert(original_map[q] == q);
      // create inverse map.
      new_map[num_qubits-q-1] = q;
  }
 
  // Simulate circuit.
  qhipster::mpi::StateBarrier();
  gettimeofday(&time, (struct timezone*)0);
  start =  time.tv_sec + time.tv_usec * 1.0e-6;
  for (int q=num_qubits-10; q<num_qubits; ++q)
  {
    psi_trivial_order.ApplyHadamard(q);
    psi_trivial_order.ApplyPauliX(q);
    psi_trivial_order.ApplyPauliY(q);
    psi_trivial_order.ApplyPauliZ(q);
  } 
  qhipster::mpi::StateBarrier();
  gettimeofday(&time, (struct timezone*)0);
  end =  time.tv_sec + time.tv_usec * 1.0e-6;
  if (myrank==0)
      std::cout << "trivial qubit order --> Simulation time = " << end-start << "\n";
  
/////////////////////////////////////////////////////////////////////////////////////////

  QubitRegister<ComplexDP> psi_inverse_order(num_qubits, "base", 0);
  psi_inverse_order.PermuteQubits(new_map, "direct");
  psi_trivial_order.ApplyPauliZ(0); // dummy operation to avoid the first MPI communication during timing
  // Simulate circuit.
  qhipster::mpi::StateBarrier();
  gettimeofday(&time, (struct timezone*)0);
  start =  time.tv_sec + time.tv_usec * 1.0e-6;
  for (int q=num_qubits-10; q<num_qubits; ++q)
  {
    psi_inverse_order.ApplyHadamard(q);
    psi_inverse_order.ApplyPauliX(q);
    psi_inverse_order.ApplyPauliY(q);
    psi_inverse_order.ApplyPauliZ(q);
  } 
  qhipster::mpi::StateBarrier();
  gettimeofday(&time, (struct timezone*)0);
  end =  time.tv_sec + time.tv_usec * 1.0e-6;
  if (myrank==0)
      std::cout << "inverse qubit order --> Simulation time = " << end-start << "\n";
    
  // Verify that the final state is the same.
  psi_inverse_order.PermuteQubits(original_map, "direct");
  double overlap_sq = std::norm(psi_trivial_order.ComputeOverlap(psi_inverse_order));
  if (myrank==0)
      std::cout << "Squared overlap of states at the end of the two simulations = " << overlap_sq << "\n\n";
         
  return 0;
}
