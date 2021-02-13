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

#include "../include/qureg.hpp"

// The scope is applying a sequence of num_gates=40 gates to the quantum register.
// The form of each gate is the same and is provided as quantum channel described
// by its chi matrix.

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  int my_rank=0, tot_num_procs=1;
  iqs::mpi::Environment env(argc, argv);
#ifdef INTELQS_HAS_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &tot_num_procs);
#endif

  // Partition the MPI environment into groups of processes. One group per pool state.
  unsigned num_pool_states = tot_num_procs;
  env.UpdateStateComm(num_pool_states);
  assert(env.GetPoolSize()  == tot_num_procs && "Wrong number of states in the pool.");
  assert(env.GetNumStates() == tot_num_procs && "Wrong number of states in the pool.");
  assert(env.IsUsefulRank() == true && "No dummy process in this IQS application.");
  my_rank = env.GetPoolRank();

  unsigned num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
  {
      num_threads = omp_get_num_threads();
  }
#endif

/// --- PARAMETERS ------------------------------------------- ///
  // number of qubits
  unsigned num_qubits = 8;
  // number of (one-qubit) gates
  unsigned num_gates = 40;
  // number of repetition of the (stochastic) noisy circuit
  unsigned num_ensemble_states = num_pool_states*20;
/// ---------------------------------------------------------- ///
  // Quantum channel corresponding to an ideal Hadamard gate:
  CM4x4<ComplexDP> chi_hadamard;
  chi_hadamard(1, 1) = {0.5, 0};
  chi_hadamard(1, 3) = {0.5, 0};
  chi_hadamard(3, 1) = {0.5, 0};
  chi_hadamard(3, 3) = {0.5, 0};
  // Verify that other entries are initialized to {0, 0}.
  assert(std::norm(chi_hadamard(0, 0))==0);
  assert(std::norm(chi_hadamard(0, 1))==0);
  assert(std::norm(chi_hadamard(0, 2))==0);
  assert(std::norm(chi_hadamard(0, 3))==0);
  assert(std::norm(chi_hadamard(1, 0))==0);
  assert(std::norm(chi_hadamard(1, 2))==0);
  assert(std::norm(chi_hadamard(2, 0))==0);
  assert(std::norm(chi_hadamard(2, 1))==0);
  assert(std::norm(chi_hadamard(2, 2))==0);
  assert(std::norm(chi_hadamard(2, 3))==0);
  assert(std::norm(chi_hadamard(3, 0))==0);
  assert(std::norm(chi_hadamard(3, 2))==0);
  // Solve eigensystem of chi and normalize its eigenvalues as probabilities.
  chi_hadamard.SolveEigenSystem();
  if (!my_rank)
      chi_hadamard.Print();
/// ---------------------------------------------------------- ///

  if(argc != 2)
  {
      fprintf(stderr, "usage: %s <num_qubits> \n", argv[0]);
      exit(1);
  }
  else
  {
      num_qubits = atoi(argv[1]);
      if (!my_rank)
          std::cout << "-- Parameters:\n"
                    << "num pool states = " << num_pool_states << "\n"
                    << "num ensemble states = " << num_ensemble_states << "\n"
                    << "num qubits = " << num_qubits << "\n"
                    << "num gates = " << num_gates << "\n";
  }

/////////////////////////////////////////////////////////////////////////////////////////
  // Noiseless circuit. 
  if (!my_rank) std::cout << "-- Run noiseless circuit\n";

  // State initialization.
  iqs::QubitRegister<ComplexDP> psi0(num_qubits);
  psi0.Initialize("base",0);

  // Circuit formed by Hadamard gates:
  unsigned i = 0, qubit;
  while (i < num_gates)
  {
      qubit = i%num_qubits;
      psi0.ApplyHadamard(qubit);
      i+=1;
  }
  // Compute the probability of qubit 0 to be in |1>.
  double probability = psi0.GetProbability(0);

/////////////////////////////////////////////////////////////////////////////////////////
  // Noisy circuit. 

  // Declare the random number generator and initialize its seed.
  iqs::RandomNumberGenerator<double> rng;
  std::size_t rng_seed = 77777;
  rng.SetSeedStreamPtrs( rng_seed );

  // State initialization.
  iqs::QubitRegister<ComplexDP> psi1(num_qubits);
  psi1.SetRngPtr(&rng);

  // Circuit formed by Hadamard gates:
  double overlap_squared_noisy = 0.;
  double probability_noisy = 0.;
  for (unsigned j=0; j<num_ensemble_states/num_pool_states; j++)
  {
      std::cout << "-- Run noisy circuit (rank=" << my_rank << ", repetition=" << j <<")\n";
      psi1.Initialize("base",0);
      i=0;
      while (i < num_gates)
      {
          qubit = i%num_qubits;
          psi1.ApplyChannel(qubit, chi_hadamard);
          i+=1;
      }
      // Compute the probability of qubit 0 to be in |1>.
      probability_noisy += psi1.GetProbability(0);
      // Overlap with ideal state.
      overlap_squared_noisy += std::norm( psi0.ComputeOverlap(psi1) );
  }
  // Incoherent average across the pool and repetitions.
  overlap_squared_noisy = env.IncoherentSumOverAllStatesOfPool<double> (overlap_squared_noisy);
  overlap_squared_noisy /= double(num_ensemble_states);
  //
  probability_noisy = env.IncoherentSumOverAllStatesOfPool<double> (probability_noisy);
  probability_noisy /= double(num_ensemble_states);

/////////////////////////////////////////////////////////////////////////////////////////

  // Print a few information on screen.
  // Computation of the overlap between the ideal state and those exposed to noise:
  if (my_rank == 0) 
      std::cout << "---- summary of simulation: \n"
                << "Number of states in the ensamble = " << num_ensemble_states << "\n"
                << "(divided in " << num_pool_states << " parallel states)\n"
                << "Overlap-squared between ideal and noisy states = "
                << overlap_squared_noisy << "\n"
                << "Probability (of qubit 0 to be in |1>) in the noiseless case = " << probability << "\n"
                << "Probability (...                 ...) with noise = " << probability_noisy << "\n\n";

  return 1;
}
