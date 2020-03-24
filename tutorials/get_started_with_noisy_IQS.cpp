//------------------------------------------------------------------------------
// Copyright (C) 2019 Intel Corporation 
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

/// @file getting_started_with_IQS.cpp
/// Tutorial on the use of Intel Quantum Simulator (IQS) for parallel noisy simulation.

#include <cassert>
#include <iostream>
#include <math.h>	// for constant M_PI
#include <vector>

// Explicitly include MPI methods when appropriate flag is defined.
#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

// Include the header file with the declaration of all classes and methods of IQS.
#include "../qureg/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

// Start of the main program (C++ language).
int main(int argc, char **argv)
{
#ifndef INTELQS_HAS_MPI
  std::cout << "\nThis noisy simulation is thought to be run with MPI.\n"
            << "To do so, please set the option '-DIqsMPI=ON' when calling CMake.\n\n"
            << "However the code will execute also without MPI.\n\n";
#endif

/////////////////////////////////////////////////////////////////////////////////////////
// Setting the MPI environment
/////////////////////////////////////////////////////////////////////////////////////////

  // Create the MPI environment, passing the same argument to all the ranks.
  qhipster::mpi::Environment env(argc, argv);

  // In IQS the noise is introduced using the Stochastic Schroedinger Equation,
  // implemented via noise gates. This means that a single (noisy) quantum circuit
  // is substituted by an ensemble of many (ideal) circuits, each having
  // additional stochastic gates and corresponding to a specific realization
  // of the stochastic term of the equation.

  // First of all, one needs to decide how many distinct states to simulate in the pool.
  // In this example, we use a single rank per state.
  int tot_num_ranks = 1;
#ifdef INTELQS_HAS_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &tot_num_ranks);
#endif
  int num_pool_states = tot_num_ranks;

  // NOTE: Why didn't we use the method 'qhipster::mpi::Environment::GetPoolSize()'?
  //       Because the initial environment cuts the number of ranks down to the closest
  //       power of 2.

  // Partition the MPI environment into groups of processes. One group per pool state.
  env.UpdateStateComm(num_pool_states);
  assert(env.GetPoolSize()  == tot_num_ranks);
  assert(env.GetNumStates() == tot_num_ranks);
  assert(env.IsUsefulRank() == true);

  // Second, one needs to decide how many states form the ensemble.
  // Here we choose the smaller multiple of the num_pool_states grater or equal to 200.
  int min_num_ensemble_states = 200;
  int num_ensemble_states = min_num_ensemble_states;
  if (min_num_ensemble_states%num_pool_states != 0 )
      num_ensemble_states += num_pool_states-(min_num_ensemble_states%num_pool_states);
  assert (num_ensemble_states%num_pool_states == 0 );

  // IQS has functions that simplify some MPI instructions. However, it is important
  // to keep trace of the current rank.
  int my_rank = env.GetPoolRank();
  {
      int world_rank = 0;
#ifdef INTELQS_HAS_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif
      assert(my_rank==world_rank);
  }

/////////////////////////////////////////////////////////////////////////////////////////
// Quantum state initialization and ideal circuit
/////////////////////////////////////////////////////////////////////////////////////////

  // Number of qubits
  int num_qubits = 8;

  // The circuit consist of three single-qubit rotations per qubit, the first
  // in the X basis, the second in the Y basis and the third in the Z basis.

  // The noisy simulation depends on the gate parallelism. Here we consider the
  // sequential application of the gates, one at a time. The order is:
  // 1. X rotations are implemented, starting from qubit 0 and in increasing qubit order;
  // 2. Y rotations are implemented, starting from qubit 0 and in increasing qubit order;
  // 3. Z rotations are implemented, starting from qubit 0 and in increasing qubit order.

  // To generate random numbers, IQS provides a wrapper around VSL random number generator.
  // If MKL is not available, a standard MT19937 generator is used.
  // We need to declare the (pseudo) random number generator...
  qhipster::RandomNumberGenerator<double> rng;
  // ... and initialize its seed:
  std::size_t rng_seed = 77777;
  rng.SetSeedStreamPtrs( rng_seed );

  // NOTE: the random number generator is able to generate three different kinds
  //       of random numbers:
  //       *local* --> different for each pool rank
  //       *state* --> common to all ranks of the same state
  //       *pool*  --> common to all ranks of the pool

  // All angles of rotations are random. However once the circuit we want to simulate
  // is determined, there is no stochasticity in the rotation angles across the ensemble.
  // The rotation angles must therefore be given as *pool* random numbers.
  std::vector<double> x_angles (num_qubits), y_angles (num_qubits), z_angles (num_qubits);
  rng.UniformRandomNumbers( x_angles.data(), x_angles.size(), 0., M_PI, "pool");
  rng.UniformRandomNumbers( y_angles.data(), y_angles.size(), 0., M_PI, "pool");
  rng.UniformRandomNumbers( z_angles.data(), z_angles.size(), 0., M_PI, "pool");
  // Ideal (i.e. noiseless) state.
  // |psi> = |00000000>
  QubitRegister<ComplexDP> psi(num_qubits);
  psi.Initialize("base",0);

  // At this point we have 5 copies of the ideal state. One for each state in the pool.
  // We enable the collection of statistics from state_id = 0 only.
  if (env.GetStateId()==0)
      psi.EnableStatistics();
  if (my_rank==0) std::cout << "---- ideal circuit \n\n";
  for (int q=0; q<num_qubits; ++q)
  {
      psi.ApplyRotationX (q, x_angles[q]);
      psi.ApplyRotationY (q, y_angles[q]);
      psi.ApplyRotationZ (q, z_angles[q]);
  }
  if (env.GetStateId()==0)
  { 
      std::cout << "When 'Statistics' is enabled, one can print to screen a few info "
                << "related to the time-cost of simulating the quantum circuit.\n";
      psi.GetStatistics();
  }

#ifdef INTELQS_HAS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Compute the probability of qubit 0 to be in |1>.
  double probability = psi.GetProbability(0);

/////////////////////////////////////////////////////////////////////////////////////////
// Quantum state evolution in presence of noise
/////////////////////////////////////////////////////////////////////////////////////////

  // State for slow decoherence.
  QubitRegister<ComplexDP> psi_slow(num_qubits);
  // One can use the same random number generator for each state or a different one.
  // Here we use the same.
  psi_slow.SetRngPtr(&rng);
  // T_1 and T_2 times for slow decoherence
  double T_1_slow=1000. , T_2_slow=500. ;
  psi_slow.SetNoiseTimescales(T_1_slow, T_2_slow);

  // State for fast decoherence.
  QubitRegister<ComplexDP> psi_fast(num_qubits);
  // Here too we use the same random number generator.
  psi_fast.SetRngPtr(&rng);
  // T_1 and T_2 times for fast decoherence
  double T_1_fast=40.   , T_2_fast=20.  ;
  psi_fast.SetNoiseTimescales(T_1_fast, T_2_fast);

  // All single-qubit rotations have the same duration:
  double duration = 1.5;
 
// ---------------- slow decoherence
  if (my_rank==0) std::cout << "\n---- slow decoherence \n\n";
  qhipster::mpi::PoolBarrier();
  double overlap_squared_slow = 0.;
  double probability_slow = 0.;
  for (int j=0; j<num_ensemble_states/num_pool_states; j++)
  {
      // Quantum circuit with explicit noise gates added to simulate noise.
      psi_slow.Initialize("base", 0);
      // Noise after state preparation:
      for (int q=0; q<num_qubits; ++q)
          psi_slow.ApplyNoiseGate (q, duration);
      // X-basis rotations.
      for (int qubit=0; qubit<num_qubits; ++qubit)
      {
          psi_slow.ApplyRotationX (qubit, x_angles[qubit]);
          // Since the gates are sequential, noise acts on all qubits after each gate.
          for (int q=0; q<num_qubits; ++q)
              psi_slow.ApplyNoiseGate (q, duration);
      }
      // Y-basis rotations.
      for (int qubit=0; qubit<num_qubits; ++qubit)
      {
          psi_slow.ApplyRotationY (qubit, y_angles[qubit]);
          // Since the gates are sequential, noise acts on all qubits after each gate.
          for (int q=0; q<num_qubits; ++q)
              psi_slow.ApplyNoiseGate (q, duration);
      }
      // Z-basis rotations.
      for (int qubit=0; qubit<num_qubits; ++qubit)
      {
          psi_slow.ApplyRotationZ (qubit, z_angles[qubit]);
          // Since the gates are sequential, noise acts on all qubits after each gate.
          for (int q=0; q<num_qubits; ++q)
              psi_slow.ApplyNoiseGate (q, duration);
      }
      // Noise before state measurement has already been implemented.

      // Compute the probability of qubit 0 to be in |1>.
      probability_slow += psi_slow.GetProbability(0);

      // Overlap with ideal state.
      overlap_squared_slow += std::norm( psi.ComputeOverlap(psi_slow) );
  }

  // Compute average per group in the pool.
  overlap_squared_slow /= (double)(num_ensemble_states/num_pool_states);
  probability_slow /= (double)(num_ensemble_states/num_pool_states);

  // Incoherent average across the pool.
  overlap_squared_slow = env.IncoherentSumOverAllStatesOfPool<double> (overlap_squared_slow);
  overlap_squared_slow /= double(qhipster::mpi::Environment::GetNumStates());
  //
  probability_slow = env.IncoherentSumOverAllStatesOfPool<double> (probability_slow);
  probability_slow /= double(qhipster::mpi::Environment::GetNumStates());

  // NOTE: For the noise model considered, noise gates can be fused with each other.
  //       In fact, the overall effect on the ensemble is the same for two consecutive
  //       noise gates (on the same qubit) as a single one with duration equal to the
  //       sum of the individual durations.
  //       Exploiting this fact can reduce the computation time considerably.
  //       In the example of this tutorial, instead of using:
  //           N + 3*N*N 
  //       noise gates one could have used just:
  //           4*N
  //       noise gates.

// ---------------- fast decoherence
  if (my_rank==0) std::cout << "---- fast decoherence \n\n";
  qhipster::mpi::PoolBarrier();
  double overlap_squared_fast = 0.;
  double probability_fast = 0.;
  for (int j=0; j<num_ensemble_states/num_pool_states; j++)
  {
      // Quantum circuit with explicit noise gates added to simulate noise.
      psi_fast.Initialize("base", 0);
      // Noise after state preparation:
      for (int q=0; q<num_qubits; ++q)
          psi_fast.ApplyNoiseGate (q, duration);
      // X-basis rotations.
      for (int qubit=0; qubit<num_qubits; ++qubit)
      {
          psi_fast.ApplyRotationX (qubit, x_angles[qubit]);
          // Since the gates are sequential, noise acts on all qubits after each gate.
          for (int q=0; q<num_qubits; ++q)
              psi_fast.ApplyNoiseGate (q, duration);
      }
      // Y-basis rotations.
      for (int qubit=0; qubit<num_qubits; ++qubit)
      {
          psi_fast.ApplyRotationY (qubit, y_angles[qubit]);
          // Since the gates are sequential, noise acts on all qubits after each gate.
          for (int q=0; q<num_qubits; ++q)
              psi_fast.ApplyNoiseGate (q, duration);
      }
      // Z-basis rotations.
      for (int qubit=0; qubit<num_qubits; ++qubit)
      {
          psi_fast.ApplyRotationZ (qubit, z_angles[qubit]);
          // Since the gates are sequential, noise acts on all qubits after each gate.
          for (int q=0; q<num_qubits; ++q)
              psi_fast.ApplyNoiseGate (q, duration);
      }
      // Noise before state measurement has already been implemented.

      // Compute the probability of qubit 0 to be in |1>.
      probability_fast += psi_fast.GetProbability(0);
      // 

      // Overlap with ideal state.
      overlap_squared_fast += std::norm( psi.ComputeOverlap(psi_fast) );
  }

  // Compute average over all ensemble.
  overlap_squared_fast = env.IncoherentSumOverAllStatesOfPool<double> (overlap_squared_fast);
  overlap_squared_fast /= double(num_ensemble_states);
  //
  probability_fast = env.IncoherentSumOverAllStatesOfPool<double> (probability_fast);
  probability_fast /= double(num_ensemble_states);

/////////////////////////////////////////////////////////////////////////////////////////
// Quantum state initialization and ideal circuit
/////////////////////////////////////////////////////////////////////////////////////////

  // Print a few information on screen.
  // Computation of the overlap between the ideal state and those exposed to noise:
  if (my_rank == 0) 
      std::cout << "---- summary of simulation: \n"
                << "Overlap-squared between ideal and 'slow decoherence' state = "
                << overlap_squared_slow << "\n"
                << "Overlap-squared between ideal and 'fast decoherence' state = "
                << overlap_squared_fast << "\n----\n"
                << "Probability in the noiseless case = " << probability << "\n"
                << "Probability with slow decoherence = " << probability_slow << "\n"
                << "Probability with fast decoherence = " << probability_fast << "\n\n";

//  double e = psi2.MaxAbsDiff(psi1);

  return 1;
}
