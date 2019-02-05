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

#include "../qureg/qureg.hpp"

// The scope is applying a sequence of num_gates=40 gates to the quantum register.
// The form of each gate is the same:
//   controlled 1-qubit operation defined by the 2x2 matrix G
// but each pair (control,target) for the involved qubit is randomly generated.

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  openqu::mpi::Environment env(argc, argv);
  myrank = env.rank();
  nprocs = openqu::mpi::Environment::size();
  if (env.is_usefull_rank() == false) return 0;
#endif
  int num_threads = 88;
  glb_affinity.set_thread_affinity(num_threads);

/// --- PARAMETERS ------------------------------------------- ///
  // number of qubits
  unsigned num_qubits = 8;
  // number of (two-qubit) gates
  unsigned num_gates=20;
  // number of repetition of the (stochastic) noisy circuit
  unsigned num_noisy_circuits=num_gates*10;
  // T_1 and T_2 times for slow decoherence
  double T_1_slow=1000. , T_2_slow=500. ;
  double T_1_fast=40.   , T_2_fast=20.  ;
  // T_1 and T_2 times for slow decoherence
/// ---------------------------------------------------------- ///

  std::size_t tmp_size = 0;
  if(argc != 2)
  {
      fprintf(stderr, "usage: %s <num_qubits> \n", argv[0]);
      exit(1);
  }
  else
  {
      num_qubits = atoi(argv[1]);
  }

  // single-qubit operation that will be implemented (in a conditioned way)
  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};

  // simplified G gate: Pauli X
/*  G(0, 0) = {0. , 0.};
  G(0, 1) = {1. , 0.};
  G(1, 0) = {1. , 0.};
  G(1, 1) = {0. , 0.}; */

  // generate random pairs for control/target qubits
  // for a total of 20 2-qubit gates
  std::vector<std::pair<unsigned, unsigned>> qpair;
  std::default_random_engine generator;
  unsigned RNG_seed = 12345;
  std::uniform_int_distribution<int> qubit1(0, num_qubits - 1);
  std::uniform_int_distribution<int> qubit2(0, num_qubits - 1);
  unsigned i = 0;
  while (i < num_gates)
  {
      unsigned q1 = qubit1(generator);
      unsigned q2 = qubit2(generator);
      if (q1 != q2)
      {
         qpair.push_back(std::make_pair(q1, q2));
         i++;
      }
  }


  // ideal state
  QubitRegister<ComplexDP> psi0(num_qubits);
  // slow decoherence
  NoisyQureg<ComplexDP> psi1(num_qubits, RNG_seed, T_1_slow, T_2_slow);
  // fast decoherence
  NoisyQureg<ComplexDP> psi2(num_qubits, RNG_seed, T_1_fast, T_2_fast);
  psi0.Initialize("base", 0);


  // ---------------- no noise
  {
      psi0.EnableStatistics();
      for (auto &p : qpair)
          psi0.ApplyControlled1QubitGate(p.first, p.second, G);
    
      for(int pos = 0; pos < num_qubits; pos++)
          psi0.Apply1QubitGate(pos, G);
      psi0.GetStatistics();
  }

 
// ---------------- slow decoherence
  std::cout << " slow decoherence \n";
  double over_sq_1 = 0.;
  for (unsigned j=0; j<num_noisy_circuits; j++)
  {
      psi1.Initialize("base", 0);
      psi1.ResetTimeForAllQubits();
      for (auto &p : qpair)
          psi1.ApplyControlled1QubitGate(p.first, p.second, G);
    
      for(int pos = 0; pos < num_qubits; pos++)
          psi1.Apply1QubitGate(pos, G);

      psi1.ApplyNoiseGatesOnAllQubits();
      over_sq_1 = over_sq_1 + std::norm( psi0.ComputeOverlap(psi1) ) ;
      if (j%100==0)
          std::cout << " j=" << j << " , logical gate count for"
                    << " psi1 =" << psi1.GetTotalExperimentalGateCount()
                    << " (they should be " << num_gates+num_qubits << ") \n";
  }
  over_sq_1 = over_sq_1/(double)num_noisy_circuits;


  // ---------------- fast decoherence
  std::cout << " fast decoherence \n";
  double over_sq_2 = 0.;
  for (unsigned j=0; j<num_noisy_circuits; j++)
  {
      psi2.Initialize("base", 0);
      psi2.ResetTimeForAllQubits();
      for (auto &p : qpair)
          psi2.ApplyControlled1QubitGate(p.first, p.second, G);
    
      for (int pos = 0; pos < num_qubits; pos++)
          psi2.Apply1QubitGate(pos, G);

      psi2.ApplyNoiseGatesOnAllQubits();
      over_sq_2 = over_sq_2 + std::norm( psi0.ComputeOverlap(psi2) ) ;
  }
  over_sq_2 = over_sq_2/(double)num_noisy_circuits;

  // ---------------- 
  // computation of the overlap between the ideal state and those exposed to noise
  if (myrank == 0) 
      std::cout << " Overlap-squared between ideal and 'slow decoherence' state = "
                << over_sq_1 << "\n"
                << " Overlap-squared between ideal and 'fast decoherence' state = "
                << over_sq_2 << "\n";

//  double e = psi2.maxabsdiff(psi1);
}
