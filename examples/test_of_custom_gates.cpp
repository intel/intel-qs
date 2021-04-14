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

#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
  iqs::mpi::Environment env(argc, argv);
  unsigned myrank = env.GetStateRank();
  unsigned nprocs = iqs::mpi::Environment::GetStateSize();

/// --- PARAMETERS ------------------------------------------- ///
  int num_qubits = 8;
  int num_gates = 1;
/// ---------------------------------------------------------- ///
  std::size_t tmp_size = 0;
  if(argc != 2)
  {
      if (!myrank)
          fprintf(stderr, "usage: %s <num_qubits> \n", argv[0]);
      exit(1);
  }
  else
  {
      int _qubits = atoi(argv[1]);
      if ((_qubits <= 0) || (_qubits > 1000000))
      {
          fprintf(stderr, "<num_qubits> was (%d) which is invalid or negative.\n", _qubits);
          exit(1);
      }
      num_qubits = (unsigned)_qubits; 
  }
  if (env.IsUsefulRank() == false) return 0;

#pragma omp parallel
#pragma omp master
  {
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_num_threads();
#endif
    fprintf(stdout, "OMP number of threads = %d \n", nthreads);
  }
  
  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};


/// --- RUN SINGLE-QUBIT GATES ------------------------------- ///

// Initialize the qubit register in state |A>, then apply the custom one-qubit gate G
// to each qubit sequentially.

  iqs::QubitRegister<ComplexDP> psi_A(num_qubits, "base", 0);
//  iqs::QubitRegister<ComplexDP> psi_A(num_qubits, "rand", -1);

  // with specialization
  psi_A.TurnOnSpecialize();
  for(int pos = 0; pos < num_qubits; pos++)
  {
      if (myrank == 0) printf(" ---------------------------------- \n");
      psi_A.EnableStatistics();
      psi_A.Apply1QubitGate(pos, G);
      psi_A.GetStatistics();
      psi_A.DisableStatistics();
  }

/// --- RUN SINGLE-QUBIT GATES ------------------------------- ///

// Initialize two quantum states: |B> and |C>=|B>.
// Then generate random pairs of control-target qubits.
// Apply both controlled-gates and one-qubit gates sequentially, first on |B> and
// then on |C>. The gates on |C> have the specialization 'on'.
// Verify that the final state |B> and |C> are equivalent.

  std::vector<std::pair<unsigned, unsigned>> qpair;
  std::default_random_engine generator;
  std::uniform_int_distribution<int> qubit1(0, num_qubits - 1);
  std::uniform_int_distribution<int> qubit2(0, num_qubits - 1);
  unsigned i = 0;
  while (i < 20)
  {
      unsigned q1 = qubit1(generator);
      unsigned q2 = qubit2(generator);
      if (q1 != q2)
      {
          qpair.push_back(std::make_pair(q1, q2));
          i++;
      }
  }

  iqs::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(7777);
  iqs::QubitRegister<ComplexDP> psi_B(num_qubits, "base", 0 );
  psi_B.SetRngPtr(&rnd_generator);
  psi_B.Initialize("rand",1);

  iqs::QubitRegister<ComplexDP> psi_C(psi_B);

  {
    // no specialization
    psi_B.EnableStatistics();
    for (auto &p : qpair) {
       psi_B.ApplyControlled1QubitGate(p.first, p.second, G);
    }
    
    for(int pos = 0; pos < num_qubits; pos++) {
       psi_B.Apply1QubitGate(pos, G);
    }
    psi_B.GetStatistics();
  }

  {
    // with specialization
    psi_C.TurnOnSpecialize();
    psi_C.EnableStatistics();
    for (auto &p : qpair) {
       psi_C.ApplyControlled1QubitGate(p.first, p.second, G);
    }

    for(int pos = 0; pos < num_qubits; pos++) {
       psi_C.Apply1QubitGate(pos, G);
    }
    psi_C.GetStatistics();
  }

  double e = psi_C.MaxAbsDiff(psi_B);
  if (myrank == 0)
      printf("Max abs difference between the entries of |B> and |C>:\n e = %lf\n", e);

}
