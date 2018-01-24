//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation 
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

#include "qureg/qureg.hpp"

// The scope is applying a sequence of Ngates=40 gates to the quantum register.
// The form of each gate is the same:
//   controlled 1-qubit operation defined by the 2x2 matrix G
// but each pair (control,target) for the involved qubit is randomly generated.

int main(int argc, char **argv)
{
  openqu::mpi::Environment env(argc, argv);
  if (env.is_usefull_rank() == false) return 0;
  int myrank = env.rank();

  unsigned Nqubits = 3;
  std::size_t tmpSize = 0;
  if(argc != 2)
  {
     fprintf(stderr, "\nusage: %s <Nqubits> \n\n", argv[0]);
     exit(1);
  }
  else
  {
     Nqubits = atoi(argv[1]);
  }

  // number of gates
  unsigned Ngates=20;
  // number of repetition of the (stochastic) noisy circuit
  unsigned Ncirc=200;

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
  std::uniform_int_distribution<int> qubit1(0, Nqubits - 1);
  std::uniform_int_distribution<int> qubit2(0, Nqubits - 1);
  unsigned i = 0;
  while (i < Ngates) {
     unsigned q1 = qubit1(generator);
     unsigned q2 = qubit2(generator);
     if (q1 != q2)
     {
        qpair.push_back(std::make_pair(q1, q2));
        i++;
     }
  }


  QbitRegister<ComplexDP> psi0(Nqubits);
  NoisyQureg<ComplexDP> psi1(Nqubits, 100., 50.);
  NoisyQureg<ComplexDP> psi2(Nqubits, 40., 20.);
  psi0.Init("base", 0);


  // ---------------- no noise
  {
    psi0.EnbStat();
    for (auto &p : qpair)
        psi0.applyControlled1QubitGate(p.first, p.second, G);
    
    for(int pos = 0; pos < Nqubits; pos++)
        psi0.apply1QubitGate(pos, G);
    psi0.GetStat();
  }

  // ---------------- slow decoherence
  double over_sq_1 = 0.;
  for (unsigned j=0; j<Ncirc; j++)
  {
    psi1.Init("base", 0);
    psi1.reset_time_for_all_qubits();
    for (auto &p : qpair)
        psi1.applyControlled1QubitGate(p.first, p.second, G);
    
    for(int pos = 0; pos < Nqubits; pos++)
        psi1.apply1QubitGate(pos, G);

    psi1.apply_noise_gates_on_all_qubits();
    over_sq_1 = over_sq_1 + std::norm( psi0.compute_overlap(psi1) ) ;
  }
  over_sq_1 = over_sq_1/(double)Ncirc;

  // ---------------- fast decoherence
  double over_sq_2 = 0.;
  for (unsigned j=0; j<Ncirc; j++)
  {
    psi2.Init("base", 0);
    psi2.reset_time_for_all_qubits();
    for (auto &p : qpair)
        psi2.applyControlled1QubitGate(p.first, p.second, G);
    
    for (int pos = 0; pos < Nqubits; pos++)
        psi2.apply1QubitGate(pos, G);

    psi2.apply_noise_gates_on_all_qubits();
    over_sq_2 = over_sq_2 + std::norm( psi0.compute_overlap(psi2) ) ;
  }
  over_sq_2 = over_sq_2/(double)Ncirc;

  // ---------------- 
  // computation of the overlap between the ideal state and those exposed to noise
  if (myrank == 0) 
      std::cout << " Overlap-squared between ideal and 'slow decoherence' state = "
                << over_sq_1 << "\n"
                << " Overlap-squared between ideal and 'fast decoherence' state = "
                << over_sq_2 << "\n";

//  double e = psi2.maxabsdiff(psi1);
}
