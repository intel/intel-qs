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

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
#ifdef INTELQS_HAS_MPI
  openqu::mpi::Environment env(argc, argv);
  if (env.is_usefull_rank() == false) return 0;
  unsigned myrank = env.rank();
#else
  unsigned myrank=0;
#endif

  unsigned num_qubits = 3;
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

  
  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};


  // generate random pairs for control qubits
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


  QubitRegister<ComplexDP> psi1(num_qubits, "rand", -1);
  QubitRegister<ComplexDP> psi2(num_qubits, "rand", -1);
  QubitRegister<ComplexDP> psi3(num_qubits, "rand", -1);
  {
    // no specialization
    psi1.EnableStatistics();
    for (auto &p : qpair) {
       psi1.ApplyControlled1QubitGate(p.first, p.second, G);
    }
    
    for(int pos = 0; pos < num_qubits; pos++) {
       psi1.Apply1QubitGate(pos, G);
    }
    psi1.GetStatistics();
  }

  {
    // with specialization
    psi2.TurnOnSpecialize();
    psi2.EnableStatistics();
    for (auto &p : qpair) {
       psi2.ApplyControlled1QubitGate(p.first, p.second, G);
    }

    for(int pos = 0; pos < num_qubits; pos++) {
       psi2.Apply1QubitGate(pos, G);
    }
    psi2.GetStatistics();
  }

  double e = psi2.maxabsdiff(psi1);
  if (myrank == 0) printf("e = %lf\n", e);
}
