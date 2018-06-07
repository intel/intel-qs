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
#ifdef NOREPA_HAS_MPI
  openqu::mpi::Environment env(argc, argv);
  if (env.is_usefull_rank() == false) return 0;
  unsigned myrank = env.rank();
#else
  unsigned myrank=0;
#endif

  int nthreads = 88;
  glb_affinity.set_thread_affinity(nthreads);

/// --- PARAMETERS ------------------------------------------- ///
  int num_qubits = 8;
  int num_gates = 1;
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

  
  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};

  QubitRegister<ComplexDP> psi1(num_qubits, "rand", -1);

  // with specialization
  psi1.TurnOnSpecialize();
  for(int pos = 0; pos < num_qubits; pos++)
  {
      if (myrank == 0) printf(" ---------------------------------- \n");
      psi1.EnableStatistics();
      psi1.Apply1QubitGate(pos, G);
      psi1.GetStatistics();
      psi1.ResetStatistics();
  }
}
