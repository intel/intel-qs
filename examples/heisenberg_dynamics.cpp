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

#include "../include/qureg.hpp"

using namespace std;

#include <complex>
#include <iomanip>
#include <iostream>


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


/* 
Dynamics of a transverse Heisenberg model. Single first-order Trotter step.
Hamiltonian:
H = -J_z \sum Z_j Z_{j+1} - g J_z \sum X_j


-----Rx---c--------c-----------------------------------------------     
          |        |                 
----------X---Rz---X---Rx---c--------c-----------------------------                       
                            |        |       
----------------------------X---Rz---X---Rx---c--------c-----------                                  
                                              |        |    
----------------------------------------------X---Rz---X---Rx------

In the code below, the initial condition is set such that spins
are rotated in X-Z plane, with the first spin 180 degrees counter
to the others.
*/



int main(int argc, char **argv)
{

/////////////////////////////////////////////////////////////////////////////////////////
// Setting the MPI environment
/////////////////////////////////////////////////////////////////////////////////////////
  
  unsigned myrank=0, nprocs=1;
  iqs::mpi::Environment env(argc, argv);
  if (env.IsUsefulRank()==false) return 0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();

  int num_qubits;

  std::size_t tmp_size = 0;
  if(argc != 2)
  {
    fprintf(stderr, "usage: %s <num_qubits> \n", argv[0]);
    exit(1);
  }
  else
  {
    num_qubits = atoi(argv[1]);
    fprintf(stderr, "Number of qubits: %i \n",num_qubits);
  }
  
  double expectval = 0.;
  
  
  iqs::QubitRegister<ComplexDP> psi(num_qubits,"base",1);
  // 1 means zeroth qubit is flipped. 0 means none are flipped.
  psi.EnableStatistics();

  // Set model parameters
  double J = 1.;
  double g = 1.0;
  double dt = 0.1;
 
  // Initial state before Trotter step:
  // Rotated slightly in X-Z plane, with first spin 180deg from others
  for(int i=0;i<num_qubits;i++)
  {
    psi.ApplyRotationY(i,3.14159/6.);
  }
  
  // Expectation value on each qubit, before simulation
  for(int i=0;i<num_qubits;i++)
  {
    expectval = 0.;
    expectval = psi.ExpectationValueZ(i);
    printf("<Z> on qubit %d: %.12f \n",i,expectval);
  }
  printf("\n");
  
  printf("Running one Trotter step.\n");
  // Apply one Trotter step of Hamiltonian
  for(int i=0;i<num_qubits-1;i++)
  {
    psi.ApplyRotationX(i,g*J*dt);
    
    psi.ApplyCPauliX(i,i+1);
    psi.ApplyRotationZ(i+1,J*dt);
    psi.ApplyCPauliX(i,i+1);
  }
  
  // Last qubit must also have an Rx rotation
  psi.ApplyRotationX(num_qubits-1,g*J*dt);
  
  // Expectation value on each qubit
  for(int i=0;i<num_qubits;i++)
  {
    expectval = 0.;
    expectval = psi.ExpectationValueZ(i);
    printf("<Z> on qubit %d: %.12f \n",i,expectval);
  }

}

