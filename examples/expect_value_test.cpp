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

// Test for qHiPSTER:
// proper implementation of the methods to compute the expectation value of Pauli strings
// on single, two or multiple qubitsi.

#include "../include/qureg.hpp"

using namespace std;

#include <complex>
#include <iomanip>	// to use: setw() in making tables
#include <iostream>	// to use: std::cout, std::cin and std::endl

/////////////////////////////////////////////////////////////////////////////////////////

// Definition of a utility macro to print only from the main rank.
#define MPIout \
if (iqs::mpi::Environment::GetStateRank()==0) std::cout

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  unsigned myrank=0, nprocs=1;
  iqs::mpi::Environment env(argc, argv);
  if (env.IsUsefulRank()==false) return 0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  if (nprocs > 1)
  {
      if (myrank==0)
          fprintf(stderr, "example to be launched with a single process\n");
      exit(1);
  }
  double expectation;

  MPIout << "------------------\n"
         << "   Single qubit   \n"
         << "------------------\n\n";

  iqs::QubitRegister<ComplexDP> psi(1,"base",1);
  psi.EnableStatistics();  
  psi.ApplyHadamard(0);

  psi.Print(" initial state |psi>=|-> : ");

  expectation = psi.ExpectationValueX(0);
  MPIout << "<psi|X|psi> = " << expectation << "\n\n";

  psi.Print(" current state should still be |psi>=|-> : ");

// using the general method
  std::vector<unsigned> qubits(1,0);
  std::vector<unsigned> observables(1,1);
  expectation = psi.ExpectationValue(qubits,observables);
  MPIout << " general method  __  <psi|X|psi> = " << expectation << "\n\n";

  psi.Print(" current state should still be |psi>=|-> : ");

  psi.ApplyPauliZ(0);
  MPIout << "\n";
  psi.Print(" current state should be Z|psi>=|+> : ");

  expectation = psi.ExpectationValueX(0);
  MPIout << " <psi|X|psi> = " << expectation << "\n";
  expectation = psi.ExpectationValue(qubits,observables);
  MPIout << " general method  __  <psi|X|psi> = " << expectation << "\n";

  std::cout << std::endl;
  psi.GetStatistics(); 

  MPIout << " goodbye \n" << std::endl;

  MPIout << "------------------\n"
         << "     4 qubits     \n"
         << "------------------\n\n";

  iqs::QubitRegister<ComplexDP> phi(4,"base",0);
  phi.ApplyPauliX(1);
  phi.ApplyHadamard(2);
  phi.ApplyPauliX(3);
  phi.ApplyHadamard(3);

  phi.Print(" initial state |phi> = |0> |1> |+> |-> = |01+-> :\n");

  qubits.assign({0,2});
  observables.assign({3,1});
  expectation = phi.ExpectationValue(qubits,observables);
  MPIout << " <phi|Z_0 X_2|phi> = " << expectation << "  <== should be  1\n";

  qubits.assign({0,2});
  observables.assign({3,2});
  expectation = phi.ExpectationValue(qubits,observables);
  MPIout << " <phi|Z_0 Y_2|phi> = " << expectation << "  <== should be  0\n";

  qubits.assign({1,2});
  observables.assign({3,1});
  expectation = phi.ExpectationValue(qubits,observables);
  MPIout << " <phi|Z_1 X_2|phi> = " << expectation << "  <== should be -1\n";

  qubits.assign({0,3});
  observables.assign({3,1});
  expectation = phi.ExpectationValue(qubits,observables);
  MPIout << " <phi|Z_0 X_3|phi> = " << expectation << "  <== should be -1\n";

  return 0;
}

