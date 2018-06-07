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

/* Test for qHiPSTER:
 * proper implementaiton of the methods to compute the expectation value of Pauli strings
 * on single, two or multiple qubitsi.
 */

#include "../qureg/qureg.hpp"

using namespace std;

#include <iostream>	// to use: std::cout, std::cin and std::endl
#include <iomanip>	// to use: setw() in making tables
#include <complex>


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  unsigned myrank=0, nprocs=1;
#ifdef NOREPA_HAS_MPI
  openqu::mpi::Environment env(argc, argv);
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
//  MPI_rank(MPI_COMM_WORLD, &myrank);
#endif

  double sum = 0.;

  std::cout << "------------------\n"
            << "   Single qubit   \n"
            << "------------------\n";

  QubitRegister<ComplexDP> psi(1,"base",1);
  psi.EnableStatistics();  
  psi.ApplyHadamard(0);

  psi.Print(" initial state |psi>=|-> : ");

  psi.ExpectationValueX(0,sum);
  std::cout << "<psi|X|psi> = " << sum << "\n";

  psi.Print(" current state should still be |psi>=|-> : ");

// using the general method
  std::vector<unsigned> qubits(1,0);
  std::vector<unsigned> observables(1,1);
  sum=0;
  psi.ExpectationValue(qubits,observables,sum);
  std::cout << " general method  __  <psi|X|psi> = " << sum << "\n";

  psi.Print(" current state should still be |psi>=|-> : ");

  psi.ApplyPauliZ(0);
  psi.Print(" current state should be Z|psi>=|+> : ");

  sum=0;
  psi.ExpectationValueX(0,sum);
  std::cout << " <psi|X|psi> = " << sum << "\n";
  sum=0;
  psi.ExpectationValue(qubits,observables,sum);
  std::cout << " general method  __  <psi|X|psi> = " << sum << "\n";

  std::cout << "\n ------ repeat the experiment 100 times ---- \n";

  for (int i=0; i<100; i++)
      psi.ExpectationValueX(0,sum);
  sum=0;
  psi.ExpectationValueX(0,sum);
  std::cout << "after 100 expectation evaluations <psi|X|psi> = " << sum << "\n";

  std::cout << std::endl;
  psi.GetStatistics(); 

  if (myrank==0) std::cout << " goodbye \n" << std::endl;

  std::cout << "------------------\n"
            << "     4 qubits     \n"
            << "------------------\n";

  QubitRegister<ComplexDP> phi(4,"base",0);
  phi.ApplyPauliX(1);
  phi.ApplyHadamard(2);
  phi.ApplyPauliX(3);
  phi.ApplyHadamard(3);

  phi.Print(" initial state |phi> = |0> |1> |+> |-> = |01+-> :\n");

  qubits.assign({0,2});
  observables.assign({3,1});
  sum=0;
  phi.ExpectationValue(qubits,observables,sum);
  std::cout << " <phi|Z_0 X_2|phi> = " << sum << "\n";

  qubits.assign({0,2});
  observables.assign({3,2});
  sum=0;
  phi.ExpectationValue(qubits,observables,sum);
  std::cout << " <phi|Z_0 Y_2|phi> = " << sum << "\n";

  qubits.assign({1,2});
  observables.assign({3,1});
  sum=0;
  phi.ExpectationValue(qubits,observables,sum);
  std::cout << " <phi|Z_1 X_2|phi> = " << sum << "\n";

  qubits.assign({0,3});
  observables.assign({3,1});
  sum=0;
  phi.ExpectationValue(qubits,observables,sum);
  std::cout << " <phi|Z_0 X_3|phi> = " << sum << "\n";

  return 0;
}

