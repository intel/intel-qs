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


// Include the header file including the declaration of the QbitRegister class
// and of its methods.
#include "qureg/qureg.hpp"


// Start of the main program (C++ language).
int main(int argc, char **argv)
{
  // Create the MPI environment, passing the same argument to all the ranks.
  openqu::mpi::Environment env(argc, argv);
  // qHiPSTER is structured so that only even number of ranks are used to store
  // and manipulate the quantum state. In case the number of ranks is not supported,
  // try to decrease it by 1 until it is.
  if (env.is_usefull_rank() == false) return 0;
  // qHiPSTER has functions that simplify some MPI instructions. However, it is important
  // to keep trace of the current rank.
  int myid = env.rank();



  // The number of qubits in the quantum register is provided by an integer variable.
  // Notice that qubits will have indices: 0, 1, 2, ..., N-1.
  int N = 4;
  std::size_t tmpSize = 0;

  //// FIRST EXAMPLE ////
  {
      if (myid == 0)
      {
          printf("\n --- FIRST EXAMPLE --- \n");
          std::cout << "The 4-qubit register is initialized in state |0000>. \n"
                    << "A sequence of gates is then performed: a Pauli X gate on qubit 0; \n"
                    << "followed by a Hadamard gate on each other qubit (namely qubits 1, 2, \n"
                    << "and 3); then a CNOT gate (control qubit is 2, target qubit is 1). \n"
                    << "Finally qubit 1 is measured in the computational basis.\n";
      }

      // Create the state of a quantum register, having N qubits.
      // The state is initialized as a computational basis state (using the keyword "base")
      // corresponding to the index 0. The index corresponds to a N-bit integer in decimal
      // representation. With N qubits there are 2^N indices, from 0 to 2^{N-1}.
      QbitRegister<ComplexDP> psi1(N, "base", 0);

      // Let us apply a Pauli X gate on qubit 0, effectively flipping it from |0> to |1>.
      psi1.applyPauliX(0);

      // Let us apply an Hadamard gate on all other qubits.
      for (unsigned q=1; q<N; ++q)
      {
          psi1.applyHadamard(q);
      }

      // Two qubit gates are applied in a similar way. For example, a C-NOT between qubit 2
      // (control qubit) and qubit 1 (target qubit):
      unsigned control_q = 2;
      unsigned target_q = 1;
      psi1.applyCPauliX( control_q , target_q );

      // To extract information from the quantum register, one can obtain the probability of
      // measuring a certain qubit, here qubit 1, in the Z basis and obtaining the outcome "-1".
      // This corresponds to the probability of qubit 1 being in the |1> state.
      unsigned measured_qubit = 1;
      double prob = 0;
      prob = psi1.getProbability( measured_qubit );

      // Print such probability to screen, only if MPI rank=0.
      // This is done to avoid each rank to write the same information to screen.
      if (myid == 0)
          printf("The probability that qubit %d is in state |1> is %g\n", measured_qubit, prob);
  }


  //// SECOND EXAMPLE ////
  {
      if (myid == 0)
      {
          printf("\n --- SECOND EXAMPLE --- \n");
          std::cout << "The 4-qubit register is initialized in a random state. \n"
                    << "A sequence of gates is then performed: a Pauli X gate on qubit 0; \n"
                    << "followed by a Hadamard gate on each other qubit (namely qubits 1, 2, \n"
                    << "and 3); then a user-defined single-qubit gate on qubit 3; then the \n"
                    << "same user-defined gate is applied to qubit 2 controlled by qubit 1. \n"
                    << "Finally the Pauli string X_0 . Z_2 . Y_3 is jointly measured.\n";
      }
      // Create the state of a quantum register, having N qubits.
      // The state is initialized as a random state (using the keyword "rand"):
      // This requires a random number generator (RNG), that we initialize just before the second
      // register. Notice that '777' plays the role of the seed to initialize the RNG.and the seed
      std::default_random_engine generator;
      QbitRegister<ComplexDP> psi2(N, "rand", 777);

      // Let us apply a X Pauli gate on qubit 0, effectively flipping it from |0> to |1>.
      psi2.applyPauliX(0);

      // Let us apply an Hadamard gate on all other qubits.
      for (unsigned q=1; q<N; ++q)
      {
          psi2.applyHadamard(q);
      }

      // One can define an arbitrary single qubit gate and apply it to the chosen qubit.
      // The quantum gate G is given by a 2x2 unitary matrix:
      TM2x2<ComplexDP> G;
      G(0, 0) = {0.592056606032915, 0.459533060553574}; 
      G(0, 1) = {-0.314948020757856, -0.582328159830658};
      G(1, 0) = {0.658235557641767, 0.070882241549507}; 
      G(1, 1) = {0.649564427121402, 0.373855203932477};

      unsigned chosen_qubit = 3;
      psi2.apply1QubitGate(chosen_qubit, G);

      // It is also possible to apply the arbitrary gate specified by G controlled on the state
      // of another qubit. G is applied only when control_qubit is in |0>.
      unsigned control_qubit = 1;
      chosen_qubit = 2;
      psi2.applyControlled1QubitGate( control_qubit, chosen_qubit, G);

      // To extract information from the quantum register, one can obtain the expectation value
      // of Pauli strings. For example, consider the Paulis tring given by:
      //     X_0 . id_1 . Z_2 , Y_3
      // Such observable is defined by the position of the non-trivial Pauli matrices:
      std::vector<unsigned> qubits_to_be_measured= {0,2,3};
      // And by the corresponding Pauli matrices (X=1, Y=2, Z=3)
      std::vector<unsigned> observables = {1,3,2};

      // The expectation value <psi2|X_0.id_1.Z_2.Y_3|psi2> is obtained via:
      double average = 0.;
      psi2.expectationValue(qubits_to_be_measured, observables, average);
      // Print the expectation value to screen.
      if (myid == 0)
      {
          printf("Expectation value <psi2| X_0 . id_1 . Z_2 . Y_3 |psi2> = %g\n\n", average);
      }
  }
}
