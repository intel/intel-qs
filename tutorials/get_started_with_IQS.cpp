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
/// Tutorial on the basic use of Intel Quantum Simulator (IQS).

// Include the header file with the declaration of all classes and methods of IQS.
#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

// Start of the main program (C++ language).
int main(int argc, char **argv)
{
#ifndef INTELQS_HAS_MPI
  std::cout << "\nThis introductory code is thought to be run with MPI.\n"
            << "To do so, please set the option '-DIqsMPI=ON' when calling CMake.\n\n"
            << "However the code will execute also without MPI.\n\n";
#endif


/////////////////////////////////////////////////////////////////////////////////////////
// Setting the MPI environment
/////////////////////////////////////////////////////////////////////////////////////////

  // Create the MPI environment, passing the same argument to all the ranks.
  qhipster::mpi::Environment env(argc, argv);
  // IQS is structured so that only 2^k ranks are used to store and manipulate
  // the quantum state. In case the number of ranks differ from a power of 2,
  // all ranks in excess are effectively excluded from the computations and called
  // 'dummy'. The dummy ranks should be terminated.
  if (env.IsUsefulRank() == false) return 0;
  // IQS has functions that simplify some MPI instructions. However, it is important
  // to keep trace of the current rank.
  int myid = env.GetStateRank();

  // NOTE: Above, we asked for the 'state rank', meaning that we are considering the MPI
  // ranks that are used to store and manipulate a single quantum state. This difference
  // is unnecessary when simulating ideal quantum circuits, but once noise is introduced
  // via the noise-gate approach the situation may differ.

/////////////////////////////////////////////////////////////////////////////////////////
// Initialize the state of the quantum register
/////////////////////////////////////////////////////////////////////////////////////////
/* IQS stores a full representation of the quantum state in the computational basis.
 * In practice, the quantum state of N qubits is represented as a complex vector with
 * 2^N components.
 *
 * Each componenet corresponds to the probability amplitude of a specific computational
 * basis state:
 *     ψ(k)=⟨k|ψ⟩
 * with the index k corresponding to the N-bit integer in decimal representation, and
 * k∈{0,1,2,…,2N−1}.
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // Allocate memory for the quantum register's state and initialize it to |0000>.
  // This can be achieved by using the codeword "base".
  int num_qubits = 4;
  QubitRegister<ComplexDP> psi (num_qubits);
  std::size_t index = 0;
  psi.Initialize("base", index);

  // The state can be initialized to a random state. To allow such initialization,
  // we need to declare a (pseudo) random number generator...
  qhipster::RandomNumberGenerator<double> rng;
  // ... and associate it (via pointer) to the quantum state.
  psi.SetRngPtr(&rng);
  // Set the seed:
  int rng_seed = 7777;
  rng.SetSeedStreamPtrs( rng_seed );

  // NOTE: the random number generator is able to generate three different kinds
  //       of random numbers:
  //       *local* --> different for each pool rank
  //       *state* --> common to all ranks of the same state
  //       *pool*  --> common to all ranks of the pool

  // Initialize the state to a random state, this can be achieved with the codeword "rand"
  // followed by 0 if we desire to use *local* random numbers (this speed up the process
  // of generating the random numbers).
  psi.Initialize("rand", 0);

/////////////////////////////////////////////////////////////////////////////////////////
// Display the quantum state
/////////////////////////////////////////////////////////////////////////////////////////
/* It is important to be able to access and visualize the quantum state.
 * IQS allows to access the single componenets of the state or to print a comprehensive
 * description.
 * What index is associated to state |1011⟩? In decimal representation one has:
 *     1011 → 1×2^0 + 0×2^1 + 1×2^2 + 1×2^3 = 1+4+8 = 13
 * therefore it corresponds to the computational basis state with index 13.
 *
 * NOTE: contrary to what is adopted in decimal notation, our binary representation
 * must be read from left to right (from least significant to most significant bit).
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // Initialize the state to |1000>.
  // The index of |1000> in decimal representation is 1.
  index = 1;
  psi.Initialize("base", index);

  // There are 2^4=16 amplitudes, corresponding to |0000>, |1000>, |0100>, |1100> ... |1111>
  std::size_t num_amplitudes = UL(1L << UL(num_qubits));
  ComplexDP amplitude;
  std::stringstream buffer;
  buffer << "\nExplicit list of all amplitudes of |1000>:\n";
  for (index=0; index<num_amplitudes; ++index)
  {
      amplitude = psi.GetGlobalAmplitude(index);
      buffer << "\tpsi(" << index << ") = <" << index << "|psi> = " << amplitude << "\n";
  }

  // Print to screen. The second argument is set to true if all ranks need to print.
  // If it is set to false, then only the main rank prints to screen.
  qhipster::mpi::Print(buffer.str(),false);

  // A complete description of the state is provided by the QubitRegister method Print().
  // One can provide a brief label to the state description.
  std::string label;
  label = "Computational basis state |1000>";
  psi.Print(label);
  if (myid==0) std::cout << std::endl;

/////////////////////////////////////////////////////////////////////////////////////////
// One-qubit gates
/////////////////////////////////////////////////////////////////////////////////////////
/* In the gate-model of quantum computation, one manipulates the quantum state
 * by means of unitary transformations acting on one or two qubits.
 * Let us apply a few of the standard one-qubit gates.
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // State is |1000>.
  // Flip qubit 1 by applying the Pauli X gate: |1000> ==> |1100>
  psi.ApplyPauliX(1);

  // Display all amplitudes.
  psi.Print("Computational basis state |1100>");
  if (myid==0) std::cout << std::endl;
    
  // Apply the Hadamard gate on qubit 0: |1100> ==> |-100> ~ |0100>-|1100>
  psi.ApplyHadamard(0);

  // Apply Pauli Z gate on qubit 1:  |-100> ==> -|-100>
  int qubit = 1;
  psi.ApplyPauliZ(qubit);

  // Apply Pauli X gate on qubit 0: -|-100> ==>  |-100>
  psi.ApplyPauliX(0);

/////////////////////////////////////////////////////////////////////////////////////////
// Two-qubit gates
/////////////////////////////////////////////////////////////////////////////////////////
/* To achieve universal quantum computation, it is enought to implement one-qubit
 * gates and a single type of two-qubit gate. The essential requirement is that such
 * two-qubit gate is able to generate entanglement. Usually the controlled-not gate
 * (CNOT in the following) is the operation of choice.
 * IQS provides built-in methods to implement a much broader variety of two-qubit gates.
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // Currently, state is |-100>.
  // Apply a CNOT(1,0): flip qubit 0 conditioned on the state of qubit 1.
  // |-100> ==> -|-100>
  int control_qubit = 1;
  int target_qubit = 0;
  psi.ApplyCPauliX(control_qubit, target_qubit);

  // The application of the previous CNOT did not create any entanglement.
  // This is achieved by exchanging the role of control and target qubits.
  // Apply a CNOT(0,1): flip qubit 1 conditioned on the state of qubit 0.
  // -|-100> ~ -|0100>+|1100> ==> -|0100>+|1000>
  control_qubit = 0;
  target_qubit = 1;
  psi.ApplyCPauliX(control_qubit, target_qubit);

  // Display all amplitudes.
  psi.Print("Entangled state ~ -|0100>+|1000>");
  if (myid==0) std::cout << std::endl;

/////////////////////////////////////////////////////////////////////////////////////////
// Custom gates
/////////////////////////////////////////////////////////////////////////////////////////
/* If IQS does not provide the gates needed in your circuit, it is possible to
 * implement custom one-qubit gates and controlled gates.
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // Define an arbitrary single qubit gate.
  // The quantum gate G is given by a 2x2 unitary matrix, here using a custom made class
  // defined inside the IQS library.
  TM2x2<ComplexDP> G;
  G(0,0) = { 0.592056606032915,   0.459533060553574};
  G(0,1) = {-0.314948020757856, - 0.582328159830658};
  G(1,0) = { 0.658235557641767, + 0.070882241549507};
  G(1,1) = { 0.649564427121402, + 0.373855203932477};

  // To verify that G is unitary, we will compute the norm of psi before and after
  // the application of G. This is a necessary but not sufficient condition for the
  // unitarity of G.
  double initial_norm, final_norm;
  initial_norm = psi.ComputeNorm();
  double accepted_error = 1e-15;
  if ( std::abs(initial_norm-1.) > accepted_error )
    if (myid==0)
        std::cout << "Even before the application of G, state psi had normalization "
                  << initial_norm << "\n";

  // Apply the custom gate G to qubit 0.
  qubit = 0;
  psi.Apply1QubitGate(qubit,G);
  final_norm = psi.ComputeNorm();
  qhipster::mpi::Print("State |psi> was multiplied by a 2x2 unitary matrix G.",false);
  if (initial_norm != final_norm && myid == 0)
      std::cout << "The application of G changed the norm of state psi: from "
                << initial_norm << " to " << final_norm << "\n\n";
  else if (myid==0)
      std::cout << "Sanity check: norm was unchanged by G.\n\n";

  // It is also possible to apply the arbitrary gate specified by G conditioned
  // on the state of another qubit. G is applied only when the control qubit is in |1>.
  control_qubit = 1;
  target_qubit = 0;
  psi.ApplyControlled1QubitGate(control_qubit, target_qubit, G);

/////////////////////////////////////////////////////////////////////////////////////////
// Single-qubit measurements
/////////////////////////////////////////////////////////////////////////////////////////
/* To extract information from the quantum register, one can obtain the probability
 * of measuring a certain qubit in the computational basis and obtaining the outcome
 * "1" (i.e. the state being in |1⟩).
 * Once the probability is known, one can draw a random number to simulate the
 * stochastic outcome of the measurement and collapse the wavefunction accordingly.
 *
 * NOTE: Computing the probability of a certain outcome does not collapse automatically
 * the wavefunction. This is helpful when the probabilities of multiple measurements
 * have to be computed without re-executing the quantum simulation.
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // Compute the probability of observing qubit 1 in state |1>.
  int measured_qubit = 1;
  double probability = psi.GetProbability(measured_qubit);
  // The expectation value of <Z1> can be computed from the above probability
  double expectation = -1*probability + 1*(1-probability);

  if (myid==0)
      std::cout << "The probability that qubit " << measured_qubit
                << " is in state |1> is " << probability
                << " \t<== should be 0.5\n"
                << "The Von Neumann measurement can be simulated by generating "
                << "the outcome at random and then collapsing the state accordingly:\n";

  // Draw random number in [0,1)
  double r;
  r = 0.66;
//  r = np.random.rand() //FIXME
  if (r < probability)
  {
      // Collapse the wavefunction according to qubit 1 being in |1>.
      if (myid==0)
          std::cout << "Simulated outcome is 1. Collapse the function accordingly.\n";
      psi.CollapseQubit(measured_qubit,true);
  }
  else
  {
      // Collapse the wavefunction according to qubit 1 being in |0>
      if (myid==0)
          std::cout << "Simulated outcome is 0. Collapse the function accordingly.\n";
      psi.CollapseQubit(measured_qubit,false);
  }

  // In both cases one needs to re-normalize the wavefunction:
  psi.Normalize();

  probability = psi.GetProbability( measured_qubit );
  if (myid==0)
      std::cout << "After collapse, the probability that qubit " << measured_qubit
                << " is in state |1> is " << probability
                << " \t<== should be " << (r<probability?"1":"0") << "\n";

/////////////////////////////////////////////////////////////////////////////////////////
// Expectation value of products of Pauli matrices
/////////////////////////////////////////////////////////////////////////////////////////
/* To extract information from the quantum register, one can obtain the expectation
 * value of Pauli strings. For example, consider the Pauli string given by:
 *     X0⊗id1⊗Z2⊗Z3 
 * such observable is defined by:
 * - the position of the non-trivial Pauli matrices, in this case {0,2,3}
 * - the corresponding Pauli matrices ( XX =1,  YY =2,  ZZ =3).
 * To facilitate the verification of the expectation value, we reinitialize
 * the quantum state to |+−01⟩ .
 *
 * We will also consider the Pauli string:
 *     X0⊗id1⊗Z2⊗Y3
 */
/////////////////////////////////////////////////////////////////////////////////////////

  // Prepare the state |+-01>.
  num_qubits = 4;
  index = 0;
  psi.Initialize("base", index);
  psi.ApplyPauliX(1);
  psi.ApplyPauliX(3);
  psi.ApplyHadamard(0);
  psi.ApplyHadamard(1);
  if (myid==0)
      std::cout << "\nRe-initialize |psi> to state |+-01> and compute "
                << "a few expectation values of Pauli strings.\n";

  // The Pauli string given by:  X_0 . id_1 . Z_2 . Z_3
  // Such observable is defined by the position of the non-trivial Pauli matrices:
  std::vector<unsigned> qubits_to_be_measured = {0,2,3};
  // And by the corresponding Pauli matrices (X=1, Y=2, Z=3)
  std::vector<unsigned> observables = {1,3,3};

  // The expectation value <psi|X_0.id_1.Z_2.Z_3|psi> is obtained via:
  expectation = psi.ExpectationValue(qubits_to_be_measured, observables, 1.);
  if (myid==0)
      std::cout << "Expectation value <psi|X_0.id_1.Z_2.Z_3|psi> = " << expectation
                << " \t<== should be -1\n";

  // The expectation value <psi|X_0.id_1.Z_2.Y_3|psi> is obtained via:
  observables = {1,3,2};
  expectation = psi.ExpectationValue(qubits_to_be_measured, observables, 1.);
  if (myid==0)
      std::cout << "Expectation value <psi|X_0.id_1.Z_2.Y_3|psi> = " << expectation
                << " \t<== should be  0\n";

/////////////////////////////////////////////////////////////////////////////////////////

  // Quantum states and MPI environment are automatically destroyed at the return.
  return 0;
}
