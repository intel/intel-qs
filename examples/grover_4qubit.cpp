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


// Include the header file including the declaration of the QubitRegister class.
#include "../include/qureg.hpp"

// Start of the main program (C++ language).
int main(int argc, char **argv)
{
  // Create the MPI environment, passing the same argument to all the ranks.
  qhipster::mpi::Environment env(argc, argv);
  // qHiPSTER is structured so that only even number of ranks are used to store
  // and manipulate the quantum state. In case the number of ranks is not supported,
  // try to decrease it by 1 until it is.
  if (env.IsUsefulRank() == false) return 0;
  // qHiPSTER has functions that simplify some MPI instructions. However, it is important
  // to keep trace of the current rank.
  int myid = env.GetRank();


  // The number of qubits in the quantum register is provided by an integer variable.
  // Notice that qubits will have indices: 0, 1, 2, ..., N-1.
  int N = 4;
  std::size_t tmpSize = 0;

   //// Grover Search EXAMPLE ////

      if (myid == 0) printf("\n --- GROVER EXAMPLE --- \n");
      int Ngrover = 4;
      // Create the state of a quantum register, having N qubits.
      // The state is initialized as a computational basis state (using the keyword "base")
      // corresponding to the index 0. The index corresponds to a N-bit integer in decimal
      // representation. With N qubits there are 2^N indices, from 0 to 2^{N-1}.
      QubitRegister<ComplexDP> psig(Ngrover, "base", 0);
	  
      psig.Print("Initial State =");

      //////////////////////
      // Initialization  ///
      //////////////////////

      // Let us apply an Hadamard gate on all other qubits.
      for (unsigned q=0; q<Ngrover; ++q)
      {
          psig.ApplyHadamard(q);
      }
      psig.Print("After Hadamard =");
	  
      /////////////////////////////////////////////////////////
      // Oracle for the state |0>_0 . |1>_1 . |0>_2 . |0>_3 ///
      /////////////////////////////////////////////////////////

      psig.ApplyPauliX(0);
      psig.ApplyPauliX(2);
      psig.ApplyPauliX(3);
		
		
      double PI_4 = M_PI/4.0;
	  
      psig.ApplyCPhaseRotation(0,3,PI_4);
      psig.ApplyCPauliX(0,1);
      psig.ApplyCPhaseRotation(1,3,-PI_4);
      psig.ApplyCPauliX(0,1);
      psig.ApplyCPhaseRotation(1,3,PI_4);
      psig.ApplyCPauliX(1,2);
      psig.ApplyCPhaseRotation(2,3,-PI_4);
      psig.ApplyCPauliX(0,2);
      psig.ApplyCPhaseRotation(2,3,PI_4);
      psig.ApplyCPauliX(1,2);
      psig.ApplyCPhaseRotation(2,3,-PI_4);
      psig.ApplyCPauliX(0,2);
      psig.ApplyCPhaseRotation(2,3,PI_4);
  
      psig.ApplyPauliX(0);
      psig.ApplyPauliX(2);
      psig.ApplyPauliX(3);
		 
      /////////////////////
      // Amplification  ///
      /////////////////////
      for (unsigned q=0; q<Ngrover; ++q)
      {
          psig.ApplyHadamard(q);
      }
		
      for (unsigned q=0; q<Ngrover; ++q)
      {
          psig.ApplyPauliX(q);
      }
		
      /////////////////
      // CCCZ Gate  ///
      /////////////////
      
      psig.ApplyCPhaseRotation(0,3,PI_4);
      psig.ApplyCPauliX(0,1);
      psig.ApplyCPhaseRotation(1,3,-PI_4);
      psig.ApplyCPauliX(0,1);
      psig.ApplyCPhaseRotation(1,3,PI_4);
      psig.ApplyCPauliX(1,2);
      psig.ApplyCPhaseRotation(2,3,-PI_4);
      psig.ApplyCPauliX(0,2);
      psig.ApplyCPhaseRotation(2,3,PI_4);
      psig.ApplyCPauliX(1,2);
      psig.ApplyCPhaseRotation(2,3,-PI_4);
      psig.ApplyCPauliX(0,2);
      psig.ApplyCPhaseRotation(2,3,PI_4);

      /////////////////////
      // End CCCZ Gate  ///
      /////////////////////

      for (unsigned q=0; q<Ngrover; ++q)
      {
          psig.ApplyPauliX(q);
      }
      for (unsigned q=0; q<Ngrover; ++q)
      {
          psig.ApplyHadamard(q);
      }

      psig.Print("Measurement =");
}
