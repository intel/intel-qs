//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

#include "qureg.hpp"
#include "highperfkernels.hpp"

/* 
 * Implements a Toffoli gate courtesy of Sonika Johri.
 *
 * Parameters:
 * 	psi - the wave function to apply the Toffoli gate to.
 * 	  qubit1 - target qubit being NOT'ed
 * 	  qubit2 - control qubit 1
 * 	  qubit3 - control qubit 2
 */
template<typename Type>
void QbitRegister<Type>::applyToffoli(unsigned const qubit1, 
		                        unsigned const qubit2, 
					  unsigned const qubit3) {

  openqu::TinyMatrix<Type, 2, 2, 32> V;
  V(0,0)={1.0/2.0,-1.0/2.0};
  V(0,1)={1.0/2.0,1.0/2.0};
  V(1,0)={1.0/2.0,1.0/2.0};
  V(1,1)={1.0/2.0,-1.0/2.0};

  openqu::TinyMatrix<Type, 2, 2, 32> V_dag;
  V_dag(0,0)={1.0/2.0,1.0/2.0};
  V_dag(0,1)={1.0/2.0,-1.0/2.0};
  V_dag(1,0)={1.0/2.0,-1.0/2.0};
  V_dag(1,1)={1.0/2.0,1.0/2.0};

  applyControlled1QubitGate(qubit2, qubit1, V);
  applyCPauliX(qubit3,qubit2);
  applyControlled1QubitGate(qubit2, qubit1, V_dag);
  applyCPauliX(qubit3,qubit2);
  applyControlled1QubitGate(qubit3, qubit1, V);

  return;
}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;

