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

/// \addtogroup qureg
/// @{

/// @file qureg_applytoffoli.cpp
/// @brief Define the @c QubitRegister method for the application of the Toffoli gate.

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply the Toffoli gate.
/// @param control_1 index of the 1st control qubit
/// @param control_2 index of the 2nd control qubit
/// @param target index of the target qubit
///
/// Implemented by decomposing the Toffoli gate in 5 two-qubit gates.
/// Courtesy of Sonika Johri.
/// If both control qubits are in state |1\>, then the target qubit is flipped.
template<typename Type>
void QubitRegister<Type>::ApplyToffoli(unsigned const control_1, 
                                       unsigned const control_2, 
                                       unsigned const target)
{
  qhipster::TinyMatrix<Type, 2, 2, 32> V;
  V(0,0)={1.0/2.0,-1.0/2.0};
  V(0,1)={1.0/2.0,1.0/2.0};
  V(1,0)={1.0/2.0,1.0/2.0};
  V(1,1)={1.0/2.0,-1.0/2.0};

  qhipster::TinyMatrix<Type, 2, 2, 32> V_dag;
  V_dag(0,0)={1.0/2.0,1.0/2.0};
  V_dag(0,1)={1.0/2.0,-1.0/2.0};
  V_dag(1,0)={1.0/2.0,-1.0/2.0};
  V_dag(1,1)={1.0/2.0,1.0/2.0};

  ApplyControlled1QubitGate( control_1, target, V );
  ApplyCPauliX( control_2, control_1 );
  ApplyControlled1QubitGate( control_1, target, V_dag );
  ApplyCPauliX( control_2, control_1 );
  ApplyControlled1QubitGate( control_2, target, V );

  return;
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
