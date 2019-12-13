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

#ifndef QAOA_EXTRA_FEATURES_H
#define QAOA_EXTRA_FEATURES_H

#include <vector>

#include "../../qureg/qureg.hpp"

namespace qaoa
{

/////////////////////////////////////////////////////////////////////////////////////////
// This function uses the IQS vector not as a quantum state but as a vector of size
// 2^num_qubits. This function initializes such vector to represent the diagonal of
// a classical Hamiltonian. Specifically, the cost function is that of the MaxCut problem
// based on the graph instance passed in input via its adjacency matrix.
// The graph has undirected edges of uniform weigth.
//
// It returns the max number of edges that can be cut.
/////////////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  int InitializeVectorAsMaxCutCostFunction(QubitRegister<Type> & diag,
                                           std::vector<int> & adjacency);
  
/////////////////////////////////////////////////////////////////////////////////////////
// Implement the QAOA operation:
//     |psi> <== exp(-i C gamma) |psi>
// for classical Hamiltonian C.
// The diagonal of C is provided as a (not-normalized) IQS vector.
/////////////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  void ImplementQaoaLayerBasedOnCostFunction(QubitRegister<Type> & psi,
                                             QubitRegister<Type> & diag,
                                             typename QubitRegister<Type>::BaseType gamma);

/////////////////////////////////////////////////////////////////////////////////////////
// Return the expectation value of the cost function.
/////////////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  typename QubitRegister<Type>::BaseType
  GetExpectationValueFromCostFunction(const QubitRegister<Type> & psi,
                                      const QubitRegister<Type> & diag);

/////////////////////////////////////////////////////////////////////////////////////////
// Return the histogram with the probability of measuring a graph coloring
// associated to a specific value of the cost function.
/////////////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  std::vector<typename QubitRegister<Type>::BaseType>
  GetHistogramFromCostFunction(const QubitRegister<Type> & psi,
                               const QubitRegister<Type> & diag, int max_value);

/////////////////////////////////////////////////////////////////////////////////////////

}		// namespace 'qaoa'

#endif		// end of guards QAOA_EXTRA_FEATURES_H
