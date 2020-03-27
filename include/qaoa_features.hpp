#ifndef QAOA_EXTRA_FEATURES_HPP
#define QAOA_EXTRA_FEATURES_HPP

#include <vector>

#include "qureg.hpp"

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
// Return the expectation value of the (cost function)^2.
/////////////////////////////////////////////////////////////////////////////////////////

  template<typename Type>
  typename QubitRegister<Type>::BaseType
  GetExpectationValueSquaredFromCostFunction(const QubitRegister<Type> & psi,
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

#endif		// end of guards QAOA_EXTRA_FEATURES_HPP
