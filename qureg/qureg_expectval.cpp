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

#define __ONLY_NORMALIZED_STATES__

#include "qureg.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_expectval.cpp
///  @brief Define the @c QbitRegister methods related to expectation values of Pauli strings.

/// @brief Compute expectation value of Pauli X for qubit over the full-register state
/// @param qubit index of the involved qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --> sum + coeff * <psi| X_qubit |psi>
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueX(unsigned qubit, BaseType &sum, BaseType coeff)
{
// compute the iexpectation value <psi|X|psi> = <psi|H.Z.H|psi>
  applyHadamard(qubit);
  sum += coeff*(1. - 2.*getProbability(qubit));
// recover the initial state |psi>
  applyHadamard(qubit);
}


#if 0
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueX(unsigned qubit, BaseType &sum)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+X)|psi>
  openqu::TinyMatrix<Type, 2, 2, 32> x_plus_2id;
  x_plus_2id(0, 1) = x_plus_2id(1, 0) = Type(1., 0.);
  x_plus_2id(0, 0) = x_plus_2id(1, 1) = Type(2., 0.);
  apply1QubitGate(qubit, x_plus_2id);

// update sum adding <psi|X|psi>
  BaseType final_norm = computenorm();
  sum += (final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 2, 2, 32> inverse;
  inverse(0, 1) = inverse(1, 0) = Type(-1./3., 0.);
  inverse(0, 0) = inverse(1, 1) = Type( 2./3., 0.);
  apply1QubitGate(qubit, inverse);
}
#endif


/// @brief Compute expectation value of Pauli Y for qubit over the full-register state
/// @param qubit index of the involved qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --> sum + coeff * <psi| Y_qubit |psi>
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueY(unsigned qubit, BaseType &sum, BaseType coeff)
{
// G is matrix from change of basis Y --> Z , such that G^dagger.Z.G=Y
  openqu::TinyMatrix<Type, 2, 2, 32> G;
  BaseType f = 1. / std::sqrt(2.);
  G(0, 0) = G(1, 0) = Type(f, 0.);
  G(0, 1) = Type(0.,-f);
  G(1, 1) = Type(0., f);
// compute the iexpectation value <psi|Y|psi> = <psi|G^dagger.Z.G|psi>
  apply1QubitGate(qubit,G);
  sum += coeff*(1. - 2.*getProbability(qubit));
// recover the initial state |psi> by applying G^dagger
  G(0, 0) = G(0, 1) = Type(f, 0.);
  G(1, 0) = Type(0., f);
  G(1, 1) = Type(0.,-f);
  apply1QubitGate(qubit,G);
}


/// @brief Compute expectation value of Pauli Z for qubit over the full-register state
/// @param qubit index of the involved qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --> sum + coeff * <psi| Z_qubit |psi>
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueZ(unsigned qubit, BaseType &sum, BaseType coeff)
{
  sum += coeff*(1. - 2.*getProbability(qubit));
}


/// @brief Compute expectation value of Pauli X.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --> sum + coeff * <psi| X_qubit.X_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueXX(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+XX)|psi>
  openqu::TinyMatrix<Type, 4, 4, 32> xx_plus_2id;
  xx_plus_2id(0, 1) = xx_plus_2id(0, 2) = xx_plus_2id(1, 0) = xx_plus_2id(1, 3) = Type(0., 0.);
  xx_plus_2id(2, 0) = xx_plus_2id(2, 3) = xx_plus_2id(3, 1) = xx_plus_2id(3, 2) = Type(0., 0.);
  xx_plus_2id(0, 3) = xx_plus_2id(1, 2) = xx_plus_2id(2, 1) = xx_plus_2id(3, 0) = Type(1., 0.);
  xx_plus_2id(0, 0) = xx_plus_2id(1, 1) = xx_plus_2id(2, 2) = xx_plus_2id(3, 3) = Type(2., 0.);
  apply2QubitGate(qubit, qubit2, xx_plus_2id);

// update sum adding <psi|XX|psi>
  BaseType final_norm = computenorm();
  sum += coeff*(final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 4, 4, 32> inverse;
  inverse(0, 1) = inverse(0, 2) = inverse(1, 0) = inverse(1, 3) = Type(0., 0.);
  inverse(2, 0) = inverse(2, 3) = inverse(3, 1) = inverse(3, 2) = Type(0., 0.);
  inverse(0, 3) = inverse(1, 2) = inverse(2, 1) = inverse(3, 0) = Type(-1./3., 0.);
  inverse(0, 0) = inverse(1, 1) = inverse(2, 2) = inverse(3, 3) = Type( 2./3., 0.);
  apply2QubitGate(qubit, qubit2, inverse);
}


/// @brief Compute expectation value of Pauli Y.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Y_qubit.X_qubit2 |psi\>
//------------------------------------------------------------------------------
// TODO : check how the basis is it is: 00-01-10-11 or 00-10-01-11 !!!!
//        this code uses the standard 00-01-10-11 despite this being opposite to backend storage convention
template <class Type>
void QbitRegister<Type>::expectationValueYX(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+XY)|psi>
  openqu::TinyMatrix<Type, 4, 4, 32> xy_plus_2id;
  xy_plus_2id(0, 1) = xy_plus_2id(0, 2) = xy_plus_2id(1, 0) = xy_plus_2id(1, 3) = Type(0., 0.);
  xy_plus_2id(2, 0) = xy_plus_2id(2, 3) = xy_plus_2id(3, 1) = xy_plus_2id(3, 2) = Type(0., 0.);
  xy_plus_2id(0, 3) = xy_plus_2id(1, 2) = Type(0.,-1.);
  xy_plus_2id(2, 1) = xy_plus_2id(3, 0) = Type(0., 1.);
  xy_plus_2id(0, 0) = xy_plus_2id(1, 1) = xy_plus_2id(2, 2) = xy_plus_2id(3, 3) = Type(2., 0.);
  apply2QubitGate(qubit, qubit2, xy_plus_2id);

// update sum adding <psi|XY|psi>
  BaseType final_norm = computenorm();
  sum += coeff*(final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 4, 4, 32> inverse;
  inverse(0, 1) = inverse(0, 2) = inverse(1, 0) = inverse(1, 3) = Type(0., 0.);
  inverse(2, 0) = inverse(2, 3) = inverse(3, 1) = inverse(3, 2) = Type(0., 0.);
  inverse(0, 3) = inverse(1, 2) = Type(0., 1./3.);
  inverse(2, 1) = inverse(3, 0) = Type(0.,-1./3.);
  inverse(0, 0) = inverse(1, 1) = inverse(2, 2) = inverse(3, 3) = Type( 2./3., 0.);
  apply2QubitGate(qubit, qubit2, inverse);
}


/// @brief Compute expectation value of Pauli Z.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Z_qubit.X_qubit2 |psi\>
//------------------------------------------------------------------------------
// TODO : check how the basis is it is: 00-01-10-11 or 00-10-01-11 !!!!
//        this code uses the standard 00-01-10-11 despite this being opposite to backend storage convention
template <class Type>
void QbitRegister<Type>::expectationValueZX(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+XZ)|psi>
  openqu::TinyMatrix<Type, 4, 4, 32> xz_plus_2id;
  xz_plus_2id(0, 2) = xz_plus_2id(0, 3) = xz_plus_2id(1, 2) = xz_plus_2id(1, 3) = Type(0., 0.);
  xz_plus_2id(2, 0) = xz_plus_2id(2, 1) = xz_plus_2id(3, 0) = xz_plus_2id(3, 1) = Type(0., 0.);
  xz_plus_2id(0, 1) = xz_plus_2id(1, 0) = Type( 1., 0.);
  xz_plus_2id(2, 3) = xz_plus_2id(3, 2) = Type(-1., 0.);
  xz_plus_2id(0, 0) = xz_plus_2id(1, 1) = xz_plus_2id(2, 2) = xz_plus_2id(3, 3) = Type(2., 0.);
  apply2QubitGate(qubit, qubit2, xz_plus_2id);

// update sum adding <psi|XZ|psi>
  BaseType final_norm = computenorm();
  sum += coeff*(final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 4, 4, 32> inverse;
  inverse(0, 2) = inverse(0, 3) = inverse(1, 2) = inverse(1, 3) = Type(0., 0.);
  inverse(2, 0) = inverse(2, 1) = inverse(3, 0) = inverse(3, 1) = Type(0., 0.);
  inverse(0, 1) = inverse(1, 0) = Type(-1./3., 0.);
  inverse(2, 3) = inverse(3, 2) = Type( 1./3., 0.);
  inverse(0, 0) = inverse(1, 1) = inverse(2, 2) = inverse(3, 3) = Type( 2./3., 0.);
  apply2QubitGate(qubit, qubit2, inverse);
}


/// @brief Compute expectation value of Pauli Y.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Y_qubit.X_qubit2 |psi\>
//------------------------------------------------------------------------------
// TODO : check how the basis is it is: 00-01-10-11 or 00-10-01-11 !!!!
//        this code uses the standard 00-01-10-11 despite this being opposite to backend storage convention
template <class Type>
void QbitRegister<Type>::expectationValueXY(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
  expectationValueYX(qubit2, qubit, sum, coeff);
}


/// @brief Compute expectation value of Pauli Y.Y for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Y_qubit.Y_qubit2 |psi\>
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueYY(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+YY)|psi>
  openqu::TinyMatrix<Type, 4, 4, 32> yy_plus_2id;
  yy_plus_2id(0, 1) = yy_plus_2id(0, 2) = yy_plus_2id(1, 0) = yy_plus_2id(1, 3) = Type(0., 0.);
  yy_plus_2id(2, 0) = yy_plus_2id(2, 3) = yy_plus_2id(3, 1) = yy_plus_2id(3, 2) = Type(0., 0.);
  yy_plus_2id(0, 3) = yy_plus_2id(3, 0) = Type(-1., 0.);
  yy_plus_2id(1, 2) = yy_plus_2id(2, 1) = Type( 1., 0.);
  yy_plus_2id(0, 0) = yy_plus_2id(1, 1) = yy_plus_2id(2, 2) = yy_plus_2id(3, 3) = Type(2., 0.);
  apply2QubitGate(qubit, qubit2, yy_plus_2id);

// update sum adding <psi|YY|psi>
  BaseType final_norm = computenorm();
  sum += coeff*(final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 4, 4, 32> inverse;
  inverse(0, 1) = inverse(0, 2) = inverse(1, 0) = inverse(1, 3) = Type(0., 0.);
  inverse(2, 0) = inverse(2, 3) = inverse(3, 1) = inverse(3, 2) = Type(0., 0.);
  inverse(0, 3) = inverse(3, 0) = Type( 1./3., 0.);
  inverse(1, 2) = inverse(2, 1) = Type(-1./3., 0.);
  inverse(0, 0) = inverse(1, 1) = inverse(2, 2) = inverse(3, 3) = Type( 2./3., 0.);
  apply2QubitGate(qubit, qubit2, inverse);
}


/// @brief Compute expectation value of Pauli Z.Y for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Z_qubit.Y_qubit2 |psi\>
//------------------------------------------------------------------------------
// TODO : check how the basis is it is: 00-01-10-11 or 00-10-01-11 !!!!
//        this code uses the standard 00-01-10-11 despite this being opposite to backend storage convention
template <class Type>
void QbitRegister<Type>::expectationValueZY(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+YZ)|psi>
  openqu::TinyMatrix<Type, 4, 4, 32> yz_plus_2id;
  yz_plus_2id(0, 2) = yz_plus_2id(0, 3) = yz_plus_2id(1, 2) = yz_plus_2id(1, 3) = Type(0., 0.);
  yz_plus_2id(2, 0) = yz_plus_2id(2, 1) = yz_plus_2id(3, 0) = yz_plus_2id(3, 1) = Type(0., 0.);
  yz_plus_2id(0, 1) = yz_plus_2id(3, 2) = Type( 0.,-1.);
  yz_plus_2id(1, 0) = yz_plus_2id(2, 3) = Type( 0., 1.);
  yz_plus_2id(0, 0) = yz_plus_2id(1, 1) = yz_plus_2id(2, 2) = yz_plus_2id(3, 3) = Type(2., 0.);
  apply2QubitGate(qubit, qubit2, yz_plus_2id);

// update sum adding <psi|YZ|psi>
  BaseType final_norm = computenorm();
  sum += coeff*(final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 4, 4, 32> inverse;
  inverse(0, 2) = inverse(0, 3) = inverse(1, 2) = inverse(1, 3) = Type(0., 0.);
  inverse(2, 0) = inverse(2, 1) = inverse(3, 0) = inverse(3, 1) = Type(0., 0.);
  inverse(0, 1) = inverse(3, 2) = Type(0., 1./3.);
  inverse(1, 0) = inverse(2, 3) = Type(0.,-1./3.);
  inverse(0, 0) = inverse(1, 1) = inverse(2, 2) = inverse(3, 3) = Type( 2./3., 0.);
  apply2QubitGate(qubit, qubit2, inverse);
}


/// @brief Compute expectation value of Pauli X.Z for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| X_qubit.Z_qubit2 |psi\>
//------------------------------------------------------------------------------
// TODO : check how the basis is it is: 00-01-10-11 or 00-10-01-11 !!!!
//        this code uses the standard 00-01-10-11 despite this being opposite to backend storage convention
template <class Type>
void QbitRegister<Type>::expectationValueXZ(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
  expectationValueZX(qubit2, qubit, sum, coeff);
}


/// @brief Compute expectation value of Pauli Y.Z for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Y_qubit.Z_qubit2 |psi\>
//------------------------------------------------------------------------------
// TODO : check how the basis is it is: 00-01-10-11 or 00-10-01-11 !!!!
//        this code uses the standard 00-01-10-11 despite this being opposite to backend storage convention
template <class Type>
void QbitRegister<Type>::expectationValueYZ(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
  expectationValueZY(qubit2, qubit, sum, coeff);
}


/// @brief Compute expectation value of Pauli Z.Z for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --\> sum + coeff * \<psi| Z_qubit.Z_qubit2 |psi\>
//------------------------------------------------------------------------------
template <class Type>
void QbitRegister<Type>::expectationValueZZ(unsigned qubit, unsigned qubit2, BaseType &sum, BaseType coeff)
{
#ifdef __ONLY_NORMALIZED_STATES__
  BaseType initial_norm = 1.;
#else
  BaseType initial_norm = computenorm();
#endif

// given initial qureg state in |psi>
// compute the non-normalized vector corresponding to (2+ZZ)|psi>
  openqu::TinyMatrix<Type, 4, 4, 32> zz_plus_2id;
  zz_plus_2id(0, 1) = zz_plus_2id(0, 2) = zz_plus_2id(0, 3) = Type(0., 0.);
  zz_plus_2id(1, 0) = zz_plus_2id(1, 2) = zz_plus_2id(1, 3) = Type(0., 0.);
  zz_plus_2id(2, 0) = zz_plus_2id(2, 1) = zz_plus_2id(2, 3) = Type(0., 0.);
  zz_plus_2id(3, 0) = zz_plus_2id(3, 1) = zz_plus_2id(3, 2) = Type(0., 0.);
  zz_plus_2id(0, 0) = zz_plus_2id(3, 3) = Type(3., 0.); 
  zz_plus_2id(1, 1) = zz_plus_2id(2, 2) = Type(1., 0.);
  apply2QubitGate(qubit, qubit2, zz_plus_2id);

// update sum adding <psi|ZZ|psi>
  BaseType final_norm = computenorm();
  sum += coeff*(final_norm*final_norm - 5.*initial_norm*initial_norm)/4.;

// undo the computation to recover the initial state |psi>
  openqu::TinyMatrix<Type, 4, 4, 32> inverse;
  inverse(0, 1) = inverse(0, 2) = inverse(0, 3) = Type(0., 0.);
  inverse(1, 0) = inverse(1, 2) = inverse(1, 3) = Type(0., 0.);
  inverse(2, 0) = inverse(2, 1) = inverse(2, 3) = Type(0., 0.);
  inverse(3, 0) = inverse(3, 1) = inverse(3, 2) = Type(0., 0.);
  inverse(0, 0) = inverse(3, 3) = Type( 1./3., 0.);
  inverse(1, 1) = inverse(2, 2) = Type(  1.  , 0.);
  apply2QubitGate(qubit, qubit2, inverse);
}


/// @brief Utility function that computes the Hamming weight of a bitstring given as integer.
//------------------------------------------------------------------------------
// Hamming weight function
std::size_t Hamming_weight(std::size_t x)
{
  std::size_t count=0;
  for (count=0; x; count++)
    x &= x-1;
  return count;
}


/// @brief Compute expectation value of a Pauli string for multiple qubits over the full-register state.
/// @param qubits vector of the involved qubit indices
/// @param observables vector of the involved pauli operators {1,2,3} -> {X,Y,Z}
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// At the end of the function, the variable 'sum' is updated via:\n
///     sum --> sum + coeff * <psi| Pauli_String |psi>\n
/// where the Pauli_String is defined by:\n
///     observables[0]_qubit[0] . observables[1]_qubit[1] ...\n
/// and each Pauli operator is encoded according to:\n
///     {1,2,3} --> {X,Y,Z}
//------------------------------------------------------------------------------
// observable:  1==PauliX , 2==PauliY , 3==PauliZ
template <class Type>
void QbitRegister<Type>::expectationValue(std::vector<unsigned> &qubits, std::vector<unsigned> &observables, BaseType &sum, BaseType coeff)
{
  assert( qubits.size() == observables.size() );

// G is matrix from change of basis Y --> Z, such that Ginv.Z.G = Y 
  openqu::TinyMatrix<Type, 2, 2, 32> G;
  BaseType f = 1. / std::sqrt(2.);
  G(0, 0) = G(1, 0) = Type(f , 0.);
  G(0, 1) = Type(0.,-f);
  G(1, 1) = Type(0., f);
// G^dagger = G^-1
  openqu::TinyMatrix<Type, 2, 2, 32> Ginv;
  Ginv(0, 0) = Ginv(0, 1) = Type(f , 0.);
  Ginv(1, 0) = Type(0., f);
  Ginv(1, 1) = Type(0.,-f);

  for (std::size_t i=0; i<qubits.size(); i++)
  {
    if (observables[i]==1)
      applyHadamard(qubits[i]);
    else if (observables[i]==2)
      apply1QubitGate(qubits[i],G);
    else if (observables[i]==3)
      continue;
    else
      assert(0);	// should never be called
  }  

// compute the expectation value
  MPI_Comm comm = openqu::mpi::Environment::comm();
  std::size_t myrank = openqu::mpi::Environment::rank();
  BaseType local_value = 0;
  std::size_t glb_start = UL(myrank) * localSize();
// integer in binary notation with 1 located at the position of the qubits
  std::size_t y=0;
  for (std::size_t i=0; i<qubits.size(); i++)
    y += 1 << qubits[i];

#pragma omp parallel
{
  std::size_t x;
  #pragma omp for reduction(+ : local_value)
  for(std::size_t i = 0; i < localSize(); i++)
  {
     x = glb_start + i;
     if ( Hamming_weight( x & y ) & 1 )	// odd number of 1 in the qubits involved
       local_value -= std::abs(state[i]) * std::abs(state[i]);
     else
       local_value += std::abs(state[i]) * std::abs(state[i]);
  }
}
  BaseType global_value;
  // MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_value, &global_value,  MPI_SUM, comm);

// update sum
  sum += coeff * global_value;

// undo the quantum gates to recover the initial state |psi>
  for (std::size_t i=0; i<qubits.size(); i++)
  {
    if (observables[i]==1)
      applyHadamard(qubits[i]);
    else if (observables[i]==2)
      apply1QubitGate(qubits[i],Ginv);
    else if (observables[i]==3)
      continue;
    else
      assert(0);	// should never be called
  }  
}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;

/// @}
