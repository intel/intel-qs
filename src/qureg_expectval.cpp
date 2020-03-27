#include "../include/qureg.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_expectval.cpp
///  @brief Define the @c QubitRegister methods for the expectation values of Pauli strings.

/// @brief Compute expectation value of Pauli X for qubit over the full-register state
/// @param qubit index of the involved qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| X_qubit |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueX(unsigned qubit, BaseType coeff)
{
// compute the expectation value <psi|X|psi> = <psi|H.Z.H|psi>
  ApplyHadamard(qubit);
  BaseType expectation = 1. - 2.*GetProbability(qubit);
// recover the initial state |psi>
  ApplyHadamard(qubit);
  return coeff * expectation;
}


/// @brief Compute expectation value of Pauli Y for qubit over the full-register state
/// @param qubit index of the involved qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Y_qubit |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueY(unsigned qubit, BaseType coeff)
{
// G is matrix from change of basis Y --> Z , such that G^dagger.Z.G=Y
  qhipster::TinyMatrix<Type, 2, 2, 32> G;
  BaseType f = 1. / std::sqrt(2.);
  G(0, 0) = G(1, 0) = Type(f, 0.);
  G(0, 1) = Type(0.,-f);
  G(1, 1) = Type(0., f);
// compute the iexpectation value <psi|Y|psi> = <psi|G^dagger.Z.G|psi>
  Apply1QubitGate(qubit,G);
  BaseType expectation = 1. - 2.*GetProbability(qubit);
// recover the initial state |psi> by applying G^dagger
  G(0, 0) = G(0, 1) = Type(f, 0.);
  G(1, 0) = Type(0., f);
  G(1, 1) = Type(0.,-f);
  Apply1QubitGate(qubit,G);

  return coeff * expectation;
}


/// @brief Compute expectation value of Pauli Z for qubit over the full-register state
/// @param qubit index of the involved qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Z_qubit |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueZ(unsigned qubit, BaseType coeff)
{
  BaseType expectation = 1. - 2.*GetProbability(qubit);
  return coeff * expectation;
}


/// @brief Utility function that computes the Hamming weight of a bitstring given as integer.
//------------------------------------------------------------------------------
// Hamming weight function
std::size_t HammingWeight(std::size_t x)
{
  std::size_t count=0;
  for (count=0; x; count++)
    x &= x-1;
  return count;
}


/// @brief Compute expectation value of a Pauli string for multiple qubits over the full-register state.
/// @param qubits vector of the involved qubit indices
/// @param observables vector of the involved pauli operators {1,2,3} -> {X,Y,Z}
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value (with |psi> possibly not normalized):\n
///     coeff * <psi| Pauli_String |psi>\n
/// where the Pauli_String is defined by:\n
///     observables[0]_qubit[0] . observables[1]_qubit[1] ...\n
/// and each Pauli operator is encoded according to:\n
///     {1,2,3} --> {X,Y,Z}
//------------------------------------------------------------------------------
// observable:  1==PauliX , 2==PauliY , 3==PauliZ
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValue(std::vector<unsigned> &qubits,
                                      std::vector<unsigned> &observables,
                                      BaseType coeff)
{
// Checks on the input values.
  assert( qubits.size() == observables.size() );
  for (unsigned j=0; j<qubits.size(); ++j)
  {
      assert(qubits[j]<num_qubits);
      assert(observables[j]>0 && observables[j]<4);
  }

// Special cases defined above for single or double Pauli matrices.
// Recall that the special case with double Pauli matrix uses Apply2QubitGate()
// that is currently available only for the single-MPI-rank case.
  unsigned nprocs = 1;
  nprocs = qhipster::mpi::Environment::GetStateSize();

  if (qubits.size()==0)
  {
      return coeff;
  }
  else if (qubits.size()==1)
  {
      if (observables[0]==1)
          return ExpectationValueX(qubits[0],coeff);
      else if (observables[0]==2)
          return ExpectationValueY(qubits[0],coeff);
      else if (observables[0]==3)
          return ExpectationValueZ(qubits[0],coeff);
  }

// G is matrix from change of basis Y --> Z, such that Ginv.Z.G = Y 
  qhipster::TinyMatrix<Type, 2, 2, 32> G;
  BaseType f = 1. / std::sqrt(2.);
  G(0, 0) = G(1, 0) = Type(f , 0.);
  G(0, 1) = Type(0.,-f);
  G(1, 1) = Type(0., f);
// G^dagger = G^-1
  qhipster::TinyMatrix<Type, 2, 2, 32> Ginv;
  Ginv(0, 0) = Ginv(0, 1) = Type(f , 0.);
  Ginv(1, 0) = Type(0., f);
  Ginv(1, 1) = Type(0.,-f);

  for (std::size_t i=0; i<qubits.size(); i++)
  {
    if (observables[i]==1)
      ApplyHadamard(qubits[i]);
    else if (observables[i]==2)
      Apply1QubitGate(qubits[i],G);
    else if (observables[i]==3)
      continue;
    else
      assert(0);	// should never be called
  }  

// compute the expectation value
  std::size_t myrank = qhipster::mpi::Environment::GetStateRank();
  BaseType local_value = 0;
  std::size_t glb_start = UL(myrank) * LocalSize();
// integer in binary notation with 1 located at the position of the qubits
  std::size_t y=0;
  for (std::size_t i=0; i<qubits.size(); i++)
  {
      y += 1 << qubits[i];
  }

#pragma omp parallel
  {
      std::size_t x;
      #pragma omp for reduction(+ : local_value)
      for(std::size_t i = 0; i < LocalSize(); i++)
      {
         x = glb_start + i;
         if ( HammingWeight( x & y ) & 1 )	// odd number of 1 in the qubits involved
           local_value -= std::norm(state[i]);
         else
           local_value += std::norm(state[i]);
      }
  }

  BaseType global_value ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  // MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, comm);
  qhipster::mpi::MPI_Allreduce_x(&local_value, &global_value, 1, MPI_SUM, comm);
#else
  global_value = local_value;
#endif

// update expectation
  BaseType expectation = global_value;

// undo the quantum gates to recover the initial state |psi>
  for (std::size_t i=0; i<qubits.size(); i++)
  {
    if (observables[i]==1)
      ApplyHadamard(qubits[i]);
    else if (observables[i]==2)
      Apply1QubitGate(qubits[i],Ginv);
    else if (observables[i]==3)
      continue;
    else
      assert(0);	// should never be called
  }  

  return coeff * expectation;
}


/////////////////////////////////////////////////////////////////////////////////////////
// The specialized funcitons are only wraps around the general case.
// Introduced for user's convenience.
/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Compute expectation value of Pauli X.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| X_qubit.X_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueXX(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  std::vector<unsigned> qubits = {qubit, qubit2};
  std::vector<unsigned> observables = {1, 1};
  return this->ExpectationValue(qubits, observables, coeff);
}


/// @brief Compute expectation value of Pauli Y.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Y_qubit.X_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueYX(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  std::vector<unsigned> qubits = {qubit, qubit2};
  std::vector<unsigned> observables = {2, 1};
  return this->ExpectationValue(qubits, observables, coeff);
}


/// @brief Compute expectation value of Pauli Z.X for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Z_qubit.X_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueZX(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  std::vector<unsigned> qubits = {qubit, qubit2};
  std::vector<unsigned> observables = {3, 1};
  return this->ExpectationValue(qubits, observables, coeff);
}


/// @brief Compute expectation value of Pauli X.Y for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| X_qubit.Y_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueXY(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  return ExpectationValueYX(qubit2, qubit, coeff);
}


/// @brief Compute expectation value of Pauli Y.Y for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Y_qubit.Y_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueYY(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  std::vector<unsigned> qubits = {qubit, qubit2};
  std::vector<unsigned> observables = {2, 2};
  return this->ExpectationValue(qubits, observables, coeff);
}


/// @brief Compute expectation value of Pauli Z.Y for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Z_qubit.Y_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueZY(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  std::vector<unsigned> qubits = {qubit, qubit2};
  std::vector<unsigned> observables = {3, 2};
  return this->ExpectationValue(qubits, observables, coeff);
}


/// @brief Compute expectation value of Pauli X.Z for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| X_qubit.Z_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueXZ(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  return ExpectationValueZX(qubit2, qubit, coeff);
}


/// @brief Compute expectation value of Pauli Y.Z for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param coeff scalar coefficient that multiplies the expectation value (optional)
///
/// Return the expectation value:\n
///     coeff * <psi| Y_qubit.Z_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueYZ(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  return ExpectationValueZY(qubit2, qubit, coeff);
}


/// @brief Compute expectation value of Pauli Z.Z for two qubits over the full-register state
/// @param qubit index of the first qubit
/// @param qubit2 index of the second qubit
/// @param sum contains the initial value to which the expectation value will be added
/// @param coeff scalar coefficient that multiplies the expectation value
///
/// Return the expectation value:\n
///     coeff * <psi| Z_qubit.Z_qubit2 |psi>
//------------------------------------------------------------------------------
template <class Type>
typename QubitRegister<Type>::BaseType
QubitRegister<Type>::ExpectationValueZZ(unsigned qubit, unsigned qubit2, BaseType coeff)
{
  std::vector<unsigned> qubits = {qubit, qubit2};
  std::vector<unsigned> observables = {3, 3};
  return this->ExpectationValue(qubits, observables, coeff);
}


template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
