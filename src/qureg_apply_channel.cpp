#include "../include/qureg.hpp"

/// \addtogroup qureg
///  @{

/// @file qureg_apply_channel.cpp
/// @brief Define the @c QubitRegister methods related to implementing simulations
/// in presence of noise using the chi matrix.

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply 1-qubit channel.
/// @param qubit the index of the involved qubit
/// @param chi the chi-matrix representing the quantum channel

template <class Type>
void QubitRegister<Type>::ApplyChannel(const unsigned qubit, CM4x4<Type> & chi)
{
  assert (rng_ptr_ != nullptr);

  // Select one of the eigenvalues of chi.
  BaseType r;
  rng_ptr_->UniformRandomNumbers(&r, 1, 0, 1, "state");
  unsigned k=0;
  while (r>chi.GetEigenCumulativeProbability(k))
  {
      ++k;
      if (k>=4)
          assert(0 && "Error: p_cum should be normalized to 1.");
  }
  // Matrix corresponding to the eigenvectors (seen as linear combinations of the Pauli matrices).
  std::vector<Type> estate = chi.GetEigenVector(k);
  // Convert it into an explicit 2x2 matrix.
  TM2x2<Type> matrix; 
  matrix(0, 0) = estate[0] + estate[3];
  matrix(0, 1) = estate[1] - Type(0., 1.) * estate[2];
  matrix(1, 0) = estate[1] + Type(0., 1.) * estate[2];
  matrix(1, 1) = estate[0] - estate[3];
  // If the eigen-value was negative, change overall sign of the channels:
  if (std::real(chi.GetEigenValue(k))<0)
      this->overall_sign_of_channels *= -1;
  // Apply the coresponding 1-qubit gate.
  QubitRegister<Type>::Apply1QubitGate(qubit, matrix);

  //FIXME: delete
#if 0
  std::cout << "k = " << k << " , matrix of the 'gate':\n"
            << "|\t" << matrix(0,0) << "\t" << matrix(0,1) << "\t|\n"
            << "|\t" << matrix(1,0) << "\t" << matrix(1,1) << "\t|\n";
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply 2-qubit channel.
/// @param qubit1 the index of the 1st involved qubit
/// @param qubit2 the index of the 2nd involved qubit
/// @param chi the chi-matrix representing the quantum channel

template <class Type>
void QubitRegister<Type>::ApplyChannel(const unsigned qubit1, const unsigned qubit2, CM16x16<Type> & chi)
{
  assert (rng_ptr_ != nullptr);

  // Select one of the eigenvalues of chi.
  BaseType r;
  rng_ptr_->UniformRandomNumbers(&r, 1, 0, 1, "state");
  unsigned k=0;
  while (r>chi.GetEigenCumulativeProbability(k))
      ++k;
  // Matrix corresponding to the eigenvectors (seen as linear combinations of the Pauli matrices).
  std::vector<Type> estate = chi.GetEigenVector(k);
  // Convert it into an explicit 4x4 matrix.
  // TODO: Verify the choice of the order of matrix elements used in Apply2QubitGate().
  //       Quantum info convention or that of IQS?
  // TODO: This is related to the order of Pauli matrixces in the basis for 2-qubit channels.
  // Here we chose the Pauli order:
  //   {id0.id1, id0.X0, id0.Y1, id0.Z1, X0.id1, ..., Z0.Z1}
  // and the Quantum Info convention:
  //   { |0>_0 . |0>_1  ,  |0>_0 . |1>_1  ,  |1>_0 . |0>_1  ,  |1>_0 . |1>_1 }
  TM4x4<Type> matrix;
  matrix(0, 0) = estate[0] + estate[3] + estate[12] + estate[15];
  matrix(1, 1) = estate[0] - estate[3] + estate[12] + estate[15];
  matrix(2, 2) = estate[0] + estate[3] - estate[12] - estate[15];
  matrix(3, 3) = estate[0] - estate[3] + estate[12] - estate[15];

  matrix(0, 1) = estate[1] + estate[13] - Type(0., 1.)*(estate[2] + estate[14]);
  matrix(0, 1) = std::conj( matrix(1, 0) );
  matrix(2, 3) = estate[1] - estate[13] - Type(0., 1.)*(estate[2] - estate[14]);
  matrix(3, 2) = std::conj( matrix(2, 3) );

  matrix(0, 2) = estate[4] + estate[7]  - Type(0., 1.)*(estate[8] + estate[11]);
  matrix(1, 3) = estate[4] - estate[7]  - Type(0., 1.)*(estate[8] - estate[11]);
  matrix(0, 3) = estate[5] - estate[10] - Type(0., 1.)*(estate[9] + estate[6]);
  matrix(1, 2) = estate[5] + estate[10] - Type(0., 1.)*(estate[9] - estate[6]);

  matrix(2, 0) = std::conj( matrix(0, 2) );
  matrix(3, 1) = std::conj( matrix(1, 3) );
  matrix(3, 0) = std::conj( matrix(0, 3) );
  matrix(2, 1) = std::conj( matrix(1, 2) );

  // Apply the coresponding 2-qubit gate.
  QubitRegister<Type>::Apply2QubitGate(qubit1, qubit2, matrix);
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

} // end namespace iqs

/// @}
