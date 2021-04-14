#include "../include/qureg.hpp"

/// \addtogroup qureg
///  @{

/// @file qureg_measure.cpp
/// @brief Define the @c QubitRegister methods related to measurement operations.

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Return 'true' if qubit is separable and in a computational state.
/// @param qubit the index of the involved qubit
/// @param tolerance the tolerance of the total probability of opposite outcome
/// @return 'true'  if bit is separable and in computational state\n
///         'false' otherwise
///
/// Function to check if one qubit is in a (separable) classical state.
template <class Type>
bool QubitRegister<Type>::IsClassicalBit(unsigned qubit, BaseType tolerance) const
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;

  assert(qubit<num_qubits);
  unsigned position = (*qubit_permutation)[qubit];
  assert(position<num_qubits);

  std::size_t delta = 1UL << position;

  if (position < M)
  {
      bool up = false, down = false;
      for (std::size_t i = 0; i < LocalSize(); i += 2 * delta)
      {
          for (std::size_t j = 0; j < delta; ++j)
          {
              up = up || (std::norm(state[i + j]) > tolerance);
              down = down || (std::norm(state[i + j + delta]) > tolerance);
              if (up && down)
                  return false;
          }
      }
  }
  else
  {
      int up = 0, down = 0;
      std::size_t src_glb_start = UL(myrank) * LocalSize();
      if (check_bit(src_glb_start, position) == 0)
      {
        down = 0;
        for (std::size_t j = 0; j < LocalSize(); ++j) 
            up = up || (std::norm(state[j]) > tolerance);
      }
      else
      {
        up = 0;
        for (std::size_t j = 0; j < LocalSize(); ++j) 
            down = down || (std::norm(state[j]) > tolerance);
      }
      // printf("[%3d] up:%d down:%d\n", myrank, up, down);
      int glb_up, glb_down;
#ifdef INTELQS_HAS_MPI
// MPI_COMM_WORLD has been changed to iqs::mpi::Environment::GetStateComm()
      MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
      MPI_Allreduce(&up, &glb_up, 1, MPI_INT, MPI_LOR, comm);
      MPI_Allreduce(&down, &glb_down, 1, MPI_INT, MPI_LOR, comm);
#else
      assert(0);	// it should never be quibt >= M without MPI
#endif
      // printf("[%3d] glb_up:%d glb_down:%d\n", myrank, glb_up, glb_down);
      if (glb_up && glb_down) return false;

      iqs::mpi::StateBarrier();
  }
  // printf("[%d] here\n", myrank);
  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Collapse the wavefunction as if qubit was measured in the computational basis.
/// @param qubit the index of the involved qubit
/// @param value boolean corresponing to: 'false'=|0> and 'true'=|1>
///
/// Depending on the measurement outcome, provided here through the boolean value
/// 'false'=|0> and 'true'=|1>, half of the amplitudes are resetted to zero.
/// Notice that the state normalization is not preserved.
template <class Type>
void QubitRegister<Type>::CollapseQubit(unsigned qubit, bool value)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;

  assert(qubit<num_qubits);
  unsigned position = (*qubit_permutation)[qubit];
  assert(position<num_qubits);

  std::size_t delta = 1UL << position;

  if (position < M)
  { 
      for (std::size_t i = value ? 0 : delta; i < LocalSize(); i += 2 * delta)
          for (std::size_t j = 0; j < delta; ++j) state[i + j] = 0.;
  }
  else
  {
      std::size_t src_glb_start = UL(myrank) * LocalSize();
      if (check_bit(src_glb_start, position) == 0 && value == true)
      {
          for (std::size_t j = 0; j < LocalSize(); ++j)
              state[j] = 0.;
      }
      else if (check_bit(src_glb_start, position) == 1 && value == false)
      {
          for (std::size_t j = 0; j < LocalSize(); ++j)
              state[j] = 0.;
      }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Return the probability of outcome '-1' when Pauli Z is measured on the qubit.
/// @param qubit the index of the involved qubit
///
/// Return the probability corresponding to the qubit being in state |1>.
/// The state is left unchanged and not collapsed.
template <class Type>
typename QubitRegister<Type>::BaseType QubitRegister<Type>::GetProbability(unsigned qubit)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;

  assert(qubit<num_qubits);
  unsigned position = (*qubit_permutation)[qubit];
  assert(position<num_qubits);

  std::size_t delta = 1UL << position;
  BaseType local_P = 0.;
  if (position < M)
  {
      // here, '0' and '1' for qubit state are within the same rank
      for (std::size_t i = delta; i < LocalSize(); i += 2 * delta)
      {
          for (std::size_t j = 0; j < delta; ++j)
              local_P += std::norm(state[i + j]);
      }
  }
  else
  {
      std::size_t src_glb_start = UL(myrank) * LocalSize();
      if (check_bit(src_glb_start, position) == 1)
      {
        for (std::size_t j = 0; j < LocalSize(); ++j)
            local_P += std::norm(state[j]);
      }
  }

  BaseType global_P;
#ifdef INTELQS_HAS_MPI
// MPI_COMM_WORLD has been changed to iqs::mpi::Environment::GetStateComm()
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  MPI_Allreduce(&local_P, &global_P, 1, MPI_DOUBLE, MPI_SUM, comm);
#else
  global_P = local_P;
#endif
  return global_P;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ?? explanation needed ??
template <class Type>
bool QubitRegister<Type>::GetClassicalValue(unsigned qubit, BaseType tolerance) const
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;

  assert(qubit<num_qubits);
  unsigned position = (*qubit_permutation)[qubit];
  assert(position<num_qubits);

  std::size_t delta = 1UL << position;
  int bit_is_zero = 0, bit_is_one = 0; 
  if (position < M)
  {
      for (std::size_t i = 0; i < LocalSize(); i += 2 * delta)
      {
          for (std::size_t j = 0; j < delta; ++j)
          {
              if (std::norm(state[i + j]) > tolerance)
              {
                  bit_is_zero = 1;
                  goto done;
              }
              if (std::norm(state[i + j + delta]) > tolerance)
              {
                  bit_is_one  = 1;
                  goto done;
              }
          }
      }
  }
  else
  {
      std::size_t src_glb_start = UL(myrank) * LocalSize();
      if (check_bit(src_glb_start, position) == 0)
      {
          for (std::size_t j = 0; j < LocalSize(); j++)
          {
              if (std::norm(state[j]) > tolerance)
              {
                  bit_is_zero = 1;
                  break;
              }
           }
      } else {
          for (std::size_t j = 0; j < LocalSize(); j++)
          {
              if (std::norm(state[j]) > tolerance)
              {
                  bit_is_one = 1;
                  break;
              }
          }
      }
  }

  done:  
  int glb_bit_is_zero, glb_bit_is_one;
#ifdef INTELQS_HAS_MPI
// MPI_COMM_WORLD has been changed to iqs::mpi::Environment::GetStateComm()
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  MPI_Allreduce(&bit_is_zero, &glb_bit_is_zero, 1, MPI_INT, MPI_LOR, comm);
  MPI_Allreduce(&bit_is_one , &glb_bit_is_one , 1, MPI_INT, MPI_LOR, comm);
#else
  glb_bit_is_zero = bit_is_zero;
  glb_bit_is_one  = bit_is_one ;
#endif

  if (glb_bit_is_zero == 1 && glb_bit_is_one == 0)
      return false;
  else if (glb_bit_is_zero == 0 && glb_bit_is_one == 1)
      return true;
  else
      assert(false);  // this should never be called

  return false;   // dummy return
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

} // end namespace iqs

/// @}
