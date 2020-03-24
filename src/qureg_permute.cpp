#include "../include/qureg.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_permute.cpp
/// @brief Define the @c QubitRegister methods to permute the index of the qubits.

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Permute(std::vector<std::size_t> permutation_new_vector)
{
  assert(num_qubits == permutation_new_vector.size());

  Permutation &permutation_old = *permutation;
  Permutation permutation_new(permutation_new_vector);

#ifndef INTELQS_HAS_MPI
  std::vector<Type> state_new(LocalSize(), 0);

  for (std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t to =
          permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
      state_new[to] = state[i];
  }
  // state = state_new;
  for (std::size_t i = 0; i < LocalSize(); i++)
      state[i] = state_new[i];
#else
  unsigned myrank = qhipster::mpi::Environment::GetStateRank();
  unsigned nprocs = qhipster::mpi::Environment::GetStateSize();
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();

  // Dummy multi-node permutation code
  std::vector<Type> glb_state(GlobalSize(), 0);
#ifdef BIGMPI
  MPIX_Allgather_x(&(state[0]), LocalSize(), MPI_DOUBLE_COMPLEX, &(glb_state[0]), LocalSize(),
                MPI_DOUBLE_COMPLEX, comm);
#else
  MPI_Allgather(&(state[0]), LocalSize(), MPI_DOUBLE_COMPLEX, &(glb_state[0]), LocalSize(),
                MPI_DOUBLE_COMPLEX, comm);
#endif //BIGMPI
  for (std::size_t i = 0; i < glb_state.size(); i++)
  {
      std::size_t glbind =
          permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
      std::size_t rank = glbind / LocalSize();
      if (rank == myrank)
      {
          std::size_t lclind = glbind - rank * LocalSize();
          assert(lclind < LocalSize());
          state[lclind] = glb_state[i];
      }
  }
#endif
  permutation_old = permutation_new;

#if 0
  // do it multinode
  // calculate displacements for other nodes
  std::vector <std::size_t> counts(nprocs, 0), displs(nprocs, 0);
  for(std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t glbind =
          permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
    std::size_t rank = glbind / LocalSize(); 
    assert(rank < nprocs);
    counts[rank]++;
  }
  // compute displacements for each rank
  for(std::size_t i = 1; i < nprocs; i++)
      displs[i] = displs[i-1] + counts[i-1]; 
 
  // fill in outgoing buffer as key value std::pairs
  std::vector<std::pair<std::size_t, Type>> tmp(LocalSize());
  for(std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t glbind =
          permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
      std::size_t rank = glbind / LocalSize();
      std::size_t lclind = glbind - rank * nprocs;
      tmp[displs[rank]
  }
#endif
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
