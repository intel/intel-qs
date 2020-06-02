/// @file qureg_permute.cpp
/// @brief Define the @c QubitRegister methods to permute the index of the qubits.

#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::Permute(std::vector<std::size_t> permutation_new_vector)
{
  assert(num_qubits == permutation_new_vector.size());

  Permutation &permutation_old = *permutation;
  Permutation permutation_new(permutation_new_vector);

  unsigned myrank = qhipster::mpi::Environment::GetStateRank();
  unsigned nprocs = qhipster::mpi::Environment::GetStateSize();

  if (nprocs==1)
  // Single-node implementation.
  {
      // state_old = state;
      std::vector<Type> state_old(LocalSize(), 0);
      for (std::size_t i = 0; i < LocalSize(); i++)
          state_old[i] = state[i];
    
      for (std::size_t i = 0; i < LocalSize(); i++)
      {
          std::size_t to = permutation_new.program2data_(permutation_old.data2program_(i));
          // FIXME delete assertion below after testing
          assert(to == permutation_old.bin2dec(permutation_new.program2data(permutation_old.data2program(i))));
          state[to] = state_old[i];
      }
  }
  else
  // Multi-node implementation.
  {
#ifndef INTELQS_HAS_MPI
      assert(0);
#else
      MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
    
      // FIXME: This is the dummy multi-node permutation code.
      //        It builds the full state locally!
      // TODO : Write the actual distributed implementation.

      // Create a global state locally, then initialize it to the current global state.
      std::vector<Type> glb_state(GlobalSize(), 0);
#ifdef BIGMPI
      MPIX_Allgather_x(&(state[0]), LocalSize(), MPI_DOUBLE_COMPLEX, &(glb_state[0]), LocalSize(),
                    MPI_DOUBLE_COMPLEX, comm);
#else
      MPI_Allgather(&(state[0]), LocalSize(), MPI_DOUBLE_COMPLEX, &(glb_state[0]), LocalSize(),
                    MPI_DOUBLE_COMPLEX, comm);
#endif //BIGMPI

      // Update the original state from its record in 'glb_state'.
      std::size_t to_lclind;
      for (std::size_t i = 0; i < glb_state.size(); i++)
      {
          std::size_t to_glbind = permutation_new.program2data_(permutation_old.data2program_(i));
          std::size_t to_rank = to_glbind / LocalSize();
          if (to_rank == myrank)
          {
              to_lclind = to_glbind - to_rank * LocalSize();
              assert(to_lclind < LocalSize());
              state[to_lclind] = glb_state[i];
          }
      }
#endif
  }

  // permutation_old is a reference to the permutation pointed by the class variable 'permutation'.
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
