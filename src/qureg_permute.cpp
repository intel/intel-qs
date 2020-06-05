/// @file qureg_permute.cpp
/// @brief Define the @c QubitRegister methods to permute the index of the qubits.

#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::Permute(std::vector<std::size_t> new_map,
                                  std::string style_of_map)
{
  assert(num_qubits == new_map.size());

  unsigned nprocs = qhipster::mpi::Environment::GetStateSize();
  if (nprocs==1)
  // Single-node implementation.
  {
      this->PermuteLocal(new_map, style_of_map);
  }
  else
  // Multi-node implementation.
  {
      Permutation &permutation_old = *permutation;
      Permutation permutation_new(new_map, style_of_map);

#ifndef INTELQS_HAS_MPI
      assert(0);
#else
      unsigned myrank = qhipster::mpi::Environment::GetStateRank();
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
      // permutation_old is a reference to the permutation pointed by the class variable 'permutation'.
      permutation_old = permutation_new;
  }


#if 0
  // do it multinode
  // calculate displacements for other nodes
  std::vector <std::size_t> counts(nprocs, 0), displs(nprocs, 0);
  for(std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t glbind =
          permutation_old.bin2dec(permutation_new.program2data(permutation_old.data2program(i)));
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
          permutation_old.bin2dec(permutation_new.program2data(permutation_old.data2program(i)));
      std::size_t rank = glbind / LocalSize();
      std::size_t lclind = glbind - rank * nprocs;
      tmp[displs[rank]
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::EmulateSwap(unsigned qubit_1, unsigned qubit_2)
{
  assert(qubit_1 < num_qubits);
  assert(qubit_2 < num_qubits);

  // Current position of program qubits 1,2.
  unsigned position_1 = (*permutation)[qubit_1];
  unsigned position_2 = (*permutation)[qubit_2];
  assert(position_1 < num_qubits);
  assert(position_2 < num_qubits);

  // Their position are exchanged in the emulation of the SWAP.
  permutation->ExchangeTwoElements(position_1, position_2);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::PermuteLocal(std::vector<std::size_t> new_map, std::string style_of_map)
{
  // Determine the inverse map.
  assert(new_map.size() == this->num_qubits);
  std::vector<std::size_t> new_inverse_map = new_map;
  if (style_of_map=="direct")
      for (std::size_t qubit = 0; qubit < new_map.size(); qubit++)
          new_inverse_map[new_map[qubit]] = qubit;
  else if (style_of_map!="inverse")
      assert(0);

  // Verify that new map mantains the current distinction between local and global qubits.
  std::vector<std::size_t> & old_inverse_map = permutation->imap;
  std::size_t M = this->num_qubits - qhipster::ilog2(qhipster::mpi::Environment::GetStateSize());
  std::vector<bool> local(new_inverse_map.size(), 0);
  for (unsigned j=0; j<M; ++j)
      local[new_inverse_map[j]] = 1;
  for (unsigned j=0; j<M; ++j)
      assert( local[old_inverse_map[j]] > 0 );
  
  // Initialize the utility vector: state_old = state;
  Permutation &permutation_old = *permutation;
  Permutation permutation_new(new_inverse_map, "inverse");
  std::vector<Type> state_old(LocalSize(), 0);
#pragma omp parallel for
  for (std::size_t i = 0; i < LocalSize(); i++)
      state_old[i] = state[i];
  
#pragma omp parallel for
  for (std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t to = permutation_new.program2data_(permutation_old.data2program_(i));
      state[to] = state_old[i];
  }

  // Update permutation:
  // permutation_old is a reference to the permutation pointed by the class variable 'permutation'.
  permutation_old = permutation_new;
}

/////////////////////////////////////////////////////////////////////////////////////////

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/////////////////////////////////////////////////////////////////////////////////////////
