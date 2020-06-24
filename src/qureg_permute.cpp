/// @file qureg_permute.cpp
/// @brief Define the @c QubitRegister methods to permute the index of the qubits.

#include "../include/qureg.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::PermuteQubits(std::vector<std::size_t> new_map, std::string style_of_map)
{
  assert(num_qubits == new_map.size());

  unsigned nprocs = qhipster::mpi::Environment::GetStateSize();

  if (nprocs==1)
  // Single-node implementation.
  {
      this->PermuteLocalQubits(new_map, style_of_map);
  }
  else
  // Multi-node implementation.
  {
      Permutation &qubit_permutation_old = *qubit_permutation;
      Permutation qubit_permutation_new(new_map, style_of_map);

      std::size_t M = this->num_qubits - qhipster::ilog2(qhipster::mpi::Environment::GetStateSize());
      std::vector<std::size_t> int_1_imap, int_2_imap;
      qubit_permutation_old.ObtainIntemediateInverseMaps(qubit_permutation_new.map, M, int_1_imap, int_2_imap);

      // The overall permutation is divided in three parts:
      // - permutation between local qubits
      // - permutation between global qubits
      // - pairwise exchange of local-global qubits
      this->PermuteLocalQubits(int_1_imap, "inverse");
      this->PermuteGlobalQubits(int_2_imap, "inverse");
      this->PermuteByLocalGlobalExchangeOfQubitPairs(new_map, style_of_map);
      // This functions already update qubit_permutation.

#ifndef INTELQS_HAS_MPI
      assert(0);
#elif 0
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
          std::size_t to_glbind = qubit_permutation_new.program2data_(qubit_permutation_old.data2program_(i));
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

#if 0
  // do it multinode
  // calculate displacements for other nodes
  std::vector <std::size_t> counts(nprocs, 0), displs(nprocs, 0);
  for(std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t glbind =
          qubit_permutation_old.bin2dec(qubit_permutation_new.program2data(qubit_permutation_old.data2program(i)));
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
          qubit_permutation_old.bin2dec(qubit_permutation_new.program2data(qubit_permutation_old.data2program(i)));
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

  // Their position are exchanged in the emulation of the SWAP.
  qubit_permutation->ExchangeTwoElements(qubit_1, qubit_2);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::PermuteLocalQubits(std::vector<std::size_t> new_map, std::string style_of_map)
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
  // and that only the local qubits are (eventually) updated.
  std::vector<std::size_t> & old_inverse_map = qubit_permutation->imap;
  std::size_t M = this->num_qubits - qhipster::ilog2(qhipster::mpi::Environment::GetStateSize());
  std::vector<bool> local(new_inverse_map.size(), 0);
  for (unsigned pos=0; pos<M; ++pos)
      local[new_inverse_map[pos]] = 1;
  for (unsigned pos=0; pos<M; ++pos)
      assert( local[old_inverse_map[pos]] > 0 );
  for (unsigned pos=M; pos<num_qubits; ++pos)
      assert( old_inverse_map[pos] == new_inverse_map[pos] );
  
  // Initialize the utility vector: state_old = state;
  Permutation &qubit_permutation_old = *qubit_permutation;
  Permutation qubit_permutation_new(new_inverse_map, "inverse");
  std::vector<Type> state_old(LocalSize(), 0);
#pragma omp parallel for
  for (std::size_t i = 0; i < LocalSize(); i++)
      state_old[i] = state[i];
  
#pragma omp parallel for
  for (std::size_t i = 0; i < LocalSize(); i++)
  {
      std::size_t to = qubit_permutation_new.program2data_(qubit_permutation_old.data2program_(i));
      state[to] = state_old[i];
  }

  // Update permutation:
  // permutation_old is a reference to the permutation pointed by the class variable 'qubit_permutation'.
  qubit_permutation_old = qubit_permutation_new;
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::PermuteGlobalQubits(std::vector<std::size_t> new_map, std::string style_of_map)
{
#ifndef INTELQS_HAS_MPI
  assert(0);
#else
  assert(new_map.size() == this->num_qubits);
  Permutation qubit_permutation_new(new_map, style_of_map);
  std::vector<std::size_t> new_direct_map = qubit_permutation_new.map;
  std::vector<std::size_t> new_inverse_map = qubit_permutation_new.imap;

  // Verify that new map mantains the current distinction between local and global qubits
  // and that only the global qubits are (eventually) updated.
  std::vector<std::size_t> old_direct_map = qubit_permutation->map;
  std::vector<std::size_t> old_inverse_map = qubit_permutation->imap;
  std::size_t M = this->num_qubits - qhipster::ilog2(qhipster::mpi::Environment::GetStateSize());
  std::vector<bool> global(new_inverse_map.size(), 0);
  for (unsigned pos=M; pos<num_qubits; ++pos)
      global[new_inverse_map[pos]] = 1;
  for (unsigned pos=M; pos<num_qubits; ++pos)
      assert( global[old_inverse_map[pos]] > 0 );
  for (unsigned pos=0; pos<M; ++pos)
      assert( old_inverse_map[pos] == new_inverse_map[pos] );
  
  // TODO: non-identity reordering of the global qubits may be implemented by reindexing MPI processes.

  // When more than two qubits change position, the content of MPI ranks are not simply
  // exchanged in paris, but one may have longer cycles like A-->B-->C-->D-->A
  // myrank sends to destination and receives from source
  std::size_t myrank = qhipster::mpi::Environment::GetStateRank();
  std::size_t source(0), destination(0);
  std::size_t glb_start = UL(myrank) * LocalSize(); // based on the position
  
  unsigned qubit;
  unsigned position;
  for (unsigned pos=M; pos<num_qubits; ++pos)
  {
      qubit = old_inverse_map[pos];
      position = new_direct_map[qubit];
      if (check_bit(glb_start, pos) == 1)
          destination += 1UL << (position-M);
      //
      qubit = new_inverse_map[pos];
      position = old_direct_map[qubit];
      if (check_bit(glb_start, pos) == 1)
          source += 1UL << (position-M);
  }
  
  MPI_Status status;
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();

  Type *tmp_state = TmpSpace();
  std::size_t lcl_size = LocalSize();
  size_t lcl_chunk = TmpSize();
  for(size_t c = 0; c < lcl_size; c += lcl_chunk)
  {
      // As tag, we use the source of the corresponding communication.
      std::size_t sendtag(myrank), recvtag(source);
      //printf("%d --> myrank = %d --> %d\n", destination, myrank, source);
      qhipster::mpi::MPI_Sendrecv_x(&(state[c])    , lcl_chunk, destination, sendtag,
                                    &(tmp_state[0]), lcl_chunk, source, recvtag,
                                    comm, &status);
      #pragma omp parallel for 
      for (std::size_t i = 0; i < lcl_chunk; i++)
          state[c+i] = tmp_state[i];
  }
  qubit_permutation->SetNewPermutationFromMap(new_map, style_of_map);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::PermuteByLocalGlobalExchangeOfQubitPairs(std::vector<std::size_t> new_map,
                                                                   std::string style_of_map)
{
  // More than one qubit pair may be exchanged.
  Permutation new_qubit_permutation(new_map, style_of_map);
  // Record if qubits has already been updated: 0 if not seen, 1 if updated.
  std::vector<unsigned> exchanged_qubits(num_qubits);
  unsigned num_pairs = 0;
  unsigned old_position, new_position, partner_qubit;
  std::size_t M = this->num_qubits - qhipster::ilog2(qhipster::mpi::Environment::GetStateSize());
  for (unsigned qubit=0; qubit<num_qubits; ++qubit)
  {
      if (exchanged_qubits[qubit]!=0)
          continue;
      old_position = (*qubit_permutation)[qubit];
      new_position = new_qubit_permutation[qubit];
      if ( new_position != old_position )
      {
          // Find its partner and verify that 1) it was not yet exchanged, 2) they form a 2-cycle
          partner_qubit = qubit_permutation->Find(new_position);
          assert(exchanged_qubits[partner_qubit]==0);
          assert(new_qubit_permutation[partner_qubit]==old_position);
          // Verify that one is local and the other global
          if (old_position<M)
              assert(new_position>=M);
          else
              assert(new_position<M);

          // Actual implementation
          ApplySwap(qubit, partner_qubit);   // move/update the data
          qubit_permutation->ExchangeTwoElements(qubit, partner_qubit);
          ++num_pairs;
          exchanged_qubits[qubit]=num_pairs;
          exchanged_qubits[partner_qubit]=num_pairs;
      }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/////////////////////////////////////////////////////////////////////////////////////////
