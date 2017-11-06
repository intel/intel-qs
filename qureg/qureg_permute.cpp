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

#include "qureg.hpp"

template <class Type>
void QbitRegister<Type>::Permute(std::vector<std::size_t> permutation_new_vec)
{

  assert(nqbits == permutation_new_vec.size());

  Permutation &permutation_old = *permutation;
  Permutation permutation_new(permutation_new_vec);

#ifndef OPENQU_HAVE_MPI
  std::vector<Type> state_new(localSize(), 0);

  for (std::size_t i = 0; i < localSize(); i++) {
    std::size_t to =
        permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
    state_new[to] = state[i];
  }
  // state = state_new;
  for (std::size_t i = 0; i < localSize(); i++) state[i] = state_new[i];
#else
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  // Dummy multi-node permutation code
  std::vector<Type> glb_state(globalSize(), 0);
  MPI_Allgather(&(state[0]), localSize(), MPI_DOUBLE_COMPLEX, &(glb_state[0]), localSize(),
                MPI_DOUBLE_COMPLEX, comm);
  for (std::size_t i = 0; i < glb_state.size(); i++) {
    std::size_t glbind =
        permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
    std::size_t rank = glbind / localSize();
    if (rank == myrank) {
      std::size_t lclind = glbind - rank * localSize();
      assert(lclind < localSize());
      state[lclind] = glb_state[i];
    }
  }
#endif
  permutation_old = permutation_new;

#if 0
  // do it multinode
  // calculate displacements for other nodes
  std::vector <std::size_t> counts(nprocs, 0), displs(nprocs, 0);
  for(std::size_t i = 0; i < localSize(); i++)
  {
    std::size_t glbind = permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
    std::size_t rank = glbind / localSize(); 
    assert(rank < nprocs);
    counts[rank]++;
  }
  // compute displacements for each rank
  for(std::size_t i = 1; i < nprocs; i++)
    displs[i] = displs[i-1] + counts[i-1]; 
 
  // fill in outgoing buffer as key value std::pairs
  std::vector<std::pair<std::size_t, Type>> tmp(localSize());
  for(std::size_t i = 0; i < localSize(); i++)
  {
    std::size_t glbind = permutation_old.bin2dec(permutation_new.perm2lin(permutation_old.lin2perm(i)));
    std::size_t rank = glbind / localSize();
    std::size_t lclind = glbind - rank * nprocs;
    tmp[displs[rank]
  }
#endif
}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;

