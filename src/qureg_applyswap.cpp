#include "../include/qureg.hpp"
#include "../include/highperfkernels.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_applyswap.cpp
/// @brief Define the @c QubitRegister methods for the application of swap-related gates.

/////////////////////////////////////////////////////////////////////////////////////////

namespace iqs {

/// @brief SWAP gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitly, the gate corresponds to the matrix:\n
///           | 1  0  0  0 |\n
///    SWAP = | 0  0  1  0 |\n
///           | 0  1  0  0 |\n
///           | 0  0  0  1 |\n
template <class Type>
void QubitRegister<Type>::ApplySwap(unsigned qubit1, unsigned qubit2)
{
  iqs::TinyMatrix<Type, 2, 2, 32> notg;
  notg(0, 0) = notg(1, 1) = {0, 0};
  notg(0, 1) = notg(1, 0) = {1, 0};
  ApplySwap_helper(qubit1, qubit2, notg);
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief iSWAP gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitly, the gate corresponds to the matrix:\n
///            | 1  0  0  0 |\n
///    iSWAP = | 0  0  i  0 |\n
///            | 0  i  0  0 |\n
///            | 0  0  0  1 |\n
template <class Type>
void QubitRegister<Type>::ApplyISwap(unsigned qubit1, unsigned qubit2)
{
  iqs::TinyMatrix<Type, 2, 2, 32> g;
  g(0, 0) = g(1, 1) = {0, 0};
  g(0, 1) = g(1, 0) = {0, 1};
  ApplySwap_helper(qubit1, qubit2, g);
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief sqrt(iSWAP) gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitly, the gate corresponds to the matrix:\n
///                            | 1  0  0  0 |\n
///    sqrt(iSWAP) = 1/sqrt(2) | 0  1  i  0 |\n
///                            | 0  i  1  0 |\n
///                            | 0  0  0  1 |\n
template <class Type>
void QubitRegister<Type>::ApplySqrtISwap(unsigned qubit1, unsigned qubit2)
{
  iqs::TinyMatrix<Type, 2, 2, 32> g;
  BaseType f = 1. / std::sqrt(2.);
  g(0, 0) = g(1, 1) = Type(f, 0);
  g(0, 1) = g(1, 0) = Type(0, f);
  ApplySwap_helper(qubit1, qubit2, g);
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Gate of the 'iSWAP' form.
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitly, the gate corresponds to the matrix:\n
///               | 1   0   0  0 |\n
///    iSWAP(m) = | 0  m00 m01 0 |\n
///               | 0  m10 m11 0 |\n
///               | 0   0   0  1 |\n
/// with m01=m10.
template <class Type>
void QubitRegister<Type>::ApplyISwapRotation(unsigned qubit1, unsigned qubit2, TM2x2<Type> const& m)
{
  assert(m(0,1)==m(1,0));
  ApplySwap_helper(qubit1, qubit2, m);
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief sqrt(sqrt(iSWAP)) gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
template <class Type>
void QubitRegister<Type>::Apply4thRootISwap( unsigned qubit1, unsigned qubit2)
{
  auto a = std::polar(.5, M_PI / 8.);
  auto b = std::polar(.5, 7. * M_PI / 8.);

  Type f0(a - b);
  Type f1(a + b);

  iqs::TinyMatrix<Type, 2, 2, 32> g;
  g(0, 0) = f0;
  g(0, 1) = f1;
  g(1, 0) = f1;
  g(1, 1) = f0;
  ApplySwap_helper(qubit1, qubit2, g);
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
bool QubitRegister<Type>::ApplySwap_helper(unsigned qubit_1, unsigned qubit_2, TM2x2<Type> const&m)
{
  // Update counter of the statistics.
  if (gate_counter != nullptr)
  {
      // IQS count the gates acting on specific program qubits.
      gate_counter->TwoQubitIncrement(qubit_1, qubit_2);
  }

  // flush fusion buffers
  if (fusion == true)
  {
      ApplyFusedGates();
  }

  assert(qubit_1 < num_qubits);
  assert(qubit_2 < num_qubits);
  unsigned position_1 = (*qubit_permutation)[qubit_1];
  assert(position_1 < num_qubits);
  unsigned position_2 = (*qubit_permutation)[qubit_2];
  assert(position_2 < num_qubits);

  // For simplicity, we choose position_1 s.t. position_1 < position_2
  // All SWAP-type gates are symmetric and therefore there is no change to the matrix m.
  // The qubits are not used anymore.
  if (position_1 > position_2)
  {
      std::swap(position_1, position_2);
      assert(m(0,1)==m(1,0));
  }

  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;

  // Consider the size of the local part of the state and of temporary buffers.
  std::size_t lcl_size_half = LocalSize() / 2UL;
  std::size_t lcl_size_quarter = lcl_size_half / 2UL;
  // At least 2 amplitudes per MPI process. Corresponding to assert(M<num_qubits-1).
  assert(lcl_size_half >= 1);
  // The buffer size is given by TmpSize(). It is unnecessary to use more then half the LocalSize().
  size_t lcl_chunk = TmpSize();
  if (lcl_chunk > lcl_size_half) 
      lcl_chunk = lcl_size_half;
  else
      assert((lcl_size_half % lcl_chunk) == 0);

  Type m00 = m(0, 0),
       m01 = m(0, 1),
       m10 = m(1, 0),
       m11 = m(1, 1);

  std::size_t delta_1 = 1 << position_1 ;
  std::size_t delta_2 = 1 << position_2;

  std::string gate_name = "TQG("+iqs::toString(position_1)+","+iqs::toString(position_2)+")::"+m.name;

  if (timer)
      timer->Start(gate_name, position_1, position_2);

  if (position_1 < M && position_2 < M)
  {
#if 0
    // Original summation
    for (std::size_t i = 0; i < LocalSize(); i += 2 * delta_2)
        for (std::size_t j = 0; j < delta_2; j += 2 * delta_1)
            for (std::size_t k = 0; k < delta_1; ++k)
            {
                std::size_t i0 = i + j + k + delta_1;
                std::size_t i1 = i + j + k + delta_2;
                Type in0 = state[i0], in1 = state[i1];
                state[i0] = m00 * in0 + m01 * in1;
                state[i1] = m10 * in0 + m11 * in1;
            }
#elif 0
    // Modified summation to be substituted by Loop_TN
    for (std::size_t i = 0; i < LocalSize(); i += 2 * delta_2)
        for (std::size_t j = i+0; j < i+delta_2; j += 2 * delta_1)
            for (std::size_t i0 = j+delta_1; i0 < j+2*delta_1; ++i0)
            {
                std::size_t i1 = i0 + (delta_2-delta_1);
                Type in0 = state[i0], in1 = state[i1];
                state[i0] = m00 * in0 + m01 * in1;
                state[i1] = m10 * in0 + m11 * in1;
            }
#else
    Loop_TN(state, 0UL, LocalSize(), 2*delta_2,
                   0UL, delta_2, 2*delta_1,
                   delta_1, 2*delta_1,
            delta_2-delta_1,  m, specialize, timer);
#endif
  }
  else
  {
        HP_DistrSwap(position_1, position_2, m);
  }

  if (timer)
      timer->Stop();

#if 0
    int itask, jtask;

    if(check_bit(glb_start, qubit2) == 0)
    {
       printf("%d ==> %d\n", myrank, myrank + (1 << (qubit2 - M)));
       iqs::mpi::Environment::RemapStateRank(myrank + (1 << (qubit2 - M)));
    }
    else
    {
       printf("%d ==> %d\n", myrank, myrank - (1 << (qubit2 - M)));
       iqs::mpi::Environment::RemapStateRank(myrank - (1 << (qubit2 - M)));
    }
    printf("here2\n");
#endif

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Arbitrary controlled two-qubit gate.
/// @param control_position position of the control qubit in the current permutation
/// @param target_position position of the target qubit in the current permutation
/// @param m 2x2 matrix corresponding to the quantum gate
template <class Type>
double QubitRegister<Type>::HP_DistrSwap(unsigned low_position, unsigned high_position,
                                         TM2x2<Type> const&m)
{
  assert(LocalSize() > 1);
#ifndef INTELQS_HAS_MPI
  assert(0);
#else
  MPI_Status status;
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  std::size_t myrank = iqs::mpi::Environment::GetStateRank();

  std::size_t M = num_qubits - iqs::ilog2(iqs::mpi::Environment::GetStateSize());
  std::size_t L = UL(low_position), H = UL(high_position);
  // Used when L < H and H >= M, while L may be < or >= M.

  unsigned itask, jtask;
  int tag1 = 1, tag2 = 2, tag3 = 3, tag4 = 4;
  std::size_t src_glb_start = UL(myrank) * LocalSize();

  // 1. allocate temp buffer
  Type *tmp_state = TmpSpace();
  std::size_t lcl_size_half = LocalSize() / 2L;
  assert(lcl_size_half <= std::numeric_limits<int>::max());

  size_t lcl_chunk = TmpSize();
  if (lcl_chunk > lcl_size_half) 
      lcl_chunk = lcl_size_half;
  else
      assert((lcl_size_half % lcl_chunk) == 0);

  double t, tnet = 0;
  Type *state0=nullptr, *state1=nullptr;

  // Case L < M <= H
  //  Steps:     1.         2.           3.              4.
  //          i    j   | i    j   |  i      j     |  i       j
  //          s1   d1  | s1   d1  |  s1     d1&s1 |  s1&d1   d1&s1
  //          s2   d2  | s2   d2  |  s2&d2  d2    |  s2&d2   d2&s2
  //          T    T   | d2   s1  |  d2&s2  s1&d1 |  T       T
  //
  // Special case when L = M-1 and H >= M
  //  Steps:     1.         2.           3.              4.
  //          i    j   | i    j   |  i      j     |  i       j
  //          s1   d1  | s1   d1  |  s1     d1&s2 |  s1      d1&s1
  //          s2   d2  | s2   d2  |  s2&d1  d2    |  s2&d2   d2
  //          T    T   | d1   s2  |  d1&s2  s2&d1 |  T       T
  // 2. i sends s2 (instead of s1) and j sends d1 (instead of d2)
  // 4. the temporary storage does not need to be communicated, it can be erased

  if (L < M) // H >= M
  {
      assert(H >= M);
      if (check_bit(src_glb_start, H) == 0)
      {
          itask = myrank;
          jtask = itask + (1UL << (H-M));
      }
      else
      {
          jtask = myrank;
          itask = jtask - (1UL << (H-M));
      }

      // Which chunk do we have to Send/Recv?
      // If chunk == lcl_size_half, then all chunks are needed, and if chunk=2^L it means
      // that M=L+1 and there is a special case hardcoded.
      // If chunk >= 2^(L+1), then all chunks are needed since every chunk include
      // global indices with L-bit=0,1
      // If chunk < 2^(L+1), then chunks only include global indices with either L-bit=0
      // or L-bit=1
    
      if ( lcl_chunk == lcl_size_half  ||  lcl_chunk >= (1UL<<(L+1)) )
      {
          for(size_t c = 0; c < lcl_size_half; c += lcl_chunk)
          {
              if (itask == myrank)  // this is itask
              {
                  // 2. src sends s1 to dst into dT
                  //    dst sends d2 to src into dT
                  t = sec();
          
                  // When L+1=M we can avoid a cycle of communication.
                  // In this special case, irank exchanges the second half (L=1, H=0) and jrank the first half (L=0, H=1).
                  std::size_t start_ind = (L != M-1) ? c : lcl_size_half+c;
                  iqs::mpi::MPI_Sendrecv_x(&(state[start_ind]), lcl_chunk, jtask, tag1,
                                                &(tmp_state[0])    , lcl_chunk, jtask, tag2,
                                                comm, &status);
                  tnet += sec() - t;
        
                  // 3. src and dst compute
                  if (L == M-1) {
                      Loop_SN(0UL, lcl_chunk, &(state[start_ind]), tmp_state, 0UL, 0UL,
                              m, specialize, timer);
                  } else {
                      Loop_DN(0UL, lcl_chunk, L, &(state[c]), tmp_state, lcl_size_half+(1UL<<L), 0UL,
                              m, specialize, timer);
                  }
        
                  // 4. src sends sT to dst into d2
                  //    dst sends dT to src into s1
                  t = sec();
                  if (L != M-1)
                  {
                      iqs::mpi::MPI_Sendrecv_x(&(tmp_state[0])    , lcl_chunk, jtask, tag3,
                                                    &(state[start_ind]), lcl_chunk, jtask, tag4,
                                                    comm, &status);
                  }
                  tnet += sec() - t;
              }
              else  // this is jtask
              {
                  // 2. src sends s1 to dst into dT
                  //    dst sends d2 to src into sT
                  t = sec();
        
                  // When L+1=M we can avoid a cycle of communication.
                  // In this special case, irank exchanges the second half (L=1, H=0) and jrank the first half (L=0, H=1).
                  std::size_t start_ind = (L != M-1) ? lcl_size_half+c : c;
                  iqs::mpi::MPI_Sendrecv_x(&(state[start_ind]), lcl_size_half, itask, tag2,
                                                &(tmp_state[0])    , lcl_size_half, itask, tag1,
                                                comm, &status);
                  tnet += sec() - t;
        
                  // 3. src and dst compute
                  if (L == M-1) {
                      Loop_SN(0UL, lcl_chunk, &(state[start_ind]), tmp_state, 0UL, 0UL,
                              m, specialize, timer);
                  } else {
                      Loop_DN(0UL, lcl_chunk, L, tmp_state, &(state[c]), (1UL<<L), 0UL,
                              m, specialize, timer);
                  }
        
                  // 4. src sends sT to dst into d2
                  //    dst sends dT to src into s1
                  t = sec();
                  if (L != M-1)
                  {
                      iqs::mpi::MPI_Sendrecv_x(&(tmp_state[0])    , lcl_size_half, itask, tag4,
                                                    &(state[start_ind]), lcl_size_half, itask, tag3,
                                                    comm, &status);
                  }
                  tnet += sec() - t;
              }
          }
      }
      else
      {
          // case: lcl_chunk < lcl_size_half && lcl_chunk <= 2^L
          // should not happen when temporary size is lcl_size_half 
          assert(0);
          // TODO: need to be implemented
      }
  }
  // Case  M <= L < H
  //  Steps:     1.         2.           3.              4.
  //          i    j   | i    j   |  i      j     |  i       j
  //          s1   d1  | s1   d1  |  s1     d1&s1 |  s1&d1   d1&s1
  //          s2   d2  | s2   d2  |  s2&d2  d2    |  s2&d2   d2&s2
  //          T    T   | d2   s1  |  d2&s2  s1&d1 |  T       T
  else // L, H >= M
  {
      assert (L>=M && H>=M);
      if (check_bit(src_glb_start, L) == 1 &&
          check_bit(src_glb_start, H) == 0)
      {
          itask = myrank;
          jtask = itask - (1UL << (L-M)) + (1UL << (H-M));
      }
      else if (check_bit(src_glb_start, L) == 0 &&
               check_bit(src_glb_start, H) == 1)
      {
          jtask = myrank;
          itask = jtask + (1UL << (L-M)) - (1UL << (H-M));
      }
      else
      {
          // early return if no rank permutation
          return true;
      }

      // Which chunk do we have to Send/Recv?
      // It is a trivial operation in which pairs of amplitudes with the same local index
      // are updated in according to m.
    
        for(size_t c = 0; c < lcl_size_half; c += lcl_chunk)
        {
          if(itask == myrank)  // this is itask
          {
              // 2. src sends s1 to dst into dT
              //    dst sends d2 to src into sT
              t = sec();
              iqs::mpi::MPI_Sendrecv_x(&(state[c])    , lcl_chunk, jtask, tag1,
                                            &(tmp_state[0]), lcl_chunk, jtask, tag2,
                                            comm, &status);
              tnet += sec() - t;
    
              // 3. src and dst compute
              Loop_SN(0UL, lcl_chunk, &(state[c]), tmp_state, lcl_size_half, 0UL, m, specialize, timer);
    

              // 4. src sends sT to dst into d2
              //    dst sends dT to src into s1
              t = sec();
              iqs::mpi::MPI_Sendrecv_x(&(tmp_state[0]), lcl_chunk, jtask, tag3,
                                            &(state[c])    , lcl_chunk, jtask, tag4,
                                            comm, &status);
              tnet += sec() - t;
          }
          else  // this is jtask
          {
              // 2. src sends s1 to dst into dT
              //    dst sends d2 to src into sT
              t = sec();
              iqs::mpi::MPI_Sendrecv_x(&(state[c+lcl_size_half]), lcl_size_half, itask, tag2,
                                            &(tmp_state[0])          , lcl_size_half, itask, tag1,
                                            comm, &status);
              tnet += sec() - t;
    
              // 3. src and dst compute
               Loop_SN(0UL, lcl_chunk, &(state[c]), tmp_state, 0UL, 0UL, m, specialize, timer);

              // 4. src sends sT to dst into d2
              //    dst sends dT to src into s1
              t = sec();
              iqs::mpi::MPI_Sendrecv_x(&(tmp_state[0])          , lcl_size_half, itask, tag4,
                                            &(state[lcl_size_half+c]), lcl_size_half, itask, tag3,
                                            comm, &status);
              tnet += sec() - t;
          }
        }
  }

  double netsize = 2.0 * sizeof(Type) * 2.0 * double(lcl_size_half), netbw = netsize / tnet;
  if (timer) timer->record_cm(tnet, netbw);
#endif

  return 0.0;
}
/////////////////////////////////////////////////////////////////////////////////////////

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

} // end namespace iqs

/// @}
