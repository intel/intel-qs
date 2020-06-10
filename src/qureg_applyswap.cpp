#include "../include/qureg.hpp"
#include "../include/highperfkernels.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_applyswap.cpp
/// @brief Define the @c QubitRegister methods for the application of swap-related gates.

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief SWAP gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitely, the gate corresponds to the matrix:\n
///           | 1  0  0  0 |\n
///    SWAP = | 0  0  1  0 |\n
///           | 0  1  0  0 |\n
///           | 0  0  0  1 |\n
template <class Type>
void QubitRegister<Type>::ApplySwap(unsigned qubit1, unsigned qubit2)
{
  qhipster::TinyMatrix<Type, 2, 2, 32> notg;
  notg(0, 0) = notg(1, 1) = {0, 0};
  notg(0, 1) = notg(1, 0) = {1, 0};
#if 0
  TODO(Use this implementation of swap till we fix tmp buffer size issue with code below) 
  TODO(    namely need to be able to use Send and Recv properly)
  TODO(    The same problem with controlled gates was already solved. Import solution.)
  unsigned b1 = qubit1, b2 = qubit2;
  ApplyCPauliX(b1, b2);
  ApplyCPauliX(b2, b1);
  ApplyCPauliX(b1, b2);
#else
  ApplySwap_helper(qubit1, qubit2, notg);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief iSWAP gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitely, the gate corresponds to the matrix:\n
///            | 1  0  0  0 |\n
///    iSWAP = | 0  0  i  0 |\n
///            | 0  i  0  0 |\n
///            | 0  0  0  1 |\n
template <class Type>
void QubitRegister<Type>::ApplyISwap(unsigned qubit1, unsigned qubit2)
{
  qhipster::TinyMatrix<Type, 2, 2, 32> g;
  g(0, 0) = g(1, 1) = {0, 0};
  g(0, 1) = g(1, 0) = {0, 1};
  ApplySwap_helper(qubit1, qubit2, g);
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief sqrt(iSWAP) gate
/// @param qubit1 index of the first qubit
/// @param qubit2 index of the second qubit
///
/// Explicitely, the gate corresponds to the matrix:\n
///                            | 1  0  0  0 |\n
///    sqrt(iSWAP) = 1/sqrt(2) | 0  1  i  0 |\n
///                            | 0  i  1  0 |\n
///                            | 0  0  0  1 |\n
template <class Type>
void QubitRegister<Type>::ApplySqrtISwap(unsigned qubit1, unsigned qubit2)
{
  qhipster::TinyMatrix<Type, 2, 2, 32> g;
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
/// Explicitely, the gate corresponds to the matrix:\n
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

  qhipster::TinyMatrix<Type, 2, 2, 32> g;
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
  TODO(Error with the tmp buffer size. The code  needs to be able to use Send and Recv properly.)
  TODO(The same problem with controlled gates was already solved. Import solution.)

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
  unsigned position_1 = (*permutation)[qubit_1];
  assert(position_1 < num_qubits);
  unsigned position_2 = (*permutation)[qubit_2];
  assert(position_2 < num_qubits);

  // For simplicity, we choose position_1 s.t. position_1 < position_2
  // All SWAP-type gates are symmetric and therefore there is no change to the matrix m.
  if (position_1 > position_2) std::swap(position_1, position_2);

  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  log2_nprocs = qhipster::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;

  // Consider the size of the local part of the state and of temporary buffers.
  std::size_t lcl_size_half = LocalSize() / 2UL;
  std::size_t lcl_size_quarter = lcl_size_half / 2UL;
  // At least 4 amplitudes per MPI process. Corresponding to assert(M<num_qubits-2).
  assert(lcl_size_quarter >= 1);
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

  int tag1 = 1,tag2 = 2;
  int itask, jtask;

  std::size_t delta_1 = 1 << position_1 ;
  std::size_t delta_2 = 1 << position_2;

  //FIXME delete next 4 lines
  unsigned &qubit1 = position_1;
  unsigned &qubit2 = position_2;
  std::size_t &delta1 = delta_1;
  std::size_t &delta2 = delta_2;

  std::string gate_name = "TQG("+qhipster::toString(position_1)+","+qhipster::toString(position_2)+")::"+m.name;

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
#ifdef INTELQS_HAS_MPI
    MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
    MPI_Status status;
    // HP_Distrpair(qubit1, qubit2);
    std::size_t src_glb_start = UL(myrank) * LocalSize();
    if (qubit1 < M)
    {
        // printf("here1\n");
        if (check_bit(src_glb_start, qubit2) == 0)
        {
            std::size_t dst_glb_start = set_bit(src_glb_start, qubit2);
            // printf("lcl_size=%lu dst_glb_start=%lu\n", LocalSize(), dst_glb_start);
            assert((dst_glb_start % LocalSize()) == 0);
            itask = myrank;
            jtask = dst_glb_start / LocalSize();
            assert(jtask > myrank);
            // printf("%2d(%3lu) ==> %2d(%3lu)\n", myrank, src_glb_start, jtask, dst_glb_start);
        } else if (check_bit(src_glb_start, qubit2) == 1) {
            std::size_t dst_glb_start = clear_bit(src_glb_start, qubit2);
            // printf("lcl_size=%lu dst_glb_start=%lu\n", LocalSize(), dst_glb_start);
            assert((dst_glb_start % LocalSize()) == 0);
            jtask = myrank;
            itask = dst_glb_start / LocalSize();
            assert(itask < myrank);
            // printf("%2d(%3lu) ==> %2d(%3lu)\n", myrank, src_glb_start, itask, dst_glb_start);
        }
    }
    else
    {
        // printf("here2\n");
        if (check_bit(src_glb_start, qubit1) == 1 &&
            check_bit(src_glb_start, qubit2) == 0)
        {
            std::size_t dst_glb_start = set_bit(clear_bit(src_glb_start, qubit1), qubit2);
            // printf("lcl_size=%lu dst_glb_start=%lu\n", LocalSize(), dst_glb_start);
            assert((dst_glb_start % LocalSize()) == 0);
            itask = myrank;
            jtask = dst_glb_start / LocalSize();
            // qhipster::mpi::Environment::RemapStateRank(jtask);
            assert(jtask > myrank);
            // printf("%2d(%3lu) ==> %2d(%3lu)\n", myrank, src_glb_start, jtask, dst_glb_start);
        }
        else if (check_bit(src_glb_start, qubit1) == 0 &&
                 check_bit(src_glb_start, qubit2) == 1)
        {
            std::size_t dst_glb_start = clear_bit(set_bit(src_glb_start, qubit1), qubit2);
            // printf("lcl_size=%lu dst_glb_start=%lu\n", LocalSize(), dst_glb_start);
            assert((dst_glb_start % LocalSize()) == 0);
            jtask = myrank;
            itask = dst_glb_start / LocalSize();
            // qhipster::mpi::Environment::RemapStateRank(itask);
            assert(itask < myrank);
            // printf("%2d(%3lu) ==> %2d(%3lu)\n", myrank, src_glb_start, itask, dst_glb_start);
        }
        else
        {
            // early return if no rank permutation
            return true;
            qhipster::mpi::Environment::RemapStateRank(myrank);
        }
    }

    // 1. allocate temp buffer
    Type *tmp_state = TmpSpace();

    Type *state0=NULL, *state1=NULL;
    double t, tnet = 0;
    if(itask == myrank)
    {  // this is itask
       // 2. src sends s1 to dst into dT
       //    dst sends d2 to src into dT
       t = sec();
   
       std::size_t start_ind = (qubit1 + 1 != M) ? 0 : lcl_size_half;
       qhipster::mpi::MPI_Sendrecv_x(&(state[start_ind]), lcl_size_half, jtask, tag1,
                                     &(tmp_state[0])    , lcl_size_half, jtask, tag2,
                                     comm, &status);

       tnet += sec() - t;

       // 3. src and dst compute
       std::size_t indsht0 = lcl_size_half, indsht1 = 0;
       state0 = (qubit1 < M - 1) ? (state + delta1): state;
       state1 = tmp_state;
       for (std::size_t j = 0; j < lcl_size_half; j += 2 * delta1)
       {
         for(std::size_t i = j; i < std::min(lcl_size_half, j + delta1); i++)
         {
             std::size_t i0 = i + indsht0;
             std::size_t i1 = i + indsht1;
#if 0
             std::swap(state0[i0], state1[i1]);
#else
             Type in0 = state0[i0], in1 = state1[i1];
             state0[i0] = m00 * in0 + m01 * in1;
             state1[i1] = m10 * in0 + m11 * in1;
#endif
         }
       }

       t = sec();
       if (qubit1 + 1 != M)
       {
         qhipster::mpi::MPI_Sendrecv_x(&(tmp_state[0]), lcl_size_half, jtask, tag1,
                                       &(state[0])    , lcl_size_half, jtask, tag2,
                                       comm, &status);
       }
       tnet += sec() - t;

    }
    else
    {  // this is jtask
       // 2. src sends s1 to dst into dT
       //    dst sends d2 to src into dT
       t = sec();
       std::size_t start_ind = (qubit1 + 1 == M) ? 0 : lcl_size_half;
       qhipster::mpi::MPI_Sendrecv_x(&(state[start_ind]), lcl_size_half, itask, tag2,
                                     &(tmp_state[0])    , lcl_size_half, itask, tag1,
                                     comm, &status);
       tnet += sec() - t;

       // 3. src and dst compute
       std::size_t indsht0 = 0, indsht1 = 0;
       state0 = (qubit1 < M - 1) ? (tmp_state + delta1) : tmp_state;
       state1 = state;
       for (std::size_t j = 0; j < lcl_size_half; j += 2 * delta1)
       {
         for(std::size_t i = j; i < std::min(lcl_size_half, j + delta1); i++)
         {
             std::size_t i0 = i + indsht0;
             std::size_t i1 = i + indsht1;
#if 0
             std::swap(state0[i0], state1[i1]);
#else
             Type in0 = state0[i0], in1 = state1[i1];
             state0[i0] = m00 * in0 + m01 * in1;
             state1[i1] = m10 * in0 + m11 * in1;
#endif
         }
       }

       t = sec();
       if(qubit1 + 1 != M) {
         qhipster::mpi::MPI_Sendrecv_x(&(tmp_state[0])        , lcl_size_half, itask, tag2,
                                       &(state[lcl_size_half]), lcl_size_half, itask, tag1,
                                       comm, &status);
       }
       tnet += sec() - t;
    }
     
#else
    assert(0);
#endif
  }

  if (timer)
      timer->Stop();

#if 0
    int itask, jtask;

    if(check_bit(glb_start, qubit2) == 0)
    {
       printf("%d ==> %d\n", myrank, myrank + (1 << (qubit2 - M)));
       qhipster::mpi::Environment::RemapStateRank(myrank + (1 << (qubit2 - M)));
    }
    else
    {
       printf("%d ==> %d\n", myrank, myrank - (1 << (qubit2 - M)));
       qhipster::mpi::Environment::RemapStateRank(myrank - (1 << (qubit2 - M)));
    }
    printf("here2\n");
#endif



  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Unnecessary method, used only to debug the ApplySwap.
// FIXME TODO: remove after distributed implementation of ApplySwap has been tested.
template <class Type>
void QubitRegister<Type>::DebugSwap(unsigned b1, unsigned b2)
{
  assert(b1 < num_qubits);
  assert(b2 < num_qubits);

  ApplyCPauliX(b1, b2);
  ApplyCPauliX(b2, b1);
  ApplyCPauliX(b1, b2);
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
