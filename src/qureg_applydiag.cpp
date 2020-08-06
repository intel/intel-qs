/// @file qureg_applydiag.cpp
/// @brief Define the @c QubitRegister methods for the application of two-qubit diagonal gates.

#include "../include/qureg.hpp"
#include "../include/highperfkernels.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// General comment.
// To distinguish between program qubits (used in the algorithm) and data qubits
// (used in the representation of the quantum state), we use the term:
// - 'position' to refer to data qubits
// - 'qubit' ro refer to program qubits
/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::ApplyDiagSimp(unsigned qubit_1, unsigned qubit_2,  TM4x4<Type> const &m)
{
  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
#endif

  Type d00 = m[0][0],
       d11 = m[1][1],
       d22 = m[2][2],
       d33 = m[3][3];

  unsigned position_1 = (*qubit_permutation)[qubit_1];
  unsigned position_2 = (*qubit_permutation)[qubit_2];
  assert(position_1 < num_qubits);
  assert(position_2 < num_qubits);

  std::size_t src_glb_start = UL(myrank) * LocalSize();
  for (std::size_t i = 0;  i < LocalSize(); i++)
  {
    if(check_bit(src_glb_start + i, position_1) == 0 &&
       check_bit(src_glb_start + i, position_2) == 0 )
        state[i] *= d00;
    else if(check_bit(src_glb_start + i, position_1) == 0 &&
            check_bit(src_glb_start + i, position_2) == 1 )
        state[i] *= d11;
    else if(check_bit(src_glb_start + i, position_1) == 1 &&
            check_bit(src_glb_start + i, position_2) == 0 )
        state[i] *= d22;
    else 
        state[i] *= d33;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
void QubitRegister<Type>::ApplyDiag(unsigned qubit_1, unsigned qubit_2,  TM4x4<Type> const &m)
{
  assert(qubit_1 < num_qubits);
  assert(qubit_2 < num_qubits);
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

  unsigned position_1 = (*qubit_permutation)[qubit_1];
  unsigned position_2 = (*qubit_permutation)[qubit_2];
  assert(position_1 < num_qubits);
  assert(position_2 < num_qubits);

  unsigned myrank=0, nprocs=1, log2_nprocs=0;
#ifdef INTELQS_HAS_MPI
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  log2_nprocs = qhipster::ilog2(nprocs);
#endif
  unsigned M = num_qubits - log2_nprocs;

  std::size_t delta1 = 1 << position_1;
  std::size_t delta2 = 1 << position_2;

  Type d00 = m[0][0],
       d11 = m[1][1],
       d22 = m[2][2],
       d33 = m[3][3];
  std::size_t src_glb_start = UL(myrank) * LocalSize();
#if 0
  bool controlled = (d00 == 1. && d11 == 1.);

  // currently disabled. controlled part is done inline inside controlled gate
  if (controlled == true)
  {
      // controlled 2-qubit diagonal gate
      if (position_1 < M && position_2 < M)
      {
          printf("here1\n");
          unsigned delta1_ = std::min(delta1, delta2);
          unsigned delta2_ = std::max(delta1, delta2);
          TM2x2<Type> m2x2;
          m2x2[0][0] = d22;
          m2x2[0][1] = 0.;
          m2x2[1][0] = 0.;
          m2x2[1][1] = d33;
          Loop_TN(state, 
                  0UL, LocalSize(), 2UL * delta2_,
                  0UL, delta2_,     2UL * delta1_,
                  delta1_, 2UL*delta1_, delta2_,
                  m2x2, NULL);
      }
      else if (position_1 >= M && position_2 < M)
      {
          printf("here2\n");
          if (check_bit(src_glb_start, position_1) != 0) {
            for (std::size_t j = 0; j < LocalSize(); j += 2 * delta2) {
              for (std::size_t k = 0; k < delta2; ++k) {
                 state[j + k                 ] *= d22;
                 state[j + k + delta2        ] *= d33;
              }
            }
          } 
      }
      else if (position_1 < M && position_2 >= M)
      {
          Type d0, d1;
          printf("here3\n");
          d1 = (check_bit(src_glb_start, position_2) == 0) ? d22 : d33;
          for (std::size_t j = 0; j < LocalSize(); j += 2 * delta1) {
            for (std::size_t k = 0; k < delta1; ++k) {
               state[j + k + delta1        ] *= d1;
            }
          }
      }
      else if (position_1 >= M && position_2 >= M)
      {
          printf("here4\n");
          if (check_bit(src_glb_start, position_1) == 1 &&
              check_bit(src_glb_start, position_2) == 0 ) {
            for (std::size_t i = 0;  i < LocalSize(); i++)
              state[i] *= d22;
          } else if (check_bit(src_glb_start, position_1) == 1 &&
                     check_bit(src_glb_start, position_2) == 1 ) {
            for (std::size_t i = 0;  i < LocalSize(); i++)
              state[i] *= d33;
          }
      }
  } else 
  #endif
  {
    // general 2-qubit diagonal gate
    if (position_1 < M && position_2 < M)
    {
        unsigned delta1_ = std::min(delta1, delta2);
        unsigned delta2_ = std::max(delta1, delta2);
        for (std::size_t i = 0; i < LocalSize(); i += 2 * delta2_)
        {
            for (std::size_t j = 0; j < delta2_; j += 2 * delta1_)
            {
                for (std::size_t k = 0; k < delta1_; ++k)
                {
                    state[i + j + k                 ]  *= d00;
                    state[i + j + k + delta2        ]  *= d11;
                    state[i + j + k + delta1        ]  *= d22;
                    state[i + j + k + delta1 + delta2] *= d33;
                }
            }
        }
    }
    else if (   position_1 >= M && position_2 <  M
             || position_1 <  M && position_2 >= M)
    {
        Type d0, d1;
        unsigned delta;
        if (position_1 >= M && position_2 < M)  {
            d0 = (check_bit(src_glb_start, position_1) == 0) ? d00 : d22;
            d1 = (check_bit(src_glb_start, position_1) == 0) ? d11 : d33;
            delta = delta2;
        }
        else
        {
            d0 = (check_bit(src_glb_start, position_2) == 0) ? d00 : d11;
            d1 = (check_bit(src_glb_start, position_2) == 0) ? d22 : d33;
            delta = delta1;
        }
        for (std::size_t j = 0; j < LocalSize(); j += 2 * delta)
        {
            for (std::size_t k = 0; k < delta; ++k)
            {
                state[j + k                 ] *= d0;
                state[j + k + delta         ] *= d1;
            }
        } 
    }
    else if (position_1 >= M && position_2 >= M)
    {
        if (check_bit(src_glb_start, position_1) == 0 &&
            check_bit(src_glb_start, position_2) == 0 ) {
            for (std::size_t i = 0;  i < LocalSize(); i++)
                state[i] *= d00;
        }
        else if (check_bit(src_glb_start, position_1) == 0 &&
                 check_bit(src_glb_start, position_2) == 1 )
        {
            for (std::size_t i = 0;  i < LocalSize(); i++)
                state[i] *= d11;
        }
        else if (check_bit(src_glb_start, position_1) == 1 &&
                 check_bit(src_glb_start, position_2) == 0 )
        {
            for (std::size_t i = 0;  i < LocalSize(); i++)
                state[i] *= d22;
        }
        else
        {
            for (std::size_t i = 0;  i < LocalSize(); i++)
                state[i] *= d33;
        }
    }
  }

}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;
