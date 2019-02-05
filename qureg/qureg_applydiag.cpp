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
#include "highperfkernels.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_applydiag.cpp
/// @brief Define the @c QubitRegister methods for the application of two-qubit diagonal gates.

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::ApplyDiagSimp(unsigned qubit1, unsigned qubit2,  TM4x4<Type> const &m)
{
  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
#endif

  Type d00 = m[0][0],
       d11 = m[1][1],
       d22 = m[2][2],
       d33 = m[3][3];

  std::size_t src_glb_start = UL(myrank) * LocalSize();
  for (std::size_t i = 0;  i < LocalSize(); i++)
  {
    if(check_bit(src_glb_start + i, qubit1) == 0 &&
       check_bit(src_glb_start + i, qubit2) == 0 )
        state[i] *= d00;
    else if(check_bit(src_glb_start + i, qubit1) == 0 &&
            check_bit(src_glb_start + i, qubit2) == 1 )
        state[i] *= d11;
    else if(check_bit(src_glb_start + i, qubit1) == 1 &&
            check_bit(src_glb_start + i, qubit2) == 0 )
        state[i] *= d22;
    else 
        state[i] *= d33;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::ApplyDiag(unsigned qubit1_, unsigned qubit2_,  TM4x4<Type> const &m)
{
 // flush fusion buffers
  if (fusion == true)
  {
      ApplyFusedGates();
  }

  assert(qubit1_ < num_qubits);
  unsigned qubit1 = (*permutation)[qubit1_];
  assert(qubit1 < num_qubits);
  assert(qubit2_ < num_qubits);
  unsigned qubit2 = (*permutation)[qubit2_];
  assert(qubit2 < num_qubits);

  unsigned myrank=0, nprocs=1, log2_nprocs=0;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  log2_nprocs = openqu::ilog2(nprocs);
#endif
  unsigned M = num_qubits - log2_nprocs;

  std::size_t delta1 = 1 << qubit1;
  std::size_t delta2 = 1 << qubit2;

  Type d00 = m[0][0],
       d11 = m[1][1],
       d22 = m[2][2],
       d33 = m[3][3];
  std::size_t src_glb_start = UL(myrank) * LocalSize();
  bool controlled = (d00 == 1. && d11 == 1.);

  
  #if 0
  // currently disabled. controlled part is done inline inside controlled gate
  if (controlled == true)
  {
      // controlled 2-qubit diagonal gate
      if (qubit1 < M && qubit2 < M)
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
      else if (qubit1 >= M && qubit2 < M)
      {
          printf("here2\n");
          if (check_bit(src_glb_start, qubit1) != 0) {
            for (std::size_t j = 0; j < LocalSize(); j += 2 * delta2) {
              for (std::size_t k = 0; k < delta2; ++k) {
                 state[j + k                 ] *= d22;
                 state[j + k + delta2        ] *= d33;
              }
            }
          } 
      }
      else if (qubit1 < M && qubit2 >= M)
      {
          Type d0, d1;
          printf("here3\n");
          d1 = (check_bit(src_glb_start, qubit2) == 0) ? d22 : d33;
          for (std::size_t j = 0; j < LocalSize(); j += 2 * delta1) {
            for (std::size_t k = 0; k < delta1; ++k) {
               state[j + k + delta1        ] *= d1;
            }
          }
      }
      else if (qubit1 >= M && qubit2 >= M)
      {
          printf("here4\n");
          if (check_bit(src_glb_start, qubit1) == 1 &&
              check_bit(src_glb_start, qubit2) == 0 ) {
            for (std::size_t i = 0;  i < LocalSize(); i++)
              state[i] *= d22;
          } else if (check_bit(src_glb_start, qubit1) == 1 &&
                     check_bit(src_glb_start, qubit2) == 1 ) {
            for (std::size_t i = 0;  i < LocalSize(); i++)
              state[i] *= d33;
          }
      }
  } else 
  #endif
  {
    // general 2-qubit diagonal gate
    if (qubit1 < M && qubit2 < M)
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
    else if (qubit1 >= M && qubit2 < M || qubit1 < M && qubit2 >= M)
    {
        Type d0, d1;
        unsigned delta;
        if (qubit1 >= M && qubit2 < M)  {
            d0 = (check_bit(src_glb_start, qubit1) == 0) ? d00 : d22;
            d1 = (check_bit(src_glb_start, qubit1) == 0) ? d11 : d33;
            delta = delta2;
        }
        else
        {
            d0 = (check_bit(src_glb_start, qubit2) == 0) ? d00 : d11;
            d1 = (check_bit(src_glb_start, qubit2) == 0) ? d22 : d33;
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
    else if (qubit1 >= M && qubit2 >= M)
    {
        if (check_bit(src_glb_start, qubit1) == 0 &&
            check_bit(src_glb_start, qubit2) == 0 ) {
            for (std::size_t i = 0;  i < LocalSize(); i++)
                state[i] *= d00;
        }
        else if (check_bit(src_glb_start, qubit1) == 0 &&
                 check_bit(src_glb_start, qubit2) == 1 )
        {
            for (std::size_t i = 0;  i < LocalSize(); i++)
                state[i] *= d11;
        }
        else if (check_bit(src_glb_start, qubit1) == 1 &&
                 check_bit(src_glb_start, qubit2) == 0 )
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

/// @}
