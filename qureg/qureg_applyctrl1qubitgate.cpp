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

/// @file qureg_applyctrl1qubitgate.cpp
/// @brief Define the @c QubitRegister methods for the application of controlled one-qubit gates.

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Arbitrary two-qubit gate.
/// @param qubit_high index of the first qubit
/// @param qubit_low index of the second qubit
/// @param m 4x4 matrix corresponding to the quantum gate
template <class Type>
double QubitRegister<Type>::HP_Distrpair(unsigned control_position, unsigned target_position,
                                         TM2x2<Type> const&m)
{
#ifdef INTELQS_HAS_MPI
  MPI_Status status;
  MPI_Comm comm = openqu::mpi::Environment::comm();
  std::size_t myrank = openqu::mpi::Environment::rank();

  assert(target_position < num_qubits);
  assert(control_position < num_qubits);
  std::size_t M = num_qubits - openqu::ilog2(openqu::mpi::Environment::size());
  std::size_t C = UL(control_position), T = UL(target_position);

  //  Steps:     1.         2.           3.              4.
  //          i    j   | i    j   |  i      j     |  i       j
  //          s1   d1  | s1   d1  |  s1     d1&s1 |  s1&d1   d1&s1
  //          s2   d2  | s2   d2  |  s2&d2  d2    |  s2&d2   d2&s2
  //          T    T   | d2   s1  |  d2&s2  s1&d1 |  T       T

  int tag1 = 1, tag2 = 2;
  std::size_t glb_start = UL(myrank) * LocalSize();
  unsigned int itask, jtask;

  if (check_bit(glb_start, T) == 0)
  {
      itask = myrank;
      jtask = itask + (1 << (T - M));
  }
  else
  {
      jtask = myrank;
      itask = jtask - (1 << (T - M));
  }

  // 1. allocate temp buffer
  Type *tmp_state = TmpSpace();
  std::size_t lcl_size_half = LocalSize() / 2L;
  assert(lcl_size_half <= std::numeric_limits<int>::max());


#if 0
  size_t lcl_chunk = 128;
  if (lcl_chunk > lcl_size_half)
      lcl_chunk = lcl_size_half;
  else
      assert((lcl_size_half % lcl_chunk) == 0);
#else
  size_t lcl_chunk = TmpSize();
  if (lcl_chunk != lcl_size_half)
  {
      fprintf(stderr, "Need to fix chunking first\n");
      assert(0);
  }
#endif

  double t, tnet = 0;
  for(size_t c = 0; c < lcl_size_half; c += lcl_chunk)
  {
    if (itask == myrank)  // this is itask
    {
        // 2. src sends s1 to dst into dT
        //    dst sends d2 to src into dT
        t = sec();
        MPI_Sendrecv_x(&(state[c]), lcl_chunk, jtask, tag1, &(tmp_state[0]),
                       lcl_chunk, jtask, tag2, comm, &status);
        tnet += sec() - t;

        // 3. src and dst compute
        if (M - C == 1) {
            Loop_SN(0L, lcl_chunk, state, tmp_state, lcl_size_half, 0L, m, specialize, timer);
        } else {
            Loop_DN((UL(1) << C), lcl_size_half, C, state, tmp_state, lcl_size_half, 0L,
                    m, specialize, timer);
        }

        t = sec();
        MPI_Sendrecv_x(&(tmp_state[0]), lcl_chunk, jtask, tag1, &(state[c]),
                       lcl_chunk, jtask, tag2, comm, &status);
        tnet += sec() - t;
    }
    else  // this is jtask
    {
        // 2. src sends s1 to dst into dT
        //    dst sends d2 to src into dT
        t = sec();
        MPI_Sendrecv_x(&(state[lcl_size_half + c]), lcl_chunk, itask, tag2,
                       &(tmp_state[0]), lcl_chunk, itask, tag1, comm,
                       &status);
        tnet += sec() - t;

        if (M - C == 1)
        {}    // this is intentional special case: nothing happens
        else
            Loop_DN((UL(1) << C), lcl_size_half, C, tmp_state, state, 0L, 0L, m, specialize, timer);

        t = sec();
        MPI_Sendrecv_x(&(tmp_state[0]), lcl_chunk, itask, tag2,
                       &(state[lcl_size_half + c]), lcl_chunk, itask, tag1,
                       comm, &status);
        tnet += sec() - t;
    }
  }

  double netsize = 2.0 * sizeof(Type) * 2.0 * D(lcl_size_half), netbw = netsize / tnet;
  if (timer) timer->record_cm(tnet, netbw);

#else
  assert(0);
#endif

  return 0.0;
}




/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Helper for the application of controlled one-qubit gates.
// Apply gate to the state vector in the range sind-eind
template <class Type>
bool QubitRegister<Type>::ApplyControlled1QubitGate_helper(unsigned control_, unsigned qubit_,
                                                          TM2x2<Type> const&m,
                                                          std::size_t sind, std::size_t eind)
{
  assert(control_ != qubit_);
  assert(control_ < num_qubits);
#if 0
  printf("New permutation: ");
  for(unsigned i = 0; i < permutation->size(); i++) printf("%u ", (*permutation)[i]);
  printf("\n");
#endif
  unsigned control = (*permutation)[control_];
  assert(control < num_qubits);
  assert(qubit_ < num_qubits);
  unsigned qubit = (*permutation)[qubit_];
  assert(qubit < num_qubits);

  std::size_t C = control, T = qubit;

  unsigned myrank=0, nprocs=1, log2_nprocs=0;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  log2_nprocs = openqu::ilog2(openqu::mpi::Environment::size());
#endif
  unsigned M = num_qubits - log2_nprocs;
  bool HasDoneWork = false;

  std::size_t src_glb_start = UL(myrank) * LocalSize();
  // check for special case of diagonal
  bool diagonal = (m[0][1].real() == 0. && m[0][1].imag() == 0. &&
                   m[1][0].real() == 0. && m[1][0].imag() == 0.);
  
  std::string gate_name = "CSQG("+openqu::toString(C)+","+openqu::toString(T)+")::"+m.name;

  if (timer) timer->Start(gate_name, C, T);

  #if 0
  // not currently used, because it messes up fusion optimization, 
  // not yet supported for diagonal gates
  if (m[0][1].real() == 0. && m[0][1].imag() == 0. &&
      m[1][0].real() == 0. && m[1][0].imag() == 0.)
  {

     Type one = Type(1., 0.);
     openqu::TinyMatrix<Type, 4, 4, 32> md;
     md(0, 0) = Type(1., 0.);
     md(1, 1) = Type(1., 0.);
     md(2, 2) = m[0][0];
     md(3, 3) = m[1][1];

     ApplyDiag(control, qubit, md);
     assert(eind - sind == LocalSize());
     assert(fusion == false);

     return true;
  }
  else 
  #endif
  {
      std::size_t cpstrideexp = C - M;

      if(C < M  && T < M)
      {
        if(C > T)
        {
            // special case when we are blocking in LLC
            // case when C stride is bigger than LLC block size
            // in this case, we only update state if we are
            // within part of the vector that has Cth bit set to 
            // one, since this is control gate
            // Otherwise, we skip computation all together
            if((C >= log2llc) && (LocalSize() > (eind - sind)))
            {
                if(check_bit(sind, C) == 1)
                {
                    Loop_DN(sind, eind, T, state, state,
                            0, 1UL<<T, m, specialize, timer);
                    HasDoneWork = true;
                }
            }
            else
            {
                Loop_TN(state, 
                        sind,  eind,        1UL<<C+1UL,
                        1UL<<C, 1UL<<C+1UL, 1UL<<T+1UL,
                        0L,     1UL<<T,     1UL<<T, m, specialize, timer);
                HasDoneWork = true;
            }
          
        }
        else
        {
            Loop_TN(state, 
                    sind,     eind,       1UL<<T+1UL,
                    0L,       1UL<<T,     1UL<<C+1UL,
                    1UL<<C,   1UL<<C+1UL, 1UL<<T, m, specialize, timer);
            HasDoneWork = true;
        }
      }
      else if (C >= M && T < M)
      {
          assert(C > T);
          if(((myrank >> cpstrideexp) % 2) != 0)
          {
              Loop_DN(sind, eind, T, state, state, 0L, (1UL << T), m, specialize, timer);
              HasDoneWork = true;
          }
      }
      else if (C >= M && T >= M)
      {
          if(((myrank >> cpstrideexp) % 2) != 0)
          {
            if (specialize && diagonal)
            {
                if (check_bit(src_glb_start, T) == 0 )
                    ScaleState(sind, eind, state, m[0][0], timer);
                else
                    ScaleState(sind, eind, state, m[1][1], timer);
            }
            else
            {
                HP_Distrpair(T, m);
                // printf("HPD 1\n");
            }
            HasDoneWork = true;
          } else {
              TODO(Way to fix problem with X and Y specializaion)
              // openqu::mpi::Environment::remaprank(myrank);
          }
      }
      else if (C < M && T >= M)
      {
          if (specialize && diagonal)
          {
              TM2x2<Type> md;
              md[0][0] = {1.0, 0};
              md[0][1] = md[1][0] = {0., 0.};
              md[1][1] = (check_bit(src_glb_start, T) == 0) ? m[0][0] : m[1][1];
              TODO(Insert Loop_SN specialization for this case)
              Loop_DN(sind, eind, C, state, state, 0, 1UL<<C, md, specialize, timer); 
          }
          else
          {
              HP_Distrpair(C, T, m);
              // printf("HPD 2\n");
          }
          HasDoneWork = true;
      }
      else
          assert(0);
    
  }
  if(timer) timer->Stop();

  return HasDoneWork;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Arbitrary two-qubit gate.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
/// @param m 2x2 matrix corresponding to the single-qubit gate (implemented if control qubit is in |1\>)
template <class Type>
void QubitRegister<Type>::ApplyControlled1QubitGate(unsigned control, unsigned qubit,
                                                    TM2x2<Type> const&m)
{
  assert(qubit < num_qubits);

  if (fusion == true)
  {
      assert((*permutation)[qubit] < num_qubits);
      if ((*permutation)[qubit] < log2llc)
      {
          std::string name = "cqg";
          fwindow.push_back(std::make_tuple(name, m, control, qubit));
          return;
      }
      else
      {
          ApplyFusedGates();
          goto L;
      }
  }
  L:
  ApplyControlled1QubitGate_helper(control, qubit, m, 0UL, LocalSize());
}



/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Controlled rotation around the X axis by an angle theta.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
/// @param theta rotation angle
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     exp( -i X theta/2 )\n
/// This convention is based on the fact that the generators
/// of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
template <class Type>
void QubitRegister<Type>::ApplyCRotationX(unsigned const control, unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> rx;
  rx(0, 1) = rx(1, 0) = Type(0, -std::sin(theta / 2.));
  rx(0, 0) = rx(1, 1) = Type(std::cos(theta / 2.), 0);
  ApplyControlled1QubitGate(control, qubit, rx);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Controlled rotation around the Y axis by an angle theta.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
/// @param theta rotation angle
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     exp( -i Y theta/2 )\n
/// This convention is based on the fact that the generators
/// of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
template <class Type>
void QubitRegister<Type>::ApplyCRotationY(unsigned const control, unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> ry;
  ry(0, 1) = Type(-std::sin(theta / 2.), 0.);
  ry(1, 0) = Type( std::sin(theta / 2.), 0.);
  ry(0, 0) = ry(1, 1) = Type(std::cos(theta / 2.), 0);
  ApplyControlled1QubitGate(control, qubit, ry);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Controlled rotation around the Z axis by an angle theta.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
/// @param theta rotation angle
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     exp( -i Z theta/2 )\n
/// This convention is based on the fact that the generators
/// of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
template <class Type>
void QubitRegister<Type>::ApplyCRotationZ(unsigned const control, unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> rz;
  rz(0, 0) = Type(std::cos(theta / 2.), -std::sin(theta / 2.));
  rz(1, 1) = Type(std::cos(theta / 2.), std::sin(theta / 2.));
  rz(0, 1) = rz(1, 0) = Type(0., 0.);
  ApplyControlled1QubitGate(control, qubit, rz);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Cotrolled X Pauli operator.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     i * exp( -i X pi/2 ) = X 
template <class Type>
void QubitRegister<Type>::ApplyCPauliX(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> px;
  px(0, 0) = Type(0., 0.);
  px(0, 1) = Type(1., 0.);
  px(1, 0) = Type(1., 0.);
  px(1, 1) = Type(0., 0.);
  ApplyControlled1QubitGate(control, qubit, px);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Cotrolled Y Pauli operator.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     i * exp( -i Y pi/2 ) = Y 
template <class Type>
void QubitRegister<Type>::ApplyCPauliY(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> py;
  py(0, 0) = Type(0., 0.);
  py(0, 1) = Type(0., -1.);
  py(1, 0) = Type(0., 1.);
  py(1, 1) = Type(0., 0.);
  ApplyControlled1QubitGate(control, qubit, py);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Cotrolled X Pauli operator.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     i * exp( -i Z pi/2 ) = Z 
template <class Type>
void QubitRegister<Type>::ApplyCPauliZ(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(-1., 0.);
  ApplyControlled1QubitGate(control, qubit, pz);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Cotrolled square root of the Z Pauli operator.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to:\n
///     sqrt(Z) 
template <class Type>
void QubitRegister<Type>::ApplyCPauliSqrtZ(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(0., 1.);
  ApplyControlled1QubitGate(control, qubit, pz);
}



/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Cotrolled hadamard gate.
/// @param control index of the control qubit
/// @param qubit index of the target qubit
///
/// Explicitely, when control qubit is in |1\>, the gate corresponds to the 2x2 matrix:\n
///     | 1/sqrt(2)   1/sqrt(2) |\n
///     | 1/sqrt(2)  -1/sqrt(2) |
template <class Type>
void QubitRegister<Type>::ApplyCHadamard(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> h;
  BaseType f = 1. / std::sqrt(2.);
  h(0, 0) = h(0, 1) = h(1, 0) = Type(f, 0.);
  h(1, 1) = Type(-f, 0.);
  ApplyControlled1QubitGate(control, qubit, h);
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
