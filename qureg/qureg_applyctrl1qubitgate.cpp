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



template <class Type>
double QbitRegister<Type>::HP_Distrpair(unsigned cpos, unsigned tpos, TM2x2<Type> const&m)
{
#ifdef OPENQU_HAVE_MPI
  MPI_Status status;
  MPI_Comm comm = openqu::mpi::Environment::comm();
  std::size_t myrank = openqu::mpi::Environment::rank();

  assert(tpos < nqbits);
  assert(cpos < nqbits);
  std::size_t M = nqbits - openqu::ilog2(openqu::mpi::Environment::size()), C = UL(cpos), T = UL(tpos);

  /*
    Steps:     1.         2.           3.              4.
            i    j   | i    j   |  i      j     |  i       j
            s1   d1  | s1   d1  |  s1     d1&s1 |  s1&d1   d1&s1
            s2   d2  | s2   d2  |  s2&d2  d2    |  s2&d2   d2&s2
            T    T   | d2   s1  |  d2&s2  s1&d1 |  T       T

  */

  int tag1 = 1, tag2 = 2;
  std::size_t glb_start = UL(myrank) * localSize();
  unsigned int itask, jtask;

  if (check_bit(glb_start, T) == 0) {
    itask = myrank;
    jtask = itask + (1 << (T - M));
  } else {
    jtask = myrank;
    itask = jtask - (1 << (T - M));
  }

  // 1. allocate temp buffer
  Type *tmpstate = tmpspace();
  std::size_t lcl_size_half = localSize() / 2L;
  assert(lcl_size_half <= std::numeric_limits<int>::max());


#if 0
  size_t lcl_chunk = 128;
  if (lcl_chunk > lcl_size_half)
    lcl_chunk = lcl_size_half;
  else
    assert((lcl_size_half % lcl_chunk) == 0);
#else
  size_t lcl_chunk = tmpSize();
  if (lcl_chunk != lcl_size_half) {
    fprintf(stderr, "Need to fix chunking first\n");
    assert(0);
  }
#endif

  double t, tnet = 0;
  for(size_t c = 0; c < lcl_size_half; c += lcl_chunk) {
    if (itask == myrank) {  // this is itask
      // 2. src sends s1 to dst into dT
      //    dst sends d2 to src into dT
      t = sec();
      MPI_Sendrecv_x(&(state[c]), lcl_chunk, jtask, tag1, &(tmpstate[0]),
                     lcl_chunk, jtask, tag2, comm, &status);
      tnet += sec() - t;

      // 3. src and dst compute
      if (M - C == 1) {
        Loop_SN(0L, lcl_chunk, state, tmpstate, lcl_size_half, 0L, m, specialize, timer);
      } else {
        Loop_DN((UL(1) << C), lcl_size_half, C, state, tmpstate, lcl_size_half, 0L,
                m, specialize, timer);
      }

      t = sec();
      MPI_Sendrecv_x(&(tmpstate[0]), lcl_chunk, jtask, tag1, &(state[c]),
                     lcl_chunk, jtask, tag2, comm, &status);
      tnet += sec() - t;

    } else {  // this is jtask
      // 2. src sends s1 to dst into dT
      //    dst sends d2 to src into dT
      t = sec();
      MPI_Sendrecv_x(&(state[lcl_size_half + c]), lcl_chunk, itask, tag2,
                     &(tmpstate[0]), lcl_chunk, itask, tag1, comm,
                     &status);
      tnet += sec() - t;

      if (M - C == 1) {
        // this is intentional special case: nothing happens
      } else {
        Loop_DN((UL(1) << C), lcl_size_half, C, tmpstate, state, 0L, 0L, m, specialize, timer);
      }

      t = sec();
      MPI_Sendrecv_x(&(tmpstate[0]), lcl_chunk, itask, tag2,
                     &(state[lcl_size_half + c]), lcl_chunk, itask, tag1,
                     comm, &status);
      tnet += sec() - t;
    }
  }

  double netsize = 2.0 * sizeof(Type) * 2.0 * D(lcl_size_half), netbw = netsize / tnet;
  if (timer) {timer->record_cm(tnet, netbw); }

#else
  assert(0);
#endif

  return 0.0;
}



// Apply gate to the state vector in the range sind-eind
template <class Type>
bool QbitRegister<Type>::applyControlled1QubitGate_helper(unsigned control_, unsigned qubit_, TM2x2<Type> const&m,
                                                    std::size_t sind, std::size_t eind)
{
  assert(control_ != qubit_);
  assert(control_ < nqbits);
#if 0
  printf("New permutation: ");
  for(unsigned i = 0; i < permutation->size(); i++) printf("%u ", (*permutation)[i]);
  printf("\n");
#endif

  unsigned control = (*permutation)[control_];
  assert(control < nqbits);
  assert(qubit_ < nqbits);
  unsigned qubit = (*permutation)[qubit_];
  assert(qubit < nqbits);
  std::size_t C = control, T = qubit;

  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  unsigned log2_nprocs = openqu::ilog2(openqu::mpi::Environment::size());
  unsigned M = nqbits - log2_nprocs;
  bool HasDoneWork = false;

  std::size_t src_glb_start = UL(myrank) * localSize();
  // check for special case of diagonal
  bool diagonal = (m[0][1].real() == 0. && m[0][1].imag() == 0. &&
                   m[1][0].real() == 0. && m[1][0].imag() == 0.);
  
  std::string gname = "CSQG("+openqu::toString(C)+","+openqu::toString(T)+")::"+m.name;
  if(timer) timer->Start(gname, C, T);

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

     applyDiag(control, qubit, md);
     assert(eind - sind == localSize());
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
        if((C >= log2llc) && (localSize() > (eind - sind)))
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
         if (specialize && diagonal) {
           if (check_bit(src_glb_start, T) == 0 ) {
             ScaleState(sind, eind, state, m[0][0], timer);
           } else {
             ScaleState(sind, eind, state, m[1][1], timer);
           }
         } else {
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
       if (specialize && diagonal) {
         TM2x2<Type> md;
         md[0][0] = {1.0, 0};
         md[0][1] = md[1][0] = {0., 0.};
         md[1][1] = (check_bit(src_glb_start, T) == 0) ? m[0][0] : m[1][1];
         TODO(Insert Loop_SN specialization for this case)
         Loop_DN(sind, eind, C, state, state, 0, 1UL<<C, md, specialize, timer); 
       } else {
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

template <class Type>
void QbitRegister<Type>::applyControlled1QubitGate(unsigned control, unsigned qubit, TM2x2<Type> const&m)
{
  assert(qubit < nqbits);

  if (fusion == true) {
    assert((*permutation)[qubit] < nqbits);
    if ((*permutation)[qubit] < log2llc) {
      std::string name = "cqg";
      fwindow.push_back(std::make_tuple(name, m, control, qubit));
      return;
    } else {
      applyFusedGates();
      goto L;
    }
  }
  L:
  applyControlled1QubitGate_helper(control, qubit, m, 0UL, localSize());
}


template <class Type>
void QbitRegister<Type>::applyCRotationX(unsigned const control, unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> rx;
  rx(0, 1) = rx(1, 0) = Type(0, -std::sin(theta / 2.));
  rx(0, 0) = rx(1, 1) = Type(std::cos(theta / 2.), 0);
  applyControlled1QubitGate(control, qubit, rx);

}

template <class Type>
void QbitRegister<Type>::applyCRotationY(unsigned const control, unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> ry;
  ry(0, 1) = Type(-std::sin(theta / 2.), 0.);
  ry(1, 0) = Type( std::sin(theta / 2.), 0.);
  ry(0, 0) = ry(1, 1) = Type(std::cos(theta / 2.), 0);
  applyControlled1QubitGate(control, qubit, ry);
}

template <class Type>
void QbitRegister<Type>::applyCRotationZ(unsigned const control, unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> rz;
  rz(0, 0) = Type(std::cos(theta / 2.), -std::sin(theta / 2.));
  rz(1, 1) = Type(std::cos(theta / 2.), std::sin(theta / 2.));
  rz(0, 1) = rz(1, 0) = Type(0., 0.);
  applyControlled1QubitGate(control, qubit, rz);
}

template <class Type>
void QbitRegister<Type>::applyCPauliX(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> px;
  px(0, 0) = Type(0., 0.);
  px(0, 1) = Type(1., 0.);
  px(1, 0) = Type(1., 0.);
  px(1, 1) = Type(0., 0.);
  applyControlled1QubitGate(control, qubit, px);
}

template <class Type>
void QbitRegister<Type>::applyCPauliY(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> py;
  py(0, 0) = Type(0., 0.);
  py(0, 1) = Type(0., -1.);
  py(1, 0) = Type(0., 1.);
  py(1, 1) = Type(0., 0.);
  applyControlled1QubitGate(control, qubit, py);
}

template <class Type>
void QbitRegister<Type>::applyCPauliZ(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(-1., 0.);
  applyControlled1QubitGate(control, qubit, pz);
}

template <class Type>
void QbitRegister<Type>::applyCPauliSqrtZ(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(0., 1.);
  applyControlled1QubitGate(control, qubit, pz);
}


template <class Type>
void QbitRegister<Type>::applyCHadamard(unsigned const control, unsigned const qubit)
{
  TM2x2<Type> h;
  BaseType f = 1. / std::sqrt(2.);
  h(0, 0) = h(0, 1) = h(1, 0) = Type(f, 0.);
  h(1, 1) = Type(-f, 0.);
  applyControlled1QubitGate(control, qubit, h);
}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;
