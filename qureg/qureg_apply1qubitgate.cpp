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

/** \addtogroup qureg
 *  @{
 */

/** @file qureg_apply1qubitgate.cpp
 *  @brief Define the @c QbitRegister methods corresponding to the application of single-qubit gates.
 */

template <class Type>
double QbitRegister<Type>::HP_Distrpair(unsigned pos, TM2x2<Type> const&m)
{
#ifdef OPENQU_HAVE_MPI
  MPI_Status status;
  MPI_Comm comm = openqu::mpi::Environment::comm();
  std::size_t myrank = openqu::mpi::Environment::rank();


  assert(pos < nqbits);
  int strideexp = pos;
  int memexp = nqbits - openqu::ilog2(openqu::mpi::Environment::size());
  int pstrideexp = strideexp - memexp;

  /*
    Steps:     1.         2.           3.              4.
            i    j   | i    j   |  i      j     |  i       j
            s1   d1  | s1   d1  |  s1     d1&s1 |  s1&d1   d1&s1
            s2   d2  | s2   d2  |  s2&d2  d2    |  s2&d2   d2&s2
            T    T   | d2   s1  |  d2&s2  s1&d1 |  T       T

  */

  int tag1 = 1, tag2 = 2;
  int tag3 = 3, tag4 = 4;
  std::size_t glb_start = UL(myrank) * localSize();

  // std::string s;
  unsigned int itask, jtask;
  if (check_bit(glb_start, UL(pos)) == 0) {
    itask = myrank;
    jtask = itask + (1 << pstrideexp);
    // s = openqu::toString(itask) + "==>" + openqu::toString(jtask);
  } else {
    jtask = myrank;
    itask = jtask - (1 << pstrideexp);
    // s = openqu::toString(jtask) + "==>" + openqu::toString(itask);
  }

  // openqu::mpi::print(s, true);

  if (specialize == true) { 
    // check for special case of diagonal
    bool Xgate = (m[0][0] == Type(0., 0.) && m[0][1] == Type(1., 0.) &&
                  m[1][0] == Type(1., 0.) && m[1][1] == Type(0., 0.));
    if (Xgate == true) {
      // printf("Xgate: remaping MPI rank %d <==> %d\n", jtask, itask);
      if (check_bit(glb_start, UL(pos)) == 0) {
        openqu::mpi::Environment::remaprank(jtask);
      } else {
        openqu::mpi::Environment::remaprank(itask);
      }
      TODO(Fix problem when coming here from controlled gate)
      openqu::mpi::barrier();
      if (timer) timer->record_cm(0., 0.);
      return 0.0;
    }
    bool Ygate = (m[0][0] == Type(0., 0.) && m[0][1] == Type(0., -1.) &&
                  m[1][0] == Type(0., 1.) && m[1][1] == Type(0., 0.));
    if (Ygate == true) {
      // printf("Ygate: remaping MPI rank\n");
      if (check_bit(glb_start, UL(pos)) == 0) {
        openqu::mpi::Environment::remaprank(jtask);
        ScaleState(0UL, localSize(), state, Type(0, 1.0), timer);
      } else {
        openqu::mpi::Environment::remaprank(itask);
        ScaleState(0UL, localSize(), state, Type(0, -1.0), timer);
      }
      openqu::mpi::barrier();
      if (timer) timer->record_cm(0., 0.);
      return 0.0;
    }

  }


  // 1. allocate temp buffer
  size_t lcl_chunk = tmpSize();
  Type *tmpstate = tmpspace();

  std::size_t lcl_size_half = localSize() / 2L;
  assert(lcl_size_half <= std::numeric_limits<int>::max());

  if (lcl_chunk > lcl_size_half) 
    lcl_chunk = lcl_size_half;
  else
    assert((lcl_size_half % lcl_chunk) == 0);
  
  double t, tnet = 0;
  for(size_t c = 0; c < lcl_size_half; c += lcl_chunk) {
    // if(!myrank) printf("c=%lu lcl_size_half=%lu lcl_chujnk=%lu\n", c, lcl_size_half, lcl_chunk);
    if (itask == myrank) {  // this is itask
      // 2. src sends s1 to dst into dT
      //    dst sends d2 to src into dT
      t = sec();
      MPI_Sendrecv_x(&(state[c]), lcl_chunk, jtask, tag1, &(tmpstate[0]),
                     lcl_chunk, jtask, tag2, comm, &status);
      tnet += sec() - t;

      // 3. src and dst compute
      Loop_SN(0L, lcl_chunk, &(state[c]), tmpstate, lcl_size_half, 0L, m, specialize, timer);

      t = sec();
      MPI_Sendrecv_x(&(tmpstate[0]), lcl_chunk, jtask, tag3, &(state[c]),
                     lcl_chunk, jtask, tag4, comm, &status);
      tnet += sec() - t;

    } else {  // this is jtask
      // 2. src sends s1 to dst into dT
      //    dst sends d2 to src into dT
      t = sec();
      MPI_Sendrecv_x(&(state[lcl_size_half + c]), lcl_chunk, itask, tag2,
                     &(tmpstate[0]), lcl_chunk, itask, tag1, comm,
                     &status);
      tnet += sec() - t;

      Loop_SN(0L, lcl_chunk, tmpstate, &(state[c]), 0L, 0L, m, specialize, timer);

      t = sec();
      MPI_Sendrecv_x(&(tmpstate[0]), lcl_chunk, itask, tag4,
                     &(state[lcl_size_half + c]), lcl_chunk, itask, tag3,
                     comm, &status);
      tnet += sec() - t;
    }
  }

  double netsize = 2.0 * sizeof(Type) * 2.0 * D(lcl_size_half), netbw = netsize / tnet;
  // printf("[%3d] size=%10lld tnet = %.3lf s netsize = %10.0lf bytes netbw = %6.2lf GB/s\n",
  //      it, sizeof(Type)*lcl_size_half, tnet, netsize, netbw / 1e9);

  if (timer) {timer->record_cm(tnet, netbw); };

#else
  assert(0);
#endif

  return 0.0;
}

template <class Type>
bool QbitRegister<Type>::apply1QubitGate_helper(unsigned qubit_,  TM2x2<Type> const&m, 
                                          std::size_t sind, std::size_t eind)
{
  assert(qubit_ < nqbits);
  unsigned qubit = (*permutation)[qubit_]; 
  assert(qubit < nqbits);

  TODO(Add diagonal special case)

  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  unsigned log2_nprocs = openqu::ilog2(openqu::mpi::Environment::size());
  unsigned M = nqbits - log2_nprocs;
  std::size_t P = qubit;

  std::size_t src_glb_start = UL(myrank) * localSize();
  // check for special case of diagonal
  bool diagonal = (m[0][1].real() == 0. && m[0][1].imag() == 0. &&
                   m[1][0].real() == 0. && m[1][0].imag() == 0.);

  std::string gname = "SQG("+openqu::toString(P)+")::"+m.name;
  if(timer) timer->Start(gname, P);

  if(P < M)
  {
     assert(eind - sind <= localSize());
     Loop_DN(sind, eind, UL(P), state, state, 0UL, (1UL << P), m, specialize, timer);
  }
  else
  {
     assert(eind - sind == localSize());
     if (specialize && diagonal) {
       if (check_bit(src_glb_start, P) == 0 ) {
         ScaleState(sind, eind, state, m[0][0], timer);
       } else {
         ScaleState(sind, eind, state, m[1][1], timer);
       }
     } else {
       HP_Distrpair(P, m);
     }
  }
  if(timer) timer->Stop();
  

  return true;
}

template <class Type>
void QbitRegister<Type>::apply1QubitGate(unsigned qubit, TM2x2<Type> const&m)
{
  if (fusion == true) {
    assert((*permutation)[qubit] < nqbits);
    if ((*permutation)[qubit] < log2llc) {
      std::string name = "sqg";
      fwindow.push_back(std::make_tuple(name, m, qubit, 0U));
      return;
    } else {
      applyFusedGates();
      goto L;
    }
  }

  L:
  apply1QubitGate_helper(qubit, m, 0UL, localSize());
}



/**
 * @brief Rotation around the X axis by an angle theta
 * @param qubit index of the involved qubit
 * @param theta rotation angle
 *
 * Explicitely, the gate corresponds to:\n
 *     exp( -i X theta/2 )\n
 * This convention is based on the fact that the generators
 * of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
 */
template <class Type>
void QbitRegister<Type>::applyRotationX(unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> rx;
  rx(0, 1) = rx(1, 0) = Type(0, -std::sin(theta / 2.));
  rx(0, 0) = rx(1, 1) = std::cos(theta / 2.);
  apply1QubitGate(qubit, rx);

}

/**
 * @brief Rotation around the Y axis by an angle theta
 * @param qubit index of the involved qubit
 * @param theta rotation angle
 *
 * Explicitely, the gate corresponds to:\n
 *     exp( -i Y theta/2 )\n
 * This convention is based on the fact that the generators
 * of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
 */
template <class Type>
void QbitRegister<Type>::applyRotationY(unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> ry;
  ry(0, 1) = Type(-std::sin(theta / 2.), 0.);
  ry(1, 0) = Type( std::sin(theta / 2.), 0.);
  ry(0, 0) = ry(1, 1) = std::cos(theta / 2.);
  apply1QubitGate(qubit, ry);
}

/**
 * @brief Rotation around the Z axis by an angle theta
 * @param qubit index of the involved qubit
 * @param theta rotation angle
 *
 * Explicitely, the gate corresponds to:\n
 *     exp( -i Z theta/2 )\n
 * This convention is based on the fact that the generators
 * of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
 */
template <class Type>
void QbitRegister<Type>::applyRotationZ(unsigned const qubit, BaseType theta)
{
  openqu::TinyMatrix<Type, 2, 2, 32> rz;
  rz(0, 0) = Type(std::cos(theta / 2.), -std::sin(theta / 2.));
  rz(1, 1) = Type(std::cos(theta / 2.), std::sin(theta / 2.));
  rz(0, 1) = rz(1, 0) = Type(0., 0.);
  apply1QubitGate(qubit, rz);
}

/**
 * @brief Apply X Pauli operator
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to:\n
 *     i * exp( -i X pi/2 ) = X 
 */
template <class Type>
void QbitRegister<Type>::applyPauliX(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> px;
  px(0, 0) = Type(0., 0.);
  px(0, 1) = Type(1., 0.);
  px(1, 0) = Type(1., 0.);
  px(1, 1) = Type(0., 0.);
  apply1QubitGate(qubit, px);
}

/**
 * @brief Apply square root of the X Pauli operator
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to:\n
 *     sqrt(X) 
 */
template <class Type>
void QbitRegister<Type>::applyPauliSqrtX(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> px;
  px(0, 0) = Type(0.5,  0.5);
  px(0, 1) = Type(0.5, -0.5);
  px(1, 0) = Type(0.5, -0.5);
  px(1, 1) = Type(0.5,  0.5);
  apply1QubitGate(qubit, px);
}

/**
 * @brief Apply Y Pauli operator
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to:\n
 *     i * exp( -i Y pi/2 ) = Y 
 */
template <class Type>
void QbitRegister<Type>::applyPauliY(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> py;
  py(0, 0) = Type(0., 0.);
  py(0, 1) = Type(0., -1.);
  py(1, 0) = Type(0., 1.);
  py(1, 1) = Type(0., 0.);
  apply1QubitGate(qubit, py);
}

/**
 * @brief Apply square root of the Y Pauli operator
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to:\n
 *     sqrt(Y) 
 */
template <class Type>
void QbitRegister<Type>::applyPauliSqrtY(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> py;
  py(0, 0) = Type(0.5,   0.5);
  py(0, 1) = Type(-0.5, -0.5);
  py(1, 0) = Type(0.5,   0.5);
  py(1, 1) = Type(0.5,  0.5);
  apply1QubitGate(qubit, py);
}

/**
 * @brief Apply Z Pauli operator
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to:\n
 *     i * exp( -i Z pi/2 ) = Z 
 */
template <class Type>
void QbitRegister<Type>::applyPauliZ(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(-1., 0.);
  apply1QubitGate(qubit, pz);
}

/**
 * @brief Apply square root of the Z Pauli operator
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to:\n
 *     sqrt(Z) 
 */
template <class Type>
void QbitRegister<Type>::applyPauliSqrtZ(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(0., 1.);
  apply1QubitGate(qubit, pz);
}


/**
 * @brief Apply Hadamard gate
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to the 2x2 matrix:\n
 *     | 1/sqrt(2)   1/sqrt(2) |\n
 *     | 1/sqrt(2)  -1/sqrt(2) |
 */
template <class Type>
void QbitRegister<Type>::applyHadamard(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> h;
  BaseType f = 1. / std::sqrt(2.);
  h(0, 0) = h(0, 1) = h(1, 0) = Type(f, 0.);
  h(1, 1) = Type(-f, 0.);
  apply1QubitGate(qubit, h);
}

/**
 * @brief Apply T gate
 * @param qubit index of the involved qubit
 *
 * Explicitely, the gate corresponds to the 2x2 matrix:\n
 *     | 1              0           |\n
 *     | 0    cos(pi/4)+i*sin(pi/4) |
 */
template <class Type>
void QbitRegister<Type>::applyT(unsigned const qubit)
{
  openqu::TinyMatrix<Type, 2, 2, 32> t;
  t(0, 0) = Type(1.0, 0.0);
  t(0, 1) = Type(0.0, 0.0);
  t(1, 0) = Type(0.0, 0.0);
  t(1, 1) = Type(cos(M_PI/4.0), sin(M_PI/4.0));
  apply1QubitGate(qubit, t);

}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;

/** @}*/
