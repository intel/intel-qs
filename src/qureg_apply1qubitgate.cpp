/// @file qureg_apply1qubitgate.cpp
/// @brief Define the @c QubitRegister methods corresponding to the application of single-qubit gates.

#include "../include/qureg.hpp"
#include "../include/highperfkernels.hpp"
#include "../include/spec_kernels.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// General comment.
// To distinguish between program qubits (used in the algorithm) and data qubits
// (used in the representation of the quantum state), we use the term:
// - 'position' to refer to data qubits
// - 'qubit' to refer to program qubits
/////////////////////////////////////////////////////////////////////////////////////////

namespace iqs {

template <class Type>
double QubitRegister<Type>::HP_Distrpair(unsigned position, TM2x2<Type> const&m, GateSpec1Q spec, BaseType angle)
{
  assert(LocalSize() > 1);
#ifndef INTELQS_HAS_MPI
  assert(0);
#else
  MPI_Status status;
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  std::size_t myrank = iqs::mpi::Environment::GetStateRank();

  assert(position < num_qubits);
  int strideexp = position;
  int memexp = num_qubits - iqs::ilog2(iqs::mpi::Environment::GetStateSize());
  int pstrideexp = strideexp - memexp;

  //  Steps:     1.         2.           3.              4.
  //          i    j   | i    j   |  i      j     |  i       j
  //          s1   d1  | s1   d1  |  s1     d1&s1 |  s1&d1   d1&s1
  //          s2   d2  | s2   d2  |  s2&d2  d2    |  s2&d2   d2&s2
  //          T    T   | d2   s1  |  d2&s2  s1&d1 |  T       T

  int tag1 = 1, tag2 = 2;
  int tag3 = 3, tag4 = 4;
  std::size_t glb_start = UL(myrank) * LocalSize();

  // std::string s;
  unsigned int itask, jtask;
  if (check_bit(glb_start, UL(position)) == 0)
  {
      itask = myrank;
      jtask = itask + (1 << pstrideexp);
      // s = iqs::toString(itask) + "==>" + iqs::toString(jtask);
  }
  else
  {
      jtask = myrank;
      itask = jtask - (1 << pstrideexp);
      // s = iqs::toString(jtask) + "==>" + iqs::toString(itask);
  }

  // iqs::mpi::StatePrint(s, true);

  if (specialize == true)
  { 
      // check for special case of diagonal
      bool Xgate = (m[0][0] == Type(0., 0.) && m[0][1] == Type(1., 0.) &&
                    m[1][0] == Type(1., 0.) && m[1][1] == Type(0., 0.));
      if (Xgate == true)
      {
          // printf("Xgate: remaping MPI rank %d <==> %d\n", jtask, itask);
          if (check_bit(glb_start, UL(position)) == 0)
              iqs::mpi::Environment::RemapStateRank(jtask);
          else
              iqs::mpi::Environment::RemapStateRank(itask);
          // TODO: Fix problem when coming here from controlled gate.
          iqs::mpi::StateBarrier();
          if (timer)
              timer->record_cm(0., 0.);
          return 0.0;
      }
      bool Ygate = (m[0][0] == Type(0., 0.) && m[0][1] == Type(0., -1.) &&
                    m[1][0] == Type(0., 1.) && m[1][1] == Type(0., 0.));
      if (Ygate == true)
      {
          // printf("Ygate: remaping MPI rank\n");
          if (check_bit(glb_start, UL(position)) == 0)
          {
              iqs::mpi::Environment::RemapStateRank(jtask);
              ScaleState(0UL, LocalSize(), state, Type(0, 1.0), timer);
          }
          else
          {
              iqs::mpi::Environment::RemapStateRank(itask);
              ScaleState(0UL, LocalSize(), state, Type(0, -1.0), timer);
          }
          iqs::mpi::StateBarrier();
          if (timer)
              timer->record_cm(0., 0.);
          return 0.0;
      }
  }


  // 1. allocate temp buffer
  size_t lcl_chunk = TmpSize();
  Type *tmp_state = TmpSpace();

  std::size_t lcl_size_half = LocalSize() / 2L;
  assert(lcl_size_half <= std::numeric_limits<size_t>::max());

  if (lcl_chunk > lcl_size_half) 
      lcl_chunk = lcl_size_half;
  else
      assert((lcl_size_half % lcl_chunk) == 0);
  
  double t, tnet = 0;
  for(std::size_t c = 0; c < lcl_size_half; c += lcl_chunk)
  {
    // if(!myrank) printf("c=%lu lcl_size_half=%lu lcl_chujnk=%lu\n", c, lcl_size_half, lcl_chunk);
    if (itask == myrank)  // this is itask
    {
        // 2. src sends s1 to dst into dT
        //    dst sends d2 to src into dT
        t = sec();
        iqs::mpi::MPI_Sendrecv_x(&(state[c])    , lcl_chunk, jtask, tag1,
                                      &(tmp_state[0]), lcl_chunk, jtask, tag2,
                                      comm, &status);
        tnet += sec() - t;

        // 3. src and dst compute
        if (specialize2 && spec != GateSpec1Q::None)
          Loop_SN(0L, lcl_chunk, &(state[c]), tmp_state, lcl_size_half, 0L, spec, timer, angle);
        else
          Loop_SN(0L, lcl_chunk, &(state[c]), tmp_state, lcl_size_half, 0L, m, specialize, timer);

        t = sec();
        iqs::mpi::MPI_Sendrecv_x(&(tmp_state[0]), lcl_chunk, jtask, tag3,
                                      &(state[c])    , lcl_chunk, jtask, tag4,
                                      comm, &status);
        tnet += sec() - t;
    }
    else  // this is jtask
    {
        // 2. src sends s1 to dst into dT
        //    dst sends d2 to src into dT
        t = sec();
        iqs::mpi::MPI_Sendrecv_x(&(state[lcl_size_half + c]), lcl_chunk, itask, tag2,
                                      &(tmp_state[0])            , lcl_chunk, itask, tag1,
                                      comm, &status);
        tnet += sec() - t;
        if (specialize2 && spec != GateSpec1Q::None)
          Loop_SN(0L, lcl_chunk, tmp_state, &(state[c]), 0L, 0L, spec, timer, angle);
        else
          Loop_SN(0L, lcl_chunk, tmp_state, &(state[c]), 0L, 0L, m, specialize, timer);

        t = sec();
        iqs::mpi::MPI_Sendrecv_x(&(tmp_state[0])            , lcl_chunk, itask, tag4,
                                      &(state[lcl_size_half + c]), lcl_chunk, itask, tag3,
                                      comm, &status);
        tnet += sec() - t;
    }
  }

  double netsize = 2.0 * sizeof(Type) * 2.0 * double(lcl_size_half), netbw = netsize / tnet;
  // printf("[%3d] size=%10lld tnet = %.3lf s netsize = %10.0lf bytes netbw = %6.2lf GB/s\n",
  //      it, sizeof(Type)*lcl_size_half, tnet, netsize, netbw / 1e9);

  if (timer) {timer->record_cm(tnet, netbw); };
#endif
  return 0.0;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
bool QubitRegister<Type>::Apply1QubitGate_helper(unsigned qubit_,  TM2x2<Type> const&m, 
                                                 std::size_t sind, std::size_t eind,
                                                 GateSpec1Q spec, BaseType angle)
{
  assert(qubit_ < num_qubits);
  unsigned position = (*qubit_permutation)[qubit_]; 
  assert(position < num_qubits);

  // TODO: Add diagonal special case.

  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
  unsigned M = num_qubits - log2_nprocs;
  std::size_t P = position;

  std::size_t src_glb_start = UL(myrank) * LocalSize();
  // check for special case of diagonal
  bool diagonal = (m[0][1].real() == 0. && m[0][1].imag() == 0. &&
                   m[1][0].real() == 0. && m[1][0].imag() == 0.);

  std::string gate_name = "SQG("+iqs::toString(P)+")::"+m.name;

  if (timer)
      timer->Start(gate_name, P);

  if (P < M)
  {
      assert(eind - sind <= LocalSize());
      // Introduce specialized kernel here
      if (!specialize2 || (spec == GateSpec1Q::None))
	      Loop_DN(sind, eind, UL(P), state, state, 0UL, (1UL << P), m, specialize, timer);
      else
	      Loop_DN(sind, eind, UL(P), state, state, 0UL, (1UL << P), spec, timer, angle);
  }
  else
  {
      assert(eind - sind == LocalSize());
      if (specialize && diagonal)
      {
          if (check_bit(src_glb_start, P) == 0 )
              ScaleState(sind, eind, state, m[0][0], timer);
          else
              ScaleState(sind, eind, state, m[1][1], timer);
      }
      else
      {
          HP_Distrpair(P, m, spec, angle);
      }
  }

  if (timer)
      timer->Stop();
  
  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Apply1QubitGate(unsigned qubit, TM2x2<Type> const&m, GateSpec1Q spec, BaseType angle)
{
  // Update counter of the statistics.
  if (gate_counter != nullptr)
  {
      // IQS count the gates acting on specific program qubits.
      gate_counter->OneQubitIncrement(qubit);
  }

  unsigned position = (*qubit_permutation)[qubit];
  assert(position < num_qubits);

  // FIXME verify that fusion is properly working even with non-identity qubit permutation.
  if (fusion == true)
  {
      if (position < log2llc)
      {
          std::string name = "sqg";
          fwindow.push_back(std::make_tuple(name, m, qubit, 0U)); // FIXME: check if using qubit (original) or position
          return;
      }
      else
      {
          ApplyFusedGates();
          goto L;
      }
  }

  L:
  Apply1QubitGate_helper(qubit, m, 0UL, LocalSize(), spec, angle);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Rotation around the X axis by an angle theta
/// @param qubit index of the involved qubit
/// @param theta rotation angle
///
/// Explicitly, the gate corresponds to:\n
///     exp( -i X theta/2 )\n
/// This convention is based on the fact that the generators
/// of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
template <class Type>
void QubitRegister<Type>::ApplyRotationX(unsigned const qubit, BaseType theta)
{
  iqs::TinyMatrix<Type, 2, 2, 32> rx;
  rx(0, 1) = rx(1, 0) = Type(0, -std::sin(theta / 2.));
  rx(0, 0) = rx(1, 1) = std::cos(theta / 2.);
  Apply1QubitGate(qubit, rx, GateSpec1Q::RotationX, theta);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Rotation around the Y axis by an angle theta
/// @param qubit index of the involved qubit
/// @param theta rotation angle
///
/// Explicitly, the gate corresponds to:\n
///     exp( -i Y theta/2 )\n
/// This convention is based on the fact that the generators
/// of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
template <class Type>
void QubitRegister<Type>::ApplyRotationY(unsigned const qubit, BaseType theta)
{
  iqs::TinyMatrix<Type, 2, 2, 32> ry;
  ry(0, 1) = Type(-std::sin(theta / 2.), 0.);
  ry(1, 0) = Type( std::sin(theta / 2.), 0.);
  ry(0, 0) = ry(1, 1) = std::cos(theta / 2.);
  Apply1QubitGate(qubit, ry, GateSpec1Q::RotationY, theta);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Rotation around the Z axis by an angle theta
/// @param qubit index of the involved qubit
/// @param theta rotation angle
///
/// Explicitly, the gate corresponds to:\n
///     exp( -i Z theta/2 )\n
/// This convention is based on the fact that the generators
/// of rotations for spin-1/2 spins are {X/2, Y/2, Z/2}.
template <class Type>
void QubitRegister<Type>::ApplyRotationZ(unsigned const qubit, BaseType theta)
{
  iqs::TinyMatrix<Type, 2, 2, 32> rz;
  rz(0, 0) = Type(std::cos(theta / 2.), -std::sin(theta / 2.));
  rz(1, 1) = Type(std::cos(theta / 2.), std::sin(theta / 2.));
  rz(0, 1) = rz(1, 0) = Type(0., 0.);
  Apply1QubitGate(qubit, rz, GateSpec1Q::RotationZ, theta);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply X Pauli operator
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to:\n
///     i * exp( -i X pi/2 ) = X 
template <class Type>
void QubitRegister<Type>::ApplyPauliX(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> px;
  px(0, 0) = Type(0., 0.);
  px(0, 1) = Type(1., 0.);
  px(1, 0) = Type(1., 0.);
  px(1, 1) = Type(0., 0.);
  Apply1QubitGate(qubit, px, GateSpec1Q::PauliX);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply square root of the X Pauli operator
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to:\n
///     sqrt(X) 
template <class Type>
void QubitRegister<Type>::ApplyPauliSqrtX(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> px;
  px(0, 0) = Type(0.5,  0.5);
  px(0, 1) = Type(0.5, -0.5);
  px(1, 0) = Type(0.5, -0.5);
  px(1, 1) = Type(0.5,  0.5);
  Apply1QubitGate(qubit, px);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply Y Pauli operator
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to:\n
///     i * exp( -i Y pi/2 ) = Y 
template <class Type>
void QubitRegister<Type>::ApplyPauliY(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> py;
  py(0, 0) = Type(0., 0.);
  py(0, 1) = Type(0., -1.);
  py(1, 0) = Type(0., 1.);
  py(1, 1) = Type(0., 0.);
  Apply1QubitGate(qubit, py, GateSpec1Q::PauliY);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply square root of the Y Pauli operator
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to:\n
///     sqrt(Y) 
template <class Type>
void QubitRegister<Type>::ApplyPauliSqrtY(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> py;
  py(0, 0) = Type(0.5,   0.5);
  py(0, 1) = Type(-0.5, -0.5);
  py(1, 0) = Type(0.5,   0.5);
  py(1, 1) = Type(0.5,  0.5);
  Apply1QubitGate(qubit, py);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply Z Pauli operator
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to:\n
///     i * exp( -i Z pi/2 ) = Z 
template <class Type>
void QubitRegister<Type>::ApplyPauliZ(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(-1., 0.);
  Apply1QubitGate(qubit, pz, GateSpec1Q::PauliZ);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply square root of the Z Pauli operator
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to:\n
///     sqrt(Z) 
template <class Type>
void QubitRegister<Type>::ApplyPauliSqrtZ(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> pz;
  pz(0, 0) = Type(1., 0.);
  pz(0, 1) = Type(0., 0.);
  pz(1, 0) = Type(0., 0.);
  pz(1, 1) = Type(0., 1.);
  Apply1QubitGate(qubit, pz);
}



/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply Hadamard gate
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to the 2x2 matrix:\n
///     | 1/sqrt(2)   1/sqrt(2) |\n
///     | 1/sqrt(2)  -1/sqrt(2) |
template <class Type>
void QubitRegister<Type>::ApplyHadamard(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> h;
  BaseType f = 1. / std::sqrt(2.);
  h(0, 0) = h(0, 1) = h(1, 0) = Type(f, 0.);
  h(1, 1) = Type(-f, 0.);
  Apply1QubitGate(qubit, h, GateSpec1Q::Hadamard);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply rotation in the XY-plane
/// @param qubit index of the involved qubit
/// @param phi angle of the rotation axis w.r.t. X axis
/// @param theta rotation angle
///
/// Explicitly, the gate corresponds to:\n
///
///     R_XY(phi, theta) = cos(theta/2) I -i sin(theta/2) (cos(phi) X + sin(phi) Y)
///
///                      = | c(t/2)         -i s(t/2) (c(p) -i s(p)) |
///                        | -i s(t/2) (c(p) +i s(p))         c(t/2) |
///
/// Or, in other format:
///
///     R_XY(phi, theta) = exp( -i theta/2 P)
/// 
/// with P = cos(phi) X + sin(phi) Y , which can be seen as a rotated Pauli matrix.
template <class Type>
void QubitRegister<Type>::ApplyRotationXY(unsigned const qubit, BaseType phi, BaseType theta)
{
  iqs::TinyMatrix<Type, 2, 2, 32> rxy;
  rxy(0, 0) = Type(std::cos(theta / 2.), 0);
  rxy(0, 1) = Type(-std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi) );
  rxy(1, 0) = Type( std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi) );
  rxy(1, 1) = Type(std::cos(theta / 2.), 0);
  Apply1QubitGate(qubit, rxy);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply T gate
/// @param qubit index of the involved qubit
///
/// Explicitly, the gate corresponds to the 2x2 matrix:\n
///     | 1              0           |\n
///     | 0    cos(pi/4)+i*sin(pi/4) |
template <class Type>
void QubitRegister<Type>::ApplyT(unsigned const qubit)
{
  iqs::TinyMatrix<Type, 2, 2, 32> t;
  t(0, 0) = Type(1.0, 0.0);
  t(0, 1) = Type(0.0, 0.0);
  t(1, 0) = Type(0.0, 0.0);
  t(1, 1) = Type(cos(M_PI/4.0), sin(M_PI/4.0));
  Apply1QubitGate(qubit, t, GateSpec1Q::T);

}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

} // end namespace iqs
