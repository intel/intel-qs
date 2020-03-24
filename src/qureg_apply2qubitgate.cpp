#include "../include/qureg.hpp"
#include "../include/highperfkernels.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_apply2qubitgate.cpp
/// @brief Define the @c QubitRegister methods corresponding to the application of two-qubit gates.

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Arbitrary two-qubit gate.
/// @param qubit_high index of the first qubit
/// @param qubit_low index of the second qubit
/// @param m 4x4 matrix corresponding to the quantum gate
template <class Type>
void QubitRegister<Type>::Apply2QubitGate(unsigned const qubit_high, unsigned const qubit_low,
                                          TM4x4<Type> const&m)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  // The general case has not been implemented yet for distributed computation.
  assert(nprocs == 1);

  std::size_t n = GlobalSize();
  assert(qhipster::isPowerOf2(n));
  assert(qubit_low < qhipster::highestBit(n));
  assert(qubit_high < qhipster::highestBit(n));
  assert(qubit_low != qubit_high);

  std::size_t dh = 1UL << qubit_high;
  std::size_t dl = 1UL << qubit_low;

  // it is easier if the qubit indices are ordered, meaning that:
  //     qubit_high > qubit_low
  // if not, we take this into account with the variable 'invert'
  bool const invert = qubit_low > qubit_high;

  // Does the compiler pull out the if-statements and thus perform the
  // branching outside of the loops when optimizing?
  for (std::size_t i = 0; i < n; i += 2 * (invert ? dl : dh))
  {
      for (std::size_t j = 0; j < (invert ? dl : dh); j += 2 * (invert ? dh : dl))
      {
          for (std::size_t k = 0; k < (invert ? dh : dl); ++k)
          {
              auto tmp = i + j + k;
              std::size_t idx[] = {tmp, tmp + dl, tmp + dh, tmp + dh + dl};

              auto uu = state[idx[0]];  // 00
              auto ud = state[idx[1]];  // 01
              auto du = state[idx[2]];  // 10
              auto dd = state[idx[3]];  // 11

              for (unsigned t = 0; t < 4; ++t)
                  state[idx[t]] = m[t][0] * uu + m[t][1] * ud + m[t][2] * du + m[t][3] * dd;
          }
      }
  }

  // Update counter of the statistics.
  if (gate_counter != nullptr)
  {
      // Verify that permutation is identity.
      assert(qubit_high == (*permutation)[qubit_high]);
      assert(qubit_low  == (*permutation)[qubit_low ]);
      // Otherwise find better location that is compatible with the permutation.
      gate_counter->TwoQubitIncrement(qubit_high, qubit_low);
  }
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
