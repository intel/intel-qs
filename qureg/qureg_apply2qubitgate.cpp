//------------------------------------------------------------------------------
// Copyright 2017 Thomas Haener
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
void QbitRegister<Type>::apply2QubitGate(unsigned const qubit_high, unsigned const qubit_low, TM4x4<Type> const&m)
{
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  assert(nprocs == 1);

  std::size_t n = globalSize();
  assert(openqu::isPowerOf2(n));
  assert(qubit_low < openqu::highestBit(n));
  assert(qubit_high < openqu::highestBit(n));
  assert(qubit_low != qubit_high);

  std::size_t dh = 1UL << qubit_high;
  std::size_t dl = 1UL << qubit_low;

  bool const inv = qubit_low > qubit_high;

  // Does the compiler pull out the if-statements and thus perform the
  // branching outside of the loops when optimizing?
  for (std::size_t i = 0; i < n; i += 2 * (inv ? dl : dh)) {
    for (std::size_t j = 0; j < (inv ? dl : dh); j += 2 * (inv ? dh : dl)) {
      for (std::size_t k = 0; k < (inv ? dh : dl); ++k) {
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
}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;
