//------------------------------------------------------------------------------
// Copyright (C) 2019 Intel Corporation 
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

#ifndef MEASURE_TEST_HPP
#define MEASURE_TEST_HPP

#include <math.h>

#include "../../qureg/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class MeasureTest : public ::testing::Test
{
 protected:

  MeasureTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 4-qubit state is distributed in more than 8 ranks.
    if (qhipster::mpi::Environment::GetStateSize() > 8)
      GTEST_SKIP();
  }

  const std::size_t num_qubits_ = 4;
  double accepted_error_ = 1e-14;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(MeasureTest, GetProbability)
{
  QubitRegister<ComplexDP> psi (num_qubits_,"base",10);
  psi.ApplyHadamard(2);
  psi.ApplyHadamard(3);
  // |psi> = |01+-> = H2.H3|"2+8">
  ASSERT_DOUBLE_EQ(psi.GetProbability(0), 0.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(1), 1.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(2), 0.5);
  ASSERT_DOUBLE_EQ(psi.GetProbability(3), 0.5);

  psi.Initialize("base",0);
  psi.ApplyRotationX(0,-M_PI/2.);
  // |psi> = |y+>|000>
  ASSERT_DOUBLE_EQ(psi.GetProbability(0), 0.5);

  psi.Initialize("base",0);
  psi.ApplyRotationX(0, M_PI/2.);
  // |psi> = |y->|000>
  ASSERT_DOUBLE_EQ(psi.GetProbability(0), 0.5);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard MEASURE_TEST_HPP
