#ifndef MEASURE_TEST_HPP
#define MEASURE_TEST_HPP

#include <math.h>

#include "../../include/qureg.hpp"

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
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 4-qubit state is distributed in more than 8 ranks.
    if (iqs::mpi::Environment::GetStateSize() > 8)
      GTEST_SKIP();
  }

  const std::size_t num_qubits_ = 4;
  double accepted_error_ = 1e-14;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(MeasureTest, GetProbability)
{
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",10);
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
