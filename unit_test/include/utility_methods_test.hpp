#ifndef UTILITY_METHODS_TEST_HPP
#define UTILITY_METHODS_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class UtilityMethodsTest : public ::testing::Test
{
 protected:

  UtilityMethodsTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 4-qubit state is distributed in more than 2^3 ranks.
    if (qhipster::mpi::Environment::GetStateSize() > 8)
      GTEST_SKIP();
  }

  const std::size_t num_qubits_ = 4;
  double accepted_error_ = 1e-15;
  const double inv_sqrt_2 = 1./std::sqrt(2.);
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(UtilityMethodsTest, DifferenceOfTwoStates)
{
  QubitRegister<ComplexDP> psi_a (num_qubits_,"base",0);
  std::size_t index = 4;
  QubitRegister<ComplexDP> psi_b (num_qubits_,"base",index);
  // |a> = |0000>
  // |b> = |0010> = |"4">

  ASSERT_DOUBLE_EQ( std::norm(psi_a.ComputeOverlap(psi_a)), 1. );
  ASSERT_DOUBLE_EQ( std::norm(psi_b.ComputeOverlap(psi_a)), 0. );
  ASSERT_DOUBLE_EQ(psi_b.MaxAbsDiff(psi_a), 1. );
  // The max L2-norm difference can be up to 4, reached for computational basis states
  // that differ by a (-1) global phase.
  // Othogonal computational basis states have max L-2 norm of 2.
  // These considerations are correct when StateSize=1, otherwise one has to correct
  // for the 'max' among all state ranks.
  if ( psi_b.LocalSize() > index )
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 2. );
  else
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 1. );

  psi_b.ApplyHadamard(2);
  // |a> = |0000>
  // |b> = |00-0>
  ASSERT_DOUBLE_EQ( std::norm(psi_b.ComputeOverlap(psi_a)), 0.5 );
  ASSERT_DOUBLE_EQ( psi_b.ComputeOverlap(psi_a).real(), inv_sqrt_2 );
  ASSERT_DOUBLE_EQ(psi_b.MaxAbsDiff(psi_a), inv_sqrt_2 );
  // Note that the only non-zero amplitudes have index 0 and 4, for both |a> and |b>.
  if ( psi_b.LocalSize() > 4 )
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), (1-inv_sqrt_2)*(1-inv_sqrt_2)+0.5 );
  else
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a),
                       std::max( (1-inv_sqrt_2)*(1-inv_sqrt_2), 0.5) );

  psi_a.ApplyPauliX(2);
  // |a> = |0010>
  // |b> = |00-0>
  ASSERT_DOUBLE_EQ( psi_b.ComputeOverlap(psi_a).real(), -inv_sqrt_2 );

  psi_a.Initialize("base",12);
  psi_b.Initialize("base",12);
  psi_b.ApplyPauliZ(3);
  // |a> =  |0011> =  |"12">
  // |b> = -|0011> = -|"12">
  ASSERT_DOUBLE_EQ( psi_b.ComputeOverlap(psi_a).real(), -1. );
  ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 4. );

  psi_a.ApplyHadamard(2);
  psi_a.ApplyHadamard(3);
  psi_b.ApplyHadamard(2);
  psi_b.ApplyHadamard(3);
  // |a> =  |00-->
  // |b> = -|00-->
  // Note that the only non-zero amplitudes have index 0,4,8,12.
  if ( psi_b.LocalSize() > 12 )
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 4. );
  else if ( psi_b.LocalSize() > 8 )
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 3. );
  else if ( psi_b.LocalSize() > 4 )
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 2. );
  else
      ASSERT_DOUBLE_EQ(psi_b.MaxL2NormDiff(psi_a), 1. );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(UtilityMethodsTest, AddOthers)
{
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard UTILITY_METHODS_TEST_HPP
