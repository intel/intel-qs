#ifndef ONE_QUBIT_REGISTER_TEST_HPP
#define ONE_QUBIT_REGISTER_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class OneQubitRegisterTest : public ::testing::Test
{
 protected:

  OneQubitRegisterTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 1-qubit state is distributed in more than 1 rank.
    // In fact the MPI version needs to allocate half-the-local-storage for communication.
    // If the local storage is a single amplitude, this cannot be further divided.
    if (iqs::mpi::Environment::GetStateSize() > 1)
      GTEST_SKIP();
  }

  // just before the 'destructor'
//  void TearDown() override {}

  const std::size_t num_qubits_ = 1;
  double accepted_error_ = 1e-15;
};
	
//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(OneQubitRegisterTest, InitializeWithDefault)
{
  ComplexDP amplitude;

  iqs::QubitRegister<ComplexDP> psi_0;
  // |psi_0> = |0>
  amplitude = psi_0.GetGlobalAmplitude(0);
  ASSERT_DOUBLE_EQ(amplitude.real(), 1.);
  ASSERT_DOUBLE_EQ(amplitude.imag(), 0.);
  amplitude = psi_0.GetGlobalAmplitude(1);
  ASSERT_DOUBLE_EQ(amplitude.real(), 0.);
  ASSERT_DOUBLE_EQ(amplitude.imag(), 0.);
}

TEST_F(OneQubitRegisterTest, InitializeInComputationalBasis)
{
  ComplexDP amplitude;

  iqs::QubitRegister<ComplexDP> psi_0 (num_qubits_,"base",0);
  // |psi_0> = |0>
  amplitude = psi_0.GetGlobalAmplitude(0);
  ASSERT_DOUBLE_EQ(amplitude.real(), 1.);
  ASSERT_DOUBLE_EQ(amplitude.imag(), 0.);
  amplitude = psi_0.GetGlobalAmplitude(1);
  ASSERT_DOUBLE_EQ(amplitude.real(), 0.);
  ASSERT_DOUBLE_EQ(amplitude.imag(), 0.);

  iqs::QubitRegister<ComplexDP> psi_1 (num_qubits_,"base",1);
  // |psi_1> = |1>
  amplitude = psi_1.GetGlobalAmplitude(0);
  ASSERT_DOUBLE_EQ(amplitude.real(), 0.);
  ASSERT_DOUBLE_EQ(amplitude.imag(), 0.);
  amplitude = psi_1.GetGlobalAmplitude(1);
  ASSERT_DOUBLE_EQ(amplitude.real(), 1.);
  ASSERT_DOUBLE_EQ(amplitude.imag(), 0.);
}

// FIXME until the RNG issues are solved.
#if 0
TEST_F(OneQubitRegisterTest, InitializeRandomly)
{
  psi_0_.Initialize("rand",1321);
  psi_1_.Initialize("rand",1157);

  ASSERT_LE( std::abs(psi_0_.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_1_.ComputeNorm()-1) , accepted_error_ );
}
#endif

TEST_F(OneQubitRegisterTest, GetCorrectProbability)
{
  iqs::QubitRegister<ComplexDP> psi_0 (num_qubits_,"base",0);
  ASSERT_DOUBLE_EQ(psi_0.GetProbability(0), 0.);
  psi_0.ApplyHadamard(0);
  // |psi_0> = |+>
  ASSERT_LE( std::abs(psi_0.GetProbability(0)-0.5) , accepted_error_ );

  iqs::QubitRegister<ComplexDP> psi_1 (num_qubits_,"base",1);
  ASSERT_DOUBLE_EQ(psi_1.GetProbability(0), 1.);
  psi_1.ApplyPauliX(0);
  // |psi_1> = |0>
  ASSERT_DOUBLE_EQ(psi_1.GetProbability(0), 0.);
  // |psi_1> = exp(-i theta X/2)|0> = cos(theta/2) |0> -i sin(theta/2)|1>
  double theta = 1.1;
  psi_1.ApplyRotationX(0,theta);
  ASSERT_LE( std::abs(psi_1.GetProbability(0)-std::sin(theta/2)*std::sin(theta/2)),
             accepted_error_ );
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard ONE_QUBIT_REGISTER_TEST_HPP
