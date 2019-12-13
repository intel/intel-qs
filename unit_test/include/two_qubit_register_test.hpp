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

#ifndef TWO_QUBIT_REGISTER_TEST_HPP
#define TWO_QUBIT_REGISTER_TEST_HPP

#include "../../qureg/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class TwoQubitRegisterTest : public ::testing::Test
{
 protected:

  TwoQubitRegisterTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 2-qubit state is distributed in more than 2 ranks.
    // In fact the MPI version needs to allocate half-the-local-storage for communication.
    // If the local storage is a single amplitude, this cannot be further divided.
    if (qhipster::mpi::Environment::GetStateSize() > 2)
      GTEST_SKIP();
  }

  const std::size_t num_qubits_ = 2;
  double accepted_error_ = 1e-15;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(TwoQubitRegisterTest, InitializeInComputationalBasis)
{
  // |psi_0> = |00> = |q0=0> x |q1=0>
  QubitRegister<ComplexDP> psi_0 (num_qubits_,"base",0);
  // |psi_1> = |10> = |q0=1> x |q1=0>
  QubitRegister<ComplexDP> psi_1 (num_qubits_,"base",1);
  // |psi_2> = |01> = |q0=0> x |q1=1>
  QubitRegister<ComplexDP> psi_2 (num_qubits_,"base",2);
  // |psi_3> = |11> = |q0=1> x |q1=1>
  QubitRegister<ComplexDP> psi_3 (num_qubits_,"base",3);

  // Test the norm.
  ASSERT_LE( std::abs(psi_0.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_1.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_2.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_3.ComputeNorm()-1) , accepted_error_ );

  ComplexDP amplitude;
  // Test the only non-zero entry.
  amplitude = psi_0.GetGlobalAmplitude(0);
  ASSERT_DOUBLE_EQ(psi_0.GetGlobalAmplitude(0).real(), 1.);
  ASSERT_DOUBLE_EQ(psi_1.GetGlobalAmplitude(1).real(), 1.);
  ASSERT_DOUBLE_EQ(psi_2.GetGlobalAmplitude(2).real(), 1.);
  ASSERT_DOUBLE_EQ(psi_3.GetGlobalAmplitude(3).real(), 1.);

  // Test the probabilities (of being in |1>) for each qubit separately.
  ASSERT_DOUBLE_EQ(psi_0.GetProbability(0), 0.);
  ASSERT_DOUBLE_EQ(psi_0.GetProbability(1), 0.);
  ASSERT_DOUBLE_EQ(psi_1.GetProbability(0), 1.);
  ASSERT_DOUBLE_EQ(psi_1.GetProbability(1), 0.);
  ASSERT_DOUBLE_EQ(psi_2.GetProbability(0), 0.);
  ASSERT_DOUBLE_EQ(psi_2.GetProbability(1), 1.);
  ASSERT_DOUBLE_EQ(psi_3.GetProbability(0), 1.);
  ASSERT_DOUBLE_EQ(psi_3.GetProbability(1), 1.);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(TwoQubitRegisterTest, InitializeRandomly)
{
  // |psi_0> = |00>
  QubitRegister<ComplexDP> psi_0 (num_qubits_,"base",0);
  // |psi_1> = |10>
  QubitRegister<ComplexDP> psi_1 (num_qubits_,"base",1);
  // random number generator
  std::size_t rng_seed = 7777;
  qhipster::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  psi_0.SetRngPtr(&rnd_generator);
  // Initilize state randomly. Same streams but consecutive numbers.
  psi_0.Initialize("rand",1);
  ASSERT_DOUBLE_EQ( psi_0.ComputeNorm(), 1. );
  //
  psi_1.SetRngPtr(&rnd_generator);
  psi_1.Initialize("rand",1);
  ASSERT_DOUBLE_EQ( psi_1.ComputeNorm(), 1. );

  // The states should be different!
#if 0
  psi_0.Print("state -0-");
  psi_1.Print("state -1-");
  std::cout << "max L2 diff = " << psi_0.MaxL2NormDiff(psi_1)
            << "  ,  max abs diff = " << psi_0.MaxAbsDiff(psi_1)
            << "  ,  overlap = " << psi_0.ComputeOverlap(psi_1)
            << "  ,  |overlap|^2 = " << std::norm(psi_0.ComputeOverlap(psi_1)) << "\n";
#endif
  // This test is not very significant for specific state pairs since the number of
  // qubits is small. In fact, the expectation value of the overlap is 1/2^num_qubits=1/16,
  // but relatively large deviations may be expected.
  // We therefore compute the average of 20 squared overlaps.
  int num_state_pairs = 40;
  double average_squared_overlaps = 0.;
  for (int j=0; j<num_state_pairs; ++j)
  {
      psi_0.Initialize("rand",1);
      psi_1.Initialize("rand",1);
      average_squared_overlaps += std::norm(psi_0.ComputeOverlap(psi_1));
  }
  average_squared_overlaps /= double(num_state_pairs);
  // Considering that the average has std. dev. proportional to 1/sqtr(num_state_pairs),
  // and accounting for 4 std. dev.
  double expected_upper_bound = 1./double(psi_0.GlobalSize())
                                + 1./std::sqrt(double(num_state_pairs)) * 3.;
  EXPECT_LT( average_squared_overlaps, expected_upper_bound );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(TwoQubitRegisterTest, InitializeRandomlyButSame)
{
  // |psi> = |00>
  QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  // random number generator
  std::size_t rng_seed = 7777;
  qhipster::RandomNumberGenerator<double> rng;
  rng.SetSeedStreamPtrs(rng_seed);
  psi.SetRngPtr(&rng);
  //
  int num_states = qhipster::mpi::Environment::GetNumStates();
  assert (num_states==1);
  psi.Initialize("rand",num_states);
  // |psi> = |rand>

  // Initilize the copy: |copy> = |psi>
  QubitRegister<ComplexDP> psi_copy (psi);
  ASSERT_DOUBLE_EQ(psi_copy.MaxAbsDiff(psi), 0 );
  ASSERT_DOUBLE_EQ(psi_copy.MaxL2NormDiff(psi), 0 );
  //
  // Reinitialize |copy> by generating its amplitudes.
  qhipster::RandomNumberGenerator<double> rng_copy;
  rng_copy.SetSeedStreamPtrs(rng_seed);
  psi_copy.SetRngPtr(&rng_copy);
  psi_copy.Initialize("rand",num_states);
  ASSERT_DOUBLE_EQ(psi_copy.MaxAbsDiff(psi), 0 );
  ASSERT_DOUBLE_EQ(psi_copy.MaxL2NormDiff(psi), 0 );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(TwoQubitRegisterTest, Hadamard)
{
  QubitRegister<ComplexDP> psi_0 (num_qubits_,"base",0);
  psi_0.ApplyHadamard(0);
  // |psi_0> = |+0> = |q0=+> x |q1=0>
  ASSERT_NEAR( psi_0.GetProbability(0), 0.5, accepted_error_ );
  psi_0.ApplyPauliZ(0);
  // |psi_0> = |-0> = |q0=-> x |q1=0>
  ASSERT_NEAR( psi_0.GetProbability(0), 0.5, accepted_error_ );
  psi_0.ApplyHadamard(0);
  // |psi_0> = |10> = |q0=1> x |q1=0>
  ASSERT_NEAR( psi_0.GetProbability(0), 1. , accepted_error_ );

  QubitRegister<ComplexDP> psi_1 (num_qubits_,"base",1);
  psi_1.ApplyHadamard(0);
  // |psi_1> = |-0> = |q0=-> x |q1=0>
  ASSERT_NEAR( psi_1.GetProbability(0), 0.5, accepted_error_ );

  QubitRegister<ComplexDP> psi_2 (num_qubits_,"base",2);
  psi_2.ApplyHadamard(1);
  // |psi_2> = |0-> = |q0=0> x |q1=->
  ComplexDP amplitude = ComplexDP(1./std::sqrt(2.), 0. );
  ASSERT_EQ(psi_2.GetGlobalAmplitude(0), amplitude);
  ASSERT_DOUBLE_EQ(psi_2.GetGlobalAmplitude(1).real(), 0.);
  ASSERT_EQ(psi_2.GetGlobalAmplitude(2),-amplitude);
  ASSERT_DOUBLE_EQ(psi_2.GetGlobalAmplitude(3).imag(), 0.);

  QubitRegister<ComplexDP> psi_3 (num_qubits_,"base",3);
  psi_3.ApplyHadamard(0);
  psi_3.ApplyHadamard(1);
  // |psi_3> = |--> = |q0=-> x |q1=->
  amplitude = ComplexDP(1./2., 0. );
  ASSERT_COMPLEX_NEAR(psi_3.GetGlobalAmplitude(0),  amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_3.GetGlobalAmplitude(1), -amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_3.GetGlobalAmplitude(2), -amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_3.GetGlobalAmplitude(3),  amplitude, accepted_error_);
  // Verify that the Hadamard is self-inverse, at least on |psi_3>.
  psi_3.ApplyHadamard(0);
  psi_3.ApplyHadamard(1);
  // |psi_3> = |11>
  ASSERT_NEAR(psi_3.GetProbability(0),1, accepted_error_);
  ASSERT_NEAR(psi_3.GetProbability(1),1, accepted_error_);
  ASSERT_NEAR(psi_3.GetGlobalAmplitude(3).imag(), 0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard TWO_QUBIT_REGISTER_TEST_HPP
