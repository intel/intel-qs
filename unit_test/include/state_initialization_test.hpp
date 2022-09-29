#ifndef STATE_INITIALIZATION_TEST_HPP
#define STATE_INITIALIZATION_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class StateInitializationTest : public ::testing::Test
{
 protected:

  StateInitializationTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 10-qubit state is distributed in more than 512 ranks.
    if (iqs::mpi::Environment::GetStateSize() > 512)
      GTEST_SKIP() << "INFO: small state distributed among too many ranks.";

    G_(0, 0) = {0.592056606032915, 0.459533060553574}; 
    G_(0, 1) = {-0.314948020757856, -0.582328159830658};
    G_(1, 0) = {0.658235557641767, 0.070882241549507}; 
    G_(1, 1) = {0.649564427121402, 0.373855203932477};
  }

  const std::size_t num_qubits_ = 10;
  TM2x2<ComplexDP> G_;
  double accepted_error_ = 1e-15;
  bool do_print_info_ = false;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(StateInitializationTest, ComputationalBasisState)
{
  // Recall that the qubits read from left to right, with the most significant
  // bit on the right (contrary to usual decimal representation).
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4+32+512);
  // |psi> = |1010010001> = |"1+4+32+512">
  ASSERT_DOUBLE_EQ( psi.ComputeNorm(), 1.);
  ASSERT_EQ(psi.GetProbability(0),1);
  ASSERT_EQ(psi.GetProbability(2),1);
  ASSERT_EQ(psi.GetProbability(num_qubits_-2),0);
  ASSERT_EQ(psi.GetProbability(num_qubits_-1),1);
  psi.ApplyPauliX(5);
  psi.ApplyPauliX(6);
  // |psi> = |1010001001> = |'1+4+64+512">
  ASSERT_EQ(psi.GetProbability(5),0);
  ASSERT_EQ(psi.GetProbability(6),1);

  psi.Initialize("base",2+4+64);
  // |psi> = |0110001000> = |2+4+64> = |"70">
  ASSERT_DOUBLE_EQ( psi.ComputeNorm(), 1.);
  ASSERT_EQ(psi.GetProbability(2),1);
  ASSERT_EQ(psi.GetProbability(num_qubits_-4),1);
  psi.ApplyPauliZ(1);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(70).real(), -1.);
  psi.ApplyPauliZ(2);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(70).real(),  1.);
  psi.ApplyPauliZ(3);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(70).real(),  1.);
  psi.ApplyPauliX(4);
  // |psi> = |0001010110> = |"86">
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(86).real(),  1.);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(86).imag(),  0.);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(StateInitializationTest, RandomState)
{
  std::size_t rng_seed = 7777;
  iqs::RandomNumberGenerator<double> rnd_generator_1;
  rnd_generator_1.SetSeedStreamPtrs(rng_seed);
  // The "rand" style cannot be used directly in the creation of the state.
  iqs::QubitRegister<ComplexDP> psi_1(num_qubits_, "base", 0);
  psi_1.SetRngPtr(&rnd_generator_1);
  // |psi_1> = |random>
  psi_1.Initialize("rand", 1);

  ASSERT_DOUBLE_EQ(psi_1.ComputeNorm(), 1);

  // Copy |psi_1> into |psi_2>
  iqs::QubitRegister<ComplexDP> psi_2(psi_1);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Random initialization, using the same RNG but continuing its stream without resetting.
  // |psi_2> = |random>
  psi_2.SetRngPtr(&rnd_generator_1);
  psi_2.Initialize("rand", 1);
#if 0
  std::cout << "max L2 diff = " << psi_2.MaxL2NormDiff(psi_1)
            << "  ,  max abs diff = " << psi_2.MaxAbsDiff(psi_1)
            << "  ,  overlap = " << psi_2.ComputeOverlap(psi_1)
            << "  ,  |overlap|^2 = " << std::norm(psi_2.ComputeOverlap(psi_1)) << "\n";
#endif
  // With very high probability (in the choice of the rng_seed) this test must pass.
  EXPECT_LT( std::norm(psi_2.ComputeOverlap(psi_1)), 0.01 );
  //
  psi_2.Initialize("rand", 1);
  EXPECT_LT( std::norm(psi_2.ComputeOverlap(psi_1)), 0.01 );
  psi_2.Initialize("rand", 1);
  EXPECT_LT( std::norm(psi_2.ComputeOverlap(psi_1)), 0.01 );
  psi_2.Initialize("rand", 1);
  EXPECT_LT( std::norm(psi_2.ComputeOverlap(psi_1)), 0.01 );
  // More precisely, one can check average of overlaps between random states.
  // here only state |psi_2> is reinitialized.
  int num_state_pairs = 40;
  double average_squared_overlaps = 0.;
  for (int j=0; j<num_state_pairs; ++j)
  {
      psi_2.Initialize("rand",1);
      average_squared_overlaps += std::norm(psi_2.ComputeOverlap(psi_1));
  }
  average_squared_overlaps /= double(num_state_pairs);
  // Considering that the average has std. dev. proportional to 1/sqtr(num_state_pairs),
  // and accounting for 4 std. dev.
  double expected_upper_bound = 1./double(psi_1.GlobalSize())
                                + 1./std::sqrt(double(num_state_pairs)) * 3.;
  EXPECT_TRUE( average_squared_overlaps < expected_upper_bound );


  // Create a new RNG initialized to the same seed.
  // if used to initialize |psi_2>, it should recreate |psi_1>.
  iqs::RandomNumberGenerator<double> rnd_generator_2;
  rnd_generator_2.SetSeedStreamPtrs(rng_seed);
  psi_2.SetRngPtr(&rnd_generator_2);
  psi_2.Initialize("rand", 1);
  // For visual inspection, if required:
  if (do_print_info_)
  {
      for (size_t j=0; j<psi_2.GlobalSize(); ++j)
          if (psi_1[j] != psi_2[j])
              std::cout << "@@ psi_1[" << j << "] = " << psi_1[j] << "\n"
                        << "@@ psi_2[" << j << "] = " << psi_2[j] << "\n";
  }

  EXPECT_NEAR(psi_2.MaxAbsDiff(psi_1), 0, accepted_error_ );
  EXPECT_NEAR(psi_2.MaxL2NormDiff(psi_1), 0, accepted_error_ );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(StateInitializationTest, RandomStateSameSeed)
{
  // random number generator
  std::size_t rng_seed = 7315;
  iqs::RandomNumberGenerator<double> rng;
  rng.SetSeedStreamPtrs(rng_seed);
  //
  int num_states = iqs::mpi::Environment::GetNumStates();
  assert (num_states==1);
  //
  // |psi> = |0000>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  psi.SetRngPtr(&rng);
  psi.Initialize("rand",num_states);
  // |psi> = |rand>

  // Initilize the copy: |copy> = |psi>
  iqs::QubitRegister<ComplexDP> psi_copy (psi);
  ASSERT_DOUBLE_EQ(psi_copy.MaxAbsDiff(psi), 0 );
  ASSERT_DOUBLE_EQ(psi_copy.MaxL2NormDiff(psi), 0 );
  //
  iqs::RandomNumberGenerator<double> rng_copy;
  rng_copy.SetSeedStreamPtrs(rng_seed);
  psi_copy.SetRngPtr(&rng_copy);
  psi_copy.Initialize("rand",num_states);

  // For visual inspection, if required:
  if (do_print_info_)
  {
      psi.Print("|random psi> = ");
      psi_copy.Print("|copy psi> = ");
  }

  ASSERT_NEAR(psi_copy.MaxAbsDiff(psi), 0., accepted_error_ );
  ASSERT_NEAR(psi_copy.MaxL2NormDiff(psi), 0., accepted_error_ );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(StateInitializationTest, SuperpositionOfAllComputationalStates)
{
  // |psi> = |++++++++++>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"++++");

  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1);

  ComplexDP amplitude = {1/std::sqrt(psi.GlobalSize()),0};
  ASSERT_COMPLEX_NEAR( psi.GetGlobalAmplitude(71),  amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR( psi.GetGlobalAmplitude(715),  amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR( psi.GetGlobalAmplitude(791),  amplitude, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(StateInitializationTest, ComputationalBasisStateWithQubitPermutation)
{
  std::size_t index = 1+4+32+256+512; // global index using trivial qubit order
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", index);
  // |psi> = |1010010011> = |"1+4+32+256+512">
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index), ComplexDP(1, 0), accepted_error_);

  // Valid only with trivial permutation
  std::size_t myrank = iqs::mpi::Environment::GetStateRank();
  std::size_t whereid = index /psi.LocalSize();
  if (whereid == myrank)
  {
      std::size_t lclind = (index % psi.LocalSize());
      ASSERT_COMPLEX_NEAR(psi[lclind], ComplexDP(1, 0), accepted_error_);
  }

  // Define and apply qubit permutation
  //                                                 |--------|
  std::vector<std::size_t>  map = {0, 1, 2, 3, 4, 5, 9, 7, 8, 6};
  psi.PermuteQubits(map, "direct");
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index), ComplexDP(1, 0), accepted_error_);
  std::size_t index_after_permutation = psi.qubit_permutation->program2data_(index);
  whereid = index_after_permutation/psi.LocalSize();
  if (whereid == myrank)
  {
      std::size_t lclind = (index_after_permutation % psi.LocalSize());
      ASSERT_COMPLEX_NEAR(psi[lclind], ComplexDP(1, 0), accepted_error_);
  }

  // Qubit permutation as above, reinitialize the state
  index = 1+64; // 2^0+2^6
  psi.Initialize ("base", index);
  //
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index), ComplexDP(1, 0), accepted_error_);
  index_after_permutation = psi.qubit_permutation->program2data_(index);
  whereid = index_after_permutation/psi.LocalSize();
  if (whereid == myrank)
  {
      std::size_t lclind = (index_after_permutation % psi.LocalSize());
      ASSERT_COMPLEX_NEAR(psi[lclind], ComplexDP(1, 0), accepted_error_);
  }
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(StateInitializationTest, DeathTest)
{
  // Skip death-tests if compiler flag NDEBUG is defined.
#ifdef NDEBUG
  GTEST_SKIP() << "INFO: test is skipped when compiler flag NDEBUG is not defined.";
#endif

  // Skip death-tests if MPI size > 1.
  if (iqs::mpi::Environment::GetStateSize() > 1)
      GTEST_SKIP();

  // To switch off the warning message about DEATH test not being thread safe.
  ::testing::FLAGS_gtest_death_test_style = "threadsafe"; 

  // |psi> = |0000000000> = |"0">
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  // Index outside the global range.
  std::size_t index;
  index = UL(1L << UL(num_qubits_));
  ASSERT_DEATH( psi.Initialize("base", index), "");
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard STATE_INITIALIZATION_TEST_HPP
