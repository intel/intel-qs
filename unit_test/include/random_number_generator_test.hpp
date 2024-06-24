#ifndef RANDOM_NUMBER_GENERATOR_TEST_HPP
#define RANDOM_NUMBER_GENERATOR_TEST_HPP

#include <iomanip>
#include <map>

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class RandomNumberGeneratorTest : public ::testing::Test
{
 protected:

  RandomNumberGeneratorTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
        GTEST_SKIP();
  }

  // just before the 'destructor'
//  void TearDown() override {}

  double accepted_error_ = 1e-15;
};
	
//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(RandomNumberGeneratorTest, BasicUse)
{
  std::size_t num_samples = 10;
  double random [num_samples];
  iqs::RandomNumberGenerator<double> rng;
  std::size_t seed = 777;
  rng.SetSeedStreamPtrs(seed);

  // Generate the samples.
  rng.UniformRandomNumbers(random, num_samples, 0, 1, "local");
  // They should not be the same number.
  ASSERT_TRUE(random[0] != random[num_samples-1]);
  ASSERT_EQ(rng.GetNumGeneratedOrSkippedLocalNumbers(), num_samples);
  ASSERT_EQ(rng.GetNumGeneratedOrSkippedStateNumbers(), 0);
  ASSERT_EQ(rng.GetNumGeneratedOrSkippedPoolNumbers() , 0 );

  // New stream, it should be a copy of the previous one since seed is the same.
  iqs::RandomNumberGenerator<double> rng_copy;
  rng_copy.SetSeedStreamPtrs(seed);
  double random_copy [num_samples];
  // Generate the samples.
  rng_copy.UniformRandomNumbers(random_copy, num_samples, 0., 1.);
  for (std::size_t j=0; j<num_samples; ++j)
      ASSERT_EQ(random[j], random_copy[j]);
  ASSERT_EQ(rng_copy.GetNumGeneratedOrSkippedLocalNumbers(), num_samples);
  ASSERT_EQ(rng_copy.GetNumGeneratedOrSkippedStateNumbers(), 0);
  ASSERT_EQ(rng_copy.GetNumGeneratedOrSkippedPoolNumbers() , 0 );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(RandomNumberGeneratorTest, Visualization)
{
  if (iqs::mpi::Environment::GetStateRank() != 0)
      GTEST_SKIP();

  std::size_t num_samples = 20*2000;
  double random [num_samples];
  iqs::RandomNumberGenerator<double> rng;
  std::size_t seed = 777;
  rng.SetSeedStreamPtrs(seed);

  // Generate the samples.
  rng.UniformRandomNumbers(random, num_samples, 0., 20., "local");
  std::map<int, int> histogram;
  for (int n = 0; n < num_samples; ++n)
  {
      ++histogram[int(random[n])];
  }
  std::cout << "Histogram of 20*2000 numbers drawn uniformly in [0,20[ "
            << "(each * corresponds to 100 draws)\n";
  for (auto p : histogram)
  {
      std::cout << std::setw(2) << p.first << ' ' << std::string(p.second / 100, '*') << '\n';
  }
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(RandomNumberGeneratorTest, DifferentSeed)
{
  std::size_t num_samples = 100;
  double random [num_samples];
  iqs::RandomNumberGenerator<double> rng;
  std::size_t seed = 717;
  rng.SetSeedStreamPtrs(seed);
  // Generate the samples.
  rng.UniformRandomNumbers(random, num_samples, 0., 1.,"local");
  ASSERT_TRUE(random[0] != random[num_samples-1]);

  // New stream, it should be different.
  iqs::RandomNumberGenerator<double> rng_2;
  rng_2.SetSeedStreamPtrs(seed+1);
  double random_2 [num_samples];
  // Generate the samples.
  rng_2.UniformRandomNumbers(random_2, num_samples, 0., 1.);
  // Verify that the two set of numbers have null intersection.
  for (std::size_t j=0; j<num_samples; ++j)
      for (std::size_t k=0; k<num_samples; ++k)
          ASSERT_NE(random[j], random_2[k]);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(RandomNumberGeneratorTest, LargeRandomVectors)
{
  // A large vector size gave errors with GCC when running on github servers.
#if defined(__ICC) || defined(__INTEL_COMPILER)
  std::size_t num_samples = (1UL << 22);
#else
  std::size_t num_samples = (1UL << 16);
#endif
  double random [num_samples];
  iqs::RandomNumberGenerator<double> rng;
  std::size_t seed = 717;
  rng.SetSeedStreamPtrs(seed);
  // Generate the samples.
  rng.UniformRandomNumbers(random, num_samples, 0., 1., "local");
  ASSERT_TRUE(random[0] != random[num_samples-1]);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(RandomNumberGeneratorTest, SkipMethod)
{
  std::size_t half_samples = 40;
  std::size_t num_samples = 2*half_samples;
  double random [num_samples];
  iqs::RandomNumberGenerator<double> rng;
  std::size_t seed = 717;
  rng.SetSeedStreamPtrs(seed);
  // Generate the samples.
  rng.UniformRandomNumbers(random, num_samples, 0., 1., "local");

  // New stream, it should corresponds to the second half of the previous stream.
  iqs::RandomNumberGenerator<double> rng_2;
  double random_2 [half_samples];
  rng_2.SetSeedStreamPtrs(seed);
  rng_2.SkipAhead(half_samples,"local");
  // Generate the samples.
  rng_2.UniformRandomNumbers(random_2, half_samples, 0., 1., "local");
  for (std::size_t j=0; j<half_samples; ++j)
      ASSERT_EQ(random[j+half_samples], random_2[j]);

  ASSERT_EQ(  rng.GetNumGeneratedOrSkippedLocalNumbers(),
            rng_2.GetNumGeneratedOrSkippedLocalNumbers() );

  // Same test, but generating Gaussian numbers.
  rng.GaussianRandomNumbers(random, num_samples, "local");
#if 0
  // Preliminary tests seems to indicate that, without MKL, one needs 2 uses of
  // the generator per uniform double and 3 usus for Gaussian double.
  // A systematic study is not required at this point.
  rng_2.SkipAhead(half_samples,"local");
#else
  rng_2.GaussianRandomNumbers(random_2, half_samples, "local");
#endif
  rng_2.GaussianRandomNumbers(random_2, half_samples, "local");
  for (std::size_t j=0; j<half_samples; ++j)
      ASSERT_EQ(random[j+half_samples], random_2[j]);
}

//////////////////////////////////////////////////////////////////////////////

#ifdef USE_MKL 

// This test is special: it is to verify that we are counting the number of generated
// random values correctly. For example, to generate a single Gaussian number, most
// methods need to consume two uniformly distributed numbers.

TEST_F(RandomNumberGeneratorTest, VslGenerator)
{
  std::size_t quarter_samples = 8;
  std::size_t num_samples = 4*quarter_samples;
  double random [num_samples];
  std::size_t seed = 717;
  int errcode;

  VSLStreamStatePtr stream;
  errcode = vslNewStream(&stream, VSL_BRNG_MT2203+0, seed);
  assert(errcode == VSL_STATUS_OK);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, num_samples, random, 0, 1);
  assert(errcode == VSL_STATUS_OK);

  // New stream with same seed that will be used to verify alternative to skip.
  double random_2 [quarter_samples];
  VSLStreamStatePtr stream_2;
  errcode = vslNewStream(&stream_2, VSL_BRNG_MT2203+0, seed);
  assert(errcode == VSL_STATUS_OK);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  for (std::size_t j=0; j<quarter_samples; ++j)
      ASSERT_EQ(random[j], random_2[j]);
  // Skip one third of the samples by generating uniform random numbers.
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  // Now generate the numbers again: they should corresponds to the third quarter of random.
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  for (std::size_t j=0; j<quarter_samples; ++j)
      ASSERT_EQ(random[j+2*quarter_samples], random_2[j]);

  // Re-initialize stream_2 with same seed.
  errcode = vslNewStream(&stream_2, VSL_BRNG_MT2203+0, seed);
  assert(errcode == VSL_STATUS_OK);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  for (std::size_t j=0; j<quarter_samples; ++j)
      ASSERT_EQ(random[j], random_2[j]);
  // Skip two quarters of the samples by generating Gaussian random numbers:
  // each Gaussian number needs two uniformly distributed numbers.
  errcode = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream_2, quarter_samples,
                          random_2, 0.,1.);
  assert(errcode == VSL_STATUS_OK);
  // Now generate the numbers again: they should corresponds to the last quarter of random.
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  for (std::size_t j=0; j<quarter_samples; ++j)
      ASSERT_EQ(random[j+3*quarter_samples], random_2[j]);

  // Re-initialize stream_2 with same seed.
  errcode = vslNewStream(&stream_2, VSL_BRNG_MT2203+0, seed);
  assert(errcode == VSL_STATUS_OK);
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  for (std::size_t j=0; j<quarter_samples; ++j)
      ASSERT_EQ(random[j], random_2[j]);
  // Skip one quarter of the samples by generating random integers in [1,5).
  int value;
  for (std::size_t j=0; j<quarter_samples; ++j)
  {
      errcode = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, 1L, &value, 1, 5);
      assert(errcode == VSL_STATUS_OK);
  }
  // Now generate the numbers again: they should corresponds to the third quarter of random.
  errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_2, quarter_samples,
                         random_2, 0., 1.);
  assert(errcode == VSL_STATUS_OK);
  for (std::size_t j=0; j<quarter_samples; ++j)
      ASSERT_EQ(random[j+2*quarter_samples], random_2[j]);
}
#endif

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard RANDOM_NUMBER_GENERATOR_TEST_HPP
