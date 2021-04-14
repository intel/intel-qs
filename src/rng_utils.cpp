#include <iostream>

#include "../include/rng_utils.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////
// Certain methods are the same both with and without MKL.
/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
RandomNumberGenerator<Type>::RandomNumberGenerator(RandomNumberGenerator * source_rng)
{
  this->SetSeedStreamPtrs( source_rng->GetSeed() );
  this->SkipAhead( source_rng->GetNumGeneratedOrSkippedPoolNumbers() , "pool"  );
  this->SkipAhead( source_rng->GetNumGeneratedOrSkippedStateNumbers(), "state" );
  this->SkipAhead( source_rng->GetNumGeneratedOrSkippedLocalNumbers(), "local" );
}

/////////////////////////////////////////////////////////////////////////////////////////
// The next two methods will be specialized for 'float' and 'double' below.
// If Type is another typename, then report the error.

template <typename Type>
void RandomNumberGenerator<Type>::UniformRandomNumbers
( Type *value, std::size_t size , Type a, Type b, std::string shared )
{
  std::cout << " ~~~~~~~~~~~~~~ wrong type for the random number generator! ~~~~~~~~~ \n";
  assert(0);
}

template <typename Type>
void RandomNumberGenerator<Type>::GaussianRandomNumbers
( Type *value, std::size_t size , std::string shared )
{
  std::cout << " ~~~~~~~~~~~~~~ wrong type for the random number generator! ~~~~~~~~~ \n";
  assert(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

#ifndef USE_MKL

/////////////////////////////////////////////////////////////////////////////////////////
// Without MKL we only define a few methods.
/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
std::mt19937 * RandomNumberGenerator<Type>::SelectGenerator (std::string shared)
{
  if (shared=="local") {
      return &_local_generator;
  } else if (shared=="state") {
      return &_state_generator;
  } else if (shared=="pool") {
      return &_pool_generator;
  } else
      assert(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
std::mt19937 * RandomNumberGenerator<Type>::SelectGeneratorAndUpdateCounter
( std::size_t size , std::string shared )
{
  if (shared=="local") {
      _num_generated_or_skipped_local_numbers += size;
      return &_local_generator;
  } else if (shared=="state") {
      _num_generated_or_skipped_state_numbers += size;
      return &_state_generator;
  } else if (shared=="pool") {
      _num_generated_or_skipped_pool_numbers  += size;
      return &_pool_generator;
  } else
      assert(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
void RandomNumberGenerator<Type>::SetSeedStreamPtrs ( std::size_t RNG_seed )
{
  _seed = RNG_seed;
  _num_generated_or_skipped_local_numbers=0;
  _num_generated_or_skipped_state_numbers=0;
  _num_generated_or_skipped_pool_numbers =0;
  // Each rank in the pool uses a different seed.
  int num_states = mpi::Environment::GetNumStates();
  int state_id   = mpi::Environment::GetStateId();
  int pool_rank  = mpi::Environment::GetPoolRank();
  // TODO: verify that the seed actually gives independent streams.
  _pool_generator.seed(RNG_seed + 0);
  _state_generator.seed(RNG_seed + 1 + state_id);
  _local_generator.seed(RNG_seed + 1 + num_states + pool_rank);
}

template <typename Type>
void RandomNumberGenerator<Type>::SkipAhead ( std::size_t num_skip, std::string shared )
{
  // The use of .discard() depends on the actual implementation of the distribution.
  // Sometimes more than one use of the std::mt19937 generator are needed, even to
  // generate a single uniformly distributed (double) number.
  // To avoid problems, one generate and discards random numbers.

#if 1
  Type values [num_skip];
  this->UniformRandomNumbers(values, num_skip, 0., 1., shared);
  // Recall that the number of skips is updated inside UniformRandomNumbers().
  assert(shared=="local" || shared=="state" || shared=="pool");
#else
  if (shared=="local") {
      _num_generated_or_skipped_local_numbers += num_skip;
      _local_generator.discard(num_skip);
  } else if (shared=="state") {
      _num_generated_or_skipped_state_numbers += num_skip;
      _state_generator.discard(num_skip);
  } else if (shared=="pool") {
      _num_generated_or_skipped_pool_numbers  += num_skip;
      _pool_generator.discard(num_skip);
  } else
      assert(0);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////
// Generate random numbers in [0,1).

template <>
void RandomNumberGenerator<float>::UniformRandomNumbers
( float *value, std::size_t size , float a, float b, std::string shared )
{
  std::mt19937 * generator;
  generator = this->SelectGeneratorAndUpdateCounter(size, shared);

  float * temp_ptr = value;
  // If [a,b) is not [0,1), translate and scale the random numbers.
  for (std::size_t i=0; i<size; ++i)
  {
      *temp_ptr = a + (b-a) * u_distribution(*generator);
      ++temp_ptr;
  }
  return;
}

template <>
void RandomNumberGenerator<double>::UniformRandomNumbers
( double *value, std::size_t size , double a, double b, std::string shared )
{
  std::mt19937 * generator;
  generator = this->SelectGeneratorAndUpdateCounter(size, shared);

  double * temp_ptr = value;
  // If [a,b) is not [0,1), translate and scale the random numbers.
  for (std::size_t i=0; i<size; ++i)
  {
      *temp_ptr = a + (b-a) * u_distribution(*generator);
      ++temp_ptr;
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Generate a random integer in [a,b).

template <typename Type>
void RandomNumberGenerator<Type>::RandomIntegersInRange
( int *value, std::size_t size, int a, int b, std::string shared )
{
  std::mt19937 * generator;
  generator = this->SelectGeneratorAndUpdateCounter(size, shared);

  Type rand;
  int * temp_ptr = value;
  for (std::size_t i=0; i<size; ++i)
  {
      rand = u_distribution(*generator);
      rand = (Type)a + rand * Type(b - a);
      // Note that there are problems with double-->int conversion:
      // - the range of double (output of floor) is larger;
      // - but not all integers are exactly representable within the double range.
      *temp_ptr = (int)std::floor(rand);
      ++temp_ptr;
  }

  return;
}
    
/////////////////////////////////////////////////////////////////////////////////////////
// Generate normal-distributed numbers (mean=0, std.dev.=1).

template <>
void RandomNumberGenerator<float>::GaussianRandomNumbers
( float *value, std::size_t size , std::string shared )
{
  std::mt19937 * generator;
  // Assuming that to generate a Gaussian number one needs 2 uniformly distributed numbers.
  generator = this->SelectGeneratorAndUpdateCounter(2*size, shared);

  float * temp_ptr = value;
  for (std::size_t i=0; i<size; ++i)
  {
       *temp_ptr = n_distribution(*generator);
       ++temp_ptr;
  }
  return;
}

template <>
void RandomNumberGenerator<double>::GaussianRandomNumbers
( double *value, std::size_t size , std::string shared )
{
  std::mt19937 * generator;
  // Assuming that to generate a Gaussian number one needs 2 uniformly distributed numbers.
  generator = this->SelectGeneratorAndUpdateCounter(2*size, shared);

  double * temp_ptr = value;
  for (std::size_t i=0; i<size; ++i)
  {
       *temp_ptr = n_distribution(*generator);
       ++temp_ptr;
  }
  return;
}

#else

/////////////////////////////////////////////////////////////////////////////////////////
// With MPI and MKL we define multiple methods (from here to the end of the class).
/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
VSLStreamStatePtr RandomNumberGenerator<Type>::SelectGenerator (std::string shared)
{
  if (shared=="local") {
      return _local_stream;
  } else if (shared=="state") {
      return _state_stream;
  } else if (shared=="pool") {
      return _pool_stream;
  } else
      assert(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
VSLStreamStatePtr RandomNumberGenerator<Type>::SelectGeneratorAndUpdateCounter
( std::size_t size , std::string shared )
{
  if (shared=="local") {
      _num_generated_or_skipped_local_numbers += size;
      return _local_stream;
  } else if (shared=="state") {
      _num_generated_or_skipped_state_numbers += size;
      return _state_stream;
  } else if (shared=="pool") {
      _num_generated_or_skipped_pool_numbers  += size;
      return _pool_stream;
  } else
      assert(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

// Set the pointers to the VSL streams of random numbers.
// TODO FIXME: only present for backward compatibility.
template <typename Type>
void RandomNumberGenerator<Type>::SetRndStreamPtrs
( VSLStreamStatePtr pool , VSLStreamStatePtr state , VSLStreamStatePtr local )
{
  _pool_stream  = pool;
  _state_stream = state;
  _local_stream = local;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Set the seed and type of pointers
template <typename Type>
void RandomNumberGenerator<Type>::SetSeedStreamPtrs( std::size_t RNG_seed )
{
  _seed = RNG_seed;
  int pool_rank  = mpi::Environment::GetPoolRank();
  int num_states = mpi::Environment::GetNumStates();
  int state_id   = mpi::Environment::GetStateId();

  // VSL provides up to 6024 Mersenne Twister pseudorandom number generators.
  assert(mpi::Environment::GetPoolSize() + num_states + 1 <=6024);

  // Pool-wide stream.
  int errcode = vslNewStream(&_pool_stream, VSL_BRNG_MT2203+0, RNG_seed);
  assert(errcode == VSL_STATUS_OK);
  // State-wide stream.
  errcode = vslNewStream(&_state_stream, VSL_BRNG_MT2203+1+state_id, RNG_seed);
  assert(errcode == VSL_STATUS_OK);
  // Local stream.
  errcode = vslNewStream(&_local_stream, VSL_BRNG_MT2203+1+num_states+pool_rank, RNG_seed);
  assert(errcode == VSL_STATUS_OK);
}

/////////////////////////////////////////////////////////////////////////////////////////

template <typename Type>
void RandomNumberGenerator<Type>::SkipAhead ( std::size_t num_skip, std::string shared )
{
  VSLStreamStatePtr stream;
  stream = this->SelectGenerator(shared);
  int errcode = vslSkipAheadStream(stream, num_skip);
  // If the BRNG does not support skip-ahead, generate the samples and discard them.
  if (errcode==VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED)
  {
      // It is unclear why, but when num_skip is > 2^18 segmentation errors appear.
      if (num_skip < (1 << 18) )
      {
          Type values [num_skip];
          this->UniformRandomNumbers(values, num_skip, 0., 1., shared);
      }
      else
      {
          Type value;
          for (size_t j=0; j<num_skip; ++j)
              this->UniformRandomNumbers(&value, 1, 0., 1., shared);
      }
  }
  else
      assert(errcode == VSL_STATUS_OK);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Generate random numbers in [0,1).

template <>
void RandomNumberGenerator<float>::UniformRandomNumbers
( float *value, std::size_t size , float a, float b, std::string shared )
{
  VSLStreamStatePtr stream;
  stream = this->SelectGeneratorAndUpdateCounter(size, shared);
  int errcode = vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, size, value, a, b);
  assert(errcode == VSL_STATUS_OK);
}

template <>
void RandomNumberGenerator<double>::UniformRandomNumbers
( double *value, std::size_t size , double a, double b, std::string shared )
{
  VSLStreamStatePtr stream;
  stream = this->SelectGeneratorAndUpdateCounter(size, shared);
  int errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, size, value, a, b);
  assert(errcode == VSL_STATUS_OK);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Generate a single random integer in [a,b).

template <typename Type>
void RandomNumberGenerator<Type>::RandomIntegersInRange
( int *value, std::size_t size, int a, int b, std::string shared )
{
  VSLStreamStatePtr stream;
  stream = this->SelectGeneratorAndUpdateCounter(size, shared);
  int errcode = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, size, value, a, b);
  assert(errcode == VSL_STATUS_OK);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Generate Gaussian numbers (std.dev.=1, mean=0).

template <>
void RandomNumberGenerator<float>::GaussianRandomNumbers
( float *value, std::size_t size , std::string shared )
{
  VSLStreamStatePtr stream;
  // Assuming that to generate a Gaussian number one needs 2 uniformly distributed numbers.
  stream = this->SelectGeneratorAndUpdateCounter(2*size, shared);
  int errcode = vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, size, value, 0.,1.);
  assert(errcode == VSL_STATUS_OK);
}

template <>
void RandomNumberGenerator<double>::GaussianRandomNumbers
( double *value, std::size_t size , std::string shared )
{
  VSLStreamStatePtr stream;
  // Assuming that to generate a Gaussian number one needs 2 uniformly distributed numbers.
  stream = this->SelectGeneratorAndUpdateCounter(2*size, shared);
  int errcode = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, size, value, 0.,1.);
  assert(errcode == VSL_STATUS_OK);
}

#endif

/////////////////////////////////////////////////////////////////////////////////////////

template class RandomNumberGenerator<float>;
template class RandomNumberGenerator<double>;

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type, typename TypeFloat>
void ShuffleFisherYates( std::vector<Type> &array,
                         RandomNumberGenerator<TypeFloat> *rnd_generator_ptr,
                         std::string shared )
{
  int ra,rb ;
  Type temp;
  int dim = array.size();

  for ( ra = dim-1; ra>0; --ra)
  {
      rnd_generator_ptr->RandomIntegersInRange (&rb, 1UL, 0, ra+1, shared);
      // exchange the components
      if (ra!=rb)
      {
          temp=array[ra];
          array[ra]=array[rb];
          array[rb]=temp;
      }
  }
}

template void ShuffleFisherYates<int,double>
  (std::vector<int>     &, RandomNumberGenerator<double> *, std::string);
template void ShuffleFisherYates<float,double>
  (std::vector<float>   &, RandomNumberGenerator<double> *, std::string);
template void ShuffleFisherYates<double,double>
  (std::vector<double>  &, RandomNumberGenerator<double> *, std::string);
template void ShuffleFisherYates<unsigned,double>
  (std::vector<unsigned>&, RandomNumberGenerator<double> *, std::string);

/////////////////////////////////////////////////////////////////////////////////////////

}	// close namespace iqs
