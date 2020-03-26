#ifndef RNG_UTILS_HPP
#define RNG_UTILS_HPP

#include <cassert>
#include <string>
#include <vector>

#ifdef USE_MKL 
#include "mkl.h"
#include "mkl_vsl.h"
#else
#include <cmath>	// std::floor
#include <random>
#endif

#include "mpi_env.hpp"

/// @file rng_utils.hpp
/// @brief Define the @c RandomNumberGenerator class and its methods.

/////////////////////////////////////////////////////////////////////////////////////////
// utility methods related to the generation of (pseudo-)random numbers.
/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

/////////////////////////////////////////////////////////////////////////////////////////

/// @class RandomNumberGenerator
/// @brief Used to generate random numbers that are local to each rank, or common to
/// the state or the complete pool.
///
/// The generation of numbers and the method to skip ahead are more efficient
/// when MKL (and in particular VSL) is used.

// Type = float, double
template <typename Type>
class  RandomNumberGenerator
{
  private:

    std::size_t _seed=0;
    std::size_t _num_generated_or_skipped_local_numbers=0;
    std::size_t _num_generated_or_skipped_state_numbers=0;
    std::size_t _num_generated_or_skipped_pool_numbers =0;

#ifdef USE_MKL
    // vsl streams for random numbers
    VSLStreamStatePtr  _pool_stream;
    VSLStreamStatePtr _state_stream;
    VSLStreamStatePtr _local_stream;
#else
//   std::default_random_engine generator;
   std::mt19937  _pool_generator;
   std::mt19937 _state_generator;
   std::mt19937 _local_generator;
   std::uniform_real_distribution<Type> u_distribution
       = std::uniform_real_distribution<Type>(0.0,1.0);	// uniform in [0,1[
   std::normal_distribution<Type> n_distribution
       = std::normal_distribution<Type>(0.0,1.0);	// mean=0, std.dev.=1
#endif

  public:

    RandomNumberGenerator() {
#ifdef USE_MKL
        _pool_stream  = nullptr;
        _state_stream = nullptr;
        _local_stream = nullptr;
#endif
	}

    ~RandomNumberGenerator() {}

    /// Initialize RNG by copying the streams of the source RNG.
    RandomNumberGenerator(RandomNumberGenerator * source_rng);

    /// Get basic quantities.
    std::size_t GetSeed ( )
    { return _seed; }

    std::size_t GetNumGeneratedOrSkippedLocalNumbers ( )
    { return _num_generated_or_skipped_local_numbers; }
 
    std::size_t GetNumGeneratedOrSkippedStateNumbers ( )
    { return _num_generated_or_skipped_state_numbers; }

    std::size_t GetNumGeneratedOrSkippedPoolNumbers ( )
    { return _num_generated_or_skipped_pool_numbers; }

    /// Set one different seed for each MPI rank (no MKL) or assign different streams (VSL).
    void SetSeedStreamPtrs ( std::size_t RNG_seed );

    /// Skip ahead.
    void SkipAhead ( std::size_t num_skip, std::string shared = "local" );

    /// Generate random numbers in [a,b):
    /// - size indicates how many numbers
    /// - shared can be: local, state, pool
    void UniformRandomNumbers( Type *value, std::size_t size = 1UL, Type a=0., Type b=1., 
                               std::string shared = "local");

    /// Generate random gaussian numbers (mean value = 0, std.dev = 1):
    /// - size indicates how many numbers
    /// - shared can be: local, state, pool
    void GaussianRandomNumbers( Type *value, std::size_t size = 1UL,
                                std::string shared = "local");

    /// Generate random integers in [a,b), default being {0,1}..
    /// - size indicates how many numbers
    /// - shared can be: local, state, pool
    void RandomIntegersInRange( int *value, std::size_t size = 1UL, int a=0, int b=2,
                               std::string shared = "local");

#ifdef USE_MKL
    /// Set the pointers to the VSL streams of random numbers.
    void SetRndStreamPtrs( VSLStreamStatePtr pool , VSLStreamStatePtr state ,
                           VSLStreamStatePtr local );

    /// Get pointer to the corresponding VSL stream.
    VSLStreamStatePtr GetPoolStreamPtr  () const {return  _pool_stream; }
    VSLStreamStatePtr GetStateStreamPtr () const {return _state_stream; }
    VSLStreamStatePtr GetLocalStreamPtr () const {return _local_stream; }
#endif

  private:

    /// Utility function.
#ifdef USE_MKL
    VSLStreamStatePtr SelectGeneratorAndUpdateCounter( std::size_t size, std::string shared);
    VSLStreamStatePtr SelectGenerator(std::string shared);
#else
    std::mt19937 * SelectGeneratorAndUpdateCounter( std::size_t size, std::string shared );
    std::mt19937 * SelectGenerator(std::string shared);
#endif

/////////////////////////////////////////////////////////////////////////////////////////

};	// end of class RandomNumberGenerator

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// Shuffle array using the Fisher and Yates shuffle algorithm.
template<typename Type, typename TypeFloat>
void ShuffleFisherYates( std::vector<Type> &array,
                         RandomNumberGenerator<TypeFloat> *rnd_generator_ptr,
                         std::string shared="local");

/////////////////////////////////////////////////////////////////////////////////////////

}	// close namespace qhipster

#endif	// close header guards RNG_UTILS_HPP
