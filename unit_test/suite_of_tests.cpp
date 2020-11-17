//////////////////////////////////////////////////////////////////////////////
/// @file suite_of_tests.cpp
/// @brief Collection of unit tests based on the 'googletest' framework.
/// 
/// It is important to keep in mind that the framework was not developed for MPI
/// applications. For this reason we have to skip a few tests under certain
/// circumstances:
/// - if the state is divided into too many rank so that a single amplitude is
///   stored per rank ==> skip test
/// - death tests are problematic when the MPI size is greater than 1 == > skip test
/// - tests related to having multiple states in the pool typically have
///   constraints like lower bound of MPI size. If not fulfilled ==> skip test
///
/// Link to the 'googletest' framework:
/// https://github.com/google/googletest
//////////////////////////////////////////////////////////////////////////////

#include <cmath>	// std::abs(), std::cos(), std::pow(), ...
#include <iostream>

// googletest
#include "gtest/gtest.h"

// Class that we want to test: QubitRegister<ComplexDP>
#include "../include/qureg.hpp"

// To compare complex values in an approximated way.
#define ASSERT_COMPLEX_NEAR(val1,val2,error)   \
ASSERT_NEAR(val1.real(),val2.real(),error);    \
ASSERT_NEAR(val1.imag(),val2.imag(),error);

// Headers with the implementation of the various tests.

#include "include/compiler_flags_test.hpp"

// Outside class QubitRegister.
#include "include/conversion_test.hpp"
#include "include/tinymatrix_test.hpp"
#include "include/random_number_generator_test.hpp"
#include "include/gate_counter_test.hpp"
#include "include/permutation_test.hpp"

// Small registers.
#include "include/one_qubit_register_test.hpp"
#include "include/two_qubit_register_test.hpp"

// State initialization and one-qubit gates.
#include "include/apply_1q_gate_test.hpp"
#include "include/single_qubit_gates_test.hpp"
#include "include/apply_swap_gate_test.hpp"

// State Initialization and Measurement (SPAM) operations.
#include "include/state_initialization_test.hpp"
#include "include/measure_test.hpp"
#include "include/expectation_values_test.hpp"

// Utility methods and distributed implementation.
#include "include/utility_methods_test.hpp"
#include "include/chunking_communication_test.hpp"

// Noisy simulations and extra features.
#include "include/qureg_permute_test.hpp"
#include "include/noisy_simulation_test.hpp"
#include "include/qaoa_features_test.hpp"

// Pool functionality of MPI environment.
#include "include/multiple_states_test.hpp"

//////////////////////////////////////////////////////////////////////////////
// main
// We want to unit_test an MPI program.

int main(int argc, char **argv)
{
  // Initialize the MPI environmentk, if MPI exists.
  qhipster::mpi::Environment env(argc, argv);
  int my_rank_id = 0;
#ifdef INTELQS_HAS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_id);
#endif 
 
  // When MPI exists, release all process associated to dummy ranks.
//  if (env.IsUsefulRank()==false)
//    return 0;

  // When MPI exists, we want the report only from the main rank.
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
  if (my_rank_id != 0)
      delete listeners.Release(listeners.default_result_printer());
  else
      std::cout << "\n\n";

  // Usual gtest main but without returning the macro directly.
  int result;
  ::testing::InitGoogleTest(&argc, argv);
  result = RUN_ALL_TESTS();

  if (my_rank_id==0)
      std::cout << "\n\n";

  // Finalize the MPI environment.
//  qhipster::mpi::PoolBarrier();
//  env.~Environment();
  return result;
}

//////////////////////////////////////////////////////////////////////////////
