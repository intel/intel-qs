#ifndef APPLY_SWAP_GATE_TEST_HPP
#define APPLY_SWAP_GATE_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class ApplySwapGateTest : public ::testing::Test
{
 protected:

  ApplySwapGateTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
        GTEST_SKIP();

    // All tests are skipped if the 4-qubit state is distributed in more than 2^3 ranks.
    // In fact the MPI version needs to allocate half-the-local-storage for communication.
    // If the local storage is a single amplitude, this cannot be further divided.
    if (qhipster::mpi::Environment::GetStateSize() > 512)
        GTEST_SKIP();
  }

  const std::size_t num_qubits_ = 10;
  double accepted_error_ = 1e-15;
  double sqrt2_ = std::sqrt(2.);
};

//////////////////////////////////////////////////////////////////////////////

// Emulation of SWAP by update the permutation.
TEST_F(ApplySwapGateTest, Emulation)
{
  ComplexDP amplitude_0 = { 0. , 0. };
  ComplexDP amplitude_1 = { 1. , 0. };

  // Initial state |0010000000>
  std::size_t index = 4;
  QubitRegister<ComplexDP> psi (num_qubits_, "base", index);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index), amplitude_1, accepted_error_);

  ASSERT_EQ( (*psi.permutation)[1], 1);
  ASSERT_EQ( (*psi.permutation)[2], 2);
  psi.EmulateSwap(1, 2);
  ASSERT_EQ( (*psi.permutation)[1], 2);
  ASSERT_EQ( (*psi.permutation)[2], 1);
  // Initial state |0100000000>
  std::size_t index_after = 2;
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index      ), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index_after), amplitude_1, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard APPLY_SWAP_GATE_TEST_HPP
