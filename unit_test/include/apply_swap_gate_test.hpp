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

// Emulation of SWAP by update the permutation.
TEST_F(ApplySwapGateTest, ComparisonWithThreeCnots)
{
  std::size_t rng_seed = 7777;
  qhipster::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  // The "rand" style cannot be used directly in the creation of the state.
  QubitRegister<ComplexDP> psi_1(num_qubits_, "base", 0);
  psi_1.SetRngPtr(&rnd_generator);
  // |psi_1> = |random>
  psi_1.Initialize("rand", 1);

  // |psi_2> = |psi_1> = |random>
  QubitRegister<ComplexDP> psi_2(psi_1);

  // Compare two implementations of the SWAP gate.
  unsigned qubit1 = 2, qubit2 = num_qubits_-1;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
  qubit1 = 4;
  qubit2 = num_qubits_-2;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // TODO: add the SWAP between (probably) distributed qubits, like n-1 and n-2.
  // Compare two implementations of the SWAP gate.
  qubit1 = num_qubits_-1;
  qubit2 = num_qubits_-2;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );
}

//////////////////////////////////////////////////////////////////////////////

// Extra unit test for distributed implementation of the SWAP gate.
TEST_F(ApplySwapGateTest, TestForDistributedImplementation)
{
  unsigned num_qubits = 20;
  std::size_t rng_seed = 12345;
  qhipster::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  // The "rand" style cannot be used directly in the creation of the state.
  QubitRegister<ComplexDP> psi_1(num_qubits, "base", 0);
  psi_1.SetRngPtr(&rnd_generator);
  // |psi_1> = |random>
  psi_1.Initialize("rand", 1);

  // |psi_2> = |psi_1> = |random>
  QubitRegister<ComplexDP> psi_2(psi_1);

  // Compare two implementations of the SWAP gate.
  unsigned qubit1 = 2, qubit2 = num_qubits-1;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
  qubit1 = 14;
  qubit2 = num_qubits-2;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
//  qubit1 = num_qubits-1; // FIXME: uncomment when global-global SWAP is available
  qubit2 = num_qubits-2;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
  qubit1 = num_qubits-3;
  qubit2 = num_qubits-2;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
  qubit1 = num_qubits-3;
  qubit2 = num_qubits-1;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
  qubit1 = num_qubits-3;
  qubit2 = 7;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );

  // Compare two implementations of the SWAP gate.
  qubit1 = 8;
  qubit2 = 13;
  psi_1.ApplySwap(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  psi_2.ApplyCPauliX(qubit2, qubit1);
  psi_2.ApplyCPauliX(qubit1, qubit2);
  // Check that the max abs difference amplitude by amplitude.
  ASSERT_DOUBLE_EQ(psi_2.MaxAbsDiff(psi_1), 0 );
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard APPLY_SWAP_GATE_TEST_HPP
