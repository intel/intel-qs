#ifndef GATE_COUNTER_TEST_HPP
#define GATE_COUNTER_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class GateCounterTest : public ::testing::Test
{
 protected:

  GateCounterTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // Hadamard gate defined via its matrix representation.
    ComplexDP amplitude = {1/std::sqrt(2), 0};
    G_(0, 0) =  amplitude; 
    G_(0, 1) =  amplitude;
    G_(1, 0) =  amplitude; 
    G_(1, 1) = -amplitude;
    // Two-qubit diagonal gate.
    D_(0, 0) = { 1,  0};
    D_(0, 0) = {-1,  0};
    D_(0, 0) = { 0,  1};
    D_(0, 0) = { 0, -1};
    D_(0, 1) = D_(0, 2) = D_(0, 3) = {0, 0}; 
    D_(1, 0) = D_(1, 2) = D_(1, 3) = {0, 0}; 
    D_(2, 0) = D_(2, 1) = D_(2, 3) = {0, 0}; 
    D_(3, 0) = D_(3, 1) = D_(3, 2) = {0, 0}; 
  }

  const std::size_t num_qubits_ = 10;
  TM2x2<ComplexDP> G_;
  TM4x4<ComplexDP> D_;
  double accepted_error_ = 1e-15;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(GateCounterTest, OnlyOneQubitGates)
{
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  // Initilize state randomly.
  std::size_t rng_seed = 977;
  iqs::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  psi.SetRngPtr(&rnd_generator);
  psi.Initialize("rand",1);
  EXPECT_NEAR( psi.ComputeNorm(), 1., accepted_error_);

  psi.EnableStatistics();

  for (int qubit=0; qubit<num_qubits_; ++qubit)
      psi.ApplyPauliX(qubit);

  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), num_qubits_);

  psi.ApplyHadamard(0);
  psi.ApplyRotationY(1, 1.35);
  psi.ApplyPauliSqrtZ(2);
  psi.ApplyT(5);

  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), num_qubits_+4);

  psi.Apply1QubitGate(7,G_);

  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), 0);
  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), psi.gate_counter->GetTotalGateCount());
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(GateCounterTest, ProvidedTwoQubitGates)
{
  int expected_counter = 0;
  // |psi> = |1010010000> = |"1+4+32">
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4+32);
  psi.EnableStatistics();
  for(int control = 0; control < num_qubits_; control++)
  {
      for(int target = control+1; target < num_qubits_; target++)
          psi.ApplyCPauliY(control, target);
  }
  expected_counter = (num_qubits_)*(num_qubits_-1);
  expected_counter /= 2;
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), expected_counter);

  psi.ApplyCHadamard(2,0);
  psi.ApplyCRotationY(5,1, 1.35);
  psi.ApplyCPauliSqrtZ(8,2);
  expected_counter +=3;
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), expected_counter);

  psi.ApplyDiag(0,1,D_);
  expected_counter +=1;
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), expected_counter);

  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), 0);
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), psi.gate_counter->GetTotalGateCount());

  // Reset statistics.
  psi.ResetStatistics();
  expected_counter = 0;

//FIXME TODO problem with SWAP gate
  psi.ApplySwap(1,5);
#if 0
//FIXME TODO proper counting
//  psi.ApplyISwap(5,num_qubits_-1);
//  psi.ApplySqrtISwap(8,4);
  expected_counter +=1;
  // Currently, the ApplySwap is implemented by applying 3 CNOTs. The gate counter is: 
//FIXME TODO proper counting
//  expected_counter +=2;
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), expected_counter);

  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), 0);
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), psi.gate_counter->GetTotalGateCount());
#endif
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(GateCounterTest, CustomTwoQubitGates)
{
  // Arbitrary two-qubit gates are implemented only when StateSize=1.
  if (iqs::mpi::Environment::GetStateSize() > 1)
      GTEST_SKIP() << "INFO: arbitrary 2-qubit gates require non-distributed states.";

  int expected_counter = 0;
  // |psi> = |1010010000> = |"1+4+32">
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4+32);
  psi.EnableStatistics();
  psi.Apply2QubitGate(7,1,D_);
  expected_counter +=1;
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), expected_counter);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(GateCounterTest, ToffoliGate)
{
  GTEST_SKIP() << "INFO: ToffoliGate is a test under development."; // FIXME: to check if action error comes from memory limits
  // |psi> = |0000000000> = |"0">
  std::cout << "Checkpoint A\n"; // FIXME: to check origin of github action error
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  // TODO: this test creates problem to the "C++ build with CMake" action in github
  //       despite it runs fine locally. It is difficult to debug since github
  //       actions cannot be manually triggered yet.
  std::cout << "Checkpoint B\n"; // FIXME: to check origin of github action error
  // Recall that the Toffoli gate is implemented by decomposing it in 5 two-qubit gates.
  psi.EnableStatistics();
  psi.ApplyToffoli(0,1,2);
  std::cout << "Checkpoint C\n"; // FIXME: to check origin of github action error
  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), 0);
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(), 5);
  psi.ApplyToffoli(5,8,3);
  ASSERT_EQ(psi.gate_counter->GetOneQubitGateCount(), 0);
  ASSERT_EQ(psi.gate_counter->GetTwoQubitGateCount(),10);
  ASSERT_EQ(psi.gate_counter->GetTotalGateCount(), 10);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard GATE_COUNTER_TEST_HPP
