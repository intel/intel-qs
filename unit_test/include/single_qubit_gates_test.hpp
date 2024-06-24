#ifndef SINGLE_QUBIT_GATES_TEST_HPP
#define SINGLE_QUBIT_GATES_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class SingleQubitGatesTest : public ::testing::Test
{
 protected:

  SingleQubitGatesTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
        GTEST_SKIP();

    // All tests are skipped if the 10-qubit state is distributed in more than 2^9-1 ranks.
    // In fact the MPI version needs to allocate half-the-local-storage for communication.
    // If the local storage is a single amplitude, this cannot be further divided.
    if (iqs::mpi::Environment::GetStateSize() > 511)
        GTEST_SKIP() << "INFO: small state distributed among too many ranks.";

    G_(0, 0) = {0.592056606032915, 0.459533060553574}; 
    G_(0, 1) = {-0.314948020757856, -0.582328159830658};
    G_(1, 0) = {0.658235557641767, 0.070882241549507}; 
    G_(1, 1) = {0.649564427121402, 0.373855203932477};
  }

  const std::size_t num_qubits_ = 10;
  TM2x2<ComplexDP> G_;
  double accepted_error_ = 1e-15;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(SingleQubitGatesTest, PauliOperators)
{
  // Recall that the qubits read from left to right, with the most significant
  // bit on the right (contrary to usual decimal representation).
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4+32+512);
  // |psi> = |1010010001> = |"1+4+32+512">
  ASSERT_DOUBLE_EQ( psi.ComputeNorm(), 1.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(0), 1.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(2), 1.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(num_qubits_-2), 0.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(num_qubits_-1), 1.);
  psi.ApplyPauliX(5);
  psi.ApplyPauliX(6);
  // |psi> = |1010001001> = |"1+4+64+512">
  ASSERT_DOUBLE_EQ(psi.GetProbability(5), 0.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(6), 1.);

  // |psi> = |0110001000> = |"2+4+64"> = |"70">
  psi.Initialize("base",2+4+64);
  ASSERT_DOUBLE_EQ( psi.ComputeNorm(), 1.);
  ASSERT_DOUBLE_EQ(psi.GetProbability(2),1);
  ASSERT_DOUBLE_EQ(psi.GetProbability(num_qubits_-4),1);
  psi.ApplyPauliZ(1);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(70).real(), -1.);
  psi.ApplyPauliZ(2);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(70).real(),  1.);
  psi.ApplyPauliZ(3);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(70).real(),  1.);
  psi.ApplyPauliX(4);
  // |psi> = |0110101000> = |"86">
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(86).real(),  1.);
  ASSERT_DOUBLE_EQ( psi.GetGlobalAmplitude(86).imag(),  0.);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(SingleQubitGatesTest, Rotations)
{
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4+32+512);
  // |psi> = |1010010001> = |"1+4+32+512">
  // Check that Rx(t)=exp(-i t X/2)
  double angle = M_PI;
  psi.ApplyRotationX(1, angle);
  psi.ApplyPauliX(1);
  // |psi> = -i|1010010001> = -i|"549">
  ASSERT_NEAR(psi.GetProbability(1), 0., accepted_error_);
  ASSERT_DOUBLE_EQ(psi.GetGlobalAmplitude(549).real(), 0.);
  ASSERT_DOUBLE_EQ(psi.GetGlobalAmplitude(549).imag(),-1.);

  angle = M_PI/2.34;
  psi.ApplyRotationX(0, angle);
  ASSERT_DOUBLE_EQ(psi.GetProbability(0), std::pow(std::cos(angle/2.),2) );
  //
  psi.ApplyRotationX(num_qubits_-2, angle);
  ASSERT_DOUBLE_EQ(psi.GetProbability(num_qubits_-2), std::pow(std::sin(angle/2.),2) );

  angle = M_PI*0.29;
  psi.ApplyHadamard(5);
  psi.ApplyRotationZ(5, angle);
  psi.ApplyHadamard(5);
  ASSERT_DOUBLE_EQ(psi.GetProbability(5), std::pow(std::cos(angle/2.),2) );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(SingleQubitGatesTest, CustomGate)
{
  // |psi> = |1010010001> = |"1+4+32+512">
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4+32+512);
  for(int qubit = 0; qubit < num_qubits_; qubit++)
  {
      psi.Apply1QubitGate(qubit, G_);
  }
  ASSERT_DOUBLE_EQ( psi.GetProbability(4), psi.GetProbability(6) );
  // The accepted error is sometimes exceeded when gcc is used.
  ASSERT_NEAR( psi.GetProbability(0), psi.GetProbability(2), accepted_error_*10 );
  ASSERT_NEAR( psi.GetProbability(5), 1.-psi.GetProbability(8), accepted_error_*10 );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(SingleQubitGatesTest, DeathTest)
{
  // Skip death-tests if compiler flag NDEBUG is defined.
#ifdef NDEBUG
  GTEST_SKIP() << "INFO: test skipped when compiler flag NDEBUG is not defined.";
#endif

  // Skip death-tests if MPI size > 1.
  if (iqs::mpi::Environment::GetStateSize() > 1)
      GTEST_SKIP();

  // |psi> = |0000000000> = |"0">
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);

  // To switch off the warning message about DEATH test not being thread safe.
  ::testing::FLAGS_gtest_death_test_style = "threadsafe"; 
  // Qubit index beyond the register size.
  int qubit = num_qubits_;
  ASSERT_DEATH( psi.ApplyHadamard(qubit), "");
  // Negative qubit index.
  qubit = -1;
  ASSERT_DEATH( psi.ApplyHadamard(qubit), "");
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard SINGLE_QUBIT_GATES_TEST_HPP
