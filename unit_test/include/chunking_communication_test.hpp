#ifndef CHUNKING_COMMUNICATION_TEST_HPP
#define CHUNKING_COMMUNICATION_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class ChunkingCommunicationTest : public ::testing::Test
{
 protected:

  ChunkingCommunicationTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 14-qubit state is distributed in more than 2^11 ranks.
    // In fact the MPI version needs to allocate eighth-the-local-storage for communication.
    // If the local storage is a eight amplitudes, this cannot be further divided.
    if (iqs::mpi::Environment::GetStateSize() > 2048)
      GTEST_SKIP() << "INFO: small state distributed among too many ranks.";
  }

  const std::size_t num_qubits_ = 14;
  double accepted_error_ = 1e-15;
  double accepted_error_loose_ = 1e-12;

  unsigned nprocs_ = iqs::mpi::Environment::GetStateSize();
  unsigned log2_nprocs_ = iqs::ilog2( iqs::floor_power_of_two(nprocs_));
  std::size_t local_size_ = std::size_t(1UL << std::size_t(num_qubits_ - log2_nprocs_));
  std::size_t tmp_spacesize_half_    = local_size_/2UL;
  std::size_t tmp_spacesize_quarter_ = local_size_/4UL;
  std::size_t tmp_spacesize_eighth_  = local_size_/8UL;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(ChunkingCommunicationTest, Initialization)
{
  // Consider 3 different sizes for the tmp_spacesize:
  // 1/2 (default with MPI), 1/4 , 1/8 of the local_size.
  iqs::QubitRegister<ComplexDP> psi_half    (num_qubits_, "base", 1, tmp_spacesize_half_   );
  iqs::QubitRegister<ComplexDP> psi_quarter (num_qubits_, "base", 10200, tmp_spacesize_quarter_);
  iqs::QubitRegister<ComplexDP> psi_eighth  (num_qubits_, "base", 5003, tmp_spacesize_eighth_ );

  // Test the norm.
  ASSERT_LE( std::abs(psi_half.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_quarter.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_eighth.ComputeNorm()-1) , accepted_error_ );

  ComplexDP amplitude;
  // Test the only non-zero entry.
  ASSERT_DOUBLE_EQ(psi_half.GetGlobalAmplitude(1).real(), 1.);
  ASSERT_DOUBLE_EQ(psi_quarter.GetGlobalAmplitude(10200).real(), 1.);
  ASSERT_DOUBLE_EQ(psi_eighth.GetGlobalAmplitude(5003).real(), 1.);

  // Apply gates such that |psi_half> = |+++...+>
  psi_half.ApplyPauliX(0);
  for (unsigned qubit=0; qubit<num_qubits_; ++qubit)
      psi_half.ApplyHadamard(qubit);
  // Verify overlap with a newly initialized state.
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "++++", 0);
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi), ComplexDP(1,0), accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ChunkingCommunicationTest, HadamardGate)
{
  // Consider 3 different sizes for the tmp_spacesize:
  // 1/2 (default with MPI), 1/4 , 1/8 of the local_size.
  iqs::QubitRegister<ComplexDP> psi_half    (num_qubits_, "base", 0, tmp_spacesize_half_   );
  iqs::QubitRegister<ComplexDP> psi_quarter (num_qubits_, "base", 0, tmp_spacesize_quarter_);
  iqs::QubitRegister<ComplexDP> psi_eighth  (num_qubits_, "base", 0, tmp_spacesize_eighth_ );

  // Test the norm.
  ASSERT_LE( std::abs(psi_half.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_quarter.ComputeNorm()-1) , accepted_error_ );
  ASSERT_LE( std::abs(psi_eighth.ComputeNorm()-1) , accepted_error_ );

  // Apply Hadamard on all qubits.
  for (unsigned qubit=0; qubit<num_qubits_; ++qubit)
  {
      psi_half.ApplyHadamard(qubit);
      psi_quarter.ApplyHadamard(qubit);
      psi_eighth.ApplyHadamard(qubit);
  }

  // Verify overlap.
  ComplexDP overlap = {1,0};
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi_quarter), overlap, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi_eighth), overlap, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_eighth.ComputeOverlap(psi_quarter), overlap, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ChunkingCommunicationTest, CustomGate)
{
  // Consider 3 different sizes for the tmp_spacesize:
  // 1/2 (default with MPI), 1/4 , 1/8 of the local_size.
  iqs::QubitRegister<ComplexDP> psi_half    (num_qubits_, "base", 0, tmp_spacesize_half_   );
  iqs::QubitRegister<ComplexDP> psi_quarter (num_qubits_, "base", 0, tmp_spacesize_quarter_);
  iqs::QubitRegister<ComplexDP> psi_eighth  (num_qubits_, "base", 0, tmp_spacesize_eighth_ );

  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};

  // Apply custom gate on all qubits.
  for (unsigned qubit=0; qubit<num_qubits_; ++qubit)
  {
      psi_half.Apply1QubitGate(qubit, G);
      psi_quarter.Apply1QubitGate(qubit, G);
      psi_eighth.Apply1QubitGate(qubit, G);
  }

  // Verify overlap.
  ComplexDP overlap = {1,0};
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi_quarter), overlap, accepted_error_loose_);
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi_eighth), overlap, accepted_error_loose_);
  ASSERT_COMPLEX_NEAR(psi_eighth.ComputeOverlap(psi_quarter), overlap, accepted_error_loose_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ChunkingCommunicationTest, CnotGate)
{
  // Consider 3 different sizes for the tmp_spacesize:
  // 1/2 (default with MPI), 1/4 , 1/8 of the local_size.
  iqs::QubitRegister<ComplexDP> psi_half    (num_qubits_, "base", 0, tmp_spacesize_half_   );
  iqs::QubitRegister<ComplexDP> psi_quarter (num_qubits_, "base", 0, tmp_spacesize_quarter_);
  iqs::QubitRegister<ComplexDP> psi_eighth  (num_qubits_, "base", 0, tmp_spacesize_eighth_ );

  // Apply CNOT on all (ordered) qubit-pairs for the highest 6 indices.
  for (unsigned q0=num_qubits_-6; q0<num_qubits_; ++q0)
  for (unsigned q1=q0+1; q1<num_qubits_; ++q1)
  {
      psi_half.ApplyCPauliX(q0, q1);
      psi_half.ApplyCPauliX(q1, q0);
      psi_quarter.ApplyCPauliX(q0, q1);
      psi_quarter.ApplyCPauliX(q1, q0);
      psi_eighth.ApplyCPauliX(q0, q1);
      psi_eighth.ApplyCPauliX(q1, q0);
  }

  // Verify overlap.
  ComplexDP overlap = {1,0};
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi_quarter), overlap, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_half.ComputeOverlap(psi_eighth), overlap, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi_eighth.ComputeOverlap(psi_quarter), overlap, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////
/*
TEST_F(ChunkingCommunicationTest, InitializeRandomlyButSame)
{
  // |psi> = |00>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  // random number generator
  std::size_t rng_seed = 7777;
  iqs::RandomNumberGenerator<double> rng;
  rng.SetSeedStreamPtrs(rng_seed);
  psi.SetRngPtr(&rng);
  //
  int num_states = iqs::mpi::Environment::GetNumStates();
  assert (num_states==1);
  psi.Initialize("rand",num_states);
  // |psi> = |rand>

  // Initilize the copy: |copy> = |psi>
  iqs::QubitRegister<ComplexDP> psi_copy (psi);
  ASSERT_DOUBLE_EQ(psi_copy.MaxAbsDiff(psi), 0 );
  ASSERT_DOUBLE_EQ(psi_copy.MaxL2NormDiff(psi), 0 );
  //
  // Reinitialize |copy> by generating its amplitudes.
  iqs::RandomNumberGenerator<double> rng_copy;
  rng_copy.SetSeedStreamPtrs(rng_seed);
  psi_copy.SetRngPtr(&rng_copy);
  psi_copy.Initialize("rand",num_states);
  ASSERT_DOUBLE_EQ(psi_copy.MaxAbsDiff(psi), 0 );
  ASSERT_DOUBLE_EQ(psi_copy.MaxL2NormDiff(psi), 0 );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ChunkingCommunicationTest, Hadamard)
{
  iqs::QubitRegister<ComplexDP> psi_0 (num_qubits_,"base",0);
  psi_0.ApplyHadamard(0);
  // |psi_0> = |+0> = |q0=+> x |q1=0>
  ASSERT_NEAR( psi_0.GetProbability(0), 0.5, accepted_error_ );
  psi_0.ApplyPauliZ(0);
  // |psi_0> = |-0> = |q0=-> x |q1=0>
  ASSERT_NEAR( psi_0.GetProbability(0), 0.5, accepted_error_ );
  psi_0.ApplyHadamard(0);
  // |psi_0> = |10> = |q0=1> x |q1=0>
  ASSERT_NEAR( psi_0.GetProbability(0), 1. , accepted_error_ );

  iqs::QubitRegister<ComplexDP> psi_1 (num_qubits_,"base",1);
  psi_1.ApplyHadamard(0);
  // |psi_1> = |-0> = |q0=-> x |q1=0>
  ASSERT_NEAR( psi_1.GetProbability(0), 0.5, accepted_error_ );

  iqs::QubitRegister<ComplexDP> psi_2 (num_qubits_,"base",2);
  psi_2.ApplyHadamard(1);
  // |psi_2> = |0-> = |q0=0> x |q1=->
  ComplexDP amplitude = ComplexDP(1./std::sqrt(2.), 0. );
  ASSERT_EQ(psi_2.GetGlobalAmplitude(0), amplitude);
  ASSERT_DOUBLE_EQ(psi_2.GetGlobalAmplitude(1).real(), 0.);
  ASSERT_EQ(psi_2.GetGlobalAmplitude(2),-amplitude);
  ASSERT_DOUBLE_EQ(psi_2.GetGlobalAmplitude(3).imag(), 0.);

  iqs::QubitRegister<ComplexDP> psi_3 (num_qubits_,"base",3);
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
*/
//////////////////////////////////////////////////////////////////////////////

#endif	// header guard CHUNKING_COMMUNICATION_TEST_HPP
