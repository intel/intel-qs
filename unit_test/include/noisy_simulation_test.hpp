#ifndef NOISY_SIMULATION_TEST_HPP
#define NOISY_SIMULATION_TEST_HPP

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class NoisySimulationTest : public ::testing::Test
{
 protected:

  NoisySimulationTest()
  {
#ifdef INTELQS_HAS_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);
    MPI_Comm_rank(MPI_COMM_WORLD, &pool_rank_id_);
#else
    num_ranks_ = 1;
    pool_rank_id_ = 0;
#endif

    // To ensure that each of two states has at least 2 amplitude per rank.
    int min_num_qubits = qhipster::floor_power_of_two( num_ranks_) + 1;
    if (num_qubits_ < min_num_qubits)
        num_qubits_ = min_num_qubits;
  }

  // Just after the 'constructor'.
  void SetUp() override
  {
      // This kind of tests makes no sense without MPI and at least 2 ranks.
      if (num_ranks_<2)
          GTEST_SKIP();
  }

  // Just before the 'destructor'.
  void TearDown() override
  {
     if (qhipster::mpi::Environment::GetNumStates() != 1)
         qhipster::mpi::Environment::UpdateStateComm(1,false);
  }

  int num_qubits_= 6;
  double T1_=6.;
  double T2_=4.;
  double accepted_error_ = 1e-15;
  int pool_rank_id_;
  int num_ranks_;
  bool do_print_info_ = true;	// Whether printing info when the PoolComm is restructured.
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(NoisySimulationTest, OneStateAtATime)
{
  ASSERT_EQ( qhipster::mpi::Environment::GetNumStates(), 1);
  // Currently pool=state, but not always pool=MPI_COMM_WORLD since pool is
  // defined only for the useful ranks.

  if (qhipster::mpi::Environment::IsUsefulRank() == false)
      return;

  QubitRegister<ComplexDP> psi (num_qubits_,"base",1+8+16+32);
  // |psi> = |100111> = |"1+8+16+32">
  psi.ApplyHadamard(0);
  psi.ApplyHadamard(1);
  // |psi> = |-+0111>

  QubitRegister<ComplexDP> noisy_psi (psi);
  ASSERT_DOUBLE_EQ( noisy_psi.ComputeOverlap(psi).real(), 1.);
  // Set the dissipation and decoherence times.
  noisy_psi.SetNoiseTimescales(T1_,T2_);
  // Noise gates require random numbers.
  std::size_t rng_seed = 7777;
  qhipster::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  noisy_psi.SetRngPtr(&rnd_generator);
  // A certain time duration is spend with all qubits idle.
  double duration=5;
  for (int qubit=0; qubit<num_qubits_; ++qubit)
      noisy_psi.ApplyNoiseGate(qubit,duration);
  ASSERT_TRUE( noisy_psi.ComputeOverlap(psi).real() < 1.-accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(NoisySimulationTest, TwoStates)
{
  ASSERT_EQ( qhipster::mpi::Environment::GetNumStates(), 1);
  // Update state commutator.
  int num_states = 2;
  qhipster::mpi::Environment::UpdateStateComm(num_states, do_print_info_);
  ASSERT_EQ( num_states, qhipster::mpi::Environment::GetNumStates() );
  if (qhipster::mpi::Environment::IsUsefulRank() == false)
      return;

  int my_state_id = qhipster::mpi::Environment::GetStateId();
  std::size_t index = my_state_id;
  QubitRegister<ComplexDP> psi (num_qubits_, "base", index);

  // The pool has two states, |0> and |1>.
  int qubit = 0;
  double incoherent_sum, probability;
  probability = psi.GetProbability(qubit);
  ASSERT_DOUBLE_EQ( double(my_state_id), probability );
  // Sum up the probabilities incoherently.
  incoherent_sum
    = qhipster::mpi::Environment::IncoherentSumOverAllStatesOfPool<double>(probability);
  ASSERT_DOUBLE_EQ( incoherent_sum, 1. );

  // Now initialize all states of the pool to |"0">.
  psi.Initialize("base",0);
  QubitRegister<ComplexDP> noisy_psi (psi);
  // Noise gates require random numbers.
  std::size_t rng_seed = 7777;
  qhipster::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  noisy_psi.SetRngPtr(&rnd_generator);
  // If purely dissipation, the population in state |0> should increase.
  psi.SetNoiseTimescales(T1_,T1_/2.);
  double duration=T1_;
  for (int q=0; q<num_qubits_; ++q)
      noisy_psi.ApplyNoiseGate(q,duration);
  probability = noisy_psi.GetProbability(qubit);
  incoherent_sum
    = qhipster::mpi::Environment::IncoherentSumOverAllStatesOfPool<double>(probability);
  double incoherent_average = incoherent_sum / double(num_states);
//qhipster::mpi::PoolPrint("~~~~ prob : " + std::to_string(probability), true);//FIXME
//qhipster::mpi::PoolPrint("~~~~ aver : " + std::to_string(incoherent_average), true);//FIXME
  ASSERT_GT( incoherent_average, accepted_error_ );
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(NoisySimulationTest, OneStatePerRank)
{
  ASSERT_EQ( qhipster::mpi::Environment::GetNumStates(), 1);
  // Update state commutator.
  int num_states = num_ranks_;
  qhipster::mpi::Environment::UpdateStateComm(num_states,do_print_info_);
  ASSERT_EQ( qhipster::mpi::Environment::GetStateRank(), 0 );
  ASSERT_EQ( qhipster::mpi::Environment::GetStateSize(), 1 );
  ASSERT_EQ( qhipster::mpi::Environment::GetPoolRank(), pool_rank_id_ );
  ASSERT_EQ( qhipster::mpi::Environment::GetPoolSize(), num_states );

  // Initialize all states in different computational basis states (if possible).
  int my_state_id = qhipster::mpi::Environment::GetStateId();
  ASSERT_EQ( my_state_id , pool_rank_id_);

  if ( qhipster::mpi::Environment::IsUsefulRank() )
  {
      QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
      ASSERT_EQ( psi.GlobalSize(), psi.LocalSize() );
      std::size_t index = my_state_id % psi.GlobalSize();
      psi.Initialize("base",index);
      // At this point: my_state_id=0 --> |0>
      //                my_state_id=k --> |k%2^n>
      int qubit = 0;
      double probability = psi.GetProbability(qubit);
      ASSERT_DOUBLE_EQ(probability, double(index%2));
      double incoherent_sum
        = qhipster::mpi::Environment::IncoherentSumOverAllStatesOfPool<double>(probability);
      if (num_states%2==0)
          ASSERT_DOUBLE_EQ( incoherent_sum, double(num_states  )/2. );
      else
          ASSERT_DOUBLE_EQ( incoherent_sum, double(num_states-1)/2. );
  }
  else
      ASSERT_TRUE(false);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard NOISY_SIMULATION_TEST_HPP
