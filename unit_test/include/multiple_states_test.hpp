#ifndef MULTIPLE_STATES_TEST_HPP
#define MULTIPLE_STATES_TEST_HPP

// This file is included only when INTELQS_HAS_MPI is defined.

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class MultipleStatesTest : public ::testing::Test
{
 protected:

  MultipleStatesTest()
  {
#ifdef INTELQS_HAS_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);
    MPI_Comm_rank(MPI_COMM_WORLD, &pool_rank_id_);
#else
    num_ranks_ = 1;
    pool_rank_id_ = 0;
#endif
    // To ensure that each of two states has at least 2 amplitude per rank.
    int min_num_qubits = qhipster::ilog2( qhipster::floor_power_of_two( num_ranks_) ) + 1;
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

  int num_qubits_= 8;
  double accepted_error_ = 1e-15;
  int pool_rank_id_;
  int num_ranks_;
  bool do_print_info_ = false;	// Whether printing info when the PoolComm is restructured.
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(MultipleStatesTest, OneStatePerRank)
{
  ASSERT_EQ( qhipster::mpi::Environment::GetNumStates(), 1);
  // Currently pool=state, but not always pool=MPI_COMM_WORLD since pool is
  // defined only for the useful ranks.

  if ( qhipster::mpi::Environment::IsUsefulRank() )
  {
    ASSERT_EQ( qhipster::mpi::Environment::GetPoolSize(),
               qhipster::mpi::Environment::GetStateSize() );
    ASSERT_EQ( qhipster::mpi::Environment::GetPoolRank(),
               qhipster::mpi::Environment::GetStateRank() );
  }
  
  // Update state commutator.
  int num_states = num_ranks_;
  qhipster::mpi::Environment::UpdateStateComm(num_states,do_print_info_);
  ASSERT_EQ( qhipster::mpi::Environment::GetStateRank(), 0 );
  ASSERT_EQ( qhipster::mpi::Environment::GetStateSize(), 1 );
  ASSERT_EQ( qhipster::mpi::Environment::GetPoolRank(), pool_rank_id_ );

  // Initialize all states in different computational basis states (if possible).
  int my_state_id = qhipster::mpi::Environment::GetStateId();
  ASSERT_EQ( my_state_id , pool_rank_id_);

  if ( qhipster::mpi::Environment::IsUsefulRank() )
  {
      QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
      ASSERT_EQ( psi.GlobalSize(), psi.LocalSize() );
      std::size_t index = my_state_id % psi.GlobalSize();
      psi.Initialize("base",index);
      ComplexDP amplitude = psi.GetGlobalAmplitude(index);
      ASSERT_EQ(amplitude.real(),1.);
      ASSERT_EQ(amplitude.imag(),0.);

      psi.ApplyPauliZ(1);
      if (my_state_id%4==2 || my_state_id%4==3)
          ASSERT_EQ( psi.GetGlobalAmplitude(index).real(),-1.);
      else
          ASSERT_EQ( psi.GetGlobalAmplitude(index).real(), 1.);

      // Apply non-diagonal gates. They never require communication.
      psi.ApplyHadamard(1);
      psi.ApplyPauliX(1);
      psi.ApplyHadamard(1);
      amplitude = psi.GetGlobalAmplitude(index);
      ASSERT_DOUBLE_EQ(amplitude.real(),1.);
      ASSERT_DOUBLE_EQ(amplitude.imag(),0.);
  }
  else
      ASSERT_TRUE(false);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(MultipleStatesTest, TwoStates)
{
  ASSERT_EQ( qhipster::mpi::Environment::GetNumStates(), 1);
  // Update state commutator.
  int num_states = 2;

  qhipster::mpi::Environment::UpdateStateComm(num_states,do_print_info_);

  ASSERT_EQ( num_states, qhipster::mpi::Environment::GetNumStates() );

  int my_state_id = qhipster::mpi::Environment::GetStateId();
  int num_useful_ranks = 0;
  int num_dummy_ranks  = 0;

  if (qhipster::mpi::Environment::IsUsefulRank() == true)
  {
      num_useful_ranks =   qhipster::mpi::Environment::GetNumStates()
                           * qhipster::mpi::Environment::GetStateSize();
#if 0
std::cout << "state_id= " << my_state_id << " over "
          << qhipster::mpi::Environment::GetNumStates()
          << "  ,  pool_size = " << qhipster::mpi::Environment::GetPoolSize() << "\n";
#endif
      ASSERT_LE(num_useful_ranks, num_ranks_ );
      ASSERT_EQ(num_useful_ranks, qhipster::mpi::Environment::GetPoolSize() );
  }
  else
  {
      num_dummy_ranks = qhipster::mpi::Environment::GetPoolSize();
      ASSERT_LE(num_dummy_ranks, num_ranks_ );
      ASSERT_LE(num_dummy_ranks, qhipster::mpi::Environment::GetPoolSize() );
  }

  // The tests are executed only if there are at least 2 MPI ranks.
  // However to use MPI types and variables, one needs to have imported <mpi.h>.
  // Use global rank=0 to sum dummy and useful ranks.
  int tag = 0;
  // Recall that GetPoolComm() is defined only for useful ranks.
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Barrier(comm);
  if (pool_rank_id_==0)
  {
      MPI_Recv(&num_dummy_ranks, 1, MPI_INT, num_ranks_-1, tag, comm, MPI_STATUS_IGNORE);
      ASSERT_EQ(num_dummy_ranks + num_useful_ranks, num_ranks_);
  }
  else if (pool_rank_id_==num_ranks_-1)
  {
      // num_dummy_ranks > 0 only if there are unuseful ranks and for a dummy rank.
      // The global rank (num_ranks_-1) must be dummy, if any rank is.
      MPI_Send(&num_dummy_ranks, 1, MPI_INT, 0, tag, comm);
  }
#endif

  // Execute simple quantum computations.
  if ( qhipster::mpi::Environment::IsUsefulRank() )
  {
      QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
      int my_state_id = qhipster::mpi::Environment::GetStateId();
      std::size_t index = my_state_id % psi.GlobalSize();
      psi.Initialize("base",index);
      // At this point:  state_id=0 --> |0>
      //                 state_id=1 --> |1>
      ComplexDP amplitude = psi.GetGlobalAmplitude(index);
      ASSERT_EQ(amplitude.real(),1.);
      ASSERT_EQ(amplitude.imag(),0.);

      psi.ApplyPauliZ(1);
      if (my_state_id%4==2 || my_state_id%4==3)
          ASSERT_EQ( psi.GetGlobalAmplitude(index).real(),-1.);
      else
          ASSERT_EQ( psi.GetGlobalAmplitude(index).real(), 1.);

      // When num_useful_ranks>2, the last qubit is global.
      psi.ApplyCPauliX(0,num_qubits_-1);
      ASSERT_DOUBLE_EQ(psi.GetProbability(num_qubits_-1),double(index));
      // Undo all the computation performed after the initialization.
      psi.ApplyHadamard(num_qubits_-1);
      psi.ApplyCPauliZ(0,num_qubits_-1);
      psi.ApplyHadamard(num_qubits_-1);
      psi.ApplyPauliZ(1);
      amplitude = psi.GetGlobalAmplitude(index);
      ASSERT_DOUBLE_EQ(amplitude.real(),1.);
      ASSERT_DOUBLE_EQ(amplitude.imag(),0.);
  }
  else
      ASSERT_TRUE( num_dummy_ranks>0 );
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard MULTIPLE_STATES_TEST_HPP
