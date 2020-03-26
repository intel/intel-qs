/// @file mpi_env.cpp
///
///  This header implements MPI support functionality

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#include <stdexcept>
#include <vector>
#endif

#include <iomanip>
#include <iostream>
#include <thread>

// In the previous release of IQS, alternative implementations of these include files
// were adapted from the OpenQu project. This is not necessary anymore. The compiler's
// option 'STANDALONE' is now deprecated since Those contributions have been removed.
#include "../include/mpi_env.hpp"
#include "../include/mpi_utils.hpp"
#include "../include/bitops.hpp"
#include "../include/conversion.hpp"

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

namespace mpi {

/////////////////////////////////////////////////////////////////////////////////////////
// Mocking methods when MPI is not present.
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef INTELQS_HAS_MPI

Environment::Environment(int&, char**&)
{ useful_rank = true; }

Environment::~Environment() {}

void Environment::UpdateStateComm(int new_num_states, bool do_print_info)
{ assert(new_num_states==1);  printf("MPI not enabled.\n"); }

int Environment::GetPoolRank() {return 0;}
int Environment::GetStateRank() {return 0;}

int Environment::GetPoolSize() {return 1;}
int Environment::GetStateSize() {return 1;}

int Environment::GetNumRanksPerNode() {return 1;}

void Environment::RemapStateRank(int newme) {}

template <class Type>
Type Environment::IncoherentSumOverAllStatesOfPool(Type local_value)
{ return local_value; }
//
template float  Environment::IncoherentSumOverAllStatesOfPool<float>  (float  );
template double Environment::IncoherentSumOverAllStatesOfPool<double> (double );

void StateBarrier() {}
void PoolBarrier() {}

void PoolPrint(std::string s, bool all)
{ StatePrint(s,all); }

void StatePrint(std::string s, bool all)
{
  int rank = Environment::GetStateRank();

  if (all) {
    if (rank == 0) {printf("[|%d>:%3d] %s\n", Environment::GetStateId(), rank,
                                               s.c_str()); fflush(stdout);}
  } else {
    if (rank == 0) std::cout << s << std::endl;
  }
}

#endif

/////////////////////////////////////////////////////////////////////////////////////////
// A few methods have the same implementation with / without MPI.
/////////////////////////////////////////////////////////////////////////////////////////

void Print(std::string s, bool all) {StatePrint(s,all);}

void Barrier() {StateBarrier();}

int Environment::GetNumNodes() {return num_nodes;}
int Environment::GetNodeId() {return my_node_id;}

int Environment::GetNumStates() {return num_states;}
int Environment::GetStateId() {return my_state_id;}

/////////////////////////////////////////////////////////////////////////////////////////

bool Environment::useful_rank = true;

int Environment::num_ranks_per_node = 1;
int Environment::num_nodes = 1;
int Environment::my_node_id = 0;
int Environment::num_states = 1;
int Environment::my_state_id = 0;

/////////////////////////////////////////////////////////////////////////////////////////
// Actual implementation with MPI.
/////////////////////////////////////////////////////////////////////////////////////////

#ifdef INTELQS_HAS_MPI

// pool_communicator and state_communicator are static variables that need to be defined.
MPI_Comm Environment::pool_communicator = MPI_COMM_WORLD;
MPI_Comm Environment::state_communicator = MPI_COMM_WORLD;

/////////////////////////////////////////////////////////////////////////////////////////

Environment::Environment(int& argc, char**& argv) : inited_(false)
{
  int flag;
  QHIPSTER_MPI_CHECK_RESULT(MPI_Initialized,(&flag))
  if (!flag) {
    QHIPSTER_MPI_CHECK_RESULT(MPI_Init,(&argc, &argv))
    inited_ = true;
#if 0
#if defined(MVAPICH2_VERSION) 
    char * mv2_string; 
    int mv2_affinity = 1; /* this is the default behavior of MVAPICH2 */ 
    if ((mv2_string = getenv("MV2_ENABLE_AFFINITY")) != NULL)
    { 
      mv2_affinity = atoi(mv2_string); 
    } 
    if (mv2_affinity!=0 /* && procid==0 FIXME */)
    { 
      printf("WARNING: You are using MVAPICH2 with affinity enabled, probably by default.\n"); 
      printf("WARNING: This will cause performance issues for ARMCI. \n"); 
      printf("WARNING: Please rerun your job with MV2_ENABLE_AFFINITY=0 \n"); 
    } 
#endif 
#endif 
  }

  UpdateStateComm(num_states);

// TODO: eliminate if unnecessary
#if 0
#if 0
  int nuseful_ranks_per_node = 10000000;
  int first_rank_on_last_node = nranks - num_ranks_per_node;
  if(useful_rank) {
    if (my_rank < first_rank_on_last_node)
      nuseful_ranks_per_node = num_ranks_per_node;
    else
      nuseful_ranks_per_node = (usefull_nranks - (nranks - num_ranks_per_node));
  }

  int threads_per_rank = (std::thread::hardware_concurrency() / 2) / nuseful_ranks_per_node;
#else
  int threads_per_rank = (std::thread::hardware_concurrency() / 2) / num_ranks_per_node;
#endif
  assert(threads_per_rank == qhipster::openmp::omp_get_set_num_threads());
#endif

  // start synching all threads
  // MPI_Ibarrier(MPI_COMM_WORLD, &synch_request);
}

/////////////////////////////////////////////////////////////////////////////////////////

Environment::~Environment()
{
  if (inited_) { 
#if 0
    if (IsUsefulRank()) {
      MPI_Wait( &synch_request, MPI_STATUS_IGNORE );
    } else {
      int flag = 0;
      while (!flag) {
        MPI_Test( &synch_request, &flag, MPI_STATUS_IGNORE );
        // sleep(1); /* this will waste as much as a minute at the end of your job */
      }
    }
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////

// NOTE: One needs to know what he/she is doing to call this method!
void Environment::UpdateStateComm(int new_num_states, bool do_print_info)
{
  // To make the implementation more resistant, we impose a few constraints:
  // 1. same number of ranks per node;
  // 2. same number of ranks per state;
  // 3. each state has a number of ranks that is a power of 2;
  // 4. remaining ranks are 'dummy';
  // We do not impose that:
  // -. each node hosts only ranks belonwing to the same state;

  // The recommended use has the ideal values:
  //   pool_size = 2^k * num_states
  // with
  //   2^k % num_ranks_per_node = 0
  // possibly with
  //   num_ranks_per_node = 2^j
  // In the latter case: num_qubits = j+k

  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int my_node_rank;

  // Create a new communicator that includes only those ranks on the same node
  // (i.e. able to have shared memory).
  MPI_Comm node_communicator;
  MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank,
                       MPI_INFO_NULL, &node_communicator );
  MPI_Comm_rank(node_communicator, &my_node_rank);
  MPI_Comm_size(node_communicator, &num_ranks_per_node);

  // Compute number of ranks per node and impose that they are the same for every node.
  assert((world_size % num_ranks_per_node) == 0);
  // Basic info on the node id and number.
  my_node_id = world_rank / num_ranks_per_node;
  num_nodes  = world_size / num_ranks_per_node;

  // Release the node_communicator.
  MPI_Comm_free(&node_communicator);

  /////

  // The ranks are divided as follows (10 ranks and 2 states):
  //   rank_id =  0  1  2  3  4  5  6  7  8  9
  //  state_id =  0  0  0  0  1  1  1  1  .  .

  num_states = new_num_states;
  int num_ranks_per_state = qhipster::floor_power_of_two(world_size/num_states);
  my_state_id = world_rank / num_ranks_per_state;
  // Reset useful_rank.
  useful_rank = true;

  MPI_Barrier(MPI_COMM_WORLD);
  if ( num_ranks_per_state * num_states == world_size )
  {
    // All ranks are useful. Update the state_communicator.
//state_communicator = MPI_COMM_WORLD;// FIXME
    MPI_Comm_split( MPI_COMM_WORLD, my_state_id, world_rank, &state_communicator );
    // The pool rank is the world commutator.
    pool_communicator = MPI_COMM_WORLD;
  }
  else
  {
    // Not all ranks are useful.
    // Update the state_communicator, possibly to a dummy one for the non-useful ranks.
    MPI_Group  world_group, new_group;

    // Record if the rank is useful or dummy.
    useful_rank = (world_rank < num_ranks_per_state * num_states);

    // Extract the original group handle.
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    // Pool communicator (only useful ranks):
    std::vector<int> useful_state_ranks;
    std::vector<int> dummy_ranks;
    int tag;
    if (useful_rank)
    {
        tag = 0;
        for(int i = 0; i < num_ranks_per_state * num_states; i++)
            useful_state_ranks.push_back(i);

        // Select the ranks that will be part of the pool_communicator.
        MPI_Group_incl( world_group, num_ranks_per_state * num_states,
                        useful_state_ranks.data(), &new_group );
        // Create the appropriate pool_communicator.
        MPI_Comm_create_group(MPI_COMM_WORLD, new_group, tag, &pool_communicator);

        // Verify that the size of the pool_communicator is the expected one. FIXME
//        MPI_Comm_size(pool_communicator, &tag);
//        assert(tag == num_ranks_per_state * num_states);
    }
    else
    {
        tag = 1;
        for(int i = num_ranks_per_state*num_states; i < world_size; i++)
            dummy_ranks.push_back(i);
        // The dummy ranks have a dummy pool communicator.
        MPI_Group_incl( world_group, dummy_ranks.size(), dummy_ranks.data(), &new_group);
        // Create the dummy pool_communicator.
        MPI_Comm_create_group(MPI_COMM_WORLD, new_group, tag, &pool_communicator);
    }

    // State communicator:
    if (useful_rank)
    {
        tag = my_state_id;
        useful_state_ranks.clear();
        for(int i = 0; i < num_ranks_per_state; i++)
            useful_state_ranks.push_back(i + num_ranks_per_state*my_state_id);

        // Select the ranks that will be part of the state_communicator.
        MPI_Group_incl( world_group, num_ranks_per_state, useful_state_ranks.data(),
                        &new_group );
        // Create the appropriate state_communicator.
        MPI_Comm_create_group(MPI_COMM_WORLD, new_group, tag, &state_communicator);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /////

  // Print information to screen, without OpenMP affinity.

  int threads_per_rank = 1;
#ifdef _OPENMP
#pragma omp parallel
  {
      threads_per_rank = omp_get_num_threads();
  }
#endif

  int my_state_rank;
  if (IsUsefulRank()==true)
      my_state_rank = GetStateRank();
  else
      my_state_rank = -1;

  std::stringstream buffer;
#ifndef NDEBUG
  // Print to screen all information, useful for debug.
  buffer    <<  "world_rank: " << std::setw(4) << qhipster::toString(world_rank) << " ,"
            <<  " state_rank: " << std::setw(4) << qhipster::toString(my_state_rank)
            <<  " (state " << std::setw(3) << qhipster::toString(my_state_id) 
            <<  " of " << std::setw(3) << qhipster::toString(num_states) << ")"
            <<  " my_node_id: " << std::setw(4) << qhipster::toString(my_node_id)
            <<  " , num_nodes: " << std::setw(4) << qhipster::toString(num_nodes)
            <<  " , ranks/node: " << std::setw(2) << qhipster::toString(num_ranks_per_node) 
            <<  " , threads/rank: " << std::setw(2) << qhipster::toString(threads_per_rank)
            <<  (useful_rank ? " --useful" : " --dummy");
#else
  // Print to screen a few information
  buffer    <<  "world_rank: " << std::setw(4) << qhipster::toString(world_rank)
            <<  "/" << world_size << " ,"
            <<  " state_rank: " << std::setw(4) << (IsUsefulRank() ? qhipster::toString(my_state_rank) : "   -")
            <<  "/" << num_ranks_per_state << " ,"
            <<  " node_id: " << std::setw(4) << qhipster::toString(my_node_id) 
            <<  "/" << qhipster::toString(num_nodes)
            <<  " , ranks/node: " << std::setw(2) << qhipster::toString(num_ranks_per_node) 
            <<  " , threads/rank: " << std::setw(2) << qhipster::toString(threads_per_rank)
            <<  (useful_rank ? " --useful" : " --dummy");
#endif
  if (do_print_info==true)
      Print(buffer.str(), MPI_COMM_WORLD);
}

/////////////////////////////////////////////////////////////////////////////////////////

int Environment::GetPoolRank()
{
  int rank;
  QHIPSTER_MPI_CHECK_RESULT(MPI_Comm_rank,(pool_communicator, &rank))
  return rank;
}
//
int Environment::GetStateRank()
{
  int rank;
  QHIPSTER_MPI_CHECK_RESULT(MPI_Comm_rank,(state_communicator, &rank))
  return rank;
}

/////////////////////////////////////////////////////////////////////////////////////////

int Environment::GetPoolSize()
{
  int size;
  QHIPSTER_MPI_CHECK_RESULT(MPI_Comm_size,(pool_communicator, &size))
  return size;
}
//
int Environment::GetStateSize()
{
  int size;
  QHIPSTER_MPI_CHECK_RESULT(MPI_Comm_size,(state_communicator, &size))
  return size;
}

/////////////////////////////////////////////////////////////////////////////////////////

int Environment::GetNumRanksPerNode()
{
  MPI_Comm node_communicator;
  int num_ranks_per_node;
  int my_rank = GetPoolRank();
  MPI_Comm_split_type( pool_communicator, MPI_COMM_TYPE_SHARED, my_rank,
                       MPI_INFO_NULL, &node_communicator );
  MPI_Comm_size(node_communicator, &num_ranks_per_node);
  MPI_Comm_free(&node_communicator);

  return num_ranks_per_node;
}

/////////////////////////////////////////////////////////////////////////////////////////

void Environment::RemapStateRank(int new_rank_id)
{
  MPI_Comm new_comm;
  MPI_Comm_split(state_communicator, 0, new_rank_id, &new_comm);
  // int r = GetStateRank();
  // printf("state rank=%d changing to %d, oldcomm=%d newcom=%d\n",
  //        r, new_rank_id, state_communicator, newcomm);
  state_communicator = new_comm;
}

/////////////////////////////////////////////////////////////////////////////////////////

template <class Type>
Type Environment::IncoherentSumOverAllStatesOfPool(Type local_value)
{
  // Only the main rank of each state has to send its value to the main global rank.
  // However both GetProbability() and ExpectationValue() actually redistribute their
  // result to each of the state ranks (using MPI_Allreduce).
  // Therefore we can simply accumulate from all ranks and divide by:
  //   (number of states) * (num ranks per state) = (size of pool communicator)
  Type global_value ;
  MPI_Comm comm = qhipster::mpi::Environment::GetPoolComm();
  MPI_Allreduce_x(&local_value, &global_value, 1, MPI_SUM, comm);
  global_value /= Type(qhipster::mpi::Environment::GetStateSize());

// If one were to compute the average instead of the sum, we need to divide by:
//    (number of states) * (num ranks per state) = (size of pool communicator)
//  global_value /= Type(qhipster::mpi::Environment::GetPoolSize());
//  global_value /= Type(qhipster::mpi::Environment::GetNumStates());
  assert( GetStateSize() * GetNumStates() == GetPoolSize() );
  return global_value;
}
//
template float  Environment::IncoherentSumOverAllStatesOfPool<float>  (float  );
template double Environment::IncoherentSumOverAllStatesOfPool<double> (double );

/////////////////////////////////////////////////////////////////////////////////////////

MPI_Comm Environment::GetPoolComm()  {return pool_communicator;}
MPI_Comm Environment::GetStateComm() {return state_communicator;}
MPI_Comm Environment::GetComm()      {return GetStateComm();}

/////////////////////////////////////////////////////////////////////////////////////////

void PoolBarrier()
{
  MPI_Comm comm = Environment::GetPoolComm(); 
  MPI_Barrier(comm);
}
//
void StateBarrier()
{
  MPI_Comm comm = Environment::GetStateComm(); 
  MPI_Barrier(comm);
}

/////////////////////////////////////////////////////////////////////////////////////////

void PoolPrint(std::string s, bool all)
{ Print(s, Environment::GetPoolComm(), all); }

#if 1
void StatePrint(std::string s, bool all)
{ Print(s, Environment::GetStateComm(), all); }
#else
void StatePrint(std::string s, bool all)
{
  int rank = Environment::GetStateRank();
  int size = Environment::GetStateSize();

  if (all)
  {
    StateBarrier();
    if (rank == 0) {printf("[|%d>:%3d] %s\n", Environment::GetStateId(), rank,
                                              s.c_str()); fflush(stdout);}
    MPI_Comm comm = Environment::GetStateComm(); 
    std::vector<char> buffer;
    if (rank == 0)
    {
      for (int i = 1; i < size; i++)
      {
        MPI_Status status;
        QHIPSTER_MPI_CHECK_RESULT((MPI_Probe,(i, i, comm, &status)));
        int cnt = 0;
        QHIPSTER_MPI_CHECK_RESULT((MPI_Get_count,(&status, MPI_CHAR, &cnt)));
        buffer.resize(cnt);
        QHIPSTER_MPI_CHECK_RESULT(
            MPI_Recv,(&buffer[0], cnt, MPI_BYTE, i, i, comm, MPI_STATUS_IGNORE));
        printf("[|%d>:%3d] %s\n", Environment::GetStateId(), i, (char*)(&buffer[0])); 
        fflush(stdout);
      }
    }
    else
    {
      QHIPSTER_MPI_CHECK_RESULT((MPI_Send,(const_cast<char*>(s.c_str()), s.size() + 1,
                                           MPI_CHAR, 0, rank, comm)));
    }
  }
  else
  {
      StateBarrier();
      if (rank == 0) std::cout << s << std::endl;
      StateBarrier();
  }
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////

void Print(std::string s, MPI_Comm comm, bool all)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  int sender_state_id = Environment::GetStateId();

  if (all)
  {
      MPI_Barrier(comm);
      if (rank == 0)
      {
          printf("[|s%d>:%3d] %s\n", sender_state_id, rank, s.c_str());
          fflush(stdout);
      }

      std::vector<char> buffer;
      if (rank == 0)
      {
          for (int i = 1; i < size; i++)
          {
              MPI_Status status;
              QHIPSTER_MPI_CHECK_RESULT(MPI_Probe,(i, i, comm, &status));
              int cnt = 0;
              QHIPSTER_MPI_CHECK_RESULT(MPI_Get_count,(&status, MPI_CHAR, &cnt));
              buffer.resize(cnt);
              QHIPSTER_MPI_CHECK_RESULT(
                  MPI_Recv,(&buffer[0], cnt, MPI_BYTE, i, i, comm, MPI_STATUS_IGNORE));
              QHIPSTER_MPI_CHECK_RESULT(
                  MPI_Recv,(&sender_state_id, 1, MPI_INT, i, i, comm, MPI_STATUS_IGNORE));
              printf("[|s%d>:%3d] %s\n", sender_state_id, i, (char*)(&buffer[0]));
              fflush(stdout);
          }
      }
      else
      {
          QHIPSTER_MPI_CHECK_RESULT(MPI_Send,(const_cast<char*>(s.c_str()), s.size() + 1,
                                              MPI_CHAR, 0, rank, comm));
          QHIPSTER_MPI_CHECK_RESULT(MPI_Send,(&sender_state_id, 1, MPI_INT, 0, rank, comm));
      }
      MPI_Barrier(comm);
  } 
  else
  {
      MPI_Barrier(comm);
      if (rank == 0) printf("%s\n", s.c_str());
      MPI_Barrier(comm);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// for the MPI implementation

/////////////////////////////////////////////////////////////////////////////////////////

}	// end namespace mpi
}	// end namespace qhipster

/////////////////////////////////////////////////////////////////////////////////////////
