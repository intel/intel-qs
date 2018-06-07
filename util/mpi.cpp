//ff(/ Copyright (C) 2016 Theoretical Physics, ETHZ Zurich

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/** @file mpi.cpp
 *
 *  This header implements MPI support functionality
 */


#if !defined(STANDALONE)
#include "openqu/config.hpp"
#include "openqu/util/mpi.hpp"
#else
#include "mpi.hpp"
//#include "openmp.hpp"
#include "bitops.hpp"
#include "conversion.hpp"
#endif

#ifdef INTELQS_HAS_MPI
#include "mpi.h"
#include <stdexcept>
#include <vector>
#endif

#include <iostream>
#include <iomanip>
#include <thread>


#ifdef INTELQS_HAS_MPI
openqu::mpi::Environment::Environment(int& argc, char**& argv) : inited_(false)
{
  int flag;
  OPENQU_MPI_CHECK_RESULT(MPI_Initialized(&flag))
  if (!flag) {
    OPENQU_MPI_CHECK_RESULT(MPI_Init(&argc, &argv))
    inited_ = true;
#if defined(MVAPICH2_VERSION) 
    char * mv2_string; 
    int mv2_affinity = 1; /* this is the default behavior of MVAPICH2 */ 
    if ((mv2_string = getenv("MV2_ENABLE_AFFINITY")) != NULL) { 
      mv2_affinity = atoi(mv2_string); 
    } 
    if (mv2_affinity!=0 /* && procid==0 FIXME */) { 
      printf("WARNING: You are using MVAPICH2 with affinity enabled, probably by default. \n"); 
      printf("WARNING: This will cause performance issues for ARMCI. \n"); 
      printf("WARNING: Please rerun your job with MV2_ENABLE_AFFINITY=0 \n"); 
    } 
#endif 
  }

  int useful_ranks_per_node = 0;

  int myrank = rank();
  int nranks = size();
  int mynoderank, nrankspernode;
  MPI_Comm nodeComm;
  MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, myrank,
                       MPI_INFO_NULL, &nodeComm );
  MPI_Comm_rank(nodeComm, &mynoderank);
  MPI_Comm_size(nodeComm, &nrankspernode);

  assert((nranks % nrankspernode) == 0);
  mynodeid = myrank / nrankspernode;
  nnodes   = nranks / nrankspernode;

  useful_rank = true;
  int usefull_nranks = openqu::floor_power_of_two(nranks);

  if (openqu::isPowerOf2(nranks) == false) {
    // fun starts
    MPI_Group  orig_group, new_group;
    MPI_Comm   new_comm;
  

    std::vector<int> usefull_ranks, dummy_ranks;
    for(int i = 0; i < usefull_nranks; i++) usefull_ranks.push_back(i);
    for(int i = usefull_nranks; i < nranks; i++) dummy_ranks.push_back(i);

    useful_rank = (myrank < usefull_ranks.size());

    /* Extract the original group handle */
    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

    if (useful_rank) {
      MPI_Group_incl(orig_group, usefull_ranks.size(), usefull_ranks.data(), &new_group);
    } else {
      MPI_Group_incl(orig_group, dummy_ranks.size(), dummy_ranks.data(), &new_group);
    }
    /* Create new new communicator and then perform collective communications */
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);


    if (useful_rank) {
      communicator = new_comm;
    }
  } else {
  }

#if 0
#if 0
  int nuseful_ranks_per_node = 10000000;
  int first_rank_on_last_node = nranks - nrankspernode;
  if(useful_rank) {
    if (myrank < first_rank_on_last_node)
      nuseful_ranks_per_node = nrankspernode;
    else
      nuseful_ranks_per_node = (usefull_nranks - (nranks - nrankspernode));
  }

  int threads_per_rank = (std::thread::hardware_concurrency() / 2) / nuseful_ranks_per_node;
#else
  int threads_per_rank = (std::thread::hardware_concurrency() / 2) / nrankspernode;
#endif
  assert(threads_per_rank == openqu::openmp::omp_get_set_num_threads());
#endif

  

//  std::string aff = openqu::openmp::init(0);
//  int threads_per_rank = openqu::openmp::omp_get_set_num_threads();

  glb_affinity.set_thread_affinity(1);
  int threads_per_rank = glb_affinity.get_num_threads();
  std::string aff_str = glb_affinity.get_affinity_string();

  std::stringstream buffer;
  buffer    <<  "glbrank: " << std::setw(6) << openqu::toString(myrank) 
            <<  " mynodeid: " << std::setw(6) << openqu::toString(mynodeid) 
            <<  " nnodes: " << std::setw(6) << openqu::toString(nnodes) 
            <<  " ranks/node: " << std::setw(6) << openqu::toString(nrankspernode) 
            <<  " threads/rank: " << std::setw(2) << openqu::toString(threads_per_rank) 
            <<  " aff: " << std::setw(50) << aff_str
            <<  (useful_rank ? " --useful" : " --dummy");
  openqu::mpi::print(buffer.str(), MPI_COMM_WORLD);
  

  // start synching all threads
  // MPI_Ibarrier(MPI_COMM_WORLD, &synch_request);

}
#else
openqu::mpi::Environment::Environment(int&, char**&) {}
#endif

openqu::mpi::Environment::~Environment()
{
#ifdef INTELQS_HAS_MPI
  if (inited_) { 
#if 0
    if (is_usefull_rank()) {
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
#endif
}

int openqu::mpi::Environment::rank()
{
#ifdef INTELQS_HAS_MPI
  int rank;
  OPENQU_MPI_CHECK_RESULT((MPI_Comm_rank(comm(), &rank)))
  return rank;
#else
  return 0;
#endif
}

int openqu::mpi::Environment::size()
{
#ifdef INTELQS_HAS_MPI
  int size;
  OPENQU_MPI_CHECK_RESULT((MPI_Comm_size(comm(), &size)))
  return size;
#else
  return 1;
#endif
}

int openqu::mpi::Environment::get_nrankspernode()
{
#ifdef INTELQS_HAS_MPI
  MPI_Comm nodeComm;
  int nrankspernode;
  int myrank = rank();
  MPI_Comm_split_type( comm(), MPI_COMM_TYPE_SHARED, myrank,
                       MPI_INFO_NULL, &nodeComm );
  // MPI_Comm_rank(nodeComm, &mynoderank);
  MPI_Comm_size(nodeComm, &nrankspernode);

  return nrankspernode;
#else
  return 1;
#endif
}

int openqu::mpi::Environment::get_nnodes()
{
   return nnodes;
}

int openqu::mpi::Environment::get_nodeid()
{
   return mynodeid;
}



void openqu::mpi::Environment::remaprank(int newme)
{
#ifdef INTELQS_HAS_MPI
  int r=rank();
  MPI_Comm oldcomm = communicator, newcomm;
  MPI_Comm_split(oldcomm, 0, newme, &newcomm);
  // printf("rank=%d changing %d to %d oldcomm=%d newcom=%d\n", r, newme, oldcomm, oldcomm, newcomm);
  communicator = newcomm;
#endif
}


#ifdef INTELQS_HAS_MPI
MPI_Comm openqu::mpi::Environment::comm()
{
  return communicator;
}
#endif


#ifdef INTELQS_HAS_MPI

openqu::mpi::Exception::Exception(int result) : result_(result)
{
  // Get the error messaage from MPI
  char txt[MPI_MAX_ERROR_STRING];
  int len;
  MPI_Error_string(result, txt, &len);
  text_ = txt;
}

openqu::mpi::Exception::~Exception() throw() {}

#endif

void openqu::mpi::barrier()
{
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm(); 
  MPI_Barrier(comm);
#endif
}

void openqu::mpi::print(std::string s, bool all)
{
  int rank = Environment::rank();
  int size = Environment::size();

  if (all) {
    barrier();
    if (rank == 0) {printf("[%3d] %s\n", rank, s.c_str()); fflush(stdout);}
#ifdef INTELQS_HAS_MPI
    MPI_Comm comm = openqu::mpi::Environment::comm(); 
    std::vector<char> buffer;
    if (rank == 0) {
      for (int i = 1; i < size; i++) {
        MPI_Status status;
        OPENQU_MPI_CHECK_RESULT((MPI_Probe(i, i, comm, &status)));
        int cnt = 0;
        OPENQU_MPI_CHECK_RESULT((MPI_Get_count(&status, MPI_CHAR, &cnt)));
        buffer.resize(cnt);
        OPENQU_MPI_CHECK_RESULT(
            (MPI_Recv(&buffer[0], cnt, MPI_BYTE, i, i, comm, MPI_STATUS_IGNORE)));
        printf("[%3d] %s\n", i, (char*)(&buffer[0])); 
        fflush(stdout);
      }
    } else
      OPENQU_MPI_CHECK_RESULT((MPI_Send(const_cast<char*>(s.c_str()), s.size() + 1, MPI_CHAR, 0,
                                        rank, comm)));
#endif  // INTELQS_HAS_MPI
  } else {
    barrier();
    if (rank == 0) std::cout << s << std::endl;
    barrier();
  }
}

void openqu::mpi::print(std::string s, MPI_Comm comm)
{
  int rank ;
  int size ;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (1) {
    if (rank == 0) printf("[%3d] %s\n", rank, s.c_str());
#ifdef INTELQS_HAS_MPI
    std::vector<char> buffer;
    if (rank == 0) {
      for (int i = 1; i < size; i++) {
        MPI_Status status;
        OPENQU_MPI_CHECK_RESULT((MPI_Probe(i, i, comm, &status)));
        int cnt = 0;
        OPENQU_MPI_CHECK_RESULT((MPI_Get_count(&status, MPI_CHAR, &cnt)));
        buffer.resize(cnt);
        OPENQU_MPI_CHECK_RESULT(
            (MPI_Recv(&buffer[0], cnt, MPI_BYTE, i, i, comm, MPI_STATUS_IGNORE)));
        printf("[%3d] %s\n", i, (char*)(&buffer[0]));
      }
    } else
      OPENQU_MPI_CHECK_RESULT((MPI_Send(const_cast<char*>(s.c_str()), s.size() + 1, MPI_CHAR, 0,
                                        rank, comm)));
#endif  // INTELQS_HAS_MPI
  } 
}


#ifdef INTELQS_HAS_MPI
MPI_Comm openqu::mpi::Environment::communicator = MPI_COMM_WORLD;
#endif
unsigned openqu::mpi::Environment::nnodes = 1;
unsigned openqu::mpi::Environment::mynodeid = 1;
