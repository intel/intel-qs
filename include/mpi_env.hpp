#ifndef IQS_MPI_ENV_HPP
#define IQS_MPI_ENV_HPP

/// \addtogroup util
/// @{

/// @file mpi_env.hpp
///
/// This header includes the MPI library if available and provides some
/// convenience functions that are implemented even if the MPI library is not provided.

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

#include <string>
#include <unistd.h>
#include <stdexcept>

#ifdef INTELQS_HAS_MPI
#include "mpi_exception.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////

namespace qhipster {

namespace mpi {

/////////////////////////////////////////////////////////////////////////////////////////
/// \class Environment
/// A trimmed down version of the BOOST::MPI environment.  Its purpose is to initialize
/// the MPI library and partition the cluster or single threaded environment for parallel
/// operations.
/// In preparation of the Monte-Carlo simulations required for the noisy implementation,
/// we provide a communicator over the ranks involved in a single MC simulation.
///
/// Specifically, we store two communicators:
/// 1. pool_communicator: spanning all the useful ranks
/// 2. state_communicator: spanning those ranks used in a single MC simulation (i.e. state)

class Environment
{
  public:
  /// Intialize the MPI Environment.
  ///
  /// It receives the same argc and argv arguments passed to the main function.
  /// If MPI is present, but has not been initialized, then MPI_Init will be called.

  Environment(int& argc, char**& argv);
  Environment();

  /// Finalize the MPI Environment
  ///
  /// If MPI is present and has been initizlized in the constructor then
  /// MPI_Finalize will be called here.

  ~Environment();
  Environment(Environment const&) = delete;
  Environment& operator=(Environment const&) = delete;

  /////////////////////////////////////
  // One needs to know what he/she is doing to call this method.
  /// Update the state and pool communicators.
  /// @pre This can only be called when all ranks are still active.
  static void UpdateStateComm (int num_states, bool do_print_info=true);
  /////////////////////////////////////

  /// Check whether the rank is useful or not.
  static bool IsUsefulRank() {return useful_rank;}

  /// The rank of the current MPI process: pool or state.
  ///
  /// The PoolRank may not corresponds to that from MPI_COMM_WORLD due to dummy ranks.
  /// The rank is 0 if MPI is not present.
  /// @pre If MPI is present, this can only be called after intializing MPI.
  static int GetPoolRank();
  static int GetStateRank();
  static int GetRank() {return GetStateRank();}

  /// Number of MPI processes.
  ///
  /// The PoolSize may not corresponds to that from MPI_COMM_WORLD due to dummy ranks.
  /// The number of processes is 1 if MPI is not present.
  /// @pre If MPI is present, this can only be called after intializing MPI.
  static int GetPoolSize();
  static int GetStateSize();
  static int GetSize() {return GetStateSize();}

  /// Get incoherent average over all states of the pool.
  ///
  /// @param local_value the address of the value stored in the local rank.
  template <class Type>
  static Type IncoherentSumOverAllStatesOfPool(Type local_value);

  static int GetNumRanksPerNode();
  static int GetNumNodes();
  static int GetNodeId();
  static int GetStateId();
  static int GetNumStates();
  static void RemapStateRank(int newme);

#ifdef INTELQS_HAS_MPI
  static MPI_Comm pool_communicator;
  static MPI_Comm state_communicator;
  static MPI_Comm GetPoolComm();
  static MPI_Comm GetStateComm();
  static MPI_Comm GetComm();
#endif


/* 
The singleton interface that is used for global management
of MPI resources.
*/
  static void Init();
  static void Init(int &argc, char**&argv);
  static void Finalize();
  static Environment *GetSharedInstance() {return shared_instance;}

 private:
// Shared helper for constructor overload
#ifdef INTELQS_HAS_MPI
  void CommonInit(int flag);
#endif

  static Environment *shared_instance;

  bool inited_;
  static bool useful_rank;
//  int num_ranks_per_state;

  static int num_ranks_per_node;

  static int num_nodes;
  static int my_node_id;
  static int num_states;
  static int my_state_id;

#ifdef INTELQS_HAS_MPI
  MPI_Request synch_request;
#endif
};

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief An MPI barrier
///
/// It waits until all MPI processes have reached this call. This function does nothing if
/// MPI is not present.
void PoolBarrier();
void StateBarrier();
void Barrier();

/// @brief Find minimum wall clock time in all processes and return it at rank 0.
/// Returns local time for other processes.
double MinTime();

/// @brief Find maximum wall clock time in all processes and return it at rank 0.
/// Returns local time for other processes.
double MaxTime();

/////////////////////////////////////////////////////////////////////////////////////////

/// @brief Print from all MPI processes
///
/// It prints a string from all processes if @c all is true or just from the master process
///with rank 0 if @c all is false. If MPI is not present it prints the string.
///
/// If @c all is set, the string is prefixed by the number of the MPI process.
///
/// @param s the string to be printed
/// @param all a flag to specify if all processes should print or just the master process
void PoolPrint(std::string s, bool all=false);
void StatePrint(std::string s, bool all=false);
void Print(std::string s, bool all=false);
#ifdef INTELQS_HAS_MPI
void Print(std::string s, MPI_Comm comm, bool all=true);
#endif



/////////////////////////////////////////////////////////////////////////////////////////

}	// end namespace mpi
}	// end namespace qhipster

/// @}*/

#endif	// header guard IQS_MPI_ENV_HPP
