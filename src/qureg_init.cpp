#ifdef USE_MKL
#include "mkl_vsl.h"
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif

#include "../include/qureg.hpp"
/// \addtogroup qureg
/// @{

/// @file qureg_init.cpp
/// @brief Define the @c QubitRegister methods to initialize the quantum register.

/////////////////////////////////////////////////////////////////////////////////////////
/// Base constructor.
template <class Type>
QubitRegister<Type>::QubitRegister()
{
  unsigned myrank=0, nprocs=1;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();

  timer = nullptr;
  gate_counter = nullptr;
  qubit_permutation = nullptr;
  imported_state = false;
  specialize = false;
  num_qubits = 1;

  fusion = false;

  Resize(1UL);
  state_storage[0] = {1., 0.};

  if (nprocs > 1) {
    fprintf(stderr,
            "nprocs > 1: seperate tmp storage from state vector, or some routines won't work\n");
    assert(0);
  }
  timer = nullptr;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Resize(std::size_t new_num_amplitudes)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  log2_nprocs = qhipster::ilog2(nprocs);

  // FIXME GG: I believe this limits the use of "resize" to adding a single qubit
  if(GlobalSize()) assert(GlobalSize() * 2UL == new_num_amplitudes);
  num_qubits = qhipster::ilog2(new_num_amplitudes);

  local_size_  = UL(1L << UL(num_qubits - log2_nprocs));
  global_size_ = UL(1L << UL(num_qubits));
  assert(LocalSize() >= 1L);

  std::size_t num_amplitudes = (nprocs == 1) ? LocalSize() : (LocalSize() + TmpSize());
  std::size_t nbytes = num_amplitudes * sizeof(state[0]);

#if defined(USE_MM_MALLOC)
  state = (Type *)_mm_malloc(nbytes, 256);
#else
  state_storage.resize(num_amplitudes);
  state = &state_storage[0];
#endif

  if (qubit_permutation) delete qubit_permutation;
  qubit_permutation = new Permutation(num_qubits);
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::AllocateAdditionalQubit()
{
  ++num_qubits;
  Resize( UL(1) << UL(num_qubits) );
}


/////////////////////////////////////////////////////////////////////////////////////////
// The second argument is not used nor returned modified.
// 
template <class Type>
void QubitRegister<Type>::Initialize(std::size_t new_num_qubits, std::size_t tmp_spacesize_)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0, num_ranks_per_node=1;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  log2_nprocs = qhipster::ilog2(nprocs);
  assert(new_num_qubits>log2_nprocs && "Too few qubits for this number of ranks");
  num_ranks_per_node = qhipster::mpi::Environment::GetNumRanksPerNode();
  unsigned M = new_num_qubits - log2_nprocs;

  assert(new_num_qubits > 0);
  local_size_  = UL(1L << UL(new_num_qubits - log2_nprocs));
  global_size_ = UL(1L << UL(new_num_qubits));

  assert(LocalSize()>1);	// Check in case one used too many ranks.
  std::size_t lcl_size_half = LocalSize() / 2L;

  // tmp_spacesize_ is used to determine the amount of memory allocated for the communication
  // scheme during the distributed implementation of quantum operations.
  // The case where tmp_spacesize_==lcl_size_half is fully implemented and efficient.
  // Sometimes it is useful to reduce the tmp_spacesize_, to be able to simulate
  // one extra qubit. The implementation of SWAP-like gates is not ready for this yet.

  size_t hard_bound_tmp_spacesize = UL(1L << UL(30));  // 4194304 = 2^22
  if (    tmp_spacesize_ == 0		// default case
       || local_size_ <= tmp_spacesize_)	// to avoid waste of memory
  {
      // if (!myrank) printf("Setting tmp storage to half the local state size\n");
      this->tmp_spacesize_ =  lcl_size_half;
  }
  else if (tmp_spacesize_ <= hard_bound_tmp_spacesize)
  {
      assert((lcl_size_half % tmp_spacesize_) == 0);
      // if (!myrank) printf("Setting tmp storage to less than half the local state size, specifically to: %lu\n",tmp_spacesize_);
      this->tmp_spacesize_ =  tmp_spacesize_;
  }
  else
      this->tmp_spacesize_ = hard_bound_tmp_spacesize;

  this->num_qubits = new_num_qubits;
  assert(LocalSize() >= 1L);

  // Set-up initial qubit permutation and state_rank permutation.
  qubit_permutation = new Permutation(new_num_qubits);

  if ( do_print_extra_info && !myrank)
      printf("Specialization is off\n");

  timer = nullptr;
  gate_counter = nullptr;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Allocate(std::size_t new_num_qubits, std::size_t tmp_spacesize_)
{
  unsigned myrank=0, nprocs=1, num_ranks_per_node=1;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  num_ranks_per_node = qhipster::mpi::Environment::GetNumRanksPerNode();

  imported_state = false;
  specialize = false;
  fusion = false;

  Initialize(new_num_qubits, tmp_spacesize_);

  std::size_t num_amplitudes = (nprocs == 1) ? LocalSize() : (LocalSize() + TmpSize());
  std::size_t nbytes = num_amplitudes * sizeof(state[0]);

  // Print some information.
  if (do_print_extra_info && !myrank)
  {
      double MB = 1024.0 * 1024.0;
      double s;
      s = D(num_ranks_per_node) * D(nbytes);
      printf("Total storage per node  = %.2lf MB \n", s / MB);
      s = D(num_ranks_per_node) * D(LocalSize()) * D(sizeof(state[0]));
      printf("      storage per state = %.2lf MB \n", s / MB);
      if (nprocs > 1)
      {
          s = D(num_ranks_per_node) * D(TmpSize()) * D(sizeof(state[0]));
          printf("      temporary storage = %.5lf MB \n", s / MB);
      }
  }

#if defined(USE_MM_MALLOC)
  state = (Type *)_mm_malloc(nbytes, 256);
#else
  state_storage.resize(num_amplitudes);
  state = &state_storage[0];
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/// Constructor <TODO: add extra info>
template <class Type>
QubitRegister<Type>::QubitRegister(std::size_t new_num_qubits, Type *state, 
                                   std::size_t tmp_spacesize_)
{
  imported_state = true;
  Initialize(new_num_qubits, tmp_spacesize_);
  this->state = state;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// Constructor, followed by state initialization.
template <class Type>
QubitRegister<Type>::QubitRegister(std::size_t new_num_qubits, 
                                   std::string style, 
                                   std::size_t base_index,
                                   std::size_t tmp_spacesize_)
{
  Allocate(new_num_qubits, tmp_spacesize_);
  Initialize(style, base_index);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the state.
///
/// The 'style' of initialization can be:
/// - 'rand': real and imag part of each amplitudes are uniformly random,
///           using either the **local** or **pool** RNG stream,
///           then state is normalized.
/// - 'base': state of the computational basis, only a non-zero amplitude.
/// - '++++': the balanced superposition of all computational basis states.
template <class Type>
void QubitRegister<Type>::Initialize(std::string style, std::size_t base_index)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  nprocs = qhipster::mpi::Environment::GetStateSize();
  log2_nprocs = qhipster::ilog2(nprocs);
  unsigned nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  {
      nthreads = omp_get_num_threads();
  }
#endif

  double t0 = time_in_seconds();

  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < lcl; i++)
      state[i] = {0, 0};

  ///////////////////////////////////////////////////////////////////////////////////////
  // Random state:
  // real and imag part of each amplitude is uniformly random,
  // using the **local** stream associated with the seed,
  // then state is normalized.
  if (style == "rand")
  {
      // The random number generator pointed by rng_ptr_ is used.
      // It must have been initialized earlier.
      assert(rng_ptr_!=nullptr);
      // The RNG must have been initialized with a seed different from 0.
      assert(rng_ptr_->GetSeed()!=0);

      // The value of base_index determines whether the random state is common among
      // all states of the pool, or unique for each distinct state_id.
      // base_index=0         : common state, use pool RNG stream
      // base_index=num_states: individual random states, use local RNG stream
      //                        (and reduced fast-forward)
      // Any other values return error.
      assert(base_index==0 || base_index==qhipster::mpi::Environment::GetNumStates() );

      // Parallel initialization using open-source parallel RNG or VRL (if MKL is used).
      // TODO: with GCC, using OpenMP code below produces a SEGMENTATION FAULT result.
      //       This happens when randomly initializing states with 20 or more qubits (no MPI)
      //       Currently we reserve the OpenMP initialization to ICPC only.
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
          std::size_t thread_id   = omp_get_thread_num();
          std::size_t num_threads = omp_get_num_threads();
#else
          std::size_t thread_id = 0;
          std::size_t num_threads = 1;
#endif
          std::size_t chunk = LocalSize() / num_threads;
          std::size_t beginning = thread_id * chunk;
          std::size_t end = (thread_id + 1) * chunk;
          if (thread_id == num_threads - 1)
              end = LocalSize();

          // Since the threads are not executed in order, we copy the rng stream per thread.
          // at this point every copy is independent of each other and from the original.
          qhipster::RandomNumberGenerator<BaseType> thread_rng ( rng_ptr_);

          // Fast forward for the thread:
          std::size_t num_skip = 2UL * beginning;

          if (base_index>0)
          {
              thread_rng.SkipAhead(num_skip,"local");
              // Directly generate numbers in the correct range of the state storage.
//FIXME
#if 0
std::stringstream buffer;
buffer << "random init, thread_id = " << thread_id << " , from " << beginning
       << " to " << end << "\n";//FIXME
if (beginning != end) printf( (buffer.str()).c_str() );
#endif
              thread_rng.UniformRandomNumbers( (BaseType *)&(state[beginning]),
                                              2UL*(end-beginning), -1., 1., "local");
          }
          else
          {
              // If the RNG stream is the pool one, then one has to skip many more numbers:
              num_skip += 2UL * myrank * LocalSize();
              thread_rng.SkipAhead(num_skip,"pool");
              // Directly generate numbers in the correct range of the state storage.
              thread_rng.UniformRandomNumbers( (BaseType *)&(state[beginning]),
                                              2UL*(end-beginning), -1., 1., "pool");
          }
      }
      // Update the main rng_prt_ by skipping the numbers already used.
      if (base_index>0)
          rng_ptr_->SkipAhead(2UL*LocalSize(),"local");
      else
          rng_ptr_->SkipAhead(2UL*GlobalSize(),"pool" );

      // Normalize the state.
      this->Normalize();
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  // Computational basis state.
  else if (style == "base")
  {
      assert(base_index < GlobalSize());

      std::size_t whereid = base_index / LocalSize();
      if (whereid == myrank)
      {
          std::size_t lclind = (base_index % LocalSize());
          state[lclind] = {1.0, 0.0};
      }
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  // Balanced superposition of all basis states.
  else if (style == "++++")
  {
      Type amplitude = {1./std::sqrt( GlobalSize() ),0.};
      std::size_t lcl = LocalSize();
#pragma omp parallel for
      for (std::size_t i = 0; i < lcl; i++)
          state[i] = amplitude;
  }

#ifdef INTELQS_HAS_MPI
  qhipster::mpi::StateBarrier();
#endif

#if 0
  double t1 = time_in_seconds();
  if (myrank == 0) {
    printf("[%u] Time to init: %lf\n", myrank, t1 - t0);
  }
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/// Constructor that copies another @c QubitRegister object.
template <class Type>
QubitRegister<Type>::QubitRegister(const QubitRegister &in)
{
  Allocate(in.num_qubits, in.TmpSize());
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
  TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < lcl; i++)
      state[i] = in.state[i];
 
  *qubit_permutation = *(in.qubit_permutation);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Specialization using the unitary matrix structure.
///
/// When turned on, avoids full matrix multiplication in some special cases
/// to improve performance.
/// Turned off by default.
template <class Type>
void QubitRegister<Type>::TurnOnSpecialize()
{
  int myrank=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  if (do_print_extra_info && !myrank)
      printf("Specialization is on\n");
  specialize = true;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::TurnOffSpecialize()
{
  unsigned myrank=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  if (do_print_extra_info && !myrank)
      printf("Specialization is off\n");
  specialize = false;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Specialization using the executed gate types.
///
/// Avoids matrix multiplication in some common gates to improve performance.
/// Turned off by default.
///
/// Supported Gate Types:
///     - 1-Qubit Gates
///         - Hadamard
///         - Rotation(X, Y, Z)
///         - Pauli(X, Y, Z)
///         - T
///     - Controlled Gates
///         - CHadamard
///         - CRotation(X, Y, Z)
///         - CPauli(X, Y, Z)
///         - CPhaseRotation
///
/// @warning May not work with gate fusion!
template <class Type>
void QubitRegister<Type>::TurnOnSpecializeV2()
{
  int myrank=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  if (do_print_extra_info && !myrank)
    printf("Specialization v2 is on\n");
  specialize2 = true;
}

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::TurnOffSpecializeV2()
{
  unsigned myrank=0;
  myrank = qhipster::mpi::Environment::GetStateRank();
  if (do_print_extra_info && !myrank)
    printf("Specialization v2 is off\n");
  specialize2 = false;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Destructor.
template <class Type>
QubitRegister<Type>::~QubitRegister()
{
#ifdef USE_MM_MALLOC
  _mm_free(state); 
#endif
  if (timer != nullptr) delete timer;
  if (gate_counter != nullptr) delete gate_counter;
  if (qubit_permutation != nullptr) delete qubit_permutation;
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;


/// @}
