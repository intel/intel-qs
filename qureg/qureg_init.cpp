//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

#include "qureg.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_init.cpp
/// @brief Define the @c QubitRegister methods to initialize the quantum register.

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
QubitRegister<Type>::QubitRegister()
{
  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
#endif

  timer = NULL;
  permutation = NULL;
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
  timer = NULL;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Resize(std::size_t new_num_amplitudes)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  log2_nprocs = openqu::ilog2(nprocs);
#endif

  // FIXME GG: I believe this limits the use of "resize" to adding a single qubit
  if(GlobalSize()) assert(GlobalSize() * 2UL == new_num_amplitudes);
  num_qubits = openqu::ilog2(new_num_amplitudes);

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

  TODO(move permutation to WaveFunctionSimulator class)
  if (permutation) delete permutation;
  permutation = new Permutation(num_qubits);
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::AllocateAdditionalQubit()
{
  ++num_qubits;
  Resize( UL(1) << UL(num_qubits) );
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Initialize(std::size_t new_num_qubits, std::size_t tmp_spacesize_)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0, num_ranks_per_node=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  log2_nprocs = openqu::ilog2(nprocs);
  num_ranks_per_node = openqu::mpi::Environment::get_nrankspernode();
#endif
  unsigned M = new_num_qubits - log2_nprocs;

  assert(new_num_qubits > 0);
  local_size_  = UL(1L << UL(new_num_qubits - log2_nprocs));
  global_size_ = UL(1L << UL(new_num_qubits));

  std::size_t lcl_size_half = LocalSize() / 2L;

  #if 0
  if (tmp_spacesize_ == 0 || local_size_ < tmp_spacesize_ )
  {
      if (!myrank) printf("Setting tmp storage to half the local state size\n");
      this->tmp_spacesize_ =  lcl_size_half;
  }
  else
  {
      this->tmp_spacesize_ =  tmp_spacesize_;
      assert((lcl_size_half % tmpSize()) == 0);
  }
  #else
  if (num_ranks_per_node <= 2)
      this->tmp_spacesize_ =  lcl_size_half;
  else
      this->tmp_spacesize_ =  (lcl_size_half > 4194304) ? 4194304 : lcl_size_half;
  assert((lcl_size_half % TmpSize()) == 0);
  #endif

  this->num_qubits = new_num_qubits;
  assert(LocalSize() >= 1L);

  // set-up initial permutation
  permutation = new Permutation(new_num_qubits);

  if (!myrank) printf("Specialization is off\n");

  timer = NULL;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::Allocate(std::size_t new_num_qubits, std::size_t tmp_spacesize_)
{
  unsigned myrank=0, nprocs=1, num_ranks_per_node=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  num_ranks_per_node = openqu::mpi::Environment::get_nrankspernode();
#endif

  imported_state = false;
  specialize = false;
  fusion = false;

  Initialize(new_num_qubits, tmp_spacesize_);

  std::size_t num_amplitudes = (nprocs == 1) ? LocalSize() : (LocalSize() + TmpSize());
  std::size_t nbytes = num_amplitudes * sizeof(state[0]);

  // print some information
  if (!myrank)
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
  state_storage.Resize(num_amplitudes);
  state = &state_storage[0];
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
QubitRegister<Type>::QubitRegister(std::size_t new_num_qubits, Type *state, 
                                   std::size_t tmp_spacesize_)
{
  imported_state = true;
  Initialize(new_num_qubits, tmp_spacesize_);
  this->state = state;
}


/////////////////////////////////////////////////////////////////////////////////////////
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
template <class Type>
void QubitRegister<Type>::Initialize(std::string style, std::size_t base_index)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  log2_nprocs = openqu::ilog2(nprocs);
  MPI_Comm comm = openqu::mpi::Environment::comm();
#endif
  unsigned nthreads = glb_affinity.get_num_threads();

  double t0 = time_in_seconds();

  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++) state[i] = {0, 0};

//// each amplitude is random , then state is normalized ////////////////////////////////
  if (style == "rand")
  {
      BaseType local_normsq = 0;

#if defined(SEQUENTIAL_INITIALIZATION)
      // sequential
      srand(base_index);
      // no need of fast-forward rng
      for (std::size_t i = 0; i < LocalSize(); i++)
      {
          state[i] = {RAND01(), RAND01()};
          local_normsq += std::abs(state[i]) * std::abs(state[i]);
      }
#endif

#if (defined(__ICC) || defined(__INTEL_COMPILER))
// --------------------- FIXME by Gian: moved to separate method with template specialization
      RandomInitialize(base_index);
#elif 0 
// --------------------- FIXME by Gian: MIT prng_engine excluded from the choices
  // Parallel initialization using open source parallel RNG
//  std::vector<sitmo::prng_engine> eng(openqu::openmp::omp_get_set_num_threads());
#pragma omp parallel reduction(+ : local_normsq)
      {
#ifdef _OPENMP
          std::size_t thread_id   = omp_get_thread_num();
          std::size_t num_threads = omp_get_num_threads();
#else
          std::size_t thread_id = 0;
          std::size_t num_threads = 1;
#endif
          std::size_t chunk = LocalSize() / num_threads;
          std::size_t beginning = thread_id * chunk, end = (thread_id + 1) * chunk;
          if (thread_id == num_threads - 1) end = LocalSize();
          // fast forward
          eng[thread_id].discard( 2 * myrank * LocalSize() + 2 * beginning );
#pragma simd reduction(+ : local_normsq)
          for (std::size_t i = beginning; i < end; i++)
          {
              BaseType r1 = (BaseType)eng[thread_id]() / D(UINT_MAX);
              BaseType r2 = (BaseType)eng[thread_id]() / D(UINT_MAX);
              state[i] = {r1, r2};
              // std::cout << "i: " << i << " state: " << state[i];
          }
      }
#else
      std::cout << " ~~~~~~~~~~~~~~~~ no random number generator !! ~~~~~~~~~~~~~~ \n";
#endif

      std::size_t lcl = LocalSize();
#pragma omp parallel for reduction(+ : local_normsq)
      for (std::size_t i = 0; i < lcl; i++)
      {
          local_normsq += std::norm(state[i]) ;
      }

      BaseType global_normsq;
#ifdef INTELQS_HAS_MPI
      // MPI_Allreduce(&local_normsq, &global_normsq, 1, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce_x(&local_normsq, &global_normsq,  MPI_SUM, comm);
#else
      global_normsq = local_normsq;
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
//   #pragma omp parallel for simd
#else
TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#endif
      for (std::size_t i = 0; i < lcl; i++)
      {
          state[i] = state[i] / std::sqrt(global_normsq);
      }

  }
//// computational basis state //////////////////////////////////////////////////////////
  else if (style == "base")
  {
      std::size_t whereid = base_index / LocalSize();
      if (whereid == myrank)
      {
          std::size_t lclind = (base_index % LocalSize());
          state[lclind] = {1.0, 0.0};
      }
  }
//// balanced superposition of all classical bitstrings /////////////////////////////////
  else if (style == "++++")
  {
      state[0] = {1./std::sqrt( GlobalSize() ),0.};
      for (std::size_t i = 1; i < LocalSize(); i++)
      {
          state[i] = state[0];
      }
  }

#ifdef INTELQS_HAS_MPI
  openqu::mpi::barrier();
#endif

#if 0
  double t1 = time_in_seconds();
  if (myrank == 0) {
    printf("[%u] Time to init: %lf\n", myrank, t1 - t0);
  }
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
// template specialization depending on the type of Type
template <typename Type>
void QubitRegister<Type>::RandomInitialize(std::size_t base_index)
{
  std::cout << " ~~~~~~~~~~~~~~~~ wrong type for state! ~~~~~~~~~~~~~~ \n";
}
//--
template <>
void QubitRegister<ComplexDP>::RandomInitialize(std::size_t base_index)
{
    unsigned myrank=0;
#ifdef INTELQS_HAS_MPI
    myrank = openqu::mpi::Environment::rank();
#endif

    // Parallel initialization using parallel MKL RNG
#pragma omp parallel
    {
#ifdef _OPENMP
      std::size_t thread_id  = omp_get_thread_num();
      std::size_t num_threads= omp_get_num_threads();
#else
      std::size_t thread_id = 0;
      std::size_t num_threads = 1;
#endif
      VSLStreamStatePtr stream;
      std::size_t chunk = LocalSize() / num_threads;
      std::size_t beginning = thread_id * chunk, end = (thread_id + 1) * chunk;
      if (thread_id == num_threads - 1) end = LocalSize();

      int errcode = vslNewStream(&stream, VSL_BRNG_MCG31, base_index);
      assert(errcode == VSL_STATUS_OK);
      std::size_t num_skip = 2UL * (myrank * LocalSize() + beginning);
      vslSkipAheadStream(stream, num_skip);
      errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2L * (end - beginning),
                             (double *)&state[beginning], 0.0, 1.0);
      assert(errcode == VSL_STATUS_OK);
    }
}
//--
template <>
void QubitRegister<ComplexSP>::RandomInitialize(std::size_t base_index)
{
    unsigned myrank=0;
#ifdef INTELQS_HAS_MPI
    myrank = openqu::mpi::Environment::rank();
#endif

    // Parallel initialization using parallel MKL RNG
#pragma omp parallel
    {
#ifdef _OPENMP
      std::size_t thread_id  = omp_get_thread_num();
      std::size_t num_threads= omp_get_num_threads();
#else
      std::size_t thread_id = 0;
      std::size_t num_threads = 1;
#endif
      VSLStreamStatePtr stream;
      std::size_t chunk = LocalSize() / num_threads;
      std::size_t beginning = thread_id * chunk, end = (thread_id + 1) * chunk;
      if (thread_id == num_threads - 1) end = LocalSize();

      int errcode = vslNewStream(&stream, VSL_BRNG_MCG31, base_index);
      assert(errcode == VSL_STATUS_OK);
      std::size_t num_skip = 2L * (myrank * LocalSize() + beginning);
      vslSkipAheadStream(stream, num_skip);
      errcode = vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2L * (end - beginning),
                        (float *)&state[beginning], 0.0, 1.0);
      assert(errcode == VSL_STATUS_OK);
    }
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
QubitRegister<Type>::QubitRegister(const QubitRegister &in)
{
  Allocate(in.num_qubits, in.TmpSize());
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
  TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#endif

  for (std::size_t i = 0; i < lcl; i++) state[i] = in.state[i];
  *permutation = *(in.permutation);
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::TurnOnSpecialize()
{
  unsigned myrank=0;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
#endif
  if (!myrank) printf("Specialization is on\n");
  specialize = true;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::TurnOffSpecialize()
{
  unsigned myrank=0;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
#endif
  if (!myrank) printf("Specialization is off\n");
  specialize = false;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
QubitRegister<Type>::~QubitRegister()
{
#if defined(USE_MM_MALLOC)
  _mm_free(state); 
#endif
  if (timer) delete timer;
  if (permutation) delete permutation;
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;


/// @}
