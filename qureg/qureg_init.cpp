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


template <class Type>
QbitRegister<Type>::QbitRegister()
{
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();

  timer = NULL;
  permutation = NULL;
  importedstate = false;
  specialize = false;
  nqbits = 1;

  fusion = false;

  resize(1UL);
  statestorage[0] = {1., 0.};

  
  if (nprocs > 1) {
    fprintf(stderr,
            "nprocs > 1: seperate tmp storage from state vector, or some routines won't work\n");
    assert(0);
  }
  timer = NULL;
}

template <class Type>
void QbitRegister<Type>::resize(std::size_t n)
{
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();

  if(globalSize()) assert(globalSize() * 2UL == n);
  nqbits = openqu::ilog2(n);

  localsize_ = UL(1L << UL(nqbits - openqu::ilog2(nprocs)));
  globalsize_ = UL(1L << UL(nqbits));
  assert(localSize() >= 1L);

  std::size_t namplitudes = (nprocs == 1) ? localSize() : (localSize() + tmpSize());
  std::size_t nbytes = namplitudes * sizeof(state[0]);


#if defined(USE_MM_MALLOC)
  state = (Type *)_mm_malloc(nbytes, 256);
#else
  statestorage.resize(namplitudes);
  state = &statestorage[0];
#endif

  TODO(move permutation to WaveFunctionSimulator class)
  if (permutation) delete permutation;
  permutation = new Permutation(nqbits);
  
}

template <class Type>
void QbitRegister<Type>::allocateQubit()
{

  ++nqbits;
  resize(UL(1) << UL(nqbits));
}

template <class Type>
void QbitRegister<Type>::Init(std::size_t nqbits, std::size_t tmpspacesize_)
{
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  unsigned log2_nprocs = openqu::ilog2(nprocs);

  assert(nqbits > openqu::ilog2(nprocs));
  localsize_ = UL(1L << UL(nqbits - openqu::ilog2(nprocs)));
  globalsize_ = UL(1L << UL(nqbits));

  std::size_t lcl_size_half = localSize() / 2L;

  #if 0
  if (tmpspacesize_ == 0 || localsize_ < tmpspacesize_ ) {
    if (!myrank) printf("Setting tmp storage to half the local state size\n");
    this->tmpspacesize_ =  lcl_size_half;
  } else {
    this->tmpspacesize_ =  tmpspacesize_;
    assert((lcl_size_half % tmpSize()) == 0);
  }
  #else
  if (openqu::mpi::Environment::get_nrankspernode() <= 2) {
    this->tmpspacesize_ =  lcl_size_half;
  } else {
    this->tmpspacesize_ =  (lcl_size_half > 4194304) ? 4194304 : lcl_size_half;
  }
  assert((lcl_size_half % tmpSize()) == 0);
  #endif


  this->nqbits = nqbits;
  assert(localSize() >= 1L);

  // set-up initial permutation
  permutation = new Permutation(nqbits);

  if (!myrank) printf("Specialization is off\n");

  timer = NULL;
}

template <class Type>
void QbitRegister<Type>::Allocate(std::size_t nqbits, std::size_t tmpspacesize_)
{
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();

  importedstate = false;
  specialize = false;
  fusion = false;

  Init(nqbits, tmpspacesize_);

  std::size_t namplitudes = (nprocs == 1) ? localSize() : (localSize() + tmpSize());
  std::size_t nbytes = namplitudes * sizeof(state[0]);

  int nrankspernode = openqu::mpi::Environment::get_nrankspernode();

  if (!myrank) {
    double MB = 1024.0 * 1024.0;
    double s;
    s = D(nrankspernode) * D(nbytes);
    printf("Total storage per node  = %.2lf MB \n", s / MB);
    s = D(nrankspernode) * D(localSize()) * D(sizeof(state[0]));
    printf("      storage per state = %.2lf MB \n", s / MB);
    if (nprocs > 1) {
       s = D(nrankspernode) * D(tmpSize()) * D(sizeof(state[0]));
       printf("      temporary storage = %.5lf MB \n", s / MB);
    }
  }

#if defined(USE_MM_MALLOC)
  state = (Type *)_mm_malloc(nbytes, 256);
#else
  statestorage.resize(namplitudes);
  state = &statestorage[0];
#endif


  openqu::mpi::barrier();
}

template <class Type>
QbitRegister<Type>::QbitRegister(std::size_t nqbits, Type *state, 
                                 std::size_t tmpspacesize_)
{
  importedstate = true;
  Init(nqbits, tmpspacesize_);
  this->state = state;
}

template <class Type>
QbitRegister<Type>::QbitRegister(std::size_t nqbits, 
                                 std::string style, 
                                 std::size_t baseind,
                                 std::size_t tmpspacesize_)
{
  Allocate(nqbits, tmpspacesize_);
  Init(style, baseind);
}

template <class Type>
void QbitRegister<Type>::Init(std::string style, std::size_t baseind)
{
  unsigned myrank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  unsigned nthreads = glb_affinity.get_num_threads();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  double t0 = time_in_seconds();

  std::size_t lcl = localSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++) state[i] = {0, 0};

  if (style == "rand") {

    BaseType local_normsq = 0;

#if defined(SEQUENTIAL_INITIALIZATION)
    // sequential
    srand(baseind);
    // fast-forward rng
    for (std::size_t i = 0; i < myrank * localSize(); i++) {
      // intentional
      Type r = {RAND01(), RAND01()};
    }
    for (std::size_t i = 0; i < localSize(); i++) {
      state[i] = {RAND01(), RAND01()};
      local_normsq += abs(state[i]) * abs(state[i]);
    }
#endif

#if (defined(__ICC) || defined(__INTEL_COMPILER))
// --------------------- FIXME by Gian: moved to separate method with template specialization
    util_rand_init(baseind);
#elif 0 
// --------------------- FIXME by Gian: MIT prng_engine excluded form the choices
  // Parallel initialization using open source parallel RNG
//  std::vector<sitmo::prng_engine> eng(openqu::openmp::omp_get_set_num_threads());
#pragma omp parallel reduction(+ : local_normsq)
  {
#ifdef _OPENMP
    std::size_t tid = omp_get_thread_num();
    std::size_t nthreads=omp_get_num_threads();
#else
    std::size_t tid = 0;
    std::size_t nthreads = 1;
#endif
    std::size_t chunk = localSize() / nthreads;
    std::size_t beg = tid * chunk, end = (tid + 1) * chunk;
    eng[tid].discard(2 * myrank * localSize() + 2 * beg);
    if (tid == nthreads - 1) end = localSize();
#pragma simd reduction(+ : local_normsq)
    for (std::size_t i = beg; i < end; i++) {
      BaseType r1 = (BaseType)eng[tid]() / D(UINT_MAX);
      BaseType r2 = (BaseType)eng[tid]() / D(UINT_MAX);
      state[i] = {r1, r2};
      // std::cout << "i: " << i << " state: " << state[i];
    }
  }
#else
  std::cout << " ~~~~~~~~~~~~~~~~ no random number generator !! ~~~~~~~~~~~~~~ \n";
#endif

    std::size_t lcl = localSize();
#pragma omp parallel for reduction(+ : local_normsq)
    for (std::size_t i = 0; i < lcl; i++) {
      local_normsq += abs(state[i]) * abs(state[i]);
    }

    BaseType global_normsq;
#ifdef OPENQU_HAVE_MPI
    // MPI_Allreduce(&local_normsq, &global_normsq, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce_x(&local_normsq, &global_normsq,  MPI_SUM, comm);
#else
    global_normsq = local_normsq;
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
// #pragma omp parallel for simd
#else
TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#endif
    for (std::size_t i = 0; i < lcl; i++) {
      state[i] = state[i] / std::sqrt(global_normsq);
    }

  } else if (style == "base") {
    std::size_t whereid = baseind / localSize();
    if (whereid == myrank) {
      std::size_t lclind = (baseind % localSize());
      state[lclind] = {1.0, 0.0};
    }
// --------------------- added by Gian: beginning
  } else if (style == "++++") {

    state[0] = {1./std::sqrt( globalSize() ),0.};
    for (std::size_t i = 1; i < localSize(); i++) {
      state[i] = state[0];
    }
// --------------------- added by Gian: end
  }

  openqu::mpi::barrier();

#if 0
  double t1 = time_in_seconds();
  if (myrank == 0) {
    printf("[%u] Time to init: %lf\n", myrank, t1 - t0);
  }
#endif
}

// --------------------- FIXME by Gian: different functions depending on the type of Type
template <typename Type>
void QbitRegister<Type>::util_rand_init(std::size_t baseind)
{
  std::cout << " ~~~~~~~~~~~~~~~~ wrong type for state! ~~~~~~~~~~~~~~ \n";
}
//--
template <>
void QbitRegister<ComplexDP>::util_rand_init(std::size_t baseind)
{
    unsigned myrank = openqu::mpi::Environment::rank();
    unsigned nprocs = openqu::mpi::Environment::size();
    // Parallel initialization using parallel MKL RNG
#pragma omp parallel
    {
#ifdef _OPENMP
      std::size_t tid = omp_get_thread_num();
      std::size_t nthreads=omp_get_num_threads();
#else
      std::size_t tid = 0;
      std::size_t nthreads = 1;
#endif
      VSLStreamStatePtr stream;
      std::size_t chunk = localSize() / nthreads;
      std::size_t beg = tid * chunk, end = (tid + 1) * chunk;
      if (tid == nthreads - 1) end = localSize();

      int errcode = vslNewStream(&stream, VSL_BRNG_MCG31, baseind);
      assert(errcode == VSL_STATUS_OK);
      std::size_t nskip = 2UL * (myrank * localSize() + beg);
      vslSkipAheadStream(stream, nskip);
      errcode = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2L * (end - beg),
                        (double *)&state[beg], 0.0, 1.0);
      assert(errcode == VSL_STATUS_OK);
    }
}
//--
template <>
void QbitRegister<ComplexSP>::util_rand_init(std::size_t baseind)
{
    unsigned myrank = openqu::mpi::Environment::rank();
    unsigned nprocs = openqu::mpi::Environment::size();
    // Parallel initialization using parallel MKL RNG
#pragma omp parallel
    {
#ifdef _OPENMP
      std::size_t tid = omp_get_thread_num();
      std::size_t nthreads=omp_get_num_threads();
#else
      std::size_t tid = 0;
      std::size_t nthreads = 1;
#endif
      VSLStreamStatePtr stream;
      std::size_t chunk = localSize() / nthreads;
      std::size_t beg = tid * chunk, end = (tid + 1) * chunk;
      if (tid == nthreads - 1) end = localSize();

      int errcode = vslNewStream(&stream, VSL_BRNG_MCG31, baseind);
      assert(errcode == VSL_STATUS_OK);
      std::size_t nskip = 2L * (myrank * localSize() + beg);
      vslSkipAheadStream(stream, nskip);
      errcode = vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2L * (end - beg),
                        (float *)&state[beg], 0.0, 1.0);
      assert(errcode == VSL_STATUS_OK);
    }
}

template <class Type>
QbitRegister<Type>::QbitRegister(const QbitRegister &in)
{
  Allocate(in.nqbits, in.tmpSize());
  std::size_t lcl = localSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
  TODO(Remember to find 'omp parallel for simd' equivalent for gcc)
#endif

  for (std::size_t i = 0; i < lcl; i++) state[i] = in.state[i];
  *permutation = *(in.permutation);
}

template <class Type>
void QbitRegister<Type>::specializeon()
{
  unsigned myrank = openqu::mpi::Environment::rank();
  if (!myrank) printf("Specialization is on\n");
  specialize = true;
}

template <class Type>
void QbitRegister<Type>::specializeoff()
{
  unsigned myrank = openqu::mpi::Environment::rank();
  if (!myrank) printf("Specialization is off\n");
  specialize = false;
}


template <class Type>
QbitRegister<Type>::~QbitRegister()
{
#if defined(USE_MM_MALLOC)
  _mm_free(state); 
#endif
  if (timer) delete timer;
  if (permutation) delete permutation;
}

template class QbitRegister<ComplexSP>;
template class QbitRegister<ComplexDP>;
