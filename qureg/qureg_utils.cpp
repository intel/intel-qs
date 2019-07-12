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

/// @file qureg_utils.cpp
///  @brief Define the @c QubitRegister methods used as basic operations.

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
bool QubitRegister<Type>::operator==(const QubitRegister &rhs)
{
  assert(rhs.GlobalSize() == GlobalSize());
  for (std::size_t i = 0; i < rhs.LocalSize(); i++)
  {
      if (state[i] != rhs.state[i])
      {
          printf("(%lf,%lf) != (%lf,%lf)\n", state[i].real(), state[i].imag(),
                                             rhs.state[i].real(), rhs.state[i].imag());
          return false;
      }
  }
  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
QubitRegister<Type>::BaseType QubitRegister<Type>::maxabsdiff(QubitRegister &x, Type sfactor)
{
  assert(LocalSize() == x.LocalSize());
  BaseType lcl_maxabsdiff = -1.0;

  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd reduction(max : lcl_maxabsdiff)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++) {
    lcl_maxabsdiff = std::max(lcl_maxabsdiff, std::abs(state[i] - sfactor*x.state[i]));
  }

  BaseType glb_maxabsdiff ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&lcl_maxabsdiff, &glb_maxabsdiff, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce_x(&lcl_maxabsdiff, &glb_maxabsdiff,  MPI_MAX, comm);
#else
  glb_maxabsdiff = lcl_maxabsdiff;
#endif

  return glb_maxabsdiff;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
QubitRegister<Type>::BaseType QubitRegister<Type>::maxl2normdiff(QubitRegister &x)
{
  assert(LocalSize() == x.LocalSize());
  BaseType lcl_diff = 0.0;
  // #pragma omp parallel for simd reduction(+:lcl_diff)
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : lcl_diff)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++)
  {
      Type r = state[i] - x.state[i];
      lcl_diff += std::norm(r);
  }

  BaseType glb_diff;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&lcl_diff, &glb_diff, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce_x(&lcl_diff, &glb_diff,  MPI_MAX, comm);
#else
   glb_diff = lcl_diff;
#endif

  return glb_diff;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Normalize the quantum state (L2 norm).
template <class Type>
void QubitRegister<Type>::Normalize() 
{
  BaseType global_norm = ComputeNorm();
  std::size_t lcl = LocalSize();
#pragma omp parallel for 
  for(std::size_t i = 0; i < lcl; i++)
  {
     state[i] = state[i] / global_norm;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Compute the norm of the state (L2 norm).
template <class Type>
QubitRegister<Type>::BaseType QubitRegister<Type>::ComputeNorm()
{
  BaseType local_normsq = 0;
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : local_normsq)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for(std::size_t i = 0; i < lcl; i++)
  {
     local_normsq += std::norm(state[i]);
  }

  BaseType global_normsq;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&local_normsq, &global_normsq, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_normsq, &global_normsq,  MPI_SUM, comm);
#else
  global_normsq = local_normsq;
#endif

  return std::sqrt(global_normsq);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Compute the overlap <psi|this state>
/// @param psi is the second state
///
/// The overlap between this state and another state |psi\>
/// is define by:\n
///     \<psi|this state\>
template <class Type>
Type QubitRegister<Type>::ComputeOverlap( QubitRegister<Type> &psi)
{
  Type local_over = Type(0.,0.);
  BaseType local_over_re = 0.;
  BaseType local_over_im = 0.;
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for private(local_over) reduction(+ : local_over_re,local_over_im)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for(std::size_t i = 0; i < lcl; i++)
  {
     local_over = std::conj(psi[i]) * state[i] ; 
     local_over_re +=  std::real( local_over );
     local_over_im +=  std::imag( local_over );
  }
  
  BaseType global_over_re(0.) , global_over_im(0.) ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&local_over_re, &global_over_re, 1, MPI_DOUBLE, MPI_SUM, comm);
  // MPI_Allreduce(&local_over_im, &global_over_im, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_over_re, &global_over_re,  MPI_SUM, comm);
  MPI_Allreduce_x(&local_over_im, &global_over_im,  MPI_SUM, comm);
#else
  global_over_re = local_over_re;
  global_over_im = local_over_im;
#endif

  return Type(global_over_re,global_over_im);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
double QubitRegister<Type>::Entropy()
{

  std::size_t lcl = LocalSize();
  double local_Hp = 0;

  if(timer) timer->Start("ENT", 0);

  double ttot = 0., ttmp1 = sec();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : local_Hp)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++)
  {
      double pj = std::norm(state[i]) ;
      if (pj != double(0.))
      {
          local_Hp -= pj * std::log(pj);
      }
  }

  double global_Hp;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&local_Hp, &global_Hp, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_Hp, &global_Hp,  MPI_SUM, comm);
#else
  global_Hp = local_Hp;
#endif
  
  ttot = sec() - ttmp1;
 
  if (timer)
  {
      double datab = D(sizeof(state[0])) * D(lcl) / ttot;
      timer->record_sn(ttot, datab / ttot);
      timer->Stop();
  }

  return global_Hp / (double)log(double(2.0));
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
std::vector<double> QubitRegister<Type>::GoogleStats()
{
  std::vector <double> stats;

  std::size_t lcl = LocalSize();
  double two2n = D(GlobalSize());
  
  double entropy = 0, avgselfinfo=0,
         m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, 
         m7 = 0, m8 = 0, m9 = 0, m10 = 0; 

  if(timer) timer->Start("ENT", 0);

  double ttot = 0., ttmp1 = sec();

#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : entropy, avgselfinfo, m2, m3, m4, m5, m6, m7, m8, m9, m10)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  #pragma simd reduction(+ : entropy, avgselfinfo, m2, m3, m4, m5, m6, m7, m8, m9, m10)
  // #pragma novector
  for (std::size_t i = 0; i < lcl; i++)
  {
    double pj = std::norm(state[i]) ;
    if (pj != double(0.))
    {
        double nl = log(pj);
        // double nl = pj*pj;
        entropy -= pj * nl;
        avgselfinfo -= nl;
    }
    double pj2  = pj *  pj,
           pj3  = pj2 * pj,
           pj4  = pj2 * pj2,
           pj5  = pj3 * pj2,
           pj6  = pj3 * pj3,
           pj7  = pj4 * pj3,
           pj8  = pj4 * pj4,
           pj9  = pj5 * pj4,
           pj10 = pj5 * pj5;
    m2  += pj2;
    m3  += pj3;
    m4  += pj4;
    m5  += pj5;
    m6  += pj6;
    m7  += pj7;
    m8  += pj8;
    m9  += pj9;
    m10 += pj10;
  }

  double global_entropy;
  double global_avgselfinfo;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  MPI_Allreduce_x(&entropy, &global_entropy,  MPI_SUM, comm);
  MPI_Allreduce_x(&avgselfinfo, &global_avgselfinfo,  MPI_SUM, comm);
#else
  global_entropy = entropy;
  global_avgselfinfo = avgselfinfo;
#endif
  global_entropy /= (double)std::log(double(2.0));
  stats.push_back(global_entropy);
  global_avgselfinfo /= (double)log(double(2.0));
  global_avgselfinfo /= two2n;
  stats.push_back(global_avgselfinfo);

  // compute moments
  std::vector <double> m = {m2, m3, m4, m5, m6, m7, m8, m9, m10},
                       factor(m.size()), 
                       global_m(m.size());
  double factorial = 1.0;
  for(auto i = 0; i < m.size(); i++)
  {
      auto k = i + 2;
      factorial *= D(k);
      factor[i] = pow(two2n, D(k - 1)) / factorial;

      m[i] *= factor[i];
#ifdef INTELQS_HAS_MPI
      MPI_Allreduce_x(&(m[i]), &(global_m[i]),  MPI_SUM, comm);
#else
      global_m[i] = m[i];
#endif
      stats.push_back(global_m[i]);
  }

  ttot = sec() - ttmp1;

  if (timer)
  {
      double datab = D(sizeof(state[0])) * D(lcl) / ttot;
      timer->record_sn(ttot, datab / ttot);
      timer->Stop();
  }

  return stats;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Prepare information about the state amplitudes.
/// @param pcum is the cumulative probability for the (possibly partial) state
///
/// Info provided in the output string:
/// index of the computational state and corresponding amplitude.
/// Partial sum of |amplitude|^2 is computed.
//------------------------------------------------------------------------------
template <class Type, class BaseType>
std::string PrintVector(Type *state, std::size_t size, std::size_t num_qubits,
                        BaseType &cumulative_probability, Permutation *permutation)
{
  std::string str;
  int rank = 0;
  for (std::size_t i = 0; i < size; i++) {
    // std::string bin = dec2bin(rank * size + i, num_qubits, false);
    std::string bin = permutation->lin2perm(rank * size + i);
    char s[4096];
    sprintf(s, "\t%-13.8lf + i * %-13.8lf   %% |%s> p=%lf\n",
            std::real(state[i]), std::imag(state[i]),
            bin.c_str(), std::norm(state[i]) );
    str = str + s;
    cumulative_probability += std::norm(state[i]);
  }
  return std::string(str);
}


/// @brief Print on screen some information about the state.
/// @param x is the message to describe the state to be printed
///
/// The pieces of information that are printed are:
/// - permutation map
/// - all amplitudes of the computational basis
/// - cumulative_probability corresponds to the state norm
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::Print(std::string x, std::vector<std::size_t> qubits)
{
  TODO(Second argument of Print() is not used!)
  BaseType cumulative_probability = 0;

  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();
  openqu::mpi::barrier();
#endif

  if (myrank == 0)
  {
      // print permutation
      assert(permutation);
      printf("permutation: %s\n", permutation->GetMapStr().c_str());
      std::string s = PrintVector<Type, BaseType>(state, LocalSize(), num_qubits,
                                                  cumulative_probability, permutation);
      printf("%s=[\n", x.c_str());
      printf("%s", s.c_str());
#ifdef INTELQS_HAS_MPI
      for (std::size_t i = 1; i < nprocs; i++)
      {
          std::size_t len;
#ifdef BIGMPI
          MPIX_Recv_x(&len, 1, MPI_LONG, i, 1000 + i, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv(&len, 1, MPI_LONG, i, 1000 + i, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          s.resize(len);
#ifdef BIGMPI
          MPIX_Recv_x((void *)(s.c_str()), len, MPI_CHAR, i, i, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv((void *)(s.c_str()), len, MPI_CHAR, i, i, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          printf("%s", s.c_str());
      }
#endif
  }
  else
  {
#ifdef INTELQS_HAS_MPI
      std::string s = PrintVector(state, LocalSize(), num_qubits, cumulative_probability, permutation);
      std::size_t len = s.length() + 1;
#ifdef BIGMPI
      MPIX_Send_x(&len, 1, MPI_LONG, 0, 1000 + myrank, comm);
      MPIX_Send_x(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, myrank, comm);
#else
      MPI_Send(&len, 1, MPI_LONG, 0, 1000 + myrank, comm);
      MPI_Send(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, myrank, comm);
#endif //BIGMPI
#endif
  }

  BaseType glb_cumulative_probability;
#ifdef INTELQS_HAS_MPI
  MPI_Reduce(&cumulative_probability, &glb_cumulative_probability, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
#else
  glb_cumulative_probability = cumulative_probability;
#endif
  if (myrank == 0)
  {
      printf("]; %% cumulative probability = %lf\n", glb_cumulative_probability);
  }

#ifdef INTELQS_HAS_MPI
  openqu::mpi::barrier();
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Only used for the MPI implementation.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::dumpbin(std::string fn)
{
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  unsigned myrank = openqu::mpi::Environment::rank();
  MPI_Status status;
  MPI_File fh;
  MPI_File_open(comm, (char *)fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  std::size_t size = LocalSize();
  assert(size < INT_MAX);
  MPI_Offset offset = size * UL(myrank * sizeof(Type));

  double t0 = sec();
  openqu::mpi::barrier();
  MPI_File_write_at(fh, offset, (void *)(&(state[0])), size, MPI_DOUBLE_COMPLEX, &status);
  openqu::mpi::barrier();
  double t1 = sec();
  MPI_File_close(&fh);
  if (myrank == 0)
  {
      double bw = D(UL(sizeof(state[0])) * size) / (t1 - t0) / 1e6;
      printf("Dumping state to %s took %lf sec (%lf MB/s)\n", fn.c_str(), t1 - t0, bw);
  }
#else
  assert(0);
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Enable the collection of statistics (time used in computation and bandwidth).
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::EnableStatistics()
{
  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
#endif

  assert(timer == NULL);
  timer = new Timer(num_qubits, myrank, nprocs);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Print the statistics (time used in computation and bandwidth) to screen.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::GetStatistics()
{
  assert(timer);
  timer->Breakdown();
  // delete timer;
  // timer = NULL;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Reset the statistics to allow a new start.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::ResetStatistics()
{
// FIXME it does not delete the timer object!
  assert(timer);
  delete timer;
  timer = NULL;
}


template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}
