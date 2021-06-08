#include <climits>

#include "../include/qureg.hpp"
#include <fstream>
/// \addtogroup qureg
/// @{

/// @file qureg_utils.cpp
/// @brief Define the @c QubitRegister methods used as basic operations.

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Overload operator 'compare'.
///
/// The permutation of both QubitRegister objects must be the same.
template <class Type>
bool QubitRegister<Type>::operator==(const QubitRegister &rhs)
{
  assert(rhs.GlobalSize() == GlobalSize());
  assert(rhs.qubit_permutation->map == qubit_permutation->map);
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
/// @brief Return the L_infinity distance between two states, psi and x.
///
/// Before the comparison, x can be multiplied by a complex factor.
/// The permutation of both QubitRegister objects must be the same.
template <class Type>
typename QubitRegister<Type>::BaseType QubitRegister<Type>::MaxAbsDiff(QubitRegister &x, Type sfactor)
{
  assert(LocalSize() == x.LocalSize());
  assert(x.qubit_permutation->map == qubit_permutation->map);
  BaseType lcl_maxabsdiff = -1.0;

  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
  #pragma omp parallel for simd reduction(max : lcl_maxabsdiff)
#else
  // TODO: Remember to find 'omp parallel for simd reduction' equivalent for gcc.
  #pragma omp parallel for reduction(max : lcl_maxabsdiff)
#endif
  for (std::size_t i = 0; i < lcl; i++)
  {
      lcl_maxabsdiff = std::max(lcl_maxabsdiff, std::abs(state[i] - sfactor*x.state[i]));
  }

  BaseType glb_maxabsdiff ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  iqs::mpi::MPI_Allreduce_x(&lcl_maxabsdiff, &glb_maxabsdiff, 1, MPI_MAX, comm);
#else
  glb_maxabsdiff = lcl_maxabsdiff;
#endif

  return glb_maxabsdiff;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Return the amplitude given its global index.
///
/// The global index is expressed in terms of the program qubits:
///   index = pq0 * 2^0 + pq1 * 2^1 + pq2 * 2^2 + ...
/// where pqk indicates the bit {0,1} corresponding to the computational
/// basis state {|0>, |1>} of program qubit k.

template <class Type>
Type QubitRegister<Type>::GetGlobalAmplitude
  (std::size_t global_index) const
{
  assert(global_index < global_size_);
  // Transform the global_index w.r.t. the data qubit order.
  global_index = qubit_permutation->program2data_(global_index);

  // Determine in what (state) rank is the amplitude and what is its local index.
  Type amplitude;
#ifdef INTELQS_HAS_MPI
  std::size_t local_index, hosting_rank;
  hosting_rank = global_index/local_size_;
  local_index = global_index % local_size_;
  // Broadcast the value to all ranks.
  amplitude = state[local_index];
  iqs::mpi::MPI_Bcast_x(&amplitude, hosting_rank, iqs::mpi::Environment::GetStateComm());
#else
  amplitude = state[global_index];
#endif
  return amplitude;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Return the maximum L2 distance between the local parts of two states.

template <class Type>
void QubitRegister<Type>::SetGlobalAmplitude(std::size_t global_index, Type value)
{
  assert(global_index < global_size_);
  // Transform the global_index w.r.t. the data qubit order.
  global_index = qubit_permutation->program2data_(global_index);

  // Determine in what (state) rank is the amplitude and set its local index.
#ifdef INTELQS_HAS_MPI
  std::size_t local_index, hosting_rank;
  hosting_rank = global_index/local_size_;
  local_index = global_index % local_size_;
  if (hosting_rank == iqs::mpi::Environment::GetStateRank())
      state[local_index] = value;
#else
  state[global_index] = value;
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Return the maximum L2 distance between the local parts of two states.
///
/// The L2 distance between the states would be not the 'max', but the 'sum' of all
/// rank-local contributions.
/// The permutation of both QubitRegister objects must be the same.
template <class Type>
typename QubitRegister<Type>::BaseType QubitRegister<Type>::MaxL2NormDiff(QubitRegister &x)
{
  assert(LocalSize() == x.LocalSize());
  assert(x.qubit_permutation->map == qubit_permutation->map);
  BaseType lcl_diff = 0.0;
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
  #pragma omp parallel for simd reduction(+ : lcl_diff)
#else
  // TODO: Remember to find 'omp parallel for simd reduction' equivalent for gcc.
  #pragma omp parallel for reduction(+ : lcl_diff)
#endif
  for (std::size_t i = 0; i < lcl; i++)
  {
      Type r = state[i] - x.state[i];
      lcl_diff += std::norm(r);
  }

  BaseType glb_diff;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  iqs::mpi::MPI_Allreduce_x(&lcl_diff, &glb_diff, 1, MPI_MAX, comm);
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
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd
#else
#pragma omp parallel for 
#endif
  for(std::size_t i = 0; i < lcl; i++)
  {
     state[i] = state[i] / global_norm;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Compute the norm of the state (L2 norm).
template <class Type>
typename QubitRegister<Type>::BaseType QubitRegister<Type>::ComputeNorm()
{
  BaseType local_normsq = 0;
  std::size_t lcl = LocalSize();
#pragma omp parallel for reduction(+ : local_normsq)
  for(std::size_t i = 0; i < lcl; i++)
  {
     local_normsq += std::norm(state[i]);
  }

  BaseType global_normsq;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  iqs::mpi::MPI_Allreduce_x(&local_normsq, &global_normsq, 1, MPI_SUM, comm);
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
/// The permutation of both QubitRegister objects must be the same.
template <class Type>
Type QubitRegister<Type>::ComputeOverlap( QubitRegister<Type> &psi)
{
  assert(LocalSize() == psi.LocalSize());
  assert(psi.qubit_permutation->map == qubit_permutation->map);
  Type local_over = Type(0.,0.);
  BaseType local_over_re = 0.;
  BaseType local_over_im = 0.;
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
  // TODO: use 'simd' keyword?
  #pragma omp parallel for private(local_over) reduction(+ : local_over_re,local_over_im)
#else
  #pragma omp parallel for private(local_over) reduction(+ : local_over_re,local_over_im)
#endif
  for(std::size_t i = 0; i < lcl; i++)
  {
     local_over = std::conj(psi.state[i]) * state[i] ; 
     local_over_re +=  std::real( local_over );
     local_over_im +=  std::imag( local_over );
  }
  
  BaseType global_over_re(0.) , global_over_im(0.) ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  iqs::mpi::MPI_Allreduce_x(&local_over_re, &global_over_re, 1, MPI_SUM, comm);
  iqs::mpi::MPI_Allreduce_x(&local_over_im, &global_over_im, 1, MPI_SUM, comm);
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
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  iqs::mpi::MPI_Allreduce_x(&local_Hp, &global_Hp, 1, MPI_SUM, comm);
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
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  iqs::mpi::MPI_Allreduce_x(&entropy, &global_entropy, 1, MPI_SUM, comm);
  iqs::mpi::MPI_Allreduce_x(&avgselfinfo, &global_avgselfinfo, 1, MPI_SUM, comm);
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
      iqs::mpi::MPI_Allreduce_x(&(m[i]), &(global_m[i]), 1, MPI_SUM, comm);
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
std::string PrintVector(Type *state, std::size_t size, std::size_t num_elements,
                        BaseType &cumulative_probability,
                        Permutation *permutation,
                        int my_data_rank)
{
  std::string str;
  for (std::size_t i = 0; i < size; i++) {
    // std::string bin = dec2bin(myrank * size + i, num_qubits, false);
    std::string bin = permutation->data2program((std::size_t)my_data_rank * size + i);
    char s[4096];
    sprintf(s, "\t%-13.8lf + i * %-13.8lf   %% |%s> p=%lf\n",
            std::real(state[i]), std::imag(state[i]),
            (const char *)bin.c_str(), std::norm(state[i]) );
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

  int my_rank = iqs::mpi::Environment::GetStateRank();
  int nprocs = iqs::mpi::Environment::GetStateSize();
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
#endif
  iqs::mpi::StateBarrier();

  int tag;
  if (my_rank == 0)
  {
      // print permutation
      assert(qubit_permutation);
      printf("qubit permutation: %s\n", qubit_permutation->GetMapStr().c_str());
      std::string s = PrintVector<Type, BaseType>(state, LocalSize(), num_qubits,
                                                  cumulative_probability, qubit_permutation, my_rank);
      printf("%s=[\n", (const char *)x.c_str());
      printf("%s", (const char *)s.c_str());
#ifdef INTELQS_HAS_MPI
      for (std::size_t i = 1; i < nprocs; i++)
      {
          std::size_t len;
          tag = 1000+i;
#ifdef BIGMPI
          MPIX_Recv_x(&len, 1, MPI_LONG, i, tag, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv(&len, 1, MPI_LONG, i, tag, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          s.resize(len);
          tag = i;
#ifdef BIGMPI
          MPIX_Recv_x((void *)(s.c_str()), len, MPI_CHAR, i, tag, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv((void *)(s.c_str()), len, MPI_CHAR, i, tag, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          printf("%s", (const char *)s.c_str());
      }
#endif
  }
  else
  {
#ifdef INTELQS_HAS_MPI
      std::string s = PrintVector(state, LocalSize(), num_qubits, cumulative_probability, qubit_permutation, my_rank);
      std::size_t len = s.length() + 1;
      tag = 1000 + my_rank;
#ifdef BIGMPI
      MPIX_Send_x(&len, 1, MPI_LONG, 0, tag, comm);
      MPIX_Send_x(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, my_rank, comm);
#else
      MPI_Send(&len, 1, MPI_LONG, 0, tag, comm);
      MPI_Send(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, my_rank, comm);
#endif //BIGMPI
#endif
  }

  BaseType glb_cumulative_probability;
#ifdef INTELQS_HAS_MPI
  MPI_Reduce(&cumulative_probability, &glb_cumulative_probability, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
#else
  glb_cumulative_probability = cumulative_probability;
#endif
  if (my_rank == 0)
  {
      assert(my_rank==0);
      printf("]; %% cumulative probability = %lf\n", (double)glb_cumulative_probability);
  }

  iqs::mpi::StateBarrier();
}


template <class Type>
std::string AmplitudesJSON(Type *state, std::size_t size, std::size_t num_elements,
                        Permutation *permutation,
                        int my_data_rank)
{
  std::string str;
  for (std::size_t i = 0; i < size; i++) {
    // std::string bin = dec2bin(myrank * size + i, num_qubits, false);
    std::string bin = permutation->data2program((std::size_t)my_data_rank * size + i);
    char s[4096];
    sprintf(s, "\t\t\"%s\" : %-13.8lf + i * %-13.8lf,\n",
            (const char *)bin.c_str(), std::real(state[i]), std::imag(state[i]) );
    str = str + s;
  }
  // remove the last trailing comma
  if (my_data_rank == iqs::mpi::Environment::GetStateSize() - 1){
    str = str.substr(0, str.size() - 2) + "\n";
  }
  return std::string(str);
}

template <class Type>
void QubitRegister<Type>::ExportAmplitudes(std::string ofname)
{
  int my_rank = iqs::mpi::Environment::GetStateRank();
  int nprocs = iqs::mpi::Environment::GetStateSize();
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
#endif
  iqs::mpi::StateBarrier();

  int tag;
  std::ofstream of_json_file;

  if (my_rank == 0)
  {
      assert(qubit_permutation);

      of_json_file.open(ofname, std::ofstream::app);
      of_json_file << "\t\"amplitudes\" :\n\t{" << std::endl;
      std::string s = AmplitudesJSON<Type>(state, LocalSize(), num_qubits,
                                                     qubit_permutation, my_rank);
      of_json_file << s.c_str();

#ifdef INTELQS_HAS_MPI
      for (std::size_t i = 1; i < nprocs; i++)
      {
          std::size_t len;
          tag = 1000+i;
#ifdef BIGMPI
          MPIX_Recv_x(&len, 1, MPI_LONG, i, tag, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv(&len, 1, MPI_LONG, i, tag, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          s.resize(len);
          tag = i;
#ifdef BIGMPI
          MPIX_Recv_x((void *)(s.c_str()), len, MPI_CHAR, i, tag, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv((void *)(s.c_str()), len, MPI_CHAR, i, tag, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          of_json_file << s.c_str();
      }
#endif
      of_json_file.close();
  }
  else
  {
#ifdef INTELQS_HAS_MPI
      std::string s = AmplitudesJSON(state, LocalSize(), num_qubits, qubit_permutation, my_rank);
      std::size_t len = s.length() + 1;
      tag = 1000 + my_rank;
#ifdef BIGMPI
      MPIX_Send_x(&len, 1, MPI_LONG, 0, tag, comm);
      MPIX_Send_x(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, my_rank, comm);
#else
      MPI_Send(&len, 1, MPI_LONG, 0, tag, comm);
      MPI_Send(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, my_rank, comm);
#endif //BIGMPI
#endif
  }

  if (my_rank == 0)
  {
      assert(my_rank==0);
      of_json_file.open(ofname, std::ofstream::app );
      of_json_file << "\t}\n";
      of_json_file.close();
  }

  iqs::mpi::StateBarrier();
}
/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Only used for the MPI implementation.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::dumpbin(std::string fn)
{
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  unsigned myrank = iqs::mpi::Environment::GetStateRank();
  MPI_Status status;
  MPI_File fh;
  MPI_File_open(comm, (char *)fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  std::size_t size = LocalSize();
  assert(size < INT_MAX);
  MPI_Offset offset = size * UL(myrank * sizeof(Type));

  double t0 = sec();
  iqs::mpi::StateBarrier();
  MPI_File_write_at(fh, offset, (void *)(&(state[0])), size, MPI_DOUBLE_COMPLEX, &status);
  iqs::mpi::StateBarrier();
  double t1 = sec();
  MPI_File_close(&fh);
  if (myrank == 0)
  {
      double bw = D(UL(sizeof(state[0])) * size) / (t1 - t0) / 1e6;
      printf("Dumping state to %s took %lf sec (%lf MB/s)\n", (const char *)fn.c_str(), (double)(t1 - t0), (double)bw);
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
  int myrank=0, nprocs=1;
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();

  assert(timer == nullptr);
  timer = new Timer(num_qubits, myrank, nprocs);
  assert(gate_counter == nullptr);
  gate_counter = new GateCounter(num_qubits);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Print the statistics (time used in computation and bandwidth) to screen.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::GetStatistics()
{
  assert(timer != nullptr);
  timer->Breakdown();
  assert(gate_counter != nullptr);
  gate_counter->Breakdown();
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Disable the statistics and delete the timer and gate counter.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::DisableStatistics()
{
  assert(timer != nullptr);
  delete timer;
  timer = nullptr;

  assert(gate_counter != nullptr);
  delete gate_counter;
  gate_counter = nullptr;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Reset the statistics to allow a new start.
// TODO: rename it DisableStatistics()?
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::ResetStatistics()
{
  this->DisableStatistics();
  this->EnableStatistics();
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

} // end namespace iqs

/// @}
