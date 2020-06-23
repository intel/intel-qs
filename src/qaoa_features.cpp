#include "../include/qaoa_features.hpp"

#include <cassert>		// pow

#ifdef INTELQS_HAS_MPI
#include <mpi.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

namespace qaoa
{

/////////////////////////////////////////////////////////////////////////////////////////

namespace utility
{

/////////////////////////////////////////////////////////////////////////////////////////

/// Function to convert a decimal number into a binary number (expressed as vector).
/// The 0-component of the vector represents the least significant bit.
/// This convention has been used in qHiPSTER but may be unusual in quantum computation
/// where qubit zero-th is usually associated to the most significant bit.
template<typename T_decimal, typename T_bit>
void ConvertToBinary( T_decimal k, std::vector<T_bit> &z )
{
  // If decimal number is too large to be converted into a binary.
  if ( k >= (T_decimal)((T_decimal)1<<(T_decimal)(z.size()) ) )
  {
      std::cout << "Too large decimal number:\n"
                << "decimal = " << k << " , but binary has " << z.size() << "bits.\n";
      assert( 0 );
  }
  for (unsigned pos=0; pos<z.size(); ++pos)
  {
      z[pos] = k%(T_decimal)(2);
      k=k/(T_decimal)(2);
  }
}

template void ConvertToBinary<std::size_t,int> (std::size_t, std::vector<int> &);
template void ConvertToBinary<std::size_t,unsigned> (std::size_t, std::vector<unsigned> &);

////////////////////////////////////////////////////////////////////////////////

/// Function to convert a binary number (expressed as vector) into a decimal number.
/// The 0-componenet of the vector represents the least significant bit (i.e. associated
/// to the factor 2^0 in power expension).
template<typename T_bit, typename T_decimal>
void ConvertToDecimal( std::vector<T_bit> &z , T_decimal &k )
{
  k=0;
  for (unsigned pos=z.size()-1; pos>0; --pos)
  {
      k += (T_decimal)(z[pos]);
      k=k*(T_decimal)(2);
  }
  k += (T_decimal)(z[0]);
}

template void ConvertToDecimal<int,std::size_t> (std::vector<int> &, std::size_t &);
template void ConvertToDecimal<unsigned,std::size_t> (std::vector<unsigned> &, std::size_t &);

////////////////////////////////////////////////////////////////////////////////

}	// end namespace 'utility'

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<typename Type>
int InitializeVectorAsMaxCutCostFunction(QubitRegister<Type> & diag,
                                         std::vector<int> & adjacency)
{
  int num_vertices = diag.NumQubits();
  int num_edges = 0;

  // A few preliminary checks on the adjacency matrix.
  // - it should have the right size:
  assert(adjacency.size() == num_vertices*num_vertices);
  // - it should have a null diagonal:
  for (int v=0; v<num_vertices; ++v)
      assert(adjacency[v*num_vertices+v]==0);
  // - each edge is counted twice
  for (int e=0; e<adjacency.size(); ++e)
      num_edges += adjacency[e];
  assert(num_edges%2==0);
  num_edges /= 2;

  // Denote (x)^T the row vector: {-1,1,1,-1,1,...,1}
  // It indicates how the vertices of the graph are colored (either +1 or -1).
  // In this case,
  //   (x)^T.ADJ.(x) = 2*(num_uncut_edges - num_cut_edges)
  // Therefore:
  //   num_cut_edges = ( num_edges - x^T.ADJ.x /2 ) /2 
 
  std::size_t myrank = qhipster::mpi::Environment::GetStateRank();
  std::size_t glb_start = UL(myrank) * diag.LocalSize();
  int max_cut = 0;

#pragma omp parallel
  {
      std::size_t x;
      std::vector<int> xbin(num_vertices);
      int cut;
      #pragma omp for reduction(max: max_cut)
      for(std::size_t i = 0; i < diag.LocalSize(); i++)
      {
         x = glb_start + i;
         // x is the global index according to the data qubits
         // Transform the global_index w.r.t. the program qubit order.
         x = diag.qubit_permutation->data2program_(x);

         // From decimal to binary vector of {0,1}.
         utility::ConvertToBinary(x,xbin);
         // From binary vector of {0,1} to binary vector of {-1,1}.
         for (int v=0; v<num_vertices; ++v)
             if (xbin[v]==0)
                 xbin[v]=-1;
         // Compute x^T.ADJ.x
         cut = 0;
         for (int v=0; v<num_vertices; ++v)
             for (int u=0; u<num_vertices; ++u)
                 cut += adjacency[v*num_vertices + u] * xbin[v] * xbin[u];
         assert(cut%2==0);
         cut = ( num_edges - cut/2);
         assert(cut%2==0);
         cut /= 2;
         diag[i] = Type(cut,0);
         if (cut>max_cut)
             max_cut = cut;
      }
  }

#ifdef INTELQS_HAS_MPI
  int lcl_max_cut = max_cut;
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  MPI_Allreduce(&lcl_max_cut, &max_cut, 1, MPI_INT, MPI_MAX, comm);
#endif

  return max_cut;
}

template int InitializeVectorAsMaxCutCostFunction<ComplexDP>
    (QubitRegister<ComplexDP> &, std::vector<int> & );
template int InitializeVectorAsMaxCutCostFunction<ComplexSP>
    (QubitRegister<ComplexSP> &, std::vector<int> & );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
typename QubitRegister<Type>::BaseType
InitializeVectorAsWeightedMaxCutCostFunction(QubitRegister<Type> & diag,
    std::vector<typename QubitRegister<Type>::BaseType> & adjacency)
{
  typedef typename QubitRegister<Type>::BaseType BaseType;

  int num_vertices = diag.NumQubits();
  int num_edges = 0;
  BaseType total_weight = 0;

  // A few preliminary checks on the adjacency matrix.
  // - it should have the right size:
  assert(adjacency.size() == num_vertices*num_vertices);
  // - it should have a null diagonal:
  for (int v=0; v<num_vertices; ++v)
      assert(adjacency[v*num_vertices+v]==0);
  // - it should be symmetric
  //   also compute the number of edges
  for (int v1=0; v1<num_vertices; ++v1)
  for (int v2=v1+1; v2<num_vertices; ++v2)
  {
      int index = v1*num_vertices + v2;
      assert(adjacency[index] == adjacency[v2*num_vertices+v1]);
      if (adjacency[index]!=0)
      {
          num_edges += 1;
          total_weight += adjacency[index];
      }
  }

  // Denote (x)^T the row vector: {-1,1,1,-1,1,...,1}
  // It indicates how the vertices of the graph are colored (either +1 or -1).
  // In this case,
  //   (x)^T.ADJ.(x) = 2*(weight_uncut_edges - weight_cut_edges)
  // Therefore:
  //   weight_cut_edges = ( total_weight - x^T.ADJ.x /2 ) /2 
 
  std::size_t myrank = qhipster::mpi::Environment::GetStateRank();
  std::size_t glb_start = UL(myrank) * diag.LocalSize();
  BaseType max_cut = 0;

#pragma omp parallel
  {
      std::size_t x;
      std::vector<int> xbin(num_vertices);
      std::vector<BaseType> xbin_basetype(num_vertices);
      BaseType cut;
      #pragma omp for reduction(max: max_cut)
      for(std::size_t i = 0; i < diag.LocalSize(); i++)
      {
         x = glb_start + i;
         // x is the global index according to the data qubits
         // Transform the global_index w.r.t. the program qubit order.
         x = diag.qubit_permutation->data2program_(x);

         // From decimal to binary vector of {0,1}.
         utility::ConvertToBinary(x,xbin);
         // From binary vector of {0,1} to binary vector of {-1,1}.
         for (int v=0; v<num_vertices; ++v)
         {
             if (xbin[v]==0)
                 xbin_basetype[v]=-1;
             else
                 xbin_basetype[v]= 1;
         }
         // Compute x^T.ADJ.x
         cut = 0;
         for (int v=0; v<num_vertices; ++v)
             for (int u=0; u<num_vertices; ++u)
                 cut += adjacency[v*num_vertices + u] * xbin_basetype[v] * xbin_basetype[u];
         cut = ( total_weight - cut/2);
         cut /= 2;
         diag[i] = Type(cut,0);
         if (cut>max_cut)
             max_cut = cut;
      }
  }

#ifdef INTELQS_HAS_MPI
  BaseType lcl_max_cut = max_cut;
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  qhipster::mpi::MPI_Allreduce_x(&lcl_max_cut, &max_cut, 1, MPI_MAX, comm);
#endif

  return max_cut;
}

template double InitializeVectorAsWeightedMaxCutCostFunction<ComplexDP>
    (QubitRegister<ComplexDP> &, std::vector<double> & );
template float InitializeVectorAsWeightedMaxCutCostFunction<ComplexSP>
    (QubitRegister<ComplexSP> &, std::vector<float> & );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
void ImplementQaoaLayerBasedOnCostFunction(QubitRegister<Type> & psi,
                                           QubitRegister<Type> & diag,
                                           typename QubitRegister<Type>::BaseType gamma)
{
  assert( psi.LocalSize( ) == diag.LocalSize( ) );
  assert( psi.GlobalSize() == diag.GlobalSize() );
  assert( psi.qubit_permutation->map == diag.qubit_permutation->map);

  // NOTE: cosine and sine for all values of the cost function could be computed once
  //       and stored in a vector.

  // exp(-i gamma H_problem)
  #pragma omp parallel for
  for (std::size_t i=0; i < psi.LocalSize(); ++i)
       psi[i] *= Type( std::cos(gamma* diag[i].real()) , -std::sin(gamma* diag[i].real()) );
}

template void ImplementQaoaLayerBasedOnCostFunction<ComplexDP>
    (QubitRegister<ComplexDP> &, QubitRegister<ComplexDP> &, double );
template void ImplementQaoaLayerBasedOnCostFunction<ComplexSP>
    (QubitRegister<ComplexSP> &, QubitRegister<ComplexSP> &, float );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
typename QubitRegister<Type>::BaseType
GetExpectationValueFromCostFunction(const QubitRegister<Type> & psi,
                                    const QubitRegister<Type> & diag)
{
  assert( psi.qubit_permutation->map == diag.qubit_permutation->map);
  // Extract basic type from IQS objects.
  typename QubitRegister<Type>::BaseType global_expectation, local_expectation = 0.;

  #pragma omp parallel for reduction(+: local_expectation)
  for ( size_t i=0 ; i < psi.LocalSize(); ++i)
  {
      local_expectation += diag[i].real() * norm(psi[i]) ;
  }

#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  qhipster::mpi::MPI_Allreduce_x(&local_expectation, &global_expectation, 1, MPI_SUM, comm);
#else
  global_expectation = local_expectation;
#endif
  return global_expectation;
}

template double GetExpectationValueFromCostFunction<ComplexDP>
    (const QubitRegister<ComplexDP> &, const QubitRegister<ComplexDP> & );
template float GetExpectationValueFromCostFunction<ComplexSP>
    (const QubitRegister<ComplexSP> &, const QubitRegister<ComplexSP> & );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
typename QubitRegister<Type>::BaseType
GetExpectationValueSquaredFromCostFunction(const QubitRegister<Type> & psi,
                                           const QubitRegister<Type> & diag)
{
  assert( psi.qubit_permutation->map == diag.qubit_permutation->map);
  // Extract basic type from IQS objects.
  typename QubitRegister<Type>::BaseType global_expectation, local_expectation = 0.;

  #pragma omp parallel for reduction(+: local_expectation)
  for ( size_t i=0 ; i < psi.LocalSize(); ++i)
  {
      local_expectation += diag[i].real() * diag[i].real()  * norm(psi[i]) ;
  }

#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  qhipster::mpi::MPI_Allreduce_x(&local_expectation, &global_expectation, 1, MPI_SUM, comm);
#else
  global_expectation = local_expectation;
#endif
  return global_expectation;
}

template double GetExpectationValueSquaredFromCostFunction<ComplexDP>
    (const QubitRegister<ComplexDP> &, const QubitRegister<ComplexDP> & );
template float  GetExpectationValueSquaredFromCostFunction<ComplexSP>
    (const QubitRegister<ComplexSP> &, const QubitRegister<ComplexSP> & );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
std::vector<typename QubitRegister<Type>::BaseType>
GetHistogramFromCostFunction( const QubitRegister<Type> & psi,
                              const QubitRegister<Type> & diag, int max_value)
{
  // Extract basic type from IQS objects.
  typedef typename QubitRegister<Type>::BaseType Basetype;

  // A few preliminary checks:
  assert( psi.LocalSize( ) == diag.LocalSize( ) );	// Vectors with equal local size.
  assert( psi.GlobalSize() == diag.GlobalSize() );	// Vectors with equal global size.
  assert( psi.qubit_permutation->map == diag.qubit_permutation->map);
  assert(max_value>0);					// The max_value must be positive.

  int my_rank = qhipster::mpi::Environment::GetStateRank();

  // Histogram of the specific MPI (state) rank.
  Basetype local_hist[max_value+1] ;	// Initialize all elements to 0 (only with C++)
  for (int n=0; n<=max_value; ++n)
      local_hist[n]=0;

  #pragma omp parallel
  {
      int cut;
      int index_bin;
      // Histogram of the specific thread.
      Basetype private_hist[max_value+1] ;
      for (int n=0; n<=max_value; ++n)
          private_hist[n]=0;

      #pragma omp for
      for ( size_t i=0 ; i < psi.LocalSize(); ++i)
      {
          cut = (int)(diag[i].real());
          assert( cut>=0 && cut <= max_value );
          // The next line can be changed to have wider bins. In general not required.
          index_bin = cut;
          private_hist[index_bin] += norm(psi[i]) ;
      }
      #pragma omp critical
      {
          for (int n=0; n<=max_value; ++n)
              local_hist[n] += private_hist[n];
      }
  }

  // Global histogram.
  std::vector<Basetype> global_hist(max_value+1,0);
#ifdef INTELQS_HAS_MPI
  // Sum local histograms into (state) global histogram.
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  qhipster::mpi::MPI_Allreduce_x(local_hist, global_hist.data(), max_value+1, MPI_SUM, comm);
#else
  for ( int n=0; n<=max_value; ++n)
      global_hist[n]=local_hist[n];
#endif

  return global_hist;
}

template std::vector<double> GetHistogramFromCostFunction<ComplexDP>
    (const QubitRegister<ComplexDP> &, const QubitRegister<ComplexDP> &, int );
template std::vector<float>  GetHistogramFromCostFunction<ComplexSP>
    (const QubitRegister<ComplexSP> &, const QubitRegister<ComplexSP> &, int );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
std::vector<typename QubitRegister<Type>::BaseType>
GetHistogramFromCostFunctionWithWeightsRounded( const QubitRegister<Type> & psi,
                                                const QubitRegister<Type> & diag,
                                                double max_value)
{
  // Extract basic type from IQS objects. 
  typedef typename QubitRegister<Type>::BaseType Basetype;

  // A few preliminary checks:
  assert( psi.LocalSize( ) == diag.LocalSize( ) );	// Vectors with equal local size.
  assert( psi.GlobalSize() == diag.GlobalSize() );	// Vectors with equal global size.
  assert( psi.qubit_permutation->map == diag.qubit_permutation->map);
  assert(max_value>0);					// The max_value must be positive.

  int my_rank = qhipster::mpi::Environment::GetStateRank();

  // Histogram of the specific MPI (state) rank.
  int num_bins = (int)(floor(max_value))+1;
  Basetype local_hist[num_bins];
  for (int n=0; n<num_bins; ++n)
      local_hist[n]=0;
  #pragma omp parallel
  {
      double cut;
      int index_bin;
      // Histogram of the specific thread.
      Basetype private_hist[num_bins] ;
      for (int n=0; n<num_bins; ++n)
          private_hist[n]=0;

      #pragma omp for
      for ( size_t i=0 ; i < psi.LocalSize(); ++i)
      {
          cut = diag[i].real();
          // single precision: machine epsilon ~ 1e-7
          // double precision: machine epsilon ~ 1e-16
          assert( cut>=-1e-7 && cut <= (max_value + 1e-7) );
          index_bin = (int)(floor(cut+1e-7));
          private_hist[index_bin] += norm(psi[i]) ;
      }
      #pragma omp critical
      {
          for (int n=0; n<num_bins; ++n)
              local_hist[n] += private_hist[n];
      }
  }

  // Global histogram.
  std::vector<Basetype> global_hist(num_bins, 0);
#ifdef INTELQS_HAS_MPI
  // Sum local histograms into (state) global histogram.
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  qhipster::mpi::MPI_Allreduce_x(local_hist, global_hist.data(), num_bins, MPI_SUM, comm);
#else
  for ( int n=0; n<num_bins; ++n)
      global_hist[n]=local_hist[n];
#endif
    
  return global_hist;
}

template std::vector<double> GetHistogramFromCostFunctionWithWeightsRounded<ComplexDP>
    (const QubitRegister<ComplexDP> &, const QubitRegister<ComplexDP> &, double );
template std::vector<float>  GetHistogramFromCostFunctionWithWeightsRounded<ComplexSP>
    (const QubitRegister<ComplexSP> &, const QubitRegister<ComplexSP> &, double );

/////////////////////////////////////////////////////////////////////////////////////////

template<typename Type>
std::vector<typename QubitRegister<Type>::BaseType>
GetHistogramFromCostFunctionWithWeightsBinned( const QubitRegister<Type> & psi,
                                               const QubitRegister<Type> & diag,
                                               double max_value, double bin_width)
{
  // Extract basic type from IQS objects. 
  typedef typename QubitRegister<Type>::BaseType Basetype;

  // A few preliminary checks:
  assert( psi.LocalSize( ) == diag.LocalSize( ) );	// Vectors with equal local size.
  assert( psi.GlobalSize() == diag.GlobalSize() );	// Vectors with equal global size.
  assert( psi.qubit_permutation->map == diag.qubit_permutation->map);
  assert(max_value>0);					// The max_value must be positive.

  int my_rank = qhipster::mpi::Environment::GetStateRank();

  // Histogram of the specific MPI (state) rank.
  int num_bins = (int)(ceil(max_value / bin_width)) + 1;
  Basetype local_hist[num_bins] ;	// Initialize all elements to 0 (only with C++)
  for (int n=0; n<num_bins; ++n)
      local_hist[n]=0;
  #pragma omp parallel
  {
      double cut;
      int index_bin;
      // Histogram of the specific thread.
      Basetype private_hist[num_bins] ;
      for (int n=0; n<num_bins; ++n)
          private_hist[n]=0;
      #pragma omp for
      for ( size_t i=0 ; i < psi.LocalSize(); ++i)
      {
          cut = diag[i].real();
          // single precision: machine epsilon ~ 1e-7
          // double precision: machine epsilon ~ 1e-16
          assert( cut>=-1e-7 && cut <= max_value + 1e-7 );
          index_bin = (int)(floor(cut / bin_width + 1e-7));
          private_hist[index_bin] += norm(psi[i]) ;
      }
      #pragma omp critical
      {
          for (int n=0; n<num_bins; ++n)
              local_hist[n] += private_hist[n];
      }
  }

  // Global histogram.
  std::vector<Basetype> global_hist(num_bins, 0);
#ifdef INTELQS_HAS_MPI
  // Sum local histograms into (state) global histogram.
  MPI_Comm comm = qhipster::mpi::Environment::GetStateComm();
  qhipster::mpi::MPI_Allreduce_x(local_hist, global_hist.data(), num_bins, MPI_SUM, comm);
#else
  for ( int n=0; n<num_bins; ++n)
      global_hist[n]=local_hist[n];
#endif
  return global_hist;
}

template std::vector<double> GetHistogramFromCostFunctionWithWeightsBinned<ComplexDP>
    (const QubitRegister<ComplexDP> &, const QubitRegister<ComplexDP> &, double, double);
template std::vector<float>  GetHistogramFromCostFunctionWithWeightsBinned<ComplexSP>
    (const QubitRegister<ComplexSP> &, const QubitRegister<ComplexSP> &, double, double);

/////////////////////////////////////////////////////////////////////////////////////////

}	// close namespace qaoa

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
