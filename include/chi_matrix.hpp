/// @file chi_matrix.hpp
/// @brief Declare the @c ChiMatrix class.

#ifndef IQS_CHI_MATRIX_HPP
#define IQS_CHI_MATRIX_HPP

#include <cassert>
#include <complex>
#include <initializer_list>
#include <iostream>
#include <vector>

#include <typeinfo>

#include <Dense>  // Eigen

#include "tinymatrix.hpp"
#include "utils.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// declaration and definition of class ChiMatrix
/////////////////////////////////////////////////////////////////////////////////////////

/// @class ChiMatrix (described in the context of IQS).
///
/// This header defines the class template @c ChiMatrix used to include noise in the simulation.
///  (1)   rho' = \sum_ij  chi_ij  sigma_i . rho . sigma_j^\dagger
/// The order of Pauli matrix in the basis is: {id, X, Y, Z}.
/// In case of multi-qubit channel, the order of Pauli matrices is:
// TODO: check previous choice for Apply2QubitGate() and verify consistency with ApplyQuantumChannel()
/// {id0.id1, id0.X0, id0.Y1, id0.Z1, X0.id1, ..., Z0.Z1}.
///
/// In addition to the explicit matrix, it includes its eigenvalues and eigenvectors.
/// Using the Dirac notation for linear algebra (and not quantum states):
///  (2)   chi = \sum_k  E_k  |E_k> <E_k|
/// where |E_k> represents the matrix  \sum_i E_k,i . \sigma_i
///
/// NOTE: The eigenvectors |E_k> are defined up to a scalar factor.
///       However eq.(2) can be used to determine the eigenvectors up to a phase.
///       We call this procedure 'standardization'.
///
/// The class also provides functionalities to sample from the eigenvalues according
/// to the probability distribution:
///  (3)   p(k) = |E_k| / sum_k |E_k|
///
/// NOTE: When the probability distribution is defined as in eq.(3), the average
///       procedure must include an extra factor (sum_k |E_k|).
///       Here we include such factor as a renormalization of the eigenvectors.
///       We call this procedure 'renormalization'.

namespace iqs {

/// @brief A small squared matrix with dimenstion fixed at compile time.
/// Stored as a 2d vector in row-major ordering.

/////////////////////////////////////////////////////////////////////////////////////////

template <class ValueType, unsigned M, unsigned align = alignof(ValueType)>
class ChiMatrix : public TinyMatrix<ValueType, M, M, align>
{
 public:
  /// The type of elements stored in the matrix.
  using value_type = ValueType;
  /// The base type of elements stored in the matrix.
  typedef typename extract_value_type<ValueType>::value_type BaseType;
  using base_type = BaseType;
  /// An integral type large enought to store the size of the matrix.
  using size_type = unsigned;

  /// The type for a row of the matrix.
  using RowType = ValueType[M];
 
  /// Default-initialize all matrix elements.
  ChiMatrix()
  {
    static_assert(M != 0, "A zero-dimensional matrix is not allowed.");
    //this->SolveEigenSystem();
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize from a C-style array of the same dimensions.
  template <class U>
  ChiMatrix(U init[M][M])
  : TinyMatrix<ValueType, M, M, align>::TinyMatrix(init) { this->SolveEigenSystem(); }

  /// Initialize from an initializer list, i.e. a compile time given matrix.
  template <class U>
  ChiMatrix(std::initializer_list<std::initializer_list<U>> const& init)
  : TinyMatrix<ValueType, M, M, align>::TinyMatrix(init) { this->SolveEigenSystem(); }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Copy from a matrix with a potentially different type and alignment.
  template <class U, unsigned alignrhs>
  ChiMatrix(ChiMatrix<U, M, alignrhs> const& rhs)
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            this->data_[i][j] = rhs(i, j);
    this->SolveEigenSystem();
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// The default copy constructor.
  ChiMatrix(ChiMatrix const&) = default;

  /// The default assignment.
  ChiMatrix& operator=(ChiMatrix const&) = default;

/////////////////////////////////////////////////////////////////////////////////////////

  /// Return non-constant pointer to data_.
  ValueType * GetPtrToData () {return &(this->data_[0][0]);};

/////////////////////////////////////////////////////////////////////////////////////////

  /// Utility function to run an example of noiseless 1-qubit gate as quantum channel.
  // The Hadamard gate can be written as a linear combination of Pauli matrices:
  //   H = 1/sqrt(2) (X+Z)
  // The ideal Hadamard channel is the unitary transformation:
  //   H[rho] = H.rho.H^dag = H.rho.H
  //          = 1/2 ( X.rho.X + X.rho.Z + Z.rho.X + Z.rho.Z )
  // The chi matrix is:
  //                | . . . . |   <-- I
  //   chi(H) = 1/2 | . 1 . 1 |   <-- X
  //                | . . . . |   <-- Y
  //                | . 1 . 1 |   <-- Z
  // Its eigevalues/eigenstats are:
  // * +1 --> {0,1,0,1}/sqrt(2) 
  // *  0 with degeneracy 3 and any set of representative orthonormal vectors like:
  //    {0,1,0,-1}/sqrt(2), {1,0,0,0}, {0,0,1,0}
  //
  // Here: dummy, generic implementation.
  // Specific implementation below (otherwise problem to initialize value_type as {1,0}.
  void EigensystemOfIdealHadamardChannel()
  {
      std::cout << "---- dummy version of EigensystemOfIdealHadamardChannel()\n";//FIXME delete
  }
  
/////////////////////////////////////////////////////////////////////////////////////////
  
  /// Decompose chi-matrix in its eigensystem (eigenvalues not normalized as probabilities).
  void SolveEigenSystem()
  {
    // Initialize
    evalues_.assign(M, 0);
    evectors_.assign(M, evalues_);
    
    // The matrix is stored in this->data_, which is just a c-style array inside TinyMatrix
    // Use eigen map
    Eigen::Map<Eigen::Matrix< ValueType,M,M,Eigen::RowMajor>> mat( &(this->data_[0][0]) );
    
    // Create and run eigensolver
    Eigen::ComplexEigenSolver< Eigen::Matrix<ValueType,M,M,Eigen::RowMajor> > ces;
    ces.compute( mat );
  
    // Transfer over results
    for(int i=0;i<M;i++)
    {
      evalues_[i] = ces.eigenvalues()[i];
      for(int j=0;j<M;j++)
      {
        evectors_[i][j] = ces.eigenvectors().col(i)[j];
      }
    }

    assert(evalues_.size()==M && "Wrong number of eigenvalues of the chi matrix.");
    assert(evectors_.size()==M && "Wrong number of eigenvector of the chi matrix.");

    // At this point it is not guaranteed that:
    //   chi = sum_k E_k |E_k><E_k|
    // since the eigenvectors |E_k> are defined up to a scalar factor.
    // However this property is required in our implementation of quantum channels.
    // The issue can be solver by computing:
    //   <E_k| chi |E_k> = G_k
    // and rescaling the eigenvectors according to:
    //   |E'_k> = sqrt(E_k/G_k) |E_k>
    // We call eigenvectors |E'_k> standardized.

    for (int k=0; k<M; ++k)
    {
      // Compute rescale factor for k-the eigenvector.
      value_type Gk(0);
      for (int i=0; i<M; ++i)
      {
        value_type Gi(0);
        for (int j=0; j<M; ++j)
        {
          Gi += this->data_[i][j] * evectors_[k][j];
        }
        Gk += std::conj(evectors_[k][i]) * Gi;
      }
      // Rescale k-th eigenvector.
      assert(std::imag(Gk)==0 && "Error: rescale factor is not real.");
      for (int i=0; i<M; ++i)
        if (std::abs(Gk)>0)
          evectors_[k][i] *= std::sqrt(evalues_[k]/Gk);
    }

    // Normalize the probabilities and renormalize the eigenvectors.
    this->NormalizeEigenProbAndRenormalizeEigenVect();
  }
 
/////////////////////////////////////////////////////////////////////////////////////////

  /// Print to screen the chi matrix and, optionally, its eigenvalues/eigenvectors.
  void Print(bool with_eigensystem = true)
  {
    std::cout << "chi_matrix :\n"; 
    for (size_type i = 0; i < this->numRows(); ++i)
    {
        for (size_type j = 0; j < this->numCols(); ++j)
            std::cout << this->data_[i][j] << "\t";
        std::cout << "\n";
    }
    if (with_eigensystem==false)
        return;
    std::cout << "eigenvalues :\n"; 
    for (size_type i = 0; i < evalues_.size(); ++i)
        std::cout << evalues_[i] << "\t";
    std::cout << "\neigenprobs :\n"; 
    for (size_type i = 0; i < eprobs_.size(); ++i)
        std::cout << eprobs_[i] << "\t";
    for (size_type i = 0; i < evectors_.size(); ++i)
    {
        std::cout << "\neigenvector " << i << " :\n"; 
        std::vector<value_type> & evector = evectors_[i];
        for (size_type j=0; j<evector.size(); ++j)
            std::cout << evector[j] << "\t";
    }
    std::cout << "\n";
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Return the k-th eigenvalue.
  value_type GetEigenValue(size_type k) const
  { return evalues_[k]; }

  /// Return the probability of the k-th eigenvalue.
  base_type GetEigenProbability(size_type k) const
  { return eprobs_[k]; }

  /// Return the cumulative probability of the k-th eigenvalue.
  base_type GetEigenCumulativeProbability(size_type k) const
  { return ecumprobs_[k]; }

  /// Return the k-th eigenvector.
  std::vector<value_type> GetEigenVector(size_type k) const
  { return evectors_[k]; }

/////////////////////////////////////////////////////////////////////////////////////////

 protected:

  /// Eigenvalues of the chi matrix (stored as vector).
  // FIXME: the type should be base_type, but it depends on the package for the eigensystem solution
  std::vector<value_type> evalues_;
  /// Eigenvectors of the chi matrix (stored as vector of vectors).
  std::vector<std::vector<value_type>> evectors_;
  /// Probability distribution on the eigenvalues.
  std::vector<base_type> eprobs_;
  /// Cumulative probability distribution on the eigenvalues.
  std::vector<base_type> ecumprobs_;

/////////////////////////////////////////////////////////////////////////////////////////

  /// Normalize 'probabilities of eigenvector'.
  void NormalizeEigenProbAndRenormalizeEigenVect()
  {
    base_type normalization = 0;
    base_type not_normalized_probability;
    eprobs_.clear();
    ecumprobs_.clear();
    for (auto iter=evalues_.begin(); iter!=evalues_.end(); ++iter)
    {
        assert(std::imag(*iter)==0 && "Eigenvalues of chi matrix must be real.");
        // The eigenvalues must be non-negative (FIXME: check this property).
        // Gian: At the moment I believe they can be negative too, look at amplitude-damping channel.
//        assert(std::real(*iter)>=0 && "Eigenvalues of chi matrix must be non-negative.");
        not_normalized_probability = std::abs(std::real(*iter));
        eprobs_.push_back(not_normalized_probability);
        normalization += not_normalized_probability;
        ecumprobs_.push_back(normalization);
    }
    if (normalization != 0 && normalization != 1)
    {
        // Renormalize the probability distribution.
        for (size_type i=0; i<eprobs_.size(); ++i)
        {
            eprobs_[i] /= normalization;
            ecumprobs_[i] /= normalization;
        }
        // Renormalize the eigen-vectors accordingly.
        // See documentation: this is done to eliminate the factor (sum_k |E_k|) coming from
        //                    the definition of the probability distribution.
        for (size_type i=0; i<evectors_.size(); ++i)
        {
            std::vector<value_type> & evector = evectors_[i];
            for (size_type j=0; j<evector.size(); ++j)
                evector[j] *= std::sqrt(normalization);
        }
    }
  }

};

/////////////////////////////////////////////////////////////////////////////////////////

////////template <class ValueType, unsigned M, unsigned align = alignof(ValueType)>
 
  // specialization
  template <> inline
  void ChiMatrix<ComplexDP, 4, 32>::EigensystemOfIdealHadamardChannel ()
  {
    unsigned M=4;
    // Verify the form of the chi matrix.
    for (size_type i = 0; i < M; ++i)
        for (size_type j = 0; j < M; ++j)
        {
            if ( (i==1 || i==3) && (j==1 || j==3) )
                assert(this->data_[i][j]==value_type(0.5,0));
            else
                assert(std::norm(this->data_[i][j])==0);
        }
    // Set eigensystem.
    std::vector<value_type> array_zeros(M, 0);
    evalues_ = array_zeros;
    evalues_[0] = {1, 0};
    value_type one = {1, 0};
    value_type zero = {0,0};
    value_type elem = {1/std::sqrt(2),0};
    evectors_.assign(M, array_zeros);
    evectors_[0] = {zero, elem, zero,  elem};
    evectors_[1] = {zero, elem, zero, -elem};
    evectors_[2] = { one, zero, zero,  zero};
    evectors_[3] = {zero, zero,  one,  zero};
    // Normalize the probabilities and renormalize the eigenvectors.
    this->NormalizeEigenProbAndRenormalizeEigenVect();
  }
  
}	// end namespace iqs

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard IQS_CHI_MATRIX_HPP


