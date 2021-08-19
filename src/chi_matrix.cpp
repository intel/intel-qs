/// @file chi_matrix.hpp
/// @brief Declare the @c ChiMatrix class.

#include <cassert>
#include <complex>
#include <initializer_list>
#include <iostream>
#include <vector>

#include <typeinfo>

//#include "Dense"  // Eigen
#include "../extern/eigen/Eigen/Dense"  // Eigen

#include "tinymatrix.hpp"
#include "utils.hpp"
#include "chi_matrix.hpp"

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

  template <class ValueType, unsigned M, unsigned align>
  void ChiMatrix<ValueType, M, align>::EigensystemOfIdealHadamardChannel()
  {
      std::cout << "---- dummy version of EigensystemOfIdealHadamardChannel()\n";//FIXME delete
  }
 
  /// Decompose chi-matrix in its eigensystem (eigenvalues not normalized as probabilities).
  template <class ValueType, unsigned M, unsigned align>
  void ChiMatrix<ValueType, M, align>::SolveEigenSystem()
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
  template <class ValueType, unsigned M, unsigned align>
  void ChiMatrix<ValueType, M, align>::Print(bool with_eigensystem)
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

  /// Normalize 'probabilities of eigenvector'.
  template <class ValueType, unsigned M, unsigned align>
  void ChiMatrix<ValueType, M, align>::NormalizeEigenProbAndRenormalizeEigenVect()
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

/////////////////////////////////////////////////////////////////////////////////////////

////////template <class ValueType, unsigned M, unsigned align = alignof(ValueType)>
 
  // specialization
  template <>
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
  
  // used specializations
  template class ChiMatrix<ComplexDP,  1, 32>;
  template class ChiMatrix<ComplexDP,  2, 32>;
  template class ChiMatrix<ComplexDP,  4, 32>;
  template class ChiMatrix<ComplexDP, 16, 32>;
  template class ChiMatrix<ComplexDP,  2 , 8>;

  template class ChiMatrix<ComplexSP,  4, 32>;
  template class ChiMatrix<ComplexSP, 16, 32>;
}	// end namespace iqs

/////////////////////////////////////////////////////////////////////////////////////////


