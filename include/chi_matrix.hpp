/// @file chi_matrix.hpp
/// @brief Declare the @c ChiMatrix class.

#ifndef IQS_CHI_MATRIX_HPP
#define IQS_CHI_MATRIX_HPP

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
/// NOTE: The eigenvalues and eigenvectors have to be explicitly computed by calling
///       method SolveEigenSystem() after initializing all values of the chi matrix.
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
  : TinyMatrix<ValueType, M, M, align>::TinyMatrix(init) { }

  /// Initialize from an initializer list, i.e. a compile time given matrix.
  template <class U>
  ChiMatrix(std::initializer_list<std::initializer_list<U>> const& init)
  : TinyMatrix<ValueType, M, M, align>::TinyMatrix(init) { }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Copy from a matrix with a potentially different type and alignment.
  template <class U, unsigned alignrhs>
  ChiMatrix(ChiMatrix<U, M, alignrhs> const& rhs)
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            this->data_[i][j] = rhs(i, j);
   // If present, copy also quantities related to the eigensystem. 
   evalues_ = rhs.GetEigenValues();
   eprobs_ = rhs.GetEigenProbabilities();
   ecumprobs_ = rhs.GetEigenCumulativeProbabilities();
   evectors_ = rhs.GetEigenVectors();
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
  void EigensystemOfIdealHadamardChannel();
  
/////////////////////////////////////////////////////////////////////////////////////////
  
  /// Decompose chi-matrix in its eigensystem (eigenvalues not normalized as probabilities).
  void SolveEigenSystem();
 
/////////////////////////////////////////////////////////////////////////////////////////

  /// Print to screen the chi matrix and, optionally, its eigenvalues/eigenvectors.
  void Print(bool with_eigensystem = true);

/////////////////////////////////////////////////////////////////////////////////////////
// NOTE: The following methods may be invoked many times in IQS user's program.
//       To avoid unnecessaary overhead, the code does not impose checks on the index k.
//       This also require that method ChiMatrix::SolveEigenSystem() has been called.

  /// Return the k-th eigenvalue.
  value_type GetEigenValue(size_type k) const
  { return evalues_[k]; }

  /// Return all eigenvalues.
  std::vector<value_type> GetEigenValues() const
  { return evalues_; }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Return the probability of the k-th eigenvalue.
  base_type GetEigenProbability(size_type k) const
  { return eprobs_[k]; }

  /// Return the probability of all eigenvalues.
  std::vector<base_type> GetEigenProbabilities() const
  { return eprobs_; }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Return the cumulative probability of the k-th eigenvalue.
  base_type GetEigenCumulativeProbability(size_type k) const
  { return ecumprobs_[k]; }

  /// Return the cumulative probability all eigenvalue.
  std::vector<base_type> GetEigenCumulativeProbabilities() const
  { return ecumprobs_; }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Return the k-th eigenvector.
  std::vector<value_type> GetEigenVector(size_type k) const
  { return evectors_[k]; }

  /// Return all eigenvector.
  std::vector<std::vector<value_type>> GetEigenVectors() const
  { return evectors_; }

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
  void NormalizeEigenProbAndRenormalizeEigenVect();

};

/////////////////////////////////////////////////////////////////////////////////////////

////////template <class ValueType, unsigned M, unsigned align = alignof(ValueType)>
 
  // specialization
//  template <> inline
//  void ChiMatrix<ComplexDP, 4, 32>::EigensystemOfIdealHadamardChannel ();
  
}	// end namespace iqs

/////////////////////////////////////////////////////////////////////////////////////////

#endif // header guards IQS_CHI_MATRIX_HPP
