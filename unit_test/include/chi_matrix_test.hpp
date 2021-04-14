#ifndef	CHI_MATRIX_TEST_HPP
#define	CHI_MATRIX_TEST_HPP


#include <complex>
#include <iostream>
//#include <math.h>

#include "../../include/chi_matrix.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class: chi matrix
//////////////////////////////////////////////////////////////////////////////

class ChiMatrixTest : public ::testing::Test
{
 protected:

  ChiMatrixTest()
  { }

  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();
  }

};

//////////////////////////////////////////////////////////////////////////////
// Utility functions for the tests.

template <class T, unsigned N>
void chimatrix_testsize()
{
  // Create an aligned and an unaligned matrix
  iqs::ChiMatrix<T, N> mat;
  iqs::ChiMatrix<T, N, 32> mata;

  ASSERT_EQ(mat.numRows(), N);
  ASSERT_EQ(mat.numCols(), N);
  ASSERT_EQ(mat.size(), N * N);

  // test assignment and construction with different aligns

  // test element-wise assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
      {
          ASSERT_EQ(mat[i][j], mat(i, j));
          mat(i, j) = 1. + i + j;
      }

  // test const versrion
  iqs::ChiMatrix<T, N> const& matc(mat);

  // test assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 1. + i + j);

  mata = mat;

  // test copy and comparison
  iqs::ChiMatrix<T, N> matb = matc;
  iqs::ChiMatrix<T, N, 32> matd = matc;

  ASSERT_EQ(matb, mat);
  ASSERT_EQ(matb, mata);
  ASSERT_EQ(matb, matc);
  ASSERT_EQ(matc, matb);
  ASSERT_EQ(matc, matc);
  ASSERT_EQ(matd, mat);
  ASSERT_EQ(matd, mata);
  ASSERT_EQ(matd, matc);
  ASSERT_EQ(matc, matd);

  if (mat.size() != 0)
      ASSERT_EQ(&mat(0, 0), mat.getPtr());

  // test assignments
  iqs::ChiMatrix<T, N> mate;
  iqs::ChiMatrix<T, N> matf;
  mate = matc;
  matf = matc;
  ASSERT_EQ(mat, mate);
  ASSERT_EQ(mat, matf);

  // test not equal
  if (mata.size())
  {
      mata(0, 0) = -100.;
      ASSERT_TRUE(mata != matc);
      ASSERT_TRUE(matc != mata);
  }
}

//////////////////////////////////////////////////////////////////////////////

template <class T>
void chimatrix_testassign()
{
  double init[2][2] = {{1., 2.}, {3., 4.}};

  iqs::ChiMatrix<T, 2> mat = init;
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 1. + 2. * i + j);

  iqs::ChiMatrix<T, 2> mat2 = {{T(1.), T(2.)}, {T(3.), T(4.)}};
  for (unsigned i = 0; i < mat2.numRows(); ++i)
      for (unsigned j = 0; j < mat2.numCols(); ++j)
          ASSERT_EQ(mat2(i, j), 1. + 2. * i + j);

  // test for equality and inequality
  ASSERT_TRUE(mat  == init);
  ASSERT_TRUE(mat  == init);
  ASSERT_TRUE(mat2 == init);
  ASSERT_TRUE(mat  == mat2);

#if 10
  ASSERT_FALSE(mat  != init);
  ASSERT_FALSE(mat2 != init);
  ASSERT_FALSE(mat  != mat2);
#endif

  mat = {{T(0.), T(1.)}, {T(1.), T(2.)}};
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 0. + i + j);
}

//////////////////////////////////////////////////////////////////////////////

template <class T>
void chimatrix_testeigen()
{
  iqs::ChiMatrix<T, 2> chimat2 = {{T(1.), T(3.)}, {T(3.), T(7.)}};
  
  chimat2.SolveEigenSystem();

  std::cout << chimat2.GetEigenVector(0)[0] << "  " << chimat2.GetEigenVector(0)[1] << "\n";
  std::cout << chimat2.GetEigenVector(1)[0] << "  " << chimat2.GetEigenVector(1)[1] << "\n";

  float precision = 1e-7;
  // Eigenvals
  ASSERT_LT( abs(chimat2.GetEigenValue(0) - -0.24264069) , precision);
  ASSERT_LT( abs(chimat2.GetEigenValue(1) -  8.24264069) , precision);
  // Eigenvec 1
  ASSERT_LT( abs(chimat2.GetEigenVector(0)[0] - -2.69121547) , precision);
  ASSERT_LT( abs(chimat2.GetEigenVector(0)[1] -  1.11473795) , precision);
  // Eigenvec 2
  ASSERT_LT( abs(chimat2.GetEigenVector(1)[0] - -1.11473795) , precision);
  ASSERT_LT( abs(chimat2.GetEigenVector(1)[1] - -2.69121547) , precision);
  
  /*
  Ground truth:
  >>> A = np.array([[1.,3],[3,7]])
  >>> la.eig(A)
  
  Eigenvalues: [-0.24264069+0.j  8.24264069+0.j]
  Eigenvectors: [-0.92387953  0.38268343]  and  [-0.38268343 -0.92387953] 

  This eigenvectors are already standardized, but one still needs to normalize them:
  (normalized) Eigenvectors: [-2.69121547  1.11473795]  and  [-1.11473795 -2.69121547] 
  */

}

//////////////////////////////////////////////////////////////////////////////
// Test macros:

/*
TEST_F(ChiMatrixTest, Float)
{
//  chimatrix_testsize<float, 1>();
  chimatrix_testsize<float, 2>();
  chimatrix_testassign<float>();
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ChiMatrixTest, Int)
{
//  chimatrix_testsize<int, 1>();
  chimatrix_testsize<int, 2>();
  chimatrix_testassign<int>();
}
*/

//////////////////////////////////////////////////////////////////////////////

TEST_F(ChiMatrixTest, ComplexDP)
{
  chimatrix_testsize<std::complex<double>, 1>();
  chimatrix_testsize<std::complex<double>, 2>();
  chimatrix_testassign<std::complex<double>>();
  chimatrix_testeigen<std::complex<double>>();
}

//////////////////////////////////////////////////////////////////////////////

#endif
