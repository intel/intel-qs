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
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();
  }

};

//////////////////////////////////////////////////////////////////////////////
// Utility functions for the tests.

template <class T, unsigned N>
void chimatrix_testsize()
{
  // Create an aligned and an unaligned matrix
  qhipster::ChiMatrix<T, N> mat;
  qhipster::ChiMatrix<T, N, 32> mata;

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
  qhipster::ChiMatrix<T, N> const& matc(mat);

  // test assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 1. + i + j);

  mata = mat;

  // test copy and comparison
  qhipster::ChiMatrix<T, N> matb = matc;
  qhipster::ChiMatrix<T, N, 32> matd = matc;

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
  qhipster::ChiMatrix<T, N> mate;
  qhipster::ChiMatrix<T, N> matf;
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

  qhipster::ChiMatrix<T, 2> mat = init;
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 1. + 2. * i + j);

  qhipster::ChiMatrix<T, 2> mat2 = {{T(1.), T(2.)}, {T(3.), T(4.)}};
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
  qhipster::ChiMatrix<T, 2> chimat2 = {{T(1.), T(2.)}, {T(3.), T(4.)}};
  
  chimat2.SolveEigenSystem();
  
  float precision = 1e-7;
  // Eigenvals
  ASSERT_LT( abs(chimat2.GetEigenValue(0) - -0.37228132) , precision);
  ASSERT_LT( abs(chimat2.GetEigenValue(1) -  5.37228132) , precision);
  // Eigenvec 1
  ASSERT_LT( abs(chimat2.GetEigenVector(0)[0] - -0.82456484) , precision);
  ASSERT_LT( abs(chimat2.GetEigenVector(0)[1] -  0.56576746) , precision);
  // Eigenvec 2
  ASSERT_LT( abs(chimat2.GetEigenVector(1)[0] - -0.41597356) , precision);
  ASSERT_LT( abs(chimat2.GetEigenVector(1)[1] - -0.90937671) , precision);
  
  /*
  Ground truth:
  >>> A = np.array([[1.,2],[3,4]])
  >>> la.eig(A)
  (array([-0.37228132+0.j,  5.37228132+0.j]),
  COLUMN eigenvectors:
  array([[-0.82456484, -0.41597356],
         [ 0.56576746, -0.90937671]]))
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
