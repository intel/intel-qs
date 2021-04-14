#ifndef	TINYMATRIX_TEST_HPP
#define	TINYMATRIX_TEST_HPP

#include <complex>
#include <iostream>

#include "../../include/tinymatrix.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class: tiny mattrix
//////////////////////////////////////////////////////////////////////////////

class TinyMatrixTest : public ::testing::Test
{
 protected:

  TinyMatrixTest()
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

template <class T, unsigned N, unsigned M>
void tinymatrix_testsize()
{
  // create an aligned and an unaligned matrix
  iqs::TinyMatrix<T, N, M> mat;
  iqs::TinyMatrix<T, N, M, 32> mata;

  ASSERT_EQ(mat.numRows(), N);
  ASSERT_EQ(mat.numCols(), M);
  ASSERT_EQ(mat.size(), N * M);

  // test assignment and construction with different aligns

  // test element-wise assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
      {
          //ASSERT_EQ(mat[i][j], mat(i, j)); //TODO: when T=float initialization is -nan
          mat(i, j) = 1. + i + j;
      }

  // test const versrion
  iqs::TinyMatrix<T, N, M> const& matc(mat);

  // test assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 1. + i + j);

  mata = mat;

  // test copy and comparison
  iqs::TinyMatrix<T, N, M> matb = matc;
  iqs::TinyMatrix<T, N, M, 32> matd = matc;

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
  iqs::TinyMatrix<T, N, M> mate;
  iqs::TinyMatrix<T, N, M> matf;
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
void tinymatrix_testassign()
{
  double init[2][2] = {{1., 2.}, {3., 4.}};

  iqs::TinyMatrix<T, 2, 2> mat = init;
  for (unsigned i = 0; i < mat.numRows(); ++i)
      for (unsigned j = 0; j < mat.numCols(); ++j)
          ASSERT_EQ(mat(i, j), 1. + 2. * i + j);

  iqs::TinyMatrix<T, 2, 2> mat2 = {{T(1.), T(2.)}, {T(3.), T(4.)}};
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
// Test macros:

TEST_F(TinyMatrixTest, Float)
{
  tinymatrix_testsize<float, 1, 1>();
  tinymatrix_testsize<float, 1, 2>();
  tinymatrix_testsize<float, 2, 2>();
  tinymatrix_testassign<float>();
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(TinyMatrixTest, Int)
{
  tinymatrix_testsize<int, 1, 1>();
  tinymatrix_testsize<int, 1, 2>();
  tinymatrix_testsize<int, 2, 2>();
  tinymatrix_testassign<int>();
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(TinyMatrixTest, ComplexDP)
{
  tinymatrix_testsize<std::complex<double>, 1, 1>();
  tinymatrix_testsize<std::complex<double>, 1, 2>();
  tinymatrix_testsize<std::complex<double>, 2, 2>();
  tinymatrix_testassign<std::complex<double>>();
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard TINYMATRIX_TEST_HPP
