// Copyright (C) 2015 Theoretical Physics, ETH Zurich

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/** @file tinymatrix_test.cpp
 *
 *  Tests for tinymatrix_test.cpp
 */

#include "tinymatrix.hpp"
#include <complex>
#include <iostream>

template <class T, unsigned N, unsigned M>
void testsize()
{
  // create an aligned and an unaligned matrix
  openqu::TinyMatrix<T, N, M> mat;
  openqu::TinyMatrix<T, N, M, 32> mata;

  assert(mat.numRows() == N);
  assert(mat.numCols() == M);
  assert(mat.size() == N * M);

  // test assignment and construction with different aligns

  // test element-wise assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
    for (unsigned j = 0; j < mat.numCols(); ++j)

      // test element-wise assignment
      for (unsigned i = 0; i < mat.numRows(); ++i)
        for (unsigned j = 0; j < mat.numCols(); ++j) {
          assert(mat[i][j] == mat(i, j));
          mat(i, j) = 1. + i + j;
        }

  // test const versrion
  openqu::TinyMatrix<T, N, M> const& matc(mat);

  // test assignment
  for (unsigned i = 0; i < mat.numRows(); ++i)
    for (unsigned j = 0; j < mat.numCols(); ++j) assert(mat(i, j) == 1. + i + j);

  mata = mat;

  // test copy and comparison
  openqu::TinyMatrix<T, N, M> matb = matc;
  openqu::TinyMatrix<T, N, M, 32> matd = matc;

  assert(matb == mat);
  assert(matb == mata);
  assert(matb == matc);
  assert(matc == matb);
  assert(matc == matc);
  assert(matd == mat);
  assert(matd == mata);
  assert(matd == matc);
  assert(matc == matd);

  if (mat.size() != 0) assert(&mat(0, 0) == mat.getPtr());

  // test assignments
  openqu::TinyMatrix<T, N, M> mate;
  openqu::TinyMatrix<T, N, M> matf;
  mate = matc;
  matf = matc;
  assert(mat == mate);
  assert(mat == matf);

  // test not equal
  if (mata.size()) {
    mata(0, 0) = -100.;
    assert(mata != matc);
    assert(matc != mata);
  }
}

template <class T>
void testassign()
{
  double init[2][2] = {{1., 2.}, {3., 4.}};

  openqu::TinyMatrix<T, 2, 2> mat = init;
  for (unsigned i = 0; i < mat.numRows(); ++i)
    for (unsigned j = 0; j < mat.numCols(); ++j) assert(mat(i, j) == 1. + 2. * i + j);

#ifndef __INTEL_COMPILER
  openqu::TinyMatrix<T, 2, 2> mat2 = {{T(1.), T(2.)}, {T(3.), T(4.)}};
  for (unsigned i = 0; i < mat2.numRows(); ++i)
    for (unsigned j = 0; j < mat2.numCols(); ++j) assert(mat2(i, j) == 1. + 2. * i + j);

  // test for equality and inequality

  assert(mat == init);
  assert(mat2 == init);
  assert(mat == mat2);

  assert(!(mat != init));
  assert(!(mat2 != init));
  assert(!(mat != mat2));

  mat = {{T(0.), T(1.)}, {T(1.), T(2.)}};
  for (unsigned i = 0; i < mat.numRows(); ++i)
    for (unsigned j = 0; j < mat.numCols(); ++j) assert(mat(i, j) == 0. + i + j);
#endif
}

template <class T>
void test()
{
  testsize<T, 1, 1>();
  testsize<T, 1, 2>();
  testsize<T, 2, 2>();
  testassign<T>();
}
int main()
{
  test<float>();
  test<float>();
  test<int>();
  test<std::complex<double>>();

  return 0;
}
