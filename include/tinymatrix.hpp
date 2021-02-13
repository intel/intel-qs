// Copyright (C) 2015 Theoretical Physics, ETH Zurich
// 
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

#ifndef IQS_TINYMATRIX_HPP
#define IQS_TINYMATRIX_HPP

#include <cassert>
#include <initializer_list>
#include <iostream>

/// \addtogroup util
/// @{

/// @file tinymatrix.hpp
///
/// This header defines the class template @c TinyMatrix, which stores a matrix
/// whose (small) dimensions are fixed at compile time

namespace iqs {

/// @brief A small matrix with dimensions fixed at compile time
///
/// The matrix is stored intenally as a two-dimensional C array, and thus in row-major
/// ordering.

/////////////////////////////////////////////////////////////////////////////////////////

template <class ValueType, unsigned M, unsigned N = M, unsigned align = alignof(ValueType)>
class TinyMatrix
{
 public:
  /// The type of elements stored in the matrix.
  using value_type = ValueType;
  /// A pointer to elements of the matrix.
  using pointer = ValueType*;
  /// A pointer o elements of a const matrix.
  using const_pointer = ValueType const*;
  /// A reference to elements of the matrix.
  using reference = ValueType&;
  /// An integral type large enought to store the size of the matrix.
  using size_type = unsigned;

  /// The type for a row of the matrix.
  using RowType = ValueType[N];

  /// Default-initialize all matrix elements.
  TinyMatrix() { static_assert(N * M != 0, "A zero-dimensional matrix is not allowed."); }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Initialize from a C-style array of the same dimensions.
  template <class U>
  TinyMatrix(U init[M][N])
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            data_[i][j] = init[i][j];
  }

  /// Initialize from an initializer list, i.e. a compile time given matrix.
  template <class U>
  TinyMatrix(std::initializer_list<std::initializer_list<U>> const& init)
  {
    size_type i = 0;
    for (auto const& line : init)
    {
        size_type j = 0;
        for (auto const& elem : line)
            data_[i][j++] = elem;
        ++i;
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Copy from a matrix with a potentially different type and alignment.
  template <class U, unsigned alignrhs>
  TinyMatrix(TinyMatrix<U, M, N, alignrhs> const& rhs)
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            data_[i][j] = rhs(i, j);
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// The default copy constructor.
  TinyMatrix(TinyMatrix const&) = default;

  /// The default assignment.
  TinyMatrix& operator=(TinyMatrix const&) = default;

/////////////////////////////////////////////////////////////////////////////////////////

  /// Assign from a matrix with a potentially different type and alignment.
  template <class U, unsigned alignrhs>
  TinyMatrix& operator=(TinyMatrix<U, M, N, alignrhs> const& rhs)
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            data_[i][j] = rhs(i, j);
    return *this;
  }

  /// Assign from a C-style array.
  template <class U>
  TinyMatrix& operator=(U const (&rhs)[M][N])
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            data_[i][j] = rhs[i][j];
    return *this;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// The number of matrix rows.
  constexpr size_type numRows() const { return M; }

  /// The number of matrix columns.
  constexpr size_type numCols() const { return N; }

  /// \brief The size of the matrix, i.e. the number of matrix elements.
  /// This is the same as number of rows times  number of columns
  constexpr size_type size() const { return N * M; }

  /// Access a matrix element of a const matrix
  ///   \param i the row index
  ///   \param j the column index
  ///   \pre i<numRows() & j<numCols()
  value_type operator()(size_type i, size_type j) const
  {
    assert(i < this->numRows() && "Row index out of range");
    assert(j < this->numCols() && "Column index out of range");
    return data_[i][j];
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Access a matrix element.
  ///   \param i the row index
  ///   \param j the column index
  ///   \pre i<numRows() & j<numCols()
  reference operator()(size_type i, size_type j)
  {
    assert(i < this->numRows() && "Row index out of range");
    assert(j < this->numCols() && "Column index out of range");
    return data_[i][j];
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Compare two matrices element-wise for equality.
  template <class U, unsigned alignrhs>
  bool operator==(TinyMatrix<U, M, N, alignrhs> const& rhs) const
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            if (data_[i][j] != rhs(i, j))
                return false;
    return true;
  }

  /// Compare two matrices element-wise for inequality.
  template <class U, unsigned alignrhs>
  bool operator!=(TinyMatrix<U, M, N, alignrhs> const& rhs) const
  {
    return !(*this == rhs);
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Compare two matrices element-wise for equality.
  template <class U>
  bool operator==(U const (&rhs)[M][N])
  {
    for (size_type i = 0; i < this->numRows(); ++i)
        for (size_type j = 0; j < this->numCols(); ++j)
            if (data_[i][j] != rhs[i][j]) return false;
    return true;
  }

  /// Compare two matrices element-wise for inequality.
  template <class U>
  bool operator!=(U const (&rhs)[M][N])
  {
    return !(*this == rhs);
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Obtain a pointer to the first element of the matrix
  const_pointer getPtr() const { return &data_[0][0]; }

  /// C-style array subscript
  ///
  /// the TinyMatrix can be indexed both using the mat(i,j) syntax or the
  /// C-style
  /// mat[i][j] syntax

  RowType& operator[](unsigned i)
  {
    assert(i < this->numRows());
    return data_[i];
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// C-style array subscript for a const matrix
  ///
  /// the TinyMatrix can be indexed both using the mat(i,j) syntax or the
  /// C-style
  /// mat[i][j] syntax

  RowType const& operator[](unsigned i) const
  {
    assert(i < this->numRows());
    return data_[i];
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Get submatrices
  ///
  /// Returns the submatrix starting at i_start, j_start of size MSub, NSub
  /// using stride i_stride, j_stride
  /// \param i_start The starting row index
  /// \þaram j_start The starting column index
  /// \param i_stride The row stride to use for accessing elements
  /// \param j_stride The column stride to use for accessing elements
  ///
  /// \tparam MSub The number of rows of the submatrix
  /// \tparam NSub The number of columns of the submatrix
  ///
  /// \pre Strides are strictly positive
  /// \pre Parameters actually represent a submatrix (no index out of bounds)

  template <unsigned MSub, unsigned NSub = MSub>
  TinyMatrix<ValueType, MSub, NSub, align> getSubMatrix(unsigned i_start = 0,
                                                        unsigned j_start = 0,
                                                        unsigned i_stride = 1,
                                                        unsigned j_stride = 1) const
  {
    assert(i_stride > 0 && j_stride > 0);
    assert((MSub - 1) * i_stride + i_start < M && (NSub - 1) * j_stride + j_start < N);

    TinyMatrix<ValueType, MSub, NSub, align> tmp;

    for (unsigned i = i_start, i_sub = 0; i_sub < MSub; i += i_stride, ++i_sub)
      for (unsigned j = j_start, j_sub = 0; j_sub < NSub; j += j_stride, ++j_sub)
        tmp(i_sub, j_sub) = (*this)(i, j);

    return tmp;
  }

/////////////////////////////////////////////////////////////////////////////////////////

 void print(std::string name) 
 {
   printf("name: %s\n", (const char *)name.c_str());
   for (size_type i = 0; i < this->numRows(); ++i) {
     for (size_type j = 0; j < this->numCols(); ++j)
     {
        std::cout << data_[i][j].real() << " + i*" <<  data_[i][j].imag() << " ";
     }
     printf("\n"); 
   }
 }

/////////////////////////////////////////////////////////////////////////////////////////

 std::string tostr() const
 {
   std::string str;
   str = "{";
   for (size_type i = 0; i < this->numRows(); ++i) {
     str += "{";
     for (size_type j = 0; j < this->numCols(); ++j)
     {
       char s[2048];
       sprintf(s, "%.3lf+%.3lf ",  data_[i][j].real(), data_[i][j].imag());
       str += s;
     }
     str += "}";
   }
   str += "{";
   return str;
 }

/////////////////////////////////////////////////////////////////////////////////////////

 std::string name;
 protected:
  alignas(align == 0 ? 8 : align) ValueType data_[M][N];
};

/////////////////////////////////////////////////////////////////////////////////////////

}	// end namespace iqs

/// @}*/

#endif	// header guard IQS_TINYMATRIX_HPP
