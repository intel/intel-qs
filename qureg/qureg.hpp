//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------
#pragma once

#include "permute.hpp"
#include "util/utils.hpp"
#include "util/timer.hpp"

#include "util/alignedallocator.hpp"
#include "util/mpi.hpp"
#include "util/bitops.hpp"
#include "util/conversion.hpp"
#include "util/tinymatrix.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <fstream>




#if defined(__ICC) || defined(__INTEL_COMPILER)
#include <mkl.h>
#endif
#include <algorithm>  // for std::swap
#include <memory>  // for std::swap
#include <numeric>
#include <functional>
#include <tuple>

template<typename T>
struct extract_value_type //lets call it extract_value_type
{
    typedef T value_type;
};

template<template<typename> class X, typename T>
struct extract_value_type<X<T>>   //specialization
{
    typedef T value_type;
};

template<class Type>
using TM2x2 = openqu::TinyMatrix<Type, 2, 2, 32>;

template<class Type>
using TM4x4 = openqu::TinyMatrix<Type, 4, 4, 32>;


template <class Type = ComplexDP>
class QbitRegister
{
 public:
  using value_type = Type;
  typedef typename extract_value_type<Type>::value_type BaseType;

  // constructors / destructors
  QbitRegister();
  QbitRegister(std::size_t nqbits, std::string style = "", 
               std::size_t baseind = 0, std::size_t tmpspacesize_ = 0);
  QbitRegister(const QbitRegister &in);
  QbitRegister(std::size_t nqbits, Type *state, std::size_t tmpspacesize_ = 0);
  ~QbitRegister();

  // allocation & initialization
  void allocateQubit();
  void Allocate(std::size_t nqbits, std::size_t tmpspacesize_);
  void Init(std::size_t nqbits, std::size_t tmpspacesize_);

  inline Type& operator[] (std::size_t ind) { return state[ind]; }
  inline Type& operator[] (std::size_t ind) const { return state[ind]; }

  std::size_t localSize() const { return localsize_; }
  std::size_t globalSize() const { return globalsize_; }
  void resize(std::size_t n);
  std::size_t size() const { return globalsize_; }
  unsigned numQubits() const { return nqbits; }
  Type *tmpspace() const { return state + localSize(); }
  size_t tmpSize() const {return tmpspacesize_;}

  inline bool check_bit(std::size_t var, std::size_t pos) const
  {
    std::size_t one = (std::size_t)1, posl = UL(pos);
    return var & (one << posl);
  }
  inline std::size_t set_bit(std::size_t var, std::size_t pos) const
  {
    std::size_t one = (std::size_t)1, posl = UL(pos);
    return var | (one << posl);
  }
  inline std::size_t clear_bit(std::size_t var, std::size_t pos) const
  {
     return (var & ~(UL(1) << UL(pos)));
  }

  void EnbStat();
  void GetStat();

  void Permute(std::vector<std::size_t> permutation_new_vec);


  // Generic gates
  // single qubit gates
  bool apply1QubitGate_helper(unsigned qubit,  TM2x2<Type> const&m,
                              std::size_t sstate_ind, std::size_t estate_ind);
  void apply1QubitGate(unsigned qubit, TM2x2<Type> const&m);
  // constrolled gates
  bool applyControlled1QubitGate_helper(unsigned control, unsigned qubit, TM2x2<Type> const&m,
                                      std::size_t sind, std::size_t eind);
  void applyControlled1QubitGate(unsigned control, unsigned qubit, TM2x2<Type> const&m);
  // swap gates
  bool applySwap_helper(unsigned qubit1, unsigned qubit2, TM2x2<Type> const&m);
  void applySwap(unsigned qubit1, unsigned qubit2);
  void applyISwap(unsigned qubit1, unsigned qubit2);
  void apply4thRootISwap(unsigned qubit1, unsigned qubit2);
  void applySqrtISwap(unsigned qubit1, unsigned qubit2);
  void applyISwapRotation(unsigned qubit1, unsigned qubit2, TM2x2<Type> const&m);
  void swap(unsigned b1, unsigned b2);
  // diagonal gates
  void applyDiagSimp(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  void applyDiag(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  void applyDiagControl(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  void applyDiagGeneral(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  // two-qubit gates
  void apply2QubitGate(unsigned const qubit_high, unsigned const qubit_low, TM4x4<Type> const&m);
  // specialized gates
  void applyRotationX(unsigned const qubit, BaseType theta);
  void applyRotationY(unsigned const qubit, BaseType theta);
  void applyRotationZ(unsigned const qubit, BaseType theta);
  void applyPauliX(unsigned const qubit);
  void applyPauliY(unsigned const qubit);
  void applyPauliZ(unsigned const qubit);
  void applyPauliSqrtX(unsigned const qubit);
  void applyPauliSqrtY(unsigned const qubit);
  void applyPauliSqrtZ(unsigned const qubit);
  void applyT(unsigned const qubit);
  void applyToffoli(unsigned const qubit1, unsigned const qubit2, unsigned const qubit3);
  void applyHadamard(unsigned const qubit);

  void applyCRotationX(unsigned const qubit, unsigned const qubit2, BaseType theta);
  void applyCRotationY(unsigned const qubit, unsigned const qubit2, BaseType theta);
  void applyCRotationZ(unsigned const qubit, unsigned const qubit2, BaseType theta);
  void applyCPauliX(unsigned const qubit1, unsigned const qubit2);
  void applyCPauliY(unsigned const qubit1, unsigned const qubit2);
  void applyCPauliZ(unsigned const qubit1, unsigned const qubit2);
  void applyCPauliSqrtZ(unsigned const qubit1, unsigned const qubit2);
  void applyCHadamard(unsigned const qubit1, unsigned const qubit2);

  // fusion  
  void fusionon(unsigned log2llc = 20);
  void fusionoff();
  bool is_fusion_enabled();
  void applyFusedGates();

  // gate specialization (experimental)
  void specializeon();
  void specializeoff();

  // measurement
  bool getClassicalValue(unsigned qubit, BaseType tolerance = 1.e-13) const;
  bool isClassicalBit(unsigned qubit, BaseType tolerance = 1.e-13) const;
  void collapseQubit(unsigned qubit, bool value);
  BaseType getProbability(unsigned qubit);

  // expectation values without state update
  void expectationValueX(unsigned const qubit, BaseType &sum, BaseType coeff=1.);
  void expectationValueY(unsigned const qubit, BaseType &sum, BaseType coeff=1.);
  void expectationValueZ(unsigned const qubit, BaseType &sum, BaseType coeff=1.);
  void expectationValueXX(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueXY(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueXZ(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueYX(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueYY(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueYZ(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueZX(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueZY(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValueZZ(unsigned const qubit, unsigned const qubit2, BaseType &sum, BaseType coeff=1.);
  void expectationValue(std::vector<unsigned> &qubits, std::vector<unsigned> &observables, BaseType &sum, BaseType coeff=1.);

  // utilities
  bool operator==(const QbitRegister &rhs);
  void Init(std::string style, std::size_t baseind);
  void util_rand_init(std::size_t baseind);
  BaseType maxabsdiff(QbitRegister &x, Type sfactor = Type(1.0, 0.));
  BaseType maxl2normdiff(QbitRegister &x);
  void dumpbin(std::string fn);
  double entropy();
  std::vector<double> googleStats();
  void normalize();
  BaseType computenorm();
  Type compute_overlap( QbitRegister<Type> &psi );

  void Print(std::string x, std::vector<std::size_t> qbits = {});

  double HP_Distrpair(unsigned pos, TM2x2<Type> const&m);
  double HP_Distrpair(unsigned control, unsigned qubit, TM2x2<Type> const&m);

  // Members
  std::size_t nqbits;
  std::vector<Type, openqu::AlignedAllocator<Type, 256>> statestorage;
  Type *state;
  Permutation *permutation;
  Timer *timer;
  std::size_t llc_watermarkbit;
  bool importedstate;
  bool specialize;

  // temporary buffer for fusion
  bool fusion;
  unsigned log2llc;
  std::vector<std::tuple<std::string, TM2x2<Type>, unsigned, unsigned>> fwindow;

 private:
  std::size_t localsize_;
  std::size_t globalsize_;
  std::size_t tmpspacesize_;
};

template <typename Type>
using BaseType = typename QbitRegister<Type>::BaseType;

//
// Derived class of QbitRegister that allows measurement of qubit gate depth.
//
#include "QbitRegisterMetric.hpp"

//
// Derived class of QbitRegister that automatically implements noise gates.
//
#include "NoisyQureg.hpp"

