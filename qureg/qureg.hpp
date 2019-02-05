//------------------------------------------------------------------------------
// Copyright (C) 2017 Intel Corporation 
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
//------------------------------------------------------------------------------

#pragma once

#include "permute.hpp"
#include "../util/utils.hpp"
#include "../util/timer.hpp"

#include "../util/alignedallocator.hpp"
#include "../util/mpi.hpp"
#include "../util/bitops.hpp"
#include "../util/conversion.hpp"
#include "../util/tinymatrix.hpp"

#include <algorithm>	// for std::swap
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>	// for std::swap
#include <numeric>
#include <tuple>
#include <vector>

#if defined(__ICC) || defined(__INTEL_COMPILER)
#include <mkl.h>
#endif

/// \addtogroup qureg
/// @{

/// @file qureg.hpp
/// @brief Declare the @c QubitRegister class.

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

/////////////////////////////////////////////////////////////////////////////////////////
// Qubitregister class declaration
/////////////////////////////////////////////////////////////////////////////////////////

template <class Type = ComplexDP>
class QubitRegister
{
  public:
    using value_type = Type;
    typedef typename extract_value_type<Type>::value_type BaseType;

  // constructors / destructors
    QubitRegister();
    QubitRegister(std::size_t num_qubits, std::string style = "", 
                 std::size_t base_index = 0, std::size_t tmp_spacesize_ = 0);
    QubitRegister(const QubitRegister &in);
    QubitRegister(std::size_t num_qubits, Type *state, std::size_t tmp_spacesize_ = 0);
    ~QubitRegister();

  // allocation & initialization
  void AllocateAdditionalQubit();
  void Allocate(std::size_t new_num_qubits, std::size_t tmp_spacesize_);
  void Initialize(std::size_t new_num_qubits, std::size_t tmp_spacesize_);

  inline Type& operator[] (std::size_t index) { return state[index]; }
  inline Type& operator[] (std::size_t index) const { return state[index]; }

  std::size_t LocalSize() const { return local_size_; }
  std::size_t GlobalSize() const { return global_size_; }

  void Resize(std::size_t new_num_amplitudes);
  std::size_t size() const { return global_size_; }
  unsigned NumQubits() const { return num_qubits; }
  Type *TmpSpace() const { return state + LocalSize(); }
  size_t TmpSize() const {return tmp_spacesize_;}

  // bit manipulation
  inline bool check_bit(std::size_t variable, std::size_t position) const
  {
    std::size_t one = (std::size_t)1, position_long = UL(position);
    return variable & (one << position_long);
  }
  inline std::size_t set_bit(std::size_t variable, std::size_t position) const
  {
    std::size_t one = (std::size_t)1, position_long = UL(position);
    return variable | (one << position_long);
  }
  inline std::size_t clear_bit(std::size_t variable, std::size_t position) const
  {
     return (variable & ~(UL(1) << UL(position)));
  }

  void EnableStatistics();
  void GetStatistics();
  void ResetStatistics();

  void Permute(std::vector<std::size_t> permutation_new_vec);


  // Generic gates
  // single qubit gates
  bool Apply1QubitGate_helper(unsigned qubit,  TM2x2<Type> const&m,
                              std::size_t sstate_ind, std::size_t estate_ind);
  void Apply1QubitGate(unsigned qubit, TM2x2<Type> const&m);
  // constrolled gates
  bool ApplyControlled1QubitGate_helper(unsigned control, unsigned qubit, TM2x2<Type> const&m,
                                        std::size_t sind, std::size_t eind);
  void ApplyControlled1QubitGate(unsigned control, unsigned qubit, TM2x2<Type> const&m);
  // swap gates
  bool ApplySwap_helper(unsigned qubit1, unsigned qubit2, TM2x2<Type> const&m);
  void ApplySwap(unsigned qubit1, unsigned qubit2);
  void ApplyISwap(unsigned qubit1, unsigned qubit2);
  void Apply4thRootISwap(unsigned qubit1, unsigned qubit2);
  void ApplySqrtISwap(unsigned qubit1, unsigned qubit2);
  void ApplyISwapRotation(unsigned qubit1, unsigned qubit2, TM2x2<Type> const&m);
  void Swap(unsigned b1, unsigned b2);
  // diagonal gates
  void ApplyDiagSimp(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  void ApplyDiag(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  void ApplyDiagControl(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  void ApplyDiagGeneral(unsigned qubit1, unsigned qubit2, TM4x4<Type> const&m);
  // two-qubit gates
  void Apply2QubitGate(unsigned const qubit_high, unsigned const qubit_low, TM4x4<Type> const&m);
  // specialized gates
  void ApplyRotationX(unsigned const qubit, BaseType theta);
  void ApplyRotationY(unsigned const qubit, BaseType theta);
  void ApplyRotationZ(unsigned const qubit, BaseType theta);
  void ApplyPauliX(unsigned const qubit);
  void ApplyPauliY(unsigned const qubit);
  void ApplyPauliZ(unsigned const qubit);
  void ApplyPauliSqrtX(unsigned const qubit);
  void ApplyPauliSqrtY(unsigned const qubit);
  void ApplyPauliSqrtZ(unsigned const qubit);
  void ApplyT(unsigned const qubit);
  void ApplyToffoli(unsigned const qubit1, unsigned const qubit2, unsigned const qubit3);
  void ApplyHadamard(unsigned const qubit);

  void ApplyCRotationX(unsigned const qubit, unsigned const qubit2, BaseType theta);
  void ApplyCRotationY(unsigned const qubit, unsigned const qubit2, BaseType theta);
  void ApplyCRotationZ(unsigned const qubit, unsigned const qubit2, BaseType theta);
  void ApplyCPauliX(unsigned const qubit1, unsigned const qubit2);
  void ApplyCPauliY(unsigned const qubit1, unsigned const qubit2);
  void ApplyCPauliZ(unsigned const qubit1, unsigned const qubit2);
  void ApplyCPauliSqrtZ(unsigned const qubit1, unsigned const qubit2);
  void ApplyCHadamard(unsigned const qubit1, unsigned const qubit2);

  // fusion  
  void TurnOnFusion(unsigned log2llc = 20);
  void TurnOffFusion();
  bool IsFusionEnabled();
  void ApplyFusedGates();

  // gate specialization (experimental)
  void TurnOnSpecialize();
  void TurnOffSpecialize();

  // measurement
  bool GetClassicalValue(unsigned qubit, BaseType tolerance = 1.e-13) const;
  bool IsClassicalBit(unsigned qubit, BaseType tolerance = 1.e-13) const;
  void CollapseQubit(unsigned qubit, bool value);
  BaseType GetProbability(unsigned qubit);

  // expectation values without state update
  void ExpectationValueX(unsigned const qubit, BaseType &sum, BaseType coeff=1.);
  void ExpectationValueY(unsigned const qubit, BaseType &sum, BaseType coeff=1.);
  void ExpectationValueZ(unsigned const qubit, BaseType &sum, BaseType coeff=1.);
  void ExpectationValueXX(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueXY(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueXZ(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueYX(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueYY(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueYZ(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueZX(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueZY(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValueZZ(unsigned const qubit, unsigned const qubit2,
                          BaseType &sum, BaseType coeff=1.);
  void ExpectationValue(std::vector<unsigned> &qubits, std::vector<unsigned> &observables,
                        BaseType &sum, BaseType coeff=1.);

  // utilities
  bool operator==(const QubitRegister &rhs);
  void Initialize(std::string style, std::size_t base_index);
  void RandomInitialize(std::size_t base_index);
  BaseType maxabsdiff(QubitRegister &x, Type sfactor = Type(1.0, 0.));
  BaseType maxl2normdiff(QubitRegister &x);
  void dumpbin(std::string fn);
  double Entropy();
  std::vector<double> GoogleStats();
  void Normalize();
  BaseType ComputeNorm();
  Type ComputeOverlap( QubitRegister<Type> &psi );

  void Print(std::string x, std::vector<std::size_t> qbits = {});

  double HP_Distrpair(unsigned pos, TM2x2<Type> const&m);
  double HP_Distrpair(unsigned control, unsigned qubit, TM2x2<Type> const&m);

  // Members
  std::size_t num_qubits;
  std::vector<Type, openqu::AlignedAllocator<Type, 256>> state_storage;
  Type *state;
  Permutation *permutation;
  Timer *timer;
  std::size_t llc_watermarkbit;
  bool imported_state;
  bool specialize;

  // temporary buffer for fusion
  bool fusion;
  unsigned log2llc;
  std::vector<std::tuple<std::string, TM2x2<Type>, unsigned, unsigned>> fwindow;

 private:
  std::size_t local_size_;
  std::size_t global_size_;
  std::size_t tmp_spacesize_;
};

template <typename Type>
using BaseType = typename QubitRegister<Type>::BaseType;

/// @}

//
// Derived class of QubitRegister that allows measurement of qubit gate depth.
//
#include "QubitRegisterMetric.hpp"
//
// Derived class of QubitRegister that automatically implements noise gates.
//
#include "NoisyQureg.hpp"

