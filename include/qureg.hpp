#pragma once

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

#ifdef USE_MKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "permute.hpp"
// utility files
#include "utils.hpp"
#include "mpi_utils.hpp"
#include "rng_utils.hpp"
#include "timer.hpp"
#include "gate_counter.hpp"
#include "alignedallocator.hpp"
#include "mpi_env.hpp"
#include "bitops.hpp"
#include "conversion.hpp"
#include "tinymatrix.hpp"

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
using TM2x2 = qhipster::TinyMatrix<Type, 2, 2, 32>;


template<class Type>
using TM4x4 = qhipster::TinyMatrix<Type, 4, 4, 32>;

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
  // The 'style' of initialization can be:
  // - 'rand': real and imag part of each amplitudes are uniformly random,
  //           using either the **local** or **pool** RNG stream,
  //           then state is normalized.
  // - 'base': state of the computational basis, only a non-zero amplitude.
  // - '++++': the balanced superposition of all computational basis states.
  void Initialize(std::string style, std::size_t base_index);

  // overload [] operator to return the amplitude stored at the local index.
  inline Type& operator[] (std::size_t index) { return state[index]; }
  inline Type& operator[] (std::size_t index) const { return state[index]; }
  // get the amplitude corresponding to a global index (with MPI broadcast).
  Type GetGlobalAmplitude(std::size_t index) const;

  std::size_t LocalSize() const { return local_size_; }
  std::size_t GlobalSize() const { return global_size_; }

  void Resize(std::size_t new_num_amplitudes);
  std::size_t size() const { return global_size_; }
  std::size_t NumQubits() const { return num_qubits; }
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
  void DisableStatistics();
  void ResetStatistics();

  void Permute(std::vector<std::size_t> permutation_new_vec);


  // Generic gates
  // single qubit gates
  bool Apply1QubitGate_helper(unsigned qubit,  TM2x2<Type> const&m,
                              std::size_t sstate_ind, std::size_t estate_ind);
  void Apply1QubitGate(unsigned qubit, TM2x2<Type> const&m);
  // constrolled gates
  bool ApplyControlled1QubitGate_helper(unsigned control_qubit, unsigned target_qubit,
                                        TM2x2<Type> const&m,
                                        std::size_t sind, std::size_t eind);
  void ApplyControlled1QubitGate(unsigned control_qubit, unsigned target_qubit,
                                 TM2x2<Type> const&m);
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

  void ApplyCRotationX(unsigned const control_qubit, unsigned const target_qubit,
                       BaseType theta);
  void ApplyCRotationY(unsigned const control_qubit, unsigned const target_qubit,
                       BaseType theta);
  void ApplyCRotationZ(unsigned const control_qubit, unsigned const target_qubit,
                       BaseType theta);
  void ApplyCPauliX(unsigned const control_qubit, unsigned const target_qubit);
  void ApplyCPauliY(unsigned const control_qubit, unsigned const target_qubit);
  void ApplyCPauliZ(unsigned const control_qubit, unsigned const target_qubit);
  void ApplyCPauliSqrtZ(unsigned const control_qubit, unsigned const target_qubit);
  void ApplyCHadamard(unsigned const control_qubit, unsigned const target_qubit);

  void ApplyCPhaseRotation(unsigned const qubit, unsigned const qubit2, BaseType theta);

  // Algorithms
  void ApplyNCU(TM2x2<Type> gate, const std::vector<std::size_t>& ctrl_indices, const std::vector<std::size_t>& aux_indices, unsigned const target);

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
  BaseType ExpectationValueX(unsigned const qubit, BaseType coeff=1.);
  BaseType ExpectationValueY(unsigned const qubit, BaseType coeff=1.);
  BaseType ExpectationValueZ(unsigned const qubit, BaseType coeff=1.);
  BaseType ExpectationValueXX(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueXY(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueXZ(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueYX(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueYY(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueYZ(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueZX(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueZY(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValueZZ(unsigned const qubit, unsigned const qubit2,
                              BaseType coeff=1.);
  BaseType ExpectationValue(std::vector<unsigned> &qubits,
                            std::vector<unsigned> &observables,
                            BaseType coeff=1.);

  // noisy simulation:
  BaseType GetT1   () {return T_1_;   }
  BaseType GetT2   () {return T_2_;   }
  BaseType GetTphi () {return T_phi_; }
  void SetNoiseTimescales(BaseType T1, BaseType T2);
  void ApplyNoiseGate(const unsigned qubit, const BaseType duration);
//FIXME DELETE  BaseType IncoherentAverageOverAllStatesOfPool (BaseType local_value);


  // Utilities:
  bool operator==(const QubitRegister &rhs);
  BaseType MaxAbsDiff(QubitRegister &x, Type sfactor = Type(1.0, 0.));
  BaseType MaxL2NormDiff(QubitRegister &x);
  void dumpbin(std::string fn);
  double Entropy();
  std::vector<double> GoogleStats();
  void Normalize();
  BaseType ComputeNorm();
  Type ComputeOverlap( QubitRegister<Type> &psi );

  void Print(std::string x, std::vector<std::size_t> qbits = {});

  double HP_Distrpair(unsigned pos, TM2x2<Type> const&m);
  double HP_Distrpair(unsigned control, unsigned qubit, TM2x2<Type> const&m);

  // related to the internal random number generator.
  qhipster::RandomNumberGenerator<BaseType> * GetRngPtr () {return rng_ptr_; }
  void ResetRngPtr () {rng_ptr_=nullptr; }
  void SetRngPtr (qhipster::RandomNumberGenerator<BaseType> * rng_ptr) {rng_ptr_=rng_ptr; }
  void SetSeedRngPtr (std::size_t seed)
  {assert(rng_ptr_); rng_ptr_->SetSeedStreamPtrs(seed); }

  // Members
  std::size_t num_qubits;
  std::vector<Type, qhipster::AlignedAllocator<Type, 256>> state_storage;
  Type *state;
  Permutation *permutation;
  Timer *timer;
  GateCounter *gate_counter;
  std::size_t llc_watermarkbit;
  bool imported_state;
  bool specialize;

  // temporary buffer for fusion
  bool fusion;
  unsigned log2llc;
  std::vector<std::tuple<std::string, TM2x2<Type>, unsigned, unsigned>> fwindow;

  // set option of printing more info.
  static void SetDoPrintExtraInfo( bool value )
  { do_print_extra_info = value; }

 private:
  std::size_t local_size_;
  std::size_t global_size_;
  std::size_t tmp_spacesize_;
  static bool do_print_extra_info;

  qhipster::RandomNumberGenerator<BaseType> * rng_ptr_ = nullptr;
  BaseType T_1_;	// T_1   given in terms of the chosen time unit
  BaseType T_2_;	// T_2   given in terms of the chosen time unit
  BaseType T_phi_;	// T_phi given in terms of the chosen time unit

  private:
    QubitRegister<Type>& operator=(const QubitRegister<Type>& src) { return *this; }
};

template <typename Type>
bool QubitRegister<Type>::do_print_extra_info = false;

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

