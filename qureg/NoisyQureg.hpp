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
#include <vector>
#include <random>

/// \addtogroup NoisyQureg
/// @{

/// @file NoisyQureg.h
///  @brief Define the class @c NoisyQureg that iautomatically includes noise gates.

using namespace std;

/// @brief Class that expand @c QitRegister states by adding noise between "logical" gates.
///
/// @param Nqubit number of qubits

template <class Type = ComplexDP>
class NoisyQureg : public QbitRegister<Type> {

  private :

    typedef typename QbitRegister<Type>::BaseType BaseType;

    std::vector<BaseType> TimeFromLastGate;		// given in terms of the chosen time unit
    BaseType TimeOneQubitLogicalGate = 1.;
    BaseType TimeTwoQubitLogicalGate = 1.;

    unsigned OneQubitLogicalGateCount=0;
    unsigned TwoQubitLogicalGateCount=0;

    // ~~~~~~~~~~ parameters characterizing the noise model ~~~~~~~~~
    BaseType T_1 = 1000;				// T_1   given in terms of the chosen time unit
    BaseType T_2 = 100;					// T_2   given in terms of the chosen time unit
    BaseType T_phi = 1./( 1./T_2 - 1./(2.*T_1) );	// T_phi given in terms of the chosen time unit

    // Pseudo random-number-generator to sample from normal distribution
    // (for the rotation angles of the noise gates)
    std::default_random_engine generator;
    std::normal_distribution<BaseType> gaussian_RNG{0.,1.} ;
    

  public :

    /// Constructor
    NoisyQureg<Type>(unsigned Nqubits , unsigned RNG_seed=12345 ,
                     BaseType T1=2000 , BaseType T2=1000 ) :
        QbitRegister<Type>(Nqubits)
    {
      T_1 = T1;
      T_2 = T2;
      TimeFromLastGate.assign(Nqubits,0.);

      // Initialization of the seed for the generation of the noise gate parameters
      generator.seed( RNG_seed );
//      srand48( RNG_seed );
    }

    /// Default destructor
    ~NoisyQureg<Type>() {};

    // Utilities
    void reset_time_for_all_qubits();
    void apply_noise_gates_on_all_qubits();

    // Noise model
    void update_noise_model(BaseType, BaseType);

    // Get stats
    unsigned GetTotalQubitLogicalGateCount();
    unsigned GetOneQubitLogicalGateCount();
    unsigned GetTwoQubitLogicalGateCount();

    // Reset decoherence time and implement noise gate
    void AddNoiseOneQubitGate(unsigned const);
    void AddNoiseTwoQubitGate(unsigned const, unsigned const);
    void NoiseGate(unsigned const);

    // Perform gates
    void apply1QubitGate(unsigned const, openqu::TinyMatrix<Type, 2, 2, 32>);
    void applyHadamard(unsigned const);
    void applyRotationX(unsigned const, BaseType);
    void applyRotationY(unsigned const, BaseType);
    void applyRotationZ(unsigned const, BaseType);
    void applyCPauliX(unsigned const, unsigned const);
    void applyControlled1QubitGate(unsigned const, unsigned const, openqu::TinyMatrix<Type, 2, 2, 32>);

};


// ----------------------------------------------------------
// ----------------------------------------------------------
// ----------------------------------------------------------

/// Reset to zero the time elapsed for each and every qubit in the register.
template <class Type>
void NoisyQureg<Type>::reset_time_for_all_qubits()
{
  unsigned Nqubits = this->nqbits;
  // increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<Nqubits; q++)
      TimeFromLastGate[q] = 0. ;
}


/// @brief Apply the noise gates on each and every qubit. Then reset to time counter.
///
/// This is useful, for example, at the end of a circuit before measuring the
/// quantities of interest: One has to apply the noise corresponding to the idle
/// evolution between the last logical gate and the final time.
/// The time from last logical gate is then resetted to zero for every qubit.
template <class Type>
void NoisyQureg<Type>::apply_noise_gates_on_all_qubits()
{
  unsigned Nqubits = this->nqbits;
  // increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<Nqubits; q++)
      NoiseGate(q);
  reset_time_for_all_qubits();
}


/// Update the T_1 and T_2 values in accordance to the new noise model.
template <class Type>
void NoisyQureg<Type>::update_noise_model(BaseType T1, BaseType T2)
{
  T_1 = T1;
  T_2 = T2;
}


/// Return the current number of (logical) single-qubit gates.
template <class Type>
unsigned NoisyQureg<Type>::GetOneQubitLogicalGateCount()
{
  return OneQubitLogicalGateCount;
}


/// Return the current number of (logical) two-qubit gates.
template <class Type>
unsigned NoisyQureg<Type>::GetTwoQubitLogicalGateCount()
{
  return TwoQubitLogicalGateCount;
}


/// Return the current number of (logical) 1- and 2-qubit gates.
template <class Type>
unsigned NoisyQureg<Type>::GetTotalQubitLogicalGateCount()
{
  return OneQubitLogicalGateCount + TwoQubitLogicalGateCount;
}


/// Include and execute the noise gate corresponding to the idle time of a single qubit.
template <class Type>
void NoisyQureg<Type>::AddNoiseOneQubitGate(unsigned const qubit)
{
  unsigned Nqubits = this->nqbits;
  // Implement the noise gate
  NoiseGate(qubit);
  // Increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<Nqubits; q++)
      TimeFromLastGate[q] += TimeOneQubitLogicalGate ;
  // Reset the time elapsed from last logical gate on the specific qubit
  TimeFromLastGate[qubit]=0.;
  // Update counter for (logical) one-qubit gates
  OneQubitLogicalGateCount++;
}


/// Include and execute the noise gate corresponding to the idle time of two qubits.
template <class Type>
void NoisyQureg<Type>::AddNoiseTwoQubitGate(unsigned const q1, unsigned const q2)
{
  unsigned Nqubits = this->nqbits;
  // Implement the noise gate
  NoiseGate(q1);
  NoiseGate(q2);
  // Increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<Nqubits; q++)
      TimeFromLastGate[q] += TimeTwoQubitLogicalGate ;
  // Reset the time elapsed from last logical gate on the specific qubits
  TimeFromLastGate[q1]=0.;
  TimeFromLastGate[q2]=0.;
  // Update counter for (logical) two-qubit gates
  TwoQubitLogicalGateCount++;
}


/// @grief Noise gate corresponding to single-qubit rotation with appropriate (stochastic) angle.
///
/// To obtain a single rotation around an arbitrary axis we use the relations:
///      | a b c |               | h-f |
///  R = | d e f |    -->    u = | c-g |    ---> abs(u) = 2 sin( 'angle' )
///      | g h i |               | d-b |    ---> u/abs(u) = rotation axis  
template <class Type>
void NoisyQureg<Type>::NoiseGate(unsigned const qubit )
{
  BaseType t = TimeFromLastGate[qubit] ;
  if (t==0) return;

  BaseType p_X , p_Y , p_Z ;
  p_X = (1. - std::exp(-t/T_1) )/4.;
  p_Y = (1. - std::exp(-t/T_1) )/4.;
  p_Z = (1. - std::exp(-t/T_2) )/2. + (1. - std::exp(-t/T_1) )/4.;
  assert( p_X>0 && p_Y>0 && p_Z>0 );

  // Computation of the standard deviations for the noise gate parameters
  BaseType s_X , s_Y , s_Z ;
  s_X = std::sqrt( -std::log(1.-p_X) );
  s_Y = std::sqrt( -std::log(1.-p_Y) );
  s_Z = std::sqrt( -std::log(1.-p_Z) );

  // Generate angle and direction of Pauli rotation for Pauli-twirld noise channel
  BaseType v_X , v_Y , v_Z;
  v_X = gaussian_RNG(generator) * s_X /2.;
  v_Y = gaussian_RNG(generator) * s_Y /2.;
  v_Z = gaussian_RNG(generator) * s_Z /2.;

  // Compose the 3-dimensional rotation: R_X(v_X).R_Y(v_Y).R_Z(v_Z)
  //       |  cos Y cos Z                          -cos Y sin Z                           sin Y        |
  //   R = |  sin X sin Y cos Z + cos X sin Z      -sin X sin Y sin Z + cos X cos Z      -sin X cos Y  |
  //       | -cos X sin Y cos Z + sin X sin Z       cos X sin Y sin Z + sin X cos Z       cos X cos Y  |
  //
  //       | cos X sin Y sin Z + sin X cos Z + sin X cos Y |
  //   u = |    sin Y + cos X sin Y cos Z - sin X sin Z    |
  //       | sin X sin Y cos Z + cos X sin Z + cos Y sin Z |
  std::vector<BaseType> u = { std::cos(v_X) * std::sin(v_Y) * std::sin(v_Z) +
                              std::sin(v_X) * std::cos(v_Z) +
                              std::sin(v_X) * std::cos(v_Y),
                                 std::sin(v_Y) +
                                 std::cos(v_X) * std::sin(v_Y) * std::cos(v_Z) -
                                 std::sin(v_X) * std::sin(v_Z),
                                    std::sin(v_X) * std::sin(v_Y) * std::cos(v_Z) +
                                    std::cos(v_X) * std::sin(v_Z) +
                                    std::cos(v_Y) * std::sin(v_Z)                   };  
  BaseType norm_u = std::sqrt( std::norm(u[0]) + std::norm(u[1]) + std::norm(u[2]) );
  std::vector<BaseType> axis = { u[0]/norm_u , u[1]/norm_u , u[2]/norm_u };
  // Compute the angle of rotation
  BaseType trace_R =  std::cos(v_Y) * std::cos(v_Z)
                    - std::sin(v_X) * std::sin(v_Y) * std::sin(v_Z)
                    + std::cos(v_X) * std::cos(v_Z)
                    + std::cos(v_X) * std::cos(v_Y) ;
  BaseType angle = std::acos( (trace_R-1.)/2. );
  if (false) std::cout << " angle of rotation = " << angle << "\n";
  // Alternative expression for the angle, from Wikipedia.
  if (false) std::cout << " angle of rotation (alternative formula) = " << std::asin( norm_u/2. ) << "\n";


  // Costruct the 1/2-spin matrix corresponding to the axis above
  //    sigma_axis
  // and use it to implement the single-qubit rotation :
  //    rot = exp( i * angle * sigma_axis )
  openqu::TinyMatrix<Type, 2, 2, 32> rot;
  BaseType s(std::sin(angle/2.)) , c(std::cos(angle/2.)) ;
  rot(0, 0) = Type(  c         ,  s*axis[2] );
  rot(0, 1) = Type(  s*axis[1] ,  s*axis[0] );
  rot(0, 1) = Type( -s*axis[1] ,  s*axis[0] );
  rot(1, 1) = Type(  c         , -s*axis[2] );

  // Apply the noise gate
  QbitRegister<Type>::apply1QubitGate(qubit,rot);
}


// ----------------------------------------------------------
// -------- list of modified gates --------------------------
// ----------------------------------------------------------


template <class Type>
void NoisyQureg<Type>::apply1QubitGate(unsigned const q,
                                       openqu::TinyMatrix<Type, 2, 2, 32> V)
{
  AddNoiseOneQubitGate(q); 
  QbitRegister<Type>::apply1QubitGate(q,V);
}

template <class Type>
void NoisyQureg<Type>::applyHadamard(unsigned const q)
{
  AddNoiseOneQubitGate(q); 
  QbitRegister<Type>::applyHadamard(q);
}

template <class Type>
void NoisyQureg<Type>::applyRotationX(unsigned const q, BaseType theta)
{
  AddNoiseOneQubitGate(q); 
  QbitRegister<Type>::applyRotationX(q,theta);
}

template <class Type>
void NoisyQureg<Type>::applyRotationY(unsigned const q, BaseType theta)
{
  AddNoiseOneQubitGate(q); 
  QbitRegister<Type>::applyRotationY(q,theta);
}

template <class Type>
void NoisyQureg<Type>::applyRotationZ(unsigned const q, BaseType theta)
{
  AddNoiseOneQubitGate(q); 
  QbitRegister<Type>::applyRotationZ(q,theta);
}

template <class Type>
void NoisyQureg<Type>::applyCPauliX(unsigned const q1, unsigned const q2)
{
  AddNoiseTwoQubitGate(q1,q2); 
  QbitRegister<Type>::applyCPauliX(q1,q2);
}

template <class Type>
void NoisyQureg<Type>::applyControlled1QubitGate(unsigned const q1, unsigned const q2,
                                                 openqu::TinyMatrix<Type, 2, 2, 32> V)
{
  AddNoiseTwoQubitGate(q1,q2); 
  QbitRegister<Type>::applyControlled1QubitGate(q1,q2,V);
}

/// @}*/
