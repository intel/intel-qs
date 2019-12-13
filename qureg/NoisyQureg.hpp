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

#ifndef NOISY_QUREG_H
#define NOISY_QUREG_H

#include <random>
#include <vector>

/// \addtogroup NoisyQureg
/// @{

/// @file NoisyQureg.h
/// @brief Define the class @c NoisyQureg that automatically includes noise gates.

using namespace std;

/// @brief Class that expand @c QubitRegister states by adding noise between "logical" gates.
///
/// @param num_qubit is the number of qubits
/// When we refer to ``experimental'' gates, it means that noise gates are excluded.
/// For the simulation to be faithful (i.e. with one-to-one correspondence with the experimental gates),
/// we include among the experimental gates both the gates for the algorithm and those to schedule it
/// according to the connectivity of the specific hardware.

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type = ComplexDP>
class NoisyQureg : public QubitRegister<Type>
{
  private :

    typedef typename QubitRegister<Type>::BaseType BaseType;

    std::vector<BaseType> time_from_last_gate;	// given in terms of the chosen time unit
    BaseType time_one_qubit_experimental_gate = 1.;
    BaseType time_two_qubit_experimental_gate = 1.;

    // matrix with number of gates between qubits (diagonal is one-qubit gates)
    std::vector<unsigned> experimental_gate_count_matrix;
    unsigned one_qubit_experimental_gate_count=0;
    unsigned two_qubit_experimental_gate_count=0;

    // ~~~~~~~~~~ parameters characterizing the noise model ~~~~~~~~~
    BaseType T_1 = 1000;		// T_1   given in terms of the chosen time unit
    BaseType T_2 = 100;			// T_2   given in terms of the chosen time unit
    BaseType T_phi = 1./( 1./T_2 - 1./(2.*T_1) );	// T_phi given in terms of the chosen time unit

    // Pseudo random-number-generator to sample from normal distribution
    // (for the rotation angles of the noise gates)
    std::default_random_engine generator;
    std::normal_distribution<BaseType> gaussian_RNG{0.,1.} ;
    

  public :

    /// Constructor
    NoisyQureg<Type>(unsigned num_qubits , unsigned RNG_seed=12345 ,
                     BaseType T1=2000 , BaseType T2=1000 ) :
        QubitRegister<Type>(num_qubits)
    {
      T_1 = T1;
      T_2 = T2;
      time_from_last_gate.assign(num_qubits,0.);
      experimental_gate_count_matrix.assign(num_qubits*num_qubits,0);

      // Initialization of the seed for the generation of the noise gate parameters
      generator.seed( RNG_seed );
//      srand48( RNG_seed );
    }

    /// Default destructor
    ~NoisyQureg<Type>() {};

    // Initialization of the state (without noise)
    void Initialize(std::string style, std::size_t base_index);

    // Utilities
    void ResetTimeForAllQubits();
    void ApplyNoiseGatesOnAllQubits();

    // Noise model
    void SetDecoherenceTime(BaseType, BaseType);
    void SetGateDurations(BaseType, BaseType);

    // Get stats
    unsigned GetTotalExperimentalGateCount();
    unsigned GetOneQubitExperimentalGateCount();
    unsigned GetTwoQubitExperimentalGateCount();
    std::vector<unsigned> GetExperimentalGateCount(unsigned q1);
    unsigned GetExperimentalGateCount(unsigned q1, unsigned q2);

    // Reset decoherence time and implement noise gate
    void AddNoiseOneQubitGate(unsigned const);
    void AddNoiseTwoQubitGate(unsigned const, unsigned const);
    void NoiseGate(unsigned const);
    // Old funcation that can be deleted and should NOT be used:
    void NoiseGate_OLD(unsigned const);

    // Perform gates
    void Apply1QubitGate(unsigned const, qhipster::TinyMatrix<Type, 2, 2, 32>);
    void ApplyHadamard(unsigned const);
    void ApplyRotationX(unsigned const, BaseType);
    void ApplyRotationY(unsigned const, BaseType);
    void ApplyRotationZ(unsigned const, BaseType);
    void ApplyCPauliX(unsigned const, unsigned const);
    void ApplyControlled1QubitGate(unsigned const, unsigned const, qhipster::TinyMatrix<Type, 2, 2, 32>);

};


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


// ------------ initialize state ---------------------------------------------------------

template <class Type>
void NoisyQureg<Type>::Initialize(std::string style, std::size_t base_index)
{
  one_qubit_experimental_gate_count = 0;
  two_qubit_experimental_gate_count = 0;
  ResetTimeForAllQubits();
  QubitRegister<Type>::Initialize(style,base_index);
}

// ------------ utils --------------------------------------------------------------------

/// Reset to zero the time elapsed for each and every qubit in the register.
template <class Type>
void NoisyQureg<Type>::ResetTimeForAllQubits()
{
  unsigned num_qubits = this->num_qubits;
  // increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<num_qubits; q++)
      time_from_last_gate[q] = 0. ;
}


/// @brief Apply the noise gates on each and every qubit. Then reset to time counter.
///
/// This is useful, for example, at the end of a circuit before measuring the
/// quantities of interest: One has to apply the noise corresponding to the idle
/// evolution between the last logical gate and the final time.
/// The time from last logical gate is then resetted to zero for every qubit.
template <class Type>
void NoisyQureg<Type>::ApplyNoiseGatesOnAllQubits()
{
  unsigned num_qubits = this->num_qubits;
  // increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<num_qubits; q++)
      NoiseGate(q);
  ResetTimeForAllQubits();
}

/// Set the decoherence time in terms of T_1 and T_2 values (in accordance to the new noise model).
template <class Type>
void NoisyQureg<Type>::SetDecoherenceTime(BaseType T1, BaseType T2)
{
  T_1 = T1;
  T_2 = T2;
}

/// Update the duration of single- and two- qubit gates.
template <class Type>
void NoisyQureg<Type>::SetGateDurations(BaseType Ts, BaseType Td)
{
  time_one_qubit_experimental_gate = Ts;
  time_two_qubit_experimental_gate = Td;
}

// ------------ count the experimental gates ---------------------------------------------

/// Return the current number of (experimental) single-qubit gates.
template <class Type>
unsigned NoisyQureg<Type>::GetOneQubitExperimentalGateCount()
{ return one_qubit_experimental_gate_count; }


/// Return the current number of (experimental) two-qubit gates.
template <class Type>
unsigned NoisyQureg<Type>::GetTwoQubitExperimentalGateCount()
{ return two_qubit_experimental_gate_count; }


/// Return the current number of (experimental) 1- and 2-qubit gates.
template <class Type>
unsigned NoisyQureg<Type>::GetTotalExperimentalGateCount()
{ return one_qubit_experimental_gate_count + two_qubit_experimental_gate_count; }

/// Return the number of (experimental) gates involving qubit q.
template <class Type>
std::vector<unsigned> NoisyQureg<Type>::GetExperimentalGateCount(unsigned q)
{
  std::vector<unsigned> SingleRowOfMatrix (experimental_gate_count_matrix.begin()+ q * this->num_qubits,
                                           experimental_gate_count_matrix.begin()+ (q+1) * this->num_qubits);
  return SingleRowOfMatrix;
}

/// Return the number of (experimental) gates involving qubits q1,q2.
template <class Type>
unsigned NoisyQureg<Type>::GetExperimentalGateCount(unsigned q1, unsigned q2)
{
  return experimental_gate_count_matrix[ q1 * this->num_qubits + q2 ];
}

// ------------ execute the noise gates --------------------------------------------------

/// Include and execute the noise gate corresponding to the idle time of a single qubit.
template <class Type>
void NoisyQureg<Type>::AddNoiseOneQubitGate(unsigned const qubit)
{
  unsigned num_qubits = this->num_qubits;
  // Implement the noise gate
  NoiseGate(qubit);
  // Increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<num_qubits; q++)
      time_from_last_gate[q] += time_one_qubit_experimental_gate ;
  // Reset the time elapsed from last logical gate on the specific qubit
  time_from_last_gate[qubit]=0.;
  // Update counter for (logical) one-qubit gates
  one_qubit_experimental_gate_count++;
  experimental_gate_count_matrix[qubit*num_qubits+qubit]++;
}


/// Include and execute the noise gate corresponding to the idle time of two qubits.
template <class Type>
void NoisyQureg<Type>::AddNoiseTwoQubitGate(unsigned const q1, unsigned const q2)
{
  unsigned num_qubits = this->num_qubits;
  // Implement the noise gate
  NoiseGate(q1);
  NoiseGate(q2);
  // Increase the idle time for all the qubits (TODO no gate parallelism is assumed)
  for (unsigned q = 0; q<num_qubits; q++)
      time_from_last_gate[q] += time_two_qubit_experimental_gate ;
  // Reset the time elapsed from last logical gate on the specific qubits
  time_from_last_gate[q1]=0.;
  time_from_last_gate[q2]=0.;
  // Update counter for (logical) two-qubit gates
  two_qubit_experimental_gate_count++;
  experimental_gate_count_matrix[q1*num_qubits+q2]++;
  experimental_gate_count_matrix[q2*num_qubits+q1]++;
}


/// @brief Noise gate corresponding to single-qubit rotation with appropriate (stochastic) angle.
///
/// Each noise gate is the product of three rotations around X,Y,Z axis (by a small angle each).
/// We compute their product before applying it to the quantum register.
template <class Type>
void NoisyQureg<Type>::NoiseGate(unsigned const qubit )
{
  BaseType t = time_from_last_gate[qubit] ;
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

  // Direct construction of the 2x2 matrix corresponding to the noise gate
  //     U_noise = exp(-i v_X X) * exp(-i v_Y Y) * exp(-i v_Z Z)
  // Helpful quantities:
  //     (A) = cos v_z -i sin v_z
  //     (B) = cos v_x * cos v_Y -i sin v_X * sin v_Y
  //     (C) = cos v_x * sin v_Y -i sin v_X * cos v_Y
  // Then :
  //               | A*B   -A'*C' |
  //     U_noise = |              |
  //               | A*C    A'*B' |

  Type A , B , C ;
  A = { std::cos(v_Z) , -std::sin(v_Z) };
  B = { std::cos(v_X)*std::cos(v_Y) , -std::sin(v_X)*std::sin(v_Y) };
  C = { std::cos(v_X)*std::sin(v_Y) , -std::sin(v_X)*std::cos(v_Y) };

  qhipster::TinyMatrix<Type, 2, 2, 32> U_noise;
  U_noise(0, 0) = A*B;
  U_noise(0, 1) = -std::conj(A)*std::conj(C);
  U_noise(1, 0) = A*C;
  U_noise(1, 1) =  std::conj(A)*std::conj(B);

  // Apply the noise gate
  QubitRegister<Type>::Apply1QubitGate(qubit,U_noise);
}


/// @brief Noise gate corresponding to single-qubit rotation with appropriate (stochastic) angle.
///
/// ** OLD OLD OLD OLD OLD OLD **
///
/// Kept for historical reasons, it shouldy be deletead.
///
/// To obtain a single rotation around an arbitrary axis we use the relations:
///      | a b c |               | h-f |
///  R = | d e f |    -->    u = | c-g |    ---> abs(u) = 2 sin( 'angle' )
///      | g h i |               | d-b |    ---> u/abs(u) = rotation axis  
template <class Type>
void NoisyQureg<Type>::NoiseGate_OLD(unsigned const qubit )
{
  BaseType t = time_from_last_gate[qubit] ;
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

  std::vector<BaseType> R = {  std::cos(v_Y)*std::cos(v_Z) ,
                              -std::cos(v_Y)*std::sin(v_Z) ,  
                               std::sin(v_Y) ,
                           std::sin(v_X)*std::sin(v_Y)*std::cos(v_Z) + std::cos(v_X)*std::sin(v_Z) ,
                          -std::sin(v_X)*std::sin(v_Y)*std::sin(v_Z) + std::cos(v_X)*std::cos(v_Z) ,
                          -std::sin(v_X) * std::cos(v_Y),
                              -std::cos(v_X)*std::sin(v_Y)*std::cos(v_Z) + std::sin(v_X)*std::sin(v_Z),
                               std::cos(v_X)*std::sin(v_Y)*std::sin(v_Z) + std::sin(v_X)*std::cos(v_Z),
                               std::cos(v_X)*std::cos(v_Y) };

  std::vector<BaseType> u = { R[3*2+1] - R[3*1+2] ,
                              R[3*0+2] - R[3*2+0] ,
                              R[3*1+0] - R[3*0+1] } ;
  BaseType norm_u = std::sqrt( std::norm(u[0]) + std::norm(u[1]) + std::norm(u[2]) );
  std::vector<BaseType> axis = { u[0]/norm_u , u[1]/norm_u , u[2]/norm_u };
  // Compute the angle of rotation
  BaseType trace_R =  R[0] + R[4] + R[8] ;
  BaseType angle = std::acos( (trace_R-1.)/2. );

  // Alternative expression for the angle, from Wikipedia:
  BaseType angle_alternative = std::asin( norm_u/2. ) ;
  // if     0<angle<Pi/2  then angle = angle_alternative
  // if  Pi/2<angle<Pi    then angle = Pi - angle_alternative

  // It may happen that "axis" still requires a "minus sign"
  if (false)
  {
      std::vector<double> R_u , R_minus_u ;
      BaseType si(std::sin(angle)) , co(std::cos(angle)) ;
      R_u = { co + axis[0]*axis[0]*(1.-co) ,
              axis[0]*axis[1]*(1.-co) - axis[2]*si ,
              axis[0]*axis[2]*(1.-co) + axis[1]*si ,
                axis[0]*axis[1]*(1.-co) + axis[2]*si ,
                co + axis[1]*axis[1]*(1.-co) ,
                axis[1]*axis[2]*(1.-co) - axis[0]*si ,
                  axis[0]*axis[2]*(1.-co) - axis[1]*si ,
                  axis[1]*axis[2]*(1.-co) + axis[0]*si ,
                  co + axis[2]*axis[2]*(1.-co) };
      R_minus_u = { co + axis[0]*axis[0]*(1.-co) ,
                    axis[0]*axis[1]*(1.-co) + axis[2]*si ,
                    axis[0]*axis[2]*(1.-co) - axis[1]*si ,
                      axis[0]*axis[1]*(1.-co) - axis[2]*si ,
                      co + axis[1]*axis[1]*(1.-co) ,
                      axis[1]*axis[2]*(1.-co) + axis[0]*si ,
                        axis[0]*axis[2]*(1.-co) + axis[1]*si ,
                        axis[1]*axis[2]*(1.-co) - axis[0]*si ,
                        co + axis[2]*axis[2]*(1.-co) };

//      print_vector(R,"rotation matrix:\n");
//      print_vector(R_u,"reconstructed with u:\n");
//      print_vector(R_minus_u,"reconstructed with -u:\n");

      // u was defined up to a sign. To resolve the ambiguity:
      BaseType R_element , Ru_element;
      if (axis[0] != 0.)
      {
          R_element = std::cos(v_X) * std::sin(v_Y) * std::sin(v_Z) + std::sin(v_X) * std::cos(v_Z) ;
          Ru_element = axis[2]*axis[1]*(1-std::cos(angle)) + axis[0]*std::sin(angle) ;
          if ( std::abs(R_element - Ru_element)> 1e-7 ) std::cout << " wrong sign from axis[0]\n";
//          else std::cout << " correct sign from axis[0]\n";
      }
      if (axis[1] != 0.)
      {
          R_element = std::sin(v_Y) ;
          Ru_element = axis[2]*axis[0]*(1-std::cos(angle)) + axis[1]*std::sin(angle) ;
          if ( std::abs(R_element - Ru_element)> 1e-7 ) std::cout << " wrong sign from axis[1]\n";
//          else std::cout << " correct sign from axis[1]\n";
      }
      if (axis[2] != 0.)
      {
          R_element = std::sin(v_X) * std::sin(v_Y) * std::cos(v_Z) + std::cos(v_X) * std::sin(v_Z) ;
          Ru_element = axis[0]*axis[1]*(1-std::cos(angle)) + axis[2]*std::sin(angle) ;
          if ( std::abs(R_element - Ru_element)> 1e-7 ) std::cout << " wrong sign from axis[2]\n";
//          else std::cout << " correct sign from axis[2]\n";
      }
  }

  // Compute the angle of rotation
  trace_R =  std::cos(v_Y) * std::cos(v_Z)
             - std::sin(v_X) * std::sin(v_Y) * std::sin(v_Z)
             + std::cos(v_X) * std::cos(v_Z)
             + std::cos(v_X) * std::cos(v_Y) ;
  angle = std::acos( (trace_R-1.)/2. );
  if (false) std::cout << " angle of rotation = " << angle << "\n";
  // Alternative expression for the angle, from Wikipedia.
  if (false) std::cout << " angle of rotation (alternative formula) = " << std::asin( norm_u/2. ) << "\n";


  // Costruct the 1/2-spin matrix corresponding to the axis above
  //    sigma_axis
  // and use it to implement the single-qubit rotation :
  //    rot = exp( i * angle * sigma_axis )
  qhipster::TinyMatrix<Type, 2, 2, 32> rot;
  BaseType s(std::sin(angle/2.)) , c(std::cos(angle/2.)) ;
  rot(0, 0) = Type(  c         ,  s*axis[2] );
  rot(0, 1) = Type(  s*axis[1] ,  s*axis[0] );
  rot(1, 0) = Type( -s*axis[1] ,  s*axis[0] );
  rot(1, 1) = Type(  c         , -s*axis[2] );

  // Apply the noise gate
  QubitRegister<Type>::Apply1QubitGate(qubit,rot);
}


// ----------------------------------------------------------
// -------- list of modified gates --------------------------
// ----------------------------------------------------------


template <class Type>
void NoisyQureg<Type>::Apply1QubitGate(unsigned const q,
                                       qhipster::TinyMatrix<Type, 2, 2, 32> V)
{
  AddNoiseOneQubitGate(q); 
  QubitRegister<Type>::Apply1QubitGate(q,V);
}

template <class Type>
void NoisyQureg<Type>::ApplyHadamard(unsigned const q)
{
  AddNoiseOneQubitGate(q); 
  QubitRegister<Type>::ApplyHadamard(q);
}

template <class Type>
void NoisyQureg<Type>::ApplyRotationX(unsigned const q, BaseType theta)
{
  AddNoiseOneQubitGate(q); 
  QubitRegister<Type>::ApplyRotationX(q,theta);
}

template <class Type>
void NoisyQureg<Type>::ApplyRotationY(unsigned const q, BaseType theta)
{
  AddNoiseOneQubitGate(q); 
  QubitRegister<Type>::ApplyRotationY(q,theta);
}

template <class Type>
void NoisyQureg<Type>::ApplyRotationZ(unsigned const q, BaseType theta)
{
  AddNoiseOneQubitGate(q); 
  QubitRegister<Type>::ApplyRotationZ(q,theta);
}

template <class Type>
void NoisyQureg<Type>::ApplyCPauliX(unsigned const q1, unsigned const q2)
{
  AddNoiseTwoQubitGate(q1,q2); 
  QubitRegister<Type>::ApplyCPauliX(q1,q2);
}

template <class Type>
void NoisyQureg<Type>::ApplyControlled1QubitGate(unsigned const q1, unsigned const q2,
                                                 qhipster::TinyMatrix<Type, 2, 2, 32> V)
{
  AddNoiseTwoQubitGate(q1,q2); 
  QubitRegister<Type>::ApplyControlled1QubitGate(q1,q2,V);
}

/// @}

#endif	// header guard NOISY_QUREG_H
