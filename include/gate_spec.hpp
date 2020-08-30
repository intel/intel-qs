/// @file gate_spec.hpp
/// @brief Declare structs used in specialization v2.

#ifndef GATE_SPEC_HPP
#define GATE_SPEC_HPP

namespace qhipster {

/// @brief Known 1-qubit gate types used internally with specialization v2
enum class GateSpec1Q { 
  Hadamard=0, 
  RotationX, RotationY, RotationZ, 
  PauliX, PauliY, PauliZ,
  T, 
  None
};

/// @brief Known 2-qubit gate types used internally with specialization v2
enum class GateSpec2Q {
  CHadamard=0,
  CRotationX, CRotationY, CRotationZ,
  CPauliX, CPauliY, CPauliZ,
  CPhase,
  None
};

GateSpec1Q ConvertSpec2to1(GateSpec2Q spec);

}

#endif
