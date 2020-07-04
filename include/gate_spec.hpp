#ifndef GATE_SPEC_HPP
#define GATE_SPEC_HPP

namespace qhipster {

enum class GateSpec1Q { 
  Hadamard=0, 
  RotationX, RotationY, RotationZ, 
  PauliX, PauliY, PauliZ,
  T, 
  None
};

enum class GateSpec2Q {
  CHadamard=0,
  CRotationX, CRotationY, CRotationZ,
  CPauliX, CPauliY, CPauliZ,
  CPhase,
  None
};

}

#endif
