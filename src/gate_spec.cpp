#include "../include/gate_spec.hpp"
using iqs::GateSpec1Q;
using iqs::GateSpec2Q;

GateSpec1Q iqs::ConvertSpec2to1(GateSpec2Q spec)
{
  switch (spec)
  {
    case GateSpec2Q::CHadamard:
      return GateSpec1Q::Hadamard;
    case GateSpec2Q::CRotationX:
      return GateSpec1Q::RotationX;
    case GateSpec2Q::CRotationY:
      return GateSpec1Q::RotationY;
    case GateSpec2Q::CRotationZ:
      return GateSpec1Q::RotationZ;
    case GateSpec2Q::CPauliX:
      return GateSpec1Q::PauliX;
    case GateSpec2Q::CPauliY:
      return GateSpec1Q::PauliY;
    case GateSpec2Q::CPauliZ:
      return GateSpec1Q::PauliZ;
    default:
      return GateSpec1Q::None;
  }
}
