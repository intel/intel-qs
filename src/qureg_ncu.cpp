#include "qureg.hpp"
#include "ncu.hpp"
#include "GateCache.hpp"
#include "mat_ops.hpp"
#include <iostream>

template <class Type>
void QubitRegister<Type>::ApplyNCU(
                TM2x2<Type> gate,
                const std::vector<std::size_t>& ctrl_indices, 
                const std::vector<std::size_t>& aux_indices, 
                unsigned const target)
{
    ncu::NCU<Type> ncu;
    ncu.initialiseMaps(ctrl_indices.size());

    auto gcache = ncu.getGateCache();

    // initial support for Pauli gates only
    std::string gate_label = "";
    if( gate == gcache.construct_pauli_x() ){
        gate_label = "X";
    }
    else if( gate == gcache.construct_pauli_y() ){
        gate_label = "Y";
    }
    else if( gate == gcache.construct_pauli_z() ){
        gate_label = "Z";
    }
    else {
        std::cerr << "Gate not currently support for NCU: " << gate.tostr() << std::endl;
        std::abort();
    }
    ncu.applyNQubitControl(*this, ctrl_indices, aux_indices, target, gate_label, gate, 0);
}

template void QubitRegister<ComplexDP>::ApplyNCU(
                TM2x2<ComplexDP> gate,
                const std::vector<std::size_t>& ctrl_indices, 
                const std::vector<std::size_t>& aux_indices, 
                unsigned const target);
template void QubitRegister<ComplexSP>::ApplyNCU(
                TM2x2<ComplexSP> gate,
                const std::vector<std::size_t>& ctrl_indices, 
                const std::vector<std::size_t>& aux_indices, 
                unsigned const target);
