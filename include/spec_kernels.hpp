#ifndef SPEC_KERNELS_HPP
#define SPEC_KERNELS_HPP

#include "qureg.hpp"
#include "gate_spec.hpp"

using qhipster::GateSpec1Q;
using qhipster::GateSpec2Q;

template< typename Type >
__attribute__((noinline))
void Loop_SN(std::size_t gstart, std::size_t gend,
             Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1,
	         GateSpec1Q spec, Timer *timer, double angle);

template< typename Type >
__attribute__((noinline))
void Loop_DN(std::size_t gstart, std::size_t gend, std::size_t pos,
             Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1,
	         GateSpec1Q spec, Timer *timer, double angle=0);

template <class Type>
__attribute__((noinline))
void Loop_TN(Type *state,
             std::size_t c11, std::size_t c12, std::size_t c13,
             std::size_t c21, std::size_t c22, std::size_t c23,
             std::size_t c31, std::size_t c32, 
             std::size_t ind_shift, GateSpec2Q spec, Timer *timer, 
             double angle=0);

#endif
