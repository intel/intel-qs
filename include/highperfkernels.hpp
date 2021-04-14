#ifndef HIGH_PERF_KERNELS_HPP
#define HIGH_PERF_KERNELS_HPP

#include "qureg.hpp"

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////

template< typename Type >
__attribute__((noinline))
void Loop_SN(std::size_t start, std::size_t end, Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1, TM2x2<Type> const&m, 
             bool specialize, Timer *timer);


template< typename Type >
__attribute__((noinline))
void Loop_DN(std::size_t gstart, std::size_t gend, std::size_t pos,
             Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1,
             TM2x2<Type> const&m,
             bool specialize, Timer *timer);

template< typename Type >
__attribute__((noinline))
void Loop_TN(Type *state, std::size_t c11, std::size_t c12,
             std::size_t c13, std::size_t c21, std::size_t c22, std::size_t c23, 
             std::size_t c31, std::size_t c32, std::size_t ind_shift, 
             TM2x2<Type> const&m, bool specialize, Timer *timer);

template< typename Type >
void ScaleState(std::size_t start, std::size_t end, Type *state,
                const Type &s, Timer *timer);

// for debugging purposes
static __attribute__((noinline)) void LOOP_START() {printf("\n");}
static __attribute__((noinline)) void LOOP_END() {printf("\n");}

/////////////////////////////////////////////////////////////////////////////////////////

} // close namespace iqs

#endif	// header guard HIGH_PERF_KERNELS_HPP
