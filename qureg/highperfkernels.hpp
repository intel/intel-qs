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

#ifndef HIGH_PERF_KERNELS_H
#define HIGH_PERF_KERNELS_H

#include "qureg.hpp"

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

#endif	// header guard HIGH_PERF_KERNELS_H
