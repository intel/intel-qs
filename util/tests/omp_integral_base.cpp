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

/**
 * @file omp_integral_base.cpp
 *
 * This file tests the implementation of OMP wrapper by computing pi using a
 * Reimann sum approach.
 *
 */
#include "openmp_affinity_noomp.hpp"
#include "openmp_affinity_corei7.hpp"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdexcept>
#include <chrono>

#define GET_PI_DEVIATION(v) (M_PI-v)

#ifndef _OPENMP
qhipster::openmp::AffinityNoOmp affinity;
#else
qhipster::openmp::AffinityCoreI7 affinity;
#endif

/**
 */
int main(int argc, char **argv) {

    double integral_start = 0.0,
           integral_end   = 1.0;
    unsigned long N = 4000000000; // number of steps.

    // Compute delta_x of the integral.
    double delta_x = (integral_end - integral_start) / N;
    double pi = 0.0;
    
    const int global_num_threads = 12;

    // Parameter set for the specific architecture.
    //omp_set_num_threads(global_num_threads);
    affinity.set_thread_affinity(global_num_threads);

    // Start performance check.
    auto start = std::chrono::steady_clock::now();
    double sum = 0.0;

#pragma omp parallel
    {
        int threadID = omp_get_thread_num();

        printf("Thread %d\n", threadID);

        // Compute the Riemann sum for the approximate integral.
        #pragma omp for reduction(+:sum)
        for(unsigned long i=0;i<N;i++) {
            double x = (i+0.5) * delta_x;
            sum = sum + (4.0 / (1.0+x*x)); 
        }
    }

    pi = sum * delta_x;


    // End performance check.
    auto end = std::chrono::steady_clock::now();
    auto diff = end-start;

    printf("Result: %.16f w/ error=%.16f Exec Time=(%lu) ms\n",\
               pi,\
               GET_PI_DEVIATION(pi),\
               std::chrono::duration_cast<std::chrono::milliseconds>(diff).count());

    return 1;
}
