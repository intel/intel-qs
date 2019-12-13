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

/// @file omp_test1.cpp
///
/// This file tests creation/destruction of an OMP session.

#include <stdexcept>
#include <stdio.h>

#include "openmp_affinity_corei7.hpp"

int main(int argc, char **argv)
{
    double integral_start = 0.0,
           integral_end   = 1.0;
    unsigned long N = 100000; // number of steps.

    // Compute delta_x of the integral.
    double delta_x = (integral_end - integral_start) / (double)N;
    double sum = 0.0;
    double x;

    // Compute the Riemann sum for the approximate integral.
    for(unsigned long i=0;i<N;i++) {
        x = (i+0.5)*delta_x;
        sum += (4.0 / (1.0+x*x));
    }

    sum *= delta_x;

    printf("Result: %.10f\n", sum);

    return 1;
}
