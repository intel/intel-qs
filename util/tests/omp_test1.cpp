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
 * @file omp_test1.cpp
 *
 * This file tests creation/destruction of an OMP session.
 *
 */
#include "openmp_affinity_corei7.hpp"
#include "openmp_affinity_noomp.hpp"
#include <stdexcept>
#include <stdio.h>


#ifndef _OPENMP
qhipster::openmp::AffinityNoOmp affinity;
#else
qhipster::openmp::AffinityCoreI7 affinity;
#endif

void print_hello_world(int _id) {
    printf("Hello World [%d]\n", _id);
}


int main(int argc, char **argv) {

    //omp_set_num_threads(8);
    affinity.set_thread_affinity(8);

#pragma omp parallel
    {
        int ID = omp_get_thread_num();

        print_hello_world(ID);
    }

    return 1;
}
