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
#include "openmp_affinity_corei7.hpp"
#include <string>

using namespace qhipster::openmp;


AffinityCoreI7::AffinityCoreI7(){ }
AffinityCoreI7::~AffinityCoreI7(){ }

void AffinityCoreI7::set_thread_affinity(int thread_count) { 

    // Core I7 processors generally have 2 HW threads per processor core. So,
    // we'll set that here.
    omp_set_dynamic(0);       // Explicitly disable dynamic teams.
    omp_set_num_threads(thread_count);

#if defined(__ICC) || defined(__INTEL_COMPILER)
    affinity_str = "KMP_AFFINITY=compact,1,0,granularity=fine";
    kmp_set_defaults(affinity_str.c_str());
#else
#pragma message("No kmp_set_defaults() function currently implemented for GNU g++")
#endif

}


unsigned AffinityCoreI7::get_num_threads( ) { 
    int th = 3;
#pragma omp parallel
#pragma omp master
    {
        th = omp_get_num_threads();
    }
    return th;
}


std::string AffinityCoreI7::get_affinity_string() {
    return affinity_str;
}

#ifdef _OPENMP
qhipster::openmp::AffinityCoreI7 glb_affinity;
#endif
