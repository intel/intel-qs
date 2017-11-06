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
#include "openmp_affinity_noomp.hpp"

using namespace qhipster::openmp;

AffinityNoOmp::AffinityNoOmp( ) { }
AffinityNoOmp::~AffinityNoOmp( ){ }

void AffinityNoOmp::set_thread_affinity(int thread_count) { 
    return; 
}

unsigned AffinityNoOmp::get_num_threads() { 
    return 1; 
}

std::string AffinityNoOmp::get_affinity_string() {
    return "";
}

#ifndef _OPENMP
qhipster::openmp::AffinityNoOmp glb_affinity;
#endif
