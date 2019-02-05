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

#pragma once

#include "openmp_affinity_v1.hpp"
#include <omp.h>

namespace qhipster {
namespace openmp {

class AffinityCoreI7 : public qhipster::openmp::IOmpAffinityV1 {
    private:
        std::string affinity_str;

    public:
        AffinityCoreI7( );
        ~AffinityCoreI7( );

    void set_thread_affinity(int);
    unsigned get_num_threads();
    std::string get_affinity_string();
};

}
}

