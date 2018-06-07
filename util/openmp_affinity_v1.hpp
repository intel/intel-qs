//------------------------------------------------------------------------------
// Copyright (C) 2017 Intel Corporation 
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

#pragma once
#include <string>

namespace qhipster {
namespace openmp {

class IOmpAffinityV1 {

    public:
        IOmpAffinityV1( ) { }
        virtual ~IOmpAffinityV1( ) { }

    // Pure virtual functions.
    virtual void set_thread_affinity(int)  = 0;
    virtual unsigned get_num_threads() = 0;
    virtual std::string get_affinity_string() = 0;
};

}
}
