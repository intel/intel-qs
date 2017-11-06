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
#include <mpi.h>
#include <exception>
#include <string>

namespace qhipster {

class MpiWrapperException : public std::exception {
    public:
        MpiWrapperException(int);
        virtual ~MpiWrapperException() throw();
        virtual const char* what() const throw() { return this->_ec_text.c_str(); }

    private:
        int _ec_value;
        std::string _ec_text;
};

}


#define QH_MPI_STATUS_CHECK(_api)\
    {\
        int _scode = _api;\
        if (MPI_SUCCESS != _scode) throw qhipster::MpiWrapperException(_scode);\
    }
