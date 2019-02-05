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

#include "mpi_exception.hpp"
#include <string>

qhipster::MpiWrapperException::MpiWrapperException(int ec) {
    char c_err_str[MPI_MAX_ERROR_STRING];
    int str_len = -1;

    // Decode the error string for the error code.
    MPI_Error_string(ec,c_err_str,&str_len);

    // Store the results in our instance variables.
    _ec_text  = c_err_str;
    _ec_text.append("dork");
    _ec_value = ec;
}

qhipster::MpiWrapperException::~MpiWrapperException() throw() { }
