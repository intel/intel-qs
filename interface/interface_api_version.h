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

/**
 * \addtogroup interface
 * @{
 */

/**
 * @file interface_api_version.h
 *
 * This header file contains functions for reporting the version of both the
 * QASM interface to qHiPSTER as well as the version of the qHiPSTER library
 * that was linked into the interface and in use.
 */

#define INTERFACE_VERSION_STRING "RC-1.1.0"

/**
 * Print out the QASM interface version string. The format will be 
 * major.minor.revision.
 * @param args Unused.
 * @return 0.
 */
unsigned long quiversion(std::string args);


/**
 * Print out the qHiPSTER library version string. The format will be 
 * major.minor.revision.
 * @param args Unused.
 * @return 0.
 */
unsigned long quversion(std::string args);

/** @} */
