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
 * @file interface_api_memory.h
 *
 * This file provides handlers for malloc/free operations on sets of qubits.
 *
 */


/**
 * Parses the number of qubits to allocate from args and attempts to create a
 * quantum register with that number of qubits.
 * @param args A string representing the number of qubits to allocate (integer 0..inf)
 * @return 0 Success\n 
 *        >0 Error/Warning with message to STDERR.
 *
 * __side effect WARNING__!! Operates on a global wavefunction variable named (@ref psi1)
 */
unsigned long qumalloc(std::string args);


/**
 * Free a qHiPSTER quantum register (wavefunction).
 * @param args Unused...can set to NULL.
 * @return 0 Success\n
 *        >0 Error/Warning with message to STDERR.
 *
 * __side effect WARNING__!! Operates on a global wavefunction variable named (@ref psi1)
 */
unsigned long qufree(std::string args);

/** @}*/
