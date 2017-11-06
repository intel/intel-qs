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
 * @file interface_api_qubitid.h
 *
 * This header file contains functions for conveniently handling the assignment
 * of qubit identifiers to qubit strings passed in as operands to QASM instructions.
 */


/**
 * query_qubit_id()
 * @param qubit_name The string name passed from the QASM instruction.
 * @return An integer representing the the qubit identifier that matches the
 *         qubit name string. If the qubit name is not found then an id is generated, 
 *         added to an internal tracking table, and the id is returned.
 */
int query_qubit_id(std::string qubit_name);

/** @}*/
