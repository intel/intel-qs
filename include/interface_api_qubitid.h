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
