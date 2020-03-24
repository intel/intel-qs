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
