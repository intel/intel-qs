#pragma once

/**
 * \addtogroup interface
 * @{
 */

/**
 * @file interface_api_qasm.h
 *
 * This file contains the interface to the QASM command handler function.  The QASM
 * operation is broken up into a command token and an args token.  The command is
 * mapped to a hash table of functions designed to handle the QASM op and then the
 * function is invoked and a result returned.
 *
 */

/**
 * Handle the execution of a QASM instruction in qHiPSTER.
 * @param op The QASM operation we wish to execute.
 * @param args The operands, if any, for the operation.
 * @return 0 Success - Op Handled\n
 *        >0 Failure - Op not handled or error.
 */
unsigned long ExecuteHandler(std::string op, std::string args);

/** @}*/
