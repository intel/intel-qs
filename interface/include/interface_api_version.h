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
