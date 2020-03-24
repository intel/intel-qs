#include <algorithm>
#include <iostream>

#include "timer.hpp"

double time_in_seconds(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

/////////////////////////////////////////////////////////////////////////////////////////

void qhipster::WhatCompileDefinitions()
{
std::cout << "Compiler flags:\n"

#ifdef INTELQS_HAS_MPI
          << "          INTELQS_HAS_MPI --> [YES]\n"
#else
          << "          INTELQS_HAS_MPI --> [ NO]\n"
#endif

#ifdef USE_MKL
          << "                  USE_MKL --> [YES]\n"
#else
          << "                  USE_MKL --> [ NO]\n"
#endif

#ifdef _OPENMP
          << "                  _OPENMP --> [YES]\n"
#else
          << "                  _OPENMP --> [ NO]\n"
#endif

#ifdef __INTEL_COMPILER
          << "         __INTEL_COMPILER --> [YES]\n"
#else
          << "         __INTEL_COMPILER --> [ NO]\n"
#endif

#ifdef __ICC
          << "                    __ICC --> [YES]\n"
#else
          << "                    __ICC --> [ NO]\n"
#endif

#ifdef USE_MM_MALLOC
          << "            USE_MM_MALLOC --> [YES]\n"
#else
          << "            USE_MM_MALLOC --> [ NO]\n"
#endif

#ifdef STANDALONE
          << "               STANDALONE --> [YES]\n"
#else
          << "               STANDALONE --> [ NO]\n"
#endif

#ifdef NDEBUG
          << "                   NDEBUG --> [YES]\n"
#else
          << "                   NDEBUG --> [ NO]\n"
#endif

#ifdef __ONLY_NORMALIZED_STATES
          << " __ONLY_NORMALIZED_STATES --> [YES]\n";
#else
          << " __ONLY_NORMALIZED_STATES --> [ NO]\n";
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////
