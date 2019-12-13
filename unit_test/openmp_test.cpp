//------------------------------------------------------------------------------
// Copyright 2019 Intel Corporation
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

/// @file openmp_test.cpp
///
/// Tests for the OpenMP settings.

#include <chrono>
#include <iostream>
#include <math.h>	// for constant M_PI
#include <stdexcept>
#include <stdio.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////

void PrintHelloWorld(int thread_id)
{
    printf("Hello World from thread [%d]\n", thread_id);
}

/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  int num_threads = 1;
#ifdef _OPENMP
  num_threads = omp_get_num_threads();

  //////////////////////////////////////////////////////////////////////////////
  // First of all, print to screen from each and every thread (in random order).
#pragma omp parallel
  {
      int id = omp_get_thread_num();
      PrintHelloWorld(id);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Second, we want to compute the delta_x of an integral.
  double integral_start = 0.0,
         integral_end   = 1.0;
  unsigned long num_intervals = 4000000000;

  // Compute delta_x of the integral.
  double delta_x = (integral_end - integral_start) / num_intervals;
  double pi = 0.0;
    
  // Start performance check.
  auto start = std::chrono::steady_clock::now();
  double sum = 0.0;

#pragma omp parallel
  {
      int id = omp_get_thread_num();
      printf("Thread [%d]\n", id);
      // Compute the Riemann sum for the approximate integral.
#pragma omp for reduction(+:sum)
      for(unsigned long i=0; i<num_intervals; i++)
      {
          double x = (i+0.5) * delta_x;
          sum = sum + (4.0 / (1.0+x*x));
      }
  }
  pi = sum * delta_x;

  // End performance check.
  auto end = std::chrono::steady_clock::now();
  auto diff = end-start;

  printf("Result: %.16f w/ error=%.16f Exec Time=(%lu) ms\n",\
         pi,\
         M_PI-pi,\
         std::chrono::duration_cast<std::chrono::milliseconds>(diff).count());
#else
  printf("This program was compiled for execution without OpenMP.\m");
#endif

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
