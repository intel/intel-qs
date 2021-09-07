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

#include "../include/qureg.hpp"

#ifdef USE_MKL
#include <mkl.h>
#ifdef INTELQS_HAS_MPI
#include <mkl_cdft.h>
#endif
#endif


/////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void qft(iqs::QubitRegister<Type> &psi)
{
  int my_rank = iqs::mpi::Environment::GetStateRank();
  int n = iqs::ilog2(psi.size());

  // main computation
  //
  //  6  -H--*----*------*--------...--*----------------x------  6
  //  5  ----R-H--|-*----|-*------...--|-*--------------|-x----  5
  //  4  ---------R-R-H--|-|-*----...--|-|-*------------|-|-x--  4
  //  3  ----------------R-R-R-H--...--|-|-|-*----------|-|-|--  3
  //  2  -------------------------...--|-|-|-|-*--------|-|-x--  2
  //  1  -------------------------...--|-|-|-|-|-*------|-x----  1
  //  0  -------------------------...--R-R-R-R-R-R-H----x------  0
  //
  // See: https://en.wikipedia.org/wiki/Quantum_Fourier_transform#Circuit_implementation
  //
  // Note that:
  // 1. the qubit indices are inverted since the quantum-information convention has
  //    qubit 0 corresponding to the highest-significance bit, while IQS convention has
  //    qubit 0 corresponding to the lowest-significance bit (as standard in computer science).
  // 2. in the phase-shift gate, the role of target and control is interchangeable.
  for (int i = n - 1; i >= 0; i--)
  {
    for (int j = n - 1; j > i; j--)
    {
      int k = j - i;
      iqs::TinyMatrix<Type, 2, 2, 32> phaseshift;
      phaseshift(0, 0) = {1, 0};
      phaseshift(0, 1) = {0, 0};
      phaseshift(1, 0) = {0, 0};
      phaseshift(1, 1) = Type(std::cos(M_PI / (double)(UL(1) << UL(k))),
                              std::sin(M_PI / (double)(UL(1) << UL(k))) );
//      if (!my_rank) std::cout << "CP(" << j << "," << i << ")\n";
      psi.ApplyControlled1QubitGate(j, i, phaseshift);
    }
//      if (!my_rank) std::cout << "H(" << i << ")\n";
    psi.ApplyHadamard(i);
  }

  // perform swapping
  for (int i = 0; i < (n / 2); i++) {
//      if (!my_rank) std::cout << "SWAP(" << i << "," << n-1-i << ")\n";
    psi.ApplySwap(i, n - 1 - i);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
static void cfft(iqs::QubitRegister<Type> &x)
{
  
  int nprocs = iqs::mpi::Environment::GetStateSize();
#ifdef USE_MKL
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = iqs::mpi::Environment::GetStateComm();
  DFTI_DESCRIPTOR_DM_HANDLE desc;
  MKL_LONG v;

  // Create descriptor for 1D FFT
  MKL_LONG status =
      (sizeof(iqs::BaseType<Type>) == 8) ? 
      DftiCreateDescriptorDM(comm, &desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, x.GlobalSize()) :
      DftiCreateDescriptorDM(comm, &desc, DFTI_SINGLE, DFTI_COMPLEX, 1, x.GlobalSize());
  assert(status==0);

  // iqs::mpi::barrier(); exit(0);
  DftiGetValueDM(desc, CDFT_LOCAL_SIZE, &v);
  assert(status==0);
  std::vector <Type> work(v);
  status = DftiSetValueDM(desc, CDFT_WORKSPACE, work.data());
  assert(status==0);
  status = DftiSetValueDM(desc, DFTI_BACKWARD_SCALE,
                          1.0 / std::sqrt((iqs::BaseType<Type>)x.GlobalSize()));
  assert(status==0);
  status = DftiCommitDescriptorDM(desc);
  assert(status==0);
  status = DftiComputeBackwardDM(desc, &(x[0]));
  assert(status==0);
  status = DftiFreeDescriptorDM(&desc);
  assert(status==0);
#else
  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status =
      (sizeof(iqs::BaseType<Type>) == 8) ? 
      DftiCreateDescriptor(&descriptor, DFTI_DOUBLE, DFTI_COMPLEX, 1, x.GlobalSize()) :
      DftiCreateDescriptor(&descriptor, DFTI_SINGLE, DFTI_COMPLEX, 1, x.GlobalSize());
  assert(status==0);
  status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_INPLACE);  
  assert(status==0);
  status = DftiSetValue(descriptor, DFTI_BACKWARD_SCALE, 1.0 / std::sqrt((double)x.size()));
  assert(status==0);
  status = DftiCommitDescriptor(descriptor);           // Finalize the descriptor
  assert(status==0);
  status = DftiComputeBackward(descriptor, &(x[0]), NULL);  // Compute the Backward FFT
  assert(status==0);
  status = DftiFreeDescriptor(&descriptor);            // Free the descriptor
  assert(status==0);
#endif
#else
  assert(nprocs == 1);
  TODO(Replace with distributed FFTW);
#if 0
// Original code was:
  std::vector<Type, Alloc> y = x;
#else
// Alternative code is:
// std::vector<Type, openqu::AlignedAllocator<Type, 256> > y = x;
  assert(x.LocalSize() == x.GlobalSize() );
  std::vector<Type> y ;
  for (std::size_t i = 0; i< x.LocalSize(); ++i)
      y.push_back(x[i]);
#endif
  int N = y.size();
  for (int k = 0; k < N; k++)
  {
    y[k] = {0, 0};
    for (int j = 0; j < N; j++)
    {
      double arg = 2.0 * M_PI * double(j) * double(k) / double(N);
      Type e = Type(std::cos(arg), std::sin(arg));
      y[k] += x[j] * e;
    }
    y[k] /= std::sqrt(double(N));
  }
#if 0
// Original code was:
  x = y;
#else
// Alternative code is:
  for (std::size_t i = 0; i< x.LocalSize(); ++i)
      x[i] = y[i];
#endif
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
  iqs::mpi::Environment env(argc, argv);
  if (env.IsUsefulRank() == false) return 0;
  int myrank = env.GetStateRank();
  if (myrank==0)
      std::cout << " ----------------------------------------------------------- \n"
                << " ---- QFT test: if MPI is enabled, more info are provided -- \n"
                << " ----------------------------------------------------------- \n";

  int num_qubits = 3;
  if (argc != 2)
  {
    if (myrank==0)
        fprintf(stderr, "usage: %s <num_qubits>\n", argv[0]);
    exit(1);
  }
  else
  {
    num_qubits = atoi(argv[1]);
  }
#ifndef USE_MKL
  if (myrank==0)
      std::cout << "Without MKL, the implementation of the classical FFT is slow\n";
#endif

  // Single precision.
  {
    using Type = ComplexSP;

    if (myrank == 0) std::cout << "\nstate initialization (single precision)\n";
    iqs::RandomNumberGenerator<float> rng_sp;
    rng_sp.SetSeedStreamPtrs(777);
    iqs::QubitRegister<Type> psi1(num_qubits, "base", 0);
    psi1.SetRngPtr(&rng_sp);
    psi1.Initialize("rand",1);
    iqs::QubitRegister<Type> psi2(psi1);
  
    if (myrank == 0) std::cout << "CFFT (single precision)\n";
    cfft<Type>(psi1);
    psi2.EnableStatistics();
    psi2.TurnOnSpecialize();

    if (myrank == 0) std::cout << "QFT  (single precision)\n";
    qft<Type>(psi2);
    psi2.GetStatistics();
    double e1 = psi2.MaxAbsDiff(psi1);
    double e2 = psi2.MaxL2NormDiff(psi1);
    if (myrank == 0)
        printf("SP::qufft error vs classical max(absdiff: %le l2normdiff: %le)\n", e1, e2);
  }

  // Double precision.
  {
    using Type = ComplexDP;

    if (myrank == 0) std::cout << "\nstate initialization (double precision)\n";
    iqs::RandomNumberGenerator<double> rng_dp;
    rng_dp.SetSeedStreamPtrs(777);
    iqs::QubitRegister<Type> psi1(num_qubits, "base", 0);
    psi1.SetRngPtr(&rng_dp);
    psi1.Initialize("rand",1);
    assert( std::abs(psi1.ComputeNorm()-1.)<1e-10);
    iqs::QubitRegister<Type> psi2(psi1);

    if (myrank == 0) std::cout << "CFFT (double precision)\n";
    cfft<Type>(psi1);
    psi2.EnableStatistics();
    psi2.TurnOffSpecialize();

    if (myrank == 0) std::cout << "QFT  (double precision)\n";
    qft<Type>(psi2);
    psi2.GetStatistics();
    double e1 = psi2.MaxAbsDiff(psi1);
    double e2 = psi2.MaxL2NormDiff(psi1);
    if (myrank == 0)
        printf("DP::qufft error vs classical max(absdiff: %le l2normdiff: %le)\n", e1, e2);
  }
             
  return 0;
}
