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

#include "../qureg/qureg.hpp"

#if (defined(__ICC) || defined(__INTEL_COMPILER))
#include <mkl.h>
#if defined(INTELQS_HAS_MPI)
#include <mkl_cdft.h>
#endif
#endif


/////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void qft(QubitRegister<Type> &psi)
{
  int n = openqu::ilog2(psi.size());

  // main computation
  for (int i = n - 1; i >= 0; i--)
  {
    for (int j = n - 1; j > i; j--)
    {
      int k = j - i;
      openqu::TinyMatrix<Type, 2, 2, 32> phaseshift;
      phaseshift(0, 0) = {1, 0};
      phaseshift(0, 1) = {0, 0};
      phaseshift(1, 0) = {0, 0};
      phaseshift(1, 1) = Type(cos(M_PI / D(UL(1) << UL(k))), sin(M_PI / D(UL(1) << UL(k))));
      psi.ApplyControlled1QubitGate(j, i, phaseshift);
    }
    psi.ApplyHadamard(i);
  }

  // perform swapping
  for (int i = 0; i < (n / 2); i++) {
    psi.ApplySwap(i, n - 1 - i);
  }
}



/////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
static void cfft(QubitRegister<Type> &x)
{
  
  int nprocs = openqu::mpi::Environment::size();
#if (defined(__ICC) || defined(__INTEL_COMPILER))
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  DFTI_DESCRIPTOR_DM_HANDLE desc;
  MKL_LONG v;

  // Create descriptor for 1D FFT
  MKL_LONG status =
      (sizeof(BaseType<Type>) == 8) ? 
      DftiCreateDescriptorDM(comm, &desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, x.GlobalSize()) :
      DftiCreateDescriptorDM(comm, &desc, DFTI_SINGLE, DFTI_COMPLEX, 1, x.GlobalSize());

  // openqu::mpi::barrier(); exit(0);
  DftiGetValueDM(desc, CDFT_LOCAL_SIZE, &v);
  std::vector <Type> work(v);
  DftiSetValueDM(desc, CDFT_WORKSPACE, work.data());

  
  DftiSetValueDM(desc, DFTI_BACKWARD_SCALE, 1.0 / std::sqrt((BaseType<Type>)x.GlobalSize()));
  DftiCommitDescriptorDM(desc);
  DftiComputeBackwardDM(desc, &(x[0]));
  DftiFreeDescriptorDM(&desc);
#else
  DFTI_DESCRIPTOR_HANDLE descriptor;
  MKL_LONG status;
  status = DftiCreateDescriptor(&descriptor, DFTI_DOUBLE, DFTI_COMPLEX, 1, x.GlobalSize());  
  status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_INPLACE);  
  status = DftiSetValue(descriptor, DFTI_BACKWARD_SCALE, 1.0 / sqrt((double)x.size()));
  status = DftiCommitDescriptor(descriptor);           // Finalize the descriptor
  status = DftiComputeBackward(descriptor, &(x[0]), NULL);  // Compute the Backward FFT
  status = DftiFreeDescriptor(&descriptor);            // Free the descriptor
#endif
#else
  assert(nprocs == 1);
  std::vector<Type, Alloc> y = x;
  TODO(Replace with distributed FFTW)
  int N = y.size();
  for (int k = 0; k < N; k++)
  {
    y[k] = {0, 0};
    for (int j = 0; j < N; j++)
    {
      double arg = 2.0 * M_PI * D(j) * D(k) / D(N);
      Type e = Type(std::cos(arg), std::sin(arg));
      y[k] += x[j] * e;
    }
    y[k] /= std::sqrt(D(N));
  }
  x = y;
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
#ifdef INTELQS_HAS_MPI
  openqu::mpi::Environment env(argc, argv);
  if (env.is_usefull_rank() == false) return 0;
  unsigned myrank = env.rank();
#else
  unsigned myrank=0;
  std::cout << " ---------------------------------------------------------- \n"
            << " ---- the QFT test works only with MPI and a single rank -- \n"
            << " ---------------------------------------------------------- \n";
  return 0;
#endif

  int num_qubits = 3;
  if (argc != 2)
  {
    fprintf(stderr, "usage: %s <num_qubits>\n", argv[0]);
    exit(1);
  }
  else
  {
    num_qubits = atoi(argv[1]);
  }

  if (myrank == 0) std::cout << "QFFT\n";

  {
  using Type = ComplexSP;
  QubitRegister<Type> psi1(num_qubits, "rand", 0), psi2(psi1);
  cfft<Type>(psi1);
  psi2.EnableStatistics();
  psi2.TurnOnSpecialize();
  qft<Type>(psi2);
  psi2.GetStatistics();
  double e1 = psi2.maxabsdiff(psi1);
  double e2 = psi2.maxl2normdiff(psi1);
  if (myrank == 0)
      printf("SP::qufft error vs classical max(absdiff: %le l2normdiff: %le)\n", e1, e2);
  }

  {
  using Type = ComplexDP;
  QubitRegister<Type> psi1(num_qubits, "rand", 0), psi2(psi1);
  cfft<Type>(psi1);
  psi2.EnableStatistics();
  psi2.TurnOffSpecialize();
  qft<Type>(psi2);
  psi2.GetStatistics();
  double e1 = psi2.maxabsdiff(psi1);
  double e2 = psi2.maxl2normdiff(psi1);
  if (myrank == 0)
      printf("DP::qufft error vs classical max(absdiff: %le l2normdiff: %le)\n", e1, e2);
  }
             
  return 0;
}
