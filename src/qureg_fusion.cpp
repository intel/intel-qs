#include "../include/qureg.hpp"

/// \addtogroup qureg
/// @{

/// @file qureg_fusion.cpp
/// @brief Define the @c QubitRegister methods related to an optimization called 'fusion'.

namespace iqs {

/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::TurnOnFusion(unsigned log2llc)
{
  unsigned myrank=0, nprocs=1, log2_nprocs=0;
#ifdef INTELQS_HAS_MPI
  myrank = iqs::mpi::Environment::GetStateRank();
  nprocs = iqs::mpi::Environment::GetStateSize();
  log2_nprocs = iqs::ilog2(nprocs);
#endif
  unsigned M = num_qubits - log2_nprocs;

  if (log2llc >= M) {
    if (!myrank) printf("Fusion is not enabled: num_qubits (%lu) is too small\n", num_qubits);
    fusion = false;
  } else {
    if (!myrank) printf("Fusion is enabled: log2llc = %u num_qubits = %lu\n", log2llc, num_qubits);
    this->log2llc = log2llc;
    fusion = true;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::TurnOffFusion()
{
  if (fwindow.size())
  {
      ApplyFusedGates();
  }
  fusion = false;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
bool QubitRegister<Type>::IsFusionEnabled()
{
  return fusion;
}


/////////////////////////////////////////////////////////////////////////////////////////
template <class Type>
void QubitRegister<Type>::ApplyFusedGates()
{
  #if 0
  std::size_t myrank = iqs::mpi::Environment::GetStateRank();
  if ( myrank==0 && fwindow.size() > 1) {
    printf("fused: ");
    for (auto &f : fwindow) {
      std::string &type = std::get<0>(f);
      unsigned qubit1 = std::get<2>(f);
      unsigned qubit2 = std::get<3>(f);
      if (type == "sqg")
        printf("[%s %u] ",  type.c_str(), qubit1);
      else
        printf("[%s %u %u] ",  type.c_str(), qubit1, qubit2);
    }
    printf("\n");
  }
  #endif

  size_t blocksize = (fwindow.size() == 1) ? LocalSize() : (1L << UL(log2llc));
  for (size_t l = 0; l < LocalSize(); l += blocksize) {
    for (auto &f : fwindow) {
      // check if work has been done
      std::string &type = std::get<0>(f);
      TM2x2<Type> &m = std::get<1>(f);
      std::size_t qubit1 = std::get<2>(f);
      std::size_t qubit2 = std::get<3>(f);
      if (type == "sqg") {
        Apply1QubitGate_helper(qubit1, m, l, l + blocksize);
      } else if (type == "cqg") {
        ApplyControlled1QubitGate_helper(qubit1, qubit2, m, l, l + blocksize);
      } else {
        assert(0); // not yet implemented
      }

    }
  }
  fwindow.resize(0);
}

template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

} // end namespace iqs

/// @}
