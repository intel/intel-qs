#ifndef APPLY_QUANTUM_CHANNEL_TEST_HPP
#define APPLY_QUANTUM_CHANNEL_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class ApplyQuantumChannel : public ::testing::Test
{
 protected:

  ApplyQuantumChannel()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
        GTEST_SKIP();

    // All tests are skipped if the 4-qubit state is distributed in more than 2^3 ranks.
    // In fact the MPI version needs to allocate half-the-local-storage for communication.
    // If the local storage is a single amplitude, this cannot be further divided.
    if (iqs::mpi::Environment::GetStateSize() > 8)
        GTEST_SKIP();

    iqs::mpi::StateBarrier();
  }

  const std::size_t num_qubits_ = 4;
  double accepted_error_ = 1e-15;
  double sqrt2_ = std::sqrt(2.);
};

//////////////////////////////////////////////////////////////////////////////

// Hadamard gate
TEST_F(ApplyQuantumChannel, IdealHadamard)
{
  CM4x4<ComplexDP> chi ;
  chi(0,0) = chi(0,1) = chi(0,2) = chi(0,3) = ComplexDP(0,0);
  chi(1,0) = chi(2,2) = ComplexDP(0,0);
  chi(1,1) = chi(1,3) = ComplexDP(0.5, 0);
  chi(2,0) = chi(2,1) = chi(2,2) = chi(2,3) = ComplexDP(0,0);
  chi(3,0) = chi(3,2) = ComplexDP(0,0);
  chi(3,1) = chi(3,3) = ComplexDP(0.5, 0);
  chi.EigensystemOfIdealHadamardChannel();
  // Initial state |0100>
  std::size_t index = 4;
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", index);
  // The application of quantum channels requires associating a RNG to psi.
  std::size_t rng_seed = 7777;
  iqs::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  psi.SetRngPtr(&rnd_generator);
  // Apply Hadamard on qubit 2, twice.
  int qubit = 2;
  psi.ApplyHadamard(qubit);
  // There is no need of averages since there is a unique non-zero eigenstate of chi.
  psi.ApplyChannel(qubit, chi);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(index), ComplexDP(1,0), accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

// Depolarizing channel:
//
//    rho' = (1-p) rho + p/3 ( X.rho.X + Y.rho.Y + Z.rho.Z )
TEST_F(ApplyQuantumChannel, DepolarizingChannel)
{
#ifndef IQS_WITH_NOISE
  GTEST_SKIP() << "INFO: Library Eigen is not used for noiseless simulations";
#else
  double p = 0.01;
  CM4x4<ComplexDP> chi;
  for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
      {
          if (i!=j)
              chi(i, j) = ComplexDP(0,0);
          else if (i==0)
              chi(i, j) = ComplexDP(1-p,0);
          else
              chi(i, j) = ComplexDP(p/3,0);
      }
  chi.SolveEigenSystem();
  // Initial state |00+1>
  std::size_t index = 1;
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", index);
  psi.ApplyHadamard(1);
  // The application of quantum channels requires associating a RNG to psi.
  std::size_t rng_seed = 7777;
  iqs::RandomNumberGenerator<double> rnd_generator;
  rnd_generator.SetSeedStreamPtrs(rng_seed);
  psi.SetRngPtr(&rnd_generator);
  // Keeping the qubits idle, apply the depolarizing channel 20 times per qubit. 
  // Take average over 100 realizations.
  int num_time_steps = 20;
  int num_ensemble_states = 100;
  std::vector<double> overlap_squared (num_time_steps);
  for (int s=0; s<num_ensemble_states; ++s)
  {
      iqs::QubitRegister<ComplexDP> psi_s(psi);
      overlap_squared[0] += std::norm( psi_s.ComputeOverlap(psi) );
      for (int t=1; t<num_time_steps; ++t)
      {
          for (int qubit=0; qubit<num_qubits_; ++qubit)
              psi.ApplyChannel(qubit, chi);
          overlap_squared[t] += std::norm( psi_s.ComputeOverlap(psi) );
      }
  }
  // Print the decaying overlap squared.
  if (iqs::mpi::Environment::GetStateRank() == 0)
  {
      std::cout << "Decay of overlap while exposed to depolarizing channel (p=" << p << "):\n";
      for (int t=0; t<num_time_steps; t+=10)
      {
          std::cout << "t=" << t << ",\t|<psi(t)|psi(0)>|^2 = "
                    << overlap_squared[t]/double(num_ensemble_states) << "\n";
      }
  }
#endif
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard APPLY_QUANTUM_CHANNEL_TEST_HPP
