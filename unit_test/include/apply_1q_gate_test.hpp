#ifndef APPLY_1Q_GATE_TEST_HPP
#define APPLY_1Q_GATE_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class Apply1QGateTest : public ::testing::Test
{
 protected:

  Apply1QGateTest()
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
        GTEST_SKIP() << "INFO: small state distributed among too many ranks.";

    std::cout << "state_rank_id = " << iqs::mpi::Environment::GetStateRank() << "\n";//FIXME delete
    iqs::mpi::StateBarrier();
  }

  const std::size_t num_qubits_ = 4;
  double accepted_error_ = 1e-15;
  // The gates involve the last qubit. If the state is |000j>, one has:
  int qubit_ = num_qubits_-1;
  std::size_t j0_ = 0, j1_ = 8; // j0_ is index when j=0, j1_ is index when j=1.
  double sqrt2_ = std::sqrt(2.);
};

//////////////////////////////////////////////////////////////////////////////
// For all 1-qubit gates we test the expected behavior:
// - on the last qubit
// - for initial states in orthogonal basis (X and Z eigensattes)
//////////////////////////////////////////////////////////////////////////////

// Hadamard gate
TEST_F(Apply1QGateTest, Hadamard)
{
  ComplexDP amplitude_0, amplitude_1;
  // Initial state |000>|0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",j0_);
  psi.ApplyHadamard(qubit_);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { 1./sqrt2_, 0 };
  amplitude_1 = { 1./sqrt2_, 0 };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|1>
  psi.Initialize("base",j1_);
  psi.ApplyHadamard(qubit_);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { 1./sqrt2_, 0 };
  amplitude_1 = {-1./sqrt2_, 0 };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // At this point we have verified that the Hadamard works properly on |0>,|1>.
  // It is enough to verify that H.H = id
  // Initial state |000>|+>
  psi.Initialize("base",j0_);
  psi.ApplyHadamard(qubit_);
  psi.ApplyHadamard(qubit_);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { 1. , 0. };
  amplitude_1 = { 0. , 0. };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|->
  psi.Initialize("base",j1_);
  psi.ApplyHadamard(qubit_);
  psi.ApplyHadamard(qubit_);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { 0. , 0. };
  amplitude_1 = { 1. , 0. };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

// exp( -i X theta/2 )
TEST_F(Apply1QGateTest, RotationX)
{
  double angle = 0.83;
  ComplexDP amplitude_0, amplitude_1;
  // Initial state |000>|0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",j0_);
  psi.ApplyRotationX(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { std::cos(angle/2.), 0                 };
  amplitude_1 = { 0                 ,-std::sin(angle/2.)};
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|1>
  psi.Initialize("base",j1_);
  psi.ApplyRotationX(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { 0                 ,-std::sin(angle/2.)};
  amplitude_1 = { std::cos(angle/2.), 0                 };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|+>
  psi.Initialize("base",j0_);
  psi.ApplyHadamard(qubit_);
  psi.ApplyRotationX(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { std::cos(angle/2.)/sqrt2_, -std::sin(angle/2.)/sqrt2_};
  amplitude_1 = { std::cos(angle/2.)/sqrt2_, -std::sin(angle/2.)/sqrt2_};
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|->
  psi.Initialize("base",j1_);
  psi.ApplyHadamard(qubit_);
  psi.ApplyRotationX(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = {  std::cos(angle/2.)/sqrt2_,  std::sin(angle/2.)/sqrt2_};
  amplitude_1 = { -std::cos(angle/2.)/sqrt2_, -std::sin(angle/2.)/sqrt2_};
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

// exp( -i Y theta/2 )
TEST_F(Apply1QGateTest, RotationY)
{
  double angle = 0.75;
  ComplexDP amplitude_0, amplitude_1;
  // Initial state |000>|0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",j0_);
  psi.ApplyRotationY(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = {  std::cos(angle/2.), 0.};
  amplitude_1 = {  std::sin(angle/2.), 0.};
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|1>
  psi.Initialize("base",j1_);
  psi.ApplyRotationY(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = { -std::sin(angle/2.), 0.};
  amplitude_1 = {  std::cos(angle/2.), 0.};
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

// exp( -i Z theta/2 )
TEST_F(Apply1QGateTest, RotationZ)
{
  double angle = 0.35;
  ComplexDP amplitude_0, amplitude_1;
  // Initial state |000>|0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",j0_);
  psi.ApplyRotationZ(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = {  std::cos(angle/2.), -std::sin(angle/2.) };
  amplitude_1 = {  0.                , 0.                    };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
  // Initial state |000>|1>
  psi.Initialize("base",j1_);
  psi.ApplyRotationZ(qubit_,angle);
  ASSERT_DOUBLE_EQ(psi.ComputeNorm(), 1.);
  amplitude_0 = {  0.                , 0.                    };
  amplitude_1 = {  std::cos(angle/2.), std::sin(angle/2.) };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j0_), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(j1_), amplitude_1, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

// Rotation in the XY-plane: exp( -i P theta/2 )
// with P = cos(phi) X + sin(phi) Y
// It can be decomposed in (time from right to left as for matrix multiplications):
//   R_XY(phi, theta) = R_Z(phi).R_X(theta).R_Z(-phi)
TEST_F(Apply1QGateTest, RotationXY)
{
  // Initial state |....> at random
  iqs::RandomNumberGenerator<double> rng;
  std::size_t seed = 777;
  rng.SetSeedStreamPtrs(seed);
  iqs::QubitRegister<ComplexDP> psi (num_qubits_);
  psi.SetRngPtr(&rng);
  psi.Initialize("rand", 0);
  iqs::QubitRegister<ComplexDP> psi_0(psi);
  // Apply rotation in the XY plane and then cancel it with two rotations.
  double phi = 0.87;
  double theta = 0.45;
  unsigned qubit = 0;
  psi.ApplyRotationXY(qubit, phi, theta);
  psi.ApplyRotationZ(qubit, -phi);
  psi.ApplyRotationX(qubit, -theta);
  psi.ApplyRotationZ(qubit,  phi);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).real() - 1) < accepted_error_);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).imag() - 0) < accepted_error_);
  ASSERT_TRUE( psi.ComputeOverlap(psi_0).real() > 1.-accepted_error_);
  // Apply rotation in the XY plane and then cancel it with two rotations.
  phi = 1.1;
  theta = 0.13;
  qubit = 3;
  psi.ApplyRotationXY(qubit, phi, theta);
  psi.ApplyRotationZ(qubit, -phi);
  psi.ApplyRotationX(qubit, -theta);
  psi.ApplyRotationZ(qubit,  phi);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).real() - 1) < accepted_error_);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).imag() - 0) < accepted_error_);
  // Special case phi=0 corresponding to a X rotation.
  phi = 0.;
  theta = 0.28;
  qubit = 1;
  psi.ApplyRotationXY(qubit, phi, theta);
  psi.ApplyRotationX(qubit, -theta);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).real() - 1) < accepted_error_);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).imag() - 0) < accepted_error_);
  // Special case phi=pi/2 corresponding to a Y rotation.
  phi = M_PI/2;
  qubit = 2;
  psi.ApplyRotationXY(qubit, phi, theta);
  psi.ApplyRotationY(qubit, -theta);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).real() - 1) < accepted_error_);
  ASSERT_TRUE( std::abs(psi.ComputeOverlap(psi_0).imag() - 0) < accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(Apply1QGateTest, CustomGate)
{

  TM2x2<ComplexDP> G;
  G(0, 0) = {0.592056606032915, 0.459533060553574}; 
  G(0, 1) = {-0.314948020757856, -0.582328159830658};
  G(1, 0) = {0.658235557641767, 0.070882241549507}; 
  G(1, 1) = {0.649564427121402, 0.373855203932477};
  // |psi> = |1010> = |"1+4">
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",1+4);
  for(int qubit = 0; qubit < num_qubits_; qubit++)
  {
      psi.Apply1QubitGate(qubit, G);
  }
  ASSERT_DOUBLE_EQ( psi.GetProbability(0), psi.GetProbability(2) );
  ASSERT_DOUBLE_EQ( psi.GetProbability(1), psi.GetProbability(3) );
  // The accepted error is sometimes exceeded when gcc is used.
  ASSERT_NEAR( psi.GetProbability(0), 1.-psi.GetProbability(1), accepted_error_*10 );
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard APPLY_1Q_GATE_TEST_HPP
