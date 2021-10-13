#ifndef EXPECTATION_VALUES_TEST_HPP
#define EXPECTATION_VALUES_TEST_HPP

#include "../../include/qureg.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class ExpectationValuesTest : public ::testing::Test
{
 protected:

  ExpectationValuesTest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();

    // All tests are skipped if the 4-qubit state is distributed in more than 8 ranks.
    if (iqs::mpi::Environment::GetStateSize() > 8)
      GTEST_SKIP() << "INFO: small state distributed among too many ranks.";
  }

  const std::size_t num_qubits_ = 6;
  double accepted_error_ = 1e-14;
  std::vector<unsigned> observables_;
  std::vector<unsigned> qubits_;
  double expectation_ = 0.;
  double coeff_ = 1.;
  double sqrt2_ = std::sqrt(2.);
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(ExpectationValuesTest, YZBaseChange)
{
  ComplexDP amplitude_0, amplitude_1;
  // In the implementation, the matrix G is used.
  // G is matrix from change of basis Y --> Z, such that Ginv.Z.G = Y 
  iqs::TinyMatrix<ComplexDP, 2, 2, 32> G;
  double f = 1. / std::sqrt(2.);
  G(0, 0) = G(1, 0) = ComplexDP(f , 0.);
  G(0, 1) = ComplexDP(0.,-f);
  G(1, 1) = ComplexDP(0., f);
  // G^dagger = G^-1
  iqs::TinyMatrix<ComplexDP, 2, 2, 32> Ginv;
  Ginv(0, 0) = Ginv(0, 1) = ComplexDP(f , 0.);
  Ginv(1, 0) = ComplexDP(0., f);
  Ginv(1, 1) = ComplexDP(0.,-f);

  // Test from Y to Z basis (using matrix G).
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",0);
  psi.ApplyRotationX(0,-M_PI/2.);
  // |psi> = |+y>|00000>
  psi.Apply1QubitGate(0,G);
  // |psi> = |0> |00000>
  amplitude_0 = {  1. ,  0. };
  amplitude_1 = {  0. ,  0. };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
  psi.Apply1QubitGate(0,Ginv);
  // |psi> = |+y> |00000>
  amplitude_0 = {  f  ,  0. };
  amplitude_1 = {  0. ,  f  };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
  //
  psi.Initialize("base",0);
  psi.ApplyRotationX(0, M_PI/2.);
  // |psi> = |-y>|00000>
  psi.Apply1QubitGate(0,G);
  // |psi> = |1> |00000>
  amplitude_0 = {  0. ,  0. };
  amplitude_1 = {  1. ,  0. };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
  psi.Apply1QubitGate(0,Ginv);
  // |psi> = |1> |00000>
  amplitude_0 = {  f  ,  0. };
  amplitude_1 = {  0. , -f  };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);

  // Test from Z to Y basis (using the inverse G^-1).
  psi.Initialize("base",0);
  psi.Apply1QubitGate(0,Ginv);
  // |psi> = |+y>|00000>
  amplitude_0 = {  f  ,  0. };
  amplitude_1 = {  0. ,  f  };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
  psi.Apply1QubitGate(0,G);
  // |psi> = |0> |00000>
  amplitude_0 = {  1. ,  0. };
  amplitude_1 = {  0. ,  0. };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
  //
  psi.Initialize("base",1);
  psi.Apply1QubitGate(0,Ginv);
  // |psi> = |-y>|00000>
  amplitude_0 = {  f  ,  0. };
  amplitude_1 = {  0. , -f  };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
  psi.Apply1QubitGate(0,G);
  // |psi> = |1> |00000>
  amplitude_0 = {  0. ,  0. };
  amplitude_1 = {  1. ,  0. };
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(0), amplitude_0, accepted_error_);
  ASSERT_COMPLEX_NEAR(psi.GetGlobalAmplitude(1), amplitude_1, accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ExpectationValuesTest, ExpectationOneQubit)
{
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",10);
  // |psi> = |010100> = |"2+8">
  psi.ApplyHadamard(2);
  psi.ApplyHadamard(3);
  psi.ApplyRotationX(4,-M_PI/2.);
  psi.ApplyRotationX(5, M_PI/2.);
  // |psi> = |01> |+-> |+y-y>
  ASSERT_NEAR(psi.ExpectationValueX(0, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueX(1, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueX(2, coeff_),  1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueX(3, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueX(4, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueX(5, coeff_),  0., accepted_error_);
  //
  ASSERT_NEAR(psi.ExpectationValueY(0, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueY(1, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueY(2, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueY(3, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueY(4, coeff_),  1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueY(5, coeff_), -1., accepted_error_);
  //
  ASSERT_NEAR(psi.ExpectationValueZ(0, coeff_),  1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZ(1, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZ(2, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZ(3, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZ(4, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZ(5, coeff_),  0., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ExpectationValuesTest, ExpectationTwoQubits)
{
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",10);
  // |psi> = |010100> = |"2+8">
  psi.ApplyHadamard(2);
  psi.ApplyHadamard(3);
  psi.ApplyRotationX(4,-M_PI/2.);
  psi.ApplyRotationX(5, M_PI/2.);
  // |psi> = |01> |+-> |+y-y>
  ASSERT_NEAR(psi.ExpectationValueXX(0, 1, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueXX(1, 2, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueXX(2, 3, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueXY(3, 4, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueXY(3, 5, coeff_),  1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueXZ(2, 1, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueXZ(2, 1, coeff_), -1., accepted_error_);
  //
  ASSERT_NEAR(psi.ExpectationValueYY(0, 1, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueYY(1, 2, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueYY(4, 5, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueYZ(4, 0, coeff_),  1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueYZ(4, 1, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueYX(5, 2, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueYX(5, 3, coeff_),  1., accepted_error_);
  //
  ASSERT_NEAR(psi.ExpectationValueZZ(0, 1, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZZ(1, 2, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZZ(2, 3, coeff_),  0., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZX(0, 2, coeff_),  1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZX(0, 3, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZY(1, 4, coeff_), -1., accepted_error_);
  ASSERT_NEAR(psi.ExpectationValueZY(1, 5, coeff_),  1., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(ExpectationValuesTest, GeneralMethod)
{
  iqs::QubitRegister<ComplexDP> psi (num_qubits_,"base",2+8);
  // |psi_0> = |010100> = |"2+8">
  ASSERT_EQ(psi.GetProbability(0),0);
  ASSERT_EQ(psi.GetProbability(1),1);
  psi.ApplyHadamard(0);
  psi.ApplyHadamard(1);
  // |psi> = |+-0100>
  ASSERT_DOUBLE_EQ(psi.GetProbability(0), 0.5);
  ASSERT_DOUBLE_EQ(psi.GetProbability(1), 0.5);
  
  // observable = X0 . X1 . Z2 . Z3
  qubits_      = {0,1,2,3};
  observables_ = {1,1,3,3};
  expectation_ = psi.ExpectationValue(qubits_, observables_, coeff_);
  ASSERT_NEAR(expectation_,  1., accepted_error_);
  // observable = X1 . X0 . Z2 . Z4
  qubits_      = {1,0,2,4};
  expectation_ = psi.ExpectationValue(qubits_, observables_, coeff_);
  ASSERT_NEAR(expectation_, -1., accepted_error_);
  // observable = X0 . Z2 . Z3 . Z5 . Z4
  qubits_      = {0,2,3,5,4};
  observables_ = {1,3,3,3,3};
  expectation_ = psi.ExpectationValue(qubits_, observables_, coeff_);
  ASSERT_NEAR(expectation_, -1., accepted_error_);
  // observable = X0 . id1 . Z2 . Z3
  qubits_      = {3,0,2};
  observables_ = {3,1,3};
  ASSERT_NEAR(psi.ExpectationValue(qubits_, observables_, coeff_), -1., accepted_error_);

  psi.Initialize("base",10);
  psi.ApplyHadamard(2);
  psi.ApplyHadamard(3);
  psi.ApplyRotationX(4,-M_PI/2.);
  psi.ApplyRotationX(5, M_PI/2.);
  // |psi> = |01> |+-> |+y-y>
  // observable = Z1 . X2 . Y5
  qubits_      = {1,2,5};
  observables_ = {3,1,2};
  ASSERT_NEAR(psi.ExpectationValue(qubits_, observables_, coeff_),  1., accepted_error_);
  // observable = Z0 . X2 . Y4 . Y5
  qubits_      = {0,2,4,5};
  observables_ = {3,1,2,2};
  ASSERT_NEAR(psi.ExpectationValue(qubits_, observables_, coeff_), -1., accepted_error_);
  // observable = Z0 . Z1 . X2 . X3 . Y4 . Y5
  qubits_      = {0,1,2,3,4,5};
  observables_ = {3,3,1,1,2,2};
  ASSERT_NEAR(psi.ExpectationValue(qubits_, observables_, coeff_), -1., accepted_error_);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard EXPECTATION_VALUES_TEST_HPP
