#ifndef QAOA_FEATURES_TEST_HPP
#define QAOA_FEATURES_TEST_HPP

#include "../../include/qureg.hpp"
#include "../../include/qaoa_features.hpp"

//////////////////////////////////////////////////////////////////////////////
// Test fixture class.

class QaoaFeaturestest : public ::testing::Test
{
 protected:

  QaoaFeaturestest()
  { }

  // just after the 'constructor'
  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
        GTEST_SKIP();

    // All tests are skipped if the 6-qubit state is distributed in more than 2^5 ranks.
    // In fact the MPI version needs to allocate half-the-local-storage for communication.
    // If the local storage is a single amplitude, this cannot be further divided.
    if (qhipster::mpi::Environment::GetStateSize() > 32)
        GTEST_SKIP();
  }

  const std::size_t num_qubits_ = 6;
  double accepted_error_ = 1e-15;
};

//////////////////////////////////////////////////////////////////////////////
// Functions developed to facilitate the simulation of QAOA circuits.
//////////////////////////////////////////////////////////////////////////////

TEST_F(QaoaFeaturestest, qaoa_maxcut)
{
  // Instance of the max-cut problem provided as adjacency matrix.
  // It is a ring of 6 vertices:
  //
  //   0--1--2
  //   |     |
  //   5--4--3
  //
  std::vector<int> adjacency = {0, 1, 0, 0, 0, 1,
                                1, 0, 1, 0, 0, 0,
                                0, 1, 0, 1, 0, 0,
                                0, 0, 1, 0, 1, 0,
                                0, 0, 0, 1, 0, 1,
                                1, 0, 0, 0, 1, 0};
  QubitRegister<ComplexDP> diag (num_qubits_,"base",0);
  int max_cut_value;
  max_cut_value = qaoa::InitializeVectorAsMaxCutCostFunction(diag,adjacency);

  // Among other properties, only two bipartition has cut=0.
  ComplexDP amplitude;
  amplitude = { 0, 0 };
  ASSERT_COMPLEX_NEAR(diag.GetGlobalAmplitude(0), amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR(diag.GetGlobalAmplitude(diag.GlobalSize()-1), amplitude,
                      accepted_error_);
  // No bipartition can cut a single edge.
  for (size_t j=0; j<diag.LocalSize(); ++j)
      ASSERT_GT( std::abs(diag[j].real()-1.), accepted_error_);

  // Perform QAOA simulation (p=1).
  QubitRegister<ComplexDP> psi  (num_qubits_,"++++",0);
  double gamma = 0.4;
  double beta  = 0.3;
  // Emulation of the layer based on the cost function: 
  qaoa::ImplementQaoaLayerBasedOnCostFunction(psi, diag, gamma);
  // Simulation of the layer based on the local transverse field: 
  for (int qubit=0; qubit<num_qubits_; ++qubit)
      psi.ApplyRotationX(qubit, beta);
  // Get average of cut value:
  double expectation = qaoa::GetExpectationValueFromCostFunction( psi, diag);
  
  // Get histogram of the cut values:
  std::vector<double> histo = qaoa::GetHistogramFromCostFunction(psi, diag, max_cut_value);
  ASSERT_EQ(histo.size(), max_cut_value+1);
  double average=0;
  for (int j=0; j<histo.size(); ++j)
      average += double(j)*histo[j];
  ASSERT_DOUBLE_EQ(expectation, average);

#if 0
  // If the permutation of diag and psi is different, then functions should fail.
  diag.EmulateSwap(0,1);
  EXPECT_DEATH( qaoa::ImplementQaoaLayerBasedOnCostFunction(psi, diag, gamma), "" );
  EXPECT_DEATH( qaoa::GetExpectationValueFromCostFunction( psi, diag), "" );
  EXPECT_DEATH( qaoa::GetHistogramFromCostFunction(psi, diag, max_cut_value), "" );
#endif
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(QaoaFeaturestest, qaoa_weighted_maxcut)
{
  // Instance of the max-cut problem provided as adjacency matrix.
  // It is a ring of 6 vertices:
  //
  //   0--1--2
  //   |     |
  //   5--4--3
  //
  // where the verical edges have weight 1.4
  std::vector<double> adjacency = {0  , 1  , 0  , 0  , 0  , 1.4,
                                   1  , 0  , 1  , 0  , 0  , 0  ,
                                   0  , 1  , 0  , 1.4, 0  , 0  ,
                                   0  , 0  , 1.4, 0  , 1  , 0  ,
                                   0  , 0  , 0  , 1  , 0  , 1  ,
                                   1.4, 0  , 0  , 0  , 1  , 0  };
  QubitRegister<ComplexDP> diag (num_qubits_,"base",0);
  double max_cut_value;
  max_cut_value = qaoa::InitializeVectorAsWeightedMaxCutCostFunction(diag,adjacency);

  // Among other properties, only two bipartition has cut=0.
  ComplexDP amplitude;
  amplitude = { 0, 0 };
  ASSERT_COMPLEX_NEAR(diag.GetGlobalAmplitude(0), amplitude, accepted_error_);
  ASSERT_COMPLEX_NEAR(diag.GetGlobalAmplitude(diag.GlobalSize()-1), amplitude,
                      accepted_error_);

  // Case in which only 2 is dis-aligned:
  // 001000 = 1*2^2
  amplitude = { 1+1.4, 0 };
  size_t index = 2*2;
  ASSERT_COMPLEX_NEAR(diag.GetGlobalAmplitude(index), amplitude, accepted_error_);

  // Case in which only 2 and 5 are dis-aligned:
  // 001001 = 1*2^2 + 1*2^5
  amplitude = { 1+1.4+1+1.4, 0 };
  index = 4+32;
  ASSERT_COMPLEX_NEAR(diag.GetGlobalAmplitude(index), amplitude, accepted_error_);

  // No bipartition can cut a single edge.
  for (size_t j=0; j<diag.LocalSize(); ++j)
      ASSERT_GT( std::abs(diag[j].real()-1.), accepted_error_);

  // Perform QAOA simulation (p=1).
  QubitRegister<ComplexDP> psi  (num_qubits_,"++++",0);
  double gamma = 0.4;
  double beta  = 0.3;
  // Emulation of the layer based on the cost function: 
  qaoa::ImplementQaoaLayerBasedOnCostFunction(psi, diag, gamma);
  // Simulation of the layer based on the local transverse field: 
  for (int qubit=0; qubit<num_qubits_; ++qubit)
      psi.ApplyRotationX(qubit, beta);
  // Get average of cut value:
  double expectation = qaoa::GetExpectationValueFromCostFunction(psi, diag);
  
  // Histogram for rounded cutvals and check if it matches expval to the tolerance.
  std::vector<double> histo = qaoa::GetHistogramFromCostFunctionWithWeightsRounded(psi, diag, max_cut_value);
  ASSERT_EQ(histo.size(), (int)(floor(max_cut_value))+1);
  double average=0;
  for (int j=0; j<histo.size(); ++j)
      average += double(j)*histo[j];
  // The expval will be within less than 1.0 of the actual since the cutvals are rounded down to nearest 1.0.
  ASSERT_TRUE( (abs(expectation - average) )<=1.0+1e-7);
    
  // Histogram for rounded cutvals and check if it matches expval to the tolerance.
  double bin_width = 0.1;
  std::vector<double> histo2 = qaoa::GetHistogramFromCostFunctionWithWeightsBinned(psi, diag, max_cut_value, bin_width);
  ASSERT_EQ(histo2.size(), (int)(ceil(max_cut_value / bin_width)) + 1);
  average = 0.0;
  for (int j=0; j<histo2.size(); ++j)
      average += double(j)*bin_width*histo2[j];
  // The expval will be within less than bin_width of the actual since the cutvals are rounded down to the bin_width.
  ASSERT_TRUE( (abs(expectation - average) )<=bin_width+1e-7);

#if 0
  // If the permutation of diag and psi is different, then functions should fail.
  diag.EmulateSwap(0,1);
  EXPECT_DEATH( qaoa::ImplementQaoaLayerBasedOnCostFunction(psi, diag, gamma), "" );
  EXPECT_DEATH( qaoa::GetExpectationValueFromCostFunction( psi, diag), "" );
  EXPECT_DEATH( qaoa::GetHistogramFromCostFunctionWithWeightsRounded(psi, diag, max_cut_value), "" );
  EXPECT_DEATH( qaoa::GetHistogramFromCostFunctionWithWeightsBinned(psi, diag, max_cut_value, bin_width), "" );
#endif
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard QAOA_FEATURES_TEST_HPP
