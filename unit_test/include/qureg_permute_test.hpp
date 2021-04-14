/// @file qureg_permute_test.hpp
/// @brief Unit test for the @c QubitRegister class.

#ifndef	QUREG_PERMUTE_TEST_HPP
#define	QUREG_PERMUTE_TEST_HPP

#include "../../include/qureg.hpp"
#include "../../include/permutation.hpp"

#include <math.h>
#include <string>
#include <vector>

//////////////////////////////////////////////////////////////////////////////
// Test fixture class: permutation for qubit ordering
//////////////////////////////////////////////////////////////////////////////

class QuregPermuteTest : public ::testing::Test
{
 protected:

  QuregPermuteTest()
  { }

  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (iqs::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();
    // Number of local qubits.
    M_ = num_qubits_ - iqs::ilog2(iqs::mpi::Environment::GetStateSize());
    
  }

  const std::size_t num_qubits_ = 20;
  std::size_t M_;
  std::vector<std::size_t>  map_;
  std::vector<std::size_t> imap_;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(QuregPermuteTest, InitialPermutation)
{
  // Initial state |0010...0>
  std::size_t index = 4;
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", index);
  unsigned position;
  // Verify that the qubit permutation in psi is the identity.
  for (unsigned qubit=0; qubit<num_qubits_; ++qubit)
  {
      position = qubit;
      EXPECT_EQ(position, (*psi.qubit_permutation)[qubit]);
  }
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(QuregPermuteTest, PermuteLocalQubits)
{
  // Initial state |0000...0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", 0);
  psi[0]=0;
  // Set a few amplitudes:
  psi.SetGlobalAmplitude(2, ComplexDP(2,0)); // 2^1
  psi.SetGlobalAmplitude(8, ComplexDP(8,0)); // 2^3
  psi.SetGlobalAmplitude(1024, ComplexDP(1024,0)); // 2^10
  psi.SetGlobalAmplitude(4096, ComplexDP(4096,0)); // 2^12

  std::size_t index = (1UL << (num_qubits_-1)) + 512;
  psi.SetGlobalAmplitude(index, ComplexDP(0,1));
  index = (1UL << (num_qubits_-1)) + 1024;
  psi.SetGlobalAmplitude(index, ComplexDP(0,2));

  index = (1UL << (num_qubits_-2)) + 512;
  psi.SetGlobalAmplitude(index, ComplexDP(0,3));
  index = (1UL << (num_qubits_-2)) + 1024;
  psi.SetGlobalAmplitude(index, ComplexDP(0,4));
  //       |--|     |--------|            |-------|
  map_  = {1, 0, 2, 6, 4, 5, 3, 7, 8, 9, 12, 11, 10, 13, 14, 15, 16, 17, 18, 19};
  if (12 > M_)
      GTEST_SKIP();
  else
  {
      psi.PermuteLocalQubits(map_, "direct");
      if (iqs::mpi::Environment::GetStateRank() == 0)
      {
          ASSERT_EQ(psi[2].real(), 0);
          ASSERT_EQ(psi[1].real(), 2);
          ASSERT_EQ(psi[8].real(), 0);
          ASSERT_EQ(psi[64].real(), 8);
          ASSERT_EQ(psi[1024].real(), 4096);
          ASSERT_EQ(psi[4096].real(), 1024);
      }
      else if (M_<num_qubits_ &&
               iqs::mpi::Environment::GetStateRank() == (1UL << (num_qubits_-1-M_)))
      {
          ASSERT_EQ(psi[512].imag(), 1);
          ASSERT_EQ(psi[1024].imag(), 0);
          ASSERT_EQ(psi[4096].imag(), 2);
      }
      else if (M_<num_qubits_-1 &&
               iqs::mpi::Environment::GetStateRank() == (1UL << (num_qubits_-2-M_)))
      {
          ASSERT_EQ(psi[512].imag(), 3);
          ASSERT_EQ(psi[1024].imag(), 0);
          ASSERT_EQ(psi[4096].imag(), 4);
      }
      else
      {
          ASSERT_EQ(psi[2].real(), 0);
          ASSERT_EQ(psi[1].real(), 0);
          ASSERT_EQ(psi[8].real(), 0);
          ASSERT_EQ(psi[64].real(), 0);
          ASSERT_EQ(psi[512].real(), 0);
          ASSERT_EQ(psi[1024].real(), 0);
          ASSERT_EQ(psi[4096].real(), 0);
      }
  }
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(QuregPermuteTest, PermuteGlobalQubits)
{
  // Initial state |0000...0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", 0);
  psi[0]=0;
  // Set a few amplitudes:
  std::vector<std::size_t> indices = {             2,              8,           1024,           4096};
  std::vector<ComplexDP>   values  = {ComplexDP(1,0), ComplexDP(2,0), ComplexDP(3,0), ComplexDP(4,0)};

  indices.push_back((1UL << (num_qubits_-1)) + 512);
  values.push_back(ComplexDP(0,1));
  indices.push_back((1UL << (num_qubits_-1)) + 1024);
  values.push_back(ComplexDP(0,2));

  indices.push_back((1UL << (num_qubits_-2)) + 512);
  values.push_back(ComplexDP(0,3));
  indices.push_back((1UL << (num_qubits_-2)) + 1024);
  values.push_back(ComplexDP(0,4));
  indices.push_back((1UL << (num_qubits_-2)) + 4096);
  values.push_back(ComplexDP(0,5));

  indices.push_back((1UL << (num_qubits_-3)) + 0);
  values.push_back(ComplexDP(0,6));
  indices.push_back((1UL << (num_qubits_-3)) + 4096);
  values.push_back(ComplexDP(0,7));

  indices.push_back((1UL << (num_qubits_-4)) + 0);
  values.push_back(ComplexDP(0,8));
  indices.push_back((1UL << (num_qubits_-4)) + 512);
  values.push_back(ComplexDP(0,9));
  indices.push_back((1UL << (num_qubits_-4)) + 4096);
  values.push_back(ComplexDP(0,10));

  ASSERT_EQ(indices.size(), values.size());
  for (unsigned i=0; i<indices.size(); ++i)
      psi.SetGlobalAmplitude(indices[i], values[i]); 

  std::size_t myrank = iqs::mpi::Environment::GetStateRank();
  std::size_t num_global_qubits = num_qubits_-M_;

  if (num_global_qubits >= 4)
  {
      //                                                              |-------|---|
      map_  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 17, 16, 18};
      psi.PermuteGlobalQubits(map_, "direct");
      for (unsigned i=0; i<indices.size(); ++i)
          ASSERT_EQ(psi.GetGlobalAmplitude(indices[i]), values[i]);
  }
  iqs::mpi::StateBarrier();

  if (num_global_qubits >= 3)
  {
      //                                                                  |-------|
      map_  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 19, 18, 17};
      psi.PermuteGlobalQubits(map_, "direct");
      for (unsigned i=0; i<indices.size(); ++i)
          ASSERT_EQ(psi.GetGlobalAmplitude(indices[i]), values[i]);
  }
  iqs::mpi::StateBarrier();

  if (num_global_qubits >= 2)
  {
      //                                                                      |---|
      map_  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 18};
      psi.PermuteGlobalQubits(map_, "direct");
      for (unsigned i=0; i<indices.size(); ++i)
          ASSERT_EQ(psi.GetGlobalAmplitude(indices[i]), values[i]);
  }
  iqs::mpi::StateBarrier();
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(QuregPermuteTest, PermuteQubits)
{
  // Initial state |0000...0>
  iqs::QubitRegister<ComplexDP> psi (num_qubits_, "base", 0);
  psi[0]=0;
  // Set a few amplitudes:                       2^1             2^3            2^10            2^12
  std::vector<std::size_t> indices = {             2,              8,           1024,           4096};
  std::vector<ComplexDP>   values  = {ComplexDP(1,0), ComplexDP(2,0), ComplexDP(3,0), ComplexDP(4,0)};

  indices.push_back( (1UL << (num_qubits_-1)) +  512);
  values.push_back(ComplexDP(0,1));
  indices.push_back( (1UL << (num_qubits_-1)) + 1024);
  values.push_back(ComplexDP(0,2));

  indices.push_back( (1UL << (num_qubits_-2)) +  512);
  values.push_back(ComplexDP(0,3));
  indices.push_back( (1UL << (num_qubits_-2)) + 1024);
  values.push_back(ComplexDP(0,4));

  ASSERT_EQ(indices.size(), values.size());
  for (unsigned i=0; i<indices.size(); ++i)
      psi.SetGlobalAmplitude(indices[i], values[i]); 

  // The map corresponds to several k-cycles.
  //    (0,4,6,16) (3,5,19) (8,10) (18,17,9) 

  //      {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}  // qubit
  map_  = {  4,  1,  2,  5,  6, 19, 16,  7, 10, 18,  8, 11, 12, 13, 14, 15,  0,  9, 17,  3}; // position

  // Code to print extra information.
#if 0
  std::vector<std::size_t> int_1_imap, int_2_imap;
  std::size_t M = num_qubits_ - iqs::ilog2(iqs::mpi::Environment::GetStateSize());
  psi.qubit_permutation->ObtainIntemediateInverseMaps(map_, M, int_1_imap, int_2_imap);
  if (iqs::mpi::Environment::GetStateRank()==0)
  {
    std::cout << "original inverse map:\n";
    for (unsigned pos=0; pos<num_qubits_; ++pos)
        std::cout << psi.qubit_permutation->Find(pos) << "\t";
    std::cout << "\nfirst intermediate imap:\n";
    for (unsigned pos=0; pos<num_qubits_; ++pos)
        std::cout << int_1_imap[pos] << "\t";
    std::cout << "\nsecond intermediate imap:\n";
    for (unsigned pos=0; pos<num_qubits_; ++pos)
        std::cout << int_2_imap[pos] << "\t";
    std::cout << "\nfinal inverse map:\n";
    Permutation new_permutation(map_, "direct");
    for (unsigned pos=0; pos<num_qubits_; ++pos)
        std::cout << new_permutation.Find(pos) << "\t";
    std::cout << "\n";
  }
#endif

  psi.PermuteQubits(map_, "direct");
  // Implemented in three steps:
  // - shuffle  local qubits only
  // - shuffle global qubits only
  // - pairwise exchange of local-global qubits

  for (unsigned i=0; i<indices.size(); ++i)
      ASSERT_EQ(psi.GetGlobalAmplitude(indices[i]), values[i]);
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard QUREG_PERMUTE_TEST_HPP
