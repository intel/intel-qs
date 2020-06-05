/// @file permutation_test.hpp
/// @brief Unit test for the @c Permutation class.

#ifndef	PERMUTATION_TEST_HPP
#define	PERMUTATION_TEST_HPP

#include "../../include/permutation.hpp"

#include <math.h>
#include <string>
#include <vector>

//////////////////////////////////////////////////////////////////////////////
// Test fixture class: permutation for qubit ordering
//////////////////////////////////////////////////////////////////////////////

class PermutationTest : public ::testing::Test
{
 protected:

  PermutationTest()
  { }

  void SetUp() override
  {
    // All tests are skipped if the rank is dummy.
    if (qhipster::mpi::Environment::IsUsefulRank() == false)
      GTEST_SKIP();
  }

  std::size_t num_bits_;
  std::vector<std::size_t> map_;
  std::vector<std::size_t> imap_;
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(PermutationTest, InitializationToIdentitylPermutation)
{
  num_bits_ = 6;
  Permutation permutation(num_bits_);
  ASSERT_EQ(permutation.num_qubits, num_bits_);

  // The permutation is initialized to the trivial one.
  map_ = {0, 1, 2, 3, 4, 5};
  for (unsigned i=0; i<num_bits_; ++i)
      ASSERT_EQ(permutation[i], map_[i]);
  // Inverse map.
  std::vector<std::size_t> imap = permutation.imap;
  for (unsigned i=0; i<num_bits_; ++i)
      ASSERT_EQ(imap[map_[i]], i);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(PermutationTest, BasicUse)
{
  num_bits_ = 6;
  Permutation permutation(num_bits_);
  ASSERT_EQ(permutation.num_qubits, num_bits_);

  // Map from program qubits to data qubits.
  map_  = {1, 2, 0, 3, 5, 4};
  imap_ = {2, 0, 1, 3, 5, 4}; // hard-coded inverse of above map_
  permutation.SetNewPermutationFromMap(map_, "direct");
  for (unsigned i=0; i<num_bits_; ++i)
      ASSERT_EQ(permutation[i], map_[i]);
  // Inverse map.
  std::vector<std::size_t> imap = permutation.imap;
  for (unsigned i=0; i<num_bits_; ++i)
  {
      ASSERT_EQ(imap[map_[i]], i       );
      ASSERT_EQ(imap[i]      , imap_[i]);
  }

  // Update the permutation using map_ as the inverse map.
  permutation.SetNewPermutationFromMap(map_, "inverse");
  for (unsigned i=0; i<num_bits_; ++i)
  {
      ASSERT_EQ(permutation.imap[i],  map_[i] );
      ASSERT_EQ(permutation.map[i] , imap_[i] );
  }
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(PermutationTest, ExchangeTwoElements)
{

  num_bits_ = 6;
  Permutation permutation(num_bits_);
  map_ = {1, 2, 0, 3, 5, 4};
  permutation.SetNewPermutationFromMap(map_, "direct");

  // Emulate a SWAP between program qubits: 2,4
  unsigned element_1 = 2;
  unsigned element_2 = 4;
  unsigned position_1 = permutation[element_1];
  unsigned position_2 = permutation[element_2];
  permutation.ExchangeTwoElements(element_1, element_2);
  ASSERT_EQ(map_[element_1], position_1);
  map_[element_1] = position_2;
  map_[element_2] = position_1;
  for (unsigned i=0; i<num_bits_; ++i)
      ASSERT_EQ(permutation[i], map_[i]);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(PermutationTest, Data2Program)
{
  num_bits_ = 3;
  Permutation permutation(num_bits_);
  map_ = {0, 1, 2};
  permutation.SetNewPermutationFromMap(map_, "direct");
  std::size_t dim = (0x1 << num_bits_);
  ASSERT_EQ(dim, std::pow(2, num_bits_));
  for (std::size_t v=0; v<dim; ++v)
      ASSERT_EQ(v, permutation.data2program_(v));

  // Basic permutation formed by a 2-cycle.
  map_ = {1, 0, 2};
  std::vector<std::size_t> expected_map = {0,2,1,3, 4,6,5,7};	// hardcoded example
  permutation.SetNewPermutationFromMap(map_, "direct");
  for (std::size_t v=0; v<dim; ++v)
  {
      //std::cout << "v=" << v << " --> " << permutation.program2data_(v) << "\n";
      ASSERT_EQ(permutation.program2data_(v), expected_map[v]);
      // When the permutation is only 1- and 2-cycles, map==imap.
      ASSERT_EQ(permutation.data2program_(v), permutation.program2data_(v));
  }

  // Permutation with a 3-cycle.
  map_ = {1, 2, 0};
  // Explicit computation:
  //   program rep  ---->   data rep
  //   '2' '1' '0'          '2' '1' '0'
  //    0   0   0            0   0   0
  //    0   0   1            0   1   0    (since program bit '0' is mapped to data bit '1')
  //    0   1   0            1   0   0    (since program bit '1' is mapped to data bit '2')
  //    0   1   1            1   1   0    (see above)
  //       ...                  ...
  expected_map = {0,2,4,6, 1,3,5,7};	// hardcoded example
  permutation.SetNewPermutationFromMap(map_, "direct");
  std::size_t u;
  //std::cout << "program representation --> data representation\n";
  for (std::size_t v=0; v<dim; ++v)
  {
      u = permutation.program2data_(v);
      //std::cout << "v="    << v << "=" << permutation.dec2bin(v, num_bits_)
      //          << " --> " << u << "=" << permutation.dec2bin(u, num_bits_) << "\n";
      ASSERT_EQ(u, expected_map[v]);
  }
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(PermutationTest, ObtainIntemediateInverseMaps)
{
  // Hard-coded example, thought for N=8 qubits of which M=5 are local.
  //
  //  map: {0,1,2,3,4,5,6,7} --> {1,5,2,3,6,0,7,4}
  // imap: {0,1,2,3,4,5,6,7} --> {5,0,2,3,7,1,4,6}
  //
  // Adding intermediate maps:
  // imap: {0,1,2,3,4,5,6,7} --> {1,0,2,3,4,5,6,7} --> {1,0,2,3,4,5,7,6} --> {5,0,2,3,7,1,4,6}
  //           original         only-local updates    only-global updates         target
  num_bits_ = 8;
  std::size_t M=5;
  map_ = {0, 1, 2, 3, 4, 5, 6, 7};
  Permutation permutation(map_, "direct");
  std::vector<std::size_t> target_map = {1, 5, 2, 3, 6, 0, 7, 4};
  std::vector<std::size_t> expected_int_1_imap = {1, 0, 2, 3, 4, 5, 6, 7};
  std::vector<std::size_t> expected_int_2_imap = {1, 0, 2, 3, 4, 5, 7, 6};

  std::vector<std::size_t> int_1_imap, int_2_imap;
  permutation.ObtainIntemediateInverseMaps(target_map, M, int_1_imap, int_2_imap);
  for (std::size_t pos=0; pos<num_bits_; ++pos)
  {
      ASSERT_EQ(int_1_imap[pos], expected_int_1_imap[pos]);
      ASSERT_EQ(int_2_imap[pos], expected_int_2_imap[pos]);
  }

  // Hard-coded example, thought for N=8 qubits of which M=5 are local.
  //
  //  map: {1,5,2,3,6,0,7,4} --> {1,6,7,4,3,2,0,5}
  // imap: {5,0,2,3,7,1,4,6} --> {6,0,5,4,3,7,1,2}
  //
  // Adding intermediate maps:
  // imap: {5,0,2,3,7,1,4,6} --> {2,0,5,7,3,1,4,6} --> {2,0,5,7,3,4,1,6} --> {6,0,5,4,3,7,1,2}
  //           original         only-local updates    only-global updates         target
  map_ = {1, 5, 2, 3, 6, 0, 7, 4};
  permutation.SetNewPermutationFromMap(map_, "direct");
  target_map = {1, 6, 7, 4, 3, 2, 0, 5};
  expected_int_1_imap = {2,0,5,7,3,1,4,6};
  expected_int_2_imap = {2,0,5,7,3,4,1,6};

  permutation.ObtainIntemediateInverseMaps(target_map, M, int_1_imap, int_2_imap);
  for (std::size_t pos=0; pos<num_bits_; ++pos)
  {
      ASSERT_EQ(int_1_imap[pos], expected_int_1_imap[pos]);
      ASSERT_EQ(int_2_imap[pos], expected_int_2_imap[pos]);
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Here we test the class Permutation using a vector state that is much simpler than
// a QubitRegister state. For this reason, we define the simple class:

namespace utest
{
  // Update decimal representation of v such that bit map[i] is moved to bit i.
  inline std::size_t perm(std::size_t v, std::size_t *map, std::size_t num_qubits)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < num_qubits; i++)
        v_ = v_ | (((v & (1 << map[i])) >> map[i]) << i);
    return v_;
  }


  class State
  {
   public:

    std::string name;
    Permutation p;
    std::vector<ComplexDP> state;
    bool do_print;
  
    // Creator of the state.
    State(Permutation p_, std::string name_, bool do_print_=true)
      : p(p_), name(name_), do_print(do_print_)
    {
        state.resize(1 << p.num_qubits);
        for (std::size_t i = 0; i < state.size(); i++)
            state[i] = {D(i % 3), D(i % 16)};
    }
 
    // Permute the order of the entries according to the Permutation pnew.
    void permute(Permutation pnew)
    {
        Permutation pold = p;
        assert(pnew.num_qubits == pold.num_qubits);
        std::vector<ComplexDP> state_new(state.size(), 0);
    
        // printf("map: %s imap: %s\n", pnew.GetMapStr().c_str(), pold.GetImapStr().c_str());
        std::vector<std::size_t> map(pnew.num_qubits, 0);
        for (std::size_t i = 0; i < pnew.num_qubits; i++)
        {
            //                   pnew.imap[i]   -->  program bit associated to new data bit 'i'
            //          pold.map[pnew.imap[i]]  -->  old data bit associated to new data bit 'i'
            // map[i] = pold.map[pnew.imap[i]]  -->  map['new data bit i'] --> associated old data bit
            map[i] = pold.map[pnew.imap[i]];
            // printf("%d ", map[i]);
        }
        // printf("\n");
    
        __int64 t0 = __rdtsc();
        double s0 = sec();
        for (std::size_t i = 0; i < state.size(); i++)
        {
            // 'i'  is the index in the old data representation.
            // 'to' is the corresponding index in the new data representation
            std::size_t to = perm(i, &(map[0]), p.num_qubits);
            assert(to == pnew.program2data_(pold.data2program_(i)));
            assert(to == pold.bin2dec(pnew.program2data(pold.data2program(i))));
            assert(to == pnew.bin2dec(pnew.program2data(pold.data2program(i))));
            state_new[to] = state[i];
        }
        __int64 t1 = __rdtsc();
        double s1 = sec();
        double bw = D(state.size()) * D(sizeof(state[0])) * 2.0 / D(s1 - s0);
        if (do_print) printf("cycles per shuffle: %.2lf bw=%.2lf GB/s\n", D(t1 - t0) / D(state.size()), bw / 1e9);
    
        p = pnew;
        state = state_new;
    }
  
    // Print the state vector.
    void print()
    {
        if (do_print==false) return;
        printf("name::%s %s\n", name.c_str(), p.GetMapStr().c_str());
        printf("data  program  state[data]\n");
        for (std::size_t i = 0; i < state.size(); i++) {
            printf("%i     %s      {%lf %lf}\n", i, p.data2program(i).c_str(), real(state[i]), imag(state[i]));
      }
    }

    // Comparison between state vectors.
    bool operator==(const State &rhs)
    {
        assert((const std::size_t)state.size() == rhs.state.size());
        for (std::size_t i = 0; i < state.size(); i++)
        {
            if (state[i] != rhs.state[i])
                return false;
        }
        return true;
    }
  };
}

//////////////////////////////////////////////////////////////////////////////
// Unit test using the ad-hoc utest::State object.

TEST_F(PermutationTest, PermutationOfSpecializedStateClass)
{
  num_bits_ = 3;

  utest::State s(Permutation({0, 1, 2}), "s", false);
  s.print();
  // Original order:
  for (std::size_t i=0; i<(1<<num_bits_); ++i)
      ASSERT_DOUBLE_EQ(s.state[i].imag(), double(i));

#if 1
  Permutation p({2, 0, 1});
//  p.Print();
  s.permute(p);
  s.print();
  // New order, explicit computation:
  //   program rep  ---->   data rep
  //   '2' '1' '0'          '2' '1' '0'
  //    0   0   0            0   0   0
  //    0   0   1            1   0   0    (since program bit '0' is mapped to data bit '2')
  //    0   1   0            0   0   1    (since program bit '1' is mapped to data bit '0')
  //    0   1   1            1   0   1    (see above)
  //       ...                  ...
  std::vector<std::size_t> expected_map  = {0,4,1,5, 2,6,3,7};	// hardcoded example
  std::vector<std::size_t> expected_imap = {0,2,4,6, 1,3,5,7};	// hardcoded example
  for (std::size_t i=0; i<(1<<num_bits_); ++i)
  {
      // From data to program.
      ASSERT_DOUBLE_EQ(s.state[i].imag(), expected_imap[i]);
      // From program to data.
      ASSERT_DOUBLE_EQ(i, s.state[expected_map[i]].imag());
  }

  s.permute(Permutation({0, 1, 2}));
  s.print();
  // Back to original order:
  for (std::size_t i=0; i<(1<<num_bits_); ++i)
      ASSERT_DOUBLE_EQ(s.state[i].imag(), double(i));

#else
#if 0
   State s1(Permutation({2, 1, 0}), "s"), s2(s1);
   s1.permute(Permutation({2,0,1}));
   s2.permute(Permutation({2,0,1}));
   assert(s1 == s2);
   printf("SUCCESS!\n");
#else
  std::vector<std::size_t> map(num_qubits, 0);
  iota(map.begin(), map.end(), 0);
  State s1(Permutation(map), "s1"), s2(s1);
  s2.permute(Permutation(map));
  s2.permute(Permutation(map));
  s2.permute(Permutation(map));
#endif
#endif
}

//////////////////////////////////////////////////////////////////////////////

#endif	// header guard PERMUTATION_TEST_HPP
