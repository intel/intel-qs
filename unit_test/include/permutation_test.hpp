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
};

//////////////////////////////////////////////////////////////////////////////
// Test macros:

TEST_F(PermutationTest, InitializationToTrivialPermutation)
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
  map_ = {1, 2, 0, 3, 5, 4};
  permutation.SetNewPermutation(map_);
  for (unsigned i=0; i<num_bits_; ++i)
      ASSERT_EQ(permutation[i], map_[i]);
  // Inverse map.
  std::vector<std::size_t> imap = permutation.imap;
  for (unsigned i=0; i<num_bits_; ++i)
      ASSERT_EQ(imap[map_[i]], i);
}

//////////////////////////////////////////////////////////////////////////////

TEST_F(PermutationTest, Lin2Perm)
{
  num_bits_ = 3;
  Permutation permutation(num_bits_);
  map_ = {0, 1, 2};
  permutation.SetNewPermutation(map_);
  std::size_t dim = (0x1 << num_bits_);
  ASSERT_EQ(dim, std::pow(2, num_bits_));
  for (std::size_t v=0; v<dim; ++v)
      ASSERT_EQ(v, permutation.lin2perm_(v));

  // When the permutation is only 1- and 2-cycles, map==imap.
  map_ = {1, 0, 2};
  permutation.SetNewPermutation(map_);
  for (std::size_t v=0; v<dim; ++v)
      std::cout << "v=" << v << " --> " << permutation.lin2perm_(v) << "\n";

  // Permutation with a 3-cycle.
  map_ = {1, 2, 0};
  permutation.SetNewPermutation(map_);
  std::size_t u;
  for (std::size_t v=0; v<dim; ++v)
  {
      u = permutation.lin2perm_(v);
      std::cout << "v="    << v << "=" << permutation.dec2bin(v, num_bits_)
                << " --> " << u << "=" << permutation.dec2bin(u, num_bits_) << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////
// Here we test the class Permutation using a vector state that is much simpler than
// a QubitRegister state. For this reason, we define the simple class:

namespace utest
{
  class State
  {
   public:
    std::string name;
    Permutation p;
    std::vector<ComplexDP> state;
  
    // Creator of the state.
    State(Permutation p_, std::string name_) : p(p_), name(name_)
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
      for (std::size_t i = 0; i < pnew.num_qubits; i++) {
        map[i] = pold.map[pnew.imap[i]];
        // printf("%d ", map[i]);
      }
      // printf("\n");
  
      __int64 t0 = __rdtsc();
      double s0 = sec();
      for (std::size_t i = 0; i < state.size(); i++) {
#if 0
        std::size_t to1 = perm(i, &(map[0]), p.num_qubits);
        std::size_t to2 = pnew.perm2lin_(pold.lin2perm_(i));
        std::size_t to = pold.bin2dec(pnew.perm2lin(pold.lin2perm(i)));
        assert(to == to1);
        assert(to == to2);
        state_new[to] = state[i];
#else
        std::size_t to_ = perm(i, &(map[0]), p.num_qubits);  // pnew.perm2lin_(pold.lin2perm_(i));
        state_new[to_] = state[i];
#endif
      }
      __int64 t1 = __rdtsc();
      double s1 = sec();
      double bw = D(state.size()) * D(sizeof(state[0])) * 2.0 / D(s1 - s0);
      printf("cycles per shuffle: %.2lf bw=%.2lf GB/s\n", D(t1 - t0) / D(state.size()), bw / 1e9);
  
      p = pnew;
      state = state_new;
    }
  
    // Print the state vector.
    void print()
    {
      printf("name::%s %s\n", name.c_str(), p.GetMapStr().c_str());
      for (std::size_t i = 0; i < state.size(); i++) {
        printf("%s {%lf %lf}\n", p.lin2perm(i).c_str(), real(state[i]), imag(state[i]));
      }
    }

    // Comparison between state vectors.
    bool operator==(const State &rhs)
    {
      assert((const std::size_t)state.size() == rhs.state.size());
      for (std::size_t i = 0; i < state.size(); i++) {
        if (state[i] != rhs.state[i]) return false;
      }
  
      return true;
    }
  };
}

//////////////////////////////////////////////////////////////////////////////
// Unit test using the ad-hoc utest::State object.

TEST_F(PermutationTest, PermutationOfAdHocState)
{
  num_bits_ = 3;
#if 1
  Permutation p({2, 0, 1});
  p.prange();
  utest::State s(Permutation({0, 1, 2}), "s");
  s.print();
  s.permute(p);
  s.print();
  s.permute(Permutation({0, 1, 2}));
  s.print();
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
