#ifndef PERMUTATION_HPP
#define PERMUTATION_HPP

#include <map>
#include <numeric>
#include <vector>

#include "utils.hpp"
#include "conversion.hpp"

inline std::size_t perm(std::size_t v, std::size_t *map, std::size_t num_qubits)
{
  std::size_t v_ = 0;
  for (std::size_t i = 0; i < num_qubits; i++)
      v_ = v_ | (((v & (1 << map[i])) >> map[i]) << i);
  return v_;
}

/////////////////////////////////////////////////////////////////////////////////////////
// declaration and definition of class Permutation
/////////////////////////////////////////////////////////////////////////////////////////

/// @class Permutation 
/// Define the order of the qubits in the representation of teh quantum state
///
/// The N-qubit quantum state |psi> is stored as a 2^N complex vector.
/// The global index i corresponds to the entry:
///   state[i] = |i0>_dq0 * |i1>_dq1 * |i2>_dq2 * ...
/// where ik is the k-th bit of i in its N-bit representation (with i0 being the
/// least significant bit) and dqk corresponds to the k-th data qubit.
///
/// By default, data qubits and program qubits correspond trivially:
///   pq(0) --> dq(0)      ,  pq(1) --> dq(1)      ,  ...
/// However this does not have to be the case. Using the Permutation object, one has:
///   pq(0) --> dq(map[0]) ,  pq(1) --> dq(map[1]) ,  ...
/// and its inverse:
///   dq(0) --> pq(imap[0]),  dq(1) --> pq(imap[1]),  ...
/// 
///
/// @var map: map
/// @var imap: inverse map

class Permutation
{
 public:
  std::vector<std::size_t> map, imap;
  std::size_t num_qubits;

/////////////////////////////////////////////////////////////////////////////////////////
// overload operator [] to address single entries of the permutation map

  unsigned operator[](std::size_t i)
  {
    assert(i <= num_qubits);
    return (unsigned) map[i];
  }

  unsigned operator[](unsigned i)
  {
    assert(i <= num_qubits);
    return map[i];
  }

  int operator[](int i)
  {
    assert(i <= num_qubits);
    return (int)map[i];
  }

  std::size_t size() {return map.size();}

/////////////////////////////////////////////////////////////////////////////////////////
// Return a desription of map or imap as a string with space separated entries.

  std::string GetMapStr()
  {
    std::string s;
    for (std::size_t i = 0; i < map.size(); i++)
        s += " " + qhipster::toString(map[i]);
    return s;
  }

  std::string GetImapStr()
  {
    std::string s;
    for (std::size_t i = 0; i < imap.size(); i++)
        s += " " + qhipster::toString(imap[i]);
    return s;
  }

/////////////////////////////////////////////////////////////////////////////////////////
// Creator of the Permutation object.

  Permutation(std::size_t num_qubits)
  {
    this->num_qubits = num_qubits;
    std::vector<std::size_t> m(num_qubits);
    // Initialize the trivial permutation: {0,1,2,... num_qubits-1}
    iota(m.begin(), m.end(), 0);
    // for(auto i:m) printf("%d ", i); printf("\n");
    // exit(0);
    SetNewPermutation(m);
  }

  Permutation(std::vector<std::size_t> m)
  {
    num_qubits = m.size();
    SetNewPermutation(m);
  }

/////////////////////////////////////////////////////////////////////////////////////////
// Find data qubit associated with program qubit 'position' .
// If it exists, Find(pos)=imap(pos).

  std::size_t Find(std::size_t position)
  {
    bool found = false;
    std::size_t i;
    for (i = 0; i < map.size(); i++)
    {
        if (map[i] == position)
         {
            found = true;
            break;
        }
    }
    assert(found);
    return i;
  }

/////////////////////////////////////////////////////////////////////////////////////////
// Set new permutation after validating the input as a valid permutation.

  void SetNewPermutation(std::vector<std::size_t> m)
  {
    map = m;
    // Check consistency of map.
    assert (map.size() == this->num_qubits);
    std::vector<bool> exist(map.size(), 0);
    for (auto &m : map)
        exist[m] = 1;
    for (auto e : exist)
        assert(e > 0);

    // Compute inverse map.
    imap = map;
    for (std::size_t i = 0; i < map.size(); i++)
        imap[map[i]] = i;
  }

/////////////////////////////////////////////////////////////////////////////////////////
// Utility functions to transform between binary, decimal and string representation.

  std::string dec2bin(std::size_t in, std::size_t num_bits)
  {
    std::string s;
    s.resize(num_bits);
    for (std::size_t i = 0; i < num_bits; i++)
    {
        s[i] = (in & 0x1) == 0 ? '0' : '1';
        in = in >> 1;
    }
    return s;
  }

  std::size_t bin2dec(std::string in)
  {
    std::size_t v = 0;
    for (std::size_t i = 0; i < in.size(); i++)
        v += UL(1 << i) * UL(in[i] - '0');
    return v;
  }

/////////////////////////////////////////////////////////////////////////////////////////
// Understanding bitwise operations:
//   v & (1 << map[i]))   ---->   zero bitstring, apart from the bit value corresponding
//                                to the data qubit associated with program qubit i
//                                the bit value is in position map[i]
//
//   [ans] >> map[i]      ---->   the above bit is moved to position 0
//
//   [ans] << i           ---->   the above bit is moved to position i
//
// Notice that when one loops over i, the bit value always ends up in a different position.
// Therefore, taking the 'bitwise or' will simply set the i-th bit once at a time.

  // Transform index v in [0,2^N[ from the data representation to the program representation.
  // FIXME: the name is poorly chosen, possibly opposite to intention
  inline std::size_t lin2perm_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < num_qubits; i++)
        v_ = v_ | (((v & (((std::size_t) 1) << map[i])) >> map[i]) << i);
    return v_;
  }

  // Transform index v in [0,2^N[ from the program representation to the data representation.
  // FIXME: the name is poorly chosen, possibly opposite to intention
  inline std::size_t perm2lin_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < num_qubits; i++)
        v_ = v_ | (((v & (((std::size_t)1) << imap[i])) >> imap[i]) << i);
    return v_;
  }

  std::string lin2perm(std::size_t v)
  {
    std::string s = dec2bin(v, num_qubits), sp(s);
    for (std::size_t i = 0; i < s.size(); i++) sp[i] = s[map[i]];
    if (0)
    {
        std::size_t v_ = 0;
        for (std::size_t i = 0; i < num_qubits; i++)
            v_ = v_ | (((v & (((std::size_t)1) << map[i])) >> map[i]) << i);
        printf("sp=%s new:%s\n", sp.c_str(), dec2bin(v_, num_qubits).c_str());
    }
    return sp;
  }

  std::string lin2perm(std::string s)
  {
    std::string sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[map[i]];
    return sp;
  }

  std::string perm2lin(std::size_t v)
  {
    std::string s = dec2bin(v, num_qubits), sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[imap[i]];
    return sp;
  }
  std::string perm2lin(std::string s)
  {
    std::string sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[imap[i]];
    return sp;
  }

  void prange()
  {
#if 0
    printf("map:  ");
    for(auto & i : map) printf("%d", i);
    printf("\n");
    printf("imap: ");
    for(auto & i : imap) printf("%d", i);
    printf("\n");
    for(std::size_t i = 0; i < (1 << num_qubits); i++)
    {
       printf("%s(%lld)\n", dec2bin(i, num_qubits).c_str(), i);
    }

    printf("\n");
    for(std::size_t i = 0; i < (1 << num_qubits); i++)
    {
       printf("%s(%lld)\n", lin2perm(i).c_str(), bin2dec(lin2perm(i)));
    }

    printf("\n");
#endif

    for (std::size_t i = 0; i < (UL(1) << num_qubits); i++)
    {
#if 0
       printf("%s ==> %s\n", lin2perm(i).c_str(), 
               perm2lin(lin2perm(i)).c_str());
#else
      printf("map(%3lu) = %3lu\n", i, bin2dec(lin2perm(i)));
#endif
    }
  }
};

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// FIXME FIXME: the code below was probably written by Misha to test the permutation class.
//              The MPI part has not been stripped from it.

#if defined(MAIN)
class State
{
 public:
  std::string name;
  Permutation p;
  std::vector<ComplexType> state;

  State(Permutation p_, std::string name_) : p(p_), name(name_)
  {
    state.resize(1 << p.num_qubits);
#pragma omp parallel for
    for (std::size_t i = 0; i < state.size(); i++) state[i] = {D(i % 3), D(i % 7)};
  }

/////////////////////////////////////////////////////////////////////////////////////////
  void permute(Permutation pnew)
  {
    Permutation pold = p;

    assert(pnew.num_qubits == pold.num_qubits);
    std::vector<ComplexType> state_new(state.size(), 0);

    // printf("map: %s imap: %s\n", pnew.GetMapStr().c_str(), pold.GetImapStr().c_str());
    std::vector<std::size_t> map(pnew.num_qubits, 0);
    for (std::size_t i = 0; i < pnew.num_qubits; i++) {
      map[i] = pold.map[pnew.imap[i]];
      // printf("%d ", map[i]);
    }
    // printf("\n");

    __int64 t0 = __rdtsc();
    double s0 = sec();
#pragma omp parallel for
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

/////////////////////////////////////////////////////////////////////////////////////////
  void print()
  {
    printf("name::%s %s\n", name.c_str(), p.GetMapStr().c_str());
    for (std::size_t i = 0; i < state.size(); i++) {
      printf("%s {%lf %lf}\n", p.lin2perm(i).c_str(), real(state[i]), imag(state[i]));
    }
  }
/////////////////////////////////////////////////////////////////////////////////////////
  bool operator==(const State &rhs)
  {
    assert((const std::size_t)state.size() == rhs.state.size());
    for (std::size_t i = 0; i < state.size(); i++) {
      if (state[i] != rhs.state[i]) return false;
    }

    return true;
  }
};


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


std::size_t main(std::size_t argc, char **argv)
{
  std::size_t num_qubits = 3, num_threads = 1;
  if (argc != 3)
  {
      fprintf(stderr, "usage: %s <num_qubits> <num_threads>\n", argv[0]);
      exit(1);
  }
  else
  {
      num_qubits = atoi(argv[1]);
      num_threads = atoi(argv[2]);
  }
  initomp(num_threads);
#if 1
  Permutation p({2, 0, 1});
  p.prange();
  State s(Permutation({0, 1, 2}), "s");
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
#endif

#endif	// header guard PERMUTATION_HPP
