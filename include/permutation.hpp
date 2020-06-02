/// @file permutation.hpp
/// @brief Declare the @c Permutation class.

#ifndef PERMUTATION_HPP
#define PERMUTATION_HPP

#include <map>
#include <numeric>
#include <vector>

#include "utils.hpp"
#include "conversion.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// declaration and definition of class Permutation
/////////////////////////////////////////////////////////////////////////////////////////

/// @class Permutation (described in the context of IQS).
/// Define the order of the qubits in the representation of the quantum state.
///
/// The N-qubit quantum state |psi> is stored as a 2^N complex vector.
/// The global index i corresponds to the entry:
///   state[i] = |i0>_dq0 * |i1>_dq1 * |i2>_dq2 * ...
/// where ik is the k-th bit of i in its N-bit representation (with i0 being the
/// least significant bit) and dqk corresponds to the k-th data qubit.
///
/// The quantum algorithm is written in terms of gates (or other quantum operations)
/// acting on 'program qubits'. Instead, the way IQS represents a quantum register state
/// is based on 'data qubits', which are in 1:1 correspondence with program qubits
/// but may be in a different order. The order of data qubits determine which qubit
/// is 'local' and which is 'global' from the point of view of teh MPI communication.
///
/// When a QubitRegister is initialized, data qubits and program qubits correspond trivially:
///   pq(0) --> dq(0)      ,  pq(1) --> dq(1)      ,  ...
/// However this can be changed by using Permutations. Specifically one has:
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

/// Overload operator [] to address single entries of the permutation map.
///
/// Map from program qubits (used in algorithm) and data qubits (used in state representation):
///   program qubit k --> data qubit map[k]
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

/// Return a desription of the direct map as a string with space-separated entries.
  std::string GetMapStr()
  {
    std::string s;
    for (std::size_t i = 0; i < map.size(); i++)
        s += " " + qhipster::toString(map[i]);
    return s;
  }

/// Return a desription of the inverse map as a string with space-separated entries.
  std::string GetImapStr()
  {
    std::string s;
    for (std::size_t i = 0; i < imap.size(); i++)
        s += " " + qhipster::toString(imap[i]);
    return s;
  }

/////////////////////////////////////////////////////////////////////////////////////////

/// Create the identity Permutation.
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

/// Create and initialize the permutation from a vector.
  Permutation(std::vector<std::size_t> m)
  {
    num_qubits = m.size();
    SetNewPermutation(m);
  }

/////////////////////////////////////////////////////////////////////////////////////////

/// Find the data qubit associated with a specific program qubit.
///
/// If it exists, Find(position)=imap(position).
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

/// Set new permutation after validating the input vector as a valid permutation.
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

/// Utility functions to transform between decimal and binary.
///
///   decimal input = i0*2^0 + i1*2^1 + ...  ---->  "i0i1i2..."
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

/// Utility functions to transform between binary and decimal.
///
///   "i0i1i2..."  ---->  decimal integer = i0*2^0 + i1*2^1 + ...
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

  /// Transform index v in [0,2^N[ from the data representation to the program representation.
  inline std::size_t data2program_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < num_qubits; i++)
        v_ = v_ | (((v & (((std::size_t) 1) << map[i])) >> map[i]) << i);
    return v_;
  }

  /// Transform index v in [0,2^N[ from the program representation to the data representation.
  inline std::size_t program2data_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < num_qubits; i++)
        v_ = v_ | (((v & (((std::size_t)1) << imap[i])) >> imap[i]) << i);
    return v_;
  }

  /// Transform index v in [0,2^N[ from data decimal index to program string representation.
  std::string data2program(std::size_t v)
  {
    std::string s = dec2bin(v, num_qubits), sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[map[i]];

    if (0)
    {
        std::size_t v_ = 0;
        for (std::size_t i = 0; i < num_qubits; i++)
            v_ = v_ | (((v & (((std::size_t)1) << map[i])) >> map[i]) << i);
        printf("sp=%s new:%s\n", sp.c_str(), dec2bin(v_, num_qubits).c_str());
    }
    return sp;
  }

  /// Transform index v in [0,2^N[ from data string to program string representation.
  std::string data2program(std::string s)
  {
    std::string sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[map[i]];
    return sp;
  }

  /// Transform index v in [0,2^N[ from program decimal index to data string representation.
  std::string program2data(std::size_t v)
  {
    std::string s = dec2bin(v, num_qubits), sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[imap[i]];
    return sp;
  }

  /// Transform index v in [0,2^N[ from program string to data string representation.
  std::string program2data(std::string s)
  {
    std::string sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[imap[i]];
    return sp;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Print permutation.
  // Previously called 'prange()'
  void Print()
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
       printf("%s(%lld)\n", date2program(i).c_str(), bin2dec(data2program(i)));
    }

    printf("\n");
#endif

    for (std::size_t i = 0; i < (UL(1) << num_qubits); i++)
    {
#if 0
       printf("%s ==> %s\n", data2program(i).c_str(), 
               program2data(data2program(i)).c_str());
#else
      printf("map(%3lu) = %3lu\n", i, bin2dec(data2program(i)));
#endif
    }
  }
};

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard PERMUTATION_HPP
