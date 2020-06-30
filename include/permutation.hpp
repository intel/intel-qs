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
/// Define the order of N objects
/// (often, the N qubits in the representation of the quantum state).
///
/// For the identity permutation, elements and their position correspond trivially:
///   element  -->  position
///      0              0
///      1              1
///      2              2
///     ...            ...
///     N-1            N-1
///
/// However this can be changed by using a Permutation object. Specifically one has:
///   element  -->  position
///      0            map(0)
///      1            map(1)
///      2            map(2)
///     ...            ...
///     N-1           map(N-1)
///
/// and its inverse:
///   position -->  element 
///      0           imap(0)
///      1           imap(1)
///      2           imap(2)
///     ...            ...
///     N-1          imap(N-1)
///
/// @var map: map
/// @var imap: inverse map

class Permutation
{
 public:
  std::vector<std::size_t> map, imap;
  std::size_t num_elements;

/////////////////////////////////////////////////////////////////////////////////////////

/// Overload operator [] to address single entries of the permutation map.
///
/// Map from elements to their positions.
/// When elements are qubits, they represent program qubits while their positions correspond
/// to data qubits. One has the correspondence:
///         element k --> position map[k]
///   program qubit k --> data qubit map[k]

  // Appropriate only for rvalue.

  unsigned operator[](std::size_t i) const
  {
    assert(i < num_elements);
    return (unsigned) map[i];
  }

  unsigned operator[](unsigned i) const
  {
    assert(i < num_elements);
    return map[i];
  }

  int operator[](int i) const
  {
    assert(i < num_elements);
    return (int)map[i];
  }

/////////////////////////////////////////////////////////////////////////////////////////

/// Return the number of elements.
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
  Permutation(std::size_t num_elements)
  {
    this->num_elements = num_elements;
    std::vector<std::size_t> m(num_elements);
    // Initialize the trivial permutation: {0,1,2,... num_elements-1}
    iota(m.begin(), m.end(), 0);
    // for(auto i:m) printf("%d ", i); printf("\n");
    // exit(0);
    SetNewPermutationFromMap(m, "direct");
  }

/// Create and initialize the permutation from a vector.
  Permutation(std::vector<std::size_t> m, std::string style_of_map="direct")
  {
    num_elements = m.size();
    SetNewPermutationFromMap(m, style_of_map);
  }

/////////////////////////////////////////////////////////////////////////////////////////

/// Find the element associated with a specific position.
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

/// Set new permutation from direct or inverse map after validating the input as a valid permutation.
  void SetNewPermutationFromMap(std::vector<std::size_t> m, std::string style_of_map="direct")
  {
    // Check consistency of m.
    assert (m.size() == this->num_elements);
    std::vector<bool> exist(m.size(), 0);
    for (auto &j : m)
        exist[j] = 1;
    for (auto e : exist)
        assert(e > 0);

    // Compute direct map from inverse map.
    if (style_of_map=="direct")
    {
        map = m;
        imap.resize(num_elements);
        for (std::size_t q = 0; q < map.size(); q++)
            imap[map[q]] = q;
    }
    else if (style_of_map=="inverse")
    {
        imap = m;
        map.resize(num_elements);
        for (std::size_t pos = 0; pos < imap.size(); pos++)
            map[imap[pos]] = pos;
    }
    else assert(0);
  }

/////////////////////////////////////////////////////////////////////////////////////////

/// Exchange two elements in the permutation (i.e. in map).
///
/// Explicitly:
///    map[e1] = p1   ---->    map[e1] = p2
///    map[e2] = p2   ---->    map[e2] = p1
/// imap is updated to reflect the change in map.
  void ExchangeTwoElements(std::size_t element_1, std::size_t element_2)
  {
    std::size_t position_1 = map[element_1];
    std::size_t position_2 = map[element_2];
    map[element_1] = position_2;
    map[element_2] = position_1;
    imap[position_1] = element_2;
    imap[position_2] = element_1;
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
// Consider the case when the elements represent the bits of a base-2 number.
// Reorder the bits according to their position in the Permutation.
// What is the new number?
/////////////////////////////////////////////////////////////////////////////////////////
// Understanding bitwise operations:
//   v & (1 << map[i]))   ---->   zero bitstring, apart from the bit value corresponding
//                                to the position associated with element i
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
    for (std::size_t i = 0; i < num_elements; i++)
        v_ = v_ | (((v & (((std::size_t) 1) << map[i])) >> map[i]) << i);
    return v_;
  }

  /// Transform index v in [0,2^N[ from the program representation to the data representation.
  inline std::size_t program2data_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < num_elements; i++)
        v_ = v_ | (((v & (((std::size_t)1) << imap[i])) >> imap[i]) << i);
    return v_;
  }

  /// Transform index v in [0,2^N[ from data decimal index to program string representation.
  std::string data2program(std::size_t v)
  {
    std::string s = dec2bin(v, num_elements), sp(s);
    for (std::size_t i = 0; i < s.size(); i++)
        sp[i] = s[map[i]];

    if (0)
    {
        std::size_t v_ = 0;
        for (std::size_t i = 0; i < num_elements; i++)
            v_ = v_ | (((v & (((std::size_t)1) << map[i])) >> map[i]) << i);
        printf("sp=%s new:%s\n", sp.c_str(), dec2bin(v_, num_elements).c_str());
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
    std::string s = dec2bin(v, num_elements), sp(s);
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

  /// Print permutation effect on indices.
  void Print()
  {
      printf("qubit permutation: %s\n", this->GetMapStr().c_str());
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Print permutation effect on indices.
  // Previously called 'prange()'
  void PrintRange()
  {
#if 0
    printf("map:  ");
    for(auto & i : map) printf("%d", i);
    printf("\n");
    printf("imap: ");
    for(auto & i : imap) printf("%d", i);
    printf("\n");
    for(std::size_t i = 0; i < (1 << num_elements); i++)
    {
       printf("%s(%lld)\n", dec2bin(i, num_elements).c_str(), i);
    }

    printf("\n");
    for(std::size_t i = 0; i < (1 << num_elements); i++)
    {
       printf("%s(%lld)\n", date2program(i).c_str(), bin2dec(data2program(i)));
    }

    printf("\n");
#endif

    for (std::size_t i = 0; i < (UL(1) << num_elements); i++)
    {
#if 0
       printf("%s ==> %s\n", data2program(i).c_str(), 
               program2data(data2program(i)).c_str());
#else
      printf("map(%3lu) = %3lu\n", i, bin2dec(data2program(i)));
#endif
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////
;;
  /// Intermediate imaps dividing the elements in two groups: those with position < or >= M.
  ///
  /// Given the current Permutation and a new target map, find two intermediate inverse maps
  /// such that:
  /// * intermediate_imap_1 is obtained from old_imap by changing only the first M values
  /// * intermediate_imap_2 is obtained from intermediate_imap_1 by changing only the last
  ///   N-M values and can be transformed in new_imap by exchanging pairs (i,j) with i<M, j>=M.
  void ObtainIntemediateInverseMaps(std::vector<std::size_t> target_map,  size_t M,
                                    std::vector<std::size_t> & int_1_imap,
                                    std::vector<std::size_t> & int_2_imap)
  {
    // Check the validity of inputs.
    assert(M <= this->num_elements);
    assert (target_map.size() == this->num_elements);
    std::vector<bool> exist(target_map.size(), 0);
    for (auto &j : target_map)
        exist[j] = 1;
    for (auto e : exist)
        assert(e > 0);

    // Intermediate map updating the first M elements of the inverse map.
    std::size_t pos_try;
    int_1_imap = imap;
    for (std::size_t pos = 0; pos<M; pos++)
    {
        pos_try = target_map[ imap[pos] ];
        while (pos_try >= M)
            pos_try = target_map[ imap[pos_try] ];
        int_1_imap[pos_try] = imap[pos];
    }

    // Intermediate map updating the last num_elements-M elements of the inverse map.
    int_2_imap = int_1_imap;
    for (std::size_t pos = M; pos<this->num_elements; pos++)
    {
        pos_try = target_map[ int_1_imap[pos] ];
        while (pos_try < M)
            pos_try = target_map[ int_1_imap[pos_try] ];
        int_2_imap[pos_try] = int_1_imap[pos];
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////

};

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard PERMUTATION_HPP
