//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------
#pragma once

#include "util/utils.hpp"
#include <numeric>
#include <map>
#include <vector>
#include "util/conversion.hpp"

inline std::size_t perm(std::size_t v, std::size_t *map, std::size_t nqbits)
{
  std::size_t v_ = 0;
  for (std::size_t i = 0; i < nqbits; i++) v_ = v_ | (((v & (1 << map[i])) >> map[i]) << i);
  return v_;
}

class Permutation
{
 public:
  std::vector<std::size_t> map, imap;
  std::size_t nqbits;

  unsigned operator[](std::size_t i)
  {
    assert(i <= nqbits);
    return (unsigned) map[i];
  }

  unsigned operator[](unsigned i)
  {
    assert(i <= nqbits);
    return map[i];
  }

  std::size_t size() {return map.size();}

  std::string GetMapStr()
  {
    std::string s;
    for (std::size_t i = 0; i < map.size(); i++) s += " " + openqu::toString(map[i]);
    return s;
  }

  std::string GetImapStr()
  {
    std::string s;
    for (std::size_t i = 0; i < imap.size(); i++) s += " " + openqu::toString(imap[i]);
    return s;
  }

  Permutation(std::size_t nqbits)
  {
    this->nqbits = nqbits;
    std::vector<std::size_t> m(nqbits);
    iota(m.begin(), m.end(), 0);
    // for(auto i:m) printf("%d ", i); printf("\n");
    // exit(0);
    SetNewPerm(m);
  }

  Permutation(std::vector<std::size_t> m)
  {
    nqbits = m.size();
    SetNewPerm(m);
  }

  std::size_t Find(std::size_t pos)
  {
    bool found = false;
    std::size_t i;
    for (i = 0; i < map.size(); i++)
      if (map[i] == pos) {
        found = true;
        break;
      }

    assert(found);
    return i;
  }

  void SetNewPerm(std::vector<std::size_t> m)
  {

    map = m;
    // check consistency of map
    std::vector<bool> exist(map.size(), 0);
    for (auto &m : map) exist[m] = 1;
    for (auto e : exist) assert(e > 0);

    // compute inverse map
    imap = map;
    for (std::size_t i = 0; i < map.size(); i++) imap[map[i]] = i;
  }

  std::string dec2bin(std::size_t in, std::size_t nbits)
  {
    std::string s;
    s.resize(nbits);
    for (std::size_t i = 0; i < nbits; i++) {
      s[i] = (in & 0x1) == 0 ? '0' : '1';
      in = in >> 1;
    }
    return s;
  }
  std::size_t bin2dec(std::string in)
  {
    std::size_t v = 0;
    for (std::size_t i = 0; i < in.size(); i++) {
      v += UL(1 << i) * UL(in[i] - '0');
    }
    return v;
  }

  inline std::size_t lin2perm_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < nqbits; i++) v_ = v_ | (((v & (1 << map[i])) >> map[i]) << i);
    return v_;
  }
  inline std::size_t perm2lin_(std::size_t v)
  {
    std::size_t v_ = 0;
    for (std::size_t i = 0; i < nqbits; i++)
      v_ = v_ | (((v & (((std::size_t)1) << imap[i])) >> imap[i]) << i);
    return v_;
  }

  std::string lin2perm(std::size_t v)
  {
    std::string s = dec2bin(v, nqbits), sp(s);
    for (std::size_t i = 0; i < s.size(); i++) sp[i] = s[map[i]];
    if (0) {
      std::size_t v_ = 0;
      for (std::size_t i = 0; i < nqbits; i++)
        v_ = v_ | (((v & (((std::size_t)1) << map[i])) >> map[i]) << i);
      printf("sp=%s new:%s\n", sp.c_str(), dec2bin(v_, nqbits).c_str());
    }
    return sp;
  }

  std::string lin2perm(std::string s)
  {
    std::string sp(s);
    for (std::size_t i = 0; i < s.size(); i++) sp[i] = s[map[i]];
    return sp;
  }

  std::string perm2lin(std::size_t v)
  {
    std::string s = dec2bin(v, nqbits), sp(s);
    for (std::size_t i = 0; i < s.size(); i++) sp[i] = s[imap[i]];
    return sp;
  }
  std::string perm2lin(std::string s)
  {
    std::string sp(s);
    for (std::size_t i = 0; i < s.size(); i++) sp[i] = s[imap[i]];
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
    for(std::size_t i = 0; i < (1 << nqbits); i++)
    {
       printf("%s(%lld)\n", dec2bin(i, nqbits).c_str(), i);
    }

    printf("\n");
    for(std::size_t i = 0; i < (1 << nqbits); i++)
    {
       printf("%s(%lld)\n", lin2perm(i).c_str(), bin2dec(lin2perm(i)));
    }

    printf("\n");
#endif

    for (std::size_t i = 0; i < (UL(1) << nqbits); i++) {
#if 0
       printf("%s ==> %s\n", lin2perm(i).c_str(), 
               perm2lin(lin2perm(i)).c_str());
#else
      printf("map(%3lu) = %3lu\n", i, bin2dec(lin2perm(i)));
#endif
    }
  }
};

#if defined(MAIN)
class State
{
 public:
  std::string name;
  Permutation p;
  std::vector<ComplexType> state;
  State(Permutation p_, std::string name_) : p(p_), name(name_)
  {
    state.resize(1 << p.nqbits);
#pragma omp parallel for
    for (std::size_t i = 0; i < state.size(); i++) state[i] = {D(i % 3), D(i % 7)};
  }
  void permute(Permutation pnew)
  {
    Permutation pold = p;

    assert(pnew.nqbits == pold.nqbits);
    std::vector<ComplexType> state_new(state.size(), 0);

    // printf("map: %s imap: %s\n", pnew.GetMapStr().c_str(), pold.GetImapStr().c_str());
    std::vector<std::size_t> map(pnew.nqbits, 0);
    for (std::size_t i = 0; i < pnew.nqbits; i++) {
      map[i] = pold.map[pnew.imap[i]];
      // printf("%d ", map[i]);
    }
    // printf("\n");

    __int64 t0 = __rdtsc();
    double s0 = sec();
#pragma omp parallel for
    for (std::size_t i = 0; i < state.size(); i++) {
#if 0
      std::size_t to1 = perm(i, &(map[0]), p.nqbits);
      std::size_t to2 = pnew.perm2lin_(pold.lin2perm_(i));
      std::size_t to = pold.bin2dec(pnew.perm2lin(pold.lin2perm(i)));
      assert(to == to1);
      assert(to == to2);
      state_new[to] = state[i];
#else
      std::size_t to_ = perm(i, &(map[0]), p.nqbits);  // pnew.perm2lin_(pold.lin2perm_(i));
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

  void print()
  {
    printf("name::%s %s\n", name.c_str(), p.GetMapStr().c_str());
    for (std::size_t i = 0; i < state.size(); i++) {
      printf("%s {%lf %lf}\n", p.lin2perm(i).c_str(), real(state[i]), imag(state[i]));
    }
  }
  bool operator==(const State &rhs)
  {
    assert((const std::size_t)state.size() == rhs.state.size());
    for (std::size_t i = 0; i < state.size(); i++) {
      if (state[i] != rhs.state[i]) return false;
    }

    return true;
  }
};

std::size_t main(std::size_t argc, char **argv)
{
  std::size_t nqbits = 3, nthreads = 1;
  if (argc != 3) {
    fprintf(stderr, "usage: %s <nqbits> <nthreads>\n", argv[0]);
    exit(1);
  } else {
    nqbits = atoi(argv[1]);
    nthreads = atoi(argv[2]);
  }
  initomp(nthreads);
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
  std::vector<std::size_t> map(nqbits, 0);
  iota(map.begin(), map.end(), 0);
  State s1(Permutation(map), "s1"), s2(s1);
  s2.permute(Permutation(map));
  s2.permute(Permutation(map));
  s2.permute(Permutation(map));
#endif
#endif
}
#endif
