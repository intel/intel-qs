// Copyright (C) 2015 Theoretical Physics, ETHZ Zurich

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "random_engine.hpp"

#include <iostream>
#include <map>
#include <random>

template <class RNG>
void print_numbers(RNG& engine, unsigned n)
{
  for (unsigned i = 0; i < n; ++i) std::cout << engine() << std::endl;
}

template <class RNG>
void distribution(RNG& engine)
{
  // perform 4 trials, each succeeds 1 in 2 times
  std::binomial_distribution<> d(4, 0.5);

  std::map<int, int> hist;
  for (int n = 0; n < 10000; ++n) {
    ++hist[d(engine)];
  }
  for (auto p : hist) {
    std::cout << p.first << ' ' << std::string(p.second / 100, '*') << '\n';
  }
}

int main()
{

  {  /// Call API
    openqu::random::engine<> engine("http://random.openqu.org/api");

    std::cout << "+++ GETTING 300 QRNG +++" << std::endl;
    print_numbers(engine, 300);

    std::cout << "+++ NORMAL DISTRIBUTION +++" << std::endl;
    distribution(engine);
  }
  {  /// Call API with different backend
    openqu::random::engine<std::mt19937> engine("http://random.openqu.org/api");

    std::cout << "+++ GETTING 10 QRNG +++" << std::endl;
    print_numbers(engine, 10);

    std::cout << "+++ NORMAL DISTRIBUTION +++" << std::endl;
    distribution(engine);
  }
  {  /// Wrong API url
    openqu::random::engine<> engine("http://wrong__url__/api");

    std::cout << "+++ GETTING 300 QRNG +++" << std::endl;
    print_numbers(engine, 300);
  }
}