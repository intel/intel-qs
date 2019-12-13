// Copyright (C) 2015 Theoretical Physics, ETH Zurich

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

#include <iostream>
#include <vector>

extern "C" {
double ddot_(const int*, const double*, const int*, const double*, const int*);
}

int main()
{
  std::vector<double> values(2, 1.);

  int N = 2;
  int one = 1;

//  double norm = ddot_(&N, &values[0], &one, &values[0], &one);
//  std::cout << "Norm: " << norm << std::endl;
  return 0;
}
