//------------------------------------------------------------------------------
// Copyright (C) 2017 Intel Corporation 
//
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
//------------------------------------------------------------------------------

#include "../qureg/qureg.hpp"

using Type = ComplexDP;
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
template<typename Type>
void GenerateGateSet(std::vector <std::pair<std::string,TM2x2<Type>>> &sqg,
                     std::vector <std::pair<std::string,TM2x2<Type>>> &cqg)
{
  TM2x2<Type> eye, X, sqrtX, Y, sqrtY, Z, sqrtZ, H, T;

  sqg.resize(0);
  cqg.resize(0);

  // create vector of single qubit gates
  eye(0,0) = {1., 0.};
  eye(0,1) = {0., 0.};
  eye(1,0) = {0., 0.};
  eye(1,1) = {1., 0.};
  eye.name = "i";

  X(0, 0) = Type(0., 0.); X(0, 1) = Type(1., 0.);
  X(1, 0) = Type(1., 0.); X(1, 1) = Type(0., 0.);
  X.name = "x";

  sqrtX(0, 0) = Type(0.5,  0.5); sqrtX(0, 1) = Type(0.5, -0.5);
  sqrtX(1, 0) = Type(0.5, -0.5); sqrtX(1, 1) = Type(0.5,  0.5);
  sqrtX.name = "x_1_2";

  Y(0, 0) = Type(0., 0.); Y(0, 1) = Type(0., -1.);
  Y(1, 0) = Type(0., 1.); Y(1, 1) = Type(0., 0.);
  Y.name = "y";

  sqrtY(0, 0) = Type(0.5,   0.5); sqrtY(0, 1) = Type(-0.5, -0.5);
  sqrtY(1, 0) = Type(0.5,   0.5); sqrtY(1, 1) = Type(0.5,  0.5);
  sqrtY.name = "y_1_2";

  Z(0, 0) = Type(1., 0.); Z(0, 1) = Type(0., 0.);
  Z(1, 0) = Type(0., 0.); Z(1, 1) = Type(-1., 0.);
  Z.name = "z";

  sqrtZ(0, 0) = Type(1., 0.); sqrtZ(0, 1) = Type(0., 0.);
  sqrtZ(1, 0) = Type(0., 0.); sqrtZ(1, 1) = Type(0., 1.);
  sqrtZ.name = "z_1_2";

  T(0, 0) = Type(1.0, 0.0); T(0, 1) = Type(0.0, 0.0);
  T(1, 0) = Type(0.0, 0.0); T(1, 1) = Type(cos(M_PI/4.0), sin(M_PI/4.0));
  T.name = "t";

  double f = 1. / std::sqrt(2.);
  H(0, 0) = H(0, 1) = H(1, 0) = Type(f, 0.);
  H(1, 1) = Type(-f, 0.);
  H.name = "h";

  sqg = {{" h ",     H    }, {" x_1_2 ", sqrtX}, {" y_1_2 ", sqrtY}, {" z_1_2 ", sqrtZ}, {" t ",     T    }};
  cqg = {{" cz ", Z},        {" ch ",    H    }};

}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  unsigned myrank;
#ifdef INTELQS_HAS_MPI
  openqu::mpi::Environment env(argc, argv);
  if (env.is_usefull_rank() == false) return 0;
  myrank = env.rank();
#endif
  
  unsigned num_qubits = 3, num_threads = 1;
  std::size_t tmp_size = 0;
  if(argc != 2)
  {
     fprintf(stderr, "usage: %s <num_qubits>\n", argv[0]);
     exit(1);
  }
  else
  {
     num_qubits = atoi(argv[1]);
  }


  std::vector <std::pair<std::string,TM2x2<Type>>> sqg, cqg;
  GenerateGateSet<Type>(sqg, cqg);

  QubitRegister<ComplexDP> psi1(num_qubits, "rand", 0, 2097152);
  QubitRegister<ComplexDP> psi2(num_qubits, "rand", 0, 2097152);

  psi2.TurnOnSpecialize();
  psi2.EnableStatistics();

  for(auto g: sqg)
  {
    for(int pos = 16; pos < num_qubits; pos++)
    {
       psi1.Apply1QubitGate(pos, g.second);
       psi2.Apply1QubitGate(pos, g.second);
       assert(psi1 == psi2);
    }
  }

  for (auto &g : cqg)
  {
    TM2x2<Type> &m = g.second;
    for (unsigned q1 = 0; q1 < num_qubits; q1++)
    {
      for (unsigned q2 = 0; q2 < num_qubits; q2++)
      {
        if (q1 != q2)
        {
          psi1.ApplyControlled1QubitGate(q1, q2, m);
          psi2.ApplyControlled1QubitGate(q1, q2, m);
          assert(psi1 == psi2);
        }
      }
    }
  }

  psi2.GetStatistics();
  
}
