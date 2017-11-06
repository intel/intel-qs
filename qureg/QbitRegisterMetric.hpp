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
#include <vector>
using namespace std;

template <class Type = ComplexDP>
class QbitRegisterMetric: public QbitRegister<Type> {
  int iTotalQubitGateCount=0;
  int iOneQubitGateCount=0;
  int iTwoQubitGateCount=0;
  std::vector<int> vParallelDepth;
  void OneQubitIncrements(int);
  void TwoQubitIncrements(int,int);

public:
  //Constructor
  QbitRegisterMetric<Type>(int iNQbits):QbitRegister<Type>(iNQbits){
    vParallelDepth.resize(iNQbits);
  }

  //Get stats
  int GetTotalQubitGateCount();
  int GetOneQubitGateCount();
  int GetTwoQubitGateCount();
  int GetParallelDepth();

  //Perform gates
  void applyHadamard(int);
  void applyRotationX(int, double);
  void applyRotationY(int, double);
  void applyRotationZ(int, double);
  void applyCPauliX(int, int);
  void applyControlled1QubitGate(int, int, openqu::TinyMatrix<Type, 2, 2, 32>);
};

template <class Type>
int QbitRegisterMetric<Type>::GetOneQubitGateCount(){
  return iOneQubitGateCount;
}

template <class Type>
int QbitRegisterMetric<Type>::GetTwoQubitGateCount(){
  return iTwoQubitGateCount;
}

template <class Type>
int QbitRegisterMetric<Type>::GetTotalQubitGateCount(){
  return iTotalQubitGateCount;
}

template <class Type>

int QbitRegisterMetric<Type>::GetParallelDepth(){
  return *std::max_element(std::begin(vParallelDepth), std::end(vParallelDepth));
}

template <class Type>
void QbitRegisterMetric<Type>::OneQubitIncrements(int q){
  iTotalQubitGateCount++;
  iOneQubitGateCount++;
  vParallelDepth[q]++;
}

template <class Type>
void QbitRegisterMetric<Type>::TwoQubitIncrements(int q1, int q2){
  iTotalQubitGateCount++;
  iTwoQubitGateCount++;
  int iNewDepth = max(vParallelDepth[q1],vParallelDepth[q2])+1;
  vParallelDepth[q1]=iNewDepth;
  vParallelDepth[q2]=iNewDepth;
}

template <class Type>
void QbitRegisterMetric<Type>::applyHadamard(int q){
  QbitRegister<Type>::applyHadamard(q);
  OneQubitIncrements(q); 
}

template <class Type>
void QbitRegisterMetric<Type>::applyRotationX(int q, double theta){
  QbitRegister<Type>::applyRotationX(q,theta);
  OneQubitIncrements(q); 
}

template <class Type>
void QbitRegisterMetric<Type>::applyRotationY(int q, double theta){
  QbitRegister<Type>::applyRotationY(q,theta);
  OneQubitIncrements(q); 
}

template <class Type>
void QbitRegisterMetric<Type>::applyRotationZ(int q, double theta){
  QbitRegister<Type>::applyRotationZ(q,theta);
  OneQubitIncrements(q); 
}

template <class Type>
void QbitRegisterMetric<Type>::applyCPauliX(int q1, int q2){
  QbitRegister<Type>::applyCPauliX(q1,q2);
  TwoQubitIncrements(q1,q2); 
}

template <class Type>
void QbitRegisterMetric<Type>::applyControlled1QubitGate(int q1, int q2, openqu::TinyMatrix<Type, 2, 2, 32> V){
  QbitRegister<Type>::applyControlled1QubitGate(q1,q2,V);
  TwoQubitIncrements(q1,q2); 
}
