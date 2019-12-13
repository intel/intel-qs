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

#ifndef IQS_GATE_COUNTER_HPP
#define IQS_GATE_COUNTER_HPP

#include <vector>

/////////////////////////////////////////////////////////////////////////////////////////

/// @file gate_counter.hpp

/////////////////////////////////////////////////////////////////////////////////////////

/// \class GateCounter
/// The GateCounter class serves two main purposes:
/// 1) To count the number of gates applied, divided by kind.
/// 2) To estimate the circuit depth if scheduled in a greedy way.
class GateCounter
{
 private:
  int num_qubits;
  int total_gate_count=0;
  int one_qubit_gate_count=0;
  int two_qubit_gate_count=0;
  std::vector<int> parallel_depth;

 public:

/////////////////////////////////////////////////////////////////////////////////////////

  // Constructor
  GateCounter (int new_num_qubits)
    : num_qubits(new_num_qubits)
  {
    parallel_depth.assign(num_qubits,0);
  }

  // Destructor
  ~GateCounter()
  { }

/////////////////////////////////////////////////////////////////////////////////////////

  void Reset()
  {
    parallel_depth.assign(num_qubits,0);
  }

/////////////////////////////////////////////////////////////////////////////////////////

  //Get stats
  int GetTotalGateCount()    { return total_gate_count; }
  int GetOneQubitGateCount() { return one_qubit_gate_count; }
  int GetTwoQubitGateCount() { return two_qubit_gate_count; }
  int GetParallelDepth()
  { return *std::max_element(std::begin(parallel_depth), std::end(parallel_depth)); }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Update the counters and depth due to the action of a one-qubit gate.
  void OneQubitIncrement (int qubit)
  {
    ++total_gate_count;
    ++one_qubit_gate_count;
    ++parallel_depth[qubit];
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Update the counters and depth due to the action of a two-qubit gate.
  void TwoQubitIncrement (int qubit_0, int qubit_1)
  {
    ++total_gate_count;
    ++two_qubit_gate_count;
    int new_depth = std::max(parallel_depth[qubit_0], parallel_depth[qubit_1]) +1;
    parallel_depth[qubit_0] = new_depth;
    parallel_depth[qubit_1] = new_depth;
  }

/////////////////////////////////////////////////////////////////////////////////////////

  /// Print the values of counters and depth.
  void Breakdown()
  {
    if (qhipster::mpi::Environment::GetStateRank() == 0)
    {
        printf("The quantum circuit is composed of %d one-qubit gates and "
               "%d two-qubitgates, for a total of %d gates.\nThe greedy depth "
               "(all gates lasting one clock cycle) is %d.\n",
               GetOneQubitGateCount(), GetTwoQubitGateCount(),
               GetTotalGateCount(), GetParallelDepth() );
    }
  }

};

/////////////////////////////////////////////////////////////////////////////////////////

#endif	// header guard IQS_GATE_COUNTER_HPP
