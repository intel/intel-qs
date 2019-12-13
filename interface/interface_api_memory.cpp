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
#include <iostream>
#include <stdexcept>

#include "../qureg/qureg.hpp"

using namespace std;


// Global variables related to Psi-function .malloc/.free routines.
using Type = ComplexDP;
extern QubitRegister<Type> *psi1;
extern bool fPsiAllocated;


unsigned long qumalloc(string args) {

    int num_qubits = 0;

    // Check for an attempt at a memory leak.
    if (fPsiAllocated) {
        cerr << ".malloc called before .free - memory leak detected."<<endl;
        return 1;
    }

    // Convert the argument to an integer number of qubits.
    try
    {
        num_qubits = stoi(args,nullptr,10);

    } catch(const invalid_argument& inv_arg) {
        cerr << ".malloc (size) - size parameter was not an integer."<<endl;
    } catch(const out_of_range& oor_arg) {
        cerr << ".malloc (size) - size parameter was too large to fit in an integer."<<endl;
    }

    // Ensure wavefunction register is in the allowed range.
    if ((num_qubits > 0) && (num_qubits <= 43)) {
        psi1 = new QubitRegister<Type>(num_qubits);

        if (psi1) {
            (*psi1)[0] = 1;
            fPsiAllocated = true;
            cout << "Allocated ["<<num_qubits<<"] qubits."<<endl;
            return 0;
        } 
    }
    else 
    {
        cerr << ".malloc (size) - allocations in the range from 1..43 only."<<endl;
    }
 

    return 1;
}


unsigned long qufree(string args) {

    if (fPsiAllocated) {
        delete psi1;
        psi1 = nullptr;
        fPsiAllocated = false;

        cout << "Qubits successfully free."<<endl;
        return 0;
    }

    cerr << "Harmless attempt to .free before .malloc. Nothing done."<<endl;

    return 1;
}
