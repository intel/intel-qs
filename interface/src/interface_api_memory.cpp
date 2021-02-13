#include <iostream>
#include <stdexcept>

#include "../../include/qureg.hpp"

using namespace std;


// Global variables related to Psi-function .malloc/.free routines.
using Type = ComplexDP;
extern iqs::QubitRegister<Type> *psi1;
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
        psi1 = new iqs::QubitRegister<Type>(num_qubits);

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
