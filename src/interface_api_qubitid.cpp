#include <unordered_map>
#include <string>

using namespace std;

/** Qubits are assigned a unique identifier that starts from this base value. */
#define QUBIT_ID_BASE 100

/**
 * @file interface_api_qubitid.cpp
 *
 * This file implements unique qubit id assignment and query operations.
 *
 * Each QASM instruction that comes into the qHiPSTER QASM interface has a
 * unique identifier string for the operand:
 *
 * 	H    qubit0
 * 	CNOT qubit0,qubit1
 * 	...
 * 
 * qHiPSTER applies qubit operations to a register of qubits that are indexed
 * from [0...(N-1)] where N is the total number of qubits. For the QASM snippet
 * above, the qubits will be assigned a unique id in the order they are
 * encountered in the QASM. To explicitly assign qubit identifiers to qHiPSTER
 * register locations you must add the (qubit) instruction to your QASM in the
 * right order as follows:
 *
 * 	qubit qubit0
 * 	qubit qubit1
 * 	H qubit0
 * 	CNOT qubit0,qubit1
 * 	...
 * 
 * This ensures that qubit0 -> qreg[0], qubit1 -> qureg[1] in the qHiPSTER
 * quantum register space.
 *
 **/

// Hash table mapping unique qubit identifiers to each qubit operand.
unordered_map<string,int> qubit_id_table;

// Variable to keep track of the next qubit to assign.
int next_qubit_id = QUBIT_ID_BASE;


// See the header file (interface_api_qubitid.h) for the description of the function.
int query_qubit_id(string qubit_name) {

   // Retrieve the qubit_id from the qubit_id_table.
   int& qubit_ref = qubit_id_table[qubit_name];

   // Insert a new qubit_id in the Hash table if we did not get back a valid one.
   if(qubit_ref < QUBIT_ID_BASE)
   {
      qubit_ref = next_qubit_id;
      next_qubit_id++;
   }

   // Translate the Hash encoded ID to a qHiPSTER API relevant ID and return it.
   //cout << "Assigned id: "<<qubit_ref<<endl;
   return (qubit_ref - QUBIT_ID_BASE); 
}

