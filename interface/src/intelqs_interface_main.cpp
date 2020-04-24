#include <iostream>
#include <fstream>
#include <unordered_map>
#include <functional>
#include <stdexcept>

#include "../../include/qureg.hpp"
#include "../include/interface_api_qasm.h"

using namespace std;

// Global variables related to Psi-function .malloc/.free routines.
using Type = ComplexDP;
QubitRegister<Type> *psi1 = nullptr;
bool fPsiAllocated = false;

int main(int argc, char*argv[])
{
    qhipster::mpi::Environment env(argc, argv);
    string line = "";
    string token = "";

    if (env.IsUsefulRank() == false) 
        return 0;

    int myid = env.GetStateRank();
    using Type = ComplexDP;

    while(true) {
        getline(cin,line);

        if(line.length() >= 1) {
            int token_end = line.find_first_of(' ');
            unsigned long result = 1;

            token = line.substr(0,token_end);
            if(!token.empty()) {
	/*
               function<long(string)> func = qufun_table[token];
               if(func) {
                  result = func(line.substr(token_end+1,line.length()));
               }
*/
	       result = ExecuteHandler(token,line.substr(token_end+1,line.length()));

               if (result > 0) {
                   cerr << "Qasm Op failed - ["<<token<<"]"<<endl;
               }
            }
        } else
          break;
    }

return 0;
}
