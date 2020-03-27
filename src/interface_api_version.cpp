#include <iostream>

#include "../include/qureg_version.hpp"
#include "../include/interface_api_version.h"


using namespace std;


unsigned long quiversion(string args) {

    cout << INTERFACE_VERSION_STRING << endl;

    return 0;
}


unsigned long quversion(string args) {

    cout << GetQhipsterVersion() << endl;

    return 0;
}
