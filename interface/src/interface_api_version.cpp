#include <iostream>

#include "../include/qureg_version.hpp"
#include "../include/interface_api_version.h"

unsigned long quiversion(std::string args) {

    std::cout << INTERFACE_VERSION_STRING << std::endl;

    return 0;
}


unsigned long quversion(std::string args) {

    std::cout << iqs::GetQhipsterVersion() << std::endl;

    return 0;
}
