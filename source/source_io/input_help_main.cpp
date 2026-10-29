#include "input_help.h"

#include <iostream>

int main()
{
    ModuleIO::ParameterHelp::generate_yaml(std::cout);
    return 0;
}
