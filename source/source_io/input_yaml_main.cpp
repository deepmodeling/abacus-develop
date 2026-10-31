// Standalone dumper of the INPUT parameter reference.
//
// Prints the same YAML as `abacus --generate-parameters-yaml`, but links only
// the Input_Item registry instead of the whole program, so the documentation
// build (.readthedocs.yaml) can produce docs/parameters.yaml without compiling
// every solver.
#include "input_help.h"

#include <iostream>

int main()
{
    ModuleIO::ParameterHelp::generate_yaml(std::cout);
    return 0;
}
