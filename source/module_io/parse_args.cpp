#include "parse_args.h"
#include <iostream>
#include <cstdlib>
#include "version.h"
#include "module_io/input.h"

namespace ModuleIO
{

void parse_args(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i) // Start from 1 to skip the program name
    {
        std::string arg = argv[i];
        if (arg == "--version" || arg == "-v" || arg == "-V")
        {
#ifdef VERSION
            const char* version = VERSION;
#else
            const char* version = "unknown";
#endif
            std::cout << "ABACUS version " << version << std::endl;
            std::exit(0);
        }
        else if (arg == "--check-input")
        {
            INPUT.check_input = true;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::exit(1);
        }
    }
}

}