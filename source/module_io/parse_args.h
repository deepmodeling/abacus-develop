#ifndef ModuleIO_PARSE_ARGS_H
#define ModuleIO_PARSE_ARGS_H
#include <iostream>
#include <cstdlib>
#include "version.h"

namespace ModuleIO
{
void parse_args(int argc, char** argv)
{
    if (argc > 1
        && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-v" || std::string(argv[1]) == "-V"))
    {
#ifdef VERSION
        const char* version = VERSION;
#else
        const char* version = "unknown";
#endif
        std::cout << "ABACUS version " << version << std::endl;
        std::exit(0);
    }
    
    return;
}
} // namespace ModuleIO

#endif