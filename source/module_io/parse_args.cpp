#include "parse_args.h"

#include <cstdlib>
#include <iostream>

#include "module_io/read_input.h"
#include "version.h"

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
            ModuleIO::ReadInput::check_mode = true;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::exit(1);
        }
    }
}

} // namespace ModuleIO