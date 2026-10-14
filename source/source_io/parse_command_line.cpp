#include "source_io/parse_command_line.h"

#include <cctype>
#include <cstdio>
#include <stdexcept>

namespace ModuleIO
{

CommandLineArgs parse_command_line(int argc, char* argv[])
{
    CommandLineArgs args;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];

        if (arg == "-p" || arg == "--parameter")
        {
            if (i + 2 >= argc)
            {
                const std::string msg = "Option " + arg
                                        + " requires <name> <value>";
                throw std::runtime_error(msg);
            }
            ++i;
            const std::string name = argv[i];
            if (name.empty())
            {
                const std::string msg = "Invalid variable name: (empty)";
                throw std::runtime_error(msg);
            }
            const unsigned char first_char
                = static_cast<unsigned char>(name[0]);
            const bool starts_with_digit = (std::isdigit(first_char) != 0);
            if (starts_with_digit)
            {
                const std::string msg = "Invalid variable name: " + name;
                throw std::runtime_error(msg);
            }
            ++i;
            const std::string value = argv[i];
            args.vars[name] = value;   // later duplicates override earlier ones
        }
        else if (arg == "-in" || arg == "--input")
        {
            if (i + 1 >= argc)
            {
                const std::string msg = "Option " + arg + " requires <file>";
                throw std::runtime_error(msg);
            }
            ++i;
            args.input_file = argv[i];
        }
        else
        {
            // Unreachable in normal flow: parse_args rejects unknown
            // flags before this function runs. Kept as a defensive check
            // so parse_command_line stays self-contained for tests.
            const std::string msg = "Unknown option: " + arg;
            throw std::runtime_error(msg);
        }
    }
    return args;
}

void print_help(const std::string& bin)
{
    // Run-control options only; the general help is owned by
    // ParameterHelp::show_general_help (via parse_args).
    std::printf(
        "Run-control options of %s:\n"
        "  -p, --parameter <name> <value>\n"
        "                         Set an INPUT variable from the command line,\n"
        "                         usable in INPUT as ${name} or $name.\n"
        "                         Repeatable; overrides INPUT 'variable'.\n"
        "  -in, --input <file>    Path to INPUT (default: INPUT)\n",
        bin.c_str());
}

} // namespace ModuleIO

