#include "source_io/parse_command_line.h"

#include <cctype>
#include <cstdio>
#include <stdexcept>

namespace ModuleIO
{

int run_control_option_size(int argc, char* argv[], int index)
{
    const std::string option = argv[index];
    int size = 0;
    std::string required;
    if (option == "-p" || option == "--parameter")
    {
        size = 3;
        required = "<name> <value>";
    }
    else if (option == "-in" || option == "--input")
    {
        size = 2;
        required = "<file>";
    }
    else if (option == "--check-input")
    {
        size = 1;
    }
    if (size > argc - index)
    {
        throw std::runtime_error("Option " + option + " requires " + required);
    }
    return size;
}

CommandLineArgs parse_command_line(int argc, char* argv[])
{
    CommandLineArgs args;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        const int option_size = run_control_option_size(argc, argv, i);

        if (option_size == 3)
        {
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
        else if (option_size == 2)
        {
            ++i;
            args.input_file = argv[i];
            if (args.input_file.empty())
            {
                throw std::runtime_error("Option " + arg + " requires a nonempty <file>");
            }
        }
        else if (option_size == 1)
        {
            // parse_args has already enabled ReadInput's check mode.
            continue;
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
    std::fprintf(stderr, "Run '%s -h' for usage.\n", bin.c_str());
}

} // namespace ModuleIO

