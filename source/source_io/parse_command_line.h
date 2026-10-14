#ifndef ModuleIO_PARSE_COMMAND_LINE_H
#define ModuleIO_PARSE_COMMAND_LINE_H

#include <map>
#include <string>

namespace ModuleIO
{

/// Run-control options parsed after ModuleIO::parse_args.
///
/// Flag ownership:
///   - Informational flags (-v/-i/-h/-s/--generate-parameters-yaml)
///     are owned by ModuleIO::parse_args and exit before this function.
///   - --check-input sets check mode in parse_args and is skipped here.
///   - -p/--parameter and -in/--input share token-count validation with
///     parse_args, then their values are parsed here, before MPI setup.
struct CommandLineArgs
{
    // -p / --parameter <name> <value>: stored variables (repeatable).
    // PART I does not apply these values to INPUT; main emits a warning.
    std::map<std::string, std::string> vars;

    // -in / --input <file> : explicit INPUT path, default "INPUT"
    std::string input_file = "INPUT";
};

/// Parses -p/--parameter and -in/--input options.
/// Throws std::runtime_error on malformed arguments (parse_args has
/// already filtered unknown flags before this function is called).
CommandLineArgs parse_command_line(int argc, char* argv[]);

/// Number of tokens occupied by one run-control option, including its values.
/// Returns zero for other options; throws std::runtime_error for missing values.
int run_control_option_size(int argc, char* argv[], int index);

/// Prints a short pointer to the canonical general help on standard error.
void print_help(const std::string& bin_name);

} // namespace ModuleIO

#endif // ModuleIO_PARSE_COMMAND_LINE_H

