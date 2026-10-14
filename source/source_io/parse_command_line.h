#ifndef ModuleIO_PARSE_COMMAND_LINE_H
#define ModuleIO_PARSE_COMMAND_LINE_H

#include <map>
#include <string>

namespace ModuleIO
{

/// Run-control options parsed after ModuleIO::parse_args.
///
/// Flag ownership:
///   - Informational flags (-v/-i/-h/-s/--generate-parameters-yaml) and
///     --check-input are owned by ModuleIO::parse_args and never reach
///     this function.
///   - -p/--parameter and -in/--input are validated (token count only)
///     and skipped by parse_args, then fully parsed here.
struct CommandLineArgs
{
    // -p / --parameter <name> <value> : INPUT variable injection (repeatable)
    // These override INPUT 'variable' definitions.
    std::map<std::string, std::string> vars;

    // -in / --input <file> : explicit INPUT path, default "INPUT"
    std::string input_file = "INPUT";
};

/// Parses -p/--parameter and -in/--input options.
/// Throws std::runtime_error on malformed arguments (parse_args has
/// already filtered unknown flags before this function is called).
CommandLineArgs parse_command_line(int argc, char* argv[]);

/// Prints usage for the run-control options handled here.
void print_help(const std::string& bin_name);

} // namespace ModuleIO

#endif // ModuleIO_PARSE_COMMAND_LINE_H

