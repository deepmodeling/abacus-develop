#include "module_io/output_parser.h"

namespace ModuleIO
{

Output_Parser::Output_Parser(Output_Interface& output) : _output(output){};

void Output_Parser::SetOutputFilename(std::string directory,
                                      std::string prefix,
                                      int is,
                                      std::string tag,
                                      std::string suffix,
                                      std::string delimiter)
{
    _output.setFilename(directory + prefix + delimiter + "SPIN" + std::to_string(is + 1) + delimiter + tag
                         + suffix);
}

void Output_Parser::SetPrecision(const int precision_in)
{
    _output.setPrecision(precision_in);
}

void Output_Parser::Write()
{
    _output.write();
}

} // namespace ModuleIO