#ifndef WRITE_OUTPUT_H
#define WRITE_OUTPUT_H

#include <string>

#include "module_io/output_interface.h"

namespace ModuleIO
{
class Output_Parser
{
  public:
    Output_Parser(Output_Interface& output);
    void SetOutputFilename(std::string directory,
                           std::string prefix,
                           int is,
                           std::string tag,
                           std::string suffix,
                           std::string delimiter);
    void SetPrecision(const int precision_in);
    void Write();

  private:
    Output_Interface& _output;
};
} // namespace ModuleIO

#endif