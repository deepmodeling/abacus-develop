#ifndef READ_INP_SCCS_H
#define READ_INP_SCCS_H
#include <string>
struct Input_para;
namespace ModuleIO
{
bool parse_solvation_model(const std::string& value, int& model, std::string& error);
bool validate_sccs_input(const Input_para& input, std::string& error);
}
#endif
