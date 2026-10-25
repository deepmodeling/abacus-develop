#ifndef READ_INP_SCCS_H
#define READ_INP_SCCS_H
#include <string>
struct Input_para;
namespace ModuleIO
{
bool parse_imp_sol(const std::string& value, int& model, std::string& error);
bool check_sccs_inp(const Input_para& input, std::string& error);
}
#endif
