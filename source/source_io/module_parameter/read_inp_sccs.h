#ifndef READ_INP_SCCS_H
#define READ_INP_SCCS_H
#include <string>
#include <vector>
struct Input_para;
namespace ModuleIO
{
bool parse_imp_sol(const std::string& value, int& model, std::string& error);
bool check_sccs_inp(const Input_para& input, std::string& error);
// Parse sccs_corespread words as finite real numbers; on failure spreads is
// unchanged and error names the offending word.
bool parse_spreads(const std::vector<std::string>& words, std::vector<double>& spreads, std::string& error);
}
#endif
