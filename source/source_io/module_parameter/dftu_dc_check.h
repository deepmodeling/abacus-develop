#ifndef DFTU_DC_CHECK_H
#define DFTU_DC_CHECK_H

#include <string>

struct Input_para;

namespace ModuleIO
{
// Called after INPUT resets; alpha_was_read distinguishes omission from 0.
std::string dftu_dc_error(const Input_para& input, bool alpha_was_read, int kpar_lcao);
double parse_dftu_alpha(const std::string& token);
} // namespace ModuleIO

#endif
