#ifndef STERNHEIMER_RUNTIME_OPTIONS_H
#define STERNHEIMER_RUNTIME_OPTIONS_H

struct Input_para;

namespace ModuleRI
{

bool sternheimer_environment_flag(const char* name, bool default_value);
int sternheimer_channel_batch_width();

// Validate the declared isolated response, not a geometry inferred from Gamma sampling.
bool validate_sternheimer_molecular_coulomb(const Input_para& input);

} // namespace ModuleRI

#endif
