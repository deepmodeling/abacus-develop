#ifndef STERNHEIMER_RUNTIME_OPTIONS_H
#define STERNHEIMER_RUNTIME_OPTIONS_H

struct Input_para;

namespace ModuleRI
{

enum class SternheimerWeakPreconditionerMode
{
    None,
    Spectral
};

bool sternheimer_environment_flag(const char* name, bool default_value);
int sternheimer_channel_batch_width();
SternheimerWeakPreconditionerMode sternheimer_weak_preconditioner_mode();
const char* sternheimer_weak_preconditioner_name(
    SternheimerWeakPreconditionerMode mode) noexcept;
double sternheimer_weak_residual_tolerance();

// Validate the declared isolated response, not a geometry inferred from Gamma sampling.
bool validate_sternheimer_molecular_coulomb(const Input_para& input);

} // namespace ModuleRI

#endif
