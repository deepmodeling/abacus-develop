#include "gint_precision_controller.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>

namespace
{
bool read_forced_precision(ModuleGint::GintRealPrecision& precision)
{
    const char* env = std::getenv("ABACUS_GINT_FORCE_CPU_REAL");
    if (env == nullptr || env[0] == '\0' || std::strcmp(env, "auto") == 0)
    {
        return false;
    }
    if (std::strcmp(env, "fp32") == 0)
    {
        precision = ModuleGint::GintRealPrecision::fp32;
        return true;
    }
    if (std::strcmp(env, "fp64") == 0)
    {
        precision = ModuleGint::GintRealPrecision::fp64;
        return true;
    }
    return false;
}
}

void GintPrecisionController::apply_runtime_override_()
{
    this->force_precision_enabled_ = read_forced_precision(this->forced_precision_);
    if (this->force_precision_enabled_)
    {
        this->current_cfg_.cpu_internal_real = this->forced_precision_;
        this->locked_fp64_ = (this->forced_precision_ == ModuleGint::GintRealPrecision::fp64);
    }
}

void GintPrecisionController::reset_for_new_scf()
{
    this->current_cfg_.cpu_internal_real = ModuleGint::GintRealPrecision::fp32;
    this->qualified_small_drho_iters_ = 0;
    this->locked_fp64_ = false;
    this->apply_runtime_override_();
}

void GintPrecisionController::update_after_iteration(int iter,
                                                     double drho,
                                                     double scf_thr,
                                                     bool conv_esolver,
                                                     bool is_restart_step)
{
    (void)conv_esolver;

    if (this->force_precision_enabled_ || this->locked_fp64_)
    {
        return;
    }

    const double switch_thr = std::max(100.0 * scf_thr, 1.0e-5);
    const bool eligible = iter >= 3 && !is_restart_step && drho <= switch_thr;

    if (eligible)
    {
        ++this->qualified_small_drho_iters_;
    }
    else
    {
        this->qualified_small_drho_iters_ = 0;
    }

    if (this->qualified_small_drho_iters_ >= 2)
    {
        this->current_cfg_.cpu_internal_real = ModuleGint::GintRealPrecision::fp64;
        this->locked_fp64_ = true;
    }
}

ModuleGint::GintExecConfig GintPrecisionController::current_config() const
{
    return this->current_cfg_;
}
