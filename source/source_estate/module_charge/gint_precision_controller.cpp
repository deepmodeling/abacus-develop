#include "gint_precision_controller.h"

#include <algorithm>

void GintPrecisionController::reset_for_new_scf()
{
    this->current_cfg_.cpu_internal_real = ModuleGint::GintRealPrecision::fp32;
    this->qualified_small_drho_iters_ = 0;
    this->locked_fp64_ = false;
}

void GintPrecisionController::update_after_iteration(int iter,
                                                     double drho,
                                                     double scf_thr,
                                                     bool conv_esolver,
                                                     bool is_restart_step)
{
    (void)conv_esolver;

    if (this->locked_fp64_)
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
