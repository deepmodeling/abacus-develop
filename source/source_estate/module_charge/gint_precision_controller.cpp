#include "gint_precision_controller.h"

#include <algorithm>

void GintPrecisionController::set_mode(const std::string& precision_mode)
{
    this->mode_ = parse_mode_(precision_mode);
}

GintPrecisionController::PrecisionMode GintPrecisionController::parse_mode_(const std::string& precision_mode)
{
    if (precision_mode == "single")
    {
        return PrecisionMode::single;
    }
    if (precision_mode == "mix")
    {
        return PrecisionMode::mix;
    }
    return PrecisionMode::double_precision;
}

void GintPrecisionController::reset_for_new_scf()
{
    switch (this->mode_)
    {
    case PrecisionMode::single:
    case PrecisionMode::mix:
        this->current_cfg_.cpu_internal_real = ModuleGint::GintRealPrecision::fp32;
        this->locked_fp64_ = false;
        break;
    case PrecisionMode::double_precision:
    default:
        this->current_cfg_.cpu_internal_real = ModuleGint::GintRealPrecision::fp64;
        this->locked_fp64_ = true;
        break;
    }
}

void GintPrecisionController::update_after_iteration(double drho, double scf_thr)
{
    if (this->locked_fp64_ || this->mode_ != PrecisionMode::mix)
    {
        return;
    }

    const double switch_thr = std::max(1000.0 * scf_thr, 1.0e-5);
    if (drho <= switch_thr)
    {
        this->current_cfg_.cpu_internal_real = ModuleGint::GintRealPrecision::fp64;
        this->locked_fp64_ = true;
    }
}

ModuleGint::GintExecConfig GintPrecisionController::current_config() const
{
    return this->current_cfg_;
}
