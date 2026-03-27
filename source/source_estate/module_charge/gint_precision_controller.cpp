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
    return PrecisionMode::double_mode;
}

void GintPrecisionController::reset_for_new_scf()
{
    switch (this->mode_)
    {
    case PrecisionMode::single:
    case PrecisionMode::mix:
        this->current_precision_ = ModuleGint::GintPrecision::fp32;
        this->locked_double_precision_ = false;
        break;
    case PrecisionMode::double_mode:
    default:
        this->current_precision_ = ModuleGint::GintPrecision::fp64;
        this->locked_double_precision_ = true;
        break;
    }
}

void GintPrecisionController::update_after_iteration(double drho, double scf_thr)
{
    if (this->locked_double_precision_ || this->mode_ != PrecisionMode::mix)
    {
        return;
    }

    const double switch_thr = std::max(1000.0 * scf_thr, 1.0e-5);
    if (drho <= switch_thr)
    {
        this->current_precision_ = ModuleGint::GintPrecision::fp64;
        this->locked_double_precision_ = true;
    }
}

ModuleGint::GintPrecision GintPrecisionController::current_precision() const
{
    return this->current_precision_;
}
