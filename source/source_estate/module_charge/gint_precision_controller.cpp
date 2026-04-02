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

bool GintPrecisionController::update_after_iteration(double drho, double scf_thr)
{
    if (this->locked_double_precision_ || this->mode_ != PrecisionMode::mix)
    {
        return false;
    }

    // Switch from fp32 to fp64 when drho is close enough to the target.
    // fp32 has ~7 significant digits (~1e-7 relative error), so we switch
    // well before that limit to let fp64 handle the final convergence.
    // The floor (kMinSwitchThreshold) prevents switching too early when
    // scf_thr is extremely tight.
    constexpr double kSwitchFactor = 1000.0;
    constexpr double kMinSwitchThreshold = 1.0e-5;
    const double switch_thr = std::max(kSwitchFactor * scf_thr, kMinSwitchThreshold);
    if (drho <= switch_thr)
    {
        this->current_precision_ = ModuleGint::GintPrecision::fp64;
        this->locked_double_precision_ = true;
        return true;
    }
    return false;
}

ModuleGint::GintPrecision GintPrecisionController::current_precision() const
{
    return this->current_precision_;
}
