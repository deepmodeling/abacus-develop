#ifndef GINT_PRECISION_CONTROLLER_H
#define GINT_PRECISION_CONTROLLER_H

#include "source_lcao/module_gint/gint_helper.h"

#include <string>

class GintPrecisionController
{
  public:
    GintPrecisionController() = default;

    void set_mode(const std::string& precision_mode);

    void reset_for_new_scf();

    /// Returns true if precision switched from fp32 to fp64 in this call.
    bool update_after_iteration(double drho, double scf_thr);

    ModuleGint::GintPrecision current_precision() const;

  private:
    enum class PrecisionMode
    {
        single,
        double_mode,
        mix
    };

    static PrecisionMode parse_mode_(const std::string& precision_mode);

    ModuleGint::GintPrecision current_precision_ = ModuleGint::GintPrecision::fp64;
    PrecisionMode mode_ = PrecisionMode::double_mode;
    bool locked_double_precision_ = true;
};

#endif
