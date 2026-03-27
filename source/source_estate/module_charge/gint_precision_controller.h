#ifndef GINT_PRECISION_CONTROLLER_H
#define GINT_PRECISION_CONTROLLER_H

#include "source_lcao/module_gint/gint_precision.h"

#include <string>

class GintPrecisionController
{
  public:
    GintPrecisionController() = default;

    void set_mode(const std::string& precision_mode);

    void reset_for_new_scf();

    void update_after_iteration(double drho, double scf_thr);

    ModuleGint::GintRealPrecision current_precision() const;

  private:
    enum class PrecisionMode
    {
        single,
        double_precision,
        mix
    };

    static PrecisionMode parse_mode_(const std::string& precision_mode);

    ModuleGint::GintRealPrecision current_precision_ = ModuleGint::GintRealPrecision::fp64;
    PrecisionMode mode_ = PrecisionMode::double_precision;
    bool locked_fp64_ = true;
};

#endif
